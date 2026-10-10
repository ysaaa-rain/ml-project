"""Frozen 6-mer encoder with explicit, tested half-open nucleotide spans."""
import json
from pathlib import Path
import numpy as np
from .audit import BRANCH, sha, save


def spans(length, k=6):
    if length < k:
        raise ValueError('Sequence shorter than token.')
    return [(i, i+k) for i in range(length-k+1)]


def token_to_base(hidden, length=81, k=6):
    if hidden.shape[-2] != length-k+1:
        raise ValueError('Special tokens must be excluded.')
    base = np.zeros((*hidden.shape[:-2], length, hidden.shape[-1]), dtype=np.float32)
    coverage = np.zeros(length, dtype=np.float32)
    for i, (start, end) in enumerate(spans(length, k)):
        base[..., start:end, :] += hidden[..., i:i+1, :]
        coverage[start:end] += 1
    return base / coverage[:, None]


def windows(base, width):
    if not 1 <= width <= base.shape[-2]:
        raise ValueError('Invalid window width.')
    return np.stack([base[..., i:i+width, :].mean(axis=-2)
                     for i in range(base.shape[-2]-width+1)], axis=-2)


def load_model(config):
    import torch
    import requests
    import importlib.util
    import sys
    torch.set_num_threads(config['threads'])
    cache = BRANCH / 'model_cache' / 'pinned'
    cache.mkdir(parents=True,exist_ok=True)
    revision = config['revision']
    code_revision = config['code_revision']
    files = {'config.json':(config['model'],revision), 'vocab.txt':(config['model'],revision),
             'model.safetensors':(config['model'],revision),
             'models.py':('neuralbioinfo/nbrg-transformers',code_revision),
             'tokenizer.py':('neuralbioinfo/nbrg-transformers',code_revision)}
    for name,(repository,rev) in files.items():
        target = cache/name
        if not target.exists():
            url = f'https://huggingface.co/{repository}/resolve/{rev}/{name}'
            print(f'downloading pinned {name}',flush=True)
            if name == 'model.safetensors':
                # Bounded HTTP ranges avoid stalled large Xet transfers on this Windows network.
                from concurrent.futures import ThreadPoolExecutor
                chunk_size = 1024*1024
                first = requests.get(url,headers={'Range':f'bytes=0-{chunk_size-1}'},timeout=(20,60))
                first.raise_for_status()
                if first.status_code != 206:
                    target.write_bytes(first.content)
                    continue
                size = int(first.headers['Content-Range'].split('/')[-1])
                def get_range(start):
                    end = min(size-1,start+chunk_size-1)
                    # Cache buster avoids a proxy serving a different cached Range.
                    for attempt in range(3):
                        try:
                            r = requests.get(url+f'?range_start={start}',headers={'Range':f'bytes={start}-{end}'},timeout=(20,60))
                            r.raise_for_status()
                            if r.status_code != 206 or r.headers.get('Content-Range') != f'bytes {start}-{end}/{size}' or len(r.content)!=end-start+1:
                                raise ValueError('Incorrect HTTP range response')
                            return r.content
                        except Exception:
                            if attempt == 2:
                                raise
                temporary = target.with_suffix('.safetensors.part')
                with temporary.open('wb') as stream, ThreadPoolExecutor(max_workers=4) as pool:
                    stream.write(first.content)
                    for idx,chunk in enumerate(pool.map(get_range,range(chunk_size,size,chunk_size)),1):
                        stream.write(chunk)
                        if idx%10 == 0:
                            stream.flush()
                            print(f'weights downloaded {idx+1} MiB',flush=True)
                if temporary.stat().st_size != size:
                    raise ValueError('Weight size mismatch')
                temporary.replace(target)
                continue
            with requests.get(url,stream=True,timeout=(20,60)) as response:
                response.raise_for_status()
                temporary = target.with_suffix(target.suffix+'.part')
                with temporary.open('wb') as stream:
                    for chunk in response.iter_content(1024*1024):
                        stream.write(chunk)
                temporary.replace(target)
    if sha(cache/'model.safetensors') != config['weights_sha256']:
        raise ValueError('Weights do not match pinned Hugging Face LFS SHA256.')
    # Import only the explicitly pinned, locally recorded official implementation.
    modules = {}
    for name in ['models','tokenizer']:
        spec = importlib.util.spec_from_file_location('emb_pinned_'+name,cache/(name+'.py'))
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        modules[name] = module
    tokenizer = modules['tokenizer'].LCATokenizer(vocab_file=str(cache/'vocab.txt'),kmer=6,shift=1)
    model_config = modules['models'].ProkBertConfig.from_pretrained(str(cache),local_files_only=True)
    model, info = modules['models'].ProkBertModel.from_pretrained(str(cache),config=model_config,
                                      local_files_only=True,output_loading_info=True)
    if info.get('missing_keys') or info.get('mismatched_keys'):
        raise RuntimeError(f'Uninitialized encoder weights: {info}')
    model.eval()
    model.requires_grad_(False)
    save(BRANCH / 'audit/model.json', {'model': config['model'], 'revision': revision,
         'code_revision': code_revision, 'loading': info,
         'files': {p.name: sha(p) for p in cache.iterdir() if p.is_file()},
         'parameters': sum(p.numel() for p in model.parameters()), 'device': 'cpu'})
    return model, tokenizer


def extract(frame, out, config, model, tokenizer):
    import torch
    out = Path(out)
    out.mkdir(exist_ok=True, parents=True)
    layers = config['layers']
    token_arrays = {l: np.lib.format.open_memmap(out / f'token_L{l}.npy', mode='w+',
                    dtype='float32', shape=(len(frame), 76, model.config.hidden_size)) for l in layers}
    batch_size = config['batch_size']
    with torch.inference_mode():
        for start in range(0, len(frame), batch_size):
            seqs = frame.sequence.iloc[start:start+batch_size].tolist()
            inputs = tokenizer(seqs, return_tensors='pt', padding=True)
            for row, seq in zip(inputs['input_ids'].tolist(), seqs):
                expected = [tokenizer.convert_tokens_to_ids(seq[i:i+6]) for i in range(76)]
                if row[1:-1] != expected or len(row) != 78 or tokenizer.unk_token_id in expected:
                    raise ValueError('Tokenizer does not match explicit 6-mer spans.')
            result = model(**inputs, output_hidden_states=True)
            for layer in layers:
                token_arrays[layer][start:start+len(seqs)] = result.hidden_states[layer][:,1:-1].numpy()
            if start % (batch_size*10) == 0:
                print(f'encoded {start}/{len(frame)}', flush=True)
    for layer, array in token_arrays.items():
        array.flush()
        base = token_to_base(array)
        np.save(out / f'base_L{layer}.npy', base)
        for width in config['windows']:
            np.save(out / f'window_L{layer}_W{width}.npy', windows(base, width))
    frame.to_csv(out / 'index.tsv', sep='\t', index=False)
    save(out / 'cache_manifest.json', {'ordered_ids': frame.sequence_id.tolist(),
          'config': config, 'files': {p.name: sha(p) for p in out.iterdir() if p.is_file()},
          'token_span': '[i,i+6), 0-based; special tokens excluded',
          'base_rule': 'average covering contextual token vectors; includes flanking sequence',
          'tss': {'one_based':61, 'zero_based':60}, 'holdout_used': False})
