"""Inventory only: deletion is performed by scoped native PowerShell after review."""
import hashlib
import json
import subprocess
from pathlib import Path

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
ROOT=BASE.parents[1]


def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    protected=json.loads((BASE/'budget_validation_20261026/summary/protected_before.json').read_text(encoding='utf8'))
    rows=[]
    for folder in sorted((BASE/'cache').iterdir()):
        if not folder.is_dir():continue
        manifest=folder/'cache_manifest.json'
        cfg=json.loads(manifest.read_text(encoding='utf8'))
        assert cfg['config']['holdout_enabled'] is False
        for path in sorted(folder.glob('*.npy')):
            resolved=path.resolve();assert resolved.is_relative_to((BASE/'cache').resolve())
            name=path.relative_to(ROOT).as_posix();assert name not in protected
            assert subprocess.run(['git','check-ignore','--quiet',name],cwd=ROOT).returncode==0
            assert not subprocess.check_output(['git','ls-files','--',name],cwd=ROOT)
            digest=sha(path);assert cfg['files'][path.name]==digest
            rows.append({'path':name,'absolute_path':str(resolved),'bytes':path.stat().st_size,
                         'sha256':digest,'action':'delete_recomputable_cache',
                         'manifest':manifest.relative_to(ROOT).as_posix()})
    assert len(rows)==84, f'Unexpected cache scope: {len(rows)}'
    result={'scope':'six historical species caches only; no source DNA, weights, per-run artifacts, logs, presentation or user changes',
            'files':rows,'file_count':len(rows),'total_bytes':sum(r['bytes'] for r in rows),
            'historical_protected_paths_intersect':0,'all_manifest_hashes_match':True,
            'restore':'fixed encoder/config and cache index in a NEW output directory, then compare original manifest hashes'}
    (HERE/'cleanup_plan.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print({'files':len(rows),'GiB':result['total_bytes']/1024**3,'status':'PLAN_VERIFIED'})


if __name__=='__main__':main()
