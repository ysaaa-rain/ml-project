"""Discovery-only fitting and fixed PWM sequence-level development evaluation."""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ALPHABET = 'ACGT'


def reverse_pwm(pwm):
    return pwm[::-1, ::-1]


def similarity(a, b, minimum_overlap=5):
    best = (-1., 0, '+')
    for strand, matrix in [('+',b), ('-',reverse_pwm(b))]:
        for offset in range(-len(matrix)+minimum_overlap, len(a)-minimum_overlap+1):
            left, right = max(0,offset), min(len(a), offset+len(matrix))
            x, y = a[left:right], matrix[left-offset:right-offset]
            if len(x) < minimum_overlap:
                continue
            # Center each column against uniform background.
            x, y = (x-.25).ravel(), (y-.25).ravel()
            denom = np.linalg.norm(x)*np.linalg.norm(y)
            score = float(x@y/denom) if denom else 0.
            if score > best[0]:
                best = score, offset, strand
    return best


def make_pwm(sequences, pseudocount=.5):
    if not sequences or len({len(s) for s in sequences}) != 1:
        raise ValueError('Aligned equal-length instances required.')
    counts = np.full((len(sequences[0]),4), pseudocount)
    for seq in sequences:
        for i, char in enumerate(seq):
            counts[i, ALPHABET.index(char)] += 1
    return counts, counts/counts.sum(axis=1,keepdims=True)


def max_scan(sequences, pwm, background):
    logodds = np.log2(pwm/background)
    width = len(pwm)
    result = []
    for seq in sequences:
        encoded = np.array([ALPHABET.index(c) for c in seq])
        forward = np.array([logodds[np.arange(width),encoded[i:i+width]].sum()
                            for i in range(len(seq)-width+1)])
        rc = np.log2(reverse_pwm(pwm)/background)
        reverse = np.array([rc[np.arange(width),encoded[i:i+width]].sum()
                            for i in range(len(seq)-width+1)])
        both = np.stack([forward,reverse])
        strand, start = np.unravel_index(np.argmax(both),both.shape)
        result.append((float(both[strand,start]),int(start), '+' if strand == 0 else '-'))
    return result


def background(sequences):
    counts = np.array([sum(s.count(c) for s in sequences) for c in ALPHABET],float)+.5
    # Strand-symmetric zero-order background for scanning both strands.
    counts = (counts+counts[::-1])/2
    return counts/counts.sum()


def bh(values):
    values = np.asarray(values,float)
    if not len(values):
        return values
    order = np.argsort(values)
    q = np.minimum.accumulate((values[order]*len(values)/np.arange(1,len(values)+1))[::-1])[::-1]
    result = np.empty_like(values)
    result[order] = np.minimum(q,1)
    return result


def group_presence(frame, presence):
    return pd.DataFrame({'group':frame.leakage_group.tolist(), 'present':presence}).groupby('group').present.max()


def test_enrichment(positive, negative):
    # Operational independent units are frozen leakage groups, not windows.
    p, n = np.asarray(positive,bool), np.asarray(negative,bool)
    table = [[int(p.sum()),int((~p).sum())],[int(n.sum()),int((~n).sum())]]
    odds, value = fisher_exact(table, alternative='greater')
    corrected = np.array(table,float)+.5
    log_or = np.log(corrected[0,0]*corrected[1,1]/(corrected[0,1]*corrected[1,0]))
    se = np.sqrt((1/corrected).sum())
    return {'positive_groups':len(p),'control_groups':len(n),'positive_hits':int(p.sum()),
            'control_hits':int(n.sum()), 'odds_ratio':float(odds), 'p':float(value),
            'or_ci_low':float(np.exp(log_or-1.96*se)), 'or_ci_high':float(np.exp(log_or+1.96*se))}


def write_meme(path, motifs, bg):
    lines = ['MEME version 4','','ALPHABET= ACGT','','strands: + -','',
             'Background letter frequencies', ' '.join(f'{c} {v:.8f}' for c,v in zip(ALPHABET,bg)),'']
    for name,pwm in motifs.items():
        lines += [f'MOTIF {name}', f'letter-probability matrix: alength= 4 w= {len(pwm)}']
        lines += [' '.join(f'{v:.9f}' for v in row) for row in pwm]
        lines += ['']
    path.write_text('\n'.join(lines),encoding='utf-8')


def read_meme(path):
    import re
    lines = path.read_text().splitlines()
    motifs = {}
    name = None
    for i,line in enumerate(lines):
        if line.startswith('MOTIF '):
            name = line.split()[1]
        if 'letter-probability matrix:' in line:
            width = int(re.search(r'w=\s*(\d+)',line).group(1))
            motifs[name] = np.array([[float(v) for v in row.split()[:4]] for row in lines[i+1:i+width+1]])
    return motifs
