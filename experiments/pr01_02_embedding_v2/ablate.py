"""Predeclared single-species diagnostic ablations, not six-species claims."""
import json
from .audit import BRANCH
from .encoder import load_model
from .run import run


def main():
    cfg = json.loads((BRANCH/'config.json').read_text(encoding='utf-8'))
    model,tokenizer = load_model(cfg)
    source = 'tjupan_bacillus_subtilis'
    for layer,width in [(3,10),(6,6),(6,8),(6,12),(6,16)]:
        run(source,20261005,'prokbert',model,tokenizer,layer,width)
    for seed in [20261005,20261006,20261007]:
        run(source,seed,'multiscale',model,tokenizer,6,16)
    run(source,20261005,'onehot',model,tokenizer,6,16)


if __name__ == '__main__':
    main()
