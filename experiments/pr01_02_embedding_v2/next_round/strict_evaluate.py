import json
import pandas as pd
from ..audit import save
from ..discovery import read_meme,background
from .benchmark import evaluate
from .common import NEXT,start,finish


def main():
    inputs=list((NEXT/'synthetic_benchmark').glob('*/data.json'))
    inputs+=list((NEXT/'synthetic_benchmark').glob('*/*/motifs.meme'))
    target,manifest=start('synthetic_strict_evaluation',json.loads((NEXT/'configs/strict_evaluation.json').read_text()),inputs)
    all_metrics=[];outputs=[]
    for path in sorted((NEXT/'synthetic_benchmark').glob('*/data.json')):
        data=json.loads(path.read_text(encoding='utf-8'));scenario=data['scenario'];seed=data['seed'];bg=background(data['calibration_control'])
        for run in sorted(path.parent.glob('*__w*')):
            method,width=run.name.split('__w');width=int(width)
            metrics,preds,matches,grammar=evaluate(read_meme(run/'motifs.meme'),scenario,data,width,bg,strict=True)
            for row in metrics:row.update(scenario=scenario['id'],seed=seed,method=method,width=width)
            all_metrics+=metrics
            dest=run/'strict';dest.mkdir(exist_ok=False)
            for name,records in [('metrics',metrics),('predictions',preds),('reference_matches',matches),('grammar',grammar)]:
                if records:pd.DataFrame(records).to_csv(dest/f'{name}.tsv',sep='\t',index=False)
            outputs+=list(dest.glob('*.tsv'))
    result=NEXT/'synthetic_benchmark/strict_metrics.tsv'
    pd.DataFrame(all_metrics).to_csv(result,sep='\t',index=False);outputs.append(result)
    finish(target,manifest,outputs,{'fits_reevaluated':132,'no_model_refit':True,'original_results_preserved':True})
    print('Strict synthetic evaluation completed')


if __name__=='__main__':main()
