from pathlib import Path
import os
ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)
"""Summarize retained FSC models and PSMC curves without guessing run settings."""
from pathlib import Path
import csv,json
out=Path('results/demography');out.mkdir(parents=True,exist_ok=True)
rows=[]
for p in sorted(Path('data/demography/fsc').glob('*.AIC')):
 with p.open() as f:r=next(csv.DictReader(f,delimiter='\t'))
 rows.append(dict(model=p.stem,**{k:float(v) for k,v in r.items()}))
rows.sort(key=lambda x:x['AIC'])
for r in rows:r['delta_AIC']=r['AIC']-rows[0]['AIC']
with (out/'FSC_retained_model_ranking.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
with Path('data/demography/fsc/model9.bestlhoods').open() as f:best=next(csv.DictReader(f,delimiter='\t'))
(out/'model9_parameters.json').write_text(json.dumps(best,indent=2)+'\n')
print('Summarized',len(rows),'retained FSC models. Best retained AIC:',rows[0]['model'])
