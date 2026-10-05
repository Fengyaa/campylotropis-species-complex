from pathlib import Path
import os
ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)
import csv
from pathlib import Path
import numpy as np
O=Path('results/local_affinity');meta=list(csv.DictReader(open(O/'samples.tsv'),delimiter='\t'));names=[r['sample'] for r in meta]
pops=sorted(set(r['population'] for r in meta[:40]));ix={p:np.array([i for i,r in enumerate(meta[:40]) if r['population']==p]) for p in pops}
group={r['population']:r['role'] for r in meta[:40]};cp=[j for j,p in enumerate(pops) if group[p]=='CPO'];cy=[j for j,p in enumerate(pops) if group[p]=='CYU'];rng=np.random.default_rng(20260926);rows=[]
for chrom,core in [('Chr8',47650001),('Chr6',36200001)]:
 z=np.load(Path('data/local_affinity')/(chrom+'.npz'));pos=z['pos'];gt=z['gt'];ok=np.ones(len(pos),bool)
 for ids in ix.values():ok&=(gt[:,ids]>=0).sum(axis=1)>=3
 windows=[('left_100kb',core-100000,core-1),('candidate_100kb',core,core+99999),('right_100kb',core+100000,core+199999)]
 if chrom=='Chr8':windows += [('MYB_gene',47707412,47710178),('association_core',47706408,47719857)]
 for label,start,end in windows:
  ids=np.flatnonzero(ok&(pos>=start)&(pos<=end));x=gt[ids];pp=pos[ids];dist=[]
  for p in pops:
   rr=x[:,ix[p]];valid=(rr[:,None,:]>=0)&(x[:,40:,None]>=0)
   d=np.abs(x[:,40:,None].astype(float)-rr[:,None,:])/2;d[~valid]=0
   avg=d.sum(2)/np.maximum(valid.sum(2),1);avg[x[:,40:]<0]=np.nan;dist.append(avg)
  ds=np.stack(dist,axis=2);bins=(pp-start)//1000;binids=[np.flatnonzero(bins==b) for b in sorted(set(bins))]
  # Equal weighting of occupied physical bins; not a confidence interval or independence claim.
  import warnings
  with warnings.catch_warnings():
   warnings.simplefilter('ignore',RuntimeWarning)
   equal=np.nanmean(np.stack([np.nanmean(ds[b],axis=0) for b in binids]),axis=0)
   draws=np.stack([np.nanmean(ds[[rng.choice(b) for b in binids]],axis=0) for _ in range(200)])
  for t in range(24):
   delta=draws[:,t,cy].mean(1)-draws[:,t,cp].mean(1)
   d=equal[t];pairs=[d[b]-d[a] for a in cp for b in cy];leave=[d[[j for j in cy if j!=omit]].mean()-d[[j for j in cp if j!=omit]].mean() for omit in range(10)]
   rows.append(dict(chrom=chrom,window=label,sample=names[40+t],population=meta[40+t]['population'],occupied_1kb_bins=len(binids),equal_bin_delta=d[cy].mean()-d[cp].mean(),equal_bin_pairs_CPO=sum(v>0 for v in pairs),equal_bin_leave_min=min(leave),equal_bin_leave_max=max(leave),random_draws=200,random_fraction_CPO=float(np.mean(delta>0)),random_delta_q025=np.quantile(delta,.025),random_delta_median=np.median(delta),random_delta_q975=np.quantile(delta,.975)))
with open(O/'bin_weighting_sensitivity.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
for chrom,label in dict.fromkeys((r['chrom'],r['window']) for r in rows):
 for pop in ['YS26','SM35','ML5']:
  rs=[r for r in rows if (r['chrom'],r['window'],r['population'])==(chrom,label,pop)]
  print(chrom,label,pop,'bins',rs[0]['occupied_1kb_bins'],'equal median',round(float(np.median([r['equal_bin_delta'] for r in rs])),4),'random CPO fraction range',min(r['random_fraction_CPO'] for r in rs),max(r['random_fraction_CPO'] for r in rs),'leave stable CPO',sum(r['equal_bin_leave_min']>0 for r in rs))
