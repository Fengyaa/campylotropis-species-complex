from pathlib import Path
import os
ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)
import csv,gzip,json
from pathlib import Path
import numpy as np
O=Path('results/local_affinity');O.mkdir(parents=True,exist_ok=True)
def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def save(p,rs):
 if not rs:return
 with open(O/p,'w') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t');w.writeheader();w.writerows(rs)
refs=read('data/metadata/balanced_panel_40.tsv')
meta=read('data/metadata/individual_Q.tsv')
targets=[r for r in meta if r['population'] in ['ML5','YS26','SM35']]
names=[r['sample'] for r in refs+targets];assert len(names)==len(set(names))
pops=[r['population'] for r in refs+targets];popnames=sorted(set(pops[:40]));groups={r['population']:r['reference_group'] for r in refs}
ix={p:np.array([i for i,x in enumerate(pops[:40]) if x==p]) for p in popnames}
save('samples.tsv',[dict(sample=n,population=pops[i],role=groups.get(pops[i],'target')) for i,n in enumerate(names)])
summary=[];pairrows=[];audit=[]
for chrom,core in [('Chr8',47650001),('Chr6',36200001)]:
 cache=Path('data/local_affinity')/(chrom+'.npz')
 z=np.load(cache);pos=z['pos'];gt=z['gt'];assert list(z['names'])==names
 # Require >=3/4 callable in EVERY reference population to keep comparison loci symmetric.
 ok=np.ones(len(pos),bool)
 for ids in ix.values():ok &= (gt[:,ids]>=0).sum(axis=1)>=3
 audit.append(dict(chrom=chrom,extracted_sites=len(pos),common_reference_QC_sites=int(ok.sum())))
 windows=[('left_100kb',core-100000,core-1),('candidate_100kb',core,core+99999),('right_100kb',core+100000,core+199999)]
 if chrom=='Chr8':windows += [('MYB_gene',47707412,47710178),('association_core',47706408,47719857)]
 windows += [('scan_'+str(j),core+j*50000,core+(j+1)*50000-1) for j in range(-10,12)]
 for label,start,end in windows:
  valid=np.flatnonzero(ok & (pos>=start)&(pos<=end))
  for mode in ['all_QC','thin1kb']:
   ids=valid
   if mode=='thin1kb':
    take=[];last=-10**12
    for k in ids:
     if pos[k]-last>=1000:take.append(k);last=pos[k]
    ids=np.array(take,dtype=int)
   x=gt[ids];pp=pos[ids]
   if not len(x):continue
   # Mean absolute allele dosage difference /2 = unphased IBS dissimilarity.
   perpop=[]
   for pop in popnames:
    rr=x[:,ix[pop]]
    d=np.abs(x[:,:,None].astype(float)-rr[:,None,:])/2
    good=(x[:,:,None]>=0)&(rr[:,None,:]>=0)
    d[~good]=0
    perpop.append(d.sum(axis=2)/np.maximum(good.sum(axis=2),1))
   ds=np.stack(perpop,axis=2)
   cp=[j for j,p in enumerate(popnames) if groups[p]=='CPO'];cy=[j for j,p in enumerate(popnames) if groups[p]=='CYU']
   for t in range(40,len(names)):
    called=x[:,t]>=0;n=int(called.sum())
    if not n:continue
    d=ds[called,t,:].mean(axis=0)
    delta=float(d[cy].mean()-d[cp].mean())
    leave=[float(d[[j for j in cy if j!=omit]].mean()-d[[j for j in cp if j!=omit]].mean()) for omit in range(len(popnames))]
    contrasts=[float(d[b]-d[a]) for a in cp for b in cy]
    summary.append(dict(chrom=chrom,window=label,start=start,end=end,mode=mode,sample=names[t],population=pops[t],retained_sites=len(x),called_sites=n,call_rate=n/len(x),distance_CPO=float(d[cp].mean()),distance_CYU=float(d[cy].mean()),delta_CYU_minus_CPO=delta,leave_one_pop_min=min(leave),leave_one_pop_max=max(leave),pair_CPO_closer=sum(v>0 for v in contrasts),pair_CYU_closer=sum(v<0 for v in contrasts),pair_ties=sum(v==0 for v in contrasts)))
    if not label.startswith('scan'):
     for a in cp:
      for b in cy:pairrows.append(dict(chrom=chrom,window=label,mode=mode,sample=names[t],population=pops[t],CPO_reference=popnames[a],CYU_reference=popnames[b],delta=float(d[b]-d[a]),called_sites=n))
   if mode=='all_QC' and label in ['left_100kb','candidate_100kb','right_100kb']:
    matrix=np.zeros((len(names),len(names)));counts=np.zeros_like(matrix)
    for i in range(len(names)):
     good=(x[:,i,None]>=0)&(x>=0);counts[i]=good.sum(axis=0)
     matrix[i]=(np.abs(x[:,i,None].astype(float)-x)*good).sum(axis=0)/np.maximum(counts[i],1)/2
    for tag,mat in [('distance',matrix),('n_shared',counts)]:
     with open(O/(chrom+'_'+label+'_'+tag+'.tsv'),'w') as f:
      f.write('sample\t'+'\t'.join(names)+'\n')
      for name,row in zip(names,mat):f.write(name+'\t'+'\t'.join(map(str,row))+'\n')
 print(chrom,'done',flush=True)
save('individual_affinity.tsv',summary);save('reference_pair_sensitivity.tsv',pairrows);save('QC_audit.tsv',audit)
agg=[]
for chrom in ['Chr8','Chr6']:
 for label in sorted(set(r['window'] for r in summary if r['chrom']==chrom)):
  for mode in ['all_QC','thin1kb']:
   for pop in ['YS26','SM35','ML5']:
    rs=[r for r in summary if (r['chrom'],r['window'],r['mode'],r['population'])==(chrom,label,mode,pop)]
    if rs:agg.append(dict(chrom=chrom,window=label,start=rs[0]['start'],end=rs[0]['end'],mode=mode,population=pop,n_individuals=len(rs),min_called_sites=min(r['called_sites'] for r in rs),median_delta=float(np.median([r['delta_CYU_minus_CPO'] for r in rs])),min_delta=min(r['delta_CYU_minus_CPO'] for r in rs),max_delta=max(r['delta_CYU_minus_CPO'] for r in rs),n_CPO_closer=sum(r['delta_CYU_minus_CPO']>0 for r in rs),n_all25_pairs_CPO=sum(r['pair_CPO_closer']==25 for r in rs),n_all25_pairs_CYU=sum(r['pair_CYU_closer']==25 for r in rs)))
save('population_summary.tsv',agg)
