"""Final pRDA/BayPass overlap. Shared 10-kb windows need not share the same SNP.
Bundled pRDA candidates passed 3.5 SD and |r|>0.5; update this input after a new scan.
"""
import csv,gzip,re
from collections import Counter
from pathlib import Path
root=Path(__file__).resolve().parents[1];out=root/'results/gea';out.mkdir(parents=True,exist_ok=True)
with gzip.open(root/'data/gea/Z3.5_r0.5.candidates.tsv.gz','rt') as f:rows=list(csv.DictReader(f,delimiter='\t'))
S=set()
for r in rows:
 m=re.fullmatch(r'HiC_scaffold(\d+):(\d+)_[^_]+',r['snp'])
 if not m or abs(float(r['correlation']))<=.5:raise ValueError('Unexpected SNP/correlation')
 S.add((int(m[1]),int(m[2])))
if len(S)!=len(rows):raise ValueError('Duplicate pRDA SNPs')
with (root/'data/gea/baypass_candidates.tsv').open() as f:
 next(f);B={(int(c.replace('HiC_scaffold','')),int(p)) for c,p in (l.split() for l in f if l.strip())}
window=lambda x:(x[0],(x[1]-1)//10000*10000+1)
sc=Counter(map(window,S));bc=Counter(map(window,B));shared=S&B;wc=Counter(map(window,shared));W=set(sc)&set(bc)
with (out/'Z3.5_r0.5.shared_10kb_windows.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['chromosome','start','end','RDA_SNPs','BayPass_SNPs','exact_shared_SNPs'])
 for c,p in sorted(W):w.writerow([f'HiC_scaffold{c}',p,p+9999,sc[c,p],bc[c,p],wc[c,p]])
with (out/'overlap_summary.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['candidate_set','RDA_SNPs','RDA_windows','BayPass_SNPs','BayPass_windows','shared_SNPs','shared_windows'])
 w.writerow(['Z3.5_r0.5',len(S),len(sc),len(B),len(bc),len(shared),len(W)])
print('pRDA SNPs:',len(S),'Shared windows:',len(W),'Shared SNPs:',len(shared))
