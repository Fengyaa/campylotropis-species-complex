"""Usage: python scripts/prepare_svd_input.py input.nex output.nex
Number a sequential 259-taxon NEXUS alignment; append population partitions.
Requires the ORIGINAL taxon order; inspect the exported name/number map before PAUP.
"""
import re,sys
from pathlib import Path
if len(sys.argv)!=3:raise SystemExit(__doc__)
source,output=map(Path,sys.argv[1:]);text=source.read_text()
m=re.search(r'(?im)^\s*MATRIX\s*$',text)
if m is None:raise ValueError('No sequential MATRIX block')
end=text.index(';',m.end());rows=[x.split() for x in text[m.end():end].splitlines() if x.strip()]
if len(rows)!=259 or any(len(x)!=2 for x in rows):raise ValueError('Expected 259 single-line sequences')
canon=lambda x:re.sub(r'[^A-Za-z0-9]','',x).lower()
if canon(rows[130][0])!='xbdrc' or canon(rows[243][0])!='xubo1407':
 raise ValueError('Outgroup positions 131/244 disagree with alignment order')
part=Path(__file__).resolve().parents[1]/'data/phylogeny/taxpartitions.nex'
output.parent.mkdir(parents=True,exist_ok=True)
if output.exists():raise FileExistsError(output)
output.write_text(text[:m.end()]+'\n'+''.join(f'{i:03d} {r[1]}\n' for i,r in enumerate(rows,1))+text[end:]+'\n'+part.read_text())
output.with_suffix('.taxon_map.tsv').write_text('number\toriginal_name\n'+''.join(f'{i:03d}\t{r[0]}\n' for i,r in enumerate(rows,1)))
print('Prepared alignment; inspect taxon map before PAUP.')
