//Number of population samples (demes)
3
//Population effective sizes (number of genes)
NCUR1
NCUR2
NCUR3
//Sample sizes
100
100
100
//Growth rates  : negative growth implies population expansion
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 M21 M31
M21 0 M32
M31 M32 0
//Migration matrix 1
0 0 0
0 0 0
0 0 0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 
2  historical event 
TDIV1 1 2 1 RES1 0 1 
TDIV2 0 2 1 RES2 0 1
/Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 7e-9 OUTEXP
