#!/usr/bin/env python3
"""Export unique ADP-binding/NBS overlap counts and corresponding GO FDR."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    nbs = {line.strip().split('.m')[0] for line in
           (ROOT / 'data/annotations/NBS_gene_ids.txt').read_text().splitlines()}
    annotations = pd.read_csv(ROOT / 'data/annotations/gene_GO.tsv.gz', sep='\t')
    adp = set(annotations.loc[annotations.ID == 'GO:0043531', 'gene_id'])
    rows = []
    for folder, contrasts in [('windows', ['fd_species', 'fd_sympatric']),
                              ('reference', ['fd_MY1', 'fd_CY32', 'fd_common'])]:
        base = ROOT / 'results' / folder
        enrichment = pd.read_csv(base / 'enrichment_all.tsv', sep='\t')
        for contrast in contrasts:
            for cutoff in [10, 5, 2.5, 1]:
                gene_file = base / f'{contrast}_{cutoff}_genes.txt'
                if gene_file.exists():
                    selected = set(gene_file.read_text().splitlines())
                    ids = selected & adp
                else:
                    # Retained enrichment tables contain the exact nonredundant hit IDs.
                    hit = enrichment[(enrichment.contrast == contrast) &
                                     (enrichment.top_percent == cutoff) &
                                     (enrichment.ID == 'GO:0043531') & (enrichment.database == 'GO')].iloc[0]
                    ids = set(str(hit.geneID).split('/')) if hit.Count else set()
                    assert ids.issubset(adp)
                hit = enrichment[(enrichment.contrast == contrast) &
                                 (enrichment.top_percent == cutoff) &
                                 (enrichment.ID == 'GO:0043531') & (enrichment.database == 'GO')].iloc[0]
                overlap = ids & nbs
                assert len(ids) == int(hit.Count)
                rows.append([contrast, cutoff, len(ids), len(overlap),
                             100 * len(overlap) / len(ids) if ids else None, hit.FDR_fixed_family])
    result = pd.DataFrame(rows, columns=['comparison', 'upper_tail_percent', 'ADP_binding_genes',
        'also_in_putative_NBS_list', 'NBS_percent_of_ADP_genes', 'ADP_GO_enrichment_FDR'])
    out = ROOT / 'results/tables'
    out.mkdir(parents=True, exist_ok=True)
    result.iloc[:8].to_csv(out / 'Table_S9.tsv', sep='\t', index=False)
    result.iloc[8:].to_csv(out / 'Table_S9_reference_sensitivity.tsv', sep='\t', index=False)
    print('Exported Table S9 and reference sensitivity table')


if __name__ == '__main__':
    main()
