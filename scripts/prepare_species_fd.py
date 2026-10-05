#!/usr/bin/env python3
"""Summarize supplied reference-based Dsuite estimates; no VCF rerun required."""
import bisect
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data/introgression'
WINDOWS = ROOT / 'data/windows/genomic_windows.tsv.gz'
OUT = ROOT / 'results/introgression'
CONFIGURATIONS = [('MY1_CYU_CPO', 'sp_set_Cyu.txt', 'MY1', 'CYU'),
                  ('CY32_CPO_CYU', 'sp_set_Cpo.txt', 'CY32', 'CPO')]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    windows = pd.read_csv(WINDOWS, sep='\t')
    assert not windows.duplicated(['chr', 'start', 'end']).any()
    lookup = {}
    for chrom, frame in windows.groupby('chr', sort=False):
        items = sorted((int(row.start), int(row.end), int(i))
                       for i, row in frame.iterrows())
        assert all(b[0] > a[1] for a, b in zip(items, items[1:]))
        lookup[chrom] = ([item[0] for item in items], items)
    population = dict(line.split() for line in
                      (ROOT / 'data/spatial/sample_population.tsv').read_text().splitlines())
    audit, manifest = {}, []
    for name, sets_name, ref, pool in CONFIGURATIONS:
        sets_path = SOURCE / sets_name
        records = [line.split() for line in sets_path.read_text().splitlines() if line.strip()]
        assert all(len(row) == 2 for row in records)
        assert len({row[0] for row in records}) == len(records)
        groups = defaultdict(list)
        for sample, group in records:
            groups[group].append(sample)
        assert set(groups['Outgroup']) == {'XB_DR_C', 'xubo1407'}
        assert len(groups[ref]) == 8 and all(population.get(s) == ref for s in groups[ref])
        assert not any(population.get(s) == ref for s in groups[pool])
        assert not any(population.get(s) == 'BM31' for s, _ in records)
        assert all(dict(records).get(s) == 'CPO' for s in ['Gln2014-8-1A', 'Gln2014-8-2A'])
        membership = pd.DataFrame([(g, population.get(s, 'OUTGROUP'), s) for s, g in records],
                                  columns=['group', 'population', 'sample'])
        membership.to_csv(OUT / f'{name}_samples.tsv', sep='\t', index=False)
        membership.groupby(['group', 'population']).size().rename('n').to_csv(
            OUT / f'{name}_population_counts.tsv', sep='\t')
        sums, counts = np.zeros(len(windows)), np.zeros(len(windows), dtype=int)
        qc = defaultdict(int)
        source = SOURCE / f'{name}_localFstats__50_25.txt.gz'
        with gzip.open(source, 'rt') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            assert {'chr', 'windowStart', 'windowEnd', 'f_d'}.issubset(reader.fieldnames)
            for row in reader:
                qc['input_rows'] += 1
                if row['chr'] == 'chr' and row['windowStart'] == 'windowStart' and row['f_d'] == 'f_d':
                    qc['repeated_header_rows'] += 1
                    continue
                if None in row or any(v is None for v in row.values()):
                    raise ValueError(f'Malformed row in {source.name}: {reader.line_num}')
                chrom, start, end, fd = row['chr'], int(row['windowStart']), int(row['windowEnd']), float(row['f_d'])
                assert chrom in lookup and end >= start >= 1
                if not np.isfinite(fd):
                    qc['nonfinite_fd'] += 1
                    continue
                if not 0 < fd < 1:
                    qc['outside_positive_unit_interval'] += 1
                    continue
                starts, items = lookup[chrom]
                j = bisect.bisect_right(starts, start) - 1
                if j < 0 or start > items[j][1]:
                    qc['unmatched_positive_rows'] += 1
                    continue
                i = items[j][2]
                sums[i] += fd
                counts[i] += 1
                qc['matched_positive_rows'] += 1
        values = np.divide(sums, counts, out=np.full(len(windows), np.nan), where=counts > 0)
        windows['fd_' + name] = values
        windows['n_raw_' + name] = counts
        qc['windows_with_positive_fd'] = int((counts > 0).sum())
        audit[name] = dict(qc, sets_group_counts={g: len(s) for g, s in groups.items()})
        for path in [source, sets_path]:
            manifest.append({'path': path.relative_to(ROOT).as_posix(),
                             'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    cols = ['fd_' + row[0] for row in CONFIGURATIONS]
    windows['n_species_fd_contrasts'] = windows[cols].notna().sum(axis=1)
    windows['species_trio_fd_recomputed'] = windows[cols].mean(axis=1)
    windows['species_fd_common_support'] = windows[cols].mean(axis=1).where(windows.n_species_fd_contrasts == 2)
    windows.to_csv(WINDOWS, sep='\t', index=False, na_rep='NA')
    paired = windows[(windows.n_species_fd_contrasts == 2) & (windows.sites >= 100)]
    audit['paired'] = {'n': len(paired), 'spearman': float(paired[cols].corr(method='spearman').iloc[0, 1]),
                       'support_counts': {str(k): int(v) for k, v in windows.n_species_fd_contrasts.value_counts().items()}}
    (OUT / 'input_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    pd.DataFrame(manifest).to_csv(OUT / 'input_manifest.tsv', sep='\t', index=False)
    print('Updated reference summaries in data/windows/genomic_windows.tsv.gz')


if __name__ == '__main__':
    main()
