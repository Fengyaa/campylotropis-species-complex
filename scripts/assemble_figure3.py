#!/usr/bin/env python3
"""Combine generated Fig. 3 panels without rasterizing vector content."""
from pathlib import Path
from pypdf import PdfReader, PdfWriter, Transformation

ROOT = Path(__file__).resolve().parents[1]
panels = ['Fig3_species_matched_gene_NBS_180mm.pdf', 'Figure_3b.pdf', 'Figure_3c.pdf']


def main():
    pages = [PdfReader(ROOT / 'results/windows' / name).pages[0] for name in panels]
    width = float(pages[0].mediabox.width)
    heights = [float(p.mediabox.height) * width / float(p.mediabox.width) for p in pages]
    writer = PdfWriter()
    canvas = writer.add_blank_page(width=width, height=sum(heights))
    y = sum(heights)
    for page, height in zip(pages, heights):
        y -= height
        canvas.merge_transformed_page(page, Transformation().scale(
            width / float(page.mediabox.width)).translate(0, y))
    canvas.compress_content_streams()
    out = ROOT / 'results/figures'
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'Figure_3.pdf').open('wb') as handle:
        writer.write(handle)


if __name__ == '__main__':
    main()
