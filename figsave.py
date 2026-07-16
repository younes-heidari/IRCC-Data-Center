"""Single place where every report figure gets written.

WHY VECTOR: these figures are consumed by pdflatex, not by a screen. A 150 dpi
PNG scaled to \\linewidth resamples to roughly 200 dpi on the page -- soft text,
fuzzy hairlines, visible stair-stepping on the plot rules. A PDF is resolution-
independent: the text stays real text (selectable, searchable, sharp at any zoom)
and the vector strokes stay crisp on paper. It is also usually SMALLER than the
equivalent high-dpi raster.

Fonts are forced to TrueType (fonttype 42). Matplotlib's default Type 3 subsets
render acceptably but are rejected by several publishers' toolchains and copy out
of the PDF as garbage, so this avoids a class of problem for free.

A PNG companion is written alongside at `PNG_DPI` -- nothing in the report uses
it, but it keeps the figures previewable in a file browser or a slide deck, which
is what they get reused for.
"""

import os

import matplotlib

# Real TrueType fonts in the PDF, not Type 3 subsets. Must be set before any
# figure is saved; importing this module is what guarantees that ordering.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

PNG_DPI = 200


def save_figure(fig, out_path, facecolor=None, also_png=True):
    """Write `fig` as a vector PDF next to `out_path`, whatever extension the
    caller passed. Returns the PDF path actually written.

    Callers still name their figures `*.png` (that is what the report's
    \\includegraphics used to reference and what reads naturally at the call
    site); the extension is swapped here so there is exactly one rule about
    figure format and it lives in one file.
    """
    base = os.path.splitext(out_path)[0]
    pdf_path = base + ".pdf"

    kw = {}
    if facecolor is not None:
        kw["facecolor"] = facecolor

    fig.savefig(pdf_path, format="pdf", **kw)
    if also_png:
        fig.savefig(base + ".png", format="png", dpi=PNG_DPI, **kw)
    return pdf_path
