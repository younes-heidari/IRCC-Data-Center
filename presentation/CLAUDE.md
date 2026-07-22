# Presentation build

Paths are relative to the repository root.

## Inputs
- Design system: presentation/DESIGN_SYSTEM.md    (authoritative — never invent values)
- Slide template: presentation/Template.html       (copy its structure exactly)
- Slide plan: presentation/PRESENTATION_PLAN.md    (slide-by-slide outline)
- Content source: report/Project_Report_full.tex   (also .txt export in report/)
- Figures: presentation/figures/                   (insert existing files, don't regenerate)
- Output: presentation/

## Rules
- Every color, font, and spacing value comes from DESIGN_SYSTEM.md.
- Reuse Template.html's slide structure; add slides, don't restyle.
- Plan slide-by-slide and get approval before generating.
