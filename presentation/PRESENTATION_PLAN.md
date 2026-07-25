# Presentation plan — IRCC 2026 Project #3 (LOCKED · built)

~25-minute talk · 20 slides · built as `Presentation Main.html` (interactive) + `Presentation Main.pdf` (print).
Evaporator ("two-zone model") and condenser ("three-zone model") slides show the
structurally different ΔP models; Martin's full form is on the evaporator slide only.

## The one-sentence story
A three-loop R-290 chiller that keeps flammable propane out of the data hall,
meets both efficiency targets (COP 3.53, IPLV 5.23), and rides Champaign's
winters to do **40 % of the year's cooling with almost no compressor work**.

## Slide flow (as built)

| # | Slide | Figure(s) |
|---|---|---|
| 1 | Title | — |
| 2 | Why data-center cooling matters (stat tiles) | `Component Figures/Datacenter.jpg` (credit: gridx.com) |
| 3 | The design task — deliverables + hard targets | — |
| 4 | Proposed design | `Component Figures/P&ID.svg` |
| 5 | Working fluid — why propane | `figures/refrigerant_comparison_cop.png` |
| 6–13 | **Component slides** — modeling + equations LEFT · P&ID spotlight + vendor photo RIGHT | see below |
| 14 | Results — full-load COP 3.53 | `figures/full_load_cop.png` |
| 15 | Results — IPLV 5.23 & 5:1 turndown | `figures/iplv_points.png` |
| 16 | Results — free cooling & seasonal | `economizer_band` + `seasonal_cop_pue` |
| 17 | Results — PUE 1.29 & compliance | `figures/compliance.png` |
| 18 | Safety — A3 charge analysis | `charge_inventory` + `charge_concentration` |
| 19 | Summary | `figures/ph_chart.png` |
| 20 | Thank you / acknowledgements | headshots + sponsor logos |

### Component slides (6–13)
Each: LEFT = modeling method + key correlations/equations; RIGHT = full P&ID
blurred except the component (orange ring) + real vendor photo.

| Slide | Component | Vendor / model | Key equation shown |
|---|---|---|---|
| 6 | Evaporator | Kelvion HP DW 500H, ≈84 plates | LMTD sizing + Han (2003) f, ΔP |
| 7 | Condenser | Kelvion HP DW 500H, ≈104 plates | 3-zone LMTD + Han condensation |
| 8 | Compressors | 2× Bitzer 4FEP-35Z + VFD | AHRI-540 10-coeff polynomial |
| 9 | EEV | Danfoss ETS / Emerson EX class | orifice ṁ = C_d·A·√(2ρΔP) |
| 10 | Oil separator | Bitzer OA1954 | suction volume flow 203.6 m³/h |
| 11 | Dry cooler | Kelvion ULF-PA106K4V-091F095 | Q̇ = ṁc_pΔT; 3+4+3 K stack |
| 12 | Economizer | Kelvion GB, ≈49 plates | A = Q̇/(U·ΔT_lm) at 4 °C OAT |
| 13 | Pumps (both) | B&G e-1510 2.5AC + 1.5AD | P = ρgQH/η |

## Rebuild the PDF after edits
```
msedge --headless --disable-gpu --print-to-pdf="<abs path>\Presentation Main.pdf" ^
  --no-pdf-header-footer --virtual-time-budget=15000 "file:///<abs path>/Presentation Main.html"
```
(Print CSS in Presentation Main.html gives one 1280×720 page per slide.)

## Notes / open items
- P&ID updated: receiver removed (matches the report's no-receiver architecture),
  economizer labeled, TCE valve repositioned. Spotlight coords re-verified.
- HX slides (evap/cond/econ) show the design flow: duty → ΔT_lm → assumed U →
  A = Q̇/(U·ΔT_lm) → plate count, plus Han (two-phase) & Martin (single-phase) ΔP.
- Title/thank-you credit "Younes Heidari, Christian Muller".
- `Glycol Pump.jpg` is 468×333 — fine at grid size, don't enlarge.
- CRAH is not vendor-selected (placeholder in report) — named in "next steps".
