# Presentation plan — IRCC 2026 Project #3

Working area for the slide deck. `Project_Report_full.pdf` and all 18 figures
(`figures/*.png`) are here. This file is the outline we iterate on.

## The one-sentence story
A three-loop R-290 chiller that keeps flammable propane out of the data hall,
meets both efficiency targets (COP 3.53, IPLV 5.23), and rides Champaign's
winters to do **40 % of the year's cooling with almost no compressor work**.

## Proposed slide flow (draft — to be reshaped around your plan)

| # | Slide | Lead visual | The point to land |
|---|---|---|---|
| 1 | Title | — | Project, team, one-line pitch |
| 2 | The problem | (photo/diagram) | 150 kW edge data center, 24/7, Champaign; A3 refrigerant can't enter the hall |
| 3 | System architecture | *P&ID (to draw)* | Three loops: water → propane → glycol → air; why three |
| 4 | Working fluid choice | `refrigerant_comparison_cop.png` | R-290: GWP 3, competitive COP; why not R-1234ze |
| 5 | The design cycle | `ph_chart.png` | p–h chart, the four state points, COP 4.67 at design |
| 6 | Why the dry cooler is huge | `temperature_cascade.png` | 10 K approach = 3 stacked steps; the price of the glycol loop |
| 7 | Meeting the targets | `compliance.png` | COP 3.53 / IPLV 5.23 / PUE 1.29 vs the mandates |
| 8 | Part-load & turndown | `iplv_points.png` + `eev_authority.png` | VFDs give the 5:1 turndown; single EEV stays controllable |
| 9 | Free cooling — the payoff | `economizer_band.png` | Load split across the economizer band; compressor collapses |
| 10 | Seasonal performance | `seasonal_cop_pue.png` or `monthly_cop_pue.png` | 40 % of cooling free; winter PUE 1.01 |
| 11 | Safety (A3) | `charge_inventory.png` + `charge_concentration.png` | 5.5–8.7 kg charge, no receiver, EN 378 ventilation |
| 12 | What's left / next steps | — | P&ID, commissioning, spec sheets, pump trim |
| 13 | Summary | `compliance.png` | Targets met, honest caveats, key design decisions |

## Figure bench (all in `figures/`)
compliance · ph_chart · temperature_cascade · hx_profiles · cop_tc_sensitivity ·
full_load_cop · iplv_points · pue_breakdown · economizer_band · monthly_cop_pue ·
seasonal_cop_pue · annual_cop_load · ambient_duration_curve ·
refrigerant_comparison_cop · charge_inventory · charge_concentration ·
eev_authority · superheat_sensitivity

## Open questions for you
- How long is the talk / how many slides is the target?
- Audience — the professor + class (technical), or broader?
- Do you want me to draft speaker notes per slide, build an HTML/PDF deck, or
  just organize the assets and outline for you to drop into PowerPoint?
- Anything in your plan I should build the flow around?

## TODO
- [ ] Confirm slide count & audience
- [ ] Lock the slide flow
- [ ] P&ID figure (the one missing visual)
- [ ] Per-slide speaker notes
