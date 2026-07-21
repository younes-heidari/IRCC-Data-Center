# Presentation assets

Everything needed to build slides, in one place. Regenerated from the code —
`report/presentation/` is a flat copy of the report figures plus the full report.

- **`Project_Report_full.pdf`** — the compiled 145-page report
- **`Project_Report_full.txt`** — the LaTeX source
- **`figures/`** — every report figure, flat, as both **`.png`** (drop straight
  into PowerPoint) and **`.pdf`** (vector, for print/large screens)

## Figure index

| File | What it shows |
|---|---|
| `compliance` | COP 3.53 / IPLV 5.23 / PUE 1.29 vs the mandatory targets |
| `ph_chart` | R-290 pressure–enthalpy chart of the design cycle |
| `temperature_cascade` | The 10 K air-to-refrigerant approach as three stacked steps |
| `hx_profiles` | Condenser (3-zone) and evaporator (2-zone) T–Q profiles |
| `cop_tc_sensitivity` | Full-load COP vs condensing temperature (why 45 °C matters) |
| `full_load_cop` | Full-load COP at AHRI vs the project design point |
| `iplv_points` | The four IPLV part-load points (A/B/C/D) |
| `pue_breakdown` | Design-day facility power split (compressor 75 %) |
| `economizer_band` | Load split + compressor power across the 4–10 °C partial band |
| `monthly_cop_pue` | Monthly COP / PUE / TUE across the Champaign year |
| `seasonal_cop_pue` | Seasonal COP / PUE / TUE |
| `annual_cop_load` | Delivered load + compressor COP across all 8,760 hours |
| `ambient_duration_curve` | Champaign TMY3 duration curve vs economizer thresholds |
| `refrigerant_comparison_cop` | R-290 vs R-1234yf vs R-1234ze(E) screening COP |
| `charge_inventory` | R-290 charge by component (5.5–8.7 kg total) |
| `charge_concentration` | Room concentration vs volume against EN 378 limits |
| `eev_authority` | EEV opening across the 5:1 turndown (one vs two valves) |
| `superheat_sensitivity` | COP vs superheat setpoint (the EN12900 artifact) |

> Note: `figures/pid/` is intentionally empty in the report tree — the P&ID
> (`fig:setup`) is the one deliverable still to be drawn.
