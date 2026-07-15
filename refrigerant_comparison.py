"""Working-fluid screening comparison (project spec subtask 1c: "Compare
different working fluid and select the most promising").

Scope: every candidate must satisfy the project's own hard constraints
(GWP < 150, ASHRAE 34 safety class A1/A2L/A3 only) and be able to run a
standard SUBCRITICAL vapor-compression cycle at this project's design
temperatures (7 C evaporating, 52 C condensing) -- which rules out CO2
(R-744, critical temp 31.1 C, would need a transcritical cycle, a
fundamentally different architecture) and ammonia (R-717, safety class
B2L -- toxic, fails the project's A1/A2L/A3-only requirement) before any
performance comparison is even run.

This is a SCREENING-level comparison, not the detailed AHRI-540 vendor-map
compressor model used elsewhere in this project (compressor.py) -- that
map only exists for the selected refrigerant/compressor pair (R-290,
Bitzer 4FEP-35Z). To compare refrigerants on an apples-to-apples basis
before a compressor is even selected, all three candidates are run through
the SAME idealized single-stage cycle (fixed isentropic efficiency,
identical superheat/subcooling/design temperatures) so any COP difference
reflects the refrigerant's own thermodynamics, not a compressor-specific
efficiency curve.
"""
import csv
import os

from CoolProp.CoolProp import PropsSI

_HERE = os.path.dirname(os.path.abspath(__file__))
REPORT_DIR = os.path.join(_HERE, "report")
FIG_DIR = os.path.join(REPORT_DIR, "figures", "refrigerant_comparison")
TABLE_DIR = os.path.join(REPORT_DIR, "tables", "refrigerant_comparison")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)

# Design conditions -- identical to the rest of this project's design point
T_EVAP_C = 7.0
T_COND_C = 52.0
SUPERHEAT_K = 8.0
SUBCOOLING_K = 5.0
ETA_ISENTROPIC = 0.70   # single fixed assumption, same for all 3 candidates
Q_DESIGN_W = 150e3

# GWP (AR5, 100-yr ITH, per ASHRAE Standard 34 GWP tables) and ASHRAE 34
# safety classification for each candidate. All three satisfy the project's
# GWP<150 and A1/A2L/A3-only constraints by construction (that's WHY these
# three and not R-32, R-454B (GWP>150), CO2 (transcritical), or NH3 (B2L)).
CANDIDATES = [
    {"name": "R-290 (propane)",     "fluid": "R290",     "gwp": 3,  "safety": "A3"},
    {"name": "R-1234yf",            "fluid": "R1234yf",  "gwp": 1,  "safety": "A2L"},
    {"name": "R-1234ze(E)",         "fluid": "R1234ze(E)", "gwp": 1, "safety": "A2L"},
]


def evaluate_candidate(fluid):
    T_evap = T_EVAP_C + 273.15
    T_cond = T_COND_C + 273.15

    p_evap = PropsSI("P", "T", T_evap, "Q", 1, fluid)
    p_cond = PropsSI("P", "T", T_cond, "Q", 1, fluid)
    T_crit = PropsSI("Tcrit", fluid) - 273.15

    T_suction = T_evap + SUPERHEAT_K
    h_suction = PropsSI("HMASS", "T", T_suction, "P", p_evap, fluid)
    s_suction = PropsSI("SMASS", "T", T_suction, "P", p_evap, fluid)

    h_disch_isentropic = PropsSI("HMASS", "P", p_cond, "SMASS", s_suction, fluid)
    w_isentropic = h_disch_isentropic - h_suction
    w_actual = w_isentropic / ETA_ISENTROPIC
    h_discharge = h_suction + w_actual
    T_discharge_C = PropsSI("T", "P", p_cond, "HMASS", h_discharge, fluid) - 273.15

    T_liquid = T_cond - SUBCOOLING_K
    h_liquid = PropsSI("HMASS", "T", T_liquid, "P", p_cond, fluid)

    q_evap_specific = h_suction - h_liquid       # refrigerating effect [J/kg]
    COP = q_evap_specific / w_actual
    m_dot = Q_DESIGN_W / q_evap_specific          # kg/s at the 150 kW design duty
    pressure_ratio = p_cond / p_evap

    rho_suction = PropsSI("D", "T", T_suction, "P", p_evap, fluid)
    vol_refrigeration_effect = q_evap_specific * rho_suction / 1000.0  # kJ/m3, suction gas basis

    return {
        "p_evap_bar": p_evap / 1e5, "p_cond_bar": p_cond / 1e5,
        "pressure_ratio": pressure_ratio, "T_crit_C": T_crit,
        "T_discharge_C": T_discharge_C, "COP": COP,
        "m_dot_kgs": m_dot, "vol_refrig_effect_kJm3": vol_refrigeration_effect,
    }


def run_comparison():
    rows = []
    for cand in CANDIDATES:
        result = evaluate_candidate(cand["fluid"])
        rows.append({**cand, **result})
    return rows


def write_csv(rows, out_path):
    fieldnames = ["name", "fluid", "gwp", "safety", "p_evap_bar", "p_cond_bar",
                  "pressure_ratio", "T_crit_C", "T_discharge_C", "COP",
                  "m_dot_kgs", "vol_refrig_effect_kJm3"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {out_path}")


def write_latex_table(rows, out_path):
    lines = []
    lines.append(r"\begin{table}[h!]")
    lines.append(r"\caption{Working-fluid screening comparison: idealized single-stage "
                  r"cycle (7\,\textdegree C evap. / 52\,\textdegree C cond., 8\,K superheat, "
                  r"5\,K subcooling, $\eta_{\mathrm{isentropic}}=0.70$ for all three candidates).}")
    lines.append(r"\label{tab:refrigerant_comparison}")
    lines.append(r"\begin{tabularx}{\linewidth}{lLLLLLLL}")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Refrigerant} & \textbf{GWP} & \textbf{Safety} & "
                  r"\textbf{$p_{evap}$ [bar]} & \textbf{$p_{cond}$ [bar]} & "
                  r"\textbf{PR [-]} & \textbf{$T_{crit}$ [\textdegree C]} & \textbf{COP [-]} \\")
    lines.append(r"\midrule")
    for r in rows:
        lines.append(
            f"{r['name']} & {r['gwp']} & {r['safety']} & {r['p_evap_bar']:.2f} & "
            f"{r['p_cond_bar']:.2f} & {r['pressure_ratio']:.2f} & {r['T_crit_C']:.1f} & "
            f"{r['COP']:.3f} \\\\"
        )
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabularx}")
    lines.append(r"\end{table}")
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Saved {out_path}")


def plot_comparison(rows, out_path):
    import numpy as np
    import matplotlib.pyplot as plt

    BLUE = "#2a78d6"
    AQUA = "#1baf7a"
    YELLOW = "#eda100"
    GREEN = "#008300"
    SURFACE = "#fcfcfb"
    PRIMARY_INK = "#0b0b0b"
    SECONDARY_INK = "#52514e"
    MUTED = "#898781"
    GRID = "#e1e0d9"
    BASELINE = "#c3c2b7"
    colors = [BLUE, YELLOW, GREEN]  # fixed categorical order, one per refrigerant identity

    names = [r["name"] for r in rows]
    cop = [r["COP"] for r in rows]
    x = np.arange(len(rows))

    fig, ax = plt.subplots(figsize=(8, 5), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.grid(True, axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(BASELINE)
    ax.spines["bottom"].set_color(BASELINE)
    ax.tick_params(colors=MUTED, labelsize=9)

    bars = ax.bar(x, cop, color=colors, width=0.55, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=9)
    ax.set_ylabel("Idealized cycle COP [-]", color=SECONDARY_INK, fontsize=10)
    ax.set_title("Working-fluid screening: COP at identical cycle conditions\n"
                  "(7C evap / 52C cond, 8K superheat, 5K subcool, "
                  f"{ETA_ISENTROPIC:.0%} isentropic efficiency, all candidates)",
                  color=PRIMARY_INK, fontsize=11, loc="left")
    for b, v in zip(bars, cop):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.2f}", ha="center",
                 va="bottom", fontsize=9, color=SECONDARY_INK)

    ax.text(0.0, -0.16, "All candidates satisfy GWP<150 and ASHRAE 34 safety class A1/A2L/A3 "
                         "by construction (screening constraint, project spec).",
            transform=ax.transAxes, fontsize=8, color=MUTED)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    rows = run_comparison()

    print(f"{'Refrigerant':<20} {'GWP':>5} {'Safety':>7} {'p_evap':>8} {'p_cond':>8} "
          f"{'PR':>6} {'Tcrit':>7} {'Tdis':>7} {'COP':>7}")
    for r in rows:
        print(f"{r['name']:<20} {r['gwp']:>5} {r['safety']:>7} {r['p_evap_bar']:>8.2f} "
              f"{r['p_cond_bar']:>8.2f} {r['pressure_ratio']:>6.2f} {r['T_crit_C']:>7.1f} "
              f"{r['T_discharge_C']:>7.1f} {r['COP']:>7.3f}")

    write_csv(rows, os.path.join(TABLE_DIR, "refrigerant_comparison.csv"))
    write_latex_table(rows, os.path.join(TABLE_DIR, "refrigerant_comparison.tex"))
    plot_comparison(rows, os.path.join(FIG_DIR, "refrigerant_comparison_cop.png"))
