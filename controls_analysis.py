"""Quantitative support for the compressor-EEV control discussion (subtask 3b).

The brief asks for "consideration of" the compressor/expansion-device
interaction, which is a written engineering discussion -- control interaction
is inherently transient, and this project's model is steady-state throughout.
A dynamic model (mass/energy storage, valve dynamics, tuned PID loops) is a
different exercise and is not attempted.

What the steady-state model CAN establish, and does here:
  1. EEV CONTROL AUTHORITY -- the valve opening required at every point of the
     5:1 load turndown. This determines whether the selected valve can actually
     modulate across the operating envelope, or saturates at the bottom of it.
  2. SUPERHEAT SETPOINT SENSITIVITY -- how much COP the superheat setpoint is
     worth, which sets how hard the EEV loop should be pushed.

Both use the same Bitzer AHRI-540 map and CoolProp properties as the rest of
the model.
"""
import os

import numpy as np
from CoolProp.CoolProp import PropsSI

from compressor import Compressor, CompressorBank
from eev import EEV

_HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(_HERE, "report")
FIG_DIR = os.path.join(REPORT, "figures", "controls")
TABLE_DIR = os.path.join(REPORT, "tables", "controls")
for d in (FIG_DIR, TABLE_DIR):
    os.makedirs(d, exist_ok=True)

R = "R290"
SH_DESIGN, SC, TO_SET = 8.0, 5.0, 10.0
APPROACH_COND_K = 10.0

# Load points: the four IPLV points plus the project's 30 kW minimum IT load.
# Each sits at its own reduced condenser-entering-air temperature, as AHRI
# 550/590 defines -- so this envelope is the real one the valve must cover.
# (percent, note, load fraction, condenser-entering-air temp). The percent is kept
# NUMERIC and the "%" sign is added per output format -- a bare "%" is a comment
# character in LaTeX and would silently swallow the rest of a table row.
ENVELOPE = [
    (100, "design", 1.00, 35.0),
    (75, "", 0.75, 26.7),
    (50, "", 0.50, 18.3),
    (25, "", 0.25, 12.8),
    (20, "30 kW min IT load", 0.20, 12.8),
]

# EEV control-authority floor: below roughly this opening a thermostatic/
# electronic valve's characteristic goes non-linear and stepper resolution and
# stiction dominate, so superheat control degrades. Not a hard cliff -- a
# design guideline.
AUTHORITY_FLOOR_PCT = 15.0


def eev_authority(n_circuits):
    """EEV opening (% of the CAREL E3V65's rated Kv) at each envelope point,
    for a given number of parallel refrigerant circuits."""
    bank = CompressorBank(n_units=2, superheat_K=SH_DESIGN, subcooling_K=SC)
    rows = []
    for pct, note, frac, t_air in ENVELOPE:
        tc = t_air + APPROACH_COND_K
        r = bank.solve_with_vfd(frac * 150e3, tc, to_setpoint_C=TO_SET)
        p_c = PropsSI("P", "T", tc + 273.15, "Q", 1, R)
        p_o = PropsSI("P", "T", r["to_C"] + 273.15, "Q", 1, R)
        h_3 = PropsSI("HMASS", "T", tc - SC + 273.15, "P", p_c, R)
        m_per = r["m_dot_total_kgh"] / 3600.0 / n_circuits
        e = EEV(m_dot_refrigerant=m_per, p_in=p_c, h_in=h_3, p_out=p_o, refrigerant=R)
        rows.append({
            "pct": pct, "note": note, "load_frac": frac, "t_air_C": t_air, "tc_C": tc,
            "to_C": r["to_C"], "m_dot_kgs": m_per, "dP_bar": (p_c - p_o) / 1e5,
            "Kv": e.Kv_required, "open_pct": e.Kv_required / e.Kv_rated * 100.0,
        })
    return rows


def superheat_sensitivity():
    """Design-point COP vs. the EEV's superheat setpoint, AT FIXED to.

    READ THIS CAREFULLY -- the trend is an artifact, not a design lever. The
    AHRI-540/EN12900 map credits ALL superheat as USEFUL (i.e. picked up in the
    evaporator, so it adds to the refrigerating effect h1-h3). Holding to fixed
    and raising the setpoint therefore raises h1 and hence COP, suggesting that
    more superheat is better.

    Physically it is the opposite. The evaporator's area is finite: boiling
    occupies ~95% of the duty and the superheat tail ~5%, but the superheat zone
    has a far poorer gas-side heat-transfer coefficient and so consumes area
    disproportionately. Demanding more superheat at fixed area steals area from
    boiling, which forces to DOWN -- and a lower to costs more COP than the
    higher h1 gains. Capturing that properly needs an area-constrained
    evaporator solve, which this function does NOT do.

    It is retained because the artifact itself is the point: it shows why the
    superheat setpoint must be argued from compressor protection and control
    stability, not read off a compressor map.
    """
    comp = Compressor()
    rows = []
    for sh in (4.0, 6.0, 8.0, 10.0, 12.0, 15.0):
        p = comp.performance(TO_SET, 45.0, sh, SC)
        rows.append({"sh_K": sh, "COP": p["COP"], "Q_kW": 2 * p["Q_w"] / 1e3})
    return rows


# ============================ FIGURES ================================
BLUE, AQUA, YELLOW, GREEN = "#2a78d6", "#1baf7a", "#eda100", "#008300"
VIOLET, RED, ORANGE = "#4a3aa7", "#e34948", "#eb6834"
SURFACE, PRIMARY, SECONDARY = "#fcfcfb", "#0b0b0b", "#52514e"
MUTED, GRID, BASELINE = "#898781", "#e1e0d9", "#c3c2b7"


def _style(ax, axis="both"):
    ax.set_facecolor(SURFACE)
    ax.grid(True, axis=axis, color=GRID, linewidth=0.8, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(BASELINE)
    ax.tick_params(colors=MUTED, labelsize=9)


def plot_authority(one, two, out_path):
    import matplotlib.pyplot as plt

    labels = [f"{r['pct']}%" + (f"\n({r['note']})" if r["note"] else "") for r in one]
    x = np.arange(len(labels))
    w = 0.36
    fig, ax = plt.subplots(figsize=(10, 5.4), facecolor=SURFACE)
    _style(ax, axis="y")
    b1 = ax.bar(x - w/2, [r["open_pct"] for r in one], w, color=BLUE, zorder=3,
                label="One EEV, full flow (as main.py models)")
    b2 = ax.bar(x + w/2, [r["open_pct"] for r in two], w, color=YELLOW, zorder=3,
                label="Two EEVs, split flow (as Ch. Compressor / BOM specify)")
    ax.axhline(AUTHORITY_FLOOR_PCT, color=RED, linestyle="--", linewidth=1.6, zorder=4)
    ax.text(len(labels) - 0.45, AUTHORITY_FLOOR_PCT + 1.5,
            f"$\\approx${AUTHORITY_FLOOR_PCT:.0f}\\% --- below this, EEV control degrades",
            ha="right", fontsize=8.5, color=RED)
    for bars in (b1, b2):
        for b in bars:
            v = b.get_height()
            ax.text(b.get_x() + b.get_width()/2, v + 1.2, f"{v:.0f}%", ha="center",
                    fontsize=8.5, color=SECONDARY)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("EEV opening [% of E3V65 rated $K_v$]", color=SECONDARY, fontsize=10)
    ax.set_title("EEV control authority across the 5:1 load turndown",
                 color=PRIMARY, fontsize=12, loc="left")
    ax.set_ylim(0, 100)
    ax.legend(fontsize=8.5, edgecolor=BASELINE, loc="upper right")
    ax.text(0.0, -0.20, "Each load point sits at its own reduced condenser-entering-air temperature "
                        "(AHRI 550/590), so this is the real envelope.\nThe two-EEV case runs out of "
                        "authority at minimum load; the one-EEV case does not. The architecture is "
                        "unresolved (see text).",
            transform=ax.transAxes, fontsize=8, color=MUTED)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_superheat(rows, out_path):
    import matplotlib.pyplot as plt

    sh = [r["sh_K"] for r in rows]
    cop = [r["COP"] for r in rows]
    fig, ax = plt.subplots(figsize=(9, 4.8), facecolor=SURFACE)
    _style(ax)
    ax.plot(sh, cop, color=BLUE, linewidth=2.2, marker="o", markersize=7, zorder=4)
    ax.axvline(SH_DESIGN, color=GREEN, linestyle=":", linewidth=1.8, zorder=3)
    ax.text(SH_DESIGN + 0.15, min(cop), "design setpoint (8 K)", rotation=90,
            va="bottom", fontsize=8.5, color=GREEN)
    for s, c in zip(sh, cop):
        ax.text(s, c + 0.006, f"{c:.3f}", ha="center", fontsize=8, color=SECONDARY)
    ax.set_xlabel("EEV superheat setpoint [K]", color=SECONDARY, fontsize=10)
    ax.set_ylabel("Design-point COP [-]  (compressor map, at FIXED $t_o$)",
                  color=SECONDARY, fontsize=10)
    ax.set_title("A trap: what the compressor map says about superheat --- and why it is wrong",
                 color=PRIMARY, fontsize=11.5, loc="left")
    span = (max(cop) - min(cop)) / np.mean(cop) * 100
    ax.annotate("map says: more superheat = better COP\n(ARTIFACT --- see caption)",
                xy=(sh[-1], cop[-1]), xytext=(-12, -38), textcoords="offset points",
                ha="right", fontsize=9, color=RED, fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=RED, lw=1.4))
    ax.text(0.0, -0.24,
            f"The map credits ALL superheat as USEFUL (EN12900 convention), so at FIXED $t_o$ a higher setpoint "
            f"raises $h_1$ and COP ({span:.1f}% across 4-15 K).\nBut evaporator area is finite: the superheat tail "
            f"has a poor gas-side coefficient and consumes area disproportionately, so demanding more\nsuperheat "
            f"steals area from boiling and forces $t_o$ DOWN --- costing more COP than the higher $h_1$ gains. "
            f"The real trend is the opposite.",
            transform=ax.transAxes, fontsize=8, color=MUTED)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Saved {out_path}")


def write_tables(one, two, sh_rows, out_dir):
    lines = [r"\begin{table}[h!]",
             r"\caption{EEV control authority across the operating envelope. Each load point "
             r"sits at its own reduced condenser-entering-air temperature per AHRI 550/590.}",
             r"\label{tab:eev_authority}",
             r"\begin{tabularx}{\linewidth}{lLLLLL}", r"\toprule",
             r"\textbf{Load} & \textbf{$T_{air}$ [\textdegree C]} & \textbf{$t_c$ [\textdegree C]} "
             r"& \textbf{$\Delta p$ [bar]} & \textbf{Opening, 1 EEV} & \textbf{Opening, 2 EEVs} \\",
             r"\midrule"]
    for a, b in zip(one, two):
        flag = r" $^{\dagger}$" if b["open_pct"] < AUTHORITY_FLOOR_PCT else ""
        lbl = f"{a['pct']}\,\%" + (f" ({a['note']})" if a['note'] else "")
        lines.append(f"{lbl} & {a['t_air_C']:.1f} & {a['tc_C']:.1f} & {a['dP_bar']:.2f} & "
                     f"{a['open_pct']:.1f}\\,\\% & {b['open_pct']:.1f}\\,\\%{flag} \\\\")
    lines += [r"\bottomrule", r"\end{tabularx}",
              r"{\small $^{\dagger}$ below the $\approx$15\,\% authority floor, where the valve's "
              r"characteristic goes non-linear and superheat control degrades.}",
              r"\end{table}"]
    p = os.path.join(out_dir, "eev_authority.tex")
    open(p, "w").write("\n".join(lines) + "\n")
    print(f"Saved {p}")

    lines = [r"\begin{table}[h!]",
             r"\caption{Design-point COP vs.\ EEV superheat setpoint, evaluated at FIXED $t_o$. "
             r"The rising trend is an artifact of the map's useful-superheat convention, not a design "
             r"lever --- see the discussion.}",
             r"\label{tab:superheat_sens}",
             r"\begin{tabularx}{\linewidth}{lLL}", r"\toprule",
             r"\textbf{Superheat setpoint} & \textbf{Design COP} & \textbf{Bank capacity at full speed} \\",
             r"\midrule"]
    for r in sh_rows:
        mark = r" \textbf{(design)}" if abs(r["sh_K"] - SH_DESIGN) < 1e-6 else ""
        lines.append(f"{r['sh_K']:.0f}\\,K{mark} & {r['COP']:.3f} & {r['Q_kW']:.1f}\\,kW \\\\")
    lines += [r"\bottomrule", r"\end{tabularx}", r"\end{table}"]
    p = os.path.join(out_dir, "superheat_sensitivity.tex")
    open(p, "w").write("\n".join(lines) + "\n")
    print(f"Saved {p}")


if __name__ == "__main__":
    one, two = eev_authority(1), eev_authority(2)
    print("EEV AUTHORITY (% of E3V65 rated Kv)")
    print(f"{'Load':<26}{'1 EEV':>9}{'2 EEVs':>9}")
    for a, b in zip(one, two):
        flag = "  <-- below authority floor" if b["open_pct"] < AUTHORITY_FLOOR_PCT else ""
        nm = f"{a['pct']}%" + (f" ({a['note']})" if a['note'] else '')
        print(f"{nm:<26}{a['open_pct']:>8.1f}%{b['open_pct']:>8.1f}%{flag}")
    print()
    sh_rows = superheat_sensitivity()
    print("SUPERHEAT SENSITIVITY (design point)")
    for r in sh_rows:
        print(f"  {r['sh_K']:>4.0f} K -> COP {r['COP']:.3f}   bank capacity {r['Q_kW']:.1f} kW")
    cops = [r["COP"] for r in sh_rows]
    print(f"  span across 4-15 K: {(max(cops)-min(cops))/np.mean(cops)*100:.1f}% of COP")

    plot_authority(one, two, os.path.join(FIG_DIR, "eev_authority.png"))
    plot_superheat(sh_rows, os.path.join(FIG_DIR, "superheat_sensitivity.png"))
    write_tables(one, two, sh_rows, TABLE_DIR)
