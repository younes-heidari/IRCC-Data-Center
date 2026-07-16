"""R-290 refrigerant charge inventory and EN 378 safety assessment.

The project brief admits an A3 refrigerant only "with charge analysis"
(subtask 1a), so this is a condition on the R-290 selection itself, not an
optional extra. It also underpins subtask 1i (safety measures).

Two parts:
  1. CHARGE INVENTORY -- component-by-component R-290 mass at the revised
     design point, built from the same geometry and CoolProp properties the
     rest of the model uses.
  2. EN 378 ASSESSMENT -- the resulting charge against R-290's flammability
     limits, the room volume an unventilated installation would need, and the
     machinery-room provisions (ventilation, detection) that make the real
     installation compliant.

Charge is reported as a RANGE, not a point value: refrigerant distribution
shifts with operating mode, and the void fractions in the two-phase
exchangers are genuinely uncertain. Sizing safety systems on the upper bound
is the defensible choice.
"""
import os

import numpy as np
from CoolProp.CoolProp import PropsSI

_HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(_HERE, "report")
FIG_SAFETY = os.path.join(REPORT, "figures", "safety")
TABLE_SAFETY = os.path.join(REPORT, "tables", "safety")
for d in (FIG_SAFETY, TABLE_SAFETY):
    os.makedirs(d, exist_ok=True)

R = "R290"

# ---------------------------------------------------------------------
# R-290 flammability data (EN 378-1 Annex E / ASHRAE 34)
# ---------------------------------------------------------------------
LFL = 0.038          # kg/m3, lower flammability limit
PRACTICAL_LIMIT = 0.008   # kg/m3, EN 378-1 practical limit for R-290 (~21% of LFL)
DETECTION_FRACTION = 0.20  # EN 378 requires detection at <= 20% LFL

# ---------------------------------------------------------------------
# Revised-design geometry (matches main.py) and state points
# ---------------------------------------------------------------------
TO, TC, SH, SC = 10.0, 45.0, 8.0, 5.0
D_H, PHI = 0.004, 1.17
B = D_H * PHI / 2               # 2.34 mm channel spacing
L_W, L_V = 0.195, 0.600         # confirmed Kelvion HP DW 500H plate
N_CP_EVAP, N_CP_COND = 38, 47   # refrigerant channels (main.py)

# Line sizes at the revised design point (main.py's select_tube output)
ID_SUCTION, ID_DISCHARGE, ID_LIQUID = 0.06287, 0.03823, 0.03388
LEN_SUCTION, LEN_DISCHARGE, LEN_LIQUID = 15.0, 15.0, 15.0   # ASSUMED, see report

RECEIVER_VOLUME = 0.025          # m3, Bitzer F252HP (25 L)
RECEIVER_FILL_RANGE = (0.20, 0.50)   # operating liquid fill fraction, assumed range


def _rho(T_C, Q=None, P=None):
    if Q is not None:
        return PropsSI("D", "T", T_C + 273.15, "Q", Q, R)
    return PropsSI("D", "T", T_C + 273.15, "P", P, R)


def inventory():
    """Component-by-component charge, as (low, high) bounds in kg."""
    p_o = PropsSI("P", "T", TO + 273.15, "Q", 1, R)
    p_c = PropsSI("P", "T", TC + 273.15, "Q", 1, R)

    V_evap = N_CP_EVAP * B * L_W * L_V
    V_cond = N_CP_COND * B * L_W * L_V
    V_suc = np.pi / 4 * ID_SUCTION**2 * LEN_SUCTION
    V_dis = np.pi / 4 * ID_DISCHARGE**2 * LEN_DISCHARGE
    V_liq = np.pi / 4 * ID_LIQUID**2 * LEN_LIQUID

    rho_l_c, rho_v_c = _rho(TC, Q=0), _rho(TC, Q=1)
    rho_l_o, rho_v_o = _rho(TO, Q=0), _rho(TO, Q=1)
    rho_liq_line = _rho(TC - SC, P=p_c)
    rho_suc = PropsSI("D", "T", TO + SH + 273.15, "P", p_o, R)
    rho_dis = _rho(TC + 20, P=p_c)   # ~superheated discharge

    # Two-phase exchangers: void fraction is the real uncertainty. Bound the
    # condenser between mostly-vapour and mostly-liquid holdup; the evaporator
    # is vapour-dominated by volume (low-density suction-side gas).
    rows = [
        ("Condenser (2-phase)", V_cond, 0.15 * rho_l_c + 0.85 * rho_v_c,
                                 0.45 * rho_l_c + 0.55 * rho_v_c),
        ("Evaporator (2-phase)", V_evap, 0.05 * rho_l_o + 0.95 * rho_v_o,
                                  0.20 * rho_l_o + 0.80 * rho_v_o),
        ("Liquid line (15 m)", V_liq, rho_liq_line, rho_liq_line),
        ("Suction line (15 m)", V_suc, rho_suc, rho_suc),
        ("Discharge line (15 m)", V_dis, rho_dis, rho_dis),
        ("Liquid receiver (25 L)", RECEIVER_VOLUME,
            RECEIVER_FILL_RANGE[0] * rho_l_c + (1 - RECEIVER_FILL_RANGE[0]) * rho_v_c,
            RECEIVER_FILL_RANGE[1] * rho_l_c + (1 - RECEIVER_FILL_RANGE[1]) * rho_v_c),
        ("Oil separator (40 L vessel)", 0.040, rho_v_c, rho_v_c),
        ("Compressors (2 x, oil-dissolved)", None, 0.4, 1.2),  # kg directly, not V*rho
    ]
    out = []
    for name, V, lo, hi in rows:
        if V is None:
            out.append((name, lo, hi, None))
        else:
            out.append((name, V * lo, V * hi, V))
    return out


def totals(inv):
    lo = sum(r[1] for r in inv)
    hi = sum(r[2] for r in inv)
    return lo, hi


def en378(m_charge):
    """EN 378 assessment for a given charge."""
    return {
        "charge": m_charge,
        # room volume needed so a FULL release stays under each threshold
        "V_practical": m_charge / PRACTICAL_LIMIT,
        "V_20pct_LFL": m_charge / (DETECTION_FRACTION * LFL),
        "V_LFL": m_charge / LFL,
        # EN 378-3 emergency ventilation for a machinery room: V = 0.014 * m^(2/3) [m3/s]
        "vent_m3s": 0.014 * m_charge ** (2 / 3),
        "vent_m3h": 0.014 * m_charge ** (2 / 3) * 3600,
        "detect_setpoint": DETECTION_FRACTION * LFL,
    }


# ============================ FIGURES ================================
BLUE, AQUA, YELLOW, GREEN = "#2a78d6", "#1baf7a", "#eda100", "#008300"
VIOLET, RED, MAGENTA, ORANGE = "#4a3aa7", "#e34948", "#e87ba4", "#eb6834"
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


def plot_inventory(inv, out_path):
    import matplotlib.pyplot as plt

    inv_s = sorted(inv, key=lambda r: -(r[1] + r[2]) / 2)
    names = [r[0] for r in inv_s]
    lo = np.array([r[1] for r in inv_s])
    hi = np.array([r[2] for r in inv_s])
    mid = (lo + hi) / 2
    y = np.arange(len(inv_s))

    fig, ax = plt.subplots(figsize=(10, 5.4), facecolor=SURFACE)
    _style(ax, axis="x")
    ax.barh(y, mid, color=BLUE, height=0.6, zorder=3)
    ax.errorbar(mid, y, xerr=[mid - lo, hi - mid], fmt="none", ecolor=SECONDARY,
                capsize=4, elinewidth=1.4, zorder=4)
    for i, (m, l, h) in enumerate(zip(mid, lo, hi)):
        ax.text(h + 0.12, i, f"{l:.2f}-{h:.2f} kg", va="center", fontsize=8.5, color=SECONDARY)
    ax.set_yticks(y); ax.set_yticklabels(names, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("R-290 charge [kg]   (bar = mid-estimate, whiskers = bounds)",
                  color=SECONDARY, fontsize=10)
    tl, th = totals(inv)
    ax.set_title(f"R-290 charge inventory at the revised design point --- total {tl:.1f} to {th:.1f} kg",
                 color=PRIMARY, fontsize=12, loc="left")
    ax.set_xlim(0, max(hi) * 1.42)
    ax.text(0.0, -0.155, "Ranges reflect two-phase void fraction (exchangers), receiver fill level, and "
                          "oil-dissolved refrigerant.\nLine lengths are an assumed 15 m each --- the project brief "
                          "scopes piping to diameters only.",
            transform=ax.transAxes, fontsize=8, color=MUTED)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_concentration(m_charge, out_path):
    """Room concentration if the entire charge is released, vs room volume,
    against the EN 378 thresholds. Shows why a machinery room (ventilation +
    detection) is required rather than optional."""
    import matplotlib.pyplot as plt

    V = np.linspace(20, 1500, 600)
    c = m_charge / V

    fig, ax = plt.subplots(figsize=(10, 5.6), facecolor=SURFACE)
    _style(ax)
    ax.plot(V, c * 1000, color=BLUE, linewidth=2.5, zorder=5,
            label=f"Concentration after full {m_charge:.1f} kg release")
    # practical limit (8) and 20% LFL (7.6) nearly coincide -- separate their labels
    for lim, col, lbl, xl, dy in (
            (LFL, RED, f"LFL = {LFL*1000:.0f} g/m$^3$ (flammable above)", 1490, 1.10),
            (PRACTICAL_LIMIT, GREEN,
             f"Practical limit {PRACTICAL_LIMIT*1000:.0f} g/m$^3$ (EN 378 occupied space)", 1490, 1.14),
            (DETECTION_FRACTION*LFL, YELLOW,
             f"20% LFL = {DETECTION_FRACTION*LFL*1000:.1f} g/m$^3$ (detection setpoint)", 400, 0.72)):
        ax.axhline(lim*1000, color=col, linestyle="--", linewidth=1.6, zorder=4)
        ax.text(xl, lim*1000*dy, lbl, ha="right", fontsize=8.5, color=SECONDARY)
        Vx = m_charge/lim
        if V[0] <= Vx <= V[-1]:
            ax.plot([Vx], [lim*1000], "o", color=col, markersize=8, zorder=6)
            ax.annotate(f"{Vx:.0f} m$^3$", (Vx, lim*1000), textcoords="offset points",
                        xytext=(10, 8), fontsize=8.5, color=col, fontweight="bold")

    V_room = 90.0
    ax.axvline(V_room, color=VIOLET, linestyle=":", linewidth=1.8, zorder=4)
    ax.plot([V_room], [m_charge/V_room*1000], "D", color=VIOLET, markersize=9, zorder=7)
    ax.annotate(f"Assumed machinery room\n{V_room:.0f} m$^3$ -> {m_charge/V_room*1000:.0f} g/m$^3$\n"
                f"= {m_charge/V_room/LFL:.1f}x LFL",
                (V_room, m_charge/V_room*1000), textcoords="offset points", xytext=(16, -6),
                fontsize=9, color=VIOLET, fontweight="bold")

    ax.set_yscale("log")
    ax.set_xlim(20, 1500); ax.set_ylim(3, 700)
    ax.set_xlabel("Room volume [m$^3$]", color=SECONDARY, fontsize=10)
    ax.set_ylabel("R-290 concentration after full release [g/m$^3$]  (log)",
                  color=SECONDARY, fontsize=10)
    ax.set_title("Why a ventilated machinery room is mandatory, not optional",
                 color=PRIMARY, fontsize=12, loc="left")
    ax.legend(loc="upper right", fontsize=9, edgecolor=BASELINE)
    ax.text(0.0, -0.145, "A realistically sized plant room sits far above the LFL on a full release. "
                          "Compliance therefore rests on EN 378\nmachinery-room provisions --- leak detection at "
                          "20% LFL plus emergency ventilation --- not on dilution by room volume.",
            transform=ax.transAxes, fontsize=8, color=MUTED)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Saved {out_path}")


def write_tables(inv, m_hi, out_dir):
    tl, th = totals(inv)
    lines = [r"\begin{table}[h!]",
             r"\caption{R-290 charge inventory at the revised design point. Ranges reflect "
             r"two-phase void fraction, receiver fill level, and oil-dissolved refrigerant.}",
             r"\label{tab:charge_inventory}",
             r"\begin{tabularx}{\linewidth}{lLLL}", r"\toprule",
             r"\textbf{Component} & \textbf{Internal volume} & \textbf{Charge (low)} & \textbf{Charge (high)} \\",
             r"\midrule"]
    for name, lo, hi, V in sorted(inv, key=lambda r: -(r[1]+r[2])/2):
        vtxt = f"{V*1000:.1f}\\,L" if V else "---"
        lines.append(f"{name} & {vtxt} & {lo:.2f}\\,kg & {hi:.2f}\\,kg \\\\")
    lines += [r"\midrule",
              rf"\textbf{{Total system charge}} & --- & \textbf{{{tl:.1f}\,kg}} & \textbf{{{th:.1f}\,kg}} \\",
              r"\bottomrule", r"\end{tabularx}", r"\end{table}"]
    p = os.path.join(out_dir, "charge_inventory.tex")
    open(p, "w").write("\n".join(lines) + "\n")
    print(f"Saved {p}")

    e = en378(m_hi)
    lines = [r"\begin{table}[h!]",
             rf"\caption{{EN 378 assessment at the upper-bound charge ({m_hi:.1f}\,kg).}}",
             r"\label{tab:en378}",
             r"\begin{tabularx}{\linewidth}{lLL}", r"\toprule",
             r"\textbf{Quantity} & \textbf{Value} & \textbf{Basis} \\", r"\midrule",
             rf"R-290 LFL & {LFL*1000:.0f}\,g/m$^3$ & EN 378-1 Annex E \\",
             rf"Practical limit & {PRACTICAL_LIMIT*1000:.0f}\,g/m$^3$ & EN 378-1, occupied space \\",
             rf"Detection setpoint (20\,\% LFL) & {e['detect_setpoint']*1000:.1f}\,g/m$^3$ & EN 378-3 \\",
             r"\midrule",
             rf"Room volume for full release $<$ practical limit & {e['V_practical']:.0f}\,m$^3$ & $m/0.008$ --- impractical \\",
             rf"Room volume for full release $<$ 20\,\% LFL & {e['V_20pct_LFL']:.0f}\,m$^3$ & $m/(0.2\cdot\mathrm{{LFL}})$ \\",
             rf"Room volume for full release $<$ LFL & {e['V_LFL']:.0f}\,m$^3$ & $m/\mathrm{{LFL}}$ --- no margin \\",
             r"\midrule",
             rf"\textbf{{Emergency ventilation airflow}} & \textbf{{{e['vent_m3h']:.0f}\,m$^3$/h}} ({e['vent_m3s']:.3f}\,m$^3$/s) & EN 378-3: $\dot V = 0.014\,m^{{2/3}}$ \\",
             r"\bottomrule", r"\end{tabularx}", r"\end{table}"]
    p = os.path.join(out_dir, "en378_assessment.tex")
    open(p, "w").write("\n".join(lines) + "\n")
    print(f"Saved {p}")


if __name__ == "__main__":
    inv = inventory()
    tl, th = totals(inv)
    print(f"{'Component':<36}{'low':>9}{'high':>9}")
    print("-" * 56)
    for name, lo, hi, V in sorted(inv, key=lambda r: -(r[1]+r[2])/2):
        print(f"{name:<36}{lo:>8.2f}kg{hi:>8.2f}kg")
    print("-" * 56)
    print(f"{'TOTAL':<36}{tl:>8.1f}kg{th:>8.1f}kg")
    print()
    e = en378(th)
    print(f"EN 378 at the upper bound ({th:.1f} kg):")
    print(f"  room for <practical limit : {e['V_practical']:.0f} m3   (impractical)")
    print(f"  room for <20% LFL         : {e['V_20pct_LFL']:.0f} m3")
    print(f"  room for <LFL             : {e['V_LFL']:.0f} m3")
    print(f"  emergency ventilation     : {e['vent_m3h']:.0f} m3/h")
    print(f"  detection setpoint        : {e['detect_setpoint']*1000:.1f} g/m3")

    plot_inventory(inv, os.path.join(FIG_SAFETY, "charge_inventory.png"))
    plot_concentration(th, os.path.join(FIG_SAFETY, "charge_concentration.png"))
    write_tables(inv, th, TABLE_SAFETY)
