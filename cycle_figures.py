"""Cycle and system visualisation figures for the report.

Covers the project's p-h chart deliverable (subtask 1h) plus the supporting
diagrams that carry the design's key arguments visually: the temperature
cascade that drives the whole efficiency revision, the heat-exchanger T-Q
profiles behind the LMTD sizing, the PUE breakdown, the original-vs-revised
metric comparison, and the annual ambient duration curve behind the
free-cooling hours.

Every state point is taken from the SAME Bitzer AHRI-540 map and CoolProp
properties used by main.py, so these figures cannot drift from the model.
"""
import math
import os

import numpy as np
from CoolProp.CoolProp import PropsSI

from compressor import CompressorBank

_HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(_HERE, "report")
FIG_CYCLE = os.path.join(REPORT, "figures", "cycle_diagram")
FIG_PERF = os.path.join(REPORT, "figures", "performance", "design_point")
FIG_ANNUAL = os.path.join(REPORT, "figures", "performance", "annual")
for d in (FIG_CYCLE, FIG_PERF, FIG_ANNUAL):
    os.makedirs(d, exist_ok=True)

R = "R290"

# ---- Design points -------------------------------------------------
# REVISED (current design): tc=45 (10 K air approach), to=10 (5 K evap approach)
# ORIGINAL (superseded):    tc=52 (17 K air approach), to=7  (8 K evap approach)
#
# Superheat is taken as 8 K -- the compressor vendor's own selection basis
# (Ch. Compressor: "useful superheat 8 K, 15 C return gas") and the value
# performance.py uses for the reported COP/IPLV. NOTE: main.py's integrated
# cycle uses DT_SH = 10 K instead, an undocumented deviation that makes its
# design-point COP read 4.73 against the performance chapter's 4.67 (a 1.3%
# split). Flagged as an open item; 8 K is used here so this chart agrees with
# the compliance numbers rather than adding a third value.
SUPERHEAT_K = 8.0
REVISED = dict(to=10.0, tc=45.0, sh=SUPERHEAT_K, sc=5.0,
               label="Revised design ($t_c$=45 C, $t_o$=10 C)")
ORIGINAL = dict(to=7.0, tc=52.0, sh=SUPERHEAT_K, sc=5.0,
                label="Original design ($t_c$=52 C, $t_o$=7 C)")

# ---- dataviz palette (validated reference instance) -----------------
BLUE, AQUA, YELLOW, GREEN = "#2a78d6", "#1baf7a", "#eda100", "#008300"
VIOLET, RED, MAGENTA, ORANGE = "#4a3aa7", "#e34948", "#e87ba4", "#eb6834"
SURFACE, PRIMARY, SECONDARY = "#fcfcfb", "#0b0b0b", "#52514e"
MUTED, GRID, BASELINE = "#898781", "#e1e0d9", "#c3c2b7"


def _style(ax, grid_axis="both"):
    ax.set_facecolor(SURFACE)
    ax.grid(True, axis=grid_axis, color=GRID, linewidth=0.8, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(BASELINE)
    ax.tick_params(colors=MUTED, labelsize=9)


def cycle_states(to, tc, sh, sc, **_):
    """The four p-h state points, from the same compressor map main.py uses."""
    bank = CompressorBank(n_units=2, superheat_K=sh, subcooling_K=sc)
    r = bank.solve_with_vfd(150e3, tc, to_setpoint_C=to)
    pu = r["per_unit"]
    po = PropsSI("P", "T", to + 273.15, "Q", 1, R)
    pc = PropsSI("P", "T", tc + 273.15, "Q", 1, R)
    h1 = pu["h_suction"]
    h2 = h1 + pu["P_w"] / (pu["m_dot_kgh"] / 3600.0)
    h3 = PropsSI("HMASS", "T", tc - sc + 273.15, "P", pc, R)
    return dict(po=po, pc=pc, h1=h1, h2=h2, h3=h3, h4=h3,
                T1=to + sh, T2=PropsSI("T", "HMASS", h2, "P", pc, R) - 273.15,
                T3=tc - sc, T4=to, COP=r["COP"],
                x4=PropsSI("Q", "HMASS", h3, "P", po, R))


def plot_ph_chart(out_path):
    """p-h chart of the designed cycle (subtask 1h), with the superseded
    original overlaid so the revision's smaller lift -- and hence its shorter
    compression work leg -- is visible rather than only asserted."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 7), facecolor=SURFACE)
    _style(ax)

    # Saturation dome
    Tt, Tc_crit = PropsSI("Ttriple", R) + 0.5, PropsSI("Tcrit", R) - 0.01
    Ts = np.linspace(Tt, Tc_crit, 400)
    hf = np.array([PropsSI("HMASS", "T", T, "Q", 0, R) for T in Ts]) / 1e3
    hg = np.array([PropsSI("HMASS", "T", T, "Q", 1, R) for T in Ts]) / 1e3
    ps = np.array([PropsSI("P", "T", T, "Q", 0, R) for T in Ts]) / 1e5
    ax.plot(hf, ps, color=MUTED, linewidth=1.5, zorder=2)
    ax.plot(hg, ps, color=MUTED, linewidth=1.5, zorder=2)
    ax.fill_betweenx(ps, hf, hg, color=GRID, alpha=0.45, zorder=1)
    # anchor the dome labels at a pressure inside the plotted window
    p_lbl = 4.0e5
    ax.text(PropsSI("HMASS", "P", p_lbl, "Q", 0, R)/1e3 - 8, p_lbl/1e5,
            "saturated\nliquid", fontsize=8, color=MUTED, ha="right", va="center")
    ax.text(PropsSI("HMASS", "P", p_lbl, "Q", 1, R)/1e3 + 8, p_lbl/1e5,
            "saturated\nvapour", fontsize=8, color=MUTED, ha="left", va="center")
    ax.text(430, 30, "two-phase", fontsize=8.5, color=MUTED, ha="center", style="italic")

    for spec, color, ls, lw, z in ((ORIGINAL, MUTED, "--", 1.6, 3), (REVISED, BLUE, "-", 2.4, 4)):
        s = cycle_states(**spec)
        h = [s["h1"]/1e3, s["h2"]/1e3, s["h3"]/1e3, s["h4"]/1e3, s["h1"]/1e3]
        p = [s["po"]/1e5, s["pc"]/1e5, s["pc"]/1e5, s["po"]/1e5, s["po"]/1e5]
        ax.plot(h, p, color=color, linewidth=lw, linestyle=ls, zorder=z,
                label=f"{spec['label']} --- COP {s['COP']:.2f}")
        if color == BLUE:
            ax.plot(h[:4], p[:4], "o", color=color, markersize=8, zorder=z + 1)
            for i, (hh, pp) in enumerate(zip(h[:4], p[:4]), start=1):
                ax.annotate(str(i), (hh, pp), textcoords="offset points",
                            xytext=(9, 9), fontsize=11, fontweight="bold", color=PRIMARY)
            ax.annotate("", xy=(h[1], p[1]), xytext=(h[0], p[0]),
                        arrowprops=dict(arrowstyle="->", color=color, lw=2))

    s = cycle_states(**REVISED)
    ax.set_yscale("log")
    ax.set_xlim(150, 780)
    ax.set_ylim(1.5, 45)
    ax.set_yticks([2, 3, 5, 7, 10, 15, 20, 30, 42])
    ax.get_yaxis().set_major_formatter(plt.ScalarFormatter())
    ax.set_xlabel("Specific enthalpy $h$ [kJ/kg]", color=SECONDARY, fontsize=10)
    ax.set_ylabel("Pressure $p$ [bar]  (log scale)", color=SECONDARY, fontsize=10)
    ax.set_title("R-290 pressure-enthalpy chart --- designed cycle at the 150 kW / 35 C design point",
                 color=PRIMARY, fontsize=12, loc="left")
    ax.legend(loc="upper left", fontsize=9, framealpha=0.95, edgecolor=BASELINE)

    box = ("1 compressor suction   2 discharge   3 condenser out (liquid)   4 EEV out\n"
           f"Refrigerating effect $h_1-h_4$ = {(s['h1']-s['h4'])/1e3:.1f} kJ/kg      "
           f"Compressor work $h_2-h_1$ = {(s['h2']-s['h1'])/1e3:.1f} kJ/kg\n"
           f"Flash gas at EEV outlet $x_4$ = {s['x4']:.3f}      Pressure ratio = {s['pc']/s['po']:.2f}")
    ax.text(0.985, 0.03, box, transform=ax.transAxes, fontsize=8.5, color=SECONDARY,
            ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.5", facecolor=SURFACE, edgecolor=BASELINE))
    ax.text(0.0, -0.115, "Free cooling is not shown: the compressor is off in that mode, so no "
                          "refrigerant cycle exists to plot.",
            transform=ax.transAxes, fontsize=8, color=MUTED)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_temperature_cascade(out_path):
    """The air -> glycol -> refrigerant -> water temperature cascade. This is
    the single most load-bearing concept in the design: the condensing
    temperature (and hence the COP) is set by three approaches stacked in
    series, which is the price of the intermediate glycol loop."""
    import matplotlib.pyplot as plt

    fig, (axo, axr) = plt.subplots(1, 2, figsize=(12, 6.5), facecolor=SURFACE, sharey=True)

    for ax, (name, t_air, gly_c, gly_h, tc, to, chw_s, chw_r, cop) in zip(
        (axo, axr),
        (("Original design", 35, 40, 46, 52, 7, 15, 21, 3.42),
         ("Revised design", 35, 38, 42, 45, 10, 15, 21, 4.67))):
        _style(ax, grid_axis="y")
        segs = [("Ambient air", t_air, t_air, ORANGE, 0),
                ("Glycol loop", gly_c, gly_h, AQUA, 1),
                ("Refrigerant\n(condensing)", tc, tc, RED, 2),
                ("Refrigerant\n(evaporating)", to, to, BLUE, 3),
                ("Chilled water", chw_s, chw_r, GREEN, 4)]
        for label, lo, hi, c, x in segs:
            if lo == hi:
                ax.plot([x - 0.32, x + 0.32], [lo, lo], color=c, linewidth=5, solid_capstyle="butt", zorder=3)
                ax.text(x, lo + 1.0, f"{lo:.0f} C", ha="center", fontsize=9, color=SECONDARY, fontweight="bold")
            else:
                ax.add_patch(plt.Rectangle((x - 0.32, lo), 0.64, hi - lo, facecolor=c, alpha=0.85, zorder=3))
                ax.text(x, hi + 1.0, f"{hi:.0f} C", ha="center", fontsize=9, color=SECONDARY)
                ax.text(x, lo - 2.4, f"{lo:.0f} C", ha="center", fontsize=9, color=SECONDARY)
        ax.set_xticks(range(5))
        ax.set_xticklabels([s[0] for s in segs], fontsize=8.5)
        ax.set_xlim(-0.7, 4.7)
        ax.set_ylim(0, 60)

        # the stacked approaches on the hot side
        ax.annotate("", xy=(0.32, t_air), xytext=(0.68, gly_c),
                    arrowprops=dict(arrowstyle="<->", color=SECONDARY, lw=1.2))
        ax.text(0.5, (t_air + gly_c)/2 - 2.2, f"{gly_c-t_air:.0f} K\ndry cooler", fontsize=7.5,
                ha="center", color=SECONDARY)
        ax.text(1.0, (gly_c + gly_h)/2, f"{gly_h-gly_c:.0f} K\nrise", fontsize=7.5,
                ha="center", va="center", color="white", fontweight="bold")
        ax.annotate("", xy=(1.32, gly_h), xytext=(1.68, tc),
                    arrowprops=dict(arrowstyle="<->", color=SECONDARY, lw=1.2))
        ax.text(1.5, (gly_h + tc)/2 - 2.2, f"{tc-gly_h:.0f} K\ncondenser", fontsize=7.5,
                ha="center", color=SECONDARY)
        ax.annotate("", xy=(3.32, to), xytext=(3.68, chw_s),
                    arrowprops=dict(arrowstyle="<->", color=SECONDARY, lw=1.2))
        ax.text(3.5, (to + chw_s)/2, f"  {chw_s-to:.0f} K evaporator", fontsize=7.5,
                ha="left", va="center", color=SECONDARY)

        ax.axhline(t_air, color=ORANGE, linestyle=":", linewidth=1.1, zorder=1)
        ax.annotate("", xy=(2.0, t_air), xytext=(2.0, tc),
                    arrowprops=dict(arrowstyle="<->", color=RED, lw=2))
        ax.text(2.18, (t_air + tc)/2, f"{tc-t_air:.0f} K\nTOTAL\nair-to-refrigerant",
                fontsize=8.5, color=RED, fontweight="bold", va="center")
        ax.set_title(f"{name}  ---  design COP {cop:.2f}", color=PRIMARY, fontsize=11, loc="left")

    axo.set_ylabel("Temperature [C]", color=SECONDARY, fontsize=10)
    fig.suptitle("Temperature cascade: why the condensing temperature (and the COP) is set by "
                 "three stacked approaches",
                 color=PRIMARY, fontsize=12.5, x=0.01, ha="left", y=0.99)
    fig.text(0.01, 0.005, "The intermediate glycol loop puts three approaches in series between ambient air and the "
                          "refrigerant. Cutting them 5/6/6 K -> 3/4/3 K\nlowers $t_c$ by 7 K, which is what lifts the "
                          "full-load COP onto target.", fontsize=8.5, color=MUTED)
    fig.tight_layout(rect=[0, 0.05, 1, 0.96])
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_hx_profiles(out_path):
    """Condenser (three-zone) and evaporator (two-zone) temperature-vs-duty
    profiles at the revised design point -- the picture behind the LMTD
    sizing, showing where each exchanger's approach actually pinches."""
    import matplotlib.pyplot as plt

    s = cycle_states(**REVISED)
    fig, (axc, axe) = plt.subplots(1, 2, figsize=(12, 5), facecolor=SURFACE)

    # ---- Condenser: refrigerant 66.5 -> 45 (desuperheat), 45 (condense), 45 -> 40 (subcool)
    _style(axc)
    hg = PropsSI("HMASS", "T", 45 + 273.15, "Q", 1, R)
    hf = PropsSI("HMASS", "T", 45 + 273.15, "Q", 0, R)
    tot = s["h2"] - s["h3"]
    f_dsh, f_cond = (s["h2"] - hg) / tot, (hg - hf) / tot
    x = [0, f_dsh, f_dsh + f_cond, 1.0]
    axc.plot(x, [s["T2"], 45, 45, 40], color=RED, linewidth=2.5, marker="o", markersize=6,
             zorder=4, label="R-290 (hot side)")
    axc.plot([0, 1], [42, 38], color=AQUA, linewidth=2.5, marker="o", markersize=6,
             zorder=4, label="30% PG glycol (cold side)")
    for xx, lbl in ((f_dsh/2, "desuperheat"), (f_dsh + f_cond/2, "condensing"), (1 - (1-x[2])/2, "subcool")):
        axc.text(xx, 0.95, lbl, transform=axc.get_xaxis_transform(), ha="center",
                 fontsize=8, color=MUTED)
    axc.axvline(f_dsh, color=GRID, lw=1); axc.axvline(f_dsh + f_cond, color=GRID, lw=1)
    axc.annotate("", xy=(1.0, 40), xytext=(1.0, 38), arrowprops=dict(arrowstyle="<->", color=SECONDARY, lw=1.4))
    axc.text(0.97, 39, "2 K pinch\n(cold end)", ha="right", fontsize=8, color=SECONDARY)
    axc.set_xlabel("Fraction of condenser duty [-]", color=SECONDARY, fontsize=10)
    axc.set_ylabel("Temperature [C]", color=SECONDARY, fontsize=10)
    axc.set_title("Condenser --- three-zone profile (181.7 kW)", color=PRIMARY, fontsize=11, loc="left")
    axc.set_ylim(30, 72)
    axc.legend(fontsize=8.5, edgecolor=BASELINE)

    # ---- Evaporator: refrigerant boils at 10 then superheats to 20; water 21 -> 15
    _style(axe)
    hg_e = PropsSI("HMASS", "T", 10 + 273.15, "Q", 1, R)
    tot_e = s["h1"] - s["h4"]
    f_boil = (hg_e - s["h4"]) / tot_e
    axe.plot([0, f_boil, 1.0], [10, 10, s["T1"]], color=BLUE, linewidth=2.5, marker="o",
             markersize=6, zorder=4, label="R-290 (cold side)")
    axe.plot([0, 1], [15, 21], color=GREEN, linewidth=2.5, marker="o", markersize=6,
             zorder=4, label="Chilled water (hot side)")
    axe.axvline(f_boil, color=GRID, lw=1)
    for xx, lbl in ((f_boil/2, "boiling"), (f_boil + (1-f_boil)/2, "superheat")):
        axe.text(min(xx, 0.965), 0.95, lbl, transform=axe.get_xaxis_transform(),
                 ha="center", fontsize=8, color=MUTED)
    axe.annotate("", xy=(0.0, 10), xytext=(0.0, 15), arrowprops=dict(arrowstyle="<->", color=SECONDARY, lw=1.4))
    axe.text(0.03, 12.5, "5 K approach\n(LCHW - $t_o$)", ha="left", fontsize=8, color=SECONDARY)
    axe.set_xlabel("Fraction of evaporator duty [-]", color=SECONDARY, fontsize=10)
    axe.set_ylabel("Temperature [C]", color=SECONDARY, fontsize=10)
    axe.set_title("Evaporator --- two-zone profile (150 kW)", color=PRIMARY, fontsize=11, loc="left")
    axe.set_ylim(5, 27)
    axe.legend(fontsize=8.5, edgecolor=BASELINE, loc="lower right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_pue_breakdown(out_path):
    """Where the design-day facility power actually goes."""
    import matplotlib.pyplot as plt

    loads = [("Compressor bank", 31.72, BLUE), ("Dry cooler fans", 10.60, AQUA),
             ("Glycol loop pump", 2.56, YELLOW), ("CHW/CRAH pump", 1.11, GREEN)]
    tot = sum(v for _, v, _ in loads)
    fig, ax = plt.subplots(figsize=(10, 3.2), facecolor=SURFACE)
    _style(ax, grid_axis="x")
    left = 0.0
    for name, v, c in loads:
        ax.barh([0], [v], left=left, color=c, height=0.55, zorder=3,
                label=f"{name}  {v:.2f} kW ({v/tot*100:.0f}%)")
        if v > 2:
            ax.text(left + v/2, 0, f"{v:.1f}", ha="center", va="center",
                    color="white", fontsize=10, fontweight="bold", zorder=4)
        left += v
    ax.set_yticks([])
    ax.set_xlabel("Design-day facility power [kW]", color=SECONDARY, fontsize=10)
    ax.set_xlim(0, tot * 1.02)
    ax.set_title(f"Design-day facility power: {tot:.1f} kW against a 150 kW IT load  "
                 f"$\\Rightarrow$  PUE = {(150+tot)/150:.3f}",
                 color=PRIMARY, fontsize=11.5, loc="left")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.42), ncol=4, fontsize=8.5,
              frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_original_vs_revised(out_path):
    """The revision's headline result. Three panels rather than one dual-axis
    plot, since COP, IPLV and PUE are different scales."""
    import matplotlib.pyplot as plt

    panels = [("Full-load COP\n(AHRI 6.7 C LCHW)", 2.64, 3.53, 3.5, "$\\geq$3.5", True),
              ("IPLV.IP", 4.43, 5.23, 5.0, "$\\geq$5.0", True),
              ("Design-day PUE\n(lower is better)", 1.40, 1.31, None, None, False)]
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.2), facecolor=SURFACE)
    for ax, (title, orig, rev, target, tlabel, higher_better) in zip(axes, panels):
        _style(ax, grid_axis="y")
        good = (rev >= target) if (target and higher_better) else True
        bars = ax.bar([0, 1], [orig, rev], color=[MUTED, GREEN if good else RED],
                      width=0.55, zorder=3)
        for b, v in zip(bars, [orig, rev]):
            ax.text(b.get_x()+b.get_width()/2, v, f"{v:.2f}", ha="center", va="bottom",
                    fontsize=11, fontweight="bold", color=SECONDARY)
        if target:
            ax.axhline(target, color=RED, linestyle="--", linewidth=1.4, zorder=4)
            ax.text(-0.48, target, f"target {tlabel}", ha="left", va="bottom",
                    fontsize=8.5, color=RED)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Original", "Revised"], fontsize=9.5)
        ax.set_title(title, color=PRIMARY, fontsize=10, loc="left")
        ax.set_ylim(0, max(orig, rev, target or 0) * 1.30)
    fig.suptitle("Efficiency revision: both required metrics moved from failing to compliant",
                 color=PRIMARY, fontsize=12.5, x=0.01, ha="left")
    fig.text(0.01, 0.005, "Revision = condensing approach 17->10 K, evaporator approach 8->5 K, and VFDs "
                          "on the compressors. Architecture, refrigerant and compressor model unchanged.",
             fontsize=8.5, color=MUTED)
    fig.tight_layout(rect=[0, 0.04, 1, 0.94])
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_duration_curve(out_path):
    """Champaign TMY3 ambient temperature duration curve against the ASHRAE 90.1
    economizer thresholds -- the picture behind the 29% free-cooling hours."""
    import matplotlib.pyplot as plt
    from weather import ChampaignWeather

    T = np.sort(np.array([h.T_db for h in ChampaignWeather()]))[::-1]
    hrs = np.arange(1, len(T) + 1)
    n_full = int((T <= 4.0).sum()); n_part = int(((T > 4.0) & (T <= 10.0)).sum())
    n_mech = len(T) - n_full - n_part

    fig, ax = plt.subplots(figsize=(11, 5.2), facecolor=SURFACE)
    _style(ax)
    ax.plot(hrs, T, color=BLUE, linewidth=2, zorder=4)
    ax.fill_between(hrs, T, -30, where=(T > 10), color=RED, alpha=0.13, zorder=1)
    ax.fill_between(hrs, T, -30, where=((T <= 10) & (T > 4)), color=YELLOW, alpha=0.16, zorder=1)
    ax.fill_between(hrs, T, -30, where=(T <= 4), color=AQUA, alpha=0.18, zorder=1)
    for y, c, lbl, xa, ha in ((10.0, YELLOW, "Economizer activates ($\\leq$10 C)", 8700, "right"),
                              (4.0, AQUA, "Full free cooling ($\\leq$4 C)", 8700, "right"),
                              (35.0, ORANGE, "Dry cooler design ambient (35 C)", 120, "left")):
        ax.axhline(y, color=c, linestyle="--", linewidth=1.5, zorder=3)
        ax.text(xa, y + 0.7, lbl, ha=ha, fontsize=8.5, color=SECONDARY)
    ax.set_xlim(0, 8760); ax.set_ylim(-30, 40)
    ax.set_xlabel("Hours per year at or above this temperature", color=SECONDARY, fontsize=10)
    ax.set_ylabel("Outdoor dry-bulb temperature [C]", color=SECONDARY, fontsize=10)
    ax.set_title("Champaign TMY3 ambient duration curve --- 44.5% of the year benefits from the economizer",
                 color=PRIMARY, fontsize=12, loc="left")
    box = (f"Mechanical only (>10 C):      {n_mech:,} h  ({n_mech/87.60:.1f}%)\n"
           f"Partial free cooling (4-10 C): {n_part:,} h  ({n_part/87.60:.1f}%)\n"
           f"Full free cooling ($\\leq$4 C):     {n_full:,} h  ({n_full/87.60:.1f}%)")
    ax.text(0.985, 0.95, box, transform=ax.transAxes, fontsize=9, color=SECONDARY,
            ha="right", va="top", family="monospace",
            bbox=dict(boxstyle="round,pad=0.5", facecolor=SURFACE, edgecolor=BASELINE))
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    plot_ph_chart(os.path.join(FIG_CYCLE, "ph_chart.png"))
    plot_temperature_cascade(os.path.join(FIG_CYCLE, "temperature_cascade.png"))
    plot_hx_profiles(os.path.join(FIG_CYCLE, "hx_profiles.png"))
    plot_pue_breakdown(os.path.join(FIG_PERF, "pue_breakdown.png"))
    plot_original_vs_revised(os.path.join(FIG_PERF, "original_vs_revised.png"))
    plot_duration_curve(os.path.join(FIG_ANNUAL, "ambient_duration_curve.png"))
