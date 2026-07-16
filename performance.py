"""System performance metrics (project spec subtask 1f): full-load COP at
both the strict AHRI 550/590-2023 rating point and the project's own design
point, IPLV.IP, and the design-day PUE -- plus the figures and LaTeX tables
that carry those numbers into the report.

The monthly/seasonal weather-driven performance (also part of subtask 1f) is
produced separately by main.py's run_monthly_seasonal_reports(); this module
covers the standardized rating-point metrics that main.py does not, and reads
main.py's monthly CSV back in to build a single seasonal-summary table.

All compressor performance comes from the SAME Bitzer 4FEP-35Z AHRI-540 map
used everywhere else (compressor.py), so these rating numbers are consistent
with the design-point cycle in main.py rather than a separate calculation.
"""
import csv
import os

from compressor import Compressor, CompressorBank

_HERE = os.path.dirname(os.path.abspath(__file__))
REPORT_DIR = os.path.join(_HERE, "report")
FIG_DIR = os.path.join(REPORT_DIR, "figures", "performance", "design_point")
TABLE_DIR = os.path.join(REPORT_DIR, "tables", "performance")
MONTHLY_CSV = os.path.join(TABLE_DIR, "monthly", "monthly_report.csv")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)

# ---------------------------------------------------------------------
# Design-basis approaches -- REVISED DESIGN (efficiency compliance, see
# the report's Performance chapter and main.py's revision header):
#   * Evaporator cold-end approach LCHW - to = 5 K (was 8 K) -- more
#     evaporator plates.
#   * Condenser air-to-refrigerant approach tc - T_air = 10 K (was 17 K):
#     dry cooler 3 K + glycol rise 4 K + condenser 3 K -- larger dry
#     cooler and condenser area, same hardware families.
#   * A VFD on the compressors replaces hot-gas bypass below the turndown
#     floor (compressor.solve_with_vfd()).
# Original-design values (8 K / 17 K, no VFD) gave COP 2.64 / IPLV 4.43.
# ---------------------------------------------------------------------
APPROACH_EVAP_K = 5.0
APPROACH_COND_K = 10.0
SUPERHEAT_K = 8.0
SUBCOOLING_K = 5.0
Q_DESIGN_W = 150e3
COP_TARGET = 3.5
IPLV_TARGET = 5.0

# Design-day PUE breakdown (from main.py's REVISED design-point run; restated
# here so the performance chapter is self-contained). These are the
# mechanical-mode design-point electrical loads at 35 C ambient / 150 kW full
# load. Compressor power reflects the VFD-trimmed operating point (to=10 C,
# tc=45 C, COP 4.67); fan power is the fan-law estimate for the re-selected
# Kelvion ULF-PA106K4V-091F095 (see main.py's P_FAN_DESIGN note -- the single
# largest remaining uncertainty in this figure).
PUE_LOADS_KW = {
    "Compressor bank": 31.72,
    "Glycol loop pump": 2.56,
    "Dry cooler fans": 10.60,
    "CHW/CRAH pump": 1.11,
}


def full_load_cop():
    """Full-load COP at (a) the strict AHRI 550/590-2023 rating point
    (6.7 C LCHW, 35 C entering air) and (b) the project's own design point
    (15 C LCHW, 35 C entering air). Both share the same 35 C entering air
    (hence tc = 52 C via the 17 K approach); they differ only in chilled-
    water temperature, which sets the evaporating temperature via the 8 K
    evaporator approach."""
    comp = Compressor()

    to_ahri = 6.7 - APPROACH_EVAP_K       # colder AHRI chilled water -> lower to
    to_design = 15.0 - APPROACH_EVAP_K    # project's warmer W17 chilled water
    tc = 35.0 + APPROACH_COND_K

    p_ahri = comp.performance(to_ahri, tc, SUPERHEAT_K, SUBCOOLING_K)
    p_design = comp.performance(to_design, tc, SUPERHEAT_K, SUBCOOLING_K)

    return {
        "ahri": {"to_C": to_ahri, "tc_C": tc, "lchw_C": 6.7,
                 "Q_kW": 2 * p_ahri["Q_w"] / 1e3, "P_kW": 2 * p_ahri["P_w"] / 1e3,
                 "COP": p_ahri["COP"]},
        "design": {"to_C": to_design, "tc_C": tc, "lchw_C": 15.0,
                   "Q_kW": 2 * p_design["Q_w"] / 1e3, "P_kW": 2 * p_design["P_w"] / 1e3,
                   "COP": p_design["COP"]},
    }


def iplv():
    """IPLV.IP per AHRI 550/590-2023 S5.2: 0.01*A + 0.42*B + 0.45*C + 0.12*D,
    where A/B/C/D are COPs at 100/75/50/25 % load, each at its own reduced
    condenser-entering-air temperature (35/26.7/18.3/12.8 C). Each air temp
    is mapped to a refrigerant condensing temperature via the same fixed
    10 K approach (design revision). Below the fixed-speed turndown floor
    the VFD modulates capacity instead of hot-gas bypassing
    (compressor.solve_with_vfd(), conservative floor-point COP)."""
    bank = CompressorBank(n_units=2, superheat_K=SUPERHEAT_K, subcooling_K=SUBCOOLING_K)
    points = [
        ("A", 1.00, 0.01, 35.0),
        ("B", 0.75, 0.42, 26.7),
        ("C", 0.50, 0.45, 18.3),
        ("D", 0.25, 0.12, 12.8),
    ]
    rows = []
    value = 0.0
    for label, load_frac, weight, t_air_c in points:
        tc = t_air_c + APPROACH_COND_K
        r = bank.solve_with_vfd(load_frac * Q_DESIGN_W, tc)
        value += weight * r["COP"]
        rows.append({
            "label": label, "load_frac": load_frac, "weight": weight,
            "t_air_C": t_air_c, "tc_C": tc, "COP": r["COP"],
            "vfd": r["status"] == "below_floor_vfd",
        })
    return {"rows": rows, "iplv": value}


def tc_sensitivity():
    """Full-load COP at the AHRI chilled-water condition (6.7 C LCHW) as the
    condensing temperature varies -- shows why the original design's 17 K
    air-to-refrigerant approach (tc = 52 C) missed the 3.5 target and how the
    revised 10 K approach (tc = 45 C) reaches it. Evaluated at the ORIGINAL
    8 K evaporator approach so the curve isolates the condensing-side effect;
    the revised design's extra margin comes from the 5 K evaporator approach
    on top of this curve."""
    comp = Compressor()
    to_ahri = 6.7 - 8.0   # original evaporator approach, isolating the tc effect
    rows = []
    for tc in [42, 45, 48, 50, 52]:
        p = comp.performance(to_ahri, tc, SUPERHEAT_K, SUBCOOLING_K)
        rows.append({"tc_C": tc, "approach_K": tc - 35.0, "COP": p["COP"]})
    return rows


def design_day_pue():
    other = sum(PUE_LOADS_KW.values())
    pue = (Q_DESIGN_W / 1e3 + other) / (Q_DESIGN_W / 1e3)
    return {"loads": dict(PUE_LOADS_KW), "P_other_kW": other, "PUE": pue}


def read_seasonal_from_monthly():
    """Groups main.py's 12 monthly rows into the 4 standard seasons and
    reports each season's average total-system COP (= 1/(PUE-1)) and PUE,
    plus how many of its months run mechanical vs free-cooling."""
    if not os.path.exists(MONTHLY_CSV):
        return None
    with open(MONTHLY_CSV) as f:
        months = list(csv.DictReader(f))

    season_of = {12: "Winter", 1: "Winter", 2: "Winter", 3: "Spring", 4: "Spring",
                 5: "Spring", 6: "Summer", 7: "Summer", 8: "Summer", 9: "Fall",
                 10: "Fall", 11: "Fall"}
    order = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    seasons = {}
    for i, row in enumerate(months):
        s = season_of[i + 1]
        pue = float(row["PUE"])
        cop_total = 1.0 / (pue - 1.0)
        seasons.setdefault(s, {"pue": [], "cop": [], "mech": 0, "free": 0, "T": []})
        seasons[s]["pue"].append(pue)
        seasons[s]["cop"].append(cop_total)
        seasons[s]["T"].append(float(row["T_air_C"]))
        if row["mode"] == "mechanical":
            seasons[s]["mech"] += 1
        else:
            seasons[s]["free"] += 1

    out = []
    for s in ["Winter", "Spring", "Summer", "Fall"]:
        d = seasons[s]
        out.append({
            "season": s,
            "T_air_C": sum(d["T"]) / len(d["T"]),
            "cop_total": sum(d["cop"]) / len(d["cop"]),
            "pue": sum(d["pue"]) / len(d["pue"]),
            "mech_months": d["mech"], "free_months": d["free"],
        })
    return out


# ============================ FIGURES ================================
_BLUE, _AQUA, _YELLOW, _GREEN, _RED = "#2a78d6", "#1baf7a", "#eda100", "#008300", "#e34948"
_SURFACE, _PRIMARY, _SECONDARY = "#fcfcfb", "#0b0b0b", "#52514e"
_MUTED, _GRID, _BASELINE = "#898781", "#e1e0d9", "#c3c2b7"


def _style(ax):
    ax.set_facecolor(_SURFACE)
    ax.grid(True, axis="y", color=_GRID, linewidth=0.8, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(_BASELINE)
    ax.spines["bottom"].set_color(_BASELINE)
    ax.tick_params(colors=_MUTED, labelsize=9)


def plot_full_load_cop(fl, out_path):
    import numpy as np
    import matplotlib.pyplot as plt

    labels = ["AHRI 550/590\nrating point\n(6.7 C LCHW)", "Project design\npoint\n(15 C LCHW)"]
    cops = [fl["ahri"]["COP"], fl["design"]["COP"]]
    colors = [_RED if c < COP_TARGET else _GREEN for c in cops]
    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(7, 5), facecolor=_SURFACE)
    _style(ax)
    bars = ax.bar(x, cops, color=colors, width=0.5, zorder=2)
    ax.axhline(COP_TARGET, color=_SECONDARY, linewidth=1.4, linestyle="--", zorder=3)
    ax.text(-0.42, COP_TARGET + 0.05, f"Target COP $\\geq$ {COP_TARGET}",
            ha="left", va="bottom", fontsize=9, color=_SECONDARY)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Full-load COP [-]", color=_SECONDARY, fontsize=10)
    ax.set_title("Full-load COP: AHRI rating point vs. project design point",
                 color=_PRIMARY, fontsize=11, loc="left")
    for b, v in zip(bars, cops):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.2f}", ha="center",
                va="bottom", fontsize=10, color=_SECONDARY)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, facecolor=_SURFACE)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_iplv(ip, out_path):
    import numpy as np
    import matplotlib.pyplot as plt

    rows = ip["rows"]
    labels = [f"{r['label']}\n{int(r['load_frac']*100)}% load\n{r['t_air_C']:.1f} C air" for r in rows]
    cops = [r["COP"] for r in rows]
    colors = [_YELLOW if r["vfd"] else _BLUE for r in rows]
    x = np.arange(len(rows))

    fig, ax = plt.subplots(figsize=(8, 5), facecolor=_SURFACE)
    _style(ax)
    bars = ax.bar(x, cops, color=colors, width=0.55, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Part-load COP [-]", color=_SECONDARY, fontsize=10)
    ax.set_title(f"IPLV.IP part-load points (weighted IPLV.IP = {ip['iplv']:.2f}, "
                 f"target $\\geq$ {IPLV_TARGET})", color=_PRIMARY, fontsize=11, loc="left")
    for b, r in zip(bars, rows):
        ax.text(b.get_x() + b.get_width() / 2, r["COP"],
                f"{r['COP']:.2f}\n(w={r['weight']})", ha="center", va="bottom",
                fontsize=8, color=_SECONDARY)
    ax.text(0.0, -0.20, "Yellow = VFD speed modulation below the fixed-speed turndown floor "
                        "(conservative floor-point COP)",
            transform=ax.transAxes, fontsize=8, color=_MUTED)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, facecolor=_SURFACE)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_tc_sensitivity(rows, out_path):
    import matplotlib.pyplot as plt

    tc = [r["tc_C"] for r in rows]
    cop = [r["COP"] for r in rows]

    fig, ax = plt.subplots(figsize=(8, 5), facecolor=_SURFACE)
    _style(ax)
    ax.grid(True, axis="both", color=_GRID, linewidth=0.8, zorder=0)
    ax.plot(tc, cop, color=_BLUE, linewidth=2, marker="o", markersize=7, zorder=3)
    ax.axhline(COP_TARGET, color=_RED, linewidth=1.4, linestyle="--", zorder=2)
    ax.text(46.5, COP_TARGET - 0.045, f"Target COP $\\geq$ {COP_TARGET}",
            ha="left", va="top", fontsize=9, color=_RED)
    ax.axvline(52, color=_MUTED, linewidth=1.0, linestyle=":", zorder=2)
    ax.text(51.85, min(cop) + 0.02, "original design (52 C)", ha="right", va="bottom",
            rotation=90, fontsize=8, color=_MUTED)
    ax.axvline(45, color=_GREEN, linewidth=1.2, linestyle=":", zorder=2)
    ax.text(44.85, min(cop) + 0.02, "revised design (45 C)", ha="right", va="bottom",
            rotation=90, fontsize=8, color=_GREEN)
    ax.set_xlabel("Condensing temperature $t_c$ [C]  (air-to-refrigerant approach = $t_c$ - 35)",
                  color=_SECONDARY, fontsize=10)
    ax.set_ylabel("Full-load COP at AHRI LCHW [-]", color=_SECONDARY, fontsize=10)
    ax.set_title("Full-load COP sensitivity to condensing temperature (AHRI 6.7 C LCHW)",
                 color=_PRIMARY, fontsize=11, loc="left")
    for x, y in zip(tc, cop):
        ax.text(x, y + 0.04, f"{y:.2f}", ha="center", va="bottom", fontsize=8, color=_SECONDARY)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, facecolor=_SURFACE)
    plt.close(fig)
    print(f"Saved {out_path}")


# ============================ TABLES =================================
def write_full_load_table(fl, out_path):
    lines = [
        r"\begin{table}[h!]",
        r"\caption{Full-load COP at the AHRI 550/590-2023 rating point and the project design point (both at 35\,\textdegree C entering air, $t_c=52$\,\textdegree C).}",
        r"\label{tab:full_load_cop}",
        r"\begin{tabularx}{\linewidth}{lLLLLL}",
        r"\toprule",
        r"\textbf{Rating basis} & \textbf{LCHW [\textdegree C]} & \textbf{$t_o$ [\textdegree C]} & \textbf{$t_c$ [\textdegree C]} & \textbf{Capacity [kW]} & \textbf{COP [-]} \\",
        r"\midrule",
        f"AHRI 550/590 rating & {fl['ahri']['lchw_C']:.1f} & {fl['ahri']['to_C']:.1f} & {fl['ahri']['tc_C']:.0f} & {fl['ahri']['Q_kW']:.1f} & {fl['ahri']['COP']:.2f} \\\\",
        f"Project design point & {fl['design']['lchw_C']:.1f} & {fl['design']['to_C']:.1f} & {fl['design']['tc_C']:.0f} & {fl['design']['Q_kW']:.1f} & {fl['design']['COP']:.2f} \\\\",
        r"\bottomrule",
        r"\end{tabularx}",
        r"\end{table}",
    ]
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Saved {out_path}")


def write_iplv_table(ip, out_path):
    lines = [
        r"\begin{table}[h!]",
        r"\caption{IPLV.IP part-load points (AHRI 550/590-2023 \S5.2). Each point is evaluated at its own reduced condenser-entering-air temperature.}",
        r"\label{tab:iplv}",
        r"\begin{tabularx}{\linewidth}{lLLLLL}",
        r"\toprule",
        r"\textbf{Point} & \textbf{Load [\%]} & \textbf{Weight} & \textbf{Entering air [\textdegree C]} & \textbf{$t_c$ [\textdegree C]} & \textbf{COP [-]} \\",
        r"\midrule",
    ]
    for r in ip["rows"]:
        note = r" $^{\dagger}$" if r["vfd"] else ""
        lines.append(f"{r['label']} & {int(r['load_frac']*100)} & {r['weight']} & "
                     f"{r['t_air_C']:.1f} & {r['tc_C']:.1f} & {r['COP']:.2f}{note} \\\\")
    lines += [
        r"\midrule",
        rf"\multicolumn{{5}}{{l}}{{\textbf{{IPLV.IP}} $= 0.01A + 0.42B + 0.45C + 0.12D$}} & \textbf{{{ip['iplv']:.2f}}} \\",
        r"\bottomrule",
        r"\end{tabularx}",
        r"{\small $^{\dagger}$ VFD speed modulation below the fixed-speed turndown floor "
        r"(COP held at the floor point's map value --- conservative, since reduced speed "
        r"lets the evaporating temperature float upward).}",
        r"\end{table}",
    ]
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Saved {out_path}")


def write_pue_table(pue, out_path):
    lines = [
        r"\begin{table}[h!]",
        r"\caption{Design-day facility PUE at the 35\,\textdegree C / 150\,kW full-load condition.}",
        r"\label{tab:pue}",
        r"\begin{tabularx}{\linewidth}{lL}",
        r"\toprule",
        r"\textbf{Electrical load} & \textbf{Power [kW]} \\",
        r"\midrule",
        r"IT load (proxy for cooling duty) & 150.00 \\",
        r"\midrule",
    ]
    for name, kw in pue["loads"].items():
        lines.append(f"{name} & {kw:.2f} \\\\")
    lines += [
        r"\midrule",
        f"Total other facility loads & {pue['P_other_kW']:.2f} \\\\",
        rf"\textbf{{PUE}} $=$ (IT $+$ other)$/$IT & \textbf{{{pue['PUE']:.3f}}} \\",
        r"\bottomrule",
        r"\end{tabularx}",
        r"\end{table}",
    ]
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Saved {out_path}")


def write_seasonal_summary_table(seasons, out_path):
    if seasons is None:
        return
    lines = [
        r"\begin{table}[h!]",
        r"\caption{Seasonal-average system performance (Champaign, IL TMY3), grouped from the 12 monthly-average runs of main.py.}",
        r"\label{tab:seasonal_perf_summary}",
        r"\begin{tabularx}{\linewidth}{lLLLL}",
        r"\toprule",
        r"\textbf{Season} & \textbf{Avg. $T_{air}$ [\textdegree C]} & \textbf{Total-system COP [-]} & \textbf{PUE [-]} & \textbf{Mechanical / free-cooling months} \\",
        r"\midrule",
    ]
    for s in seasons:
        lines.append(f"{s['season']} & {s['T_air_C']:.1f} & {s['cop_total']:.1f} & "
                     f"{s['pue']:.2f} & {s['mech_months']} / {s['free_months']} \\\\")
    lines += [
        r"\bottomrule",
        r"\end{tabularx}",
        r"\end{table}",
    ]
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    fl = full_load_cop()
    ip = iplv()
    sens = tc_sensitivity()
    pue = design_day_pue()
    seasons = read_seasonal_from_monthly()

    print("=== FULL-LOAD COP ===")
    for k in ("ahri", "design"):
        d = fl[k]
        print(f"  {k:<8}: LCHW={d['lchw_C']:.1f}C to={d['to_C']:.1f}C tc={d['tc_C']:.0f}C  "
              f"Q={d['Q_kW']:.1f}kW  COP={d['COP']:.3f}")
    print(f"\n=== IPLV.IP = {ip['iplv']:.2f}  (target >= {IPLV_TARGET}) ===")
    for r in ip["rows"]:
        print(f"  {r['label']}: {int(r['load_frac']*100)}% @ {r['t_air_C']:.1f}C air -> "
              f"tc={r['tc_C']:.1f}C  COP={r['COP']:.2f}{'  [vfd]' if r['vfd'] else ''}")
    print(f"\n=== DESIGN-DAY PUE = {pue['PUE']:.3f} ===")

    plot_full_load_cop(fl, os.path.join(FIG_DIR, "full_load_cop.png"))
    plot_iplv(ip, os.path.join(FIG_DIR, "iplv_points.png"))
    plot_tc_sensitivity(sens, os.path.join(FIG_DIR, "cop_tc_sensitivity.png"))
    write_full_load_table(fl, os.path.join(TABLE_DIR, "full_load_cop.tex"))
    write_iplv_table(ip, os.path.join(TABLE_DIR, "iplv.tex"))
    write_pue_table(pue, os.path.join(TABLE_DIR, "design_day_pue.tex"))
    write_seasonal_summary_table(seasons, os.path.join(TABLE_DIR, "seasonal_perf_summary.tex"))
