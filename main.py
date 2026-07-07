"""Integration driver: chains the four refrigerant-loop components
(compressor bank -> condenser -> EEV -> evaporator) around ONE shared
mass flow, taken from the compressor bank's own solved operating point
-- instead of each component file guessing its own placeholder m_dot
independently.

The economizer is NOT part of this chain: it is a parallel, refrigerant-
free waterside loop (glycol <-> chilled water) that only runs when the
compressor is off, so it has no shared m_dot with the cycle below.
"""
import math
import time

from CoolProp.CoolProp import PropsSI

from compressor import CompressorBank
from condenser import Condenser
from dry_cooler import DryCooler
from dry_cooler_pump import GlycolLoopPump
from eev import EEV
from evaporator import Evaporator
from water_side_pump import CHWPump, TRIM_CURVES_1_5AD_1750RPM
from weather import ChampaignWeather

R = "R290"
T_c = 52 + 273.15         # condensing temperature [K] (raised from 45 C so the condenser's
                          # saturation/subcooled approach stays positive against the dry
                          # cooler's vendor-confirmed 40/46 C glycol loop -- see dry_cooler.py)
DT_SH = 10                # suction superheat assumed by the compressor [K]
SUBCOOLING = 5            # condenser design subcooling [K] (matches eev.py)
Q_TARGET = 150e3          # design cooling duty [W]
N_UNITS = 2               # project's dual 75 kW parallel-circuit design (compressor.py)

# Annual-simulation constants (used by run_annual_simulation() below)
T_AIR_DESIGN_C = 35.0
APPROACH_K = (T_c - 273.15) - T_AIR_DESIGN_C  # 17 K, fixed condensing approach (same
                                               # simplification already used in compressor.py's
                                               # IPLV.IP calc) -- holds tc-T_air constant across
                                               # all ambient conditions rather than re-solving the
                                               # dry cooler's eps-NTU model at each point
GLYCOL_APPROACH_K = 6.0    # tc -> condenser glycol-leaving-temp approach (52->46 C at design)
FREE_COOLING_THRESHOLD_C = 4.0  # ASHRAE 90.1 full free-cooling threshold (economizer.py)

p_c = PropsSI("P", "T", T_c, "Q", 1, R)

# ---------------------------------------------------------------------
# 1) Compressor bank -- solves for the evaporating temperature (and unit
#    staging) that hits the Q_TARGET design duty at this condensing
#    temperature, with hot-gas bypass covering anything below the bank's
#    turndown floor. This is what actually closes the 40%-shortfall gap
#    from the earlier single-fixed-compressor version of this file: T_o
#    is now DERIVED from the real two-unit vendor map, not assumed.
# ---------------------------------------------------------------------
bank = CompressorBank(n_units=N_UNITS, superheat_K=DT_SH, subcooling_K=SUBCOOLING,
                       to_bounds_C=(-10.0, 15.0))
bank_result = bank.solve_with_bypass(Q_TARGET, T_c - 273.15)

T_o = bank_result["to_C"] + 273.15
p_o = PropsSI("P", "T", T_o, "Q", 1, R)

m_dot = bank_result["m_dot_total_kgh"] / 3600.0   # TOTAL flow across all active units
per_unit = bank_result["per_unit"]                # single-unit state (identical parallel units)

# The vendor map only reports aggregate Q/P/m_dot, not a discharge state --
# approximate h_2 from suction enthalpy + specific work (P_w/m_dot), PER
# UNIT (identical parallel units share the same discharge state regardless
# of how many are active). ASSUMPTION: P_w is electrical input power
# (includes motor/mechanical losses), not all of which ends up in the
# refrigerant stream, so this somewhat overstates the true discharge
# enthalpy rise -- flagged pending a proper discharge-temperature
# correlation.
h_2 = per_unit["h_suction"] + per_unit["P_w"] / (per_unit["m_dot_kgh"] / 3600.0)
T_2 = PropsSI("T", "HMASS", h_2, "P", p_c, R)

# ---------------------------------------------------------------------
# 2) Condenser -- same corrugation-family placeholder geometry as the
#    evaporator (D_h, Lambda, beta, L), pending a real datasheet.
# ---------------------------------------------------------------------
D_h, Lambda, beta, L = 0.004, 0.005, 30, 0.5
b, L_w = D_h * 1.17 / 2, 0.2
N_cp_cond = 38
A_flow_cond = N_cp_cond * b * L_w
N_cp_glycol = 24

condenser = Condenser(m_dot_refrigerant=m_dot, p_in=p_c, refrigerant=R,
                       h_in=h_2, subcooling=SUBCOOLING,
                       T_glycol_in=30 + 273.15, T_glycol_out=35 + 273.15,
                       D_h=D_h, A_flow=A_flow_cond, L=L, beta=beta, Lambda=Lambda,
                       N_cp=N_cp_cond, b=b, L_w=L_w, N_cp_glycol=N_cp_glycol)

# ---------------------------------------------------------------------
# 3) EEV -- isenthalpic throttle from the condenser outlet to p_o
# ---------------------------------------------------------------------
eev = EEV(m_dot_refrigerant=m_dot, p_in=condenser.p_out, h_in=condenser.h_out,
          p_out=p_o, refrigerant=R)

# ---------------------------------------------------------------------
# 4) Evaporator -- Q is DERIVED from m_dot and the compressor's assumed
#    suction superheat, not assumed as a fixed 150 kW target. This is
#    what actually closes the loop: whatever the compressor's geometry
#    delivers is what the evaporator is asked to produce.
# ---------------------------------------------------------------------
h_out_target = PropsSI('HMASS', 'T', T_o + DT_SH, 'P', p_o, R)   # = compressor suction state
Q_actual = m_dot * (h_out_target - eev.h_out)

N_cp_evap = 38
A_flow_evap = N_cp_evap * b * L_w
N_cp_water = 24

evaporator = Evaporator(m_dot_refrigerant=m_dot, p_in=p_o, Q=Q_actual, refrigerant=R,
                         T_water_in=21 + 273.15, T_water_out=15 + 273.15, h_in=eev.h_out,
                         D_h=D_h, A_flow=A_flow_evap, L=L, beta=beta, Lambda=Lambda,
                         N_cp=N_cp_evap, b=b, L_w=L_w, N_cp_water=N_cp_water)

# ---------------------------------------------------------------------
# 5) Refrigerant-line piping sizing (project spec Table 1 velocity
#    guidelines: gaseous low-P 4.5-20 m/s, gaseous high-P 10-18 m/s,
#    liquid <1.5 m/s). Required ID is sized for a design velocity within
#    each band, then rounded UP to the nearest standard ACR copper tube.
#
#    NOTE: this sizing does NOT feed a pressure drop back into p_o/p_c or
#    any other cycle state. That's intentional, not an oversight -- the
#    project brief explicitly says to neglect piping pressure losses
#    (heat-exchanger losses are NOT to be neglected, and aren't: see the
#    Han/Martin correlations throughout this file's HX components). This
#    section is a standalone sizing/velocity check sitting alongside the
#    cycle, per spec.
# ---------------------------------------------------------------------
# Standard ACR (Type L) hard-drawn copper tube, actual ID [m]. Only the
# 1-5/8" and 2-5/8" entries are CONFIRMED against the project's own tube
# table; the rest are typical ACR dimensions included for selection
# generality -- ASSUMED, double-check against the actual catalog in use
# if a different size gets selected here.
ACR_TUBE_ID_M = {
    '3/8"':   0.00874,
    '1/2"':   0.01199,
    '5/8"':   0.01508,
    '3/4"':   0.01831,
    '7/8"':   0.02120,
    '1-1/8"': 0.02769,
    '1-3/8"': 0.03388,
    '1-5/8"': 0.03823,   # CONFIRMED
    '2-1/8"': 0.05087,
    '2-5/8"': 0.06287,   # CONFIRMED
    '3-1/8"': 0.07480,
    '3-5/8"': 0.08733,
}


def select_tube(m_dot_line, rho, v_design):
    """Required ID for a target design velocity, then the smallest
    standard ACR tube whose actual ID still clears it, and the resulting
    actual velocity in that tube."""
    ID_required = math.sqrt(4 * m_dot_line / (rho * v_design * math.pi))
    for size, ID in sorted(ACR_TUBE_ID_M.items(), key=lambda kv: kv[1]):
        if ID >= ID_required:
            v_actual = m_dot_line / (rho * math.pi / 4 * ID ** 2)
            return size, ID, ID_required, v_actual
    raise ValueError(f"No standard ACR tube large enough for ID_required={ID_required*1000:.1f} mm")


# Suction line: compressor suction state (T_o + DT_SH, p_o), gaseous low-pressure side
rho_suction = PropsSI('D', 'T', T_o + DT_SH, 'P', p_o, R)
tube_suction, ID_suction, ID_req_suction, v_suction = select_tube(m_dot, rho_suction, v_design=15.0)

# Discharge (hot gas) line: compressor discharge state (T_2, p_c), gaseous high-pressure side
rho_discharge = PropsSI('D', 'HMASS', h_2, 'P', p_c, R)
tube_discharge, ID_discharge, ID_req_discharge, v_discharge = select_tube(m_dot, rho_discharge, v_design=14.0)

# Liquid line: condenser subcooled outlet (T_out, p_out), liquid <1.5 m/s
rho_liquid = PropsSI('D', 'T', condenser.T_out, 'P', condenser.p_out, R)
tube_liquid, ID_liquid, ID_req_liquid, v_liquid = select_tube(m_dot, rho_liquid, v_design=1.45)

# ---------------------------------------------------------------------
# 6) Design-day PUE: (IT load + other facility loads) / IT load. IT load
#    is taken as numerically equal to the design cooling duty (150 kW) --
#    a common simplifying proxy in the absence of a separate data-hall IT
#    power model. "Other facility loads" = compressor bank + glycol loop
#    pump + dry cooler fans (vendor-confirmed rated value, mechanical-
#    mode design point) + CHW pump. The glycol pump uses the LIVE
#    condenser glycol flow/temps; the CHW pump uses the LIVE evaporator
#    water flow/pressure drop -- everything else (dry cooler design
#    point, free-cooling reference point, CRAH dP, piping) reuses the
#    same design-basis values as those components' own chapters.
# ---------------------------------------------------------------------
glycol_pump = GlycolLoopPump(
    m_dot_mech=condenser.m_dot_glycol, T_glycol_hot_mech=condenser.T_glycol_out,
    T_glycol_cold_mech=condenser.T_glycol_in,
    m_dot_free=6.52, T_glycol_hot_free=15 + 273.15, T_glycol_cold_free=9 + 273.15,
    dP_drycooler_mech_ref=48.26e3, dP_econ_free_ref=10.43e3,
    D_h=D_h, L=L, beta=beta, b=b, L_w=L_w, N_cp_glycol=N_cp_glycol,
    D_pipe=0.0779, pipe_roughness=0.045e-3, L_eq_total=65.0,
)

chw_pump = CHWPump(
    m_dot_water=evaporator.m_dot_water, dp_hx_kpa=evaporator.delta_p_water / 1e3,
    dp_crah_kpa=50.0, trim_curves_gpm_ft=TRIM_CURVES_1_5AD_1750RPM,
)
chw_pump.select_trim()

P_FAN_DESIGN = 10.74e3   # W, vendor-confirmed dry cooler fan power at the mechanical-mode
                          # design point (Kelvion ULF-PA104Y4V-096Z100, dry_cooler.py)

P_compressor = bank_result["P_total_w"]
P_glycol_pump = glycol_pump.P_elec_mech
# eta_pump=0.60 matches water_side_pump.py's own contour-line reading; eta_motor=0.90
# matches the glycol pump's motor-efficiency assumption (no CHW motor-efficiency
# stage built into CHWPump itself, so applied here explicitly)
P_chw_pump = chw_pump.shaft_power_w(eta_pump=0.60) / 0.90
P_other = P_compressor + P_glycol_pump + P_FAN_DESIGN + P_chw_pump
PUE = (Q_TARGET + P_other) / Q_TARGET

# ---------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------
print("=== SHARED CYCLE MASS FLOW ===")
print(f"m_dot (from compressor bank) : {m_dot:.3f} kg/s ({m_dot*3600:.0f} kg/h)")
print()
print("=== COMPRESSOR BANK ===")
print(f"Units active : {bank_result['n_active']} / {N_UNITS}  (status: {bank_result['status']})")
print(f"Evap temp (solved) : {bank_result['to_C']:.2f} C")
print(f"Discharge : {T_2-273.15:.1f} C @ {p_c/1000:.0f} kPa  (approximated, no map discharge state)")
print(f"Work      : {bank_result['P_total_w']/1e3:.2f} kW  (electrical input, per vendor map)")
print(f"Bank COP  : {bank_result['COP']:.2f}  (catalog-predicted duty: {bank_result['Q_total_w']/1e3:.1f} kW)")
print()
print("=== CONDENSER ===")
print(f"Duty        : {condenser.Q/1e3:.2f} kW")
print(f"Subcool     : {condenser.subcool:.2f} K @ {condenser.p_out/1000:.0f} kPa")
print(f"Delta_p     : {condenser.delta_p/1e3:.2f} kPa  (refrigerant side)")
print(f"Delta_p     : {condenser.delta_p_glycol/1e3:.2f} kPa  (glycol side, Martin 1996)")
print()
print("=== EEV ===")
print(f"Quality out : {eev.x_out:.3f}")
print(f"Kv required : {eev.Kv_required:.3f} m^3/h  (Loading vs E3V65: {eev.Kv_required/eev.Kv_rated*100:.1f}%)")
print()
print("=== EVAPORATOR ===")
print(f"Superheat  : {evaporator.superheat:.2f} K @ {evaporator.p_out/1000:.0f} kPa")
print(f"Delta_p    : {evaporator.delta_p/1e3:.2f} kPa  (refrigerant side)")
print(f"Delta_p    : {evaporator.delta_p_water/1e3:.2f} kPa  (water side, Martin 1996)")
print(f"Water flow : {evaporator.m_dot_water:.3f} kg/s")
print()
print("=== REFRIGERANT PIPING (Table 1 velocity guidelines) ===")
print(f"Suction   : {rho_suction:.2f} kg/m3  req.ID {ID_req_suction*1000:.1f} mm -> "
      f"{tube_suction} ({ID_suction*1000:.2f} mm ID)  actual v = {v_suction:.1f} m/s")
print(f"Discharge : {rho_discharge:.2f} kg/m3  req.ID {ID_req_discharge*1000:.1f} mm -> "
      f"{tube_discharge} ({ID_discharge*1000:.2f} mm ID)  actual v = {v_discharge:.1f} m/s")
print(f"Liquid    : {rho_liquid:.2f} kg/m3  req.ID {ID_req_liquid*1000:.1f} mm -> "
      f"{tube_liquid} ({ID_liquid*1000:.2f} mm ID)  actual v = {v_liquid:.2f} m/s")
print()
print("=== CYCLE CLOSURE CHECK ===")
print(f"Cooling duty delivered : {evaporator.Q/1e3:.1f} kW   (design target: {Q_TARGET/1e3:.0f} kW)")
print(f"Shortfall vs target    : {(1 - evaporator.Q/Q_TARGET)*100:.1f} %")
print(f"Energy balance check   : Q_cond ({condenser.Q/1e3:.1f} kW) "
      f"vs Q_evap + W_comp ({(evaporator.Q + bank_result['P_total_w'])/1e3:.1f} kW)")
print(f"Map cross-check        : catalog-predicted Q ({bank_result['Q_total_w']/1e3:.1f} kW) "
      f"vs cycle-derived Q ({evaporator.Q/1e3:.1f} kW)")
print()
print("=== DESIGN-DAY PUE ===")
print(f"IT load (proxy)      : {Q_TARGET/1e3:.1f} kW")
print(f"Compressor bank      : {P_compressor/1e3:.2f} kW")
print(f"Glycol loop pump     : {P_glycol_pump/1e3:.2f} kW")
print(f"Dry cooler fans      : {P_FAN_DESIGN/1e3:.2f} kW  (vendor-confirmed, design point)")
print(f"CHW/CRAH pump        : {P_chw_pump/1e3:.2f} kW")
print(f"Other facility loads : {P_other/1e3:.2f} kW")
print(f"PUE = (IT + other)/IT = {PUE:.3f}")


# =======================================================================
# ANNUAL SIMULATION (first pass): loops the real 8760-hour Champaign TMY3
# weather file, classifying each hour as mechanical (OAT > 4 C) or full
# free-cooling (OAT <= 4 C). SIMPLIFICATION: the partial-free-cooling band
# (4-10 C, per economizer.py's own regime table) is folded into mechanical
# mode for this first pass -- no partial-load economizer blend is modeled
# yet (same open item already flagged in the LaTeX report). This produces
# RAW per-hour results only (mode, compressor power, fan power, delivered
# duty); aggregating into annual energy / IPLV / PUE numbers is the next
# step, not attempted here.
# =======================================================================
def mechanical_hour(bank_obj, T_air_C, Q_target=Q_TARGET):
    """Lightweight per-hour solve for a mechanical-mode hour: floats the
    condensing temperature with ambient via the design point's fixed
    17 K approach (same simplification already used in compressor.py's
    IPLV.IP calculation), then solves the compressor bank for the
    required duty. Skips re-instantiating the full condenser/EEV/
    evaporator chain -- not needed for an annual energy sum, only for the
    single design-point validation already done above."""
    tc_C = T_air_C + APPROACH_K
    result = bank_obj.solve_with_bypass(Q_target, tc_C)
    return {
        "mode": "mechanical", "T_air_C": T_air_C, "tc_C": tc_C,
        "P_compressor_w": result["P_total_w"], "Q_delivered_w": result["Q_total_w"],
        "status": result["status"],
    }


def free_cooling_hour(T_air_C, Q_target=Q_TARGET):
    """Per-hour solve for a full-free-cooling hour (OAT <= 4 C): compressor
    off, economizer covers the full duty at its own fixed design point
    (9->15 C glycol, 21->15 C water) -- economizer.py's duty doesn't vary
    with ambient as long as the dry cooler can deliver cold-enough glycol,
    which is checked separately via the dry cooler's predict_off_design()
    in run_annual_simulation() below."""
    return {
        "mode": "free_cooling", "T_air_C": T_air_C, "tc_C": None,
        "P_compressor_w": 0.0, "Q_delivered_w": Q_target, "status": "ok",
    }


def run_annual_simulation():
    """Runs mechanical_hour()/free_cooling_hour() plus the dry cooler's
    off-design fan solve across every hour of the real embedded Champaign
    TMY3 weather file. Returns the list of per-hour result dicts (not yet
    aggregated -- that's the next step)."""
    wx = ChampaignWeather()
    bank_obj = CompressorBank(n_units=N_UNITS, superheat_K=DT_SH, subcooling_K=SUBCOOLING,
                               to_bounds_C=(-10.0, 15.0))
    dry_cooler_annual = DryCooler(
        Q_design=213.23e3, T_glycol_hot_in=46 + 273.15, T_glycol_cold_out=40 + 273.15,
        T_air_in_design=T_AIR_DESIGN_C + 273.15, dT_air_design=5.56,
    )
    P_FAN_RATED = 10.74e3  # W, vendor-confirmed (Kelvion ULF-PA104Y4V-096Z100)

    results = []
    t0 = time.time()
    for hr in wx:
        T_air_C = hr.T_db
        if T_air_C <= FREE_COOLING_THRESHOLD_C:
            r = free_cooling_hour(T_air_C)
            T_glycol_hot_in_hour = 15 + 273.15
        else:
            r = mechanical_hour(bank_obj, T_air_C)
            T_glycol_hot_in_hour = (r["tc_C"] - GLYCOL_APPROACH_K) + 273.15

        Q_cond_hour = r["Q_delivered_w"] + r["P_compressor_w"]
        fan_result = dry_cooler_annual.predict_off_design(
            Q_target=Q_cond_hour, T_glycol_hot_in=T_glycol_hot_in_hour, T_air_in=T_air_C + 273.15,
        )
        r["fan_speed_frac"] = fan_result["fan_speed_frac"]
        r["P_fan_w"] = dry_cooler_annual.fan_power(fan_result["fan_speed_frac"], P_FAN_RATED)
        r["hour_of_year"] = hr.hour_of_year
        results.append(r)

    elapsed = time.time() - t0
    n_mech = sum(1 for r in results if r["mode"] == "mechanical")
    n_free = len(results) - n_mech
    print(f"\nAnnual simulation: {len(results)} hours processed in {elapsed:.1f} s")
    print(f"  Mechanical-mode hours (incl. unmodeled partial-free-cooling band) : "
          f"{n_mech} ({n_mech/len(results)*100:.1f}%)")
    print(f"  Full free-cooling hours (OAT <= {FREE_COOLING_THRESHOLD_C:.0f} C)             : "
          f"{n_free} ({n_free/len(results)*100:.1f}%)")
    n_over_speed = sum(1 for r in results if r["fan_speed_frac"] > 1.0)
    if n_over_speed:
        print(f"  WARNING: {n_over_speed} hour(s) require >100% fan speed to hit duty -- "
              f"dry cooler capacity-ceiling hours, flagged for follow-up.")
    return results


def plot_annual_cop_load(results, save_path="annual_cop_load.png"):
    """Saves a two-panel PNG (delivered load, compressor COP) across all
    8760 hours, hour-of-year on a shared x-axis. Two stacked panels
    instead of one dual-axis plot, since load [kW] and COP [-] are
    different units/scales -- a dual y-axis chart is never the right call.

    Load is delivered duty (essentially flat at the 150 kW design target
    year-round, per the project's own "near-constant IT load" assumption
    -- see Project3_Description.txt); COP = Q_delivered/P_compressor is
    only defined for mechanical-mode hours (compressor running) and is
    left as a gap during free-cooling hours (compressor off), which are
    shaded in the background on both panels so the gaps read as
    intentional, not missing data.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    hours = np.array([r["hour_of_year"] for r in results])
    load_kw = np.array([r["Q_delivered_w"] for r in results]) / 1e3
    is_free = np.array([r["mode"] == "free_cooling" for r in results])
    cop = np.array([
        r["Q_delivered_w"] / r["P_compressor_w"] if r["P_compressor_w"] > 0 else np.nan
        for r in results
    ])

    # Colors: validated categorical palette (dataviz skill reference palette)
    BLUE = "#2a78d6"
    AQUA = "#1baf7a"
    SURFACE = "#fcfcfb"
    PRIMARY_INK = "#0b0b0b"
    SECONDARY_INK = "#52514e"
    MUTED = "#898781"
    GRID = "#e1e0d9"
    BASELINE = "#c3c2b7"

    # Approximate month tick positions (365-day nominal TMY3 year, matches weather.py's 8760 rows)
    month_starts = [0, 744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016]
    month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                     "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    # Contiguous free-cooling hour ranges, for the shared background shading
    bands = []
    band_start = None
    for h, f in zip(hours, is_free):
        if f and band_start is None:
            band_start = h
        elif not f and band_start is not None:
            bands.append((band_start, h))
            band_start = None
    if band_start is not None:
        bands.append((band_start, hours[-1]))

    fig, (ax_load, ax_cop) = plt.subplots(2, 1, figsize=(11, 6), sharex=True, facecolor=SURFACE)

    for ax in (ax_load, ax_cop):
        ax.set_facecolor(SURFACE)
        ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color(BASELINE)
        ax.spines["bottom"].set_color(BASELINE)
        ax.tick_params(colors=MUTED, labelsize=8)
        for lo, hi in bands:
            ax.axvspan(lo, hi, color=AQUA, alpha=0.10, linewidth=0, zorder=0)

    ax_load.plot(hours, load_kw, color=BLUE, linewidth=1.2, zorder=2)
    ax_load.set_ylabel("Cooling load [kW]", color=SECONDARY_INK, fontsize=9)
    ax_load.set_title("Delivered cooling load and compressor COP across the year "
                       "(Champaign, IL TMY3)", color=PRIMARY_INK, fontsize=11, loc="left")

    ax_cop.plot(hours, cop, color=BLUE, linewidth=1.0, zorder=2)
    ax_cop.set_ylabel("COP [-]", color=SECONDARY_INK, fontsize=9)
    ax_cop.set_xlabel("Month", color=SECONDARY_INK, fontsize=9)
    ax_cop.set_xlim(0, 8760)
    ax_cop.set_xticks(month_starts)
    ax_cop.set_xticklabels(month_labels)

    ax_cop.text(0.0, -0.34, "Aqua shading = free-cooling hours (compressor off, COP undefined)",
                transform=ax_cop.transAxes, fontsize=8, color=MUTED)

    fig.tight_layout()
    fig.savefig(save_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"\nSaved annual COP/load chart -> {save_path}")


print()
print("=" * 72)
print("ANNUAL SIMULATION -- first pass (raw per-hour results, no aggregation yet)")
print("=" * 72)
annual_results = run_annual_simulation()
plot_annual_cop_load(annual_results)
