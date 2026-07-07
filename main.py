"""Integration driver: chains the four refrigerant-loop components
(compressor bank -> condenser -> EEV -> evaporator) around ONE shared
mass flow, taken from the compressor bank's own solved operating point
-- instead of each component file guessing its own placeholder m_dot
independently.

The economizer is NOT part of this chain: it is a parallel, refrigerant-
free waterside loop (glycol <-> chilled water) that only runs when the
compressor is off, so it has no shared m_dot with the cycle below.
"""
import csv
import math
import time

import numpy as np
from CoolProp.CoolProp import PropsSI

from compressor import CompressorBank
from condenser import Condenser
from dry_cooler import DryCooler
from dry_cooler_pump import GlycolLoopPump
from economizer import Economizer
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
ECON_L, ECON_L_W = 0.70, 0.25    # economizer's own larger frame (economizer.py's design point)

p_c = PropsSI("P", "T", T_c, "Q", 1, R)

# Shared dry cooler instance (same vendor-validated design point everywhere it's used --
# predict_off_design() doesn't mutate instance state, so reusing one object across the
# design-point report, the annual loop, and the monthly/seasonal report is safe).
DRY_COOLER = DryCooler(
    Q_design=213.23e3, T_glycol_hot_in=46 + 273.15, T_glycol_cold_out=40 + 273.15,
    T_air_in_design=T_AIR_DESIGN_C + 273.15, dT_air_design=5.56,
)

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


def select_tube(m_dot_line, rho, v_min, v_max, v_target):
    """Picks the standard ACR tube whose resulting velocity best fits
    inside the Table 1 guideline band [v_min, v_max] (closest to
    v_target within the band). Tube IDs are discrete, so a single
    design-velocity threshold (the old approach) can leave the actual
    velocity outside the guideline band once mass flow changes -- this
    checks the achieved band membership directly instead. Falls back to
    the standard size that minimizes the band violation if no single
    size satisfies both bounds at once (possible when the required ID
    falls in a gap between standard sizes)."""
    ID_target = math.sqrt(4 * m_dot_line / (rho * v_target * math.pi))
    candidates = [(size, ID, m_dot_line / (rho * math.pi / 4 * ID ** 2))
                  for size, ID in ACR_TUBE_ID_M.items()]

    in_band = [c for c in candidates if v_min <= c[2] <= v_max]
    if in_band:
        size, ID, v_actual = min(in_band, key=lambda c: abs(c[2] - v_target))
    else:
        def violation(c):
            v = c[2]
            return (v_min - v) if v < v_min else (v - v_max)
        size, ID, v_actual = min(candidates, key=violation)

    return size, ID, ID_target, v_actual


# Suction line: compressor suction state (T_o + DT_SH, p_o), gaseous low-pressure side
# Table 1 guideline band: 4.5-20 m/s
rho_suction = PropsSI('D', 'T', T_o + DT_SH, 'P', p_o, R)
tube_suction, ID_suction, ID_req_suction, v_suction = select_tube(
    m_dot, rho_suction, v_min=4.5, v_max=20.0, v_target=15.0)

# Discharge (hot gas) line: compressor discharge state (T_2, p_c), gaseous high-pressure side
# Table 1 guideline band: 10-18 m/s
rho_discharge = PropsSI('D', 'HMASS', h_2, 'P', p_c, R)
tube_discharge, ID_discharge, ID_req_discharge, v_discharge = select_tube(
    m_dot, rho_discharge, v_min=10.0, v_max=18.0, v_target=14.0)

# Liquid line: condenser subcooled outlet (T_out, p_out) -- Table 1 guideline: <1.5 m/s
rho_liquid = PropsSI('D', 'T', condenser.T_out, 'P', condenser.p_out, R)
tube_liquid, ID_liquid, ID_req_liquid, v_liquid = select_tube(
    m_dot, rho_liquid, v_min=0.0, v_max=1.5, v_target=1.2)

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
        fan_result = DRY_COOLER.predict_off_design(
            Q_target=Q_cond_hour, T_glycol_hot_in=T_glycol_hot_in_hour, T_air_in=T_air_C + 273.15,
        )
        r["fan_speed_frac"] = fan_result["fan_speed_frac"]
        r["P_fan_w"] = DRY_COOLER.fan_power(fan_result["fan_speed_frac"], P_FAN_RATED)
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


# =======================================================================
# MONTHLY / SEASONAL DETAILED REPORT: runs the FULL cycle (same physics as
# the single design-point section at the top of this file -- compressor
# bank, condenser, EEV, evaporator, refrigerant piping sizing, dry cooler
# fan, both pumps, PUE) at each month's and each season's AVERAGE dry-bulb
# temperature, instead of the lightweight compressor-only solve used for
# the 8760-hour annual loop above. 16 total runs (12 months + 4 seasons),
# so the extra per-run detail (full Condenser/EEV/Evaporator instances,
# 200-point Han/Martin integration, etc.) is cheap here.
#
# Reuses the module-level glycol_pump/chw_pump built by the design-point
# section: their trim/selection depends only on the two FIXED reference
# duty points (mechanical 9.07 kg/s @ 40/46 C, free-cooling 6.52 kg/s @
# 9/15 C, and the evaporator's own worst-case-branch water-side dP for the
# CHW pump) -- none of that changes with ambient, so P_elec_mech/
# P_elec_free/P_chw_pump are constants pulled from the SAME already-built
# instances rather than re-selected per period.
#
# INHERITED SIMPLIFICATIONS (same ones already flagged elsewhere in this
# file, not expanded here): the condenser's glycol temperatures stay fixed
# at 30/35 C regardless of month (the known, still-unresolved mismatch
# against the dry cooler's 40/46 C design point); the 4-10 C partial-
# free-cooling band is folded into full mechanical mode (same threshold
# as run_annual_simulation()).
# =======================================================================
def run_full_cycle(T_air_C, label=""):
    """Runs the full detailed cycle at one ambient dry-bulb temperature.
    Returns a flat dict covering every mass flow, pressure drop,
    temperature, and power figure needed for the monthly/seasonal report;
    fields that don't apply to the active mode are left as None."""
    row = {
        "label": label, "T_air_C": T_air_C, "mode": None,
        "m_dot_refrigerant_kgs": None, "n_active_compressors": None,
        "to_C": None, "tc_C": None, "T_discharge_C": None,
        "P_compressor_kW": None, "COP_compressor": None,
        "Q_condenser_kW": None, "condenser_subcool_K": None,
        "condenser_dP_ref_kPa": None, "condenser_dP_glycol_kPa": None,
        "eev_quality": None, "eev_Kv_m3h": None,
        "evaporator_superheat_K": None, "evaporator_dP_ref_kPa": None,
        "evaporator_dP_water_kPa": None,
        "econ_dP_water_kPa": None, "econ_dP_glycol_kPa": None,
        "econ_UA_kWK": None, "econ_LMTD_K": None,
        "tube_suction": None, "v_suction_ms": None,
        "tube_discharge": None, "v_discharge_ms": None,
        "tube_liquid": None, "v_liquid_ms": None,
        "m_dot_water_kgs": None, "m_dot_glycol_kgs": None,
        "P_glycol_pump_kW": None, "P_chw_pump_kW": None,
        "fan_speed_pct": None, "P_fan_kW": None,
        "Q_delivered_kW": None, "P_other_kW": None, "PUE": None,
        "COP_sys_fc": None,
    }

    if T_air_C <= FREE_COOLING_THRESHOLD_C:
        row["mode"] = "free_cooling"

        econ = Economizer(
            Q=Q_TARGET, T_water_in=21 + 273.15, T_water_out=15 + 273.15,
            T_glycol_in=9 + 273.15, T_glycol_out=15 + 273.15,
            D_h=D_h, L=ECON_L, beta=beta, Lambda=Lambda, b=b, L_w=ECON_L_W,
            N_cp_water=24, N_cp_glycol=24, U=3500.0,
        )

        fan_result = DRY_COOLER.predict_off_design(
            Q_target=Q_TARGET, T_glycol_hot_in=15 + 273.15, T_air_in=T_air_C + 273.15,
        )
        P_fan = DRY_COOLER.fan_power(fan_result["fan_speed_frac"], P_FAN_DESIGN)

        P_glycol_pump_row = glycol_pump.P_elec_free
        P_other_row = P_glycol_pump_row + P_fan + P_chw_pump  # compressor off
        PUE_row = (Q_TARGET + P_other_row) / Q_TARGET
        # COP_sys,fc = Q_cool / (W_fan + W_pump) -- compressor excluded since
        # it's off in free-cooling mode; only the dry-cooler fan and the two
        # circulation pumps (glycol + CHW) are actually doing work.
        COP_sys_fc = Q_TARGET / P_other_row

        row.update({
            "econ_dP_water_kPa": econ.delta_p_water / 1e3, "econ_dP_glycol_kPa": econ.delta_p_glycol / 1e3,
            "econ_UA_kWK": econ.UA / 1e3, "econ_LMTD_K": econ.LMTD,
            "m_dot_water_kgs": econ.m_dot_water, "m_dot_glycol_kgs": econ.m_dot_glycol,
            "P_compressor_kW": 0.0,
            "P_glycol_pump_kW": P_glycol_pump_row / 1e3, "P_chw_pump_kW": P_chw_pump / 1e3,
            "fan_speed_pct": fan_result["fan_speed_frac"] * 100, "P_fan_kW": P_fan / 1e3,
            "Q_delivered_kW": Q_TARGET / 1e3, "P_other_kW": P_other_row / 1e3, "PUE": PUE_row,
            "COP_sys_fc": COP_sys_fc,
        })
        return row

    # ---- MECHANICAL MODE ----
    row["mode"] = "mechanical"
    tc_C_p = T_air_C + APPROACH_K
    p_c_p = PropsSI("P", "T", tc_C_p + 273.15, "Q", 1, R)

    result = bank.solve_with_bypass(Q_TARGET, tc_C_p)
    T_o_p = result["to_C"] + 273.15
    p_o_p = PropsSI("P", "T", T_o_p, "Q", 1, R)
    m_dot_p = result["m_dot_total_kgh"] / 3600.0
    per_unit_p = result["per_unit"]

    h_2_p = per_unit_p["h_suction"] + per_unit_p["P_w"] / (per_unit_p["m_dot_kgh"] / 3600.0)
    T_2_p = PropsSI("T", "HMASS", h_2_p, "P", p_c_p, R)

    condenser_p = Condenser(m_dot_refrigerant=m_dot_p, p_in=p_c_p, refrigerant=R,
                             h_in=h_2_p, subcooling=SUBCOOLING,
                             T_glycol_in=30 + 273.15, T_glycol_out=35 + 273.15,
                             D_h=D_h, A_flow=A_flow_cond, L=L, beta=beta, Lambda=Lambda,
                             N_cp=N_cp_cond, b=b, L_w=L_w, N_cp_glycol=N_cp_glycol)

    eev_p = EEV(m_dot_refrigerant=m_dot_p, p_in=condenser_p.p_out, h_in=condenser_p.h_out,
                p_out=p_o_p, refrigerant=R)

    h_out_target_p = PropsSI('HMASS', 'T', T_o_p + DT_SH, 'P', p_o_p, R)
    Q_actual_p = m_dot_p * (h_out_target_p - eev_p.h_out)

    evaporator_p = Evaporator(m_dot_refrigerant=m_dot_p, p_in=p_o_p, Q=Q_actual_p, refrigerant=R,
                               T_water_in=21 + 273.15, T_water_out=15 + 273.15, h_in=eev_p.h_out,
                               D_h=D_h, A_flow=A_flow_evap, L=L, beta=beta, Lambda=Lambda,
                               N_cp=N_cp_evap, b=b, L_w=L_w, N_cp_water=N_cp_water)

    rho_suction_p = PropsSI('D', 'T', T_o_p + DT_SH, 'P', p_o_p, R)
    tube_s, ID_s, IDreq_s, v_s = select_tube(m_dot_p, rho_suction_p, v_min=4.5, v_max=20.0, v_target=15.0)
    rho_discharge_p = PropsSI('D', 'HMASS', h_2_p, 'P', p_c_p, R)
    tube_d, ID_d, IDreq_d, v_d = select_tube(m_dot_p, rho_discharge_p, v_min=10.0, v_max=18.0, v_target=14.0)
    rho_liquid_p = PropsSI('D', 'T', condenser_p.T_out, 'P', condenser_p.p_out, R)
    tube_l, ID_l, IDreq_l, v_l = select_tube(m_dot_p, rho_liquid_p, v_min=0.0, v_max=1.5, v_target=1.2)

    T_glycol_hot_in_p = (tc_C_p - GLYCOL_APPROACH_K) + 273.15
    fan_result = DRY_COOLER.predict_off_design(
        Q_target=condenser_p.Q, T_glycol_hot_in=T_glycol_hot_in_p, T_air_in=T_air_C + 273.15,
    )
    P_fan = DRY_COOLER.fan_power(fan_result["fan_speed_frac"], P_FAN_DESIGN)

    P_compressor_p = result["P_total_w"]
    P_glycol_pump_p = glycol_pump.P_elec_mech
    P_other_p = P_compressor_p + P_glycol_pump_p + P_fan + P_chw_pump
    Q_delivered_p = evaporator_p.Q
    PUE_p = (Q_delivered_p + P_other_p) / Q_delivered_p

    row.update({
        "m_dot_refrigerant_kgs": m_dot_p, "n_active_compressors": result["n_active"],
        "to_C": result["to_C"], "tc_C": tc_C_p, "T_discharge_C": T_2_p - 273.15,
        "P_compressor_kW": P_compressor_p / 1e3, "COP_compressor": result["COP"],
        "Q_condenser_kW": condenser_p.Q / 1e3, "condenser_subcool_K": condenser_p.subcool,
        "condenser_dP_ref_kPa": condenser_p.delta_p / 1e3,
        "condenser_dP_glycol_kPa": condenser_p.delta_p_glycol / 1e3,
        "eev_quality": eev_p.x_out, "eev_Kv_m3h": eev_p.Kv_required,
        "evaporator_superheat_K": evaporator_p.superheat, "evaporator_dP_ref_kPa": evaporator_p.delta_p / 1e3,
        "evaporator_dP_water_kPa": evaporator_p.delta_p_water / 1e3,
        "tube_suction": tube_s, "v_suction_ms": v_s,
        "tube_discharge": tube_d, "v_discharge_ms": v_d,
        "tube_liquid": tube_l, "v_liquid_ms": v_l,
        "m_dot_water_kgs": evaporator_p.m_dot_water, "m_dot_glycol_kgs": condenser_p.m_dot_glycol,
        "P_glycol_pump_kW": P_glycol_pump_p / 1e3, "P_chw_pump_kW": P_chw_pump / 1e3,
        "fan_speed_pct": fan_result["fan_speed_frac"] * 100, "P_fan_kW": P_fan / 1e3,
        "Q_delivered_kW": Q_delivered_p / 1e3, "P_other_kW": P_other_p / 1e3, "PUE": PUE_p,
    })
    return row


def compute_monthly_seasonal_averages():
    """Groups the real 8760-hour Champaign TMY3 dry-bulb data by calendar
    month (1-12) and by standard meteorological season (Winter=Dec/Jan/Feb,
    Spring=Mar/Apr/May, Summer=Jun/Jul/Aug, Fall=Sep/Oct/Nov). Returns two
    dicts of {label: average T_db [C]}."""
    wx = ChampaignWeather()
    month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    monthly_avg = {
        month_names[m - 1]: float(wx.T_db[wx.month == m].mean())
        for m in range(1, 13)
    }

    season_months = {
        "Winter": (12, 1, 2), "Spring": (3, 4, 5),
        "Summer": (6, 7, 8), "Fall": (9, 10, 11),
    }
    seasonal_avg = {
        season: float(wx.T_db[np.isin(wx.month, months)].mean())
        for season, months in season_months.items()
    }
    return monthly_avg, seasonal_avg


def run_monthly_seasonal_reports():
    """Runs run_full_cycle() at each month's and each season's average
    dry-bulb temperature and saves two CSV files -- monthly_report.csv
    (12 rows) and seasonal_report.csv (4 rows) -- covering every mass
    flow, pressure drop, temperature, and power figure from each run."""
    monthly_avg, seasonal_avg = compute_monthly_seasonal_averages()

    monthly_rows = [run_full_cycle(T_air_C, label=month) for month, T_air_C in monthly_avg.items()]
    seasonal_rows = [run_full_cycle(T_air_C, label=season) for season, T_air_C in seasonal_avg.items()]

    fieldnames = list(monthly_rows[0].keys())

    with open("monthly_report.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(monthly_rows)

    with open("seasonal_report.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(seasonal_rows)

    print(f"\nSaved monthly_report.csv ({len(monthly_rows)} rows) "
          f"and seasonal_report.csv ({len(seasonal_rows)} rows)")
    return monthly_rows, seasonal_rows


print()
print("=" * 72)
print("ANNUAL SIMULATION -- first pass (raw per-hour results, no aggregation yet)")
print("=" * 72)
annual_results = run_annual_simulation()
plot_annual_cop_load(annual_results)

print()
print("=" * 72)
print("MONTHLY / SEASONAL DETAILED REPORT")
print("=" * 72)
monthly_rows, seasonal_rows = run_monthly_seasonal_reports()
