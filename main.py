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
import os
import time

import numpy as np
from CoolProp.CoolProp import PropsSI

from compressor import CompressorBank
from condenser import Condenser
from dry_cooler import DryCooler
from dry_cooler_pump import GlycolLoopPump
from economizer import Economizer, rate_economizer
from eev import EEV
from evaporator import Evaporator
from water_side_pump import CHWPump, TRIM_CURVES_1_5AD_1750RPM
from weather import ChampaignWeather
from figsave import save_figure

# ---------------------------------------------------------------------
# Report output layout -- every figure/table this file produces lands
# directly in the matching report/ subfolder (mirrors the project spec's
# own deliverable structure: 1f seasonal performance, not an arbitrary
# "which script produced this" split), so nothing needs to be manually
# moved before it's \includegraphics{}'d / \input{}'d into the LaTeX report.
# ---------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
REPORT_DIR = os.path.join(_HERE, "report")
FIG_PERF_DESIGN = os.path.join(REPORT_DIR, "figures", "performance", "design_point")
FIG_PERF_MONTHLY = os.path.join(REPORT_DIR, "figures", "performance", "monthly")
FIG_PERF_SEASONAL = os.path.join(REPORT_DIR, "figures", "performance", "seasonal")
FIG_PERF_ANNUAL = os.path.join(REPORT_DIR, "figures", "performance", "annual")
FIG_ECON = os.path.join(REPORT_DIR, "figures", "economizer")
TABLE_PERF_MONTHLY = os.path.join(REPORT_DIR, "tables", "performance", "monthly")
TABLE_PERF_SEASONAL = os.path.join(REPORT_DIR, "tables", "performance", "seasonal")
for _d in (FIG_PERF_DESIGN, FIG_PERF_MONTHLY, FIG_PERF_SEASONAL, FIG_PERF_ANNUAL,
           FIG_ECON, TABLE_PERF_MONTHLY, TABLE_PERF_SEASONAL):
    os.makedirs(_d, exist_ok=True)

R = "R290"
# =======================================================================
# DESIGN REVISION (efficiency compliance -- see report Ch. Performance):
# the original design condensed at 52 C (17 K air-to-refrigerant approach:
# dry cooler 5 K + glycol rise 6 K + condenser 6 K) and evaporated 8 K
# below the leaving chilled water. That design missed both efficiency
# targets (AHRI full-load COP 2.64 < 3.5, IPLV.IP 4.43 < 5.0). Revised:
#   * condensing approach 17 -> 10 K (dry cooler 3 K + glycol rise 4 K +
#     condenser 3 K), i.e. tc = 45 C at the 35 C design ambient -- larger
#     dry cooler + condenser area, same hardware families;
#   * evaporator approach 8 -> 5 K, i.e. to = 10 C at 15 C LCHW -- more
#     evaporator plates;
#   * a VFD on the compressors: at full load the bank over-delivers at the
#     evaporator-dictated to = 10 C, so the drive trims speed to deliver
#     exactly 150 kW there instead of letting to float down (fixed-speed
#     behaviour), and below the old turndown floor it modulates speed
#     instead of hot-gas bypassing (compressor.solve_with_vfd()).
# Result: AHRI full-load COP 3.53 >= 3.5, IPLV.IP 5.23 >= 5.0.
# =======================================================================
T_c = 45 + 273.15         # condensing temperature [K] at the 35 C design ambient
                          # (10 K total air-to-refrigerant approach, design revision)
DT_SH = 8                 # suction superheat [K] -- the compressor vendor's own selection
                          # basis (Ch. Compressor: "useful superheat 8 K, 15 C return gas").
                          # Was 10 K, an undocumented deviation that made this file's
                          # design-point COP read 4.73 against performance.py's 4.67;
                          # aligned to 8 K so the report quotes one design-point COP.
SUBCOOLING = 5            # condenser design subcooling [K] (matches eev.py)
Q_TARGET = 150e3          # design cooling duty [W]
N_UNITS = 2               # project's dual 75 kW parallel-circuit design (compressor.py)
TO_SETPOINT_C = 10.0      # evaporating temp the evaporator's installed area sustains at
                          # full load: LCHW 15 C - 5 K design approach (design revision)

# Annual-simulation constants (used by run_annual_simulation() below)
T_AIR_DESIGN_C = 35.0
APPROACH_K = (T_c - 273.15) - T_AIR_DESIGN_C  # 10 K, fixed condensing approach (same
                                               # simplification already used in compressor.py's
                                               # IPLV.IP calc) -- holds tc-T_air constant across
                                               # all ambient conditions rather than re-solving the
                                               # dry cooler's eps-NTU model at each point
GLYCOL_APPROACH_K = 3.0    # tc -> condenser glycol-leaving-temp approach (45->42 C at design)
GLYCOL_RANGE_K = 4.0       # condenser glycol supply/return range (38->42 C at design,
                           # matches DRY_COOLER's revised 38/42 C design point below)
CP_GLYCOL_NOMINAL = 3835.0 # J/kgK, 30% PG at ~40 C -- for the mechanical-mode loop
                           # glycol flow = Q_cond / (cp x range); used only to hand the
                           # dry cooler the real loop flow (not its rating-point flow)
FREE_COOLING_THRESHOLD_C = 4.0  # ASHRAE 90.1 full free-cooling threshold (economizer.py)
ECON_L, ECON_L_W = 0.70, 0.25    # economizer's own larger frame (economizer.py's design point)

# ASHRAE 90.1 s6.5.1 economizer ACTIVATION threshold. Between this and the full
# free-cooling threshold above lies the INTEGRATED (partial) free-cooling band:
# the economizer pre-cools the CHW return by whatever the ambient allows and the
# chiller trims the remainder (Ch. Economizer, "Full vs. Partial Economizer
# Operation"). 1,358 hrs/yr in Champaign -- 15.5% of the year.
PARTIAL_FC_THRESHOLD_C = 10.0
# Dry-cooler approach on the free-cooling glycol circuit: OAT + 5 K is the
# economizer chapter's own stated basis (4 C OAT -> 9 C glycol at the sizing
# point, "the dry cooler delivers glycol between ~9 and 15 C" across the band).
DRY_COOLER_APPROACH_FC_K = 5.0

# Lowest speed the drive can actually hold. Below it the compressor CYCLES
# on/off against the min-speed operating point rather than modulating.
#
# This matters at the cold end of the partial band: with the economizer
# already carrying ~96% of the load at 4.5 C, the chiller is asked for ~6 kW,
# which is ~4% speed -- no drive does that. The power reported there is still
# right, though: ideal cycling runs at min speed for a duty fraction
# Q_req/Q_at_min, so average power = Q_req/COP -- algebraically identical to
# what the modulation model returns. What is NOT modelled is real cycling
# LOSS (start-up transients, off-cycle migration). That is bounded by the
# few kW of compressor duty involved in this band, so it cannot move PUE
# materially, but the mechanism should be reported as cycling, not modulation.
VFD_MIN_SPEED_FRAC = 0.30

p_c = PropsSI("P", "T", T_c, "Q", 1, R)

# Shared dry cooler instance -- RE-SELECTED for the revised design point.
#
# The original Kelvion ULF-PA104Y4V-096Z100 (727.56 MBH at 46/40 C glycol,
# 35 C air) cannot serve the revised 42/38 C glycol condition: its rating
# is at an 11 K glycol-to-air driving dT, while 42/38 C gives only 7 K, so
# the same coil delivers ~146 kW against the ~182 kW required (it needed
# ~150% fan speed -- the source of the 17 fan-ceiling hours flagged in the
# annual run). Capacity is a number AT A CONDITION, not a property of the box.
#
# Re-selected via Kelvion Select RT at the revised duty:
#   KELVION ULF-PA106K4V-091F095
#     670.36 MBH = 196.5 kW  @ glycol 42 -> 38 C, air 35 C  (8% margin on the
#                                                            ~182 kW required)
#     surface 20,977 ft2 = 1,949 m2  (2.1x the original's 938 m2 -- physically
#          a much larger coil, despite the LOWER headline MBH, which is quoted
#          at the harder 3 K-approach condition)
#     air 93,782 cfm, outlet 101.7 F  -> dT_air = 3.72 K
#     6 x 0.910 m EC fans @ 950/1050 rpm (10% speed headroom, unlike the
#          "Z110" variants which run 1100/1100 flat out with no reserve)
#     60 dB(A) at 33 ft -- matches the original unit's sound level, the
#          criterion the original selection was made on (see Ch. Dry Cooler)
#     glycol-side dP 8 psi = 55.16 kPa (up from 7 psi; drives the glycol pump
#          frame re-check, already an open item)
DRY_COOLER = DryCooler(
    Q_design=196.5e3, T_glycol_hot_in=42 + 273.15, T_glycol_cold_out=38 + 273.15,
    T_air_in_design=T_AIR_DESIGN_C + 273.15, dT_air_design=3.72,
)

# ---------------------------------------------------------------------
# 1) Compressor bank -- VFD operating model (design revision): the
#    evaporator's installed area fixes to = TO_SETPOINT_C (10 C) at full
#    load, and since the two-unit bank at full speed over-delivers there
#    (~201 kW at to=10/tc=45), the VFD trims speed (~75%) so the bank
#    delivers exactly Q_TARGET at that setpoint -- rather than letting
#    the evaporating temperature float down as a fixed-speed bank must
#    (which is worth ~1 COP point at this design point). Staging still
#    picks the fewest units that can reach the load at the setpoint.
# ---------------------------------------------------------------------
bank = CompressorBank(n_units=N_UNITS, superheat_K=DT_SH, subcooling_K=SUBCOOLING,
                       to_bounds_C=(-10.0, 15.0))
bank_result = bank.solve_with_vfd(Q_TARGET, T_c - 273.15, to_setpoint_C=TO_SETPOINT_C)

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
# 2) Condenser -- same corrugation family as the evaporator (D_h, Lambda,
#    beta), sharing the Kelvion HP DW 500H plate.
#
#    PLATE DIMENSIONS ARE NOW VENDOR-CONFIRMED (superseding the earlier
#    L=0.5/L_w=0.2 placeholders): the HP DW 500H plate is 195 x 600 mm
#    (0.117 m2 of heat-transfer area, 120-plate frame). L_w and L below are
#    those real dimensions. NOTE: 600 mm is the plate's stated heat-transfer
#    length; the true port-to-port flow length L_v is marginally shorter, so
#    using it here is slightly conservative for pressure drop (dP ~ L_v).
#    D_h/Lambda/beta remain the corrugation-family placeholders -- Kelvion
#    does not publish them per model.
# ---------------------------------------------------------------------
D_h, Lambda, beta, L = 0.004, 0.005, 30, 0.6
b, L_w = D_h * 1.17 / 2, 0.195
# N_cp=47 (not 38, the evaporator's value) -- the condenser carries a higher
# refrigerant mass flow than the evaporator, so it needs more channels to
# hold the same G_c~25 kg/m2s design target inside Han's validated 13-34
# kg/m2s range (Project_Report_full.txt Sec 3.11: N_cp ~= 0.553/(25*4.68e-4) ~= 47).
N_cp_cond = 47
A_flow_cond = N_cp_cond * b * L_w
N_cp_glycol = 24

# Glycol boundary condition floats with the condensing temperature (fixed
# 6 K approach to tc, fixed 6 K supply/return range) instead of a stale
# hardcoded value -- this now matches DRY_COOLER's own 40/46 C design
# point exactly (52 C tc -> 46/40 C glycol) instead of contradicting it.
T_glycol_out = (T_c - 273.15 - GLYCOL_APPROACH_K) + 273.15
T_glycol_in = T_glycol_out - GLYCOL_RANGE_K

condenser = Condenser(m_dot_refrigerant=m_dot, p_in=p_c, refrigerant=R,
                       h_in=h_2, subcooling=SUBCOOLING,
                       T_glycol_in=T_glycol_in, T_glycol_out=T_glycol_out,
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
    dP_drycooler_mech_ref=55.16e3, dP_econ_free_ref=10.43e3,   # 8 psi, re-selected unit
    D_h=D_h, L=L, beta=beta, b=b, L_w=L_w, N_cp_glycol=N_cp_glycol,
    D_pipe=0.0779, pipe_roughness=0.045e-3, L_eq_total=65.0,
)

chw_pump = CHWPump(
    m_dot_water=evaporator.m_dot_water, dp_hx_kpa=evaporator.delta_p_water / 1e3,
    dp_crah_kpa=50.0, trim_curves_gpm_ft=TRIM_CURVES_1_5AD_1750RPM,
)
chw_pump.select_trim()

# Dry cooler RATED fan power, for the RE-SELECTED Kelvion ULF-PA106K4V-091F095
# (6 x 0.910 m EC fans) at its 196.5 kW / 100%-air-flow selection point.
# ESTIMATE, not vendor-confirmed: the Select RT results table does not report fan
# power, so this is scaled from a same-family confirmed 2.685 kW/fan
# (0.960 m @ 1000 rpm) by the fan laws, P ~ N^3 * D^5:
#     6 x 2.685 * (950/1000)^3 * (0.910/0.960)^5 = 10.6 kW
# CONFIRM from the datasheet -- this is the single largest uncertainty in the PUE.
#
# This is the power at the coil's full SELECTION duty (196.5 kW). The plant only
# rejects ~182 kW at the 35 C design day, and the unit has 8% capacity margin, so
# the fans do NOT run at full design flow there -- the off-design solve below finds
# the actual speed (86.8%) and the cube law gives the actual design-day fan power
# (~6.9 kW). Using the flat 10.6 kW would charge the PUE for capacity the plant
# never uses.
P_FAN_RATED = 10.6e3
P_FAN_DESIGN = P_FAN_RATED   # kept for the annual sim's fan_power() reference

# Design-day fan power: solve the fan speed that rejects the actual 182 kW condenser
# duty at 35 C, then apply the cube law -- the same machinery the annual sim uses,
# so the design-day PUE is now consistent with the monthly/seasonal numbers.
_fan_design = DRY_COOLER.predict_off_design(
    Q_target=condenser.Q, T_glycol_hot_in=condenser.T_glycol_out, T_air_in=T_AIR_DESIGN_C + 273.15,
    m_dot_glycol=condenser.m_dot_glycol)   # the actual loop flow (11.64), not the unit rating flow
P_fan_designday = DRY_COOLER.fan_power(_fan_design["fan_speed_frac"], P_FAN_RATED)

P_compressor = bank_result["P_total_w"]
P_glycol_pump = glycol_pump.P_elec_mech
# eta_pump=0.60 matches water_side_pump.py's own contour-line reading; eta_motor=0.90
# matches the glycol pump's motor-efficiency assumption (no CHW motor-efficiency
# stage built into CHWPump itself, so applied here explicitly)
P_chw_pump = chw_pump.shaft_power_w(eta_pump=0.60) / 0.90
P_other = P_compressor + P_glycol_pump + P_fan_designday + P_chw_pump
PUE = (Q_TARGET + P_other) / Q_TARGET

# TUE (Total-power Usage Effectiveness, Patterson et al. 2013):
#     TUE = ITUE x PUE,   ITUE = (total IT power) / (compute power)
# PUE stops at the IT rack and is blind to losses INSIDE the IT equipment
# (server PSU, VRMs, on-board fans), so two facilities with identical PUE
# can differ in real compute delivered per watt. TUE closes that gap by
# referencing every facility watt to actual compute power.
#
# ASSUMPTION, not derived: this project models the 150 kW as a black-box IT
# load -- there is no server-level power model anywhere in it, and the
# project brief specifies none. ITUE = 1.20 is a mid-range literature value
# for modern volume servers (~94% PSU efficiency, ~5% on-board fans, plus
# VRM/conversion losses). Reported band: 1.10 (best-in-class 80PLUS
# Titanium, minimal fan power) to 1.30 (older or fan-heavy hardware).
#
# Consequence worth being honest about: with ITUE held constant, TUE is
# just PUE scaled by 1.20. It re-ranks nothing in this report and adds no
# new physics -- its value is that it states the IT-side overhead the PUE
# convention hides, and it shows how much of the total facility power the
# cooling plant is actually responsible for. A real TUE study would need
# measured server data.
ITUE = 1.20
ITUE_BAND = (1.10, 1.30)
TUE = ITUE * PUE

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
print("=== DESIGN-DAY PUE / TUE ===")
print(f"IT load (proxy)      : {Q_TARGET/1e3:.1f} kW")
print(f"Compressor bank      : {P_compressor/1e3:.2f} kW")
print(f"Glycol loop pump     : {P_glycol_pump/1e3:.2f} kW")
print(f"Dry cooler fans      : {P_fan_designday/1e3:.2f} kW  (off-design solve at 35 C, "
      f"{_fan_design['fan_speed_frac']*100:.0f}% of {P_FAN_RATED/1e3:.1f} kW rated)")
print(f"CHW/CRAH pump        : {P_chw_pump/1e3:.2f} kW")
print(f"Other facility loads : {P_other/1e3:.2f} kW")
print(f"PUE = (IT + other)/IT = {PUE:.3f}")
print(f"TUE = ITUE x PUE      = {TUE:.3f}   (ITUE = {ITUE:.2f} ASSUMED, band "
      f"{ITUE_BAND[0]:.2f}-{ITUE_BAND[1]:.2f} -> TUE {ITUE_BAND[0]*PUE:.3f}-{ITUE_BAND[1]*PUE:.3f})")


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
    10 K approach (same simplification already used in compressor.py's
    IPLV.IP calculation), then runs the VFD operating model at the
    evaporator-dictated to = 10 C setpoint (design revision). Skips
    re-instantiating the full condenser/EEV/evaporator chain -- not
    needed for an annual energy sum, only for the single design-point
    validation already done above."""
    tc_C = T_air_C + APPROACH_K
    result = bank_obj.solve_with_vfd(Q_target, tc_C, to_setpoint_C=TO_SETPOINT_C)
    return {
        "mode": "mechanical", "T_air_C": T_air_C, "tc_C": tc_C,
        "P_compressor_w": result["P_total_w"], "Q_delivered_w": result["Q_total_w"],
        "Q_econ_w": 0.0, "Q_mech_w": result["Q_total_w"],
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
        "Q_econ_w": Q_target, "Q_mech_w": 0.0,
    }


# The economizer AS INSTALLED, built once at its sizing point. Only its fixed
# hardware UA and the two stream heat-capacity rates are reused below -- the
# partial-band solve rates this same exchanger off-design, it does not re-size
# it. Balanced by design (6 K glycol rise mirrors the 6 K CHW span), so
# C_water ~= C_glycol and Cr ~= 1.
ECON_DESIGN = Economizer(
    Q=Q_TARGET, T_water_in=21 + 273.15, T_water_out=15 + 273.15,
    T_glycol_in=9 + 273.15, T_glycol_out=15 + 273.15,
    D_h=D_h, L=ECON_L, beta=beta, Lambda=Lambda, b=b, L_w=ECON_L_W,
    N_cp_water=24, N_cp_glycol=24, U=3500.0,
)


def econ_duty_at(T_air_C):
    """Duty the INSTALLED economizer delivers against ambient, by rating its
    fixed UA off-design (epsilon-NTU). Glycol arrives at OAT + 5 K, the
    economizer chapter's own basis; the CHW return is the fixed 21 C.

    Sanity anchor: at OAT = 4 C this returns exactly the 150 kW the LMTD
    sizing model was built for, so the rating and sizing models agree at
    their shared point and the partial band joins the full-free-cooling
    band continuously -- no discontinuity at the 4 C threshold.
    """
    T_glycol_in = (T_air_C + DRY_COOLER_APPROACH_FC_K) + 273.15
    Q, eps, NTU, Cr = rate_economizer(
        ECON_DESIGN.UA, ECON_DESIGN.m_dot_water, ECON_DESIGN.cp_water, 21 + 273.15,
        ECON_DESIGN.m_dot_glycol, ECON_DESIGN.cp_glycol, T_glycol_in,
    )
    return min(Q, Q_TARGET), T_glycol_in - 273.15


def partial_free_cooling_hour(bank_obj, T_air_C, Q_target=Q_TARGET):
    """Per-hour solve for an INTEGRATED (partial) free-cooling hour,
    4 < OAT <= 10 C -- the ASHRAE 90.1 s6.5.1 activation band.

    Water-side series, exactly as Ch. Economizer specifies: the CHW return
    (21 C) passes through the economizer, which pre-cools it as far as the
    ambient allows, and the chiller trims what's left down to the 15 C
    supply. So the compressor sees only Q_mech = 150 kW - Q_econ, at the
    same to = 10 C setpoint and the same tc = OAT + 10 K basis as a
    mechanical hour.

    The chiller's lift is genuinely low here (OAT <= 10 C -> tc <= 20 C) and
    its duty is small, so these hours are cheap -- which is the entire point
    of the 90.1 activation requirement.
    """
    Q_econ, T_glycol_C = econ_duty_at(T_air_C)
    Q_mech = Q_target - Q_econ

    if Q_mech <= 0.0:      # economizer alone covers it -> compressor off
        r = free_cooling_hour(T_air_C, Q_target)
        r["mode"] = "partial_free_cooling"
        r["T_glycol_C"] = T_glycol_C
        return r

    tc_C = T_air_C + APPROACH_K
    result = bank_obj.solve_with_vfd(Q_mech, tc_C, to_setpoint_C=TO_SETPOINT_C)
    speed_frac = result.get("speed_frac")
    return {
        "mode": "partial_free_cooling", "T_air_C": T_air_C, "tc_C": tc_C,
        "P_compressor_w": result["P_total_w"],
        "Q_delivered_w": Q_target,          # the hall still gets its full 150 kW
        "Q_econ_w": Q_econ, "Q_mech_w": Q_mech,
        "T_glycol_C": T_glycol_C,
        "speed_frac": speed_frac,
        "cycling": speed_frac is not None and speed_frac < VFD_MIN_SPEED_FRAC,
        "status": result["status"],
    }


def run_annual_simulation():
    """Runs mechanical_hour()/free_cooling_hour() plus the dry cooler's
    off-design fan solve across every hour of the real embedded Champaign
    TMY3 weather file. Returns the list of per-hour result dicts (not yet
    aggregated -- that's the next step)."""
    wx = ChampaignWeather()
    bank_obj = CompressorBank(n_units=N_UNITS, superheat_K=DT_SH, subcooling_K=SUBCOOLING,
                               to_bounds_C=(-10.0, 15.0))
    P_FAN_RATED = P_FAN_DESIGN  # re-selected Kelvion ULF-PA106K4V-091F095 (see above)

    results = []
    t0 = time.time()
    C_glycol_free = ECON_DESIGN.m_dot_glycol * ECON_DESIGN.cp_glycol   # [W/K]

    for hr in wx:
        T_air_C = hr.T_db
        if T_air_C <= FREE_COOLING_THRESHOLD_C:
            r = free_cooling_hour(T_air_C)
            T_glycol_hot_in_hour = 15 + 273.15
            m_dot_glycol_hour = None   # free cooling: economizer-loop flow, open item
        elif T_air_C <= PARTIAL_FC_THRESHOLD_C:
            r = partial_free_cooling_hour(bank_obj, T_air_C)
            # Common glycol loop at its free-cooling flow: it leaves the dry
            # cooler cold at OAT + 5 K, picks up BOTH the economizer's duty and
            # the (small) condenser rejection, and returns mixed. Closing the
            # loop on energy gives the dry cooler's hot inlet directly.
            Q_reject_hour = r["Q_delivered_w"] + r["P_compressor_w"]
            T_glycol_hot_in_hour = ((r["T_glycol_C"] + Q_reject_hour / C_glycol_free)
                                    + 273.15)
            m_dot_glycol_hour = None   # free/partial: economizer-loop flow, open item; use default
        else:
            r = mechanical_hour(bank_obj, T_air_C)
            T_glycol_hot_in_hour = (r["tc_C"] - GLYCOL_APPROACH_K) + 273.15
            # Mechanical loop flow = condenser duty / (cp x 4 K range), the real flow
            # the loop carries (not the coil's rating-point flow).
            m_dot_glycol_hour = (r["Q_delivered_w"] + r["P_compressor_w"]) / (
                CP_GLYCOL_NOMINAL * GLYCOL_RANGE_K)

        # The dry cooler rejects everything the plant absorbs -- the full IT load
        # plus whatever compressor work was needed. True in all three regimes,
        # since Q_delivered is the 150 kW hall duty however it was produced.
        Q_cond_hour = r["Q_delivered_w"] + r["P_compressor_w"]
        fan_result = DRY_COOLER.predict_off_design(
            Q_target=Q_cond_hour, T_glycol_hot_in=T_glycol_hot_in_hour, T_air_in=T_air_C + 273.15,
            m_dot_glycol=m_dot_glycol_hour,
        )
        r["fan_speed_frac"] = fan_result["fan_speed_frac"]
        r["P_fan_w"] = DRY_COOLER.fan_power(fan_result["fan_speed_frac"], P_FAN_RATED)
        r["hour_of_year"] = hr.hour_of_year
        results.append(r)

    elapsed = time.time() - t0
    n = len(results)
    n_mech = sum(1 for r in results if r["mode"] == "mechanical")
    n_part = sum(1 for r in results if r["mode"] == "partial_free_cooling")
    n_free = sum(1 for r in results if r["mode"] == "free_cooling")
    econ_kwh = sum(r.get("Q_econ_w", 0.0) for r in results) / 1e3
    load_kwh = sum(r["Q_delivered_w"] for r in results) / 1e3
    print(f"\nAnnual simulation: {n} hours processed in {elapsed:.1f} s")
    print(f"  Mechanical only (OAT > {PARTIAL_FC_THRESHOLD_C:.0f} C)                 : "
          f"{n_mech} ({n_mech/n*100:.1f}%)")
    print(f"  Integrated/partial free cooling ({FREE_COOLING_THRESHOLD_C:.0f} < OAT <= "
          f"{PARTIAL_FC_THRESHOLD_C:.0f} C) : {n_part} ({n_part/n*100:.1f}%)")
    print(f"  Full free cooling (OAT <= {FREE_COOLING_THRESHOLD_C:.0f} C)               : "
          f"{n_free} ({n_free/n*100:.1f}%)")
    print(f"  Cooling energy delivered by the economizer : {econ_kwh/1e3:.0f} MWh/yr of "
          f"{load_kwh/1e3:.0f} MWh/yr ({econ_kwh/load_kwh*100:.1f}%)")
    n_cyc = sum(1 for r in results if r.get("cycling"))
    if n_cyc:
        print(f"  NOTE: {n_cyc} partial-band hour(s) ({n_cyc/n*100:.1f}%) need less duty than "
              f"{VFD_MIN_SPEED_FRAC*100:.0f}% drive speed -> chiller CYCLES rather than "
              f"modulates.\n        Average power is unaffected (ideal cycling = modulation); "
              f"cycling losses are not modelled.")
    n_over_speed = sum(1 for r in results if r["fan_speed_frac"] > 1.0)
    if n_over_speed:
        print(f"  WARNING: {n_over_speed} hour(s) require >100% fan speed to hit duty -- "
              f"dry cooler capacity-ceiling hours, flagged for follow-up.")
    return results


def plot_annual_cop_load(results, save_path=None):
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

    if save_path is None:
        save_path = os.path.join(FIG_PERF_ANNUAL, "annual_cop_load.png")

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
    save_figure(fig, save_path, facecolor=SURFACE)
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
# All three ASHRAE 90.1 economizer regimes are modelled: mechanical
# (OAT > 10 C), integrated/partial free cooling (4 < OAT <= 10 C, where the
# economizer pre-cools the CHW return and the chiller trims the remainder),
# and full free cooling (OAT <= 4 C, compressor off) -- same thresholds as
# run_annual_simulation().
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
        # Integrated/partial free-cooling split (4 < OAT <= 10 C); None in the
        # other two regimes, where the load is carried entirely by one side.
        "Q_econ_kW": None, "Q_mech_kW": None,
        "econ_glycol_supply_C": None, "CHW_econ_outlet_C": None,
        "tube_suction": None, "v_suction_ms": None,
        "tube_discharge": None, "v_discharge_ms": None,
        "tube_liquid": None, "v_liquid_ms": None,
        "m_dot_water_kgs": None, "m_dot_glycol_kgs": None,
        "P_glycol_pump_kW": None, "P_chw_pump_kW": None,
        "fan_speed_pct": None, "P_fan_kW": None,
        "Q_delivered_kW": None, "P_other_kW": None, "PUE": None, "TUE": None,
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
            "TUE": ITUE * PUE_row,
            "COP_sys_fc": COP_sys_fc,
        })
        return row

    # ---- MECHANICAL / INTEGRATED-PARTIAL-FREE-COOLING MODE ----
    # Both regimes run the same chiller chain below; the partial band differs
    # only in that the economizer has already removed part of the load before
    # the water reaches the evaporator, so the compressor is asked for less
    # duty against a lower entering-water temperature. Set that up first.
    Q_econ_row = 0.0
    T_water_in_C = 21.0
    if T_air_C <= PARTIAL_FC_THRESHOLD_C:
        row["mode"] = "partial_free_cooling"
        Q_econ_row, T_glycol_C_row = econ_duty_at(T_air_C)
        # Water-side series: the economizer pre-cools the 21 C return, the
        # chiller trims the rest to the 15 C supply (Ch. Economizer).
        T_water_in_C = 21.0 - Q_econ_row / (ECON_DESIGN.m_dot_water * ECON_DESIGN.cp_water)
        row.update({
            "econ_UA_kWK": ECON_DESIGN.UA / 1e3,
            "econ_dP_water_kPa": ECON_DESIGN.delta_p_water / 1e3,
            "econ_dP_glycol_kPa": ECON_DESIGN.delta_p_glycol / 1e3,
            "Q_econ_kW": Q_econ_row / 1e3,
            "econ_glycol_supply_C": T_glycol_C_row,
            "CHW_econ_outlet_C": T_water_in_C,
        })
    else:
        row["mode"] = "mechanical"

    Q_duty = Q_TARGET - Q_econ_row     # what the compressor is actually asked for
    row["Q_mech_kW"] = Q_duty / 1e3

    tc_C_p = T_air_C + APPROACH_K
    p_c_p = PropsSI("P", "T", tc_C_p + 273.15, "Q", 1, R)

    result = bank.solve_with_vfd(Q_duty, tc_C_p, to_setpoint_C=TO_SETPOINT_C)
    T_o_p = result["to_C"] + 273.15
    p_o_p = PropsSI("P", "T", T_o_p, "Q", 1, R)
    m_dot_p = result["m_dot_total_kgh"] / 3600.0
    per_unit_p = result["per_unit"]

    h_2_p = per_unit_p["h_suction"] + per_unit_p["P_w"] / (per_unit_p["m_dot_kgh"] / 3600.0)
    T_2_p = PropsSI("T", "HMASS", h_2_p, "P", p_c_p, R)

    # Glycol boundary condition floats with tc_C_p (same fixed 3 K approach/
    # 4 K range as the design-point section above) instead of the stale
    # fixed 30/35 C that could fall below tc_C_p at low ambient -- see
    # T_glycol_hot_in_p below, which is now just this same T_glycol_out_p.
    T_glycol_out_p = (tc_C_p - GLYCOL_APPROACH_K) + 273.15
    T_glycol_in_p = T_glycol_out_p - GLYCOL_RANGE_K

    condenser_p = Condenser(m_dot_refrigerant=m_dot_p, p_in=p_c_p, refrigerant=R,
                             h_in=h_2_p, subcooling=SUBCOOLING,
                             T_glycol_in=T_glycol_in_p, T_glycol_out=T_glycol_out_p,
                             D_h=D_h, A_flow=A_flow_cond, L=L, beta=beta, Lambda=Lambda,
                             N_cp=N_cp_cond, b=b, L_w=L_w, N_cp_glycol=N_cp_glycol)

    eev_p = EEV(m_dot_refrigerant=m_dot_p, p_in=condenser_p.p_out, h_in=condenser_p.h_out,
                p_out=p_o_p, refrigerant=R)

    h_out_target_p = PropsSI('HMASS', 'T', T_o_p + DT_SH, 'P', p_o_p, R)
    Q_actual_p = m_dot_p * (h_out_target_p - eev_p.h_out)

    evaporator_p = Evaporator(m_dot_refrigerant=m_dot_p, p_in=p_o_p, Q=Q_actual_p, refrigerant=R,
                               T_water_in=T_water_in_C + 273.15, T_water_out=15 + 273.15,
                               h_in=eev_p.h_out,
                               D_h=D_h, A_flow=A_flow_evap, L=L, beta=beta, Lambda=Lambda,
                               N_cp=N_cp_evap, b=b, L_w=L_w, N_cp_water=N_cp_water)

    rho_suction_p = PropsSI('D', 'T', T_o_p + DT_SH, 'P', p_o_p, R)
    tube_s, ID_s, IDreq_s, v_s = select_tube(m_dot_p, rho_suction_p, v_min=4.5, v_max=20.0, v_target=15.0)
    rho_discharge_p = PropsSI('D', 'HMASS', h_2_p, 'P', p_c_p, R)
    tube_d, ID_d, IDreq_d, v_d = select_tube(m_dot_p, rho_discharge_p, v_min=10.0, v_max=18.0, v_target=14.0)
    rho_liquid_p = PropsSI('D', 'T', condenser_p.T_out, 'P', condenser_p.p_out, R)
    tube_l, ID_l, IDreq_l, v_l = select_tube(m_dot_p, rho_liquid_p, v_min=0.0, v_max=1.5, v_target=1.2)

    # The dry cooler rejects the condenser duty AND, in the partial band, the
    # economizer's duty as well -- one loop, one coil, both heat sources. In
    # that band the loop runs at its free-cooling flow and leaves the coil at
    # OAT + 5 K, so closing the loop on energy gives the mixed return.
    Q_reject_p = condenser_p.Q + Q_econ_row
    if Q_econ_row > 0.0:
        C_glycol_free = ECON_DESIGN.m_dot_glycol * ECON_DESIGN.cp_glycol
        T_glycol_hot_p = (T_glycol_C_row + Q_reject_p / C_glycol_free) + 273.15
        P_glycol_pump_p = glycol_pump.P_elec_free
    else:
        T_glycol_hot_p = T_glycol_out_p
        P_glycol_pump_p = glycol_pump.P_elec_mech

    # Pure mechanical mode: the loop carries the condenser's own glycol flow.
    # Partial band: the glycol serves both economizer and condenser at the
    # free-cooling flow (open item), so fall back to the model default there.
    m_dot_glycol_p = condenser_p.m_dot_glycol if Q_econ_row == 0 else None
    fan_result = DRY_COOLER.predict_off_design(
        Q_target=Q_reject_p, T_glycol_hot_in=T_glycol_hot_p, T_air_in=T_air_C + 273.15,
        m_dot_glycol=m_dot_glycol_p,
    )
    P_fan = DRY_COOLER.fan_power(fan_result["fan_speed_frac"], P_FAN_DESIGN)

    P_compressor_p = result["P_total_w"]
    P_other_p = P_compressor_p + P_glycol_pump_p + P_fan + P_chw_pump
    # What the data hall actually receives: the economizer's share plus the
    # chiller's. In the mechanical band Q_econ_row is 0 and this is unchanged.
    Q_delivered_p = evaporator_p.Q + Q_econ_row
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
        "TUE": ITUE * PUE_p,
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

    monthly_csv_path = os.path.join(TABLE_PERF_MONTHLY, "monthly_report.csv")
    seasonal_csv_path = os.path.join(TABLE_PERF_SEASONAL, "seasonal_report.csv")

    with open(monthly_csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(monthly_rows)

    with open(seasonal_csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(seasonal_rows)

    print(f"\nSaved {monthly_csv_path} ({len(monthly_rows)} rows) "
          f"and {seasonal_csv_path} ({len(seasonal_rows)} rows)")
    return monthly_rows, seasonal_rows


def _row_cop_total(row):
    """Unified system COP = Q_delivered / P_other, valid in either mode:
    algebraically identical to 1/(PUE-1) since PUE = (Q+P_other)/Q, so this
    single formula reproduces both COP_compressor-derived (mechanical) and
    COP_sys_fc (free-cooling) without needing an if/else on mode."""
    return 1.0 / (row["PUE"] - 1.0)


def plot_economizer_band(save_path=None):
    """How the 150 kW load splits between the economizer and the chiller across
    the whole economizer envelope, and what the compressor costs there.

    This is the ASHRAE 90.1 s6.5.1 compliance evidence in one picture: the
    economizer carries load continuously from the 10 C activation threshold
    down to 100% at the 4 C full-free-cooling threshold, with no step at
    either boundary.
    """
    import matplotlib.pyplot as plt

    if save_path is None:
        save_path = os.path.join(FIG_ECON, "economizer_band.png")

    BLUE, AQUA, YELLOW, RED = "#2a78d6", "#1baf7a", "#eda100", "#e34948"
    SURFACE, PRIMARY_INK, SECONDARY_INK = "#fcfcfb", "#0b0b0b", "#52514e"
    MUTED, GRID, BASELINE = "#898781", "#e1e0d9", "#c3c2b7"

    bank_obj = CompressorBank(n_units=N_UNITS, superheat_K=DT_SH, subcooling_K=SUBCOOLING,
                               to_bounds_C=(-10.0, 15.0))
    oat = np.linspace(0.0, 14.0, 141)
    q_econ, q_mech, p_comp = [], [], []
    for t in oat:
        if t <= FREE_COOLING_THRESHOLD_C:
            qe, qm, pc = Q_TARGET, 0.0, 0.0
        elif t <= PARTIAL_FC_THRESHOLD_C:
            r = partial_free_cooling_hour(bank_obj, float(t))
            qe, qm, pc = r["Q_econ_w"], r["Q_mech_w"], r["P_compressor_w"]
        else:
            r = mechanical_hour(bank_obj, float(t))
            qe, qm, pc = 0.0, r["Q_delivered_w"], r["P_compressor_w"]
        q_econ.append(qe / 1e3); q_mech.append(qm / 1e3); p_comp.append(pc / 1e3)

    fig, (ax_q, ax_p) = plt.subplots(2, 1, figsize=(9, 6.4), facecolor=SURFACE,
                                     sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    for ax in (ax_q, ax_p):
        ax.set_facecolor(SURFACE)
        ax.grid(True, axis="y", color=GRID, linewidth=0.8, zorder=0)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(BASELINE)
        ax.tick_params(colors=MUTED, labelsize=8)
        for thr in (FREE_COOLING_THRESHOLD_C, PARTIAL_FC_THRESHOLD_C):
            ax.axvline(thr, color=RED, linestyle="--", linewidth=1.2, zorder=1)

    ax_q.stackplot(oat, q_econ, q_mech, colors=[AQUA, BLUE], zorder=2,
                   labels=["Economizer (free)", "Chiller (compressor)"])
    ax_q.set_ylabel("Cooling delivered [kW]", color=SECONDARY_INK, fontsize=9)
    ax_q.set_title("How the 150 kW load splits across the ASHRAE 90.1 economizer envelope",
                   color=PRIMARY_INK, fontsize=11, loc="left")
    ax_q.legend(fontsize=8.5, edgecolor=BASELINE, loc="center left")
    ax_q.set_ylim(0, 165)
    ax_q.text(2.0, 157, "Full free cooling\ncompressor off", fontsize=7.5, color=MUTED, ha="center")
    ax_q.text(7.0, 157, "Integrated / partial\neconomizer + chiller", fontsize=7.5, color=MUTED, ha="center")
    ax_q.text(12.2, 157, "Mechanical only", fontsize=7.5, color=MUTED, ha="center")

    ax_p.fill_between(oat, p_comp, color=YELLOW, zorder=2, alpha=0.9)
    ax_p.set_ylabel("Compressor power [kW]", color=SECONDARY_INK, fontsize=9)
    ax_p.set_xlabel("Outdoor air dry-bulb temperature [°C]", color=SECONDARY_INK, fontsize=9)
    ax_p.text(10.15, max(p_comp) * 0.55, f"{PARTIAL_FC_THRESHOLD_C:.0f} °C\nactivation",
              fontsize=7.5, color=RED)
    ax_p.text(3.85, max(p_comp) * 0.55, f"{FREE_COOLING_THRESHOLD_C:.0f} °C\nfull free cooling",
              fontsize=7.5, color=RED, ha="right")
    ax_p.text(0.0, -0.46, "At 4 °C the curves are continuous -- the economizer's rated duty there is "
                          "exactly the 150 kW it was sized for,\nso the partial band joins full free "
                          "cooling with no step. At 10 °C there IS a step: the economizer is still "
                          "delivering\n75 kW (half the load) when the fixed 90.1 changeover switches "
                          "it off, so compressor power jumps 6.4 -> 13.0 kW.\nThat step is left money "
                          "on the table, not physics -- see the floating-changeover item in "
                          "Suggestions for Improvement.",
              transform=ax_p.transAxes, fontsize=7.5, color=MUTED)

    fig.tight_layout()
    save_figure(fig, save_path, facecolor=SURFACE)
    plt.close(fig)
    print(f"Saved {save_path}")


def plot_period_performance(rows, save_path, title):
    """Saves a two-panel bar chart (COP, PUE) across the given rows (12
    months or 4 seasons), bars colored by mode -- mechanical vs free-
    cooling -- matching the annual chart's blue/aqua convention. Two
    stacked panels instead of one dual-axis plot, since COP [-] and PUE
    [-] use very different numeric ranges (COP up to ~80 in free-cooling
    vs PUE always close to 1) -- a dual y-axis chart is never the right
    call. Direct value labels are used instead of a hover layer, since
    this is a static PNG for a printed/PDF report, not an interactive
    chart."""
    import matplotlib.pyplot as plt

    BLUE = "#2a78d6"
    AQUA = "#1baf7a"
    YELLOW = "#eda100"
    SURFACE = "#fcfcfb"
    PRIMARY_INK = "#0b0b0b"
    SECONDARY_INK = "#52514e"
    MUTED = "#898781"
    GRID = "#e1e0d9"
    BASELINE = "#c3c2b7"

    labels = [r["label"] for r in rows]
    # Three regimes, three colors -- the partial band is neither free nor fully
    # mechanical and shouldn't be colored as either.
    MODE_COLOR = {"free_cooling": AQUA, "partial_free_cooling": YELLOW,
                  "mechanical": BLUE}
    colors = [MODE_COLOR[r["mode"]] for r in rows]
    cop = [_row_cop_total(r) for r in rows]
    pue = [r["PUE"] for r in rows]
    tue = [r["TUE"] for r in rows]
    x = np.arange(len(rows))

    fig, (ax_cop, ax_pue) = plt.subplots(2, 1, figsize=(9, 6), facecolor=SURFACE)

    for ax in (ax_cop, ax_pue):
        ax.set_facecolor(SURFACE)
        ax.grid(True, axis="y", color=GRID, linewidth=0.8, zorder=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color(BASELINE)
        ax.spines["bottom"].set_color(BASELINE)
        ax.tick_params(colors=MUTED, labelsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8)

    bars_cop = ax_cop.bar(x, cop, color=colors, width=0.6, zorder=2)
    ax_cop.set_ylabel("Total system COP [-]", color=SECONDARY_INK, fontsize=9)
    ax_cop.set_title(title, color=PRIMARY_INK, fontsize=11, loc="left")
    for b, v in zip(bars_cop, cop):
        ax_cop.text(b.get_x() + b.get_width() / 2, v, f"{v:.1f}", ha="center",
                    va="bottom", fontsize=7, color=SECONDARY_INK)

    # PUE and TUE share a panel: same units, comparable range, and TUE is a
    # fixed 1.20x of PUE, so the pair reads as one story. Mode stays encoded in
    # the bar COLOR; the metric is encoded by FILL (solid PUE vs outlined TUE)
    # so the two encodings don't collide.
    bars_pue = ax_pue.bar(x - 0.19, pue, color=colors, width=0.36, zorder=2)
    bars_tue = ax_pue.bar(x + 0.19, tue, facecolor="none", edgecolor=colors,
                          linewidth=1.4, hatch="///", width=0.36, zorder=2)
    ax_pue.set_ylabel("PUE / TUE [-]", color=SECONDARY_INK, fontsize=9)
    for bars, series, fmt in ((bars_pue, pue, "{:.2f}"), (bars_tue, tue, "{:.2f}")):
        for b, v in zip(bars, series):
            ax_pue.text(b.get_x() + b.get_width() / 2, v, fmt.format(v), ha="center",
                        va="bottom", fontsize=6.5, color=SECONDARY_INK)

    ax_pue.text(0.0, -0.34, "Aqua = full free cooling (compressor off)   Yellow = integrated/partial "
                            "free cooling   Blue = mechanical      Solid = PUE   Hatched = TUE\n"
                            f"TUE = ITUE x PUE with ITUE = {ITUE:.2f} ASSUMED (no server power "
                            f"model in this project); band {ITUE_BAND[0]:.2f}-{ITUE_BAND[1]:.2f}.",
                transform=ax_pue.transAxes, fontsize=7.5, color=MUTED)

    fig.tight_layout()
    save_figure(fig, save_path, facecolor=SURFACE)
    plt.close(fig)
    print(f"Saved {save_path}")


def _latex_escape(val):
    if isinstance(val, float):
        return f"{val:.2f}"
    s = str(val)
    return s.replace("_", "\\_") if s else "--"


def write_summary_latex_table(rows, out_path, caption, label):
    """Writes a compact tabularx snippet (label, T_air, mode, COP, PUE,
    compressor power, delivered duty) meant to be \\input{} directly into
    the report, rather than hand-copying numbers from the CSV -- so the
    report table always matches whatever main.py last computed."""
    cols = [
        ("label", "Period"), ("T_air_C", "$T_{air}$ [C]"), ("mode", "Mode"),
        ("P_compressor_kW", "$P_{comp}$ [kW]"), ("Q_delivered_kW", "$Q$ [kW]"),
        ("PUE", "PUE [-]"), ("TUE", "TUE [-]"),
    ]
    lines = []
    lines.append(r"\begin{table}[h!]")
    lines.append(rf"\caption{{{caption}}}")
    lines.append(rf"\label{{{label}}}")
    lines.append(r"\begin{tabularx}{\linewidth}{l" + "L" * (len(cols) - 1) + "}")
    lines.append(r"\toprule")
    lines.append(" & ".join(f"\\textbf{{{h}}}" for _, h in cols) + r" \\")
    lines.append(r"\midrule")
    for row in rows:
        cop_total = _row_cop_total(row)
        vals = []
        for key, _ in cols:
            if key == "PUE":
                vals.append(f"{row['PUE']:.2f} (COP {cop_total:.1f})")
            elif key == "TUE":
                vals.append(f"{row['TUE']:.2f}")
            elif key == "T_air_C":
                vals.append(f"{row['T_air_C']:.1f}")
            elif key in ("P_compressor_kW", "Q_delivered_kW"):
                v = row[key]
                vals.append(f"{v:.2f}" if v not in (None, "") else "--")
            else:
                vals.append(_latex_escape(row[key]))
        lines.append(" & ".join(vals) + r" \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabularx}")
    lines.append(r"\end{table}")

    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Saved {out_path}")


print()
print("=" * 72)
print("ANNUAL SIMULATION -- first pass (raw per-hour results, no aggregation yet)")
print("=" * 72)
annual_results = run_annual_simulation()
plot_annual_cop_load(annual_results)
plot_economizer_band()

print()
print("=" * 72)
print("MONTHLY / SEASONAL DETAILED REPORT")
print("=" * 72)
monthly_rows, seasonal_rows = run_monthly_seasonal_reports()

plot_period_performance(
    monthly_rows, os.path.join(FIG_PERF_MONTHLY, "monthly_cop_pue.png"),
    "Monthly-average total system COP, PUE and TUE (Champaign, IL TMY3)",
)
plot_period_performance(
    seasonal_rows, os.path.join(FIG_PERF_SEASONAL, "seasonal_cop_pue.png"),
    "Seasonal-average total system COP, PUE and TUE (Champaign, IL TMY3)",
)
_TUE_NOTE = (rf" TUE $=$ ITUE $\times$ PUE with ITUE $= {ITUE:.2f}$ \textbf{{assumed}} "
             rf"(band {ITUE_BAND[0]:.2f}--{ITUE_BAND[1]:.2f}); this project has no "
             rf"server-level power model, so TUE is PUE scaled by a constant.")
write_summary_latex_table(
    monthly_rows, os.path.join(TABLE_PERF_MONTHLY, "monthly_summary.tex"),
    "Monthly-average system performance summary." + _TUE_NOTE, "tab:monthly_summary",
)
write_summary_latex_table(
    seasonal_rows, os.path.join(TABLE_PERF_SEASONAL, "seasonal_summary.tex"),
    "Seasonal-average system performance summary." + _TUE_NOTE, "tab:seasonal_summary",
)
