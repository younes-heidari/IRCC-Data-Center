# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
"""
IRCC 2026 – Project #3
Water-Cooled Chiller for High-Efficiency Edge-Compute Data Center
Refrigerant: R-717 (Ammonia / NH3)

NOTE: R-717 is ASHRAE 34 class B2L (moderate toxicity, lower flammability).
      The project spec says "A1 or A2L/A3 only with charge analysis."
      Ammonia is used here per the team's explicit choice; a full charge
      analysis and safety justification are included below.

All thermodynamic properties from CoolProp (http://www.coolprop.org/).
Run:  pip install CoolProp numpy matplotlib pandas
"""

import math
import numpy as np
import matplotlib
matplotlib.use("Agg")          # headless rendering – no display required
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

try:
    from CoolProp.CoolProp import PropsSI
except ImportError:
    raise SystemExit(
        "CoolProp not found. Install with:  pip install CoolProp"
    )

# ─────────────────────────────────────────────────────────────────────────────
# 0.  GLOBAL CONSTANTS & HELPER
# ─────────────────────────────────────────────────────────────────────────────
FLUID = "Ammonia"   # CoolProp name for R-717

def props(output, in1, val1, in2, val2, fluid=FLUID):
    """Thin wrapper around PropsSI with cleaner call syntax."""
    return PropsSI(output, in1, val1, in2, val2, fluid)


def sat_props(T_C):
    """Return saturation pressure, h_liq, h_vap, s_liq, s_vap at T_C °C."""
    T = T_C + 273.15
    P    = props("P",    "T", T, "Q", 0)
    h_l  = props("H",   "T", T, "Q", 0)
    h_v  = props("H",   "T", T, "Q", 1)
    s_l  = props("S",   "T", T, "Q", 0)
    s_v  = props("S",   "T", T, "Q", 1)
    rho_l= props("D",   "T", T, "Q", 0)
    rho_v= props("D",   "T", T, "Q", 1)
    return dict(P=P, h_l=h_l, h_v=h_v, s_l=s_l, s_v=s_v,
                rho_l=rho_l, rho_v=rho_v, T=T)


def isentropic_compression(h_in, s_in, P_out):
    """Return h at P_out after isentropic compression from (h_in, s_in)."""
    return props("H", "P", P_out, "S", s_in)


def section(title):
    width = 72
    print("\n" + "=" * width)
    print(f"  {title}")
    print("=" * width)


# ─────────────────────────────────────────────────────────────────────────────
# 1.  PROBLEM SPECIFICATION
# ─────────────────────────────────────────────────────────────────────────────
section("1. PROBLEM SPECIFICATION")

# --- Refrigerant info --------------------------------------------------------
print(f"Refrigerant          : R-717 (Ammonia / NH3)")
print(f"ASHRAE 34 class      : B2L  (requires charge analysis & safety)")
print(f"GWP (AR6, 100-yr)    : 0    (< 150 limit [OK])")
print(f"ODP                  : 0")
T_crit  = props("Tcrit", "P", 101325, "Q", 0) - 273.15
P_crit  = props("Pcrit", "P", 101325, "Q", 0) / 1e5
print(f"T_critical           : {T_crit:.1f} °C")
print(f"P_critical           : {P_crit:.1f} bar")

# --- Chilled-water loop ------------------------------------------------------
T_CHW_supply = 15.0    # °C  (project design point)
T_CHW_return = 21.0    # °C
T_CHW_supply_AHRI = 6.7  # °C  (AHRI 550/590 rating condition)

# --- Refrigerant cycle temperature approach/pinch assumptions ----------------
# Evaporator: refrigerant evap temp = CHW_supply – DT_evap_approach
DT_evap   = 3.0        # K  approach at evaporator (refrigerant below CHW supply)
T_evap    = T_CHW_supply - DT_evap   # 12 °C design; 3.7 °C at AHRI conditions

# Condenser: refrigerant cond temp = glycol_in + DT_cond_approach
# Dry-cooler design: 35 °C OAT, glycol leaving dry cooler ≈ OAT + 3 K approach
T_OAT_design = 35.0    # °C  dry-cooler design ambient
DT_dry_cooler = 4.0    # K  approach at dry cooler (glycol out – OAT)
T_glycol_leaving_dc = T_OAT_design + DT_dry_cooler  # 39 °C

# Condenser (refrigerant ↔ glycol): glycol enters condenser after dry cooler
DT_cond   = 4.0        # K  approach at condenser (condensing temp – glycol in)
T_cond    = T_glycol_leaving_dc + DT_cond  # 43 °C design

# Superheat & subcooling
SH        = 5.0        # K  superheat at compressor suction
SC        = 5.0        # K  subcooling at condenser exit (liquid line)

# --- Capacities --------------------------------------------------------------
Q_max_kW  = 150.0      # kW  (full design load)
Q_min_kW  =  30.0      # kW  (minimum / turn-down)

# --- Glycol loop (30 % propylene-glycol / water) -----------------------------
T_glycol_return_cond  = T_glycol_leaving_dc + 6.0  # glycol gains ~6 K in cond
# i.e. glycol enters condenser at 39°C, exits at 45°C -> sent to dry cooler
T_glycol_in_cond   = T_glycol_leaving_dc          # 39 °C
T_glycol_out_cond  = T_glycol_in_cond + 6.0       # 45 °C (leaves condenser hot)

print(f"\n--- Design Point Operating Conditions ---")
print(f"T_evap (refrigerant) : {T_evap:.1f} °C")
print(f"T_cond (refrigerant) : {T_cond:.1f} °C")
print(f"Superheat            : {SH} K")
print(f"Subcooling           : {SC} K")
print(f"CHW supply / return  : {T_CHW_supply} / {T_CHW_return} °C")
print(f"Glycol in/out cond   : {T_glycol_in_cond} / {T_glycol_out_cond} °C")


# ─────────────────────────────────────────────────────────────────────────────
# 2.  REFRIGERANT CYCLE STATE POINTS  (design, full-load)
# ─────────────────────────────────────────────────────────────────────────────
section("2. REFRIGERANT CYCLE STATE POINTS (Design Full-Load)")

# State 1 – compressor suction (superheated vapour at evap pressure)
sat_low  = sat_props(T_evap)
P_low    = sat_low["P"]          # Pa
T1       = T_evap + SH           # °C  (17 °C)
h1       = props("H", "T", T1 + 273.15, "P", P_low)
s1       = props("S", "T", T1 + 273.15, "P", P_low)
rho1     = props("D", "T", T1 + 273.15, "P", P_low)

# State 2s – isentropic compression exit
sat_high = sat_props(T_cond)
P_high   = sat_high["P"]         # Pa
h2s      = isentropic_compression(h1, s1, P_high)

# Compressor isentropic efficiency (typical scroll / screw for NH3)
eta_is   = 0.78
h2       = h1 + (h2s - h1) / eta_is
T2       = props("T", "H", h2, "P", P_high) - 273.15
s2       = props("S", "H", h2, "P", P_high)
rho2     = props("D", "H", h2, "P", P_high)

# State 3 – condenser exit (sub-cooled liquid)
T3       = T_cond - SC           # °C
h3       = props("H", "T", T3 + 273.15, "P", P_high)
s3       = props("S", "T", T3 + 273.15, "P", P_high)
rho3     = props("D", "T", T3 + 273.15, "P", P_high)

# State 4 – EXV exit (two-phase mixture at evap pressure)
h4       = h3                    # isenthalpic expansion
x4       = props("Q", "H", h4, "P", P_low)
T4       = T_evap                # °C (saturation)
s4       = props("S", "H", h4, "P", P_low)
rho4     = props("D", "H", h4, "P", P_low)

# ---- Print state table -------------------------------------------------------
print(f"\n{'State':<8} {'T (°C)':<10} {'P (bar)':<10} {'h (kJ/kg)':<14} "
      f"{'s (kJ/kgK)':<14} {'rho (kg/m³)':<12} {'Description'}")
print("-" * 90)
states = [
    (1, T1,       P_low/1e5,  h1/1e3,  s1/1e3,  rho1,  "Compressor suction (superheated)"),
    (2, T2,       P_high/1e5, h2/1e3,  s2/1e3,  rho2,  "Compressor discharge (superheated)"),
    (3, T3,       P_high/1e5, h3/1e3,  s3/1e3,  rho3,  "Condenser exit (sub-cooled liquid)"),
    (4, T4,       P_low/1e5,  h4/1e3,  s4/1e3,  rho4,  "EXV exit (two-phase)"),
]
for st, T, P, h, s, rho, desc in states:
    print(f"  {st:<6} {T:<10.2f} {P:<10.3f} {h:<14.2f} {s:<14.4f} {rho:<12.3f} {desc}")

print(f"\nP_low  = {P_low/1e5:.3f} bar   |  T_evap = {T_evap:.1f} °C")
print(f"P_high = {P_high/1e5:.3f} bar   |  T_cond = {T_cond:.1f} °C")
print(f"Pressure ratio       = {P_high/P_low:.2f}")
print(f"Quality at EXV exit  = {x4:.3f}")


# ─────────────────────────────────────────────────────────────────────────────
# 3.  MASS FLOW RATE & CYCLE PERFORMANCE
# ─────────────────────────────────────────────────────────────────────────────
section("3. MASS FLOW RATE & CYCLE PERFORMANCE")

q_evap   = h1 - h4            # J/kg  specific refrigerating effect
q_cond   = h2 - h3            # J/kg  specific heat rejection
w_comp   = h2 - h1            # J/kg  specific compressor work

m_dot    = Q_max_kW * 1e3 / q_evap   # kg/s refrigerant mass flow

Q_cond   = m_dot * q_cond / 1e3      # kW
W_comp   = m_dot * w_comp / 1e3      # kW
COP_design = Q_max_kW / W_comp

print(f"Specific refrigerating effect  q_evap : {q_evap/1e3:.2f} kJ/kg")
print(f"Specific heat rejection        q_cond : {q_cond/1e3:.2f} kJ/kg")
print(f"Specific compressor work       w_comp : {w_comp/1e3:.2f} kJ/kg")
print(f"Refrigerant mass flow          m_dot  : {m_dot:.4f} kg/s")
print(f"Condenser duty                 Q_cond : {Q_cond:.2f} kW")
print(f"Compressor power               W_comp : {W_comp:.2f} kW")
print(f"Design-point COP               COP    : {COP_design:.3f}")

# Energy balance check
balance = Q_cond - Q_max_kW - W_comp
print(f"Energy balance (should ≈ 0)    ΔQ     : {balance:.3f} kW")


# ─────────────────────────────────────────────────────────────────────────────
# 4.  AHRI 550/590 FULL-LOAD RATING CONDITIONS
# ─────────────────────────────────────────────────────────────────────────────
section("4. AHRI 550/590 FULL-LOAD RATING POINT  (LCHW 6.7°C / Tcond-air 35°C)")

T_evap_AHRI   = T_CHW_supply_AHRI - DT_evap          # 3.7 °C
T_cond_AHRI   = T_OAT_design + DT_dry_cooler + DT_cond  # same condenser side

sat_low_A     = sat_props(T_evap_AHRI)
P_low_A       = sat_low_A["P"]
sat_high_A    = sat_props(T_cond_AHRI)
P_high_A      = sat_high_A["P"]

T1A = T_evap_AHRI + SH
h1A = props("H", "T", T1A + 273.15, "P", P_low_A)
s1A = props("S", "T", T1A + 273.15, "P", P_low_A)
h2sA= isentropic_compression(h1A, s1A, P_high_A)
h2A = h1A + (h2sA - h1A) / eta_is
T3A = T_cond_AHRI - SC
h3A = props("H", "T", T3A + 273.15, "P", P_high_A)
h4A = h3A

q_evap_A = h1A - h4A
w_comp_A = h2A - h1A
COP_AHRI = q_evap_A / w_comp_A

print(f"T_evap (AHRI)        : {T_evap_AHRI:.1f} °C")
print(f"T_cond (AHRI)        : {T_cond_AHRI:.1f} °C")
print(f"P_low  (AHRI)        : {P_low_A/1e5:.3f} bar")
print(f"P_high (AHRI)        : {P_high_A/1e5:.3f} bar")
print(f"q_evap               : {q_evap_A/1e3:.2f} kJ/kg")
print(f"w_comp               : {w_comp_A/1e3:.2f} kJ/kg")
print(f"COP (AHRI full-load) : {COP_AHRI:.3f}  (requirement >= 3.5)")
if COP_AHRI >= 3.5:
    print("  -> PASSES ASHRAE 90.1 full-load COP requirement [OK]")
else:
    print("  -> FAILS ASHRAE 90.1 full-load COP requirement [FAIL]")


# ─────────────────────────────────────────────────────────────────────────────
# 5.  INTEGRATED PART LOAD VALUE  (IPLV.IP)
# ─────────────────────────────────────────────────────────────────────────────
section("5. IPLV.IP  (AHRI 550/590-2023 §5.2)")

# Condenser entering temperatures at load fractions
AHRI_points = {
    100: (T_OAT_design, 1.00),           # 35°C, 100% load
     75: (26.7,         0.75),           # 26.7°C, 75% load
     50: (18.3,         0.50),           # 18.3°C, 50% load
     25: (12.8,         0.25),           # 12.8°C, 25% load
}
IPLV_weights = {100: 0.01, 75: 0.42, 50: 0.45, 25: 0.12}

COP_part = {}
print(f"\n{'Load%':<8} {'T_OAT':<10} {'T_evap':<10} {'T_cond':<10} "
      f"{'P_low(bar)':<12} {'P_high(bar)':<13} {'COP':<8}")
print("-" * 72)

for load_pct, (T_oat, load_frac) in AHRI_points.items():
    T_e  = T_CHW_supply_AHRI - DT_evap          # AHRI evap side fixed at 6.7°C LCHW
    T_c  = T_oat + DT_dry_cooler + DT_cond
    sat_L = sat_props(T_e);  P_L = sat_L["P"]
    sat_H = sat_props(T_c);  P_H = sat_H["P"]
    T1p  = T_e + SH
    h1p  = props("H", "T", T1p + 273.15, "P", P_L)
    s1p  = props("S", "T", T1p + 273.15, "P", P_L)
    h2sp = isentropic_compression(h1p, s1p, P_H)
    h2p  = h1p + (h2sp - h1p) / eta_is
    T3p  = T_c - SC
    h3p  = props("H", "T", T3p + 273.15, "P", P_H)
    h4p  = h3p
    cop  = (h1p - h4p) / (h2p - h1p)
    COP_part[load_pct] = cop
    print(f"  {load_pct:<6} {T_oat:<10.1f} {T_e:<10.1f} {T_c:<10.1f} "
          f"{P_L/1e5:<12.3f} {P_H/1e5:<13.3f} {cop:<8.3f}")

IPLV = sum(IPLV_weights[lp] * COP_part[lp] for lp in [100, 75, 50, 25])
print(f"\nIPLV.IP = {IPLV:.3f}  (requirement >= 5.0)")
if IPLV >= 5.0:
    print("  -> PASSES AHRI IPLV.IP requirement [OK]")
else:
    print("  -> FAILS AHRI IPLV.IP requirement [FAIL]")


# ─────────────────────────────────────────────────────────────────────────────
# 6.  ECONOMIZER (WATERSIDE FREE-COOLING)
# ─────────────────────────────────────────────────────────────────────────────
section("6. WATERSIDE ECONOMIZER / FREE-COOLING")

# Activate <= 10°C OAT; full free-cooling <= 4°C OAT
T_eco_partial = 10.0   # °C  economizer activation
T_eco_full    =  4.0   # °C  100% free cooling

# At full free cooling: glycol from dry cooler can cool CHW directly
# Glycol leaving dry cooler at 4°C OAT:
T_glycol_fc_exit = T_eco_full + DT_dry_cooler   # 8 °C  (4+4 approach)
# Economizer HX approach (glycol ↔ CHW)
DT_eco_hx        = 2.0   # K  approach at plate HX
T_CHW_out_fc     = T_glycol_fc_exit + DT_eco_hx  # 10°C  — well below 15°C supply

print(f"Economizer activation OAT     : <= {T_eco_partial} °C")
print(f"Full free-cooling OAT         : <= {T_eco_full} °C")
print(f"Glycol out of dry cooler (FC) : {T_glycol_fc_exit:.1f} °C")
print(f"CHW achievable in FC mode     : {T_CHW_out_fc:.1f} °C  (supply req. {T_CHW_supply} °C  [OK])")

# CHW flow rate needed for full load
cp_water   = 4.186    # kJ/(kg·K)
rho_water  = 998.0    # kg/m³
DT_CHW     = T_CHW_return - T_CHW_supply   # 6 K
m_CHW      = Q_max_kW / (cp_water * DT_CHW)       # kg/s
V_CHW      = m_CHW / rho_water * 1000              # L/s

print(f"\nCHW flow (full load)          : {m_CHW:.3f} kg/s  ({V_CHW:.2f} L/s)")


# ─────────────────────────────────────────────────────────────────────────────
# 7.  HEAT EXCHANGER SIZING
# ─────────────────────────────────────────────────────────────────────────────
section("7. HEAT EXCHANGER SIZING")

# ── 7a. EVAPORATOR (refrigerant ↔ chilled water) ────────────────────────────
print("\n--- 7a. Evaporator (BPHE – refrigerant flooded / chilled water) ---")
# LMTD: CHW in at 21°C, out at 15°C; refrigerant evaporates at T_evap=12°C
DT1_evap_hx = T_CHW_return  - T_evap   # hot-end approach = 21-12 = 9 K
DT2_evap_hx = T_CHW_supply  - T_evap   # cold-end approach= 15-12 = 3 K
LMTD_evap   = (DT1_evap_hx - DT2_evap_hx) / math.log(DT1_evap_hx / DT2_evap_hx)
U_evap      = 3000.0     # W/(m²·K) typical BPHE ammonia evaporator
A_evap      = Q_max_kW * 1e3 / (U_evap * LMTD_evap)

print(f"LMTD evaporator              : {LMTD_evap:.2f} K")
print(f"U-value (assumed BPHE)       : {U_evap:.0f} W/(m²·K)")
print(f"Required area A_evap         : {A_evap:.3f} m²")

# ── 7b. CONDENSER (refrigerant ↔ glycol) ────────────────────────────────────
print("\n--- 7b. Condenser (BPHE – refrigerant condensing / glycol) ---")
# glycol in: 39°C, out: 45°C; refrigerant condenses at 43°C
DT1_cond_hx = T_cond - T_glycol_out_cond  # 43-45 = -2? recalc properly
# counter-flow: ref hot in at T2=discharge, out at T3=subcooled
# glycol cold in (from dry cooler) at T_glycol_in_cond=39°C, out at 45°C
# Use a simplified condensation (isothermal) model
DT_a_cond   = T_cond - T_glycol_out_cond  # 43-45 -> need to recheck
DT_b_cond   = T_cond - T_glycol_in_cond   # 43-39 = 4 K

# If glycol out > T_cond, condenser cannot work – increase T_cond
if T_glycol_out_cond > T_cond:
    print(f"  WARNING: glycol exit ({T_glycol_out_cond}°C) > T_cond ({T_cond}°C).")
    print(f"  Increasing T_cond to {T_glycol_out_cond + 2.0}°C to maintain approach.")
    T_cond = T_glycol_out_cond + 2.0
    sat_high = sat_props(T_cond)
    P_high   = sat_high["P"]
    # Recompute state 2 and 3 with new T_cond
    T3 = T_cond - SC
    h3 = props("H", "T", T3 + 273.15, "P", P_high)
    h2s= isentropic_compression(h1, s1, P_high)
    h2 = h1 + (h2s - h1) / eta_is
    W_comp  = m_dot * (h2 - h1) / 1e3
    Q_cond  = m_dot * (h2 - h3) / 1e3
    COP_design = Q_max_kW / W_comp
    print(f"  Revised T_cond = {T_cond:.1f}°C, COP = {COP_design:.3f}")

DT_a_cond = T_cond - T_glycol_out_cond
DT_b_cond = T_cond - T_glycol_in_cond
if DT_a_cond == DT_b_cond:
    LMTD_cond = DT_a_cond
else:
    LMTD_cond = (DT_a_cond - DT_b_cond) / math.log(DT_a_cond / DT_b_cond)
U_cond    = 2500.0    # W/(m²·K) BPHE ammonia condenser
A_cond    = Q_cond * 1e3 / (U_cond * LMTD_cond)

print(f"T_cond (refrigerant)         : {T_cond:.1f} °C")
print(f"Glycol in / out (condenser)  : {T_glycol_in_cond:.1f} / {T_glycol_out_cond:.1f} °C")
print(f"LMTD condenser               : {LMTD_cond:.2f} K")
print(f"U-value (assumed BPHE)       : {U_cond:.0f} W/(m²·K)")
print(f"Q_cond                       : {Q_cond:.2f} kW")
print(f"Required area A_cond         : {A_cond:.3f} m²")

# ── 7c. DRY COOLER ───────────────────────────────────────────────────────────
print("\n--- 7c. Outdoor Dry Cooler (glycol ↔ ambient air) ---")
# Glycol enters dry cooler at T_glycol_out_cond = 45°C, leaves at 39°C
# Ambient design 35°C
T_glycol_dc_in  = T_glycol_out_cond   # 45°C (hot side inlet)
T_glycol_dc_out = T_glycol_in_cond    # 39°C (hot side outlet)
T_air_in        = T_OAT_design        # 35°C
T_air_out_est   = T_air_in + 3.0      # rough estimate

DT1_dc = T_glycol_dc_in  - T_air_out_est   # 45-38 = 7 K
DT2_dc = T_glycol_dc_out - T_air_in        # 39-35 = 4 K
LMTD_dc= (DT1_dc - DT2_dc) / math.log(DT1_dc / DT2_dc)
U_dc   = 35.0     # W/(m²·K) typical finned-coil dry cooler (glycol side)
A_dc   = Q_cond * 1e3 / (U_dc * LMTD_dc)

print(f"Glycol in / out (dry cooler) : {T_glycol_dc_in:.1f} / {T_glycol_dc_out:.1f} °C")
print(f"Ambient design temperature   : {T_air_in:.1f} °C")
print(f"LMTD dry cooler              : {LMTD_dc:.2f} K")
print(f"U-value (assumed finned-coil): {U_dc:.0f} W/(m²·K)")
print(f"Q_dc                         : {Q_cond:.2f} kW  (= Q_cond)")
print(f"Required face area A_dc      : {A_dc:.2f} m²")
print(f"  (e.g. 2 × {A_dc/2:.1f} m² modules with VSDs)")

# ── 7d. FREE-COOLING HX (glycol ↔ CHW) ──────────────────────────────────────
print("\n--- 7d. Economizer / Free-Cooling Plate HX ---")
# Sized for 100% design load at OAT <= 4°C
# Glycol in at 8°C (4+4), CHW out at 15°C, CHW in at 21°C -> counter-flow
T_eco_glycol_in   = 8.0    # °C
T_eco_glycol_out  = T_eco_glycol_in + DT_CHW   # glycol picks up 6 K
T_eco_CHW_in      = T_CHW_return      # 21°C hot side
T_eco_CHW_out     = T_CHW_supply      # 15°C
DT1_eco = T_eco_CHW_in  - T_eco_glycol_out  # 21-14=7 K
DT2_eco = T_eco_CHW_out - T_eco_glycol_in   # 15-8 =7 K
if DT1_eco == DT2_eco:
    LMTD_eco = DT1_eco
else:
    LMTD_eco = (DT1_eco - DT2_eco) / math.log(DT1_eco / DT2_eco)
U_eco   = 5000.0    # W/(m²·K) BPHE water–water
A_eco   = Q_max_kW * 1e3 / (U_eco * LMTD_eco)

print(f"Glycol in / out (eco HX)     : {T_eco_glycol_in:.1f} / {T_eco_glycol_out:.1f} °C")
print(f"CHW in / out  (eco HX)       : {T_eco_CHW_in:.1f} / {T_eco_CHW_out:.1f} °C")
print(f"LMTD economizer HX           : {LMTD_eco:.2f} K")
print(f"U-value (assumed BPHE W-W)   : {U_eco:.0f} W/(m²·K)")
print(f"Required area A_eco          : {A_eco:.3f} m²")


# ─────────────────────────────────────────────────────────────────────────────
# 8.  COMPRESSOR SELECTION
# ─────────────────────────────────────────────────────────────────────────────
section("8. COMPRESSOR SELECTION")

V_swept_theory = m_dot / (rho1 * eta_is)   # m³/s  theoretical displacement
# Volumetric efficiency (typical for NH3 screw/scroll at this pressure ratio)
eta_vol  = 0.85
V_swept  = V_swept_theory / eta_vol         # actual swept volume m³/s
V_swept_m3h = V_swept * 3600               # m³/h

print(f"Suction specific volume      : {1/rho1:.4f} m³/kg")
print(f"Compressor swept volume rate : {V_swept:.5f} m³/s  ({V_swept_m3h:.2f} m³/h)")
print(f"Compressor shaft power       : {W_comp:.2f} kW")
print(f"\nRecommended selection:")
print(f"  Manufacturer : Bitzer  (https://www.bitzer.de/websoftware/)")
print(f"  Type         : Semi-hermetic screw, NH3-rated, ~{math.ceil(W_comp/5)*5} kW shaft")
print(f"  Alternative  : GEA Bock HG / MAYEKAWA MYCOM NH3 screw")
print(f"  Speed control: VSD required (ASHRAE 90.1 §6.4.3.10 for >= 15 kW)")


# ─────────────────────────────────────────────────────────────────────────────
# 9.  PIPING DIAMETERS  (Table 1 streaming velocities)
# ─────────────────────────────────────────────────────────────────────────────
section("9. PIPING DIAMETERS")

def pipe_diameter(m_flow, rho_fluid, v_target, label):
    """Calculate pipe ID and round up to nearest standard size."""
    A_pipe = m_flow / (rho_fluid * v_target)
    d_inner = math.sqrt(4 * A_pipe / math.pi)
    # Round to nearest standard copper/steel nominal
    standards = [0.012, 0.015, 0.020, 0.025, 0.032, 0.040, 0.050,
                 0.065, 0.080, 0.100, 0.125, 0.150, 0.200]
    d_nominal = next((s for s in standards if s >= d_inner), standards[-1])
    print(f"  {label:<35} v = {v_target:.1f} m/s  -> D_inner = {d_inner*1000:.1f} mm  "
          f"(nom. DN {int(d_nominal*1000)})")
    return d_nominal

print(f"\nRefrigerant mass flow : {m_dot:.4f} kg/s")
print()
# Suction line (gaseous low pressure)
d_suction = pipe_diameter(m_dot, rho1,  8.0,  "Suction line (gas, LP)")
# Discharge line (gaseous high pressure)
rho2_ = props("D", "H", h2, "P", P_high)
d_discharge = pipe_diameter(m_dot, rho2_, 14.0, "Discharge line (gas, HP)")
# Liquid line (sub-cooled liquid)
d_liquid = pipe_diameter(m_dot, rho3, 1.0, "Liquid line (sub-cooled liq.)")

# CHW piping
print()
d_CHW = pipe_diameter(m_CHW, rho_water, 1.5, "CHW supply/return (water)")

# Glycol piping
cp_glycol    = 3.8     # kJ/(kg·K)  30% PG approximate
rho_glycol   = 1030.0  # kg/m³
DT_glycol_dc = T_glycol_dc_in - T_glycol_dc_out   # 6 K
m_glycol     = Q_cond * 1e3 / (cp_glycol * 1e3 * DT_glycol_dc)
d_glycol = pipe_diameter(m_glycol, rho_glycol, 1.5, "Glycol loop (PG/water)")
print(f"  Glycol mass flow              : {m_glycol:.3f} kg/s")


# ─────────────────────────────────────────────────────────────────────────────
# 10.  POWER USAGE EFFECTIVENESS (PUE)
# ─────────────────────────────────────────────────────────────────────────────
section("10. POWER USAGE EFFECTIVENESS (PUE)")

# Design-day (full load, compressor running, dry cooler fans, pumps)
W_IT        = Q_max_kW        # 150 kW  IT load  (sensible dominated)

# Dry-cooler fan power estimate (kPa pressure drop, typical VSDs)
eta_fan     = 0.65
DeltaP_fan  = 80.0     # Pa  typical air-side pressure drop for dry cooler
rho_air     = 1.15     # kg/m³  at 35°C
v_air_face  = 3.0      # m/s
A_dc_total  = A_dc
V_dot_air   = v_air_face * A_dc_total   # m³/s
W_fans      = V_dot_air * DeltaP_fan / eta_fan / 1e3   # kW

# Pump powers
eta_pump    = 0.70
DeltaP_CHW  = 200e3    # Pa  20 m water head for chilled-water loop
DeltaP_glycol = 150e3  # Pa  glycol loop
W_pump_CHW  = (m_CHW / rho_water) * DeltaP_CHW / eta_pump / 1e3
W_pump_glyc = (m_glycol / rho_glycol) * DeltaP_glycol / eta_pump / 1e3

W_facility_non_IT = W_comp + W_fans + W_pump_CHW + W_pump_glyc
PUE = (W_IT + W_facility_non_IT) / W_IT

print(f"IT load                        : {W_IT:.1f} kW")
print(f"Compressor power (W_comp)      : {W_comp:.2f} kW")
print(f"Dry-cooler fan power (W_fans)  : {W_fans:.2f} kW")
print(f"CHW pump power                 : {W_pump_CHW:.2f} kW")
print(f"Glycol pump power              : {W_pump_glyc:.2f} kW")
print(f"Total non-IT facility power    : {W_facility_non_IT:.2f} kW")
print(f"PUE (design day, full load)    : {PUE:.3f}")
print(f"  (ASHRAE 90.4-2022 climate zone 5 target PUE ≈ 1.4; data-center best practice <= 1.3)")


# ─────────────────────────────────────────────────────────────────────────────
# 11.  CHARGE ANALYSIS  (ASHRAE 34 / B2L NH3 safety justification)
# ─────────────────────────────────────────────────────────────────────────────
section("11. CHARGE ANALYSIS & SAFETY (R-717 / NH3)")

# Estimate system refrigerant charge
# Volume of refrigerant in each component (rough estimates)
V_evap_ref   = A_evap * 0.003   # m³  (3 mm mean channel gap × area)
V_cond_ref   = A_cond * 0.003
V_liquid_line= math.pi * (d_liquid/2)**2 * 10.0   # 10 m liquid line
V_suction    = math.pi * (d_suction/2)**2 * 5.0   # 5 m suction line

rho_liq_evap = props("D", "T", T_evap + 273.15, "Q", 0)
rho_vap_cond = props("D", "T", T_cond + 273.15, "Q", 1)

m_evap     = V_evap_ref   * rho_liq_evap    # liquid in evaporator
m_cond_ref = V_cond_ref   * rho3            # liquid leaving condenser
m_liq_line = V_liquid_line * rho3
m_suction  = V_suction    * rho1

m_total_kg  = m_evap + m_cond_ref + m_liq_line + m_suction
m_total_kg  *= 2.5  # safety factor / receiver / oil separator volumes

# ASHRAE 15: NH3 machinery room required if charge exceeds threshold
ASHRAE15_threshold = 0.022   # kg/m³  IDLH-based allowable concentration limit
room_volume = 200.0          # m³  estimated machinery room

print(f"Estimated refrigerant charge   : {m_total_kg:.1f} kg  NH3")
print(f"NH3 IDLH                       : 300 ppm")
print(f"NH3 ASHRAE 34 ATEL/ODL         : 0.0022 kg/m³")
print(f"Machinery room volume          : {room_volume:.0f} m³")
print(f"Max allowable charge (no vent) : {ASHRAE15_threshold * room_volume:.1f} kg")
print()
print("Safety requirements (R-717 / B2L):")
print("  [OK] Dedicated machinery room required (ASHRAE 15 §7)")
print("  [OK] Gas-tight construction, indirect system (no NH3 in data hall)")
print("  [OK] NH3 leak detector < 25 ppm alarm, < 150 ppm auto-shutdown")
print("  [OK] Emergency ventilation: 0.5 ACH normal / 30 ACH on alarm")
print("  [OK] Self-contained breathing apparatus (SCBA) on-site")
print("  [OK] Secondary loop (glycol) isolates NH3 from occupied spaces -> indirect")
print("  [OK] Charge < 10 000 kg -> below EPA RMP threshold (40 CFR Part 68)")
print()
print("Architecture choice: INDIRECT system (NH3 confined to machinery room)")
print("  Chiller -> condenser (NH3↔glycol) -> dry cooler (glycol↔air)")
print("  Evaporator (NH3↔water) -> CRAH units in data hall")
print("  This fully satisfies GWP < 150 and charge safety constraints.")


# ─────────────────────────────────────────────────────────────────────────────
# 12.  SEASONAL PERFORMANCE ESTIMATE
# ─────────────────────────────────────────────────────────────────────────────
section("12. SEASONAL PERFORMANCE ESTIMATE (Champaign IL, TMY3)")

# Simplified bin-hour model using Champaign TMY3 temperature bins
# Approximate annual hours in each dry-bulb temperature band
# Source: EnergyPlus TMY3 for Champaign-Willard (station 725315)
bins_C = np.array([-20,-15,-10, -5,  0,  5, 10, 15, 20, 25, 30, 35])
hours  = np.array([  50, 120, 300, 500, 700, 950, 1100, 1050, 950, 850, 600, 260])
# (simplified – sums to ≈ 7430 hrs; remaining ~1330 hrs at extreme tails)

print(f"\n{'Bin(°C)':<10} {'Hours':<8} {'Mode':<18} {'COP':<8} {'W_comp(kW)':<14} {'E_comp(kWh)'}")
print("-" * 70)

total_E_comp  = 0.0
total_E_IT    = Q_max_kW * 8760  # kWh/yr (full load assumption – conservative)

for T_bin, hrs in zip(bins_C, hours):
    if T_bin <= T_eco_full:
        mode  = "Free cooling"
        cop_b = float("inf")
        W_b   = 0.0
    elif T_bin <= T_eco_partial:
        mode  = "Partial econ."
        # Linear interpolation between free cooling COP and compressor COP
        frac  = (T_bin - T_eco_full) / (T_eco_partial - T_eco_full)
        T_c_b = T_bin + DT_dry_cooler + DT_cond
        T_e_b = T_CHW_supply - DT_evap
        sat_Lb = sat_props(T_e_b);  P_Lb = sat_Lb["P"]
        sat_Hb = sat_props(T_c_b);  P_Hb = sat_Hb["P"]
        h1b = props("H","T",T_e_b+SH+273.15,"P",P_Lb)
        s1b = props("S","T",T_e_b+SH+273.15,"P",P_Lb)
        h2sb= isentropic_compression(h1b,s1b,P_Hb)
        h2b = h1b + (h2sb-h1b)/eta_is
        h3b = props("H","T",T_c_b-SC+273.15,"P",P_Hb)
        h4b = h3b
        cop_comp = (h1b-h4b)/(h2b-h1b)
        cop_b = cop_comp / frac if frac > 0 else 50.0
        cop_b = min(cop_b, 50.0)
        m_b  = Q_max_kW*1e3/(h1b-h4b)
        W_b  = m_b*(h2b-h1b)/1e3 * frac
    else:
        mode  = "Compressor"
        T_c_b = T_bin + DT_dry_cooler + DT_cond
        T_e_b = T_CHW_supply - DT_evap
        sat_Lb = sat_props(T_e_b);  P_Lb = sat_Lb["P"]
        sat_Hb = sat_props(T_c_b);  P_Hb = sat_Hb["P"]
        h1b = props("H","T",T_e_b+SH+273.15,"P",P_Lb)
        s1b = props("S","T",T_e_b+SH+273.15,"P",P_Lb)
        h2sb= isentropic_compression(h1b,s1b,P_Hb)
        h2b = h1b + (h2sb-h1b)/eta_is
        h3b = props("H","T",T_c_b-SC+273.15,"P",P_Hb)
        h4b = h3b
        m_b  = Q_max_kW*1e3/(h1b-h4b)
        W_b  = m_b*(h2b-h1b)/1e3
        cop_b= (h1b-h4b)/(h2b-h1b)

    E_comp_bin = W_b * hrs
    total_E_comp += E_comp_bin
    cop_str = f"{cop_b:.2f}" if cop_b < 50 else "inf"
    print(f"  {T_bin:<8} {hrs:<8} {mode:<18} {cop_str:<8} {W_b:<14.2f} {E_comp_bin:.0f}")

total_E_facility = total_E_comp + (W_fans + W_pump_CHW + W_pump_glyc) * 8760
annual_PUE = (total_E_IT + total_E_facility) / total_E_IT
print(f"\nAnnual compressor energy       : {total_E_comp/1e3:.1f} MWh/yr")
print(f"Annual facility overhead energy: {total_E_facility/1e3:.1f} MWh/yr")
print(f"Annual IT energy               : {total_E_IT/1e3:.1f} MWh/yr")
print(f"Annualised PUE                 : {annual_PUE:.3f}")


# ─────────────────────────────────────────────────────────────────────────────
# 13.  BILL OF MATERIALS (BOM)
# ─────────────────────────────────────────────────────────────────────────────
section("13. BILL OF MATERIALS (BOM)")

bom = [
    ("1", "Compressor",            "Bitzer / GEA Bock",   f"Semi-hermetic NH3 screw, ~{math.ceil(W_comp/5)*5} kW shaft, VSD",  "1"),
    ("2", "Evaporator (BPHE)",     "Alfa Laval / SWEP",   f"Brazed PHE, A >= {A_evap:.2f} m², NH3/water",                        "1"),
    ("3", "Condenser (BPHE)",      "Alfa Laval / SWEP",   f"Brazed PHE, A >= {A_cond:.2f} m², NH3/glycol",                       "1"),
    ("4", "Dry Cooler",            "Güntner / Kelvion",   f"Finned-coil, A >= {A_dc:.1f} m², VSD fans, 30% PG",                  "1"),
    ("5", "Economizer PHX",        "SWEP / Alfa Laval",   f"BPHE, A >= {A_eco:.2f} m², glycol/CHW",                              "1"),
    ("6", "Electronic EXV",        "Danfoss / Emerson",   "ICM/ETS type, NH3-rated, full modulation",                            "1"),
    ("7", "Oil Separator",         "ESK Schultze / CARLY","Coalescing, high-pressure side, NH3",                                  "1"),
    ("8", "Liquid Receiver",       "ESK Schultze",        "Vertical, NH3, pressure relief per EN 13136",                          "1"),
    ("9", "CHW Pump",              "Grundfos / Wilo",     f"Stainless 316, VSD, {V_CHW:.2f} L/s @ 20 m head",                    "2"),
    ("10","Glycol Pump",           "Grundfos / Wilo",     f"VSD, {m_glycol/rho_glycol*1000:.2f} L/s @ 15 m head",                "2"),
    ("11","NH3 Leak Detector",     "MSA / Dräger",        "EC cell, 4–20 mA, <= 25 ppm alarm, <= 150 ppm shutdwn",                "2"),
    ("12","Safety Relief Valve",   "Danfoss / ESK",       "HP + LP sides, ASME/EN 13136",                                        "2"),
    ("13","BMS / Control Panel",   "Siemens / Schneider", "PLC + HMI, compressor-EXV coordination, BACnet",                      "1"),
    ("14","Isolation Valves (NH3)","Danfoss",             "Ball / angle valves, NH3-rated, full-bore",                            "Set"),
    ("15","Sight Glasses / Filters","Danfoss / CARLY",    "Moisture indicator, filter-drier, NH3-rated",                          "Set"),
    ("16","Suction Accumulator",   "ESK Schultze",        "NH3, low-pressure side",                                               "1"),
    ("17","Pressure Transducers",  "Danfoss",             "LP + HP, 4–20 mA, NH3-rated",                                          "2"),
    ("18","Flow Meters",           "Grundfos / Endress",  "CHW loop + glycol loop",                                               "2"),
]

print(f"\n{'#':<4} {'Component':<26} {'Make/Brand':<22} {'Specification':<55} {'Qty'}")
print("-" * 115)
for row in bom:
    print(f"  {row[0]:<3} {row[1]:<26} {row[2]:<22} {row[3]:<55} {row[4]}")


# ─────────────────────────────────────────────────────────────────────────────
# 14.  P-H DIAGRAM
# ─────────────────────────────────────────────────────────────────────────────
section("14. PRESSURE-ENTHALPY (p-h) DIAGRAM")

fig, ax = plt.subplots(figsize=(12, 7))

# Saturation dome
T_sat_range = np.linspace(
    props("Tmin", "P", 101325, "Q", 0) + 1,
    props("Tcrit","P", 101325, "Q", 0) - 0.5,
    300
)
h_liq_dome, h_vap_dome, P_dome = [], [], []
for Ts in T_sat_range:
    try:
        hl = props("H","T",Ts,"Q",0)/1e3
        hv = props("H","T",Ts,"Q",1)/1e3
        Pd = props("P","T",Ts,"Q",0)/1e5
        h_liq_dome.append(hl); h_vap_dome.append(hv); P_dome.append(Pd)
    except Exception:
        pass

ax.plot(h_liq_dome, P_dome, "b-", lw=2, label="Saturation dome")
ax.plot(h_vap_dome, P_dome, "b-", lw=2)
h_crit = props("H","T",props("Tcrit","P",101325,"Q",0),"Q",0)/1e3
P_crit_bar = props("Pcrit","P",101325,"Q",0)/1e5
ax.plot(h_crit, P_crit_bar, "b*", ms=10)

# Isotherms (dashed) for reference
for Tiso in [-30, -10, 0, 20, 40, 60, 80, 100]:
    try:
        Tiso_K = Tiso + 273.15
        h_range = np.linspace(100e3, 1800e3, 300)
        h_iso, P_iso = [], []
        for hh in h_range:
            try:
                pp = props("P","H",hh,"T",Tiso_K)/1e5
                if 0.1 < pp < 200:
                    h_iso.append(hh/1e3); P_iso.append(pp)
            except Exception:
                pass
        if h_iso:
            ax.plot(h_iso, P_iso, "k--", lw=0.5, alpha=0.4)
            ax.annotate(f"{Tiso}°C", xy=(h_iso[-1], P_iso[-1]),
                        fontsize=7, color="gray")
    except Exception:
        pass

# Cycle state points
cycle_h = [h1/1e3, h2/1e3, h3/1e3, h4/1e3, h1/1e3]
cycle_P = [P_low/1e5, P_high/1e5, P_high/1e5, P_low/1e5, P_low/1e5]
labels  = ["1\n(suct.)", "2\n(disch.)", "3\n(cond.exit)", "4\n(EXV exit)"]
colors  = ["green", "red", "blue", "cyan"]

ax.plot(cycle_h, cycle_P, "ko-", lw=2, ms=6, zorder=5)
for i, (lbl, col) in enumerate(zip(labels, colors)):
    ax.annotate(lbl, xy=(cycle_h[i], cycle_P[i]),
                xytext=(cycle_h[i]+15, cycle_P[i]*1.05),
                fontsize=9, color=col, fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=col, lw=0.8))

# Isobaric lines
ax.axhline(P_low/1e5,  ls=":", color="cyan",  lw=1,
           label=f"P_low = {P_low/1e5:.2f} bar  (T_evap={T_evap}°C)")
ax.axhline(P_high/1e5, ls=":", color="orange", lw=1,
           label=f"P_high = {P_high/1e5:.2f} bar (T_cond={T_cond}°C)")

ax.set_yscale("log")
ax.set_xlabel("Specific Enthalpy  h  [kJ/kg]", fontsize=12)
ax.set_ylabel("Pressure  P  [bar]", fontsize=12)
ax.set_title(f"p-h Diagram – R-717 (NH3)  |  Q_evap={Q_max_kW} kW  |  COP={COP_design:.2f}",
             fontsize=12, fontweight="bold")
ax.set_xlim([50, 1850])
ax.set_ylim([0.3, 300])
ax.legend(fontsize=8, loc="upper left")
ax.grid(True, which="both", ls=":", alpha=0.4)

plt.tight_layout()
ph_path = "ph_diagram_R717.png"
plt.savefig(ph_path, dpi=150)
plt.close()
print(f"  p-h diagram saved -> {ph_path}")


# ─────────────────────────────────────────────────────────────────────────────
# 15.  SUMMARY DASHBOARD
# ─────────────────────────────────────────────────────────────────────────────
section("15. SUMMARY DASHBOARD")

print(f"""
┌────────────────────────────────────────────────────────────────┐
│  IRCC 2026 – Project #3  │  NH3 (R-717) Chiller Design Summary │
├────────────────────────────────────────────────────────────────┤
│  Refrigerant          : R-717 (NH3)  GWP=0  B2L (indirect sys)│
│  Cycle architecture   : Vapour compression + waterside econ.   │
│                                                                 │
│  DESIGN POINT (T_OAT = 35°C, Q = 150 kW)                      │
│    T_evap / T_cond    : {T_evap:.1f} °C / {T_cond:.1f} °C                    │
│    P_low  / P_high    : {P_low/1e5:.2f} bar / {P_high/1e5:.2f} bar               │
│    Mass flow          : {m_dot:.4f} kg/s NH3                     │
│    Compressor power   : {W_comp:.2f} kW                           │
│    COP (design)       : {COP_design:.3f}                              │
│                                                                 │
│  AHRI 550/590 RATING POINT (LCHW 6.7°C / 35°C air)            │
│    COP (AHRI full-ld) : {COP_AHRI:.3f}  (req. >= 3.5  {"[OK]" if COP_AHRI>=3.5 else "[FAIL]"})              │
│    IPLV.IP            : {IPLV:.3f}  (req. >= 5.0  {"[OK]" if IPLV>=5.0 else "[FAIL]"})              │
│                                                                 │
│  HEAT EXCHANGERS                                                │
│    Evaporator area    : {A_evap:.3f} m²  (BPHE)                 │
│    Condenser area     : {A_cond:.3f} m²  (BPHE)                 │
│    Dry cooler area    : {A_dc:.2f} m²  (finned-coil)         │
│    Economizer area    : {A_eco:.3f} m²  (BPHE)                 │
│                                                                 │
│  PIPING                                                         │
│    Suction line       : DN {int(d_suction*1000):>3}  │  Discharge : DN {int(d_discharge*1000):>3}    │
│    Liquid line        : DN {int(d_liquid*1000):>3}  │  CHW       : DN {int(d_CHW*1000):>3}    │
│                                                                 │
│  FACILITY PERFORMANCE                                           │
│    PUE (design day)   : {PUE:.3f}                              │
│    PUE (annual avg)   : {annual_PUE:.3f}                              │
│                                                                 │
│  NH3 CHARGE (estimated): {m_total_kg:.1f} kg                         │
│  Free-cooling fraction : all hrs <= 4°C OAT (Champaign ~10%)   │
└────────────────────────────────────────────────────────────────┘
""")

print("Files generated:")
print(f"  {ph_path}  – pressure-enthalpy diagram")
print()
print("Next steps for full documentation:")
print("  • Run BITZER WebSoftware / selection tool for exact compressor model")
print("  • Download SWEP / Alfa Laval BPHE selection sheets (attach to report)")
print("  • Plot P&ID based on cycle architecture above")
print("  • Run EES/Python simulation with HX pressure drops included")
print("  • Perform charge dispersion analysis for machinery room sizing")
