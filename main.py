"""Integration driver: chains the four refrigerant-loop components
(compressor -> condenser -> EEV -> evaporator) around ONE shared mass
flow, taken from the compressor's own geometry -- instead of each
component file guessing its own placeholder m_dot independently.

The economizer is NOT part of this chain: it is a parallel, refrigerant-
free waterside loop (glycol <-> chilled water) that only runs when the
compressor is off, so it has no shared m_dot with the cycle below.
"""
from CoolProp.CoolProp import PropsSI

from compressor import Compressor
from condenser import Condenser
from eev import EEV
from evaporator import Evaporator

R = "R290"
T_o = 5 + 273.15          # evaporating temperature [K]
T_c = 52 + 273.15         # condensing temperature [K] (raised from 45 C so the condenser's
                          # saturation/subcooled approach stays positive against the dry
                          # cooler's vendor-confirmed 40/46 C glycol loop -- see dry_cooler.py)
DT_SH = 10                # suction superheat assumed by the compressor [K]
SUBCOOLING = 5            # condenser design subcooling [K] (matches eev.py)

p_o = PropsSI("P", "T", T_o, "Q", 1, R)
p_c = PropsSI("P", "T", T_c, "Q", 1, R)

# ---------------------------------------------------------------------
# 1) Compressor -- sets the ONE shared mass flow for the whole cycle
# ---------------------------------------------------------------------
compressor = Compressor(eta_is=0.8, eta_vol=0.959, zeta_v=0.1, temp_in=T_o + DT_SH,
                         d=0.075, z=4, s=0.055, n=1750,
                         pressure_in=p_o, pressure_out=p_c, refrigerant=R)
m_dot = compressor.m_dot

# ---------------------------------------------------------------------
# 2) Condenser -- same corrugation-family placeholder geometry as the
#    evaporator (D_h, Lambda, beta, L), pending a real datasheet.
# ---------------------------------------------------------------------
D_h, Lambda, beta, L = 0.004, 0.005, 30, 0.5
b, L_w = D_h * 1.17 / 2, 0.2
N_cp_cond = 38
A_flow_cond = N_cp_cond * b * L_w
N_cp_glycol = 24

condenser = Condenser(m_dot_refrigerant=m_dot, p_in=compressor.p_2, refrigerant=R,
                       h_in=compressor.h_2, subcooling=SUBCOOLING,
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
# Report
# ---------------------------------------------------------------------
Q_TARGET = 150e3

print("=== SHARED CYCLE MASS FLOW ===")
print(f"m_dot (from compressor) : {m_dot:.3f} kg/s ({m_dot*3600:.0f} kg/h)")
print()
print("=== COMPRESSOR ===")
print(f"Discharge : {compressor.T_2-273.15:.1f} C @ {compressor.p_2/1000:.0f} kPa")
print(f"Work      : {compressor.W_v/1e3:.2f} kW")
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
print("=== CYCLE CLOSURE CHECK ===")
print(f"Cooling duty delivered : {evaporator.Q/1e3:.1f} kW   (design target: {Q_TARGET/1e3:.0f} kW)")
print(f"Shortfall vs target    : {(1 - evaporator.Q/Q_TARGET)*100:.1f} %")
print(f"Energy balance check   : Q_cond ({condenser.Q/1e3:.1f} kW) "
      f"vs Q_evap + W_comp ({(evaporator.Q + compressor.W_v)/1e3:.1f} kW)")
