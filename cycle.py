from compressor import RecipCompressor
from CoolProp.CoolProp import PropsSI

# Given Values
R = "R290"
T_o = 5 # Evaporating temperature in °C
T_c = 45 # Condensing temperature in °C
DT_SH = 10 # Superheat in K - EEV

# Adjust temperatures for EOS to Kelvin
T_o += 273.15 # K
T_c += 273.15 # K

# Retrieve saturation pressure for set temperature levels
p_o = PropsSI("P", "T", T_o, "Q", 1, R)
p_c = PropsSI("P", "T", T_c, "Q", 1, R)

# Setup for test case - BITZER 4GEP-30Z-40P
recip = RecipCompressor(zeta_v=0, temp_in=T_o+DT_SH,
                    # d=0.075, z=4, s=0.055,
                    displacement=102,
                    # eta_vol=0.959,
                    mass_flow=1072,
                    # eta_is = 0.73,
                    n=1750, pressure_in=p_o, pressure_out=p_c, refrigerant=R,
                    power_input=20400)

# Print relevant data
print(f"Displacement: {recip.V_dot_disp * 3600:.1f} m^3/h")
print(f"Volume flow inlet: {recip.V_dot_1 * 3600:.1f} m^3/h")
print(f"Mass flow: {recip.m_dot * 3600:.1f} kg/h")
print(f"Pressure ratio: {recip.pressure_ratio:.2f} -")
print(f"Discharge Temp.: {recip.T_2-273.15:.2f} °C\n")

print(f"Compressor work: {recip.W_v / 1000:.2f} kW")
print(f"Power input: {recip.P_el/1000:.2f} kW\n")

print(f"Volumetric efficiency: {recip.eta_vol:.4f} -")
print(f"Isentropic efficiency: {recip.eta_is:.2f} -")
