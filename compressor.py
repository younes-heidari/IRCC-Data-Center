import numpy as np
from CoolProp.CoolProp import PropsSI

class Compressor:

    def __init__(self, eta_is, eta_vol, zeta_v,temp_in, d, z, s, n, pressure_in, pressure_out, refrigerant):
        """Calculate the thermodynamic properties of a chosen compressor."""
        # General properties
        self.R = refrigerant
        self.w_v = None
        self.W_v = None
        self.m_dot = None
        self.V_dot_disp = None
        self.V_cyl = None
        self.is_eff = eta_is
        self.vol_eff = eta_vol
        self.d = d
        self.z = z
        self.stroke = s
        self.n = n
        self.zeta_v = zeta_v

        # Inlet State (1)
        self.T_1 = temp_in
        self.p_1 = pressure_in
        self.h_1 = None
        self.s_1 = None
        self.rho_1 = None
        self.V_dot_1 = None

        # Outlet State (2)
        self.T_2 = None
        self.p_2 = pressure_out
        self.h_2is = None
        self.h_2 = None

        self.pressure_ratio = self.p_2 / self.p_1

        # Calculation methods for completing properties and EOS
        self.calc_displacement()
        self.calc_inlet_state()
        self.calc_outlet_state()

    def calc_inlet_state(self):
        """Calculate thermodynamic properties of the inlet state based on the set temperature and pressure of the inlet state. Mass flow is calculated using a set volumetric efficiency."""
        self.rho_1 = PropsSI('D', 'T', self.T_1, 'P', self.p_1, self.R)
        self.h_1 = PropsSI('HMASS', 'T', self.T_1, 'P', self.p_1, self.R)
        self.s_1 = PropsSI('SMASS', 'T', self.T_1, 'P', self.p_1, self.R)
        self.V_dot_1 = self.V_dot_disp * self.vol_eff
        self.m_dot = self.V_dot_1 * self.rho_1

    def calc_displacement(self):
        """Calculate displacement by usage of geometric values and speed."""
        self.V_cyl = np.pi / 4 * self.d**2
        self.V_dot_disp = self.V_cyl * self.stroke * self.z * self.n / 60  # conversion from m^3/h to m^3/s

    def calc_outlet_state(self):
        """Calculate thermodynamic properties for outlet of the compressor (by using an isentropic efficiency). For the compressor work, a factor zeta_v is considered to take the heat loss into account."""
        self.h_2is = PropsSI('HMASS', 'SMASS', self.s_1, 'P', self.p_2, self.R)
        self.h_2 = self.h_1 + ((self.h_2is - self.h_1) * (1 + self.zeta_v)) / self.is_eff
        self.T_2 = PropsSI('T', 'HMASS', self.h_2, 'P', self.p_2, self.R)
        self.w_v = (self.h_2-self.h_1)/(1 + self.zeta_v)
        self.W_v = self.w_v * self.m_dot

# Given Values
R = "R290"
T_o = 5 # Evaporating temperature in °C
T_c = 45 # Condensing temperature in °C
DT_SH = 10 # Superheat in K

# Adjust temperatures for EOS to Kelvin
T_o += 273.15 # K
T_c += 273.15 # K

# Retrieve saturation pressure for set temperature levels
p_o = PropsSI("P", "T", T_o, "Q", 1, R)
p_c = PropsSI("P", "T", T_c, "Q", 1, R)

# Setup for test case
recip = Compressor(eta_is=0.8, eta_vol=0.959, zeta_v=0.1, temp_in=T_o+DT_SH,
                   d=0.075, z=4, s=0.055, n=1750, # only for recip (see comment below)
                   pressure_in=p_o, pressure_out=p_c, refrigerant=R)

# Print relevant data
print(f"Displacement: {recip.V_dot_disp * 3600:.1f} m^3/h")
print(f"Volume flow inlet: {recip.V_dot_1 * 3600:.1f} m^3/h")
# print(f"Mass flow: {recip.m_dot:.4f} kg/s")
print(f"Mass flow: {recip.m_dot * 3600:.1f} kg/h")
print(f"Pressure ratio: {recip.pressure_ratio:.2f} -")
print(f"Discharge Temp.: {recip.T_2-273.15:.2f} °C")
print(f"Specific compressor work: {recip.w_v/1000:.2f} kJ/kg")
print(f"Compressor work: {recip.W_v / 1000:.2f} kW")

### Remarks: Currently, some measures regarding the cylinder are set above (diameter, stroke, cylinder amount and speed).
### This is only applicable for a recip compressor. Otherwise we can also just use the displacement or create a subclass
### for the different types of compressors (scroll, screw, recip - centrifugal is probably not needed)

### Probably I'll switch from isentropic efficiency to compressor work (or power input) and from volumetric efficiency
### to mass flow as these will most likely be set by the compressor supplier

### zeta_v is a best guess by the script of Prof. Bradshaw (we'll have to set a value based on our compressor choice