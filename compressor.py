import numpy as np
from CoolProp.CoolProp import PropsSI

class Compressor:

    def __init__(self, eta_is, eta_vol, zeta_v,temp_in, d, z, s, n, pressure_in, pressure_out, refrigerant):
        self.W_v = None
        self.w_v = None
        self.T_2 = None
        self.h_2 = None
        self.h_2is = None
        self.s_1 = None
        self.h_1 = None
        self.m_dot = None
        self.rho_1 = None
        self.V_dot_disp = None
        self.V_cyl = None
        self.V_ dot_1 = None
        self.R = refrigerant
        self.is_eff = eta_is
        self.vol_eff = eta_vol
        self.T_1 = temp_in
        self.d = d
        self.z = z
        self.stroke = s
        self.n = n
        self.zeta_v = zeta_v

        self.p_1 = pressure_in
        self.p_2 = pressure_out

        self.pressure_ratio = self.p_2 / self.p_1
        self.calc_displacement()
        self.calc_inlet_state()
        self.calc_outlet_state()

    def calc_inlet_state(self):
        self.rho_1 = PropsSI('D', 'T', self.T_1, 'P', self.p_1, self.R)
        self.h_1 = PropsSI('HMASS', 'T', self.T_1, 'P', self.p_1, self.R)
        self.s_1 = PropsSI('SMASS', 'T', self.T_1, 'P', self.p_1, self.R)
        self.V_dot_1 = self.V_dot_disp * self.vol_eff
        self.m_dot = self.V_dot_1 * self.rho_1

    def calc_displacement(self):
        self.V_cyl = np.pi / 4 * self.d**2
        self.V_dot_disp = self.V_cyl * self.stroke * self.z * self.n / 60  # conversion from m^3/h to m^3/s

    def calc_outlet_state(self):
        self.h_2is = PropsSI('HMASS', 'SMASS', self.s_1, 'P', self.p_2, self.R)
        self.h_2 = self.h_1 + ((self.h_2is - self.h_1) * (1 + self.zeta_v)) / self.is_eff
        self.T_2 = PropsSI('T', 'HMASS', self.h_2, 'P', self.p_2, self.R)
        self.w_v = (self.h_2-self.h_1)/(1 + self.zeta_v)
        self.W_v = self.w_v * self.m_dot


screw = Compressor(eta_is=0.8, eta_vol=0.959, zeta_v=0.1, temp_in=288.15, d=0.075, z=4, s=0.055, n=1750, pressure_in=5.5e5, pressure_out=12.2e5, refrigerant="R290")
print("Displacement:", round(screw.V_dot_disp*3600, 1), "m^3/h")
print("Volume flow inlet:", round(screw.V_dot_1*3600, 1), "m^3/h")
print("Mass flow:", round(screw.m_dot,4), "kg/s")
print("Mass flow:", round(screw.m_dot*3600,1), "kg/h")
print("Pressure ratio:", round(screw.pressure_ratio,2),"-")
print("Discharge Temp.:", round(screw.T_2,2),"K")
print("Specific compressor work:", round(screw.w_v,2), "J/kg")
print("Compressor work:", round(screw.W_v/1000,2), "kW")
