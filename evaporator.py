import numpy as np
from CoolProp.CoolProp import PropsSI


class Evaporator:
    """Models the chilled-water evaporator: boils refrigerant using heat
    from the chilled water loop, refrigerant leaves superheated.

    INPUTS (from the cycle side -- provided by teammates/operating conditions):
      m_dot_refrigerant : refrigerant mass flow [kg/s]      (from compressor)
      p_in              : evaporating pressure [Pa]          (operating condition)
      Q                 : cooling load / heat duty [W]       (design spec, e.g. 150 kW)
      refrigerant       : fluid name [str]                   (e.g. "R290")
      T_water_in        : chilled-water return temp [K]      (spec: 21 C)
      T_water_out       : chilled-water supply temp [K]      (spec: 15 C)
      h_in              : refrigerant inlet enthalpy [J/kg]  (from EEV / condenser outlet)

    INPUTS (heat-exchanger geometry -- my responsibility, from datasheet):
      D_h, A_flow, L, beta, Lambda, N_cp  (see geometry block below)

    OUTPUTS (for teammates' component models):
      h_out, p_out, T_out, superheat  -> compressor suction state
      delta_p                         -> evaporator refrigerant-side pressure loss
      m_dot_water, cp_water           -> chilled-water pump sizing
      UA, LMTD                        -> plate-HX selection / sizing
    """

    def __init__(self, m_dot_refrigerant, p_in, Q, refrigerant,
                 T_water_in, T_water_out, h_in,
                 D_h, A_flow, L, beta, Lambda, N_cp):
        # General properties
        self.R = refrigerant
        self.m_dot_refrigerant = m_dot_refrigerant
        self.Q = Q

        # ==================================================================
        # >>> BRAZED-PLATE HEAT-EXCHANGER GEOMETRY <<<
        # Edit these with the real datasheet values (SWEP / Alfa Laval /
        # Kelvion) once a plate is selected. They drive the Han et al. (2003)
        # pressure-drop correlation; typical ranges + Han validity in [].
        # ------------------------------------------------------------------
        self.D_h = D_h          # hydraulic diameter [m]       (~3-5 mm; = 2b/phi)
        self.A_flow = A_flow    # refrigerant flow area [m^2]  (= N_cp * b * L_w)
        self.L = L              # plate flow length L_v [m]    (~0.3-0.6 m)
        self.beta = beta        # chevron angle [deg]          (from horizontal; Han: 20-45)
        self.Lambda = Lambda    # corrugation pitch p_co [m]   (~5-7 mm)
        self.N_cp = N_cp        # number of refrigerant channels [-]  (Han Eq. 12 factor)
        # ==================================================================

        # Refrigerant inlet (from EEV, two-phase)
        self.p_in = p_in
        self.h_in = h_in
        self.T_sat_in = None

        # Refrigerant outlet (to compressor, superheated)
        self.h_out = self.h_in + self.Q / self.m_dot_refrigerant
        self.p_out = None
        self.T_out = None
        self.superheat = None

        # Water side
        self.T_water_in = T_water_in
        self.T_water_out = T_water_out
        self.m_dot_water = None
        self.cp_water = None

        # Sizing results
        self.LMTD = None
        self.UA = None

        self.delta_p = None

        self.calc_refrigerant_inlet_state()
        self.calc_pressure_drop()
        self.calc_refrigerant_outlet_state()
        self.calc_water_mass_flow()
        self.calc_lmtd()
        self.calc_ua()

    # Refrigerant temp/enthalpy at inlet pressure (saturation-fixed since two-phase)
    def calc_refrigerant_inlet_state(self):
        self.T_sat_in = PropsSI('T', 'P', self.p_in, 'Q', 0, self.R)

    def calc_pressure_drop(self):
        G = self.m_dot_refrigerant / self.A_flow

        # Saturated properties at inlet pressure
        h_f = PropsSI('HMASS', 'P', self.p_in, 'Q', 0, self.R)
        h_g = PropsSI('HMASS', 'P', self.p_in, 'Q', 1, self.R)
        rho_l = PropsSI('D', 'P', self.p_in, 'Q', 0, self.R)
        rho_v = PropsSI('D', 'P', self.p_in, 'Q', 1, self.R)
        mu_l = PropsSI('V', 'P', self.p_in, 'Q', 0, self.R)

        x_in = (self.h_in - h_f) / (h_g - h_f)

        # Han, Lee, Kim (2003) Eq. 17-19: friction factor geometry parameters
        beta_rad = np.radians(self.beta)
        Ge3 = 64710 * (self.Lambda / self.D_h)**(-5.27) * (np.pi / 2 - beta_rad)**(-3.03)
        Ge4 = -1.314 * (self.Lambda / self.D_h)**(-0.62) * (np.pi / 2 - beta_rad)**(-0.47)

        # Boiling zone: Eq. 12 integrated over quality
        L_boiling = self.L * (h_g - self.h_in) / (self.h_out - self.h_in)
        x = np.linspace(x_in, 1.0, 200)

        # Eq. 11: equivalent mass flux; Eq. 20: equivalent Reynolds number
        G_eq = G * ((1 - x) + x * (rho_l / rho_v)**0.5)
        Re_eq = G_eq * self.D_h / mu_l
        f_tp = Ge3 * Re_eq**Ge4
        # Eq. 12: dP_fr = f * (L_v * N / D_h) * G_Eq^2 / rho_f.
        # N is HAN'S TEST-RIG refrigerant channel count (2), a fixed calibration
        # constant the friction factor was regressed against -- NOT the design
        # channel count self.N_cp. The design count enters only through
        # G = m_dot / A_flow (the per-channel mass flux); using self.N_cp here
        # over-predicts the boiling drop by N_cp/2. (Verified against Han's data
        # in condenser.py; see the condensation-paper companion correlation.)
        dpdz_tp = f_tp * G_eq**2 / (rho_l * self.D_h)

        N_HAN_CAL = 2
        delta_p_boiling = N_HAN_CAL * np.trapezoid(dpdz_tp, x) * L_boiling / (1.0 - x_in)

        # Superheating zone: single-phase vapor, Fanning friction factor.
        # Smooth-channel (Blasius) estimate -- this zone is short so its
        # contribution is small. No N_cp here: this is a physical single-channel
        # f, unlike Han's two-phase f which is calibrated to require N_cp.
        L_superheating = self.L - L_boiling
        h_avg_sh = (h_g + self.h_out) / 2
        rho_sh = PropsSI('D', 'HMASS', h_avg_sh, 'P', self.p_in, self.R)
        mu_sh = PropsSI('V', 'HMASS', h_avg_sh, 'P', self.p_in, self.R)
        Re_sh = G * self.D_h / mu_sh
        f_sh = 16 / Re_sh if Re_sh < 2000 else 0.079 * Re_sh**(-0.25)
        delta_p_superheating = 2 * f_sh * G**2 * L_superheating / (rho_sh * self.D_h)

        self.delta_p = delta_p_boiling + delta_p_superheating
        self.p_out = self.p_in - self.delta_p

    def calc_refrigerant_outlet_state(self):
        self.T_out = PropsSI('T', 'HMASS', self.h_out, 'P', self.p_out, self.R)
        T_sat_out = PropsSI('T', 'P', self.p_out, 'Q', 1, self.R)
        self.superheat = self.T_out - T_sat_out

    # Water mass flow needed to absorb Q while cooling 21 -> 15 C
    def calc_water_mass_flow(self):
        T_water_avg = (self.T_water_in + self.T_water_out) / 2
        self.cp_water = round(PropsSI('CPMASS', 'T', T_water_avg, 'P', 101325, 'Water'), 1)
        delta_T_water = self.T_water_in - self.T_water_out
        self.m_dot_water = self.Q / (self.cp_water * delta_T_water)

    def calc_lmtd(self):
        # Counter-flow: water in vs refrigerant out, water out vs refrigerant in
        delta_T_1 = self.T_water_in - self.T_out
        delta_T_2 = self.T_water_out - self.T_sat_in
        self.LMTD = (delta_T_1 - delta_T_2) / np.log(delta_T_1 / delta_T_2)

    def calc_ua(self):
        self.UA = self.Q / self.LMTD