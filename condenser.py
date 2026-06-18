import numpy as np
from CoolProp.CoolProp import PropsSI


class Condenser:
    """Models the refrigerant-to-glycol condenser: rejects heat from the
    compressor-discharge refrigerant to the intermediate propylene-glycol
    loop (which carries it to the outdoor dry cooler). Refrigerant enters
    superheated, is desuperheated, condensed, and leaves subcooled.

    INPUTS (from the cycle side -- provided by teammates/operating conditions):
      m_dot_refrigerant : refrigerant mass flow [kg/s]        (from compressor)
      p_in              : condensing pressure [Pa]            (= compressor discharge p)
      refrigerant       : fluid name [str]                    (e.g. "R290")
      h_in              : refrigerant inlet enthalpy [J/kg]   (from compressor discharge)
      subcooling        : design subcooling at outlet [K]     (my choice -- to justify)
      T_glycol_in       : glycol supply temp from dry cooler [K]
      T_glycol_out      : glycol return temp to dry cooler [K]
      glycol            : secondary-fluid name [str]          (CoolProp INCOMP brine)

    INPUTS (heat-exchanger geometry -- my responsibility, from datasheet):
      D_h, A_flow, L, beta, Lambda, N_cp  (see geometry block below)

    OUTPUTS (for teammates' component models):
      Q                               -> condenser heat duty (should ~ Q_evap + W_comp)
      h_out, p_out, T_out, subcool    -> EEV / receiver inlet state
      delta_p                         -> condenser refrigerant-side pressure loss
      m_dot_glycol, cp_glycol         -> dry-cooler-loop pump sizing
      UA, LMTD                        -> plate-HX selection / sizing
    """

    def __init__(self, m_dot_refrigerant, p_in, refrigerant, h_in, subcooling,
                 T_glycol_in, T_glycol_out,
                 D_h, A_flow, L, beta, Lambda, N_cp,
                 glycol="INCOMP::MPG[0.30]"):
        # General properties
        self.R = refrigerant
        self.m_dot_refrigerant = m_dot_refrigerant
        self.Q = None           # heat duty [W] -- computed (subcooling is the design input)

        # ==================================================================
        # >>> BRAZED-PLATE HEAT-EXCHANGER GEOMETRY <<<
        # Edit these with the real datasheet values (SWEP / Alfa Laval /
        # Kelvion) once a plate is selected. They drive the Han et al. (2003)
        # condensation pressure-drop correlation; typical ranges in [].
        # ------------------------------------------------------------------
        self.D_h = D_h          # hydraulic diameter [m]       (~3-5 mm; = 2b/phi)
        self.A_flow = A_flow    # refrigerant flow area [m^2]  (= N_cp * b * L_w)
        self.L = L              # plate flow length L_v [m]    (~0.3-0.6 m)
        self.beta = beta        # chevron angle [deg]          (from horizontal; Han: 20-45)
        self.Lambda = Lambda    # corrugation pitch p_co [m]   (~5-7 mm)
        self.N_cp = N_cp        # number of refrigerant channels [-]  (Han pressure-drop factor)
        # ==================================================================

        # Refrigerant inlet (from compressor, superheated vapor)
        self.p_in = p_in
        self.h_in = h_in
        self.T_in = None
        self.T_sat = None       # condensing (saturation) temperature at p_in

        # Refrigerant outlet (to EEV / receiver, subcooled liquid)
        self.subcooling = subcooling   # design subcooling [K] -- to be justified
        self.h_out = None
        self.p_out = None
        self.T_out = None
        self.subcool = None            # actual subcooling at outlet pressure [K]

        # Glycol side (intermediate loop to the outdoor dry cooler)
        self.glycol = glycol
        self.T_glycol_in = T_glycol_in
        self.T_glycol_out = T_glycol_out
        self.m_dot_glycol = None
        self.cp_glycol = None

        # Sizing results
        self.LMTD = None
        self.UA = None

        self.delta_p = None

        self.calc_refrigerant_inlet_state()
        self.calc_heat_duty()
        self.calc_pressure_drop()
        self.calc_refrigerant_outlet_state()
        self.calc_glycol_mass_flow()
        self.calc_lmtd()
        self.calc_ua()

    # Superheated inlet temperature and condensing (saturation) temperature at p_in
    def calc_refrigerant_inlet_state(self):
        self.T_in = PropsSI('T', 'HMASS', self.h_in, 'P', self.p_in, self.R)
        self.T_sat = PropsSI('T', 'P', self.p_in, 'Q', 1, self.R)

    # Outlet enthalpy is fixed by the design subcooling; the duty then follows.
    def calc_heat_duty(self):
        # Subcooling is referenced to the bubble point at the condensing pressure.
        T_bubble = PropsSI('T', 'P', self.p_in, 'Q', 0, self.R)
        T_out_design = T_bubble - self.subcooling
        # Liquid enthalpy is ~pressure-independent, so evaluate at p_in.
        self.h_out = PropsSI('HMASS', 'T', T_out_design, 'P', self.p_in, self.R)
        self.Q = self.m_dot_refrigerant * (self.h_in - self.h_out)

    def calc_pressure_drop(self):
        G = self.m_dot_refrigerant / self.A_flow

        # Saturated properties at condensing pressure
        h_f = PropsSI('HMASS', 'P', self.p_in, 'Q', 0, self.R)
        h_g = PropsSI('HMASS', 'P', self.p_in, 'Q', 1, self.R)
        rho_l = PropsSI('D', 'P', self.p_in, 'Q', 0, self.R)
        rho_v = PropsSI('D', 'P', self.p_in, 'Q', 1, self.R)
        mu_l = PropsSI('V', 'P', self.p_in, 'Q', 0, self.R)

        # Zone lengths split by enthalpy fraction (h_in > h_g > h_f > h_out):
        #   desuperheating (vapor) -> condensing (two-phase) -> subcooling (liquid)
        dh_total = self.h_in - self.h_out
        L_desuperheat = self.L * (self.h_in - h_g) / dh_total
        L_condense = self.L * (h_g - h_f) / dh_total
        L_subcool = self.L * (h_f - self.h_out) / dh_total

        # --- Condensing zone: Han, Lee, Kim (2003) condensation correlation ---
        # Han, Lee, Kim, "The characteristics of condensation in brazed plate
        # heat exchangers with different chevron angles," J. Korean Phys. Soc.
        # 43(1):66-73 (2003). f = Ge3 * Re_eq^Ge4 (Eqs. 23-25); equivalent mass
        # flux G_eq (Eq. 27, with G = G_c the per-channel flux) and equivalent
        # Reynolds number Re_eq (Eq. 26). Valid for Re_eq ~ 300-4000.
        beta_rad = np.radians(self.beta)
        Ge3 = 3521.1 * (self.Lambda / self.D_h)**4.17 * (np.pi / 2 - beta_rad)**(-7.75)
        Ge4 = -1.024 * (self.Lambda / self.D_h)**0.0925 * (np.pi / 2 - beta_rad)**(-1.3)

        x = np.linspace(0.0, 1.0, 200)
        G_eq = G * ((1 - x) + x * (rho_l / rho_v)**0.5)   # Eq. 27
        Re_eq = G_eq * self.D_h / mu_l                    # Eq. 26
        f_tp = Ge3 * Re_eq**Ge4                           # Eqs. 23-25
        # Eq. 12: dP_fr = f * (L_v * N / D_h) * G_eq^2 / rho_f, integrated over
        # quality. N is HAN'S TEST-RIG refrigerant channel count (2), the value
        # the friction factor was regressed against -- it is a fixed calibration
        # constant, NOT our design channel count. Verified against the paper's
        # Fig. 5 (~4-7 kPa @ G_c=34): N=2 matches, dropping it under-predicts 2x.
        # The design channel count must NOT appear here: across parallel single-
        # pass channels the port-to-port drop is independent of channel number;
        # it enters only through G_c = m_dot / A_flow.
        N_HAN_CAL = 2
        dpdz_tp = N_HAN_CAL * f_tp * G_eq**2 / (rho_l * self.D_h)
        delta_p_condense = np.trapezoid(dpdz_tp, x) * L_condense

        # --- Desuperheating zone: single-phase vapor (Fanning friction) ---
        h_avg_v = (self.h_in + h_g) / 2
        rho_v_sh = PropsSI('D', 'HMASS', h_avg_v, 'P', self.p_in, self.R)
        mu_v_sh = PropsSI('V', 'HMASS', h_avg_v, 'P', self.p_in, self.R)
        Re_v = G * self.D_h / mu_v_sh
        f_v = 16 / Re_v if Re_v < 2000 else 0.079 * Re_v**(-0.25)
        delta_p_desuperheat = 2 * f_v * G**2 * L_desuperheat / (rho_v_sh * self.D_h)

        # --- Subcooling zone: single-phase liquid (Fanning friction) ---
        h_avg_l = (h_f + self.h_out) / 2
        rho_l_sc = PropsSI('D', 'HMASS', h_avg_l, 'P', self.p_in, self.R)
        mu_l_sc = PropsSI('V', 'HMASS', h_avg_l, 'P', self.p_in, self.R)
        Re_l = G * self.D_h / mu_l_sc
        f_l = 16 / Re_l if Re_l < 2000 else 0.079 * Re_l**(-0.25)
        delta_p_subcool = 2 * f_l * G**2 * L_subcool / (rho_l_sc * self.D_h)

        self.delta_p = delta_p_desuperheat + delta_p_condense + delta_p_subcool
        self.p_out = self.p_in - self.delta_p

    def calc_refrigerant_outlet_state(self):
        self.T_out = PropsSI('T', 'HMASS', self.h_out, 'P', self.p_out, self.R)
        T_sat_out = PropsSI('T', 'P', self.p_out, 'Q', 0, self.R)
        self.subcool = T_sat_out - self.T_out

    # Glycol mass flow needed to carry Q while heating T_glycol_in -> T_glycol_out
    def calc_glycol_mass_flow(self):
        T_glycol_avg = (self.T_glycol_in + self.T_glycol_out) / 2
        self.cp_glycol = round(PropsSI('CPMASS', 'T', T_glycol_avg, 'P', 101325, self.glycol), 1)
        delta_T_glycol = self.T_glycol_out - self.T_glycol_in
        self.m_dot_glycol = self.Q / (self.cp_glycol * delta_T_glycol)

    def calc_lmtd(self):
        # Counter-flow: hot refrigerant in vs glycol out, refrigerant out vs glycol in
        delta_T_1 = self.T_in - self.T_glycol_out
        delta_T_2 = self.T_out - self.T_glycol_in
        self.LMTD = (delta_T_1 - delta_T_2) / np.log(delta_T_1 / delta_T_2)

    def calc_ua(self):
        self.UA = self.Q / self.LMTD