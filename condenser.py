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
      D_h, A_flow, L, beta, Lambda, N_cp        (refrigerant side)
      b, L_w, N_cp_glycol                       (glycol side; see geometry block below)

    OUTPUTS (for teammates' component models):
      Q                               -> condenser heat duty (should ~ Q_evap + W_comp)
      h_out, p_out, T_out, subcool    -> EEV / receiver inlet state
      delta_p                         -> condenser refrigerant-side pressure loss
      delta_p_glycol                  -> condenser glycol-side pressure loss (Martin 1996)
      m_dot_glycol, cp_glycol         -> dry-cooler-loop pump sizing
      UA, LMTD                        -> plate-HX selection / sizing
    """

    def __init__(self, m_dot_refrigerant, p_in, refrigerant, h_in, subcooling,
                 T_glycol_in, T_glycol_out,
                 D_h, A_flow, L, beta, Lambda, N_cp,
                 b, L_w, N_cp_glycol,
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

        # Glycol-side channel geometry -- same corrugation family (D_h, beta,
        # Lambda) as the refrigerant side, drives the Martin (1996) single-
        # phase friction correlation (same one used in economizer.py).
        self.b = b                      # mean channel gap [m]        (= D_h*phi/2)
        self.L_w = L_w                  # plate width [m]
        self.N_cp_glycol = N_cp_glycol  # number of glycol channels [-]
        self.A_flow_glycol = N_cp_glycol * b * L_w
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
        self.delta_p_glycol = None

        self.calc_refrigerant_inlet_state()
        self.calc_heat_duty()
        self.calc_pressure_drop()
        self.calc_refrigerant_outlet_state()
        self.calc_glycol_mass_flow()
        self.calc_glycol_pressure_drop()
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

    # -- Single-phase friction: MARTIN (1996) generalised Leveque correlation --
    # Same correlation as economizer.py (VDI Heat Atlas, Sec. N6); the glycol
    # side here is single-phase liquid, so Han's two-phase correlation (used
    # for the refrigerant side above) does not apply.
    def _martin_dp(self, m_dot, A_flow, fluid, T_avg):
        rho = PropsSI('D', 'T', T_avg, 'P', 101325, fluid)
        mu = PropsSI('V', 'T', T_avg, 'P', 101325, fluid)
        G = m_dot / A_flow
        Re = G * self.D_h / mu
        phi = np.radians(self.beta)

        if Re < 2000:
            f0 = 16.0 / Re
            f1 = 149.0 / Re + 0.9625
        else:
            f0 = (1.56 * np.log(Re) - 3.0) ** (-2)
            f1 = 9.75 / Re ** 0.289

        inv_sqrt_f = (np.cos(phi) /
                      np.sqrt(0.045 * np.tan(phi) + 0.09 * np.sin(phi) + f0 / np.cos(phi))
                      + (1 - np.cos(phi)) / np.sqrt(3.8 * f1))
        f = 1.0 / inv_sqrt_f ** 2

        delta_p = 2 * f * G ** 2 * self.L / (rho * self.D_h)
        return delta_p, Re, f, G, rho, mu

    def calc_glycol_pressure_drop(self):
        T_glycol_avg = (self.T_glycol_in + self.T_glycol_out) / 2
        (self.delta_p_glycol, self.Re_glycol, self.f_glycol,
         self.G_glycol, self.rho_glycol, self.mu_glycol) = \
            self._martin_dp(self.m_dot_glycol, self.A_flow_glycol, self.glycol, T_glycol_avg)

    def calc_lmtd(self):
        # Counter-flow: hot refrigerant in vs glycol out, refrigerant out vs glycol in
        delta_T_1 = self.T_in - self.T_glycol_out
        delta_T_2 = self.T_out - self.T_glycol_in
        self.LMTD = (delta_T_1 - delta_T_2) / np.log(delta_T_1 / delta_T_2)

    def calc_ua(self):
        self.UA = self.Q / self.LMTD


if __name__ == "__main__":
    # Cycle-side operating conditions (R290, matches main.py's design point)
    R = "R290"
    T_o = 5 + 273.15    # evaporating temperature [K]
    T_c = 52 + 273.15   # condensing temperature [K] (raised from 45 C to keep a positive
                        # approach against the dry cooler's 40/46 C glycol loop)
    DT_SH = 10           # suction superheat assumed by the compressor [K]
    SUBCOOLING = 5        # condenser design subcooling [K]

    p_o = PropsSI("P", "T", T_o, "Q", 1, R)
    p_c = PropsSI("P", "T", T_c, "Q", 1, R)
    h_in = PropsSI("HMASS", "T", T_c + 20, "P", p_c, R)  # placeholder superheated discharge
    m_dot = 0.53  # refrigerant mass flow [kg/s] -- PLACEHOLDER (matches evaporator.py's demo)

    # ==================================================================
    # >>> GEOMETRY -- PLACEHOLDER, pending real SWEP/Kelvion datasheet <<<
    D_h    = 0.004    # [m]   4mm -- locked: pco/D_h=1.25 consistent w/ Han's tested geometry
    Lambda = 0.005    # [m]   5mm
    beta   = 30       # [deg]
    L      = 0.5      # [m]   ASSUMED, mid of typical 0.3-0.6m range
    N_cp   = 38       # ASSUMED to target G~25 kg/m2s -- RECHECK vs real m_dot_refrigerant
    A_flow = N_cp * (D_h * 1.17 / 2) * 0.2   # = 0.0178 m^2 (L_w=0.2m ASSUMED)

    b            = D_h * 1.17 / 2   # mean channel gap [m]     (same corrugation family)
    L_w          = 0.2              # plate width [m]          (matches refrigerant side)
    N_cp_glycol  = 24                # ASSUMED, same order as economizer.py's glycol side
    # ==================================================================

    cond = Condenser(m_dot_refrigerant=m_dot, p_in=p_c, refrigerant=R, h_in=h_in,
                      subcooling=SUBCOOLING, T_glycol_in=30 + 273.15, T_glycol_out=35 + 273.15,
                      D_h=D_h, A_flow=A_flow, L=L, beta=beta, Lambda=Lambda, N_cp=N_cp,
                      b=b, L_w=L_w, N_cp_glycol=N_cp_glycol)

    print(f"Duty           : {cond.Q/1e3:.2f} kW")
    print(f"Subcool        : {cond.subcool:.2f} K @ {cond.p_out/1000:.0f} kPa")
    print(f"Pressure drop  : {cond.delta_p/1000:.2f} kPa  (refrigerant side)")
    print(f"Glycol flow    : {cond.m_dot_glycol:.3f} kg/s  (cp = {cond.cp_glycol} J/kg-K)")
    print(f"LMTD / UA      : {cond.LMTD:.2f} K / {cond.UA:.0f} W/K")
    print("--- Martin (1996) glycol-side pressure drop ---")
    print(f"Glycol : G={cond.G_glycol:5.1f} kg/m2s  Re={cond.Re_glycol:6.0f}  "
          f"f={cond.f_glycol:.4f}  dP={cond.delta_p_glycol/1e3:5.2f} kPa  "
          f"v={cond.G_glycol/cond.rho_glycol:.2f} m/s")