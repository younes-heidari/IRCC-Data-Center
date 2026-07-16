import numpy as np
from CoolProp.CoolProp import PropsSI


def rate_economizer(UA, m_dot_water, cp_water, T_water_in,
                    m_dot_glycol, cp_glycol, T_glycol_in):
    """RATING model (epsilon-NTU, counter-flow): given INSTALLED hardware (a
    fixed UA) and the two inlet conditions, return the duty the exchanger
    actually delivers.

    This is the inverse of the Economizer class above/below it, which is a
    SIZING model: that one is handed a duty and all four temperatures and
    back-solves the UA needed. Sizing answers "how big must it be?" at the one
    design condition; it cannot answer "what does the exchanger I already
    bought do at 7 C ambient?" -- both outlet temperatures are unknown there,
    so LMTD has nothing to stand on. Hence epsilon-NTU.

    Used for the INTEGRATED (PARTIAL) free-cooling band, 4 < OAT <= 10 C
    (ASHRAE 90.1 s6.5.1 activation threshold): the economizer pre-cools the
    returning chilled water by whatever the ambient allows, and the chiller
    trims the remainder.

    Returns (Q_w, eps, NTU, Cr).
    """
    C_water = m_dot_water * cp_water          # [W/K]
    C_glycol = m_dot_glycol * cp_glycol
    C_min, C_max = min(C_water, C_glycol), max(C_water, C_glycol)
    Cr = C_min / C_max
    NTU = UA / C_min

    # Counter-flow effectiveness. The Cr -> 1 branch is not a corner case here:
    # this exchanger is deliberately balanced (6 K glycol rise mirrors the 6 K
    # CHW span), so C_water ~= C_glycol and the general formula goes 0/0.
    if abs(1.0 - Cr) < 1e-6:
        eps = NTU / (1.0 + NTU)
    else:
        x = np.exp(-NTU * (1.0 - Cr))
        eps = (1.0 - x) / (1.0 - Cr * x)

    Q = eps * C_min * (T_water_in - T_glycol_in)
    return max(Q, 0.0), eps, NTU, Cr


class Economizer:
    """Models the waterside free-cooling economizer: a glycol-to-chilled-water
    plate heat exchanger placed in PARALLEL with the chiller. When the outdoor
    air is cold enough, the dry cooler chills the intermediate propylene-glycol
    loop directly and this economizer transfers that "coolth" to the chilled
    water, cooling it 21 C -> 15 C with the compressor OFF (free cooling).

    KEY DIFFERENCE vs. the evaporator / condenser models
    ----------------------------------------------------
    Both streams here are SINGLE-PHASE liquids (water on one side, glycol on the
    other) -- there is NO refrigerant and NO phase change. The Han, Lee & Kim
    (2003) two-phase correlation used in evaporator.py / condenser.py therefore
    does NOT apply. The refrigerant-side pressure drop is replaced by the
    single-phase chevron-plate friction correlation of MARTIN (1996) (VDI Heat
    Atlas, Sec. N6), evaluated INDEPENDENTLY on each stream.

    INPUTS (cycle / operating conditions):
      Q             : free-cooling design duty [W]          (= 100% design load, 150 kW)
      T_water_in    : chilled-water return temp [K]         (spec: 21 C)
      T_water_out   : chilled-water supply temp [K]         (spec: 15 C)
      T_glycol_in   : glycol supply from dry cooler [K]     (cold; my design choice)
      T_glycol_out  : glycol return to dry cooler [K]       (warmed; my design choice)
      U             : assumed overall HTC [W/m^2 K]         (glycol/water BPHE, ~3000-4000)

    INPUTS (plate geometry -- from datasheet / same corrugation family as evap+cond):
      D_h, L, beta, Lambda, b, L_w, N_cp_water, N_cp_glycol  (see geometry block)

    OUTPUTS (for the rest of the design):
      m_dot_water, cp_water           -> CHW pump sizing (shared with evaporator loop)
      m_dot_glycol, cp_glycol         -> dry-cooler-loop pump sizing (free-cooling mode)
      LMTD, UA, A_required, N_plates  -> plate-HX selection / sizing
      delta_p_water, delta_p_glycol   -> both-side pressure loss (pump heads)
    """

    def __init__(self, Q, T_water_in, T_water_out, T_glycol_in, T_glycol_out,
                 D_h, L, beta, Lambda, b, L_w, N_cp_water, N_cp_glycol,
                 U=3500.0, glycol="INCOMP::MPG[0.30]", water="Water"):
        # General
        self.Q = Q
        self.glycol = glycol
        self.water = water
        self.U = U

        # ==================================================================
        # >>> CHEVRON-PLATE GEOMETRY <<<
        # Same corrugation family as the evaporator/condenser (D_h, beta,
        # Lambda) so plate technology is consistent, but on a LARGER, single-
        # wall frame (see report -- DW / 45 bar not required for a glycol/water
        # duty). Drives the Martin (1996) single-phase friction correlation.
        # ------------------------------------------------------------------
        self.D_h = D_h              # hydraulic diameter [m]           (~3-5 mm; = 2b/phi)
        self.L = L                  # plate flow length L_v [m]
        self.beta = beta            # chevron angle [deg]  (used directly in Martin, see note)
        self.Lambda = Lambda        # corrugation pitch p_co [m]       (carried for reference)
        self.b = b                  # mean channel gap [m]             (= D_h*phi/2)
        self.L_w = L_w              # plate width [m]
        self.N_cp_water = N_cp_water        # number of water channels  [-]
        self.N_cp_glycol = N_cp_glycol      # number of glycol channels [-]
        # ==================================================================

        # Water side (hot stream, being cooled 21 -> 15)
        self.T_water_in = T_water_in
        self.T_water_out = T_water_out
        self.m_dot_water = None
        self.cp_water = None
        self.A_flow_water = N_cp_water * b * L_w

        # Glycol side (cold stream, being warmed by the CHW)
        self.T_glycol_in = T_glycol_in
        self.T_glycol_out = T_glycol_out
        self.m_dot_glycol = None
        self.cp_glycol = None
        self.A_flow_glycol = N_cp_glycol * b * L_w

        # Sizing results
        self.LMTD = None
        self.UA = None
        self.A_required = None
        self.N_plates = None

        # Pressure drops (one per stream)
        self.delta_p_water = None
        self.delta_p_glycol = None

        self.calc_mass_flows()
        self.calc_lmtd()
        self.calc_ua()
        self.calc_area()
        self.calc_pressure_drops()

    # -- Mass flows from the duty and each stream's temperature change ----------
    def calc_mass_flows(self):
        T_w_avg = (self.T_water_in + self.T_water_out) / 2
        self.cp_water = round(PropsSI('CPMASS', 'T', T_w_avg, 'P', 101325, self.water), 1)
        self.m_dot_water = self.Q / (self.cp_water * (self.T_water_in - self.T_water_out))

        T_g_avg = (self.T_glycol_in + self.T_glycol_out) / 2
        self.cp_glycol = round(PropsSI('CPMASS', 'T', T_g_avg, 'P', 101325, self.glycol), 1)
        self.m_dot_glycol = self.Q / (self.cp_glycol * (self.T_glycol_out - self.T_glycol_in))

    # -- Counter-flow LMTD (both streams sensible, sloped on both sides) --------
    def calc_lmtd(self):
        # Counter-flow: water_in meets glycol_out ; water_out meets glycol_in
        delta_T_1 = self.T_water_in - self.T_glycol_out
        delta_T_2 = self.T_water_out - self.T_glycol_in
        if abs(delta_T_1 - delta_T_2) < 1e-6:
            self.LMTD = delta_T_1
        else:
            self.LMTD = (delta_T_1 - delta_T_2) / np.log(delta_T_1 / delta_T_2)

    def calc_ua(self):
        self.UA = self.Q / self.LMTD

    def calc_area(self):
        self.A_required = self.UA / self.U

    # -- Single-phase friction: MARTIN (1996) generalised Leveque correlation ---
    def _martin_dp(self, m_dot, A_flow, fluid, T_avg):
        """Fanning-type friction factor for a chevron channel (Martin 1996):

            1/sqrt(f) = cos(phi) / sqrt(0.045*tan(phi)+0.09*sin(phi)+f0/cos(phi))
                        + (1-cos(phi)) / sqrt(3.8*f1)

        with, per Re range,
            Re < 2000 : f0 = 16/Re            , f1 = 149/Re + 0.9625
            Re >= 2000: f0 = (1.56*ln Re-3)^-2, f1 = 9.75/Re^0.289
        and dP = 2*f*(L/D_h)*G^2/rho  (Fanning form; identical to the single-
        phase form already used in evaporator.py / condenser.py).

        Angle convention: phi is measured from the longitudinal (flow) axis so
        that a LARGER angle => harder plate => more friction. This is the same
        sense in which Han's beta enters the two-phase correlation, so the
        project's beta = 30 deg is used DIRECTLY (no 90-beta conversion). Flip
        to (90-beta) here if the datasheet defines beta from the horizontal.
        """
        rho = PropsSI('D', 'T', T_avg, 'P', 101325, fluid)
        mu = PropsSI('V', 'T', T_avg, 'P', 101325, fluid)
        G = m_dot / A_flow                       # per-channel mass flux [kg/m^2 s]
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

    def calc_pressure_drops(self):
        T_w_avg = (self.T_water_in + self.T_water_out) / 2
        T_g_avg = (self.T_glycol_in + self.T_glycol_out) / 2
        (self.delta_p_water, self.Re_water, self.f_water,
         self.G_water, self.rho_water, self.mu_water) = \
            self._martin_dp(self.m_dot_water, self.A_flow_water, self.water, T_w_avg)
        (self.delta_p_glycol, self.Re_glycol, self.f_glycol,
         self.G_glycol, self.rho_glycol, self.mu_glycol) = \
            self._martin_dp(self.m_dot_glycol, self.A_flow_glycol, self.glycol, T_g_avg)

    # -- Plate count estimate (A x B outer-footprint proxy, same as report) -----
    def set_plate_count(self, area_per_plate, margin=1.20):
        self.N_plates = int(np.ceil(self.A_required * margin / area_per_plate))
        return self.N_plates


if __name__ == "__main__":
    # ------------------------------------------------------------------ #
    # FREE-COOLING DESIGN POINT (100% load, worst-case OAT for full free  #
    # cooling = 4 C dry-bulb; see ASHRAE 90.1 s6.5.1 economizer sizing).  #
    # ------------------------------------------------------------------ #
    Q = 150e3                              # 100% design load [W]

    # Chilled water: same spec as the evaporator loop
    T_water_in = 21 + 273.15               # CHW return  [K]
    T_water_out = 15 + 273.15              # CHW supply  [K]

    # Glycol from the dry cooler: OAT 4 C + 5 K dry-cooler approach -> 9 C supply,
    # 6 K rise across the economizer -> 15 C return (matches CHW 6 K span).
    T_glycol_in = 9 + 273.15               # glycol supply (cold)  [K]
    T_glycol_out = 15 + 273.15             # glycol return (warm)  [K]

    # Geometry -- same corrugation family as evap/cond, larger single-wall frame
    D_h = 0.004      # [m]  4 mm hydraulic diameter (= 2b/phi, phi=1.17)
    phi = 1.17
    b = D_h * phi / 2        # 2.34 mm channel gap
    Lambda = 0.005   # [m]  5 mm corrugation pitch
    beta = 30        # [deg] chevron angle (soft/medium plate)
    L = 0.70         # [m]  plate flow length (larger frame)
    L_w = 0.25       # [m]  plate width (larger frame)
    N_cp_water = 24
    N_cp_glycol = 24
    U = 3500.0       # [W/m^2 K] assumed overall HTC (glycol/water BPHE)

    econ = Economizer(Q=Q, T_water_in=T_water_in, T_water_out=T_water_out,
                      T_glycol_in=T_glycol_in, T_glycol_out=T_glycol_out,
                      D_h=D_h, L=L, beta=beta, Lambda=Lambda, b=b, L_w=L_w,
                      N_cp_water=N_cp_water, N_cp_glycol=N_cp_glycol, U=U)

    area_per_plate = 0.175   # [m^2] outer A x B footprint of the selected larger plate
    econ.set_plate_count(area_per_plate, margin=1.20)

    print("=== FREE-COOLING ECONOMIZER (glycol -> chilled water, single-phase) ===")
    print(f"Duty                 : {econ.Q/1e3:.1f} kW")
    print(f"Water   21->15 C     : m_dot = {econ.m_dot_water:.2f} kg/s "
          f"({econ.m_dot_water*3600/999:.1f} m3/h), cp = {econ.cp_water} J/kgK")
    print(f"Glycol   9->15 C     : m_dot = {econ.m_dot_glycol:.2f} kg/s "
          f"({econ.m_dot_glycol*3600/econ.rho_glycol:.1f} m3/h), cp = {econ.cp_glycol} J/kgK")
    print(f"LMTD                 : {econ.LMTD:.2f} K")
    print(f"UA                   : {econ.UA/1e3:.2f} kW/K")
    print(f"A_required (U={U:.0f}) : {econ.A_required:.2f} m^2  (bare)")
    print(f"Plate count (x1.2)   : {econ.N_plates} plates @ {area_per_plate} m2/plate")
    print("--- Martin (1996) single-phase pressure drop ---")
    print(f"Water  : G={econ.G_water:5.0f} kg/m2s  Re={econ.Re_water:6.0f}  "
          f"f={econ.f_water:.4f}  dP={econ.delta_p_water/1e3:5.2f} kPa  "
          f"v={econ.G_water/econ.rho_water:.2f} m/s")
    print(f"Glycol : G={econ.G_glycol:5.0f} kg/m2s  Re={econ.Re_glycol:6.0f}  "
          f"f={econ.f_glycol:.4f}  dP={econ.delta_p_glycol/1e3:5.2f} kPa  "
          f"v={econ.G_glycol/econ.rho_glycol:.2f} m/s")
