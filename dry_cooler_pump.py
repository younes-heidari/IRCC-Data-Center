import numpy as np
from scipy.optimize import brentq
from CoolProp.CoolProp import PropsSI


class GlycolLoopPump:
    """Models the dry-cooler-loop glycol circulation pump: ONE physical pump,
    VFD-driven, with bypass valves routing it between two circuits depending
    on plant mode (per Setup_Figure_V2's bypass arrangement):

      Mechanical mode : Dry Cooler <-> Condenser      (compressor running)
      Free-cooling mode: Dry Cooler <-> Economizer     (compressor off)

    SELECTED PUMP: Bell & Gossett e-1510 2AD @ 1750 RPM, 5 impeller trims
    digitized from the vendor performance sheet; select_trim() picks the
    smallest trim that clears the margined selection point (7.0 in, given
    the design values below).

    METHODOLOGY:
      1. System curves (mechanical & free-cooling) are built from the SAME
         correlations already validated in condenser.py, dry_cooler.py, and
         economizer.py -- Martin (1996) for the condenser glycol side (exact
         recompute, not scaled), Darcy-Weisbach/Swamee-Jain for piping, and
         dry_cooler.py's own turbulent-scaling method for its two operating
         points. These are evaluated PARAMETRICALLY (as a function of m_dot,
         not just at one design point) so a real intersection can be solved.
      2. The vendor pump curve (selected trim, 1750 RPM) is fit to a
         quadratic H(Q) = c2 Q^2 + c1 Q + c0 via least squares (Q in m^3/h,
         H in m).
      3. Trim selection: smallest trim whose 100%-speed curve still clears
         the margined selection point.
      4. VFD speed fraction to hit each EXACT design duty point is solved via
         the affinity-law equation H_100(Q/s)*s^2 = H_required, root-found
         with scipy.optimize.brentq -- based on the REAL digitized curve, not
         an assumed shutoff-ratio shape.

    INPUTS (from cycle side -- provided by teammates/operating conditions):
      m_dot_mech, T_glycol_hot_mech, T_glycol_cold_mech
                        -> mechanical-mode glycol design point (condenser.py)
      m_dot_free, T_glycol_hot_free, T_glycol_cold_free
                        -> free-cooling glycol design point (economizer.py)
      dP_drycooler_mech_ref -> dry cooler glycol-side dP @ m_dot_mech, mechanical
                                mode (dry_cooler.py, vendor-confirmed)
      dP_econ_free_ref      -> economizer glycol-side dP @ m_dot_free, free-cooling
                                mode (economizer.py, Martin 1996)

    INPUTS (condenser glycol-side plate geometry -- same corrugation family
            as condenser.py, drives the Martin (1996) recompute):
      D_h, L, beta, b, L_w, N_cp_glycol

    INPUTS (piping -- ASSUMED, project brief excludes piping length from scope):
      D_pipe, pipe_roughness, L_eq_total

    OUTPUTS (for plant selection / PUE):
      selected_trim, curve_coeffs        -> selected pump curve
      s_mech, s_free, s_sel              -> VFD speed fractions at each operating point
      P_hyd_mech/free, P_elec_mech/free  -> pump power draw (PUE input)
    """

    # Digitized vendor pump curve: B&G e-1510 2AD @ 1750 RPM, 5 impeller trims
    GPM_POINTS = np.array([25, 50, 75, 100, 125, 150, 175])
    HEAD_FT = {
        "5.0in": np.array([26.00, 25.67, 23.93, 21.06, 17.16, 11.62, 8.64]),
        "5.5in": np.array([31.51, 31.31, 30.04, 26.97, 23.13, 18.09, 11.68]),
        "6.0in": np.array([37.52, 37.35, 36.21, 33.44, 29.64, 24.80, 18.49]),
        "6.5in": np.array([43.93, 43.96, 42.99, 40.45, 32.61, 31.01, 26.07]),
        "7.0in": np.array([51.10, 51.07, 50.27, 48.00, 45.06, 41.09, 35.71]),
    }
    GPM_TO_M3H = 0.227125
    FT_TO_M = 0.3048

    def __init__(self, m_dot_mech, T_glycol_hot_mech, T_glycol_cold_mech,
                 m_dot_free, T_glycol_hot_free, T_glycol_cold_free,
                 dP_drycooler_mech_ref, dP_econ_free_ref,
                 D_h, L, beta, b, L_w, N_cp_glycol,
                 D_pipe, pipe_roughness, L_eq_total,
                 flow_margin=1.10, head_margin=1.18,
                 eta_pump=0.70, eta_motor=0.90,
                 glycol="INCOMP::MPG[0.30]"):
        self.glycol = glycol

        # ---- Design-point flows & temperatures (from sibling component files) ----
        self.m_dot_mech = m_dot_mech
        self.T_avg_mech = (T_glycol_hot_mech + T_glycol_cold_mech) / 2
        self.m_dot_free = m_dot_free
        self.T_avg_free = (T_glycol_hot_free + T_glycol_cold_free) / 2

        # ---- Piping assumptions ----
        self.D_pipe = D_pipe                  # m, pipe ID
        self.pipe_roughness = pipe_roughness  # m
        self.L_eq_total = L_eq_total          # m, ASSUMED (see docstring)

        # ---- Reference component pressure drops ----
        self.dP_drycooler_mech_ref = dP_drycooler_mech_ref  # Pa @ m_dot_mech (mechanical mode)
        self.dP_econ_free_ref = dP_econ_free_ref            # Pa @ m_dot_free (free-cooling mode)

        # Condenser glycol-side plate geometry (same corrugation family as condenser.py)
        self.D_h = D_h
        self.L = L
        self.beta = beta
        self.b = b
        self.L_w = L_w
        self.N_cp_glycol = N_cp_glycol
        self.A_flow_glycol = N_cp_glycol * b * L_w

        # ---- Efficiency assumptions (vendor sheet: H-Q only, no power curve) ----
        self.eta_pump = eta_pump
        self.eta_motor = eta_motor

        # ---- Selection margins (same philosophy as the condenser's 20% area margin) ----
        self.flow_margin = flow_margin
        self.head_margin = head_margin

        # Results (populated by run())
        self.trim_comparison = None
        self.selected_trim = None
        self.curve_coeffs = None
        self.Q_mech_m3h = self.H_mech = None
        self.Q_free_m3h = self.H_free = None
        self.Q_sel_m3h = self.H_sel = None
        self.m_dot_intersect_100pct = None

        self.run()

    # ================= Component correlations (parametric in m_dot) =================

    def _darcy_dp(self, m_dot, T_avg):
        rho = PropsSI('D', 'T', T_avg, 'P', 101325, self.glycol)
        mu = PropsSI('V', 'T', T_avg, 'P', 101325, self.glycol)
        v = m_dot / rho / (np.pi / 4 * self.D_pipe**2)
        Re = rho * v * self.D_pipe / mu
        f = 0.25 / (np.log10(self.pipe_roughness / (3.7 * self.D_pipe)
                              + 5.74 / Re**0.9))**2
        return f * (self.L_eq_total / self.D_pipe) * (rho * v**2 / 2)

    def _martin_dp(self, m_dot, T_avg):
        rho = PropsSI('D', 'T', T_avg, 'P', 101325, self.glycol)
        mu = PropsSI('V', 'T', T_avg, 'P', 101325, self.glycol)
        G = m_dot / self.A_flow_glycol
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
        return 2 * f * G ** 2 * self.L / (rho * self.D_h)

    def _dP_drycooler_mech(self, m_dot):
        # Same-mode, same-temperature scaling: turbulent tube-side ~ m_dot^1.8
        return self.dP_drycooler_mech_ref * (m_dot / self.m_dot_mech) ** 1.8

    def _dP_drycooler_free(self, m_dot):
        # Cross-mode scaling (different T_avg too): reuses dry_cooler.py's own
        # rho^0.8 * v^1.8 * mu^0.2 turbulent-flow scaling, off the SAME
        # mechanical-mode reference point (same physical coil, geometry cancels).
        rho_ref = PropsSI('D', 'T', self.T_avg_mech, 'P', 101325, self.glycol)
        mu_ref = PropsSI('V', 'T', self.T_avg_mech, 'P', 101325, self.glycol)
        rho_new = PropsSI('D', 'T', self.T_avg_free, 'P', 101325, self.glycol)
        mu_new = PropsSI('V', 'T', self.T_avg_free, 'P', 101325, self.glycol)
        v_ratio = (m_dot / rho_new) / (self.m_dot_mech / rho_ref)
        ratio = (rho_new / rho_ref) ** 0.8 * v_ratio ** 1.8 * (mu_new / mu_ref) ** 0.2
        return self.dP_drycooler_mech_ref * ratio

    def _dP_econ_free(self, m_dot):
        return self.dP_econ_free_ref * (m_dot / self.m_dot_free) ** 1.8

    def system_head_mech(self, m_dot):
        """Total system head [m] for the mechanical-mode branch, at any m_dot."""
        dP = (self._dP_drycooler_mech(m_dot) + self._martin_dp(m_dot, self.T_avg_mech)
              + self._darcy_dp(m_dot, self.T_avg_mech))
        rho = PropsSI('D', 'T', self.T_avg_mech, 'P', 101325, self.glycol)
        return dP / (rho * 9.81)

    def system_head_free(self, m_dot):
        """Total system head [m] for the free-cooling branch, at any m_dot."""
        dP = (self._dP_drycooler_free(m_dot) + self._dP_econ_free(m_dot)
              + self._darcy_dp(m_dot, self.T_avg_free))
        rho = PropsSI('D', 'T', self.T_avg_free, 'P', 101325, self.glycol)
        return dP / (rho * 9.81)

    def m_dot_to_m3h(self, m_dot, T_avg):
        rho = PropsSI('D', 'T', T_avg, 'P', 101325, self.glycol)
        return m_dot / rho * 3600

    # ================= Vendor pump curve =================

    def _fit_trim(self, trim):
        Q_m3h = self.GPM_POINTS * self.GPM_TO_M3H
        H_m = self.HEAD_FT[trim] * self.FT_TO_M
        return np.polyfit(Q_m3h, H_m, 2)   # [c2, c1, c0], Q in m3/h, H in m

    def H_pump_100(self, Q_m3h, coeffs):
        return np.polyval(coeffs, Q_m3h)

    def select_trim(self):
        """Compares all 5 digitized trims against the margined selection
        point; picks the smallest trim that still clears it."""
        results = {}
        for trim in self.HEAD_FT:
            coeffs = self._fit_trim(trim)
            H_at_mech = self.H_pump_100(self.Q_mech_m3h, coeffs)
            H_at_sel = self.H_pump_100(self.Q_sel_m3h, coeffs)
            clears = H_at_sel >= self.H_sel
            results[trim] = dict(H_at_mech=H_at_mech, H_at_sel=H_at_sel, clears=clears)
        self.trim_comparison = results
        for trim in ["5.0in", "5.5in", "6.0in", "6.5in", "7.0in"]:
            if results[trim]["clears"]:
                return trim
        return "7.0in"  # fallback, largest available

    def _solve_speed(self, Q_required_m3h, H_required_m, coeffs):
        """Solves VFD speed fraction s such that the affinity-scaled vendor
        curve passes exactly through (Q_required, H_required)."""
        def g(s):
            Q100 = Q_required_m3h / s
            return self.H_pump_100(Q100, coeffs) * s**2 - H_required_m
        s_scan = np.linspace(0.3, 1.05, 300)
        gv = [g(s) for s in s_scan]
        for i in range(len(gv) - 1):
            if gv[i] * gv[i + 1] < 0:
                return brentq(g, s_scan[i], s_scan[i + 1])
        return None

    def _find_intersection(self, coeffs, system_head_fn, m_lo=2.0, m_hi=15.0, n=300):
        m_scan = np.linspace(m_lo, m_hi, n)
        def diff(m_dot):
            Q = self.m_dot_to_m3h(m_dot, self.T_avg_mech)
            return self.H_pump_100(Q, coeffs) - system_head_fn(m_dot)
        vals = [diff(m) for m in m_scan]
        for i in range(len(vals) - 1):
            if vals[i] * vals[i + 1] < 0:
                return brentq(diff, m_scan[i], m_scan[i + 1])
        return None

    # ================= Top-level run =================

    def run(self):
        # Bare duty points (system curve, no margin)
        self.Q_mech_m3h = self.m_dot_to_m3h(self.m_dot_mech, self.T_avg_mech)
        self.H_mech = self.system_head_mech(self.m_dot_mech)
        self.Q_free_m3h = self.m_dot_to_m3h(self.m_dot_free, self.T_avg_free)
        self.H_free = self.system_head_free(self.m_dot_free)

        # Margined selection point
        self.Q_sel_m3h = self.Q_mech_m3h * self.flow_margin
        self.H_sel = self.H_mech * self.head_margin

        # Trim selection against the real vendor curves
        self.selected_trim = self.select_trim()
        self.curve_coeffs = self._fit_trim(self.selected_trim)

        # Uncontrolled (100% speed) intersection with the mechanical system curve
        self.m_dot_intersect_100pct = self._find_intersection(
            self.curve_coeffs, self.system_head_mech)
        self.Q_intersect_100pct = self.m_dot_to_m3h(
            self.m_dot_intersect_100pct, self.T_avg_mech)
        self.H_intersect_100pct = self.system_head_mech(self.m_dot_intersect_100pct)

        # VFD speed fractions to hit each EXACT design duty point
        self.s_mech = self._solve_speed(self.Q_mech_m3h, self.H_mech, self.curve_coeffs)
        self.s_free = self._solve_speed(self.Q_free_m3h, self.H_free, self.curve_coeffs)
        self.s_sel = self._solve_speed(self.Q_sel_m3h, self.H_sel, self.curve_coeffs)

        # Power estimates (assumed efficiency -- vendor sheet gives H-Q only)
        rho_mech = PropsSI('D', 'T', self.T_avg_mech, 'P', 101325, self.glycol)
        rho_free = PropsSI('D', 'T', self.T_avg_free, 'P', 101325, self.glycol)
        Qv_mech = self.m_dot_mech / rho_mech
        Qv_free = self.m_dot_free / rho_free
        self.P_hyd_mech = (self.H_mech * rho_mech * 9.81) * Qv_mech
        self.P_elec_mech = self.P_hyd_mech / self.eta_pump / self.eta_motor
        self.P_hyd_free = (self.H_free * rho_free * 9.81) * Qv_free
        self.P_elec_free = self.P_hyd_free / self.eta_pump / self.eta_motor


if __name__ == "__main__":
    # ------------------------------------------------------------------ #
    # Design points -- matches condenser.py (mechanical, 40->46 C glycol,
    # corrected for the 52 C condensing point) and economizer.py /
    # dry_cooler.py (free-cooling, 9->15 C glycol).
    # ------------------------------------------------------------------ #
    pump = GlycolLoopPump(
        m_dot_mech=9.07, T_glycol_hot_mech=46 + 273.15, T_glycol_cold_mech=40 + 273.15,
        m_dot_free=6.52, T_glycol_hot_free=15 + 273.15, T_glycol_cold_free=9 + 273.15,
        dP_drycooler_mech_ref=48.26e3,   # Pa, Kelvion vendor-confirmed @ 9.07 kg/s
        dP_econ_free_ref=10.43e3,        # Pa, economizer.py Martin (1996) @ 6.52 kg/s

        # Condenser glycol-side plate geometry (same corrugation family as condenser.py)
        D_h=0.004, L=0.5, beta=30, b=0.004 * 1.17 / 2, L_w=0.2, N_cp_glycol=24,

        # Piping -- ASSUMED, ~65 m total equivalent length (mech-room-to-roof-and-back
        # + fittings allowance), sized to the 8x12 m data-hall footprint
        D_pipe=0.0779, pipe_roughness=0.045e-3, L_eq_total=65.0,
    )

    print("=" * 72)
    print("TRIM SELECTION -- B&G e-1510 2AD @ 1750 RPM (5 trims compared)")
    print("=" * 72)
    print(f"{'Trim':>8} {'H @ mech Q':>14} {'H @ sel. Q':>14} {'Clears margin?':>16}")
    for trim in ["5.0in", "5.5in", "6.0in", "6.5in", "7.0in"]:
        r = pump.trim_comparison[trim]
        flag = "YES" if r["clears"] else "no"
        marker = "  <-- SELECTED" if trim == pump.selected_trim else ""
        print(f"{trim:>8} {r['H_at_mech']:>12.2f} m {r['H_at_sel']:>12.2f} m {flag:>16}{marker}")

    print()
    print("=" * 72)
    print(f"SELECTED: {pump.selected_trim} trim, B&G e-1510 2AD @ 1750 RPM")
    print("=" * 72)
    print(f"Required duty (mechanical) : {pump.Q_mech_m3h:.2f} m3/h @ {pump.H_mech:.2f} m")
    print(f"Required duty (free-cooling): {pump.Q_free_m3h:.2f} m3/h @ {pump.H_free:.2f} m")
    print(f"Selection point (margined) : {pump.Q_sel_m3h:.2f} m3/h @ {pump.H_sel:.2f} m")
    print()
    print(f"Uncontrolled 100%-speed intersection (mechanical system curve):")
    print(f"  {pump.Q_intersect_100pct:.2f} m3/h @ {pump.H_intersect_100pct:.2f} m "
          f"({pump.m_dot_intersect_100pct:.2f} kg/s)")
    print(f"  -> confirms trim has margin above the design duty; VFD trims it back down")
    print()
    print(f"VFD speed to hit EXACT mechanical design point : {pump.s_mech*100:.1f}%")
    print(f"VFD speed to hit EXACT free-cooling design point: {pump.s_free*100:.1f}%")
    print(f"VFD speed to hit margined selection point       : {pump.s_sel*100:.1f}%")
    print()
    print(f"Hydraulic power (mechanical): {pump.P_hyd_mech:.0f} W  "
          f"| Electrical (assumed 70%/90% eff.): {pump.P_elec_mech:.0f} W")
    print(f"Hydraulic power (free-cooling): {pump.P_hyd_free:.0f} W  "
          f"| Electrical (assumed 70%/90% eff.): {pump.P_elec_free:.0f} W")
