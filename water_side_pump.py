import numpy as np
from CoolProp.CoolProp import PropsSI

# Digitized from Bell & Gossett Series e-1510 curve sheet (B-261G), 1.5AD
# frame, 1750 RPM composite chart, standard (non "-es") curve. Points read
# directly off the vendor chart in its native units (GPM, ft) -- same
# digitization approach as dry_cooler_pump.py's 2AD/7.0in fit. Module-level
# so other drivers (e.g. main.py's PUE calculation) can reuse it directly.
TRIM_CURVES_1_5AD_1750RPM = {
    "5.0": [(0, 24), (25, 23), (50, 22), (75, 19), (90, 15), (100, 12)],
    "5.5": [(0, 29), (25, 28.5), (50, 27), (75, 24), (100, 18), (110, 13)],
    "6.0": [(0, 34), (25, 33.5), (50, 32), (75, 29), (100, 24), (120, 16)],
    "6.5": [(0, 39), (25, 38.5), (50, 37), (75, 34), (100, 29), (125, 21), (135, 15)],
    "7.0": [(0, 44), (25, 43), (50, 41), (75, 38), (100, 33.5), (125, 26), (150, 15)],
}


class CHWPump:
    """Models the chilled-water (CRAH/data-hall side) circulation pump.

    ROLE IN THE SYSTEM
    -------------------
    This pump circulates chilled water between the mechanical room (evaporator
    in chiller mode / economizer in free-cooling mode) and the data-hall
    CRAH/in-row units, in the SAME fashion the glycol pump (dry_cooler_pump.py)
    circulates propylene glycol between the dry cooler and the
    condenser/economizer. It is a SEPARATE pump on a SEPARATE (chilled-water)
    loop -- it does not touch the glycol side at all.

    Like the glycol loop, the CHW loop has two parallel branches depending on
    mode:
        Mechanical (chiller) mode : CRAH <-> evaporator  (evaporator.py)
        Free-cooling mode         : CRAH <-> economizer   (economizer.py)
    Both branches move the SAME chilled-water flow (5.97 kg/s, 21->15 C) since
    both are sized for the full 150 kW duty; the branch resistances differ
    only in the HX term (evaporator: 7.73 kPa vs. economizer: 6.67 kPa).

    SYSTEM CURVE COMPONENTS (kPa, at design flow)
    -----------------------------------------------------------------
      Evaporator / economizer (worst-case branch = evaporator)   7.73 kPa
      CRAH coil                              -- ASSUMED PLACEHOLDER, 50.0 kPa
                                                 (no CRAH unit is specified by
                                                 the project brief or vendor
                                                 datasheet; typical air-cooled
                                                 CRAH coils run 40-70 kPa at
                                                 this flow -- FLAG pending
                                                 vendor confirmation)
      Piping / fittings / elevation           -- Darcy-Weisbach (Swamee-Jain),
                                                 same method as the glycol
                                                 loop's piping term
                                                 (dry_cooler_pump.py), sized to
                                                 the SAME 65 m equivalent
                                                 length / 8x12 m footprint
                                                 assumption. Project brief says
                                                 to NEGLECT CHW piping losses --
                                                 this term is retained anyway as
                                                 a deliberate OVER-DESIGN margin
                                                 per project decision, but is
                                                 now a real calculation instead
                                                 of a flat placeholder (pass
                                                 dp_piping_kpa explicitly to
                                                 override)
      -----------------------------------------------------------------
      TOTAL (design basis used here)                          ~82.7 kPa

    VENDOR CURVE / TRIM SELECTION METHOD (same method as dry_cooler_pump.py)
    ----------------------------------------------------------------------
    Vendor performance data (Bell & Gossett Series e-1510, document B-261G)
    were digitized from the curve sheet in US customary units (GPM, ft) --
    the native units of the vendor chart -- across several impeller trims,
    and each trim fit to a quadratic H(Q) = c2*Q^2 + c1*Q + c0 (least squares).
    A trim is viable if its 100%-speed curve produces at least a margined
    head (design head + 10%) AT the design flow; the smallest viable trim is
    selected, and a VFD then trims the speed back down (affinity laws) from
    the actual 100%-speed system-curve intersection to hit the exact design
    operating point.

    NOTE ON FRAME SIZE: the glycol pump (dry_cooler_pump.py) uses the 2AD
    frame (7.0in trim) at its much larger duty (~156 GPM @ ~38 ft). At this
    pump's much smaller duty (~95 GPM @ ~28 ft), the 2AD casing would run far
    to the left of its best-efficiency point (~40-50% eta, per the vendor's
    own contour lines) -- an oversized-pump efficiency penalty that works
    directly against the project's PUE target. The smaller 1.5AD frame
    (composite envelope: 66-151 GPM, 16-50 ft, per the B-261G selection
    chart) both brackets this duty point AND keeps the pump running much
    closer to its own BEP, so it is selected here even though it is a
    different frame than the glycol pump.

    INPUTS (from cycle side -- provided by teammates/operating conditions):
      m_dot_water         : CHW mass flow rate [kg/s]  (from evaporator.py /
                             economizer.py -- both agree at 5.97 kg/s)
      dp_hx_kpa           : worst-case HX branch water-side dP [kPa]
                             (evaporator: 7.73 kPa > economizer: 6.67 kPa)
      dp_crah_kpa         : ASSUMED CRAH coil dP [kPa]           (flagged)
      trim_curves_gpm_ft  : dict {trim_label: [(Q_gpm, H_ft), ...]} digitized
                             points per impeller trim, native vendor units
      dp_piping_kpa       : piping/fittings dP [kPa]. If None (default),
                             computed via Darcy-Weisbach from D_pipe/
                             pipe_roughness/L_eq_total; pass a number to
                             override with a different value.
      D_pipe, pipe_roughness, L_eq_total : CHW piping geometry -- ASSUMED,
                             sized for ~1.8 m/s at design flow (same hydronic
                             velocity guideline as the glycol loop), same 65 m
                             equivalent length / 8x12 m footprint assumption.
      speed_rpm           : synchronous speed of the digitized chart [RPM]

    OUTPUTS (for plant selection / PUE):
      selected_trim                    -> selected pump impeller trim
      Q_full_gpm, H_full_ft            -> 100%-speed system-curve intersection
      speed_ratio                      -> VFD speed fraction to hit design duty
      hydraulic_power_w, shaft_power_w -> pump power draw (PUE input)
    """

    # Physical constants
    RHO_WATER = 998.0      # [kg/m^3] at CHW loop average temp (~18 C)
    G = 9.81               # [m/s^2]
    GPM_PER_M3S = 15850.3  # 1 m^3/s = 15850.3 US GPM
    FT_PER_M = 3.28084

    def __init__(self, m_dot_water, dp_hx_kpa, dp_crah_kpa, trim_curves_gpm_ft,
                 dp_piping_kpa=None, D_pipe=0.0627, pipe_roughness=0.045e-3,
                 L_eq_total=65.0, speed_rpm=1750.0, water="Water"):
        self.m_dot_water = m_dot_water
        self.dp_hx_kpa = dp_hx_kpa
        self.dp_crah_kpa = dp_crah_kpa
        self.D_pipe = D_pipe
        self.pipe_roughness = pipe_roughness
        self.L_eq_total = L_eq_total
        self.water = water
        self.speed_rpm = speed_rpm

        # Piping loss: real Darcy-Weisbach/Swamee-Jain calc (same method and
        # 65 m equivalent length as the glycol loop's piping term in
        # dry_cooler_pump.py) -- replaces the earlier flat 25 kPa placeholder.
        self.dp_piping_kpa = dp_piping_kpa if dp_piping_kpa is not None else self._darcy_dp_kpa()

        # ---- Design duty point, SI -> US customary (vendor chart units) ----
        self.dp_total_kpa = dp_hx_kpa + dp_crah_kpa + self.dp_piping_kpa
        self.Q_m3s = m_dot_water / self.RHO_WATER
        self.Q_gpm = self.Q_m3s * self.GPM_PER_M3S
        self.H_m = (self.dp_total_kpa * 1000.0) / (self.RHO_WATER * self.G)
        self.H_ft = self.H_m * self.FT_PER_M

        # ---- System curve constant: H_sys(Q) = k * Q^2 (Q in GPM, H in ft) --
        # Anchored at the design duty point (matches the quadratic-loss
        # assumption used throughout the report's other pressure-drop terms).
        self.k_sys = self.H_ft / self.Q_gpm ** 2

        # ---- Fit each vendor trim to a quadratic, same method as dry_cooler_pump.py ----
        self.trim_fits = {}
        for trim, pts in trim_curves_gpm_ft.items():
            Q = np.array([p[0] for p in pts])
            H = np.array([p[1] for p in pts])
            self.trim_fits[trim] = np.polyfit(Q, H, 2)  # (c2, c1, c0)

        self.selected_trim = None
        self.Q_full_gpm = None
        self.H_full_ft = None
        self.speed_ratio = None

    def _darcy_dp_kpa(self, T_avg=291.15):
        """CHW piping pressure drop via Darcy-Weisbach with the Swamee-Jain
        explicit friction factor -- same method as dry_cooler_pump.py's
        glycol-loop piping term, evaluated at the ~18 C CHW loop average
        temperature (21/15 C design in/out)."""
        rho = PropsSI('D', 'T', T_avg, 'P', 101325, self.water)
        mu = PropsSI('V', 'T', T_avg, 'P', 101325, self.water)
        v = self.m_dot_water / rho / (np.pi / 4 * self.D_pipe ** 2)
        Re = rho * v * self.D_pipe / mu
        f = 0.25 / (np.log10(self.pipe_roughness / (3.7 * self.D_pipe)
                              + 5.74 / Re ** 0.9)) ** 2
        dP = f * (self.L_eq_total / self.D_pipe) * (rho * v ** 2 / 2)
        return dP / 1000.0

    def system_curve(self, Q_gpm):
        """H_sys(Q) = k * Q^2, quadratic loss assumption [ft]."""
        return self.k_sys * Q_gpm ** 2

    def trim_curve(self, trim, Q_gpm):
        """Evaluate a fitted vendor trim's head at a given flow [ft]."""
        c2, c1, c0 = self.trim_fits[trim]
        return c2 * Q_gpm ** 2 + c1 * Q_gpm + c0

    def solve_intersection(self, trim):
        """Solve the 100%-speed intersection of a trim curve with the system
        curve: c2*Q^2 + c1*Q + c0 = k*Q^2  ->  (c2-k)*Q^2 + c1*Q + c0 = 0.
        Returns (Q_gpm, H_ft) of the physically valid (positive-flow) root.
        """
        c2, c1, c0 = self.trim_fits[trim]
        a = c2 - self.k_sys
        b = c1
        c = c0
        disc = b ** 2 - 4 * a * c
        if disc < 0:
            return None
        roots = [(-b - np.sqrt(disc)) / (2 * a), (-b + np.sqrt(disc)) / (2 * a)]
        valid = [r for r in roots if r > 0]
        if not valid:
            return None
        Q_int = min(valid)  # smallest positive root = first crossing
        H_int = self.k_sys * Q_int ** 2
        return Q_int, H_int

    def select_trim(self, margin_head=1.10):
        """Screen all trims directly AT THE DESIGN FLOW: a trim is viable if
        its 100%-speed curve produces at least margin_head x H_design at
        Q_design (the excess head is then trimmed away by the VFD). This is
        the direct physical acceptance test -- equivalent to checking that
        the trim curve lies above the system curve at the design flow, with
        a head margin for fouling/aging/uncertainty.

        The smallest trim that clears is selected, then the actual 100%-speed
        system-curve intersection is solved as the rigorous confirmation
        (this is also the point the VFD trims back FROM).
        """
        H_marg = self.H_ft * margin_head

        results = {}
        for trim in sorted(self.trim_fits, key=float):
            H_at_design_Q = self.trim_curve(trim, self.Q_gpm)
            clears = H_at_design_Q >= H_marg
            results[trim] = (H_at_design_Q, clears)

        candidates = [t for t, (H, ok) in results.items() if ok]
        self.selected_trim = min(candidates, key=float) if candidates else \
            max(results, key=lambda t: results[t][0])

        intersection = self.solve_intersection(self.selected_trim)
        if intersection is not None:
            self.Q_full_gpm, self.H_full_ft = intersection
            self.speed_ratio = self.Q_gpm / self.Q_full_gpm

        return self.selected_trim, results, H_marg

    def hydraulic_power_w(self, Q_gpm=None, H_ft=None):
        """Hydraulic (water) power at a given duty [W]. Defaults to the
        design duty point."""
        Q = self.Q_gpm if Q_gpm is None else Q_gpm
        H = self.H_ft if H_ft is None else H_ft
        Q_m3s = Q / self.GPM_PER_M3S
        H_m = H / self.FT_PER_M
        return self.RHO_WATER * self.G * Q_m3s * H_m

    def shaft_power_w(self, eta_pump, Q_gpm=None, H_ft=None):
        """Shaft (brake) power at a given duty and pump efficiency [W]."""
        return self.hydraulic_power_w(Q_gpm, H_ft) / eta_pump


if __name__ == "__main__":
    # ------------------------------------------------------------------ #
    # DESIGN POINT: CHW flow from evaporator.py / economizer.py (agree at  #
    # 5.97 kg/s, 21->15 C). Worst-case HX branch = evaporator (7.73 kPa). #
    # CRAH coil is still an ASSUMED placeholder, pending vendor data.     #
    # Piping is now a real Darcy-Weisbach calc (was a flat 25 kPa         #
    # placeholder), sized to the same 65 m equivalent length / 8x12 m     #
    # footprint assumption as the glycol loop -- retained per project     #
    # decision even though the brief says to neglect CHW piping losses.  #
    # ------------------------------------------------------------------ #
    m_dot_water = 5.97          # [kg/s]  from evaporator.py / economizer.py
    dp_hx_kpa = 7.73             # [kPa]   evaporator water-side (Martin 1996), worst-case branch
    dp_crah_kpa = 50.0           # [kPa]   ASSUMED -- no vendor CRAH unit specified (FLAG)

    pump = CHWPump(m_dot_water=m_dot_water, dp_hx_kpa=dp_hx_kpa,
                    dp_crah_kpa=dp_crah_kpa,
                    trim_curves_gpm_ft=TRIM_CURVES_1_5AD_1750RPM,
                    speed_rpm=1750.0)

    print("=== CHW (DATA-CENTER / CRAH-SIDE) CIRCULATION PUMP ===")
    print(f"Design flow      : {pump.Q_m3s*3600:.2f} m3/h ({pump.Q_gpm:.1f} GPM), "
          f"m_dot = {m_dot_water:.2f} kg/s")
    print(f"System curve dP  : HX {dp_hx_kpa:.2f} + CRAH {dp_crah_kpa:.1f} (assumed) "
          f"+ piping {pump.dp_piping_kpa:.2f} (Darcy-Weisbach) = {pump.dp_total_kpa:.2f} kPa")
    print(f"Design head      : {pump.H_m:.2f} m ({pump.H_ft:.2f} ft)")
    print(f"System curve k   : {pump.k_sys:.6e}  (H_sys = k*Q_gpm^2, ft)")

    selected, screening, H_marg = pump.select_trim(margin_head=1.10)
    print(f"\nMargined design head (10% on head, at design flow {pump.Q_gpm:.1f} GPM): {H_marg:.2f} ft")
    print("--- Trim screening (1.5AD, 1750 RPM) ---")
    for trim, (H_at_Q, clears) in sorted(screening.items(), key=lambda kv: float(kv[0])):
        print(f"  {trim} in : H@Qdesign = {H_at_Q:5.2f} ft   clears (w/ 10% margin)? {'YES' if clears else 'no'}")
    print(f"\nSelected trim: {selected} in")

    print(f"\n100%-speed (1750 RPM) system-curve intersection: "
          f"{pump.Q_full_gpm:.2f} GPM @ {pump.H_full_ft:.2f} ft")
    print(f"  ({pump.Q_full_gpm/4.4029:.2f} m3/h @ {pump.H_full_ft/3.28084:.2f} m)")
    print(f"VFD speed ratio to reach design duty : {pump.speed_ratio*100:.1f}% "
          f"= {pump.speed_ratio*pump.speed_rpm:.0f} RPM")

    P_hyd = pump.hydraulic_power_w()
    print(f"\nHydraulic power @ design duty : {P_hyd:.1f} W")
    eta_est = 0.60   # digitized from B-261G contour lines near the operating point (~55-63%)
    P_shaft = pump.shaft_power_w(eta_est)
    print(f"Estimated pump efficiency (digitized, off-BEP)  : {eta_est*100:.0f}%")
    print(f"Shaft power @ design duty      : {P_shaft:.1f} W = {P_shaft/745.7:.3f} HP")
    print("  -> nearest standard motor size: 1.5 HP")

    # Low-load / free-cooling turndown check (flow-proportional control assumption)
    Q_min_gpm = pump.Q_gpm * (30.0 / 150.0)
    H_min_ft = pump.system_curve(Q_min_gpm)
    print(f"\nAt 30 kW turndown (20% flow, {Q_min_gpm:.1f} GPM): "
          f"system head = {H_min_ft:.2f} ft -- OPEN ITEM: check CRAH/coil minimum "
          f"flow requirement and pump minimum continuous flow (thrust/overheat risk) "
          f"at this very low head/flow combination.")
