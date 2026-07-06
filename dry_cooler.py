import numpy as np
from scipy.optimize import brentq
from CoolProp.CoolProp import PropsSI


class DryCooler:
    """Models the outdoor air-cooled dry cooler: rejects heat from the
    intermediate propylene-glycol (PG) loop to ambient air through a
    finned-tube coil with variable-speed axial fans. This is the glycol
    loop's final heat-rejection step -- glycol enters hot (from the
    refrigerant condenser in mechanical mode, or straight from the
    economizer in free-cooling mode) and leaves cooled.

    WHY A DIFFERENT METHOD THAN evaporator.py / condenser.py / economizer.py
    --------------------------------------------------------------------
    Those three are plate heat exchangers: single-pass counter-flow, with
    refrigerant-side (Han, Lee & Kim 2003) or single-phase liquid-side
    (Martin 1996) chevron-plate friction correlations. The dry cooler is a
    CROSS-FLOW, finned-tube, gas/liquid exchanger -- air crosses the tube
    bundle perpendicular to the glycol flow over several rows, not
    counter-current -- so a single LMTD is not appropriate, and neither
    plate correlation applies to air over finned round tubes. Instead this
    model uses:
      - Effectiveness-NTU, cross-flow with both fluids unmixed (Incropera),
        for the thermal duty.
      - Fan affinity laws (V_dot ~ N, dP ~ N^2, P_fan ~ N^3) for the
        variable-speed-fan power, which is what actually shows up in PUE.
      - Darcy-Weisbach (circular tube) for the glycol-side pressure drop.
      - Air-side pressure drop is NOT modeled here -- finned-bundle air
        friction is proprietary fin-geometry data normally taken straight
        from the manufacturer's selection software (Kelvion, Guntner);
        treat it as an open item pending a vendor selection, same as the
        BPHE geometry parameters are flagged ASSUMED elsewhere.

    TWO USES OF THE SAME COIL (one UA, two operating points)
    ---------------------------------------------------------
    1) size_design_point(): at the MECHANICAL-MODE design condition
       (max duty, 35 C DB ambient) choose an air-side temperature rise and
       back out the UA / air flow / fan size the coil needs -- this mirrors
       how a dry cooler is actually selected against a vendor catalog.
    2) predict_off_design(): holding that UA (scaled for reduced air flow),
       solve FORWARD for the fan speed needed to hit a target duty at a
       different ambient temperature and glycol inlet temp. Used to check
       the free-cooling condition (150 kW duty, 4 C DB, glycol 15->9 C)
       required by ASHRAE 90.1 S6.5.1, and to build the IPLV/PUE part-load
       fan-power curve.

    INPUTS (design point):
      Q_design           : design heat-rejection duty [W]        (~190 kW, = Q_cond)
      T_glycol_hot_in     : glycol entering (hot), design [K]     (from condenser, 46 C)
      T_glycol_cold_out   : glycol leaving (cold), design target [K] (my choice, 40 C)
      T_air_in_design     : design ambient dry-bulb [K]           (35 C, project spec)
      dT_air_design       : chosen air-side temperature rise [K]  (my choice, ~8-12 K)
      glycol              : secondary-fluid name [str]            (CoolProp INCOMP brine)

    OUTPUTS (from size_design_point):
      UA_design, NTU_design, eps_design, Cr_design
      m_dot_air_design, m_dot_glycol_design   -> fan / pump sizing

    OUTPUTS (from predict_off_design):
      m_dot_air_req, fan_speed_frac, fan_power    -> VSD fan operating point
      T_glycol_cold_out_actual                    -> achieved leaving glycol temp
    """

    def __init__(self, Q_design, T_glycol_hot_in, T_glycol_cold_out,
                 T_air_in_design, dT_air_design,
                 glycol="INCOMP::MPG[0.30]"):
        self.glycol = glycol

        # Design-point duty and terminal temperatures
        self.Q_design = Q_design
        self.T_glycol_hot_in = T_glycol_hot_in
        self.T_glycol_cold_out = T_glycol_cold_out
        self.T_air_in_design = T_air_in_design
        self.dT_air_design = dT_air_design
        self.T_air_out_design = T_air_in_design + dT_air_design

        # Populated by size_design_point()
        self.m_dot_glycol_design = None
        self.cp_glycol_design = None
        self.m_dot_air_design = None
        self.cp_air_design = None
        self.C_glycol_design = None
        self.C_air_design = None
        self.Cmin_design = None
        self.Cmax_design = None
        self.Cr_design = None
        self.eps_design = None
        self.NTU_design = None
        self.UA_design = None
        self.air_side_UA_exponent = 0.6   # h_air ~ (m_dot_air)^0.6 (typical finned-tube Re-power law)

        self.size_design_point()

    # ------------------------------------------------------------------ #
    # Effectiveness-NTU relation, cross-flow, both fluids unmixed        #
    # (Incropera, Fundamentals of Heat and Mass Transfer)                #
    # ------------------------------------------------------------------ #
    @staticmethod
    def _eps_crossflow(NTU, Cr):
        if Cr < 1e-6:
            return 1 - np.exp(-NTU)
        return 1 - np.exp((1.0 / Cr) * NTU**0.22 * (np.exp(-Cr * NTU**0.78) - 1.0))

    # ------------------------------------------------------------------ #
    # DESIGN POINT: choose dT_air_design -> get m_dot_air, then UA        #
    # ------------------------------------------------------------------ #
    def size_design_point(self):
        # Glycol side: mass flow fixed by the duty and the chosen glycol dT
        T_gly_avg = (self.T_glycol_hot_in + self.T_glycol_cold_out) / 2
        self.cp_glycol_design = round(PropsSI('CPMASS', 'T', T_gly_avg, 'P', 101325, self.glycol), 1)
        self.m_dot_glycol_design = self.Q_design / (
            self.cp_glycol_design * (self.T_glycol_hot_in - self.T_glycol_cold_out))

        # Air side: mass flow fixed by the duty and the CHOSEN air-side rise
        T_air_avg = (self.T_air_in_design + self.T_air_out_design) / 2
        self.cp_air_design = PropsSI('CPMASS', 'T', T_air_avg, 'P', 101325, 'Air')
        self.m_dot_air_design = self.Q_design / (self.cp_air_design * self.dT_air_design)

        # Capacity rates
        self.C_glycol_design = self.m_dot_glycol_design * self.cp_glycol_design
        self.C_air_design = self.m_dot_air_design * self.cp_air_design
        self.Cmin_design = min(self.C_glycol_design, self.C_air_design)
        self.Cmax_design = max(self.C_glycol_design, self.C_air_design)
        self.Cr_design = self.Cmin_design / self.Cmax_design

        # Required effectiveness from the actual terminal temperatures
        self.eps_design = self.Q_design / (self.Cmin_design * (self.T_glycol_hot_in - self.T_air_in_design))
        if not (0 < self.eps_design < 1):
            raise ValueError(f"eps_design = {self.eps_design:.3f} is outside (0,1) -- "
                              "check design temperatures / dT_air_design choice.")

        # Invert eps-NTU for NTU (transcendental -> root-find)
        self.NTU_design = brentq(lambda NTU: self._eps_crossflow(NTU, self.Cr_design) - self.eps_design,
                                  1e-4, 50.0)
        self.UA_design = self.NTU_design * self.Cmin_design

    # ------------------------------------------------------------------ #
    # OFF-DESIGN: same coil (UA scaled for air flow), solve for the       #
    # air flow (fan speed) needed to hit a target duty at new T_air_in /  #
    # T_glycol_hot_in. Glycol flow held at the design value (fixed-speed  #
    # glycol pump assumption -- flag if a VSD glycol pump is used).       #
    # ------------------------------------------------------------------ #
    def _UA_at_air_flow(self, m_dot_air):
        # Air-side resistance assumed to dominate (typical for air-cooled
        # coils, h_air << h_liquid), so UA scales with air flow via the
        # usual Re-power law for finned tube banks; glycol-side resistance
        # held constant. ASSUMED simplification -- flag for refinement once
        # a real air-side/glycol-side UA split is available from the vendor.
        return self.UA_design * (m_dot_air / self.m_dot_air_design) ** self.air_side_UA_exponent

    def predict_off_design(self, Q_target, T_glycol_hot_in, T_air_in,
                            m_dot_air_bounds=None):
        m_dot_glycol = self.m_dot_glycol_design  # fixed-speed glycol pump assumption
        cp_glycol = PropsSI('CPMASS', 'T', T_glycol_hot_in, 'P', 101325, self.glycol)
        C_glycol = m_dot_glycol * cp_glycol

        def Q_predicted(m_dot_air):
            cp_air = PropsSI('CPMASS', 'T', T_air_in + 5, 'P', 101325, 'Air')  # approx avg
            C_air = m_dot_air * cp_air
            Cmin, Cmax = min(C_glycol, C_air), max(C_glycol, C_air)
            Cr = Cmin / Cmax
            UA = self._UA_at_air_flow(m_dot_air)
            NTU = UA / Cmin
            eps = self._eps_crossflow(NTU, Cr)
            return eps * Cmin * (T_glycol_hot_in - T_air_in)

        if m_dot_air_bounds is None:
            m_dot_air_bounds = (0.05 * self.m_dot_air_design, 1.5 * self.m_dot_air_design)

        m_dot_air_req = brentq(lambda m: Q_predicted(m) - Q_target, *m_dot_air_bounds)

        fan_speed_frac = m_dot_air_req / self.m_dot_air_design   # V_dot ~ N (approx, const air density)
        T_glycol_cold_out_actual = T_glycol_hot_in - Q_target / C_glycol

        return {
            "m_dot_air_req": m_dot_air_req,
            "fan_speed_frac": fan_speed_frac,
            "T_glycol_cold_out_actual": T_glycol_cold_out_actual,
        }

    # ------------------------------------------------------------------ #
    # Fan power: affinity law, P ~ N^3, off a rated design fan power      #
    # ------------------------------------------------------------------ #
    @staticmethod
    def fan_power(fan_speed_frac, P_fan_rated):
        return P_fan_rated * fan_speed_frac ** 3

    # ------------------------------------------------------------------ #
    # Glycol-side pressure drop, SCALED from a real vendor-reported value #
    # at a reference condition. For the SAME physical coil, Darcy-Weisbach #
    # gives dP ~ rho^0.8 * v^1.8 * mu^0.2 (Blasius, turbulent); since v ~  #
    # m_dot/rho for a fixed flow area, the unknown tube geometry (D, L,    #
    # N_circuits) cancels entirely in the ratio between two conditions on  #
    # the SAME coil. This lets a single real vendor data point (e.g. the   #
    # Kelvion-confirmed design-point dP) be scaled to a different          #
    # operating condition WITHOUT guessing tube length/diameter -- unlike  #
    # glycol_pressure_drop() below, which needs that geometry explicitly   #
    # and is only as good as the guessed L_tube.                          #
    # ------------------------------------------------------------------ #
    def scale_pressure_drop(self, dP_ref, T_ref_avg, m_dot_ref, T_new_avg, m_dot_new):
        rho_ref = PropsSI('D', 'T', T_ref_avg, 'P', 101325, self.glycol)
        mu_ref = PropsSI('V', 'T', T_ref_avg, 'P', 101325, self.glycol)
        rho_new = PropsSI('D', 'T', T_new_avg, 'P', 101325, self.glycol)
        mu_new = PropsSI('V', 'T', T_new_avg, 'P', 101325, self.glycol)

        v_ratio = (m_dot_new * rho_ref) / (m_dot_ref * rho_new)   # v ~ m_dot/rho, flow area cancels
        dP_ratio = (rho_new / rho_ref) ** 0.8 * v_ratio ** 1.8 * (mu_new / mu_ref) ** 0.2
        dP_new = dP_ref * dP_ratio
        return dP_new, dP_ratio

    # ------------------------------------------------------------------ #
    # Glycol-side pressure drop: single-phase, circular tube,             #
    # Darcy-Weisbach with Blasius/laminar friction factor (same style as   #
    # the single-phase zones in condenser.py, just a round-tube D_h).     #
    # NOTE: needs an assumed tube run length L_tube -- prefer              #
    # scale_pressure_drop() above whenever a real reference dP exists,    #
    # since that method doesn't depend on the unknown geometry at all.    #
    # ------------------------------------------------------------------ #
    def glycol_pressure_drop(self, m_dot_glycol, D_tube, L_tube, N_circuits, T_avg):
        rho = PropsSI('D', 'T', T_avg, 'P', 101325, self.glycol)
        mu = PropsSI('V', 'T', T_avg, 'P', 101325, self.glycol)
        A_tube = np.pi / 4 * D_tube ** 2
        m_dot_per_circuit = m_dot_glycol / N_circuits
        v = m_dot_per_circuit / (rho * A_tube)
        Re = rho * v * D_tube / mu
        f_darcy = 64 / Re if Re < 2300 else 0.184 * Re ** (-0.2)   # Blasius, smooth tube
        delta_p = f_darcy * (L_tube / D_tube) * 0.5 * rho * v ** 2
        return delta_p, Re, v


if __name__ == "__main__":
    # ================================================================== #
    # MECHANICAL-MODE DESIGN POINT -- REAL VENDOR SELECTION               #
    # Kelvion ULF-PA104Y4V-096Z100 (Kelvion Select RT, Jul 2026)          #
    #   Capacity  213.23 kW  (12.2% margin over the 190 kW requirement -- #
    #             catalog step; nearest available "Z100" speed tap)      #
    #   Glycol    46.0 -> 40.0 C, 30% propylene glycol, dP = 48.26 kPa    #
    #   Ambient   35.0 C DB, 213 m altitude (Champaign)                  #
    #   Air       outlet 40.56 C -> REAL air-side rise = 5.56 K          #
    #             (replaces the earlier assumed 10 K)                    #
    #   Coil      1x2x34 (68 total circuits), 4 passes, 938.5 m2 surface #
    #   Fans      4x EC, 37.8 in dia; selected op. point 1000 rpm,       #
    #             2.685 kW/fan = 10.74 kW total (nameplate max 1100 rpm, #
    #             17.60 kW total) -- replaces the earlier ASSUMED value  #
    # ================================================================== #
    dc = DryCooler(
        Q_design=213.23e3,
        T_glycol_hot_in=46 + 273.15,
        T_glycol_cold_out=40 + 273.15,
        T_air_in_design=35 + 273.15,
        dT_air_design=5.56,          # REAL, from vendor air outlet temp (was ASSUMED 10 K)
    )

    print("=== DRY COOLER -- MECHANICAL-MODE DESIGN POINT (Kelvion-validated) ===")
    print(f"Glycol flow   : {dc.m_dot_glycol_design:.2f} kg/s")
    print(f"Air flow      : {dc.m_dot_air_design:.2f} kg/s")
    print(f"Cmin / Cmax   : {dc.Cmin_design/1e3:.2f} / {dc.Cmax_design/1e3:.2f} kW/K  (Cr = {dc.Cr_design:.3f})")
    print(f"eps / NTU     : {dc.eps_design:.3f} / {dc.NTU_design:.3f}")
    print(f"UA_design     : {dc.UA_design/1e3:.2f} kW/K")

    # --- Validation against the real Kelvion datasheet ---
    rho_air_approx = 1.105   # kg/m3 at ~38 C avg, 213 m altitude
    rho_glycol_approx = 1030  # kg/m3, 30% PG at ~43 C avg
    m_dot_air_m3h = dc.m_dot_air_design / rho_air_approx * 3600
    m_dot_glycol_Ls = dc.m_dot_glycol_design / rho_glycol_approx * 1000
    print("\n--- Model vs. vendor datasheet (validation) ---")
    print(f"Air flow    : model {m_dot_air_m3h:,.0f} m3/h  vs  vendor 122,071 m3/h  "
          f"({(m_dot_air_m3h/122071-1)*100:+.1f}%)")
    print(f"Glycol flow : model {m_dot_glycol_Ls:.2f} L/s  vs  vendor 8.884 L/s  "
          f"({(m_dot_glycol_Ls/8.884-1)*100:+.1f}%)")

    # REAL fan power at the design-point operating speed (1000 rpm, per datasheet)
    P_fan_rated = 10.74e3  # W, REAL -- Kelvion selected-speed total (4 fans, was ASSUMED 4 kW)
    print(f"\nFan power @ design point (vendor-reported): {P_fan_rated/1e3:.2f} kW total "
          f"(4 fans @ 1000 rpm; nameplate max 17.60 kW @ 1100 rpm)")

    # ================================================================== #
    # FREE-COOLING CHECK (matches economizer.py design point):            #
    # Q = 150 kW, ambient 4 C DB, glycol 15 C (hot, return) -> target      #
    # 9 C (cold, supply). NOTE: Kelvion Select RT was only run at the     #
    # mechanical-mode point above -- this off-design prediction is still  #
    # the MODEL's estimate, not vendor-confirmed. Re-run the selection at #
    # this condition (or use its performance-check mode) to validate it,  #
    # same as was done for the mechanical-mode point.                     #
    # ================================================================== #
    result = dc.predict_off_design(
        Q_target=150e3,
        T_glycol_hot_in=15 + 273.15,
        T_air_in=4 + 273.15,
    )
    fan_pwr = dc.fan_power(result["fan_speed_frac"], P_fan_rated)

    print("\n=== FREE-COOLING CHECK (150 kW @ 4 C DB) -- model estimate, not yet vendor-confirmed ===")
    print(f"Air flow required   : {result['m_dot_air_req']:.2f} kg/s "
          f"({result['fan_speed_frac']*100:.0f}% of design flow)")
    print(f"Glycol leaving temp : {result['T_glycol_cold_out_actual']-273.15:.2f} C  "
          f"(target 9.0 C)")
    print(f"Fan power @ this speed : {fan_pwr/1e3:.2f} kW  "
          f"(vs. {P_fan_rated/1e3:.2f} kW rated -- affinity-law saving)")

    # Glycol-side pressure drop -- design point is REAL (Kelvion, 48.26 kPa).
    # Free-cooling point is SCALED from that real value using
    # scale_pressure_drop() -- geometry-independent, so no guessed tube
    # length is needed (unlike the earlier glycol_pressure_drop() attempt,
    # which was off by 75% because L_tube had to be guessed).
    dP_design_real = 48.26e3  # Pa, Kelvion datasheet (7 psi)
    T_design_avg = (dc.T_glycol_hot_in + dc.T_glycol_cold_out) / 2
    T_free_avg = (15 + 273.15 + 9 + 273.15) / 2

    # Free-cooling glycol mass flow (fixed-speed pump assumption, same as
    # predict_off_design(): C_glycol computed there from Q_target/dT)
    cp_glycol_free = PropsSI('CPMASS', 'T', T_free_avg, 'P', 101325, dc.glycol)
    m_dot_glycol_free = 150e3 / (cp_glycol_free * (15 - 9))

    dP_free_scaled, dP_ratio = dc.scale_pressure_drop(
        dP_ref=dP_design_real, T_ref_avg=T_design_avg, m_dot_ref=dc.m_dot_glycol_design,
        T_new_avg=T_free_avg, m_dot_new=m_dot_glycol_free,
    )
    print(f"\nGlycol-side dP, design point (Kelvion-confirmed) : {dP_design_real/1e3:.2f} kPa")
    print(f"Glycol-side dP, free-cooling (scaled from above)  : {dP_free_scaled/1e3:.2f} kPa  "
          f"(ratio {dP_ratio:.3f}; m_dot drops {dc.m_dot_glycol_design:.2f}->{m_dot_glycol_free:.2f} kg/s, "
          f"but viscosity rises ~2.8x at 12C avg -- net effect is still lower dP than design point)")
    print("(Uses the SAME physical coil's geometry implicitly -- D_tube, L_tube, N_circuits\n"
          " all cancel in the ratio, so no guessed tube length is needed for this estimate.\n"
          " Caveat: Re at the free-cooling point is close to the laminar/turbulent transition,\n"
          " so this is a reasonable extrapolation of the Blasius correlation, not a rigorous one.)")

### Remarks:
### - Mechanical-mode design point is now VENDOR-CONFIRMED: Kelvion ULF-PA104Y4V-096Z100
###   (Kelvion Select RT), selected at Q=213.23 kW (12.2% margin over the 190 kW
###   requirement), 30% PG, 46/40 C glycol, 35 C DB ambient, 213 m altitude.
### - Model validation: feeding the REAL duty and REAL air-side rise (5.56 K, not
###   the earlier assumed 10 K) into size_design_point() predicts air flow within
###   +1.7% and glycol flow within -0.9% of the vendor's reported values -- a
###   good check on the effectiveness-NTU cross-flow approach used here.
### - Fan power at the design point is now REAL (10.74 kW total, 4x EC fans at
###   their selected 1000 rpm operating point; nameplate max 17.60 kW at 1100 rpm),
###   replacing the earlier ASSUMED 4-14 kW guesses.
### - Glycol-side pressure drop: design point uses the REAL Kelvion value
###   (48.26 kPa). Free-cooling point (~32 kPa) is SCALED from that real
###   value via scale_pressure_drop() -- a geometry-independent Darcy-Weisbach
###   ratio (dP ~ rho^0.8 * v^1.8 * mu^0.2) that cancels the unknown tube
###   length/diameter entirely, so it doesn't repeat the earlier 75% guess
###   error. Caveat: Re at the free-cooling point (~2,400) sits close to the
###   laminar/turbulent transition, so treat this as a reasonable
###   extrapolation, not a rigorous result -- flag it as such in the report.
### - glycol_pressure_drop() (guessed-geometry Darcy-Weisbach) is kept as a
###   utility method for cases with NO real reference point to scale from,
###   but is no longer used in the main demo now that a real anchor exists.
### - The free-cooling THERMAL check (150 kW @ 4 C fan speed/power) is still
###   the MODEL's own estimate -- Kelvion Select RT was only run at the
###   mechanical-mode point. Re-run the selection at that condition to get a
###   vendor-confirmed number for the thermal side too.
### - air_side_UA_exponent (0.6) remains a modeling assumption for off-design
###   UA scaling; not addressed by this datasheet.
