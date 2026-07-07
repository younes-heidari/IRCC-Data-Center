import numpy as np
from scipy.optimize import brentq
from CoolProp.CoolProp import PropsSI


# ============================================================================
# AHRI STANDARD 540 / EN12900 10-COEFFICIENT COMPRESSOR MAP
# ============================================================================
#
#   y = c1 + c2*to + c3*tc + c4*to^2 + c5*to*tc + c6*tc^2
#         + c7*to^3 + c8*tc*to^2 + c9*to*tc^2 + c10*tc^3
#
#   to = evaporating dew point [degC], tc = condensing dew point [degC]
#
# Vendor: Bitzer 4FEP-35Z (semi-hermetic reciprocating), refrigerant R-290.
# Source: R290_System_Modeling_Parameters.txt (project data, this chapter).
#
# IMPORTANT -- RATING BASIS: these coefficients reproduce compressor
# performance at EN12900 STANDARD rating conditions: suction (return) gas
# temperature fixed at 20 degC ABSOLUTE (not a fixed superheat -- the actual
# superheat this represents changes with `to`), and 0 K liquid subcooling
# (saturated liquid at tc feeding the expansion valve). These are NOT the
# actual application conditions and must be corrected -- see
# `Compressor.performance()` below.
# ============================================================================

COEFF_Q_W = (114372.76126639, 3903.60053484, -1121.94016238, 45.36731498,
             -33.97146035, 3.40438432, 0.16833364, -0.33798336,
             0.00989812, -0.02557815)          # Cooling capacity [W]

COEFF_P_W = (7147.23542138, -239.01196212, 255.72496495, -7.86795916,
             8.64966259, 1.61953582, -0.04119267, 0.06978450,
             0.01025775, -0.02036026)          # Power input [W]

COEFF_M_KGH = (1009.98550735, 36.10432058, -4.33944694, 0.49390814,
               -0.09966336, 0.03776835, 0.00286517, -0.00108556,
               -0.00002386, -0.00038152)       # Mass flow [kg/h]

EN12900_RETURN_GAS_TEMP_C = 20.0   # fixed ABSOLUTE suction gas temp, per EN12900
EN12900_SUBCOOLING_K = 0.0         # saturated liquid feed, per EN12900


def _ahri540_poly(to, tc, c):
    """Evaluate the AHRI-540 / EN12900 10-coefficient polynomial."""
    c1, c2, c3, c4, c5, c6, c7, c8, c9, c10 = c
    return (c1 + c2 * to + c3 * tc + c4 * to**2 + c5 * to * tc + c6 * tc**2
            + c7 * to**3 + c8 * tc * to**2 + c9 * to * tc**2 + c10 * tc**3)


class Compressor:
    """Single Bitzer 4FEP-35Z semi-hermetic reciprocating compressor,
    R-290, modeled from its AHRI-540/EN12900 10-coefficient performance map.

    Two performance levels are available:
      * `catalog_performance()`  -- raw polynomial output, AT THE EN12900
        RATING BASIS (20 degC absolute return gas, 0 K subcooling). This is
        what the manufacturer's map coefficients directly reproduce.
      * `performance()`          -- CORRECTED to the actual application
        suction superheat and liquid subcooling, via a suction-density-ratio
        correction (standard industry practice for adapting a fixed rating
        condition to actual operating superheat/subcooling). This is the one
        that should be used everywhere else in the system model.
    """

    REFRIGERANT = "R290"

    def __init__(self, refrigerant="R290"):
        self.R = refrigerant

    # ------------------------------------------------------------------ #
    def catalog_performance(self, to_C, tc_C):
        """Raw AHRI-540 map output at the EN12900 rating basis (20 degC
        absolute return gas, 0 K subcooling). Returns dict with Q [W],
        P [W], m_dot [kg/h]."""
        Q_w = _ahri540_poly(to_C, tc_C, COEFF_Q_W)
        P_w = _ahri540_poly(to_C, tc_C, COEFF_P_W)
        m_kgh = _ahri540_poly(to_C, tc_C, COEFF_M_KGH)
        return {"Q_w": Q_w, "P_w": P_w, "m_dot_kgh": m_kgh}

    # ------------------------------------------------------------------ #
    def _liquid_enthalpy(self, tc_C, subcooling_K, p_cond):
        """Enthalpy of the liquid feeding the expansion valve, at a given
        subcooling below tc_C. Uses quality-based lookup at exactly 0 K
        subcooling to avoid a CoolProp ambiguity at the saturation point."""
        if subcooling_K <= 1e-6:
            return PropsSI("HMASS", "T", tc_C + 273.15, "Q", 0, self.R)
        return PropsSI("HMASS", "T", tc_C - subcooling_K + 273.15, "P", p_cond, self.R)

    # ------------------------------------------------------------------ #
    def performance(self, to_C, tc_C, superheat_K=8.0, subcooling_K=0.0):
        """Corrected performance at the ACTUAL application suction
        superheat and liquid subcooling.

        Correction method (density-ratio, standard industry practice):
          1. Mass flow scales with the ratio of actual-to-rated suction gas
             density (same physical displacement + volumetric efficiency,
             so a denser suction gas simply packs more mass per stroke).
          2. Capacity is recomputed directly from the corrected mass flow
             and the ACTUAL refrigeration effect (h_suction - h_liquid),
             which is the more rigorous option -- it also captures the
             (usually small) enthalpy shift from a different suction
             temperature and/or subcooling, not just the density change.
          3. Power input is scaled by the same density ratio as mass flow
             (specific compression work is only weakly sensitive to suction
             superheat over small ranges; this is the standard simplified
             correction used in manufacturer application guides).
        """
        cat = self.catalog_performance(to_C, tc_C)

        p_evap = PropsSI("P", "T", to_C + 273.15, "Q", 1, self.R)
        p_cond = PropsSI("P", "T", tc_C + 273.15, "Q", 1, self.R)

        # --- Rated (EN12900) suction & liquid states ---------------------
        T_suction_rated = EN12900_RETURN_GAS_TEMP_C + 273.15
        rho_rated = PropsSI("D", "T", T_suction_rated, "P", p_evap, self.R)
        h_suction_rated = PropsSI("HMASS", "T", T_suction_rated, "P", p_evap, self.R)
        h_liquid_rated = self._liquid_enthalpy(tc_C, EN12900_SUBCOOLING_K, p_cond)

        # --- Actual application suction & liquid states -------------------
        T_suction_actual = to_C + superheat_K + 273.15
        rho_actual = PropsSI("D", "T", T_suction_actual, "P", p_evap, self.R)
        h_suction_actual = PropsSI("HMASS", "T", T_suction_actual, "P", p_evap, self.R)
        h_liquid_actual = self._liquid_enthalpy(tc_C, subcooling_K, p_cond)

        density_ratio = rho_actual / rho_rated

        m_dot_actual_kgh = cat["m_dot_kgh"] * density_ratio
        q_ref_actual = h_suction_actual - h_liquid_actual          # J/kg
        Q_actual_w = (m_dot_actual_kgh / 3600.0) * q_ref_actual
        P_actual_w = cat["P_w"] * density_ratio

        return {
            "Q_w": Q_actual_w,
            "P_w": P_actual_w,
            "m_dot_kgh": m_dot_actual_kgh,
            "COP": Q_actual_w / P_actual_w,
            "density_ratio": density_ratio,
            "h_suction": h_suction_actual,
            "catalog": cat,
        }


class CompressorBank:
    """Bank of N identical Compressor units in parallel on a common suction/
    discharge manifold (two units, per the project's dual 75 kW circuit
    design). Both units always share the same to/tc, since they are
    physically on the same lines.

    Provides the load-staging + evaporating-temperature "float" logic needed
    to run the system at an arbitrary required duty and condensing
    temperature (e.g. an IPLV point or an annual-simulation hour): given
    Q_required and tc, solve for the to at which the ACTIVE compressor(s)
    exactly deliver Q_required, staging a second unit on if one unit alone
    cannot reach the load within a sensible evaporating-temperature range.
    """

    def __init__(self, n_units=2, superheat_K=8.0, subcooling_K=0.0,
                 to_bounds_C=(-10.0, 15.0)):
        self.units = [Compressor() for _ in range(n_units)]
        self.n_units = n_units
        self.superheat_K = superheat_K
        self.subcooling_K = subcooling_K
        self.to_bounds_C = to_bounds_C

    def _bank_capacity(self, to_C, tc_C, n_active):
        """Total cooling capacity [W] of n_active identical units at (to, tc)."""
        perf = self.units[0].performance(to_C, tc_C, self.superheat_K, self.subcooling_K)
        return n_active * perf["Q_w"]

    def _solve_to(self, Q_required_w, tc_C, n_active):
        """Root-find the evaporating temperature at which n_active units
        deliver exactly Q_required_w, within to_bounds_C. Returns
        (to_solved, status) where status is 'ok', 'below_floor' (required
        load is below what n_active units can reach even at the coldest
        allowed to -- a genuine turndown-floor limit, not a solver failure),
        or 'above_ceiling' (required load exceeds n_active units' capacity
        even at the warmest allowed to)."""
        lo, hi = self.to_bounds_C

        def f(to_C):
            return self._bank_capacity(to_C, tc_C, n_active) - Q_required_w

        f_lo, f_hi = f(lo), f(hi)
        if f_lo > 0:
            return lo, "below_floor"     # even the coldest to over-delivers
        if f_hi < 0:
            return hi, "above_ceiling"   # even the warmest to under-delivers
        return brentq(f, lo, hi), "ok"

    def solve(self, Q_required_w, tc_C):
        """Stage compressors and solve for the operating evaporating
        temperature at a required total duty and condensing temperature.

        Staging rule: try the smallest n_active (starting at 1) that can
        reach Q_required_w within to_bounds_C ('ok' status); use the
        largest n_active (all units) if none of the smaller stagings can
        reach it. If even all units together cannot reach it, this is a
        genuine turndown-floor (load too low) or capacity-ceiling (load too
        high) limit of this compressor bank -- flagged via `status` in the
        returned dict rather than silently forcing a number.
        """
        status = None
        for n_active in range(1, self.n_units + 1):
            to_solved, status = self._solve_to(Q_required_w, tc_C, n_active)
            if status != "above_ceiling":
                # Either 'ok' (this staging reaches the load), or
                # 'below_floor' (fewer/equal units already over-deliver even
                # at the coldest to -- adding MORE units only raises the
                # floor further, so this is already the closest staging;
                # stop here rather than escalating).
                break
        # If the loop exits on 'above_ceiling' even at n_active == n_units,
        # this is the full bank's capacity ceiling -- the closest it can get.

        perf = self.units[0].performance(to_solved, tc_C, self.superheat_K, self.subcooling_K)
        Q_total = n_active * perf["Q_w"]
        P_total = n_active * perf["P_w"]
        m_dot_total = n_active * perf["m_dot_kgh"]

        return {
            "n_active": n_active,
            "to_C": to_solved,
            "tc_C": tc_C,
            "Q_total_w": Q_total,
            "P_total_w": P_total,
            "m_dot_total_kgh": m_dot_total,
            "COP": Q_total / P_total,
            "status": status,
            "per_unit": perf,
        }

    def solve_with_bypass(self, Q_required_w, tc_C):
        """Same as solve(), but when the required load sits BELOW the
        bank's turndown floor ('below_floor'), models hot-gas bypass:
        the compressor(s) keep running at the floor operating point (they
        cannot unload below it), and a bypass valve recirculates enough
        hot discharge gas back to suction that the NET duty delivered to
        the evaporator load equals Q_required_w exactly.

        Physically, bypassing does not reduce the compressor's own mass
        flow or power draw -- it just throws away part of the capacity as
        recirculated hot gas -- so P stays at the floor's own value while
        Q drops to the target, which is exactly why COP collapses in this
        regime (visible in the returned 'COP' field).

        Returns the same dict shape as solve(), with 'status' relabeled
        'below_floor_bypass' when this correction was applied so callers
        can distinguish it from a normal 'ok' solve.
        """
        result = self.solve(Q_required_w, tc_C)
        if result["status"] != "below_floor":
            return result

        Q_floor = result["Q_total_w"]
        P_floor = result["P_total_w"]
        bypass_fraction = 1.0 - Q_required_w / Q_floor

        result = dict(result)  # shallow copy, don't mutate solve()'s return
        result["Q_total_w"] = Q_required_w          # net duty delivered to the load
        result["P_total_w"] = P_floor                # compressor power unchanged
        result["COP"] = Q_required_w / P_floor
        result["bypass_fraction"] = bypass_fraction
        result["status"] = "below_floor_bypass"
        return result


if __name__ == "__main__":
    # ------------------------------------------------------------------ #
    # DESIGN-POINT CHECK: reproduce the vendor selection file's stated
    # rated performance at to=7 degC, tc=50 degC, then apply the density-
    # ratio correction for the actual 8 K useful superheat (15 degC return
    # gas) vs. the EN12900 rating basis (20 degC return gas).
    # ------------------------------------------------------------------ #
    print("=== SINGLE-UNIT DESIGN POINT: Bitzer 4FEP-35Z, to=7C, tc=50C ===")
    comp = Compressor()
    to_design, tc_design = 7.0, 50.0

    cat = comp.catalog_performance(to_design, tc_design)
    print(f"Catalog (EN12900 rating basis, 20C return gas, 0K subcool):")
    print(f"  Q = {cat['Q_w']/1000:.2f} kW   P = {cat['P_w']/1000:.2f} kW   "
          f"m_dot = {cat['m_dot_kgh']:.1f} kg/h   COP = {cat['Q_w']/cat['P_w']:.2f}")

    perf = comp.performance(to_design, tc_design, superheat_K=8.0, subcooling_K=0.0)
    print(f"\nCorrected (actual 8K useful superheat, 0K subcool):")
    print(f"  Q = {perf['Q_w']/1000:.2f} kW   P = {perf['P_w']/1000:.2f} kW   "
          f"m_dot = {perf['m_dot_kgh']:.1f} kg/h   COP = {perf['COP']:.2f}")
    print(f"  (density ratio applied: {perf['density_ratio']:.4f})")

    print("\n*** RECONCILED: main.py evaluates this same map at the condenser")
    print("*** chapter's established design point (tc=52C, 5K subcooling) instead")
    print("*** of this demo's 50C/0K vendor-selection point -- the map itself is a")
    print("*** continuous polynomial in (to, tc) and is valid at either condition.")

    # ------------------------------------------------------------------ #
    # TWO-COMPRESSOR BANK: stage + solve evaporating temperature across
    # the four IPLV-style load points, at a representative condensing
    # temperature, WITH hot-gas-bypass covering any point below the
    # bank's turndown floor.
    # ------------------------------------------------------------------ #
    print("\n\n=== TWO-COMPRESSOR BANK: staging across load points (tc=50C) ===")
    bank = CompressorBank(n_units=2, superheat_K=8.0, subcooling_K=0.0)
    tc_fixed = 50.0
    for Q_kw in [150.0, 112.5, 75.0, 37.5, 30.0]:
        result = bank.solve_with_bypass(Q_kw * 1000.0, tc_fixed)
        flag = "" if result["status"] == "ok" else f"  <-- {result['status'].upper()}"
        print(f"Q_required = {Q_kw:6.1f} kW  ->  n_active = {result['n_active']}  "
              f"to = {result['to_C']:6.2f} C  Q_achieved = {result['Q_total_w']/1000:6.2f} kW  "
              f"P_total = {result['P_total_w']/1000:6.2f} kW  COP = {result['COP']:.2f}{flag}")

    print("\n*** FINDING: this compressor bank has a genuine LOW-LOAD TURNDOWN FLOOR.")
    print("*** Even a single unit at the coldest allowed evaporating temperature")
    print("*** (to=-10C) still delivers ~42 kW -- above the 30-37.5 kW minimum IT")
    print("*** load. Below that floor, solve_with_bypass() models hot-gas bypass:")
    print("*** the compressor keeps running at the floor point and a bypass valve")
    print("*** throttles the NET delivered duty down to the target, at the cost of")
    print("*** COP (power stays at the floor's value while useful Q drops) -- a real")
    print("*** but inefficient way to cover this range pending a better strategy")
    print("*** (cylinder unloading, a VFD, a smaller trim compressor, or handing")
    print("*** low-load hours to the economizer instead).")

    # ------------------------------------------------------------------ #
    # IPLV.IP -- AHRI 550/590-2023 S5.2: weighted COP at 100/75/50/25%
    # load, each evaluated at ITS OWN reduced condenser-entering-air
    # temperature (35/26.7/18.3/12.8 C). This plant condenses against an
    # intermediate glycol loop (not air directly), so each AHRI air temp
    # is converted to a refrigerant condensing temperature via the SAME
    # fixed approach established at the 100% design point (tc=52C at
    # T_air=35C design -> 17 K approach), held constant across all four
    # points. SIMPLIFICATION: a rigorous treatment would re-solve the dry
    # cooler's eps-NTU model (dry_cooler.py) at each point's reduced air
    # flow/temperature instead of assuming a fixed approach -- flagged as
    # an open item, not attempted here.
    # ------------------------------------------------------------------ #
    print("\n\n=== IPLV.IP (AHRI 550/590-2023 S5.2) ===")
    T_AIR_DESIGN_C = 35.0
    TC_DESIGN_C = 52.0
    APPROACH_K = TC_DESIGN_C - T_AIR_DESIGN_C  # = 17 K, held fixed across all 4 points

    iplv_points = [
        ("A", 1.00, 0.01, 35.0),
        ("B", 0.75, 0.42, 26.7),
        ("C", 0.50, 0.45, 18.3),
        ("D", 0.25, 0.12, 12.8),
    ]

    iplv = 0.0
    cops = {}
    for label, load_frac, weight, t_air_c in iplv_points:
        tc_point = t_air_c + APPROACH_K
        Q_point = load_frac * 150e3
        r = bank.solve_with_bypass(Q_point, tc_point)
        cops[label] = r["COP"]
        iplv += weight * r["COP"]
        flag = "" if r["status"] == "ok" else f"  <-- {r['status'].upper()}"
        print(f"  {label} ({load_frac*100:5.1f}% load, T_air={t_air_c:5.1f} C -> "
              f"tc={tc_point:5.1f} C): COP = {r['COP']:.2f}  (weight {weight}){flag}")

    print(f"\nIPLV.IP = 0.01*A + 0.42*B + 0.45*C + 0.12*D = {iplv:.2f}")
    print("(Project target: IPLV.IP >= 5.0)")
