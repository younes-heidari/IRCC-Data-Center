import numpy as np
from CoolProp.CoolProp import PropsSI


class EEV:
    """Models the electronic expansion valve (EEV): throttles subcooled liquid
    from the condenser outlet down to the evaporating pressure. Throttling is
    adiabatic with no work, so it is isenthalpic (h_out = h_in).

    INPUTS:
      m_dot_refrigerant : refrigerant mass flow [kg/s]
      p_in              : inlet pressure [Pa]
      h_in              : inlet enthalpy [J/kg]
      p_out             : evaporating pressure [Pa]
      refrigerant       : fluid name [str]

    OUTPUTS:
      h_out, T_out, x_out, delta_p, Kv_required
    """

    def __init__(self, m_dot_refrigerant, p_in, h_in, p_out, refrigerant, Cd=0.65):
        # General properties
        self.R = refrigerant
        self.m_dot_refrigerant = m_dot_refrigerant
        self.Cd = Cd

        # CAREL E3V65 rated flow coefficient
        self.Kv_rated = 1.087  # m^3/h

        # Inlet (subcooled liquid from the condenser)
        self.p_in = p_in
        self.h_in = h_in

        # Outlet (two-phase)
        self.p_out = p_out
        self.h_out = h_in
        self.delta_p = p_in - p_out

        # Sizing outputs
        self.Kv_required = None
        self.Cv_required = None

        self.calc_outlet_state()
        self.calc_sizing()

    # Isenthalpic flash to the evaporating pressure
    def calc_outlet_state(self):
        self.T_out = PropsSI('T', 'P', self.p_out, 'Q', 0, self.R)
        h_f = PropsSI('HMASS', 'P', self.p_out, 'Q', 0, self.R)
        h_g = PropsSI('HMASS', 'P', self.p_out, 'Q', 1, self.R)
        self.x_out = (self.h_out - h_f) / (h_g - h_f)

    # Calculate commercial Flow Coefficient (Kv) for valve selection
    def calc_sizing(self):
        # 1. Get liquid density at inlet
        rho_in = PropsSI('D', 'HMASS', self.h_in, 'P', self.p_in, self.R)

        # 2. Commercial Flow Coefficient (Kv) Calculation
        # Kv is the volumetric flow of water [m^3/h] at a 1 bar pressure drop.
        v_dot_m3h = (self.m_dot_refrigerant / rho_in) * 3600
        delta_p_bar = self.delta_p / 100000.0
        sg = rho_in / 1000.0  # Specific gravity relative to water

        # 3. Calculate required Kv and Cv
        self.Kv_required = v_dot_m3h * np.sqrt(sg / delta_p_bar)
        self.Cv_required = self.Kv_required / 0.865  # Conversion from Kv to Cv


if __name__ == "__main__":
    # Test conditions matching design point
    R = "R290"
    p_cond = PropsSI("P", "T", 45 + 273.15, "Q", 1, R)
    p_evap = PropsSI("P", "T", 5 + 273.15, "Q", 1, R)
    T_in = 45 - 5 + 273.15  # 5 K subcooling
    h_in = PropsSI('HMASS', 'T', T_in, 'P', p_cond, R)
    m_dot = 0.50  # Target mass flow

    valve = EEV(m_dot_refrigerant=m_dot, p_in=p_cond, h_in=h_in, p_out=p_evap, refrigerant=R)

    print(f"=== CAREL E3V65 Validation ===")
    print(f"Required Kv    : {valve.Kv_required:.3f} m^3/h")
    print(f"Valve Rated Kv : {valve.Kv_rated:.3f} m^3/h")

    # Calculate utilization percentage
    utilization = (valve.Kv_required / valve.Kv_rated) * 100
    print(f"Valve Loading  : {utilization:.1f}%")

    if utilization > 90:
        print("WARNING: Valve is nearing maximum capacity.")
    elif utilization < 50:
        print("WARNING: Valve may be oversized for this design point.")
    else:
        print("SUCCESS: Valve selection is optimal.")
