import numpy as np
from CoolProp.CoolProp import PropsSI


class EEV:
    """Models the electronic expansion valve (EEV): throttles subcooled liquid
    from the condenser outlet down to the evaporating pressure. Throttling is
    adiabatic with no work, so it is isenthalpic (h_out = h_in); the refrigerant
    flashes to a low-quality two-phase state at the evaporator inlet. This closes
    the cycle: condenser -> EEV -> evaporator.

    INPUTS (from the cycle side -- provided by teammates/operating conditions):
      m_dot_refrigerant : refrigerant mass flow [kg/s]      (from compressor)
      p_in              : inlet pressure [Pa]               (= condenser outlet p_out)
      h_in              : inlet enthalpy [J/kg]             (= condenser outlet h_out)
      p_out             : outlet / evaporating pressure [Pa] (= evaporator p_in)
      refrigerant       : fluid name [str]                  (e.g. "R290")
      Cd                : valve discharge coefficient [-]   (datasheet; ~0.6-0.7)

    OUTPUTS (for teammates' component models):
      h_out                -> evaporator inlet enthalpy (= h_in, isenthalpic)
      T_out                -> evaporating temperature (saturation at p_out)
      x_out                -> evaporator inlet quality (flash-gas fraction)
      delta_p              -> pressure drop taken across the valve
      A_orifice, d_orifice -> EEV orifice sizing estimate (for selection)

    NOTE: in operation the EEV modulates this orifice to hold a target evaporator
    superheat (the compressor/EEV control interaction). This is the steady design
    point; the control behaviour is a separate discussion item.
    """

    def __init__(self, m_dot_refrigerant, p_in, h_in, p_out, refrigerant, Cd=0.65):
        # General properties
        self.R = refrigerant
        self.m_dot_refrigerant = m_dot_refrigerant
        self.Cd = Cd

        # Inlet (subcooled liquid from the condenser)
        self.p_in = p_in
        self.h_in = h_in
        self.T_in = None
        self.subcooling_in = None

        # Outlet (two-phase, to the evaporator) -- isenthalpic throttling
        self.p_out = p_out
        self.h_out = h_in
        self.T_out = None
        self.x_out = None

        self.delta_p = p_in - p_out

        # Orifice sizing
        self.A_orifice = None
        self.d_orifice = None

        self.calc_inlet_state()
        self.calc_outlet_state()
        self.calc_orifice()

    # Inlet temperature and how subcooled the liquid is at the valve inlet
    def calc_inlet_state(self):
        self.T_in = PropsSI('T', 'HMASS', self.h_in, 'P', self.p_in, self.R)
        T_bubble = PropsSI('T', 'P', self.p_in, 'Q', 0, self.R)
        self.subcooling_in = T_bubble - self.T_in

    # Isenthalpic flash to the evaporating pressure -> two-phase outlet state
    def calc_outlet_state(self):
        self.T_out = PropsSI('T', 'P', self.p_out, 'Q', 0, self.R)   # saturation temp
        h_f = PropsSI('HMASS', 'P', self.p_out, 'Q', 0, self.R)
        h_g = PropsSI('HMASS', 'P', self.p_out, 'Q', 1, self.R)
        self.x_out = (self.h_out - h_f) / (h_g - h_f)

    # First-cut orifice sizing from an incompressible-liquid orifice balance:
    #   m_dot = Cd * A * sqrt(2 * rho * delta_p)   (rho at the liquid inlet state).
    # Real EEVs flash across the port (two-phase), so manufacturers rate by
    # Kv / capacity tables; this gives a flow area / port diameter for selection.
    def calc_orifice(self):
        rho_in = PropsSI('D', 'HMASS', self.h_in, 'P', self.p_in, self.R)
        self.A_orifice = self.m_dot_refrigerant / (self.Cd * np.sqrt(2 * rho_in * self.delta_p))
        self.d_orifice = np.sqrt(4 * self.A_orifice / np.pi)