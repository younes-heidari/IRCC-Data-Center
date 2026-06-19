import numpy as np
from CoolProp.CoolProp import PropsSI

class Compressor:

    def __init__(self, temp_in, n, pressure_in, pressure_out, refrigerant,
                 displacement=None, eta_vol=0.9, zeta_v=0, eta_is=0.7, power_input=None, eta_m=0.92, mass_flow=None):
        """Calculate the thermodynamic properties of a chosen compressor.
        Enter either volumetric efficiency (eta_vol) or mass flow. Otherwise, eta_vol is set to standard value.
        Enter either isentropic efficiency (eta_is) or electrical power input (power_input)

        :param temp_in: inlet temperature in Kelvin.
        :param n: Rotational speed in rpm.
        :param pressure_in: Inlet pressure in Pa.
        :param pressure_out: Outlet pressure in Pa.
        :param refrigerant: Refrigerant string for coolprop ("RXXX").
        :param displacement: Displacement in m^3/h.
        :param eta_vol: Volumetric efficiency. Default value: 0.9
        :param zeta_v:  Heat transfer ratio (q_v/w_v). q_v > 0 means heat loss.
                        Default value: 0 (adiabatic). Typical value from literature: -0.07
        :param eta_is: Isentropic efficiency. Default value: 0.7. Mutually exclusive with power_input.
        :param power_input: Electrical power input in W. Mutually exclusive with eta_is.
        :param eta_m: Motor efficiency, only used with power_input. Default value: 0.92.
        :param mass_flow: Mass flow in kg/h. Mutually exclusive with eta_vol.

        Notes:
            - Enter either volumetric efficiency (eta_vol) or mass flow. Otherwise, eta_vol is set to default value.
            - Enter either isentropic efficiency (eta_is) or electrical power input (power_input)
        """

        # Validation: enter either eta_vol or mass flow to calculate the other
        if eta_is !=0.7 and power_input is not None:
            raise ValueError("Either 'eta_is' or 'power_input' must be specified, not both.")
        if eta_vol != 0.9 and mass_flow is not None:
            raise ValueError("Either 'eta_vol' or 'mass_flow' must be specified, not both.")
        if not displacement:
            raise ValueError("'displacement' must be specified. Please use a subclass for geometric properties")

        # General properties
        self.R = refrigerant
        self.n = n
        self.zeta_v = zeta_v
        self.eta_m = eta_m
        self.V_dot_disp = displacement / 3600 # convert from unit m^3/h to m^3/s

        # Efficiency / Power input
        self.eta_is = eta_is if power_input is None else None
        self.P_el = power_input
        self.eta_vol = eta_vol if mass_flow is None else None
        self.m_dot = mass_flow / 3600  if mass_flow else None # convert from unit kg/h to kg/s

        # Derived quantities (calculated below)
        self.w_v = None
        self.W_v = None

        # Inlet State (1)
        self.T_1 = temp_in
        self.p_1 = pressure_in
        self.h_1 = None
        self.s_1 = None
        self.rho_1 = None
        self.V_dot_1 = None

        # Outlet State (2)
        self.T_2 = None
        self.p_2 = pressure_out
        self.h_2is = None
        self.h_2 = None

        self.pressure_ratio = self.p_2 / self.p_1

        # Calculation methods for completing properties and EOS
        self.calc_inlet_state()
        self.calc_outlet_state()

    def calc_inlet_state(self):
        """Calculate thermodynamic properties of the inlet state based on the set temperature and pressure of the inlet state. Mass flow is calculated using a set volumetric efficiency."""
        self.rho_1 = PropsSI('D', 'T', self.T_1, 'P', self.p_1, self.R)
        self.h_1 = PropsSI('HMASS', 'T', self.T_1, 'P', self.p_1, self.R)
        self.s_1 = PropsSI('SMASS', 'T', self.T_1, 'P', self.p_1, self.R)

        if self.m_dot is not None:
            # mass flow given - calculate eta_vol
            self.V_dot_1 = self.m_dot / self.rho_1
            self.eta_vol = self.V_dot_1 / self.V_dot_disp
        else:
            # eta_vol given - calculate mass flow
            self.V_dot_1 = self.V_dot_disp * self.eta_vol
            self.m_dot = self.V_dot_1 * self.rho_1


    def calc_outlet_state(self):
        """Calculate thermodynamic properties for outlet of the compressor.

        Uses isentropic efficiency (eta_is) and heat loss factor (zeta_v = q_v / w_v)
        to determine the actual outlet state. A positive zeta_v means heat is added
        from the environment. A negative heat loss means heat transfer to the environment,
        therefore lowering the discharge temperature compared to adiabatic.
        """

        self.h_2is = PropsSI('HMASS', 'SMASS', self.s_1, 'P', self.p_2, self.R)

        if self.eta_is is not None:
            # Calculate specific compressor work, accounting for heat loss with zeta_v.
            self.w_v = (self.h_2is - self.h_1) / (self.eta_is * (1 + self.zeta_v))
            self.W_v = self.w_v * self.m_dot
            self.P_el = self.W_v / self.eta_m
        else:
            # Calculate compressor work of power input and motor efficiency.
            self.W_v = self.P_el * self.eta_m
            self.w_v = self.W_v / self.m_dot
            # Calculate isentropic efficiency of isentropic compression compared to actual power draw
            self.eta_is = (self.h_2is - self.h_1) * self.eta_m / self.w_v

        # Actual outlet enthalpy: h_2 = h_1 + w_v * (1 + zeta_v)
        self.h_2 = self.h_1 + self.w_v * (1 + self.zeta_v)
        self.T_2 = PropsSI('T', 'HMASS', self.h_2, 'P', self.p_2, self.R)


class RecipCompressor(Compressor):

    def __init__(self, temp_in, n, pressure_in, pressure_out, refrigerant,
                 d=None, z=None, s=None, displacement=None,
                 eta_vol=0.9, zeta_v=0.0, eta_is=0.7, power_input=None, eta_m=0.92, mass_flow=None):
        """
        Calculate the thermodynamic properties of a chosen reciprocating compressor.

        :param d: Diameter of the cylinder in m.
        :param z: Amount of cylinders
        :param s: Stroke of the cylinders in m.

        All other parameters are passed to class Compressor. See Compressor documentation.

        Notes:
            - Enter either volumetric efficiency (eta_vol) or mass flow.
            - Enter either displacement or d, n and s.
        """

        geo_given = all(v is not None for v in [d, z, s])

        if not displacement and not geo_given:
            raise ValueError("Either displacement or d, z and s must be specified.")
        elif displacement and all([d, z, s]):
            raise ValueError("Either displacement or d, z and s must be specified, not both.")

        if geo_given:
            self.d = d
            self.z = z
            self.stroke = s

            # Calculate displacement, pass it to the class inheritance
            volume_cyl = np.pi / 4 * d ** 2
            displacement = volume_cyl * s * z * n * 60 # convert from rpm to rph - m^3/min to m^3/h
        else:
            self.d = None
            self.z = None
            self.stroke = None

        super().__init__(temp_in=temp_in, n=n, pressure_in=pressure_in, pressure_out=pressure_out,
                         refrigerant=refrigerant, displacement=displacement,
                         eta_vol=eta_vol, zeta_v=zeta_v, eta_is=eta_is,
                         power_input=power_input, eta_m=eta_m,mass_flow=mass_flow)


### zeta_v is a guess value by the script of Prof. Bradshaw (we'll have to set a value based on our compressor choice