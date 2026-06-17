# IRCC-Data-Center

Vapor-compression chiller design for a ~150 kW edge-compute data center
(IRCC 2026, Project #3). Component models written in Python using
[CoolProp](http://www.coolprop.org/) for fluid properties.

## Components

| file | component | what it does |
|---|---|---|
| `compressor.py` | reciprocating compressor | displacement → mass flow, work, discharge state |
| `evaporator`    | chilled-water brazed-plate evaporator | heat duty, pressure drop, UA |

## Usage

Each component is an importable class — instantiate it with the inputs below
and read the results off the returned object:

```python
from evaporator import Evaporator
evap = Evaporator(m_dot_refrigerant=..., p_in=..., Q=150e3, refrigerant="R290",
                  T_water_in=..., T_water_out=..., h_in=...,
                  D_h=..., A_flow=..., L=..., beta=..., Lambda=..., N_cp=...)
print(evap.UA, evap.delta_p, evap.T_out)
```

## Evaporator interface

**Cycle-side inputs** (provided by teammates / operating conditions):

| param | meaning | source |
|---|---|---|
| `m_dot_refrigerant` | refrigerant mass flow [kg/s] | compressor |
| `p_in` | evaporating pressure [Pa] | operating condition |
| `Q` | cooling load [W] | design spec (150 kW) |
| `refrigerant` | fluid name | working-fluid choice |
| `T_water_in` / `T_water_out` | chilled water 21 / 15 °C [K] | spec |
| `h_in` | refrigerant inlet enthalpy [J/kg] | EEV / condenser outlet |

**Geometry inputs** (plate datasheet — owner's responsibility):
`D_h`, `A_flow`, `L`, `beta`, `Lambda`, `N_cp`

**Outputs** (for other component models):

| output | use |
|---|---|
| `h_out`, `p_out`, `T_out`, `superheat` | compressor suction state |
| `delta_p` | evaporator refrigerant-side pressure loss |
| `m_dot_water`, `cp_water` | chilled-water pump sizing |
| `UA`, `LMTD` | plate-HX selection / sizing |

## Pressure-drop model

Two-zone refrigerant-side Δp:

- **Boiling zone** — Han, Lee & Kim (2003) two-phase friction correlation
  (Eqs. 11, 12, 17–20), integrated over vapor quality, with the `N_cp`
  channel factor.
- **Superheat zone** — single-phase Blasius friction (short zone, small
  contribution).

## Open items / caveats

- **Han validity** — calibrated on R410A/R22, β 20–45°, G 13–34 kg/m²·s.
  We use R290, so confirm the mass flux lands in range and justify the
  extrapolation. (The run prints `G` for this check.)
- **`N_cp` factor** — Han Eq. 12 is strictly self-consistent only at the rig
  value `N_cp = 2`; sanity-check Δp lands ~10–30 kPa with real geometry.
- **Compressor sizing** — its displacement must be sized so the geometry-based
  mass flow matches the load-based design flow (the run prints both).
- **Superheat / subcooling** — design choices (currently 8 K / 2 K) to justify.
