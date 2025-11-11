# vehicle_definitions.py

from models.common_params import EngineParams

# --- Propellant Density Constants (kg/m^3) ---
# Using values from Akin (ENAE 791), Page 4 [cite: 1233, 1234]
DENSITY_LOX = 1140.0
DENSITY_LH2 = 71.0
DENSITY_RP1 = 820.0  # From Akin (ENAE 791), Page 7 [cite: 1280]
# Standard gravity for Isp <-> Ve conversion
G0 = 9.80665

def get_le5_engine() -> EngineParams:
    """
    Defines a test case for the LE-5 engine.

    This is a LOX/LH2 Gas Generator (GG) cycle engine.

    Reference:
    Data primarily from Tizón & Román, 2017, Table 7.
    """
    return EngineParams(
        thrust_vac_N=103_000.0,
        isp_vac_s=450.0,  # Note: Isp not in Table 7, using typical value from Akin [cite: 1214]
        chamber_pressure_Pa=3.65e6,
        propellant_type="LOX/LH2",
        cycle_type="GG",
        mixture_ratio=5.5,
        expansion_ratio=140.0,
        fuel_density=DENSITY_LH2,
        oxidizer_density=DENSITY_LOX,
        num_chambers=1
    )


def get_ssme_engine() -> EngineParams:
    """
    Defines a test case for the SSME (RS-25) engine.

    This is a LOX/LH2 Staged Combustion (SC) cycle engine.

    Reference:
    Data primarily from Tizón & Román, 2017, Table 7.
    """
    return EngineParams(
        thrust_vac_N=2_280_000.0,
        isp_vac_s=452.0,  # Known public value
        chamber_pressure_Pa=2.04e7,
        propellant_type="LOX/LH2",
        cycle_type="SC",
        mixture_ratio=6.0,
        expansion_ratio=77.5,
        fuel_density=DENSITY_LH2,
        oxidizer_density=DENSITY_LOX,
        num_chambers=1
    )


def get_rl10a_engine() -> EngineParams:
    """
    Defines a test case for the RL10A-3-A engine.

    This is a LOX/LH2 Expander (EX) cycle engine.

    Reference:
    Data primarily from Tizón & Román, 2017, Table 7.
    """
    return EngineParams(
        thrust_vac_N=73_400.0,
        isp_vac_s=444.0,  # Known public value
        chamber_pressure_Pa=3.20e6,
        propellant_type="LOX/LH2",
        cycle_type="EX",
        mixture_ratio=5.0,
        expansion_ratio=61.1,
        fuel_density=DENSITY_LH2,
        oxidizer_density=DENSITY_LOX,
        num_chambers=1
    )


def get_rd120_engine() -> EngineParams:
    """
    Defines a test case for the RD-120 engine.

    This is a LOX/RP-1 Staged Combustion (SC) cycle engine.

    Reference:
    Data from Zandbergen, 2015, Table 3 [cite: 96]
    """
    return EngineParams(
        thrust_vac_N=833_565.0,
        isp_vac_s=350.0,  # Typical for LOX/RP-1 SC
        chamber_pressure_Pa=16.0e6,  # Known public value
        propellant_type="LOX/RP1",
        cycle_type="SC",
        mixture_ratio=2.6,  # Known public value
        expansion_ratio=106.0,  # Known public value
        fuel_density=DENSITY_RP1,
        oxidizer_density=DENSITY_LOX,
        num_chambers=1
    )