# models/zandbergen_engine.py

from .common_params import EngineParams
from typing import Dict


def estimate_engine_mass(params: EngineParams) -> Dict:
    """
    Estimates engine mass based on thrust and propellant type.
    Returns a dictionary with total mass.

    Implements formulas M-C1 (Hydro-lox) and M-S1 (Storable/Kero-lox).
    Thrust (THRUST) must be in Newtons.

    Reference:
    Simple mass and size estimation relationships... (Zandbergen, 2015), Table 4.


    Formulas:
    - M-C1 (Hydro-lox): M(kg) = 0.00514 * THRUST(N)^0.92068
    - M-S1 (Storable/Kero-lox): M(kg) = 1.104E-3 * THRUST(N) + 27.702

    Args:
        params (EngineParams): The engine's parameter object.

    Returns:
        Dict: Dictionary with 'total_mass_kg' and 'components_kg'.
    """
    thrust_N = params.thrust_vac_N

    if params.propellant_type == "LOX/LH2":
        # Formula M-C1: M = 0.00514 * THRUST^0.92068
        mass = 0.00514 * (thrust_N ** 0.92068)

    elif params.propellant_type in ["LOX/RP1", "LOX/LCH4", "Storable"]:
        # Formula M-S1: M = 1.104E-3 * THRUST + 27.702
        # We assume LCH4 (methane) follows the kero-lox/storable trend.
        mass = (1.104e-3 * thrust_N) + 27.702

    else:
        raise ValueError(f"Propellant type {params.propellant_type} not supported by Zandbergen M-C1/M-S1 model")

    return {
        "total_mass_kg": mass,
        "components_kg": {
            "N/A (Monolithic Model)": mass
        }
    }


def estimate_engine_length(params: EngineParams) -> float:
    """
    Estimates overall engine length based on thrust and propellant type.

    Implements formulas L-C1 (Hydro-lox) and L-S1 (Storable/Kero-lox).
    Thrust (THRUST) must be in Newtons.

    Reference:
    Simple mass and size estimation relationships... (Zandbergen, 2015), Table 5.
    [cite_start][cite: 887]

    Formulas:
    - L-C1 (Hydro-lox): L(m) = 0.1667 * THRUST(N)^0.2238
    - L-S1 (Storable/Kero-lox): L(m) = 0.1362 * THRUST(N)^0.2279

    Args:
        params (EngineParams): The engine's parameter object.

    Returns:
        float: Estimated engine length in meters.
    """
    thrust_N = params.thrust_vac_N

    if params.propellant_type == "LOX/LH2":
        # [cite_start]Formula L-C1: L = 0.1667 * THRUST^0.2238 [cite: 887]
        length = 0.1667 * (thrust_N ** 0.2238)

    elif params.propellant_type in ["LOX/RP1", "LOX/LCH4", "Storable"]:
        # [cite_start]Formula L-S1: L = 0.1362 * THRUST^0.2279 [cite: 887]
        length = 0.1362 * (thrust_N ** 0.2279)

    else:
        raise ValueError(f"Propellant type {params.propellant_type} not supported by Zandbergen L-C1/L-S1 model")

    return length


def estimate_engine_diameter(params: EngineParams) -> float:
    """
    Estimates overall engine envelope diameter based on thrust and propellant type.

    Implements formulas D-C1 (Hydro-lox) and D-S1 (Storable/Kero-lox).
    Thrust (THRUST) must be in Newtons.

    Reference:
    Simple mass and size estimation relationships... (Zandbergen, 2015), Table 5.
    [cite_start][cite: 887]

    Formulas:
    - D-C1 (Hydro-lox): D(m) = 0.1503 * THRUST(N)^0.19
    - D-S1 (Storable/Kero-lox): D(m) = 0.0455 * THRUST(N)^0.2745

    Args:
        params (EngineParams): The engine's parameter object.

    Returns:
        float: Estimated engine envelope diameter in meters.
    """
    thrust_N = params.thrust_vac_N

    if params.propellant_type == "LOX/LH2":
        # [cite_start]Formula D-C1: D = 0.1503 * THRUST^0.19 [cite: 887]
        diameter = 0.1503 * (thrust_N ** 0.19)

    elif params.propellant_type in ["LOX/RP1", "LOX/LCH4", "Storable"]:
        # [cite_start]Formula D-S1: D = 0.0455 * THRUST^0.2745 [cite: 887]
        diameter = 0.0455 * (thrust_N ** 0.2745)

    else:
        raise ValueError(f"Propellant type {params.propellant_type} not supported by Zandbergen D-C1/D-S1 model")

    return diameter