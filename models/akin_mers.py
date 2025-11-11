# models/akin_mers.py

import math
from .common_params import EngineParams, PropellantType
from typing import Literal, Dict

def estimate_engine_mass(thrust_N: float, expansion_ratio: float) -> Dict:
    """
    Estimates liquid pump-fed rocket engine mass based on thrust and expansion ratio.
    Returns a dictionary with total mass and component (term) breakdown.

    Reference:
    Mass Estimating Relations (Akin, ENAE 791), Page 25.
    Formula: M_Rocket_Engine(kg) = 7.81e-4*T(N) + 3.37e-5*T(N)*sqrt(Ae/At) + 59

    Args:
        thrust_N (float): Engine thrust in Newtons.
        expansion_ratio (float): Nozzle expansion ratio (Ae/At).

    Returns:
        Dict: Dictionary with 'total_mass_kg' and 'components_kg'.
    """
    term1 = 7.81e-4 * thrust_N
    term2 = 3.37e-5 * thrust_N * (expansion_ratio ** 0.5)
    term3 = 59.0

    total_mass = term1 + term2 + term3

    components = {
        "thrust_term (7.81e-4*T)": term1,
        "expansion_term (3.37e-5*T*sqrt(Ae/At))": term2,
        "constant_term (59)": term3
    }

    return {
        "total_mass_kg": total_mass,
        "components_kg": components
    }


def estimate_thrust_structure_mass(total_thrust_N: float) -> float:
    """
    Estimates thrust structure mass based on total vehicle thrust.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 25[cite: 3170].
    Formula: M_Thrust_Structure(kg) = 2.55e-4*T(N)

    Args:
        total_thrust_N (float): Total thrust of all engines supported by the structure, in Newtons.

    Returns:
        float: Estimated thrust structure mass in kg.
    """
    return 2.55e-4 * total_thrust_N


def estimate_propellant_tank_mass(volume_m3: float, propellant: Literal["LH2", "LOX", "RP1"]) -> float:
    """
    Estimates propellant tank mass based on propellant volume.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 6 [cite: 2885-2888].
    Formulas:
        M_LH2_Tank(kg) = 9.09 * V_LH2(m^3)
        M_Other_Tank(kg) = 12.16 * V_prop(m^3)

    [cite_start]Note: Page 7 [cite: 2894-2899] provides mass-based MERs, but volume-based is more general.

    Args:
        volume_m3 (float): Volume of the propellant in cubic meters.
        propellant (Literal["LH2", "LOX", "RP1"]): Type of propellant.

    Returns:
        float: Estimated tank mass in kg.
    """
    if propellant == "LH2":
        # [cite_start]M_LH2_Tank(kg) = 9.09 * V_LH2(m^3) [cite: 2886]
        return 9.09 * volume_m3
    else:
        # [cite_start]M_Tank(kg) = 12.16 * V_prop(m^3) [cite: 2888]
        # [cite_start]This is used for LOX and RP-1 per the regression plot[cite: 2858].
        return 12.16 * volume_m3


def tank_surface_area_m2(propellant_mass: float, propellant_density: float):
    """
    Calculate the estimated area of the tank surface based on propellant mass and density.
    """
    # Need area to find LOX tank insulation mass - assume a sphere
    volume = propellant_mass / propellant_density
    sphere_radius = (volume / (4 * math.pi / 3)) ** (1/3)
    return 4 * math.pi * (sphere_radius ** 2)


def estimate_cryo_insulation_mass(tank_surface_area_m2: float, propellant: Literal["LH2", "LOX"]) -> float:
    """
    Estimates cryogenic insulation mass based on tank surface area.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 8 [cite: 2904-2907].
    Formulas:
        M_LH2_Insulation(kg) = 2.88 * A_tank(m^2)
        M_LOX_Insulation(kg) = 1.123 * A_tank(m^2)

    Args:
        tank_surface_area_m2 (float): Surface area of the tank in square meters.
        propellant (Literal["LH2", "LOX"]): Type of cryogenic propellant.

    Returns:
        float: Estimated insulation mass in kg.
    """
    if propellant == "LH2":
        # [cite_start]M_LH2_Insulation <kg> = 2.88 * A_tank <kg/m^2> [cite: 2904-2905]
        return 2.88 * tank_surface_area_m2
    elif propellant == "LOX":
        # [cite_start]M_LOX_Insulation <kg> = 1.123 * A_tank <kg/m^2> [cite: 2906-2907]
        return 1.123 * tank_surface_area_m2
    elif propellant == "Methane":
        # estimated, needs double-check
        return 1.0 * tank_surface_area_m2
    else:
        return 0.0  # No insulation needed for non-cryo


def estimate_fairing_mass(fairing_surface_area_m2: float) -> float:
    """
    Estimates fairing/shroud mass based on surface area.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 20[cite: 3087].
    Formula: M_fairing(kg) = 4.95 * (A_fairing(m^2))^1.15

    Args:
        fairing_surface_area_m2 (float): Surface area of the fairing in square meters.

    Returns:
        float: Estimated fairing mass in kg.
    """
    return 4.95 * (fairing_surface_area_m2 ** 1.15)


def estimate_avionics_mass(vehicle_gross_mass_kg: float) -> float:
    """
    Estimates avionics mass based on vehicle gross mass.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 20[cite: 3091].
    Formula: M_avionics(kg) = 10 * (M_o(kg))^0.361

    Args:
        vehicle_gross_mass_kg (float): Vehicle gross mass (M_o) in kg.

    Returns:
        float: Estimated avionics mass in kg.
    """
    return 10.0 * (vehicle_gross_mass_kg ** 0.361)


def estimate_wiring_mass(vehicle_gross_mass_kg: float, vehicle_length_m: float) -> float:
    """
    Estimates wiring mass based on vehicle gross mass and length.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 20[cite: 3092].
    Formula: M_wiring(kg) = 1.058 * sqrt(M_o(kg)) * l^0.25

    Args:
        vehicle_gross_mass_kg (float): Vehicle gross mass (M_o) in kg.
        vehicle_length_m (float): Vehicle total length (l) in meters.

    Returns:
        float: Estimated wiring mass in kg.
    """
    return 1.058 * (vehicle_gross_mass_kg ** 0.5) * (vehicle_length_m ** 0.25)


def estimate_gimbal_mass(engine_thrust_N: float, chamber_pressure_Pa: float) -> float:
    """
    Estimates gimbal mass for a single engine.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 26[cite: 3176].
    Formula: M_Gimbals(kg) = 237.8 * [T(N) / P_c(Pa)]^0.9375

    [cite_start]Note: The reference slide [cite: 3176] uses P_0(Pa), but the example calculation
    [cite_start]on page 28 [cite: 3198-3199] [cite_start]and mass summary on page 30 [cite: 3213] confirm
    this refers to chamber pressure (P_c) and that the mass is per engine.
    (6 engines * 13.5 kg/eng ~= 81 kg total).

    Args:
        engine_thrust_N (float): Thrust of a single engine in Newtons.
        chamber_pressure_Pa (float): Engine chamber pressure in Pascals.

    Returns:
        float: Estimated gimbal mass in kg.
    """
    if chamber_pressure_Pa <= 0:
        return 0.0

    ratio = engine_thrust_N / chamber_pressure_Pa
    return 237.8 * (ratio ** 0.9375)