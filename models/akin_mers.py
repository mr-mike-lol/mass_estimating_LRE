# models/akin_mers.py

import math
from .common_params import EngineParams, PropellantType
from typing import Literal, Dict


def estimate_engine_mass(thrust_N: float, expansion_ratio: float) -> Dict:
    """
    Estimates liquid pump-fed rocket engine mass based on thrust and expansion ratio.
    Returns a dictionary with total mass and component (term) breakdown.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 25 [cite: 357-360].
    [cite_start]Formula: M_Rocket_Engine(kg) = 7.81e-4*T(N) + 3.37e-5*T(N)*sqrt(Ae/At) + 59 [cite: 358, 360]

    Args:
        thrust_N (float): Engine thrust in Newtons.
        expansion_ratio (float): Nozzle expansion ratio (Ae/At).

    Returns:
        Dict: Dictionary with 'total_mass_kg' and 'components_kg'.
    """
    # [cite_start]Term 1 from formula [cite: 358]
    term1 = 7.81e-4 * thrust_N
    # [cite_start]Term 2 from formula [cite: 360]
    term2 = 3.37e-5 * thrust_N * (expansion_ratio ** 0.5)
    # [cite_start]Term 3 from formula [cite: 360]
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
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 25 [cite: 362-363].
    [cite_start]Formula: M_Thrust_Structure(kg) = 2.55e-4*T(N) [cite: 363]

    Args:
        total_thrust_N (float): Total thrust of all engines supported by the structure, in Newtons.

    Returns:
        float: Estimated thrust structure mass in kg.
    """
    # [cite_start]Formula from source [cite: 363]
    return 2.55e-4 * total_thrust_N


def estimate_propellant_tank_mass(volume_m3: float, propellant: Literal["LH2", "LOX", "RP1"]) -> float:
    """
    Estimates propellant tank mass based on propellant volume.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 6 [cite: 77-81].
    [cite_start]Formulas derived from regression plot on Page 5 [cite: 51-73].
    [cite_start]M_LH2_Tank(kg) = 9.09 * V_LH2(m^3) [cite: 79]
    [cite_start]M_Other_Tank(kg) = 12.16 * V_prop(m^3) (for LOX, RP1) [cite: 81]

    Args:
        volume_m3 (float): Volume of the propellant in cubic meters.
        propellant (Literal["LH2", "LOX", "RP1"]): Type of propellant.

    Returns:
        float: Estimated tank mass in kg.
    """
    if propellant == "LH2":
        # [cite_start]M_LH2_Tank(kg) = 9.09 * V_LH2(m^3) [cite: 79]
        return 9.09 * volume_m3
    else:
        # [cite_start]M_Tank(kg) = 12.16 * V_prop(m^3) [cite: 81]
        # [cite_start]This is used for LOX and RP-1 per the regression plot[cite: 51, 56, 71].
        return 12.16 * volume_m3


def estimate_propellant_tank_mass_from_mass(propellant_mass_kg: float,
                                            propellant: Literal["LH2", "LOX", "RP1"]) -> float:
    """
    Estimates propellant tank mass based on propellant mass (alternative MER).

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 7 [cite: 86-92].
    Formulas:
    [cite_start]M_LH2_Tank(kg) = 0.128 * M_LH2(kg) [cite: 88]
    [cite_start]M_LOX_Tank(kg) = 0.0107 * M_LOX(kg) [cite: 90]
    [cite_start]M_RP1_Tank(kg) = 0.0148 * M_RP1(kg) [cite: 92]

    Args:
        propellant_mass_kg (float): Mass of the propellant in kg.
        propellant (Literal["LH2", "LOX", "RP1"]): Type of propellant.

    Returns:
        float: Estimated tank mass in kg.
    """
    if propellant == "LH2":
        # [cite_start]M_LH2_Tank(kg) = 0.128 * M_LH2(kg) [cite: 88]
        return 0.128 * propellant_mass_kg
    elif propellant == "LOX":
        # [cite_start]M_LOX_Tank(kg) = 0.0107 * M_LOX(kg) [cite: 90]
        return 0.0107 * propellant_mass_kg
    elif propellant == "RP1":
        # [cite_start]M_RP1_Tank(kg) = 0.0148 * M_RP1(kg) [cite: 92]
        return 0.0148 * propellant_mass_kg
    else:
        # Fallback for other types, though not specified in this MER
        return 0.01 * propellant_mass_kg


def estimate_cryo_insulation_mass(tank_surface_area_m2: float, propellant: Literal["LH2", "LOX"]) -> float:
    """
    Estimates cryogenic insulation mass based on tank surface area.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 8 [cite: 96-100].
    Formulas:
    [cite_start]M_LH2_Insulation(kg) = 2.88 * A_tank(m^2) [cite: 97, 100]
    [cite_start]M_LOX_Insulation(kg) = 1.123 * A_tank(m^2) [cite: 99, 100]

    Args:
        tank_surface_area_m2 (float): Surface area of the tank in square meters.
        propellant (Literal["LH2", "LOX"]): Type of cryogenic propellant.

    Returns:
        float: Estimated insulation mass in kg.
    """
    if propellant == "LH2":
        # [cite_start]M_LH2_Insulation <kg> = 2.88 <kg/m^2> * A_tank [cite: 97, 100]
        return 2.88 * tank_surface_area_m2
    elif propellant == "LOX":
        # [cite_start]M_LOX_Insulation <kg> = 1.123 <kg/m^2> * A_tank [cite: 99, 100]
        return 1.123 * tank_surface_area_m2
    elif propellant == "Methane":
        # estimated, needs double-check
        return 1.0 * tank_surface_area_m2
    else:
        return 0.0  # No insulation needed for non-cryo


def estimate_pressurized_gas_tank_mass(volume_m3: float, tank_type: Literal["COPV", "Titanium"]) -> float:
    """
    Estimates mass for high-pressure gas tanks (e.g., pressurant tanks).

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 13 [cite: 164-168].
    [cite_start]Based on regression plot Page 12 [cite: 132-161].
    Formulas:
    [cite_start]M_COPV_Tank(kg) = 115.3 * V(m^3) + 3 [cite: 166]
    [cite_start]M_Titanium_Tank(kg) = 299.8 * V(m^3) + 2 [cite: 168]

    Args:
        volume_m3 (float): Volume of the contents in cubic meters.
        tank_type (Literal["COPV", "Titanium"]): Type of high-pressure tank.

    Returns:
        float: Estimated tank mass in kg.
    """
    if tank_type == "COPV":
        # [cite_start]M_COPV Tank(kg) = 115.3 * V_contents(m^3) + 3 [cite: 166]
        return 115.3 * volume_m3 + 3.0
    elif tank_type == "Titanium":
        # [cite_start]M_Titanium_Tank(kg) = 299.8 * V_contents(m^3) + 2 [cite: 168]
        return 299.8 * volume_m3 + 2.0
    else:
        # Default to COPV if type is unknown
        return 115.3 * volume_m3 + 3.0


def estimate_small_liquid_tank_mass(volume_m3: float, tank_type: Literal["Bare", "PMD", "Diaphragm"]) -> float:
    """
    Estimates mass for smaller storable liquid tanks.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 15 [cite: 212-218].
    [cite_start]Based on regression plot Page 14 [cite: 172-210].
    Formulas:
    [cite_start]M_Bare_Tank(kg) = 27.34 * V(m^3) + 2 [cite: 214]
    [cite_start]M_PMD_Tank(kg) = 34.69 * V(m^3) + 3 [cite: 216]
    [cite_start]M_Diaphragm_Tank(kg) = 71.17 * V(m^3) + 3 [cite: 218]

    Args:
        volume_m3 (float): Volume of the contents in cubic meters.
        tank_type (Literal["Bare", "PMD", "Diaphragm"]): Type of small tank.

    Returns:
        float: Estimated tank mass in kg.
    """
    if tank_type == "Bare":
        # [cite_start]M_Bare Tank(kg) = 27.34 * V_contents(m^3) + 2 [cite: 214]
        return 27.34 * volume_m3 + 2.0
    elif tank_type == "PMD":
        # [cite_start]M_PMD Tank(kg) = 34.69 * V_contents(m^3) + 3 [cite: 216]
        return 34.69 * volume_m3 + 3.0
    elif tank_type == "Diaphragm":
        # [cite_start]M_Diaphragm Tank(kg) = 71.17 * V_contents(m^3) + 3 [cite: 218]
        return 71.17 * volume_m3 + 3.0
    else:
        # Default to PMD as a reasonable intermediate
        return 34.69 * volume_m3 + 3.0


def estimate_fairing_mass(fairing_surface_area_m2: float) -> float:
    """
    Estimates fairing/shroud mass based on surface area.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 20 [cite: 279-280].
    [cite_start]Formula: M_fairing(kg) = 4.95 * (A_fairing(m^2))^1.15 [cite: 280]

    Args:
        fairing_surface_area_m2 (float): Surface area of the fairing in square meters.

    Returns:
        float: Estimated fairing mass in kg.
    """
    # [cite_start]Formula from source [cite: 280]
    return 4.95 * (fairing_surface_area_m2 ** 1.15)


def estimate_avionics_mass(vehicle_gross_mass_kg: float) -> float:
    """
    Estimates avionics mass based on vehicle gross mass.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 20[cite: 281, 284].
    [cite_start]Formula: M_avionics(kg) = 10 * (M_o(kg))^0.361 [cite: 284]

    Args:
        vehicle_gross_mass_kg (float): Vehicle gross mass (M_o) in kg.

    Returns:
        float: Estimated avionics mass in kg.
    """
    # [cite_start]Formula from source [cite: 284]
    return 10.0 * (vehicle_gross_mass_kg ** 0.361)


def estimate_wiring_mass(vehicle_gross_mass_kg: float, vehicle_length_m: float) -> float:
    """
    Estimates wiring mass based on vehicle gross mass and length.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 20[cite: 283, 285].
    [cite_start]Formula: M_wiring(kg) = 1.058 * sqrt(M_o(kg)) * l^0.25 [cite: 285]

    Args:
        vehicle_gross_mass_kg (float): Vehicle gross mass (M_o) in kg.
        vehicle_length_m (float): Vehicle total length (l) in meters.

    Returns:
        float: Estimated wiring mass in kg.
    """
    # [cite_start]Formula from source [cite: 285]
    return 1.058 * (vehicle_gross_mass_kg ** 0.5) * (vehicle_length_m ** 0.25)


def estimate_srm_casing_mass(propellant_mass_kg: float) -> float:
    """
    Estimates Solid Rocket Motor (SRM) casing mass from propellant mass.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 25[cite: 359, 361].
    [cite_start]Formula: M_Motor_Casing = 0.135 * M_propellants [cite: 361]

    Args:
        propellant_mass_kg (float): Mass of the solid propellant in kg.

    Returns:
        float: Estimated motor casing mass in kg.
    """
    # [cite_start]Formula from source [cite: 361]
    return 0.135 * propellant_mass_kg


def estimate_gimbal_mass(engine_thrust_N: float, chamber_pressure_Pa: float) -> float:
    """
    Estimates gimbal mass for a single engine.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 26 [cite: 368-369].
    [cite_start]Formula: M_Gimbals(kg) = 237.8 * [T(N) / P_c(Pa)]^0.9375 [cite: 369]

    [cite_start]Note: The reference slide [cite: 369] uses P_0(Pa), but the example calculation
    on [cite_start]page 28 [cite: 391] (using Pc=6.897e6 Pa) and mass summary on [cite_start]page 30 [cite: 406]
    (Total Gimbals = 81 kg for 6 engines) confirm this refers to
    chamber pressure (P_c) and that the mass is per engine.
    (6 engines * 13.5 kg/eng ~= 81 kg total).

    Args:
        engine_thrust_N (float): Thrust of a single engine in Newtons.
        chamber_pressure_Pa (float): Engine chamber pressure in Pascals.

    Returns:
        float: Estimated gimbal mass in kg.
    """
    if chamber_pressure_Pa <= 0:
        return 0.0

    # [cite_start]Ratio T(N) / P_0(Pa) from formula [cite: 369]
    ratio = engine_thrust_N / chamber_pressure_Pa
    # [cite_start]Full formula from source [cite: 369]
    return 237.8 * (ratio ** 0.9375)


def estimate_gimbal_torque(engine_thrust_N: float, chamber_pressure_Pa: float) -> float:
    """
    Estimates gimbal torque for a single engine. (Note: Not a mass relation).

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 26 [cite: 370-371].
    [cite_start]Formula: Tau_Gimbals(N*m) = 990,000 * [T(N) / P_c(Pa)]^1.25 [cite: 371]

    Args:
        engine_thrust_N (float): Thrust of a single engine in Newtons.
        chamber_pressure_Pa (float): Engine chamber pressure in Pascals.

    Returns:
        float: Estimated gimbal torque in N*m.
    """
    if chamber_pressure_Pa <= 0:
        return 0.0

    # [cite_start]Ratio T(N) / P_0(Pa) from formula [cite: 371]
    ratio = engine_thrust_N / chamber_pressure_Pa
    # [cite_start]Full formula from source [cite: 371]
    return 990_000 * (ratio ** 1.25)


def calculate_sphere_area(propellant_mass: float, propellant_density: float):
    """
    Calculate the estimated area of the tank surface based on propellant mass and density.
    This is a helper function based on the logic from the article.

    Reference:
    [cite_start]Logic derived from examples on Page 9 [cite: 107-109] [cite_start]and Page 10 [cite: 117-119].
    [cite_start]Assumes a spherical tank: V = M / rho [cite: 107, 117]
    [cite_start]r = (V / (4pi/3))^(1/3) [cite: 108, 118]
    [cite_start]A = 4 * pi * r^2 [cite: 109, 119]

    Args:
        propellant_mass (float): Mass of the propellant in kg.
        propellant_density (float): Density of the propellant in kg/m^3.

    Returns:
        float: Estimated spherical surface area in m^2.
    """
    if propellant_density == 0:
        return 0.0
    # [cite_start]V = M / rho [cite: 107, 117]
    volume = propellant_mass / propellant_density
    # [cite_start]r = (V / (4pi/3))^(1/3) [cite: 108, 118]
    sphere_radius = (volume / (4 * math.pi / 3)) ** (1 / 3)
    # [cite_start]A = 4 * pi * r^2 [cite: 109, 119]
    return 4 * math.pi * (sphere_radius ** 2)


def calculate_cone_area(radius: float, height: float) -> float:
    """
    Calculates the surface area of a cone (excluding the base).
    A = pi * r * sqrt(r^2 + h^2)

    Reference:
    Mass Estimating Relations (Akin, ENAE 791), Page 21[cite: 290].

    Args:
        radius (float): Base radius of the cone (r).
        height (float): Height of the cone (h).

    Returns:
        float: Surface area in m^2.
    """
    return math.pi * radius * math.sqrt(radius**2 + height**2)

def calculate_cylinder_area(radius: float, height: float) -> float:
    """
    Calculates the surface area of a cylinder's side wall.
    A = 2 * pi * r * h

    Reference:
    Mass Estimating Relations (Akin, ENAE 791), Page 21[cite: 295].

    Args:
        radius (float): Radius of the cylinder (r).
        height (float): Height of the cylinder (h).

    Returns:
        float: Surface area in m^2.
    """
    return 2 * math.pi * radius * height

def calculate_frustum_area(radius1: float, radius2: float, height: float) -> float:
    """
    Calculates the surface area of a cone frustum.
    A = pi * (r1 + r2) * sqrt((r1 - r2)^2 + h^2)

    Reference:
    Mass Estimating Relations (Akin, ENAE 791), Page 21[cite: 292].

    Args:
        radius1 (float): Radius of the first base (r1).
        radius2 (float): Radius of the second base (r2).
        height (float): Height of the frustum (h).

    Returns:
        float: Surface area in m^2.
    """
    return math.pi * (radius1 + radius2) * math.sqrt((radius1 - radius2)**2 + height**2)