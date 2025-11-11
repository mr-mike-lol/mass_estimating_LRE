# models/akin_mers.py

import math
from typing import Literal, Dict, Tuple, Any
from vehicle_definitions import DENSITY_RP1, DENSITY_LH2, DENSITY_LOX, G0


def estimate_engine_mass(thrust_N: float, expansion_ratio: float) -> float:
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

    return total_mass


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
    This is the MER used in the SSTO example calcs.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 7 [cite: 86-92].
    Formulas:
    [cite_start]M_LH2_Tank(kg) = 0.128 * M_LH2(kg) [cite: 88, 116]
    [cite_start]M_LOX_Tank(kg) = 0.0107 * M_LOX(kg) [cite: 90, 106]
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
        # Raise error other types, though not specified in this MER
        raise ValueError(f"Unknown propellant type for mass-based tank MER: {propellant}")


def estimate_cryo_insulation_mass(tank_surface_area_m2: float, propellant: Literal["LH2", "LOX"]) -> float:
    """
    Estimates cryogenic insulation mass based on tank surface area.

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 8 [cite: 96-100].
    Formulas:
    [cite_start]M_LH2_Insulation(kg) = 2.88 * A_tank(m^2) [cite: 97, 100, 120]
    [cite_start]M_LOX_Insulation(kg) = 1.123 * A_tank(m^2) [cite: 99, 100, 110]

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
        # estimated, needs double-check; cz it's cryogenic
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
    if fairing_surface_area_m2 <= 0:
        return 0.0
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
    if vehicle_gross_mass_kg <= 0 or vehicle_length_m <= 0:
        return 0.0
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


# --- Geometry Helper Functions ---
# These functions implement the geometry formulas used in the SSTO example
# and can be used for the "optional" fairing calculations.

def _calculate_sphere_geom(propellant_mass: float, propellant_density: float) -> Tuple[float, float]:
    """
    (Internal) Calculates the radius and surface area of a spherical tank.
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
        Tuple[float, float]: (radius_m, area_m2)
    """
    if propellant_density == 0 or propellant_mass == 0:
        return 0.0, 0.0
    # [cite_start]V = M / rho [cite: 107, 117]
    volume = propellant_mass / propellant_density
    # [cite_start]r = (V / (4pi/3))^(1/3) [cite: 108, 118]
    sphere_radius = (volume / (4 * math.pi / 3)) ** (1 / 3)
    # [cite_start]A = 4 * pi * r^2 [cite: 109, 119]
    area = 4 * math.pi * (sphere_radius ** 2)
    return sphere_radius, area


def calculate_sphere_area(propellant_mass: float, propellant_density: float) -> float:
    """
    Public helper to calculate just the surface area of a spherical tank.
    See _calculate_sphere_geom for details.
    """
    _, area = _calculate_sphere_geom(propellant_mass, propellant_density)
    return area


def calculate_cone_area(radius: float, height: float) -> float:
    """
    Calculates the surface area of a cone (excluding the base).
    A = pi * r * sqrt(r^2 + h^2)

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 21[cite: 290].

    Args:
        radius (float): Base radius of the cone (r).
        height (float): Height of the cone (h).

    Returns:
        float: Surface area in m^2.
    """
    if radius <= 0 or height <= 0:
        return 0.0
    return math.pi * radius * math.sqrt(radius ** 2 + height ** 2)


def calculate_cylinder_area(radius: float, height: float) -> float:
    """
    Calculates the surface area of a cylinder's side wall.
    A = 2 * pi * r * h

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 21[cite: 295].

    Args:
        radius (float): Radius of the cylinder (r).
        height (float): Height of the cylinder (h).

    Returns:
        float: Surface area in m^2.
    """
    if radius <= 0 or height <= 0:
        return 0.0
    return 2 * math.pi * radius * height


def calculate_frustum_area(radius1: float, radius2: float, height: float) -> float:
    """
    Calculates the surface area of a cone frustum.
    A = pi * (r1 + r2) * sqrt((r1 - r2)^2 + h^2)

    Reference:
    [cite_start]Mass Estimating Relations (Akin, ENAE 791), Page 21[cite: 292].

    Args:
        radius1 (float): Radius of the first base (r1).
        radius2 (float): Radius of the second base (r2).
        height (float): Height of the frustum (h).

    Returns:
        float: Surface area in m^2.
    """
    if height <= 0:
        return 0.0
    return math.pi * (radius1 + radius2) * math.sqrt((radius1 - radius2) ** 2 + height ** 2)


# --- New Example Runner Function ---

def get_akin_ssto_default_params() -> Dict[str, Any]:
    """
    Returns the default configuration for the SSTO 1st Pass example
    from the "Mass Estimating Relations" PDF.
    """
    return {
        # [cite_start]Page 3 [cite: 24]
        "payload_kg": 5000.0,
        # [cite_start]Page 3 [cite: 23]
        "delta_v_ms": 9200.0,
        # [cite_start]Page 3 [cite: 26]
        "isp_s": 430.0,
        # [cite_start]Page 3 [cite: 28]
        "initial_delta": 0.08,
        # [cite_start]Page 4 [cite: 40]
        "mixture_ratio": 6.0,
        # [cite_start]Page 27 [cite: 375]
        "initial_twr": 1.3,
        # [cite_start]Page 27 [cite: 377]
        "num_engines": 6,
        # [cite_start]Page 27 [cite: 381]
        "chamber_pressure_Pa": 6.897e6,  # 1000 psi
        # [cite_start]Page 27 [cite: 383]
        "expansion_ratio": 30.0,
        # [cite_start]Page 22 [cite: 313]
        "payload_fairing_height_m": 7.0,
        # [cite_start]Page 22 [cite: 316]
        "intertank_fairing_height_m": 7.0,
        # [cite_start]Page 22 [cite: 319]
        "aft_fairing_height_m": 7.0,
    }


def run_akin_ssto_example(params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Runs the 1st Pass SSTO mass budget analysis based on the
    "Mass Estimating Relations" (Akin, ENAE 791) PDF.

    This function follows the logic from Page 3 to Page 30.

    Args:
        params (Dict[str, Any]): A dictionary of parameters,
                                 see `get_akin_ssto_default_params` for keys.

    Returns:
        Dict[str, Any]: A dictionary containing the detailed mass budget.
    """
    # Unpack parameters
    p = params

    # --- Step 1: Vehicle-Level 1st Pass (Page 3) ---
    # [cite_start]Page 3 [cite: 27]
    Ve = p['isp_s'] * G0
    # [cite_start]Page 3 [cite: 31]
    r = math.exp(-p['delta_v_ms'] / Ve)
    # [cite_start]Page 3 [cite: 32]
    lambda_payload = r - p['initial_delta']
    if lambda_payload <= 0:
        raise ValueError("Infeasible mission: delta (inert fraction) is too high for this delta-V.")
    # [cite_start]Page 3 [cite: 33]
    M_o = p['payload_kg'] / lambda_payload
    # [cite_start]Page 3 [cite: 34]
    M_i_initial_guess = p['initial_delta'] * M_o
    # [cite_start]Page 3 [cite: 35]
    M_p = M_o * (1 - r)

    # --- Step 2: Propellant Calcs (Page 4) ---
    # [cite_start]Page 4 [cite: 40, 42]
    M_lh2 = M_p / (1 + p['mixture_ratio'])
    # [cite_start]Page 4 [cite: 40, 43]
    M_lox = M_p * p['mixture_ratio'] / (1 + p['mixture_ratio'])

    # --- Step 3: Tank Mass (Page 7, 9, 10) ---
    # [cite_start]Page 7 [cite: 90][cite_start], Page 9 [cite: 106]
    M_lox_tank = estimate_propellant_tank_mass_from_mass(M_lox, "LOX")
    # [cite_start]Page 7 [cite: 88][cite_start], Page 10 [cite: 116]
    M_lh2_tank = estimate_propellant_tank_mass_from_mass(M_lh2, "LH2")

    # --- Step 4: Insulation Mass (Page 8-10) ---
    # [cite_start]Page 9 [cite: 108][cite_start], Page 10 [cite: 118]
    r_lox_tank, A_lox_tank = _calculate_sphere_geom(M_lox, DENSITY_LOX)
    r_lh2_tank, A_lh2_tank = _calculate_sphere_geom(M_lh2, DENSITY_LH2)
    # [cite_start]Page 8 [cite: 99][cite_start], Page 9 [cite: 110]
    M_lox_insulation = estimate_cryo_insulation_mass(A_lox_tank, "LOX")
    # [cite_start]Page 8 [cite: 97][cite_start], Page 10 [cite: 120]
    M_lh2_insulation = estimate_cryo_insulation_mass(A_lh2_tank, "LH2")

    # --- Step 5: Fairing Mass (Page 20-23) ---
    # [cite_start]Page 21-22 [cite: 290, 313, 322]
    A_payload_fairing = calculate_cone_area(r_lox_tank, p['payload_fairing_height_m'])
    # [cite_start]Page 21-22 [cite: 292, 316, 322, 323]
    A_intertank_fairing = calculate_frustum_area(r_lh2_tank, r_lox_tank, p['intertank_fairing_height_m'])
    # [cite_start]Page 21-22 [cite: 295, 319, 323]
    A_aft_fairing = calculate_cylinder_area(r_lh2_tank, p['aft_fairing_height_m'])
    # [cite_start]Page 20 [cite: 280][cite_start], Page 23 [cite: 340]
    M_payload_fairing = estimate_fairing_mass(A_payload_fairing)
    # [cite_start]Page 20 [cite: 280][cite_start], Page 23 [cite: 342]
    M_intertank_fairing = estimate_fairing_mass(A_intertank_fairing)
    # [cite_start]Page 20 [cite: 280][cite_start], Page 23 [cite: 344]
    M_aft_fairing = estimate_fairing_mass(A_aft_fairing)

    # --- Step 6: Avionics & Wiring (Page 20, 24) ---
    # [cite_start]Page 24 [cite: 353] - Approx. length based on 7+7+7 fairing sections
    vehicle_length_l = p['payload_fairing_height_m'] + p['intertank_fairing_height_m'] + p['aft_fairing_height_m']
    # [cite_start]Page 20 [cite: 284][cite_start], Page 24 [cite: 351]
    M_avionics = estimate_avionics_mass(M_o)
    # [cite_start]Page 20 [cite: 285][cite_start], Page 24 [cite: 353]
    M_wiring = estimate_wiring_mass(M_o, vehicle_length_l)

    # --- Step 7: Propulsion (Page 25-28) ---
    # [cite_start]Page 27-28 [cite: 375, 390]
    Total_Thrust_T = M_o * G0 * p['initial_twr']
    T_per_engine = Total_Thrust_T / p['num_engines']
    # [cite_start]Page 25 [cite: 358, 360][cite_start], Page 28 [cite: 392]
    M_per_engine = estimate_engine_mass(T_per_engine, p['expansion_ratio'])
    M_engines_total = M_per_engine * p['num_engines']
    # [cite_start]Page 25 [cite: 363][cite_start], Page 28 [cite: 394]
    M_thrust_structure_per_eng = estimate_thrust_structure_mass(T_per_engine)
    M_thrust_structure_total = M_thrust_structure_per_eng * p['num_engines']
    # [cite_start]Page 26 [cite: 369]
    M_gimbal_per_engine = estimate_gimbal_mass(T_per_engine, p['chamber_pressure_Pa'])
    M_gimbals_total = M_gimbal_per_engine * p['num_engines']

    # --- Step 8: Mass Summary (Page 30) ---
    components = {
        "LOX Tank": M_lox_tank,
        "LH2 Tank": M_lh2_tank,
        "LOX Insulation": M_lox_insulation,
        "LH2 Insulation": M_lh2_insulation,
        "Payload Fairing": M_payload_fairing,
        "Intertank Fairing": M_intertank_fairing,
        "Aft Fairing": M_aft_fairing,
        "Engines": M_engines_total,
        "Thrust Structure": M_thrust_structure_total,
        "Gimbals": M_gimbals_total,
        "Avionics": M_avionics,
        "Wiring": M_wiring,
    }
    # [cite_start]Page 30 [cite: 406]
    M_i_calculated = sum(components.values())
    margin = (M_i_initial_guess - M_i_calculated) / M_i_initial_guess

    results = {
        "params": p,
        "initial_calcs": {
            "M_o_kg": M_o,
            "M_i_initial_guess_kg": M_i_initial_guess,
            "M_propellant_kg": M_p,
            "M_payload_kg": p['payload_kg'],
            "mass_ratio_r": r,
            "payload_fraction_lambda": lambda_payload,
        },
        "propulsion": {
            "total_thrust_N": Total_Thrust_T,
            "thrust_per_engine_N": T_per_engine,
            "mass_per_engine_kg": M_per_engine,
            "gimbal_mass_per_engine_kg": M_gimbal_per_engine,
            "thrust_structure_mass_per_engine_kg": M_thrust_structure_per_eng,
        },
        "geometry": {
            "lox_tank_radius_m": r_lox_tank,
            "lox_tank_area_m2": A_lox_tank,
            "lh2_tank_radius_m": r_lh2_tank,
            "lh2_tank_area_m2": A_lh2_tank,
            "vehicle_length_approx_m": vehicle_length_l,
        },
        "mass_budget": {
            "components_kg": components,
            "total_inert_mass_calculated_kg": M_i_calculated,
            "total_inert_mass_initial_guess_kg": M_i_initial_guess,
            "design_margin_percent": margin * 100.0,
        }
    }
    return results


def _print_ssto_results(results: Dict[str, Any]):
    """Helper to format and print the SSTO analysis results."""

    print("=" * 70)
    print("🚀 AKIN SSTO 1st PASS ANALYSIS RESULTS")
    print(f"(Based on ENAE 791, 1st Pass, Pages 3-30)")
    print("=" * 70)

    print("\n--- Initial Vehicle Sizing (from Page 3) ---")
    ic = results['initial_calcs']
    print(f"  Gross Mass (M_o):         {ic['M_o_kg']:12,.1f} kg")
    print(f"  Propellant Mass (M_p):    {ic['M_propellant_kg']:12,.1f} kg")
    print(f"  Initial Inert Guess (M_i):{ic['M_i_initial_guess_kg']:12,.1f} kg")
    print(f"  Payload Mass (M_l):       {ic['M_payload_kg']:12,.1f} kg")

    print("\n--- Propulsion System (from Page 27-28) ---")
    pr = results['propulsion']
    p = results['params']
    print(f"  Num. Engines:             {p['num_engines']}")
    print(f"  Total Thrust:             {pr['total_thrust_N'] / 1e6:12.2f} MN")
    print(f"  Thrust per Engine:        {pr['thrust_per_engine_N'] / 1e3:12.1f} kN")
    print(f"  Mass per Engine:          {pr['mass_per_engine_kg']:12.1f} kg")

    print("\n--- Calculated Mass Budget (from Page 30) ---")
    mb = results['mass_budget']
    components = mb['components_kg']

    # [cite_start]Replicating table from Page 30 [cite: 406]
    print(f"  {'Component':<18} | {'Calculated (kg)':>15} | {'PDF Ref (kg)':>12}")
    print(f"  {'-' * 18: <18} | {'-' * 15:>15} | {'-' * 12:>12}")
    # Helper to map component names to PDF values
    pdf_ref = {
        "LOX Tank": 1245, "LH2 Tank": 2482, "LOX Insulation": 119,
        "LH2 Insulation": 586, "Payload Fairing": 645, "Intertank Fairing": 1626,
        "Aft Fairing": 1905, "Engines": 2236, "Thrust Structure": 497,
        "Gimbals": 81, "Avionics": 744, "Wiring": 886
    }
    for name, mass in components.items():
        ref_val_str = f"{pdf_ref.get(name, 0):,.0f}"
        print(f"  {name:<18} | {mass:15,.1f} | {ref_val_str:>12}")

    print(f"  {'-' * 18: <18} | {'-' * 15:>15} | {'-' * 12:>12}")
    print(f"  {'Total Inert Mass':<18} | {mb['total_inert_mass_calculated_kg']:15,.1f} | {13052:>12}")
    print(f"  {'Initial Guess':<18} | {mb['total_inert_mass_initial_guess_kg']:15,.1f} | {12240:>12}")

    print("\n--- FINAL DESIGN MARGIN ---")
    print(f"  Calculated Margin:   {mb['design_margin_percent']:15.2f} %")
    # [cite_start]Page 30 [cite: 406]
    print(f"  PDF Reference Margin: {-6.22:15.2f} %")
    print("=" * 70)


if __name__ == "__main__":
    """
    This block allows the file to be run directly (e.g., `python models/akin_mers.py`)
    to execute the SSTO 1st Pass analysis using the default parameters
    from the PDF.
    """

    print("Running SSTO 1st Pass Example directly from akin_mers.py...")

    # 1. Get the specific config from the PDF
    default_params = get_akin_ssto_default_params()

    # 2. Run the analysis
    try:
        results = run_akin_ssto_example(default_params)

        # 3. Print the formatted results
        _print_ssto_results(results)

    except Exception as e:
        print(f"\nAn error occurred during the analysis: {e}")