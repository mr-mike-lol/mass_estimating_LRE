# models/mota_schlingloff.py

from .common_params import EngineParams, PropellantType

# --- Unit Conversion Constants ---
N_PER_KN = 1000.0
PA_PER_BAR = 100_000.0


# --- Private Component Functions ---

def _estimate_turbopump_mass(thrust_kN: float, pc_bar: float,
                             c_propellant: float, c_tp: float) -> float:
    """
    Estimates Turbopump (m_tp) mass.

    Reference:
    Modeling and Analysis of a LOX/Ethanol... (Mota et al., 2018), Eq. 21 text.
    Citing Schlingloff (2005).
    Formula: m_tp = C_propellant * C_tp * (F_kN * p_c_bar)^0.71
    [cite: 230]
    """
    return c_propellant * c_tp * ((thrust_kN * pc_bar) ** 0.71)


def _estimate_valve_mass(thrust_kN: float, pc_bar: float) -> float:
    """
    Estimates Valve (m_valve) mass.

    Reference:
    Modeling and Analysis of a LOX/Ethanol... (Mota et al., 2018), Eq. 21 text.
    Citing Schlingloff (2005).
    Formula: m_valve = 0.02 * (F_kN * p_c_bar)^0.71
    [cite: 230]
    """
    return 0.02 * ((thrust_kN * pc_bar) ** 0.71)


def _estimate_injector_mass(thrust_kN: float) -> float:
    """
    Estimates Injector (m_inj) mass.

    Reference:
    Modeling and Analysis of a LOX/Ethanol... (Mota et al., 2018), Eq. 21 text.
    Citing Schlingloff (2005).
    Formula: m_inj = 0.25 * (F_kN)^0.85
    [cite: 230]
    """
    return 0.25 * (thrust_kN ** 0.85)


def _estimate_chamber_mass(thrust_kN: float) -> float:
    """
    Estimates Combustion Chamber (m_cc) mass.

    Reference:
    Modeling and Analysis of a LOX/Ethanol... (Mota et al., 2018), Eq. 21 text.
    Citing Schlingloff (2005).
    Formula: m_cc = 0.75 * (F_kN)^0.85
    [cite: 230]
    """
    return 0.75 * (thrust_kN ** 0.85)


def _estimate_nozzle_mass(thrust_kN: float, pc_bar: float, c_nozzle: float) -> float:
    """
    Estimates Nozzle/Cooling (m_nc) mass.

    Reference:
    Modeling and Analysis of a LOX/Ethanol... (Mota et al., 2018), Eq. 21 text.
    Citing Schlingloff (2005).
    Formula: m_nc = F_kN * (0.00225 * C_nozzle + (0.225 - 0.075 * C_nozzle)) / p_c_bar
    [cite: 230]
    """
    if pc_bar == 0:
        return 0.0

    term1 = 0.00225 * c_nozzle
    term2 = 0.225 - (0.075 * c_nozzle)

    return thrust_kN * (term1 + term2) / pc_bar


# --- Public Interface Function ---

def estimate_total_engine_mass(params: EngineParams,
                               has_boost_pumps: bool = False,
                               is_regen_cooled: bool = True) -> dict:
    """
    Estimates total engine mass using the Schlingloff (2005) component model.

    This function implements the component models described in Mota et al. (2018)
    [cite: 228, 230, 231, 232], which are then summed and multiplied by a correction factor.

    The model requires specific units:
    - Thrust in KiloNewtons (kN)
    - Chamber Pressure in Bar

    This function returns a dictionary containing the total estimated mass
    and the mass breakdown of the components.

    Reference:
    Modeling and Analysis of a LOX/Ethanol... (Mota et al., 2018), p. 6-7.
    [cite: 216-248]

    Args:
        params (EngineParams): The engine's parameter object.
        has_boost_pumps (bool): Flag for boost pumps (like on SSME).
                                Affects C_tp. Defaults to False.
        is_regen_cooled (bool): Flag for regenerative cooling.
                                Affects C_nozzle. Defaults to True.

    Returns:
        dict: A dictionary containing 'total_mass_kg' and a 'components' dict.
    """

    # --- 1. Convert units to those required by the model ---
    thrust_kN = params.thrust_vac_N / N_PER_KN
    pc_bar = params.chamber_pressure_Pa / PA_PER_BAR

    # --- 2. Determine model constants based on inputs ---

    # C_propellant: 0.19 (high energetic), 0.11 (low energetic) [cite: 231]
    if params.propellant_type == "LOX/LH2":
        c_propellant = 0.19
    else:
        c_propellant = 0.11  # For LOX/RP1, LOX/LCH4, Storable

    # C_tp: 0.5 (boost-pumps), 1.0 (no boost-pumps) [cite: 231]
    c_tp = 0.5 if has_boost_pumps else 1.0

    # C_nozzle: 1.0 (regenerative), 0.0 (dump cooling) [cite: 232]
    c_nozzle = 1.0 if is_regen_cooled else 0.0

    # --- 3. Calculate mass for each component ---
    components = {}
    components['turbopump'] = _estimate_turbopump_mass(thrust_kN, pc_bar, c_propellant, c_tp)
    components['valves'] = _estimate_valve_mass(thrust_kN, pc_bar)
    components['injector'] = _estimate_injector_mass(thrust_kN)
    components['chamber'] = _estimate_chamber_mass(thrust_kN)
    components['nozzle_cooling'] = _estimate_nozzle_mass(thrust_kN, pc_bar, c_nozzle)

    # --- 4. Sum components and apply correction factor ---

    # This model uses the original Schlingloff (2005) factor
    # The paper also explores a hybrid model (Eq. 23), but that
    # requires a different m_tp (Eq. 22) [cite: 237] for which we lack inputs (Pump Power).
    # We are implementing the self-contained Schlingloff model (Eq. 21).

    correction_factor = 1.34  # From Eq. 21

    component_sum = sum(components.values())
    total_mass = correction_factor * component_sum

    return {
        "total_mass_kg": total_mass,
        "components_kg": components,
        "notes": {
            "model": "Mota/Schlingloff (2005)",
            "correction_factor": correction_factor,
            "c_propellant": c_propellant,
            "c_tp": c_tp,
            "c_nozzle": c_nozzle
        }
    }