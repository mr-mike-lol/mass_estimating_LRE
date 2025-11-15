# models/mota_schlingloff.py
"""
Implements the Mota/Schlingloff (2005) analytical-statistical mass model
for liquid rocket engines, as presented in Mota et al. (2018).

Approach Overview:
This model estimates the total dry mass of a pump-fed liquid rocket
engine by breaking it down into five key components:
1.  Turbopump (m_tp)
2.  Valves (m_valve)
3.  Injector (m_inj)
4.  Combustion Chamber (m_cc)
5.  Nozzle & Cooling (m_nc)

Each component's mass is estimated using a separate empirical formula
derived from key engine parameters (Thrust, Chamber Pressure) and
technology coefficients. The sum of these component masses
is then multiplied by a final correction factor (1.34) to arrive
at the total engine mass.

This approach is described in the paper as an "analytical/statistical model"
because it "considers not only statistical data but also physical
relationships".

Pros:
+   **Detailed:** More detailed than a simple top-level regression (like
    Zandbergen), as it is "sufficiently detailed when the influence
    of the engine parameters on the engine mass... are aim of study".
+   **Technology-Aware:** The model accounts for key technology choices
    through coefficients for:
    * Propellant energy (high vs. low)
    * Presence of boost pumps
    * Cooling type (regenerative vs. dump)
+   **Component Breakdown:** Provides a mass estimate for individual
    components, not just a total number.

Cons:
-   **Empirical:** The component formulas are still regressions
    based on "historical and empirical data".
-   **Known Gaps:** The original Schlingloff (2005) model (Eq. 21)
    does not explicitly account for the engine's mixture ratio (O/F).
-   **Unit-Specific:** The formulas are rigid and require specific
    units (kN and bar) to function correctly.
"""

from models.common_params import EngineParams, PropellantType
from models.base import BaseEngineModel, ModelResult
from typing import Dict

# For standalone running
import inspect
import vehicle_definitions

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
    if (thrust_kN * pc_bar) < 0:
        return 0.0  # Avoid math domain error with negative inputs
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
    if (thrust_kN * pc_bar) < 0:
        return 0.0  # Avoid math domain error
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
    if thrust_kN < 0:
        return 0.0
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
    if thrust_kN < 0:
        return 0.0
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


# --- Public Interface Class ---

class MotaSchlingloffModel(BaseEngineModel):
    """
    Implements the Mota/Schlingloff (2005) analytical-statistical mass model
    as a class inheriting from BaseEngineModel.

    This model estimates the total dry mass of a pump-fed liquid rocket
    engine by breaking it down into five key components.

    Technology flags (e.g., boost pumps, cooling type) are set during
    instantiation.
    """

    def __init__(self,
                 has_boost_pumps: bool = False,
                 is_regen_cooled: bool = True):
        """
        Initializes the model with technology-specific flags.

        Args:
            has_boost_pumps (bool): Flag for boost pumps (like on SSME).
                                    Affects C_tp. Defaults to False.
            is_regen_cooled (bool): Flag for regenerative cooling.
                                    Affects C_nozzle. Defaults to True.
        """
        self.has_boost_pumps = has_boost_pumps
        self.is_regen_cooled = is_regen_cooled

    @property
    def model_name(self) -> str:
        """Returns the unique, human-readable name of the model."""
        return "Mota/Schlingloff (2005)"

    def estimate_mass(self, params: EngineParams) -> ModelResult:
        """
        Estimates total engine mass using the Schlingloff (2005) component model.

        This function implements the component models described in Mota et al. (2018)
        [cite: 228, 230, 231, 232], which are then summed and multiplied by a
        correction factor.

        The model requires specific units:
        - Thrust in KiloNewtons (kN)
        - Chamber Pressure in Bar

        Reference:
        Modeling and Analysis of a LOX/Ethanol... (Mota et al., 2018), p. 6-7.
        [cite: 216-248]

        Args:
            params (EngineParams): The engine's parameter object.

        Returns:
            ModelResult: A TypedDict containing 'total_mass_kg' and a 'components_kg' dict.

        Raises:
            ValueError: If chamber pressure is not a positive number.
        """

        # --- 1. Convert units and validate inputs ---
        thrust_kN = params.thrust_vac_N / N_PER_KN
        pc_bar = params.chamber_pressure_Pa / PA_PER_BAR

        if pc_bar <= 0:
            raise ValueError("Chamber pressure must be > 0 for Mota/Schlingloff calculation.")

        # --- 2. Determine model constants based on inputs ---

        # C_propellant: 0.19 (high energetic), 0.11 (low energetic) [cite: 231]
        if params.propellant_type == "LOX/LH2":
            c_propellant = 0.19
        else:
            c_propellant = 0.11  # For LOX/RP1, LOX/LCH4, Storable

        # C_tp: 0.5 (boost-pumps), 1.0 (no boost-pumps) [cite: 231]
        c_tp = 0.5 if self.has_boost_pumps else 1.0

        # C_nozzle: 1.0 (regenerative), 0.0 (dump cooling) [cite: 232]
        c_nozzle = 1.0 if self.is_regen_cooled else 0.0

        # --- 3. Calculate mass for each component ---
        components: Dict[str, float] = {}
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

        # --- 5. Return standardized ModelResult ---
        return ModelResult(
            total_mass_kg=total_mass,
            components_kg=components,
            notes={
                "correction_factor": correction_factor,
                "c_propellant": c_propellant,
                "c_tp": c_tp,
                "c_nozzle": c_nozzle,
                "has_boost_pumps": self.has_boost_pumps,
                "is_regen_cooled": self.is_regen_cooled
            }
        )


if __name__ == "__main__":
    """
    Standalone runner for the Mota/Schlingloff model.

    Allows testing a specific engine or all engines defined
    in vehicle_definitions.py.
    """

    # --- 1. CONFIGURE TEST RUN ---
    # Set this variable to a specific engine key (e.g., "ssme", "rd120") to test ONLY that one engine.
    # Set to None to test ALL available engines.
    SPECIFIC_ENGINE_KEY_TO_TEST = "default"
    # Example to test only SSME:
    # SPECIFIC_ENGINE_KEY_TO_TEST = "ssme"
    # Example to test all:
    # SPECIFIC_ENGINE_KEY_TO_TEST = None

    # --- 2. Find all available engine getter functions ---
    available_engines = {}
    for name, func in inspect.getmembers(vehicle_definitions, inspect.isfunction):
        if name.startswith("get_") and name.endswith("_engine"):
            # Store the key (e.g., 'ssme') and the function itself
            key = name.replace("get_", "").replace("_engine", "")
            available_engines[key] = func

    print(f"[__main__] Found {len(available_engines)} engines: {list(available_engines.keys())}")
    print("-" * 70)

    # --- 3. Define which engines to test (based on config) ---
    keys_to_test = []
    if SPECIFIC_ENGINE_KEY_TO_TEST:
        # --- IF (Specific Engine) ---
        if SPECIFIC_ENGINE_KEY_TO_TEST in available_engines:
            print(f"[__main__] Mode: Testing SPECIFIC engine: '{SPECIFIC_ENGINE_KEY_TO_TEST}'")
            keys_to_test = [SPECIFIC_ENGINE_KEY_TO_TEST]
        else:
            # Error: The specified key doesn't exist
            print(f"[__main__] ERROR: Specific engine key '{SPECIFIC_ENGINE_KEY_TO_TEST}' not found.")
            print(f"[__main__] Available keys are: {list(available_engines.keys())}")
    else:
        # --- ELSE (All Engines) ---
        print(f"[__main__] Mode: Testing ALL {len(available_engines)} engines.")
        keys_to_test = list(available_engines.keys())

    # --- 4. Loop and run analysis for each engine ---
    if not keys_to_test:
        if SPECIFIC_ENGINE_KEY_TO_TEST:
            # This branch is reached if a specific key was given but not found
            print("[__main__] No engines to test.")
        else:
            # This branch is reached if no engines were found at all
            print("[__main__] No engine getter functions found in vehicle_definitions.py.")

    for engine_key in keys_to_test:
        print(f"\n[__main__] Loading parameters for '{engine_key}'...")
        engine_getter = available_engines[engine_key]

        try:
            engine_params = engine_getter()

            # --- Mota/Schlingloff Specific Logic ---
            # We intelligently set the 'has_boost_pumps' flag based on
            # the engine cycle. Staged Combustion (SC) engines
            # (like SSME) often have boost pumps.
            has_boost_pumps = (engine_params.cycle_type == "SC")

            print(f"[__main__] Configuring model with: has_boost_pumps={has_boost_pumps}")

            # Instantiate the model with this configuration
            model = MotaSchlingloffModel(
                has_boost_pumps=has_boost_pumps,
                is_regen_cooled=True  # Assume regen cooled for all tests
            )

            # Call the INHERITED run_single_engine_analysis method
            # This method is defined in BaseEngineModel and handles
            # the calculation, timing, and pretty-printing.
            model.run_single_engine_analysis(engine_params)

        except Exception as e:
            print(f"\n[__main__] ERROR processing engine '{engine_key}': {e}")
            # Optional: uncomment to raise the full error
            # raise e

    print("-" * 70)
    print("[__main__] Standalone analysis complete.")