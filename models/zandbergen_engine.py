# models/zandbergen_engine.py
r"""
Implements the simple mass and size estimation relationships for pump-fed
liquid rocket engines as described by B.T.C. Zandbergen (2015).

---
Approach Overview
---
This module provides Mass Estimating Relationships (MERs) based on statistical
[cite_start]regression analysis of historical rocket engine data[cite: 23, 46]. The models
are derived from a database of over 45 engines, covering a thrust range from
[cite_start]15 kN to 8 MN[cite: 24, 33].

The primary goal is to offer simple, "back of the envelope" approximations
for mass and envelope dimensions (length, diameter) suitable for the very
[cite_start]early stages of launch vehicle conceptual design[cite: 31].

[cite_start]The models are separated into two main propellant classes[cite: 145]:
1.  **Hydro-lox ("C")**: Cryogenic LOX/LH2 engines.
2.  **Storable/Kero-lox ("S")**: Storable propellants (e.g., NTO/A-50) and
    [cite_start]semi-cryogenic propellants (e.g., LOX/RP-1, LOX/LCH4)[cite: 145, 443].

---
Pros & Cons
---
**Pros:**
* **Simplicity:** The basic models require only thrust and propellant class
    [cite_start]as inputs[cite: 25].
* [cite_start]**Speed:** Ideal for rapid trade studies during conceptual design[cite: 31].
* **Known Uncertainty:** The paper quantifies the Relative Standard Error (RSE)
    for each model, providing a clear indication of its expected inaccuracy
    (typically between $\pm10\%$ and $\pm25\%$ for mass) [cite_start][cite: 27, 152, 464].
* **Flexibility:** Provides both simple (thrust-only) models and more complex
    multi-variable models that account for cycle type, chamber pressure,
    [cite_start]expansion ratio, and number of chambers[cite: 26, 152, 237].

**Cons:**
* **Statistical, Not Physical:** These are not physics-based models. They only
    reflect historical trends and cannot capture novel design features
    [cite_start]or technology improvements not present in the database[cite: 53].
* **Database Limitations:** The models are based on data from engines
    [cite_start]primarily from the pre-2000 era[cite: 137]. Their applicability to
    very modern designs (e.g., full-flow staged combustion, 3D-printed
    engines) is not guaranteed.
* **Extrapolation Risk:** The author explicitly warns that extrapolating
    the formulas far outside the original database's thrust range
    (15 kN - 8 MN) [cite_start]can introduce large errors[cite: 446].
* **Data Ambiguity:** The source data itself has limitations, such as
    unclear definitions of what "engine mass" includes (e.g., thrust
    [cite_start]frame, actuators)[cite: 84, 85, 456].

Reference:
Zandbergen, B.T.C. (2015). Simple mass and size estimation relationships
of pump fed rocket engines for launch vehicle conceptual design.
6th European Conference for Aeronautics and Space Sciences (EUCASS).

"""

from models.common_params import EngineParams
from models.base import BaseEngineModel, ModelResult
from typing import Dict, Literal, get_args
import inspect
import vehicle_definitions

# --- Type definitions for method selection ---
MassMethod = Literal[
    "default", "M-C1", "M-C2", "M-C3",  # Hydro-lox
    "M-S1", "M-S2"  # Storable/Kero-lox
]
LengthMethod = Literal[
    "default", "L-C1", "L-C2",  # Hydro-lox
    "L-S1", "L-S2"  # Storable/Kero-lox
]
DiameterMethod = Literal[
    "default", "D-C1", "D-C2",  # Hydro-lox
    "D-S1", "D-S2", "D-S3"  # Storable/Kero-lox
]

# Helper constant for unit conversion
PA_TO_BAR = 1.0 / 100_000.0


class ZandbergenEngineModel(BaseEngineModel):
    """
    Implements the Zandbergen (2015) engine mass model as a class.

    This class wraps the statistical regression formulas, allowing the
    selection of a specific formula (e.g., "M-C1", "M-S2") via the
    constructor.
    """

    def __init__(self, method: MassMethod = "default"):
        """
        Initializes the model with a specific calculation method.

        Args:
            method (MassMethod): The specific regression formula ID
                from Table 4 (e.g., "M-C1", "M-S1").
        """
        self.method = method
        self._method_used: str = ""  # Internal storage for the method name

    @property
    def model_name(self) -> str:
        """Returns the unique, human-readable name of the model."""
        # The _method_used is set during estimation, so we use it
        # to provide a more descriptive name.
        if self._method_used:
            return f"Zandbergen (2015) - {self._method_used}"
        return f"Zandbergen (2015) - {self.method}"

    def estimate_mass(self, params: EngineParams) -> ModelResult:
        """
        Estimates engine mass based on thrust and propellant type
        using the method specified in the constructor.

        Implements formulas from Table 4.
        Thrust (THRUST) must be in Newtons.

        Reference:
        Simple mass and size estimation relationships... (Zandbergen, 2015), Table 4.

        Args:
            params (EngineParams): The engine's parameter object.

        Returns:
            ModelResult: A TypedDict containing the total mass.

        Raises:
            ValueError: If the method is not valid for the propellant type
                        or if required parameters (e.g., Ae/At) are missing.
        """
        thrust_N = params.thrust_vac_N
        prop_type = params.propellant_type
        mass = 0.0
        method_to_run = self.method
        method_used = ""

        if prop_type == "LOX/LH2":
            # --- Hydro-lox Engines ---
            if method_to_run == "default":
                method_to_run = "M-C1"  # Default for LOX/LH2

            if method_to_run == "M-C1":
                # Formula M-C1: M(kg) = 0.00514 * THRUST(N)^0.92068
                mass = 0.00514 * (thrust_N ** 0.92068)
                method_used = "M-C1"

            elif method_to_run == "M-C2":
                # Formula M-C2: M(kg) = 1.866E-10*THRUST^2 + 0.00130*THRUST + 77.4
                mass = (1.866e-10 * (thrust_N ** 2)) + (0.00130 * thrust_N) + 77.4
                method_used = "M-C2"

            elif method_to_run == "M-C3":
                # Formula M-C3: M(kg) = 1.091E-4*THRUST^1.185 * 1.108^CYCLE + 104.0
                # CYCLE = 0 for staged combustion (SC), 1 for all other cycles
                cycle_param = 0.0 if params.cycle_type == "SC" else 1.0
                mass = (1.091e-4 * (thrust_N ** 1.185) * (1.108 ** cycle_param)) + 104.0
                method_used = "M-C3"

            else:
                raise ValueError(
                    f"Method '{method_to_run}' not valid for LOX/LH2. "
                    f"Use 'M-C1', 'M-C2', or 'M-C3'."
                )

        elif prop_type in ["LOX/RP1", "LOX/LCH4", "Storable"]:
            # --- Storable and Kero-lox Engines ---
            if method_to_run == "default":
                method_to_run = "M-S1"  # Default for this class

            if method_to_run == "M-S1":
                # Formula M-S1: M(kg) = 1.104E-3 * THRUST(N) + 27.702
                mass = (1.104e-3 * thrust_N) + 27.702
                method_used = "M-S1"

            elif method_to_run == "M-S2":
                # Formula M-S2: M=(1.079E-3*THRUST+53.0165)*NT^0.0369*EXP^-0.00184
                nt = params.num_chambers
                exp = params.expansion_ratio
                if exp <= 0:
                    raise ValueError("Expansion ratio must be > 0 for M-S2 calculation")
                if nt <= 0:
                    raise ValueError("Number of chambers must be > 0 for M-S2 calculation")

                mass = (1.079e-3 * thrust_N + 53.0165) * (nt ** 0.0369) * (exp ** -0.00184)
                method_used = "M-S2"

            else:
                raise ValueError(
                    f"Method '{method_to_run}' not valid for {prop_type}. "
                    f"Use 'M-S1' or 'M-S2'."
                )

        else:
            raise ValueError(
                f"Propellant type {prop_type} not supported by Zandbergen M-C/M-S models"
            )

        # Store the method that was actually used (resolving "default")
        self._method_used = method_used

        # Add applicability note based on prop type
        if prop_type == "LOX/LH2":
            app_note = "Simple regression (Thrust + prop class). Monolithic for LOX/LH2 engines"
        else:
            app_note = "Simple regression (Thrust + prop class). Monolithic for Storable, LOX/RP1, or LOX/LCH4 engines"

        # Zandbergen is monolithic (no components)
        return ModelResult(
            total_mass_kg=mass,
            components_kg={
                "Monolithic Mass": mass
            },
            notes={
                "method_used": method_used,
                "Applicability": app_note
            }
        )


# --- STANDALONE SIZE ESTIMATION FUNCTIONS ---


def estimate_engine_length(params: EngineParams, method: LengthMethod = "default") -> float:
    """
    Estimates overall engine length based on thrust and propellant type.

    Implements formulas from Table 5.
    Thrust (THRUST) must be in Newtons.
    Pressure (PRES) must be in bar.

    Available methods:
    - Hydro-lox ("LOX/LH2"):
        - "L-C1" (default): Thrust-based
        - "L-C2": Thrust, expansion_ratio (EXP), and chamber_pressure (PRES)
    - Storable/Kero-lox ("LOX/RP1", "LOX/LCH4", "Storable"):
        - "L-S1" (default): Thrust-based
        - "L-S2": Thrust, num_chambers (NT), and expansion_ratio (EXP)

    Reference:
    Simple mass and size estimation relationships... (Zandbergen, 2015), Table 5.

    Args:
        params (EngineParams): The engine's parameter object.
        method (LengthMethod): The specific regression formula ID from Table 5.
            Defaults to "default", which uses the simplest relation (L-C1 or L-S1).

    Returns:
        float: Estimated engine length in meters.
    """
    thrust_N = params.thrust_vac_N
    prop_type = params.propellant_type
    length = 0.0

    if prop_type == "LOX/LH2":
        # --- Hydro-lox Engines ---
        if method == "default" or method == "L-C1":
            # Formula L-C1: L(m) = 0.1667 * THRUST(N)^0.2238
            length = 0.1667 * (thrust_N ** 0.2238)

        elif method == "L-C2":
            # Formula L-C2: L=0.0615*THRUST^0.3091 * EXP^0.1037 * PRES^-0.1319
            exp = params.expansion_ratio
            pres_bar = params.chamber_pressure_Pa * PA_TO_BAR
            if exp <= 0:
                raise ValueError("Expansion ratio must be > 0 for L-C2 calculation")
            if pres_bar <= 0:
                raise ValueError("Chamber pressure must be > 0 for L-C2 calculation")

            length = 0.0615 * (thrust_N ** 0.3091) * (exp ** 0.1037) * (pres_bar ** -0.1319)

        else:
            raise ValueError(f"Method '{method}' not valid for LOX/LH2. Use {get_args(LengthMethod)}.")

    elif prop_type in ["LOX/RP1", "LOX/LCH4", "Storable"]:
        # --- Storable and Kero-lox Engines ---
        if method == "default" or method == "L-S1":
            # Formula L-S1: L(m) = 0.1362 * THRUST(N)^0.2279
            length = 0.1362 * (thrust_N ** 0.2279)

        elif method == "L-S2":
            # Formula L-S2: L=0.0879*THRUST^0.2570 * NT^-0.3815 * EXP^0.0485
            nt = params.num_chambers
            exp = params.expansion_ratio
            if nt <= 0:
                raise ValueError("Number of chambers must be > 0 for L-S2 calculation")
            if exp <= 0:
                raise ValueError("Expansion ratio must be > 0 for L-S2 calculation")

            length = 0.0879 * (thrust_N ** 0.2570) * (nt ** -0.3815) * (exp ** 0.0485)

        else:
            raise ValueError(f"Method '{method}' not valid for {prop_type}. Use {get_args(LengthMethod)}.")

    else:
        raise ValueError(f"Propellant type {prop_type} not supported by Zandbergen L-C/L-S models")

    return length


def estimate_engine_diameter(params: EngineParams, method: DiameterMethod = "default") -> float:
    """
    Estimates overall engine envelope diameter based on thrust and propellant type.

    Implements formulas from Table 5.
    Thrust (THRUST) must be in Newtons.
    Pressure (PRES) must be in bar.

    Available methods:
    - Hydro-lox ("LOX/LH2"):
        - "D-C1" (default): Thrust-based
        - "D-C2": Thrust, expansion_ratio (EXP), and chamber_pressure (PRES)
    - Storable/Kero-lox ("LOX/RP1", "LOX/LCH4", "Storable"):
        - "D-S1" (default): Thrust-based
        - "D-S2": Thrust and num_chambers (NT)
        - "D-S3": Thrust, NT, expansion_ratio (EXP), and chamber_pressure (PRES)

    Reference:
    Simple mass and size estimation relationships... (Zandbergen, 2015), Table 5.

    Args:
        params (EngineParams): The engine's parameter object.
        method (DiameterMethod): The specific regression formula ID from Table 5.
            Defaults to "default", which uses the simplest relation (D-C1 or D-S1).

    Returns:
        float: Estimated engine envelope diameter in meters.
    """
    thrust_N = params.thrust_vac_N
    prop_type = params.propellant_type
    diameter = 0.0

    if prop_type == "LOX/LH2":
        # --- Hydro-lox Engines ---
        if method == "default" or method == "D-C1":
            # Formula D-C1: D(m) = 0.1503 * THRUST(N)^0.1915
            diameter = 0.1503 * (thrust_N ** 0.1915)

        elif method == "D-C2":
            # Formula D-C2: D=3.85E-3*THRUST^0.4856 * EXP^0.4361 * PRES^-0.4760
            exp = params.expansion_ratio
            pres_bar = params.chamber_pressure_Pa * PA_TO_BAR
            if exp <= 0:
                raise ValueError("Expansion ratio must be > 0 for D-C2 calculation")
            if pres_bar <= 0:
                raise ValueError("Chamber pressure must be > 0 for D-C2 calculation")

            diameter = 3.85e-3 * (thrust_N ** 0.4856) * (exp ** 0.4361) * (pres_bar ** -0.4760)

        else:
            raise ValueError(f"Method '{method}' not valid for LOX/LH2. Use {get_args(DiameterMethod)}.")

    elif prop_type in ["LOX/RP1", "LOX/LCH4", "Storable"]:
        # --- Storable and Kero-lox Engines ---
        if method == "default" or method == "D-S1":
            # Formula D-S1: D(m) = 0.0455 * THRUST(N)^0.2745
            diameter = 0.0455 * (thrust_N ** 0.2745)

        elif method == "D-S2":
            # Formula D-S2: D=0.0415*THRUST^0.2788 * NT^0.1181
            nt = params.num_chambers
            if nt <= 0:
                raise ValueError("Number of chambers must be > 0 for D-S2 calculation")

            diameter = 0.0415 * (thrust_N ** 0.2788) * (nt ** 0.1181)

        elif method == "D-S3":
            # Formula D-S3: D=0.0413*THRUST^0.3110 * NT^0.2111 * EXP^0.1349 * PRES^-0.2030
            nt = params.num_chambers
            exp = params.expansion_ratio
            pres_bar = params.chamber_pressure_Pa * PA_TO_BAR
            if nt <= 0:
                raise ValueError("Number of chambers must be > 0 for D-S3 calculation")
            if exp <= 0:
                raise ValueError("Expansion ratio must be > 0 for D-S3 calculation")
            if pres_bar <= 0:
                raise ValueError("Chamber pressure must be > 0 for D-S3 calculation")

            diameter = 0.0413 * (thrust_N ** 0.3110) * (nt ** 0.2111) * (exp ** 0.1349) * (pres_bar ** -0.2030)

        else:
            raise ValueError(f"Method '{method}' not valid for {prop_type}. Use {get_args(DiameterMethod)}.")

    else:
        raise ValueError(f"Propellant type {prop_type} not supported by Zandbergen D-C/D-S models")

    return diameter


if __name__ == "__main__":
    """
    Standalone runner for the Zandbergen model.

    Allows testing a specific engine or all engines defined
    in vehicle_definitions.py.
    """

    # --- 1. CONFIGURE TEST RUN ---
    # Set this variable to a specific engine key (e.g., "ssme", "rd120") to test ONLY that one engine.
    # Set to None to test ALL available engines.
    SPECIFIC_ENGINE_KEY_TO_TEST = "default"
    # Example to test only SSME:
    # SPECIFIC_ENGINE_KEY_TO_TEST = "ssme"

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
        if not SPECIFIC_ENGINE_KEY_TO_TEST:
            print("[__main__] No engine getter functions found in vehicle_definitions.py.")
        else:
            print("[__main__] No analysis run due to configuration error.")

    for engine_key in keys_to_test:
        # We already validated the key in step 3, so we just get it
        print(f"\n[__main__] Loading parameters for '{engine_key}'...")
        engine_getter = available_engines[engine_key]

        try:
            engine_params = engine_getter()

            # Instantiate and run all applicable models for this engine

            print(f"\n--- Running Mass Models for {engine_key} ---")
            if engine_params.propellant_type == "LOX/LH2":
                methods: list[MassMethod] = ["M-C1", "M-C2", "M-C3"]
            else:
                methods: list[MassMethod] = ["M-S1", "M-S2"]

            for method in methods:
                model = ZandbergenEngineModel(method=method)
                model.run_single_engine_analysis(engine_params)

            # --- Run Size Estimations (which are still functions) ---
            print(f"\n--- Running Size Models for {engine_key} ---")

            # --- Length Estimates (m) ---
            print("\n--- Length Estimations (m) ---")
            if engine_params.propellant_type == "LOX/LH2":
                len_methods: list[LengthMethod] = ["L-C1", "L-C2"]
            else:
                len_methods: list[LengthMethod] = ["L-S1", "L-S2"]

            for method in len_methods:
                try:
                    res = estimate_engine_length(engine_params, method=method)
                    print(f"  {method:<6}: {res:12.2f} m")
                except Exception as e:
                    print(f"  {method:<6}: ERROR ({e})")

            # --- Diameter Estimates (m) ---
            print("\n--- Diameter Estimations (m) ---")
            if engine_params.propellant_type == "LOX/LH2":
                diam_methods: list[DiameterMethod] = ["D-C1", "D-C2"]
            else:
                diam_methods: list[DiameterMethod] = ["D-S1", "D-S2", "D-S3"]

            for method in diam_methods:
                try:
                    res = estimate_engine_diameter(engine_params, method=method)
                    print(f"  {method:<6}: {res:12.2f} m")
                except Exception as e:
                    print(f"  {method:<6}: ERROR ({e})")

            print("-" * 70)


        except Exception as e:
            print(f"\n[__main__] ERROR processing engine '{engine_key}': {e}")

    print("-" * 70)
    if keys_to_test:
        print(f"[__main__] Standalone analysis complete for {len(keys_to_test)} engine(s).")
    else:
        print("[__main__] No analysis was run.")