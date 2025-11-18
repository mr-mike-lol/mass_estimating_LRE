# models/base.py

"""
Defines the abstract base class (interface) for all engine and stage mass models.

This module establishes the "Strategy" design pattern, allowing different
models (strategies) to be used interchangeably by the main application.
Each model must adhere to the contract defined by its Base Class.
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, TypedDict, Tuple
from models.common_params import EngineParams, StageParams


# --- Engine Model Interface ---

class ModelResult(TypedDict):
    """
    A standardized dictionary structure for all model estimation results.

    Attributes:
        total_mass_kg: The final, total estimated mass.
        components_kg: A dictionary of component names and their masses.
                       For monolithic models, this may contain a single
                       entry or be empty.
        notes: A dictionary for any extra metadata, such as the
               specific sub-model used (e.g., "M-C1") or reference engine.
    """
    total_mass_kg: float
    components_kg: Dict[str, float]
    notes: Dict[str, Any]


class BaseEngineModel(ABC):
    """
    Abstract Base Class for an ENGINE mass estimation model.

    This class defines the common interface that all concrete engine
    models must implement.
    """

    @property
    @abstractmethod
    def model_name(self) -> str:
        """
        Returns the unique, human-readable name of the model.

        Example:
            "Tizon/RemA (2017)"
        """
        pass

    @abstractmethod
    def estimate_mass(self, params: EngineParams) -> ModelResult:
        """
        Estimates the engine mass based on input parameters.

        Args:
            params (EngineParams): The standardized dataclass containing
                                   all necessary inputs for the model.

        Returns:
            ModelResult: A TypedDict containing the total mass, component
                         breakdown, and any relevant notes.

        Raises:
            ValueError: If the provided EngineParams are invalid or
                        unsupported by the specific model (e.g., wrong
                        cycle type, non-physical values).
        """
        pass

    # --- Concrete Helper Methods ---

    def run_single_engine_analysis(self, params: EngineParams):
        """
        Runs a standardized analysis and prints a formatted report.

        This method is implemented on the base class and uses the
        abstract 'model_name' and 'estimate_mass' methods.
        It provides a consistent test-bench runner for all subclasses.

        Args:
            params (EngineParams): The engine parameters to analyze.
        """
        print("=" * 70)
        print(f"Running Single Engine Analysis for: {self.model_name}")
        print(f"  Thrust: {params.thrust_vac_N / 1_000_000:.3f} MN")
        print(f"  Prop:   {params.propellant_type}")
        print(f"  Cycle:  {params.cycle_type}")
        print(f"  Pc:     {params.chamber_pressure_Pa / 1e6:.2f} MPa")
        print(f"  Ae/At:  {params.expansion_ratio:.1f}")
        print("=" * 70)

        try:
            # Call the abstract method, which is implemented by the subclass
            result = self.estimate_mass(params)

            total = result['total_mass_kg']

            print("  --- Results ---")

            # Handle both component and monolithic results
            components = result.get('components_kg', {})
            if components and total > 0:
                for component, mass in components.items():
                    percent = (mass / total * 100)
                    print(f"    {component:<30}: {mass:10,.1f} kg ({percent:5.1f}%)")
            elif total > 0:
                # Monolithic model with no component breakdown
                print(f"    {'Total (Monolithic)':<30}: {total:10,.1f} kg (100.0%)")
            else:
                print("    No component breakdown provided.")

            print(f"    {'-' * 30}: {'-' * 10} -- {'-' * 7}")
            print(f"    {'TOTAL MASS':<30}: {total:10,.1f} kg (100.0%)")

            # Print notes if they exist
            notes = result.get('notes', {})
            if notes:
                print("\n  --- Model Notes ---")
                for key, value in notes.items():
                    print(f"    {key}: {value}")

        except ValueError as e:
            # Handle known modeling errors (e.g., wrong cycle)
            print(f"  ERROR: {e}")
        except Exception as e:
            # Handle unexpected programming errors
            print(f"  An unexpected error occurred: {e}")

        print("-" * 70)
        print("\n")


# --- Stage Model Interface ---

class StageModelResult(TypedDict):
    """
    A standardized dictionary structure for FRACTION-based results (e.g., Kibbey).

    Attributes:
        propellant_mass_fraction: (lambda) Stage propellant / Stage total mass
        total_inert_fraction: (f_i_total) Stage inert / Stage propellant mass
        component_fractions: A nested dict with the 3 inert components
        notes: A dictionary for any extra metadata.
    """
    propellant_mass_fraction: float
    total_inert_fraction: float
    component_fractions: Dict[str, float]
    notes: Dict[str, Any]


class StageMassBudget(TypedDict):
    """
    A standardized dictionary for a full, ABSOLUTE MASS budget (in kg).

    Attributes:
        gross_mass_kg: Total vehicle mass (M_o).
        propellant_mass_kg: Total propellant mass (M_p).
        payload_mass_kg: Payload mass (M_pay).
        total_inert_mass_kg: Sum of all dry components (M_inert).
        components_kg: Dictionary of specific component masses (Tanks, Engines, etc.).
        notes: Metadata (including design margin, iteration data).
    """
    gross_mass_kg: float
    propellant_mass_kg: float
    payload_mass_kg: float
    total_inert_mass_kg: float
    components_kg: Dict[str, float]
    notes: Dict[str, Any]


class BaseStageModel(ABC):
    """
    Abstract Base Class for a STAGE mass estimation model.

    This class implements a "Template Method" pattern:
    1. `calculate_full_stage_mass_budget` (Concrete): The skeleton of the calculation.
    2. `_calculate_...` (Abstract): The specific implementation steps (Akin, Tizon, etc.).
    3. `calculate_propellant_volumes` (Concrete): Shared physics helper.
    """

    @property
    @abstractmethod
    def model_name(self) -> str:
        """Returns the unique, human-readable name of the model."""
        pass

    # --- 1. CONCRETE SHARED METHODS (Physics/Geometry) ---

    def calculate_propellant_volumes(self,
                                     params: EngineParams,
                                     propellant_mass_kg: float,
                                     ullage_factor: float = 1.0
                                     ) -> Tuple[float, float]:
        """
        Calculates Oxidizer and Fuel volumes based on mixture ratio and densities.

        This is a shared physics helper used by both Akin and Tizon models.

        Formulas:
            M_ox = M_prop * (MR / (1+MR))
            M_fuel = M_prop * (1 / (1+MR))
            V = M / rho

        Args:
            params (EngineParams): Engine parameters (contains densities and O/F ratio).
            propellant_mass_kg (float): Total propellant mass.
            ullage_factor (float, optional): Multiplier for tank volume (e.g., 1.03).
                                             Defaults to 1.0.

        Returns:
            Tuple[float, float]: (vol_oxidizer_m3, vol_fuel_m3) including ullage.
        """
        if propellant_mass_kg <= 0:
            return 0.0, 0.0

        mr = params.mixture_ratio
        # Calculate component masses
        m_ox = propellant_mass_kg * (mr / (1.0 + mr))
        m_fuel = propellant_mass_kg * (1.0 / (1.0 + mr))

        # Calculate net volumes
        if params.oxidizer_density <= 0 or params.fuel_density <= 0:
            return 0.0, 0.0

        v_ox_net = m_ox / params.oxidizer_density
        v_fuel_net = m_fuel / params.fuel_density

        # Apply ullage factor
        return (v_ox_net * ullage_factor, v_fuel_net * ullage_factor)

    # --- 2. ABSTRACT METHODS (The Implementation Details) ---

    @abstractmethod
    def _calculate_initial_sizing(self,
                                  params: EngineParams,
                                  stage: StageParams
                                  ) -> Tuple[float, float, float]:
        """
        Performs the 1st pass sizing iteration (Rocket Equation).

        Akin iterates on Mass Fraction (lambda).
        Tizon iterates on Inert Fraction (delta) or uses fixed propellant load.

        Returns:
            Tuple[float, float, float]: (M_gross, M_propellant, M_inert_initial_guess)
        """
        pass

    @abstractmethod
    def _calculate_tank_mass(self,
                             params: EngineParams,
                             stage: StageParams
                             ) -> float:
        """
        Estimates the mass of the propellant tanks.

        Implementations:
        - Akin: Uses volume-based regression (M = a * V).
        - Tizon: Uses pressure vessel physics (Stress, Thickness, Density).
        """
        pass

    @abstractmethod
    def _calculate_propulsion_system_mass(self,
                                          params: EngineParams,
                                          stage: StageParams
                                          ) -> float:
        """
        Estimates the mass of the full propulsion system.

        Implementations:
        - Akin: Sum of Engine MER + Thrust Structure MER + Gimbal MER.
        - Tizon: Sum of Engine Model (scaled) * Num_Engines.
        """
        pass

    @abstractmethod
    def _calculate_structural_mass(self,
                                   params: EngineParams,
                                   stage: StageParams
                                   ) -> float:
        """
        Estimates the mass of structural elements (Fairings, Intertanks, Adapters).
        """
        pass

    @abstractmethod
    def _calculate_avionics_and_power_mass(self,
                                           params: EngineParams,
                                           stage: StageParams
                                           ) -> float:
        """
        Estimates the mass of avionics, wiring, and power systems.
        Can return 0.0 if the specific model does not account for this.
        """
        pass

    @abstractmethod
    def _calculate_insulation_mass(self,
                                   params: EngineParams,
                                   stage: StageParams
                                   ) -> float:
        """
        Estimates the mass of cryogenic insulation.
        Can return 0.0 if the model does not account for this.
        """
        pass

    # --- 3. TEMPLATE METHOD (The Skeleton) ---

    def calculate_full_stage_mass_budget(self,
                                         params: EngineParams,
                                         stage: StageParams
                                         ) -> StageMassBudget:
        """
        Executes the full stage mass estimation process.

        This "Template Method" defines the algorithm's skeleton:
        1. Sizing (Gross/Propellant Mass).
        2. Component Estimation (calling abstract methods).
        3. Summation and Margin Calculation.

        Args:
            params (EngineParams): Engine configuration.
            stage (StageParams): Stage configuration.

        Returns:
            StageMassBudget: The complete mass breakdown in kg.
        """
        components: Dict[str, float] = {}
        notes: Dict[str, Any] = {"model_name": self.model_name}

        try:
            # --- Step 1: Initial Sizing ---
            M_o, M_p, M_i_guess = self._calculate_initial_sizing(params, stage)

            # Update the mutable StageParams object for subsequent steps
            stage.vehicle_gross_mass_kg = M_o
            stage.propellant_mass_kg = M_p
            stage.stage_inert_mass_kg = M_i_guess  # Temporary storage of guess

            notes['M_o_guess_kg'] = M_o
            notes['M_p_guess_kg'] = M_p

            # --- Step 2: Component Estimation ---
            # Concrete child classes implement these specific logic steps
            components['Tanks'] = self._calculate_tank_mass(params, stage)
            components['Propulsion'] = self._calculate_propulsion_system_mass(params, stage)
            components['Structure'] = self._calculate_structural_mass(params, stage)
            components['Avionics/Power'] = self._calculate_avionics_and_power_mass(params, stage)
            components['Insulation'] = self._calculate_insulation_mass(params, stage)

            # --- Step 3: Summation & Summary ---
            total_inert_calculated = sum(components.values())

            # Calculate margin: (Guess - Calculated) / Guess
            if M_i_guess > 0:
                margin = (M_i_guess - total_inert_calculated) / M_i_guess
            else:
                margin = 0.0

            notes['design_margin_percent'] = margin * 100.0

            return StageMassBudget(
                gross_mass_kg=M_o,
                propellant_mass_kg=M_p,
                payload_mass_kg=stage.payload_mass_kg,
                total_inert_mass_kg=total_inert_calculated,
                components_kg=components,
                notes=notes
            )

        except Exception as e:
            # Graceful error handling
            return StageMassBudget(
                gross_mass_kg=0.0,
                propellant_mass_kg=0.0,
                payload_mass_kg=stage.payload_mass_kg,
                total_inert_mass_kg=0.0,
                components_kg={},
                notes={"ERROR": f"Calculation failed: {str(e)}"}
            )

    # --- 4. LEGACY SUPPORT ---

    @abstractmethod
    def estimate_stage_fractions(self,
                                 params: EngineParams,
                                 stage_params: StageParams) -> StageModelResult:
        """
        Legacy interface for fraction-based models (e.g., Kibbey).
        Must be implemented, but can simply wrap calculate_full_stage_mass_budget
        if converting absolute mass to fractions.
        """
        pass

    def run_single_stage_analysis(self,
                                  params: EngineParams,
                                  stage_params: StageParams):
        """
        Runs a simple reporting routine.
        Defaults to calling the legacy estimate_stage_fractions for backward compatibility.
        """
        print("=" * 70)
        print(f"Running Single Stage Analysis for: {self.model_name}")
        print(f"  Engine Prop: {params.propellant_type}")
        print(f"  Target dV:   {stage_params.delta_v_ms} m/s")
        print("=" * 70)

        try:
            # Default behavior: Run the fraction model
            result = self.estimate_stage_fractions(params, stage_params)

            lambda_frac = result['propellant_mass_fraction']
            total_f_i = result['total_inert_fraction']

            print("  --- Results (Fractions) ---")
            print(f"    Propellant Fraction (lambda): {lambda_frac:.4f}")
            print(f"    Inert Fraction (f_i):         {total_f_i:.4f}")

            comps = result.get('component_fractions', {})
            if comps:
                print("\n  --- Component Breakdown (f_i) ---")
                for k, v in comps.items():
                    print(f"    {k:<20}: {v:.5f}")

        except Exception as e:
            print(f"  ERROR: {e}")

        print("-" * 70)
        print("\n")