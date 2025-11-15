# models/base.py

"""
Defines the abstract base class (interface) for all engine and stage mass models.

This module establishes the "Strategy" design pattern, allowing different
models (strategies) to be used interchangeably by the main application.
Each model must adhere to the contract defined by its Base Class.
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, TypedDict
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
    A standardized dictionary structure for all STAGE model estimation results.

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


class BaseStageModel(ABC):
    """
    Abstract Base Class for a STAGE mass estimation model.

    This class defines the common interface that all concrete stage
    models must implement.
    """

    @property
    @abstractmethod
    def model_name(self) -> str:
        """
        Returns the unique, human-readable name of the model.

        Example:
            "Kibbey (2015) Stage Model"
        """
        pass

    @abstractmethod
    def estimate_stage_fractions(self,
                                 params: EngineParams,
                                 stage_params: StageParams) -> StageModelResult:
        """
        Estimates the stage mass fractions based on engine and stage params.

        Args:
            params (EngineParams): The standardized engine dataclass.
            stage_params (StageParams): The standardized stage dataclass.

        Returns:
            StageModelResult: A TypedDict containing the mass fractions
                              and any relevant notes.

        Raises:
            ValueError: If the provided params are invalid or
                        unsupported by the specific model.
        """
        pass

    def run_single_stage_analysis(self,
                                  params: EngineParams,
                                  stage_params: StageParams):
        """
        Runs a standardized analysis and prints a formatted report
        for a stage model.

        Args:
            params (EngineParams): The engine parameters to analyze.
            stage_params (StageParams): The stage parameters to analyze.
        """
        print("=" * 70)
        print(f"Running Single Stage Analysis for: {self.model_name}")
        print(f"  Engine Prop: {params.propellant_type}")
        print(f"  Engine Cycle: {params.cycle_type}")
        print(f"  Target dV:   {stage_params.delta_v_ms} m/s")
        print(f"  Target T/W:  {stage_params.initial_twr}")
        print("=" * 70)

        try:
            # Call the abstract method, which is implemented by the subclass
            result = self.estimate_stage_fractions(params, stage_params)

            lambda_frac = result['propellant_mass_fraction']
            inert_frac = 1.0 - lambda_frac
            total_f_i = result['total_inert_fraction']

            print("  --- Results ---")
            print(f"    Estimated Propellant Fraction (lambda): {lambda_frac:.4f}")
            print(f"    Estimated Inert Fraction (1-lambda):    {inert_frac:.4f}")
            print(f"    (f_i_total, inert-to-propellant-mass):  {total_f_i:.4f}")

            print("\n  Inert Fraction Breakdown (f_i):")
            components = result.get('component_fractions', {})
            if components and total_f_i > 0:
                for comp, val in components.items():
                    percent = (val / total_f_i * 100)
                    print(f"    {comp:<10}: {val:8.5f} ({percent:5.1f}%)")
                print(f"    {'TOTAL':<10}: {total_f_i:8.5f} (100.0%)")
            else:
                print("    No component breakdown provided.")

            # Print notes if they exist
            notes = result.get('notes', {})
            if notes:
                print("\n  --- Model Notes ---")
                for key, value in notes.items():
                    print(f"    {key}: {value}")

        except ValueError as e:
            print(f"  ERROR: {e}")
        except Exception as e:
            print(f"  An unexpected error occurred: {e}")

        print("-" * 70)
        print("\n")