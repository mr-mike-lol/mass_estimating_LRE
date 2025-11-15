# main.py

"""
Main execution file for the rocket engine and stage mass model comparison.

This file orchestrates the analysis by:
1. Defining a set of test engines (from vehicle_definitions).
2. Defining a list of engine and stage model classes to run.
3. Looping through each test engine and running it against each
   applicable model.
4. Running standalone, full-stage analyses (like the Akin SSTO example).
"""

from typing import List, Tuple, Dict, Any
from models.common_params import EngineParams, StageParams, CycleType

from models.akin_mers import (
    AkinPropulsionModel, run_akin_ssto_example, print_ssto_results
)
from models.mota_schlingloff import MotaSchlingloffModel
from models.zandbergen_engine import ZandbergenEngineModel, MassMethod
from models.tizon_rema import TizonEngineModel
from models.kibbey_stage import KibbeyStageModel

# Import the base classes
from models.base import BaseEngineModel, BaseStageModel

# Import the engine/stage definitions
from vehicle_definitions import (
    get_le5_engine, get_ssme_engine, get_rl10a_engine,
    get_rd120_engine, default_rocket_params
)


def get_applicable_engine_models(params: EngineParams) -> List[BaseEngineModel]:
    """
    Instantiates and returns a list of all engine models
    applicable to the given engine's parameters.

    This acts as a "factory" that configures the models (e.g., setting
    boost pump flags or cycle types) based on the engine being tested.

    Args:
        params (EngineParams): The engine to be tested.

    Returns:
        List[BaseEngineModel]: A list of configured model instances.
    """
    models: List[BaseEngineModel] = []

    # 1. Akin (always applicable)
    models.append(AkinPropulsionModel())

    # 2. Zandbergen (add all sub-models for the correct prop type)
    if params.propellant_type == "LOX/LH2":
        models.append(ZandbergenEngineModel(method="M-C1"))
        models.append(ZandbergenEngineModel(method="M-C2"))
        models.append(ZandbergenEngineModel(method="M-C3"))
    else:
        # For LOX/RP1, LOX/LCH4, Storable
        models.append(ZandbergenEngineModel(method="M-S1"))
        models.append(ZandbergenEngineModel(method="M-S2"))

    # 3. Tizon (only for the matching cycle type)
    # This try/except block automatically handles applicability.
    try:
        # TizonEngineModel constructor raises ValueError if cycle is unsupported
        models.append(
            TizonEngineModel(cycle_type=params.cycle_type, method="historical")
        )
        models.append(
            TizonEngineModel(cycle_type=params.cycle_type, method="design")
        )
    except ValueError:
        # This engine's cycle (e.g., from RD-120) is not supported
        # by Tizon (which only has 'SC', 'GG', 'EX' for LOX/LH2), so we skip it.
        pass

    # 4. Mota/Schlingloff (always applicable, configures boost pumps)
    # We configure the model based on the engine's cycle type.
    has_boost_pumps = (params.cycle_type == "SC")
    models.append(MotaSchlingloffModel(has_boost_pumps=has_boost_pumps))

    return models


def get_applicable_stage_models() -> List[BaseStageModel]:
    """
    Instantiates and returns a list of all available stage models.

    Returns:
        List[BaseStageModel]: A list of configured model instances.
    """
    models: List[BaseStageModel] = []

    # 1. Kibbey (only one we have refactored so far)
    models.append(KibbeyStageModel())

    # In the future, we could add:
    # models.append(AkinSSTOModel())

    return models


if __name__ == "__main__":

    # --- 1. Define test case engines ---
    # Store engines in a dictionary for easy, repeatable access
    # Keys are simple identifiers, values are (Name String, Engine Object)
    test_engines: Dict[str, Tuple[str, EngineParams]] = {
        "le5": ("LE-5 (LOX/LH2 GG)", get_le5_engine()),
        "ssme": ("SSME (LOX/LH2 SC)", get_ssme_engine()),
        "rl10": ("RL10A (LOX/LH2 EX)", get_rl10a_engine()),
        "rd120": ("RD-120 (LOX/RP1 SC)", get_rd120_engine()),
    }

    # --- 2. Run Engine Mass Comparison ---
    print("||||" + "=" * 70 + "||||")
    print("||||   PART 1: ENGINE MASS MODEL COMPARISON")
    print("||||" + "=" * 70 + "||||")

    for engine_key, (engine_name, engine_params) in test_engines.items():
        print("\n" + "#" * 70)
        print(f"###   RUNNING ALL MODELS FOR: {engine_name}")
        print("#" * 70)

        # Get the list of models configured for this specific engine
        engine_models_to_run = get_applicable_engine_models(engine_params)

        if not engine_models_to_run:
            print("No applicable models found for this engine.")
            continue

        for model in engine_models_to_run:
            # Each model now runs its *own* standardized analysis
            # method, which we inherited from BaseEngineModel.
            model.run_single_engine_analysis(engine_params)

    # --- 3. Run Stage Model Example ---
    print("||||" + "=" * 70 + "||||")
    print("||||   PART 2: STAGE MASS MODEL EXAMPLE (Kibbey)")
    print("||||" + "=" * 70 + "||||")

    stage_models_to_run = get_applicable_stage_models()

    # We need a set of stage/mission parameters to run against.
    # We'll use the one from our default custom rocket.
    _custom_engine, custom_stage_params = default_rocket_params()

    # Test 1: RD-120 (LOX/RP1, most relevant for Kibbey's Atlas V ref)
    _rd120_name, rd120_params = test_engines["rd120"]
    print(f"\n--- Running Stage Models for RD-120 (LOX/RP1) ---")
    for model in stage_models_to_run:
        model.run_single_stage_analysis(rd120_params, custom_stage_params)

    # Test 2: SSME (LOX/LH2, extrapolation test)
    _ssme_name, ssme_params = test_engines["ssme"]
    print(f"\n--- Running Stage Models for SSME (LOX/LH2) ---")
    for model in stage_models_to_run:
        model.run_single_stage_analysis(ssme_params, custom_stage_params)

    # --- 4. Run Akin SSTO Full Example ---
    print("||||" + "=" * 70 + "||||")
    print("||||   PART 3: AKIN MERS SSTO EXAMPLE")
    print("||||" + "=" * 70 + "||||")

    # This analysis is more complex than a single model call,
    # so we run its specific function directly.
    # We have not refactored this part into a BaseStageModel yet.

    # Scenario 1: Run default mission but override with SSME engine
    print("\n--- Running Akin SSTO with SSME Engine ---")
    _ssme_name_ssto, ssme_params_ssto = test_engines["ssme"]

    # Get the default stage params
    _default_engine, stage_params_ssto = default_rocket_params()
    # Override the engine in the stage object
    stage_params_ssto.engine = ssme_params_ssto

    try:
        results_ssme = run_akin_ssto_example(
            ssme_params_ssto, stage_params_ssto
        )
        print_ssto_results(results_ssme, show_pdf_ref=False)
    except Exception as e:
        print(f"ERROR running Akin SSTO (SSME): {e}\n")

    # Scenario 2: Run the fully custom scenario
    print("\n--- Running Akin SSTO with Custom Params ---")
    custom_engine, custom_stage = default_rocket_params()
    try:
        custom_results = run_akin_ssto_example(custom_engine, custom_stage)
        print_ssto_results(custom_results, show_pdf_ref=False)
    except Exception as e:
        print(f"ERROR running Akin SSTO (Custom): {e}\n")