# main.py

from models.common_params import EngineParams
from models import akin_mers, mota_schlingloff
from models import zandbergen_engine, tizon_rema, kibbey_stage
from models.akin_mers import run_akin_ssto_example, print_ssto_results
from vehicle_definitions import (
    get_le5_engine, get_ssme_engine, get_rl10a_engine, get_rd120_engine, default_rocket_params
)


def run_engine_mass_comparison(engine_name: str, engine_params: EngineParams):
    """
    Runs all applicable ENGINE MASS ESTIMATION models
    for a given engine configuration and prints the results,
    including component breakdown where available.
    """
    print("=" * 70)
    print(f"   RUNNING COMPARISON FOR: {engine_name}")
    print(f"   (Params: {engine_params.thrust_vac_N / 1_000_000:.3f} MN, "
          f"P_c: {engine_params.chamber_pressure_Pa / 1e6:.2f} MPa, "
          f"Ae/At: {engine_params.expansion_ratio:.1f})")
    print("=" * 70)

    # Dictionary to store all results
    results = {}

    # --- 1. Akin (1991) Model ---
    try:
        # Calculate all propulsion-related MERs from Akin
        # that only depend on engine parameters.

        # Formula: M_Rocket_Engine(kg) = 7.81e-4*T(N) + 3.37e-5*T(N)*sqrt(Ae/At) + 59
        m_engine = akin_mers.estimate_engine_mass(
            engine_params.thrust_vac_N,
            engine_params.expansion_ratio
        )

        # Formula: M_Thrust_Structure(kg) = 2.55e-4*T(N)
        m_thrust_structure = akin_mers.estimate_thrust_structure_mass(
            engine_params.thrust_vac_N
        )

        # Formula: M_Gimbals(kg) = 237.8 * [T(N) / P_c(Pa)]^0.9375
        m_gimbals = akin_mers.estimate_gimbal_mass(
            engine_params.thrust_vac_N,
            engine_params.chamber_pressure_Pa
        )

        # Sum them up for a total "Akin Propulsion System" mass
        total_propulsion_mass = m_engine + m_thrust_structure + m_gimbals

        results["Akin (1991)"] = {
            "total_mass_kg": total_propulsion_mass,
            "components_kg": {
                "Engine (MER)": m_engine,
                "Thrust Structure (MER)": m_thrust_structure,
                "Gimbals (MER)": m_gimbals
            },
            "note": "Note: Propulsion system components (MERs). Ignores cycle."
        }
    except Exception as e:
        results["Akin (1991)"] = {"error": str(e)}

    # --- 2. Zandbergen (2015) Model ---
    try:
        # Run default (M-C1/M-S1)
        res_default = zandbergen_engine.estimate_engine_mass(engine_params)
        model_name_default = f"Zandbergen ({res_default['method_used']})"
        results[model_name_default] = {
            "total_mass_kg": res_default["total_mass_kg"],
            "note": "Note: Simple regression (Thrust + prop class). Monolithic."
        }

        # Run more complex models if applicable
        if engine_params.propellant_type == "LOX/LH2":
            # M-C3 (with cycle type)
            res_mc3 = zandbergen_engine.estimate_engine_mass(engine_params, method="M-C3")
            model_name_mc3 = f"Zandbergen ({res_mc3['method_used']})"
            results[model_name_mc3] = {
                "total_mass_kg": res_mc3["total_mass_kg"],
                "note": "Note: Hydro-lox regression (uses cycle type)."
            }

        elif engine_params.propellant_type in ["LOX/RP1", "LOX/LCH4", "Storable"]:
            # M-S2 (with NT, Expansion Ratio)
            res_ms2 = zandbergen_engine.estimate_engine_mass(engine_params, method="M-S2")
            model_name_ms2 = f"Zandbergen ({res_ms2['method_used']})"
            results[model_name_ms2] = {
                "total_mass_kg": res_ms2["total_mass_kg"],
                "note": "Note: Kero-lox regression (uses NT, Ae/At)."
            }

    except Exception as e:
        results["Zandbergen (2015)"] = {"error": str(e)}

    # --- 3. Tizón/RemA (2017) Model ---
    try:
        tizon_model = tizon_rema.TizonEngineModel(
            engine_params.cycle_type,
            method="historical"
        )
        tizon_result = tizon_model.estimate_total_mass(engine_params)
        ref_name = tizon_model.ref_engine['name']
        tizon_result["note"] = f"Note: *Relative* model. Reference: ({ref_name})."
        results["Tizón/RemA (2017)"] = tizon_result
    except ValueError as e:
        results["Tizón/RemA (2017)"] = {"error": f"Not applicable ({e})"}
    except Exception as e:
        results["Tizón/RemA (2017)"] = {"error": str(e)}

    # --- 4. Mota/Schlingloff (2005) Model ---
    try:
        has_boost_pumps = (engine_params.cycle_type == "SC")
        mota_result = mota_schlingloff.estimate_total_engine_mass(
            engine_params,
            has_boost_pumps=has_boost_pumps
        )
        mota_result["note"] = "Note: Analytical-statistical model, needs Pc, O/F, etc."
        results["Mota/Schlingloff (2005)"] = mota_result
    except Exception as e:
        results["Mota/Schlingloff (2005)"] = {"error": str(e)}

    # --- Print Results ---
    print("---  Engine Mass Estimation Results ---")

    # First, print the total mass summary table
    print("\n--- Total Mass Summary ---")
    for model_name, res in results.items():
        if "error" in res:
            mass_str = f"ERROR: {res['error']}"
            note = ""
        else:
            mass_str = f"{res.get('total_mass_kg', 0.0):,.1f} kg"
            note = res.get('note', '')

        print(f"{model_name:<28} | {mass_str:<20} | {note}")

    # Then, print component breakdowns where available
    print("\n--- Component Breakdown (where available) ---")
    for model_name, res in results.items():
        if "components_kg" in res:
            print(f"  --- {model_name} ---")
            total = res['total_mass_kg']
            for component, mass in res['components_kg'].items():
                percent = (mass / total * 100) if total > 0 else 0
                print(f"    {component:<30}: {mass:10,.1f} kg ({percent:5.1f}%)")
            print(f"    {'TOTAL':<30}: {total:10,.1f} kg (100.0%)")
            print()  # Empty line for separation

    print("\n")


def run_stage_model_example(engine_name: str, engine_params: EngineParams):
    """
    Runs an example of the STAGE MASS ESTIMATION model (Kibbey).
    """
    print("=" * 70)
    print(f"   RUNNING STAGE MODEL EXAMPLE (Kibbey, 2015)")
    print(f"   Using engine: {engine_name}")
    print(f"   Reference Stage: Atlas V (LOX/RP1)")
    print("=" * 70)

    # Kibbey's model requires mission params R (mass ratio) and Psi (T/W)
    # We use the Atlas V reference values for demonstration

    # Take Atlas V target params as an example
    example_mission_R = 6.62
    example_mission_Psi = 1.28

    try:
        kibbey_model = kibbey_stage.KibbeyStageModel()

        stage_fractions = kibbey_model.estimate_propellant_mass_fraction(
            engine_params,
            mission_mass_ratio=example_mission_R,
            mission_T_W=example_mission_Psi
        )

        lambda_frac = stage_fractions['propellant_mass_fraction']
        inert_frac = 1.0 - lambda_frac

        print(f"---   Kibbey Stage Model Results ---")

        # Add a note about relevance
        if engine_params.propellant_type == "LOX/RP1":
            print("(!) NOTE: This is the most relevant test, as both the engine (LOX/RP1)")
            print("    and the reference (Atlas V, LOX/RP1) use the same prop class.")
        else:
            print(f"(!) NOTE: This is an extrapolation. The model is scaling the LOX/RP1")
            print(f"    reference stage to work with {engine_params.propellant_type} engine params.")

        print(f"\nFor a stage with R={example_mission_R} and T/W={example_mission_Psi}:")
        print(f"  Estimated Propellant Fraction (lambda): {lambda_frac:.4f}")
        print(f"  Estimated Inert Fraction (1-lambda):    {inert_frac:.4f}")
        print(f"  (f_i_total, inert-to-propellant-mass):  {stage_fractions['total_inert_fraction']:.4f}")

        print("\n  Inert Fraction Breakdown (f_i):")
        total_f_i = stage_fractions['total_inert_fraction']
        for comp, val in stage_fractions['component_fractions'].items():
            percent = (val / total_f_i * 100) if total_f_i > 0 else 0
            print(f"    {comp:<10}: {val:8.5f} ({percent:5.1f}%)")
        print(f"    {'TOTAL':<10}: {total_f_i:8.5f} (100.0%)")


    except Exception as e:
        print(f"ERROR running Kibbey model: {e}")
    print("\n")


def run_akin_ssto_analysis_with_engine(engine_name: str, engine_params: EngineParams):
    """
    Runs the Akin SSTO 1st Pass analysis, but OVERRIDES
    the default engine parameters with those from a specific engine.
    """
    try:
        # 1. Get the default *mission* config
        # We only need the 'stage' object, the 'engine' will be replaced
        _default_engine, stage_params = default_rocket_params()

        # 2. OVERRIDE engine-specific params
        print(f"... Using {engine_name} engine data...")
        # The new engine is passed directly to the function.
        # We also need to update the stage's internal engine reference.
        stage_params.engine = engine_params

        # 3. Run the analysis
        results = run_akin_ssto_example(engine_params, stage_params)

        # 4. Print the results
        print_ssto_results(results, show_pdf_ref=False)

    except Exception as e:
        print(f"ERROR running Akin SSTO example: {e}")
    print("\n")


if __name__ == "__main__":

    # --- 1. Define test case engines ---
    # Store engines in a dictionary for easy, repeatable access
    # Keys are simple identifiers, values are (Name String, Engine Object)
    test_engines = {
        "le5": ("LE-5 (LOX/LH2 GG)", get_le5_engine()),
        "ssme": ("SSME (LOX/LH2 SC)", get_ssme_engine()),
        "rl10": ("RL10A (LOX/LH2 EX)", get_rl10a_engine()),
        "rd120": ("RD-120 (LOX/RP1 SC)", get_rd120_engine()),
    }

    # --- Define a separate *scenario* for a custom run ---
    # default_rocket_params() returns (EngineParams, StageParams)
    custom_engine, custom_stage = default_rocket_params()
    custom_scenario_name = "Custom Run"

    # --- 2. Run Engine Mass Comparison ---
    print("||||" + "=" * 70 + "||||")
    print("||||   PART 1: ENGINE MASS MODEL COMPARISON")
    print("||||" + "=" * 70 + "||||")

    # Loop through all defined engines
    # This will now work correctly
    for engine_name, engine_params in test_engines.values():
        run_engine_mass_comparison(engine_name, engine_params)

    # --- 3. Run Stage Model Example ---
    print("||||" + "=" * 70 + "||||")
    print("||||   PART 2: STAGE MASS MODEL EXAMPLE (Kibbey)")
    print("||||" + "=" * 70 + "||||")

    # Kibbey's model is calibrated on Atlas V (LOX/RP1).
    # We run the RD-120 (also LOX/RP1) as the most relevant test.
    engine_name_rd120, params_rd120 = test_engines["rd120"]
    run_stage_model_example(engine_name_rd120, params_rd120)

    # Also run for SSME to see the extrapolation
    engine_name_ssme, params_ssme = test_engines["ssme"]
    run_stage_model_example(engine_name_ssme, params_ssme)

    # --- 4. Run Akin SSTO Full Example ---
    print("||||" + "=" * 70 + "||||")
    print("||||   PART 3: AKIN MERS SSTO EXAMPLE")
    print("||||" + "=" * 70 + "||||")

    # Scenario 1: Run default mission but override with SSME engine
    run_akin_ssto_analysis_with_engine(engine_name_ssme, params_ssme)

    # Scenario 2: Run the fully custom scenario
    print("\n" + "=" * 70)
    print(f"RUNNING AKIN SSTO CUSTOM SCENARIO ({custom_scenario_name})")
    print("=" * 70)

    # Correctly call the core analysis function
    # with the custom engine and custom stage
    try:
        custom_results = run_akin_ssto_example(custom_engine, custom_stage)
        # Print the results (using the helper from akin_mers)
        print_ssto_results(custom_results, show_pdf_ref=False)
    except Exception as e:
        print(f"ERROR running Akin Custom SSTO example: {e}")
    print("\n")