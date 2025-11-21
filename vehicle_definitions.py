# vehicle_definitions.py

from models.common_params import EngineParams, StageParams, TwoStageConfig
from typing import Dict, Any, Tuple
from models.common_params import (
    DENSITY_RP1, DENSITY_LH2, DENSITY_LOX, DENSITY_LCH4, G0
)

def get_le5_engine() -> EngineParams:
    """
    Defines a test case for the LE-5 engine.

    This is a LOX/LH2 Gas Generator (GG) cycle engine.

    Reference:
    Data primarily from Tizón & Román, 2017, Table 7.
    """
    return EngineParams(
        thrust_vac_N=103_000.0,
        isp_vac_s=450.0,  # Note: Isp not in Table 7, using typical value from Akin
        chamber_pressure_Pa=3.65e6,
        propellant_type="LOX/LH2",
        cycle_type="GG",
        mixture_ratio=5.5,
        expansion_ratio=140.0,
        fuel_density=DENSITY_LH2,
        oxidizer_density=DENSITY_LOX,
        num_chambers=1,
        safety_factor=1.1  # Tizon 2017, p. 15 (non-crewed)
    )


def get_ssme_engine() -> EngineParams:
    """
    Defines a test case for the SSME (RS-25) engine.

    This is a LOX/LH2 Staged Combustion (SC) cycle engine.

    Reference:
    Data primarily from Tizón & Román, 2017, Table 7.
    """
    return EngineParams(
        thrust_vac_N=2_280_000.0,
        isp_vac_s=452.0,  # Known public value
        chamber_pressure_Pa=2.04e7,
        propellant_type="LOX/LH2",
        cycle_type="SC",
        mixture_ratio=6.0,
        expansion_ratio=77.5,
        fuel_density=DENSITY_LH2,
        oxidizer_density=DENSITY_LOX,
        num_chambers=1,
        safety_factor=1.2  # Tizon 2017, p. 15 (crewed mission)
    )


def get_rl10a_engine() -> EngineParams:
    """
    Defines a test case for the RL10A-3-A engine.

    This is a LOX/LH2 Expander (EX) cycle engine.

    Reference:
    Data primarily from Tizón & Román, 2017, Table 7.
    """
    return EngineParams(
        thrust_vac_N=73_400.0,
        isp_vac_s=444.0,  # Known public value
        chamber_pressure_Pa=3.20e6,
        propellant_type="LOX/LH2",
        cycle_type="EX",
        mixture_ratio=5.0,
        expansion_ratio=61.1,
        fuel_density=DENSITY_LH2,
        oxidizer_density=DENSITY_LOX,
        num_chambers=1,
        safety_factor=1.1  # Tizon 2017, p. 15 (non-crewed)
    )


def get_rd120_engine() -> EngineParams:
    """
    Defines a test case for the RD-120 engine.

    This is a LOX/RP-1 Staged Combustion (SC) cycle engine.

    Reference:
    Data from Zandbergen, 2015, Table 3 [cite: 96]
    """
    return EngineParams(
        thrust_vac_N=833_565.0,
        isp_vac_s=350.0,  # Typical for LOX/RP-1 SC
        chamber_pressure_Pa=16.0e6,  # Known public value
        propellant_type="LOX/RP1",
        cycle_type="SC",
        mixture_ratio=2.6,  # Known public value
        expansion_ratio=106.0,  # Known public value
        fuel_density=DENSITY_RP1,
        oxidizer_density=DENSITY_LOX,
        num_chambers=1,
        safety_factor=1.1  # Assuming non-crewed
    )


def get_akin_ssto_default_params_1st() -> Tuple[EngineParams, StageParams]:
    """
    Returns the default configuration for the SSTO 1st Pass example
    from the "Mass Estimating Relations" PDF (Pages 3-30).

    This pass uses:
    - initial_delta = 0.08
    - tank_geometry = "Sphere"
    - dV = 9200 m/s

    Returns:
        Tuple[EngineParams, StageParams]: A tuple containing the
        default engine and stage/mission parameters.
    """

    # Engine params from Page 3, 4, 27
    engine = EngineParams(
        thrust_vac_N=0.0,  # Calculated based on M_o and TWR
        isp_vac_s=430.0,
        chamber_pressure_Pa=6.897e6,  # 1000 psi
        propellant_type="LOX/LH2",
        cycle_type="SC",  # Assumed, not specified in PDF 1st pass
        mixture_ratio=6.0,
        expansion_ratio=30.0,
        fuel_density=DENSITY_LH2,
        oxidizer_density=DENSITY_LOX,
        safety_factor=1.2 # Assuming crewed SSTO
    )

    # Mission & Stage params from Page 3, 22, 27
    stage = StageParams(
        # --- Core Mission & Payload Inputs ---
        payload_mass_kg=5000.0,
        delta_v_ms=9200.0,

        # --- 1st Pass Sizing Inputs (Akin) ---
        initial_delta=0.08,
        initial_twr=1.3,
        num_engines=6,

        # --- Geometry Inputs (Akin) ---
        tank_geometry="Sphere",
        vehicle_diameter_m=0.0,  # 1st pass uses sphere radius, so diameter input is 0
        payload_fairing_height_m=7.0,
        intertank_fairing_height_m=7.0,
        aft_fairing_height_m=7.0,

        # --- Calculated/Placeholder Values ---
        engine=None,  # Set to None initially, assigned 'engine' obj below
        propellant_mass_kg=0.0,
        vehicle_gross_mass_kg=0.0,
        vehicle_length_m=0.0,
        stage_inert_mass_kg=0.0
    )

    stage.engine = engine
    return engine, stage


def get_akin_ssto_default_params_2nd() -> Tuple[EngineParams, StageParams]:
    """
    Returns the configuration for the SSTO 2nd Pass example,
    based on the results of the 1st Pass (Page 31-32).

    This pass uses:
    - initial_delta = 0.08
    - tank_geometry = "Cylinder"
    - vehicle_diameter_m =4 m
    - Fairing heights adjusted (Page 31)

    Returns:
        Tuple[EngineParams, StageParams]: A tuple containing the
        updated engine and stage/mission parameters.
    """
    # Mission & Stage params from Page 3, 22, 27
    engine, _ = get_akin_ssto_default_params_1st()

    # Mission & Stage params from Page 3, 22, 27
    stage = StageParams(
        # --- Core Mission & Payload Inputs ---
        payload_mass_kg=5000.0,
        delta_v_ms=9200.0,

        # --- 1st Pass Sizing Inputs (Akin) ---
        initial_delta=0.08,
        initial_twr=1.3,
        num_engines=6,

        # --- Geometry Inputs (Akin) ---
        tank_geometry="Cylinder",
        vehicle_diameter_m=4.0,
        payload_fairing_height_m=7.0,
        intertank_fairing_height_m=4.0,
        aft_fairing_height_m=5.0,

        # --- Calculated/Placeholder Values ---
        engine=None,  # Set to None initially, assigned 'engine' obj below
        propellant_mass_kg=0.0,
        vehicle_gross_mass_kg=0.0,
        vehicle_length_m=0.0,
        stage_inert_mass_kg=0.0
    )

    stage.engine = engine
    return engine, stage


def get_akin_ssto_default_params_3rd() -> Tuple[EngineParams, StageParams]:
    """
    Returns the configuration for the SSTO 3rd Pass example,
    based on the results of the 2nd Pass (Page 33-34).

    This pass uses:
    - initial_delta = 0.08323 (Calculated: 12716 kg / 157240 kg)
    - tank_geometry = "Cylinder"
    - vehicle_diameter_m = 4.2 m

    Returns:
        Tuple[EngineParams, StageParams]: A tuple containing the
        updated engine and stage/mission parameters.
    """
    engine, _ = get_akin_ssto_default_params_1st()

    # Mission & Stage params from Page 3, 22, 27
    stage = StageParams(
        # --- Core Mission & Payload Inputs ---
        payload_mass_kg=5000.0,
        delta_v_ms=9200.0,

        # --- 1st Pass Sizing Inputs (Akin) ---
        initial_delta=0.08323,
        initial_twr=1.3,
        num_engines=6,

        # --- Geometry Inputs (Akin) ---
        tank_geometry="Cylinder",
        vehicle_diameter_m=4.2,
        payload_fairing_height_m=7.0,
        intertank_fairing_height_m=4.25,
        aft_fairing_height_m=5.25,

        # --- Calculated/Placeholder Values ---
        engine=None,  # Set to None initially, assigned 'engine' obj below
        propellant_mass_kg=0.0,
        vehicle_gross_mass_kg=0.0,
        vehicle_length_m=0.0,
        stage_inert_mass_kg=0.0
    )

    return engine, stage


def get_default_engine() -> EngineParams:
    """
    Creates a set of default parameters for a *custom* rocket engine.
    """
    # 1. Get a known, high-performance engine
    engine = EngineParams(
        thrust_vac_N=1_000_000,
        isp_vac_s=320.0,
        chamber_pressure_Pa=10e6,  # 1000 psi
        propellant_type="LOX/LCH4",
        cycle_type="SC",
        mixture_ratio=6.0,
        expansion_ratio=30.0,
        fuel_density=DENSITY_LCH4,
        oxidizer_density=DENSITY_LOX,
        safety_factor=1.1  # Tizon 2017, p. 15 (non-crewed)
    )
    return engine


def default_rocket_params() -> Tuple[EngineParams, StageParams]:
    """
    Creates a set of default parameters for a *custom* rocket analysis,
    using the SSME as the engine and a larger payload.

    This is for use in `main.py` for custom analysis runs.

    Returns:
        Tuple[EngineParams, StageParams]: A tuple containing the
        SSME engine and a custom stage/mission.
    """
    # 1. Get a known, high-performance engine
    engine = get_default_engine()

    # 2. Define a new, custom mission
    stage = StageParams(
        # --- Core Mission & Payload Inputs ---
        payload_mass_kg=500.0,
        delta_v_ms=7000.0,

        # --- 1st Pass Sizing Inputs (Akin) ---
        initial_delta=0.08323,
        initial_twr=1.25,
        num_engines=4,

        # --- Geometry Inputs (Akin) ---
        tank_geometry="Cylinder",
        vehicle_diameter_m=1.6,
        payload_fairing_height_m=7.0,
        intertank_fairing_height_m=2,
        aft_fairing_height_m=3,

        # --- Calculated/Placeholder Values ---
        engine=None,  # Set to None initially, assigned 'engine' obj below
        propellant_mass_kg=0.0,
        vehicle_gross_mass_kg=0.0,
        vehicle_length_m=0.0,
        stage_inert_mass_kg=0.0
    )

    stage.engine = engine

    return engine, stage


def get_two_stage_example_params() -> TwoStageConfig:
    """
        Defines parameters for a 2-stage LOX/LCH4 vehicle (2m diameter).
        Goal: 500kg to LEO (~9000 m/s dV total).
        Total Length approx 20m (S1 ~13m, S2 ~7m).
    """

    # --- Stage 2 (Upper) ---
    # Engine: Vacuum optimized Methalox
    s2_engine = EngineParams(
        thrust_vac_N=30_000,  # ~30 kN (Small upper stage engine)
        isp_vac_s=355.0,  # Realistic for Methalox Vacuum (e.g. Raptor Vac is ~378, but small engines are less efficient)
        chamber_pressure_Pa=6e6,
        expansion_ratio=80.0,
        mixture_ratio=3.4,  # Standard Methalox
        propellant_type="LOX/LCH4",
        cycle_type="GG",  # Simpler cycle for small stage
        oxidizer_density=DENSITY_LOX,
        fuel_density=DENSITY_LCH4,
        safety_factor=1.1  # Safety factors
    )

    s2_params = StageParams(
        engine=s2_engine,
        propellant_mass_kg=0.0,
        vehicle_gross_mass_kg=0.0,
        vehicle_length_m=0.0,
        stage_inert_mass_kg=0.0,

        payload_mass_kg=500.0,  # Payload
        delta_v_ms=4650.0,  # Upper stage dV

        initial_twr=0.7,  # Vacuum stage can be < 1
        num_engines=1,

        # Geometry
        vehicle_diameter_m=1.8,
        tank_geometry="Cylinder",

        # Length Tuning (~6-7m total)
        # Fairing covers payload, so it adds length at top
        payload_fairing_height_m=3.5,
        intertank_fairing_height_m=0.1,  # Common bulkhead or short intertank
        aft_fairing_height_m=2,  # Engine section

        # Sizing Guess
        # Small stages have bad mass fractions. 0.16 = 16% inert mass.
        initial_delta=0.148
    )

    # --- Stage 1 ---
    # Engine: Sea-level optimized Methalox
    s1_engine = EngineParams(
        thrust_vac_N=350_000,  # ~350 kN per engine (x4 = 1.4MN)
        isp_vac_s=325.0,  # Lower for SL engine
        chamber_pressure_Pa=10e6,
        expansion_ratio=30.0,
        mixture_ratio=3.4,
        propellant_type="LOX/LCH4",
        cycle_type="GG",
        oxidizer_density=DENSITY_LOX,
        fuel_density=DENSITY_LCH4,
        safety_factor=1.1
    )

    s1_params = StageParams(
        engine=s1_engine,
        propellant_mass_kg=0.0,
        vehicle_gross_mass_kg=0.0,
        vehicle_length_m=0.0,
        stage_inert_mass_kg=0.0,
        payload_mass_kg=0.0,  # Overwritten by S2 mass
        delta_v_ms=4650.0,  # Booster dV

        initial_twr=1.35,  # Liftoff TWR
        num_engines=4,  # Cluster of 4

        # Geometry
        vehicle_diameter_m=1.8,  # Same diameter
        tank_geometry="Cylinder",

        # Length Tuning (~13-14m)
        # Payload fairing is on S2, so here it is 0 or just adapter
        payload_fairing_height_m=0.5,
        intertank_fairing_height_m=0.1,  # Between tanks
        aft_fairing_height_m=2,  # Thrust structure/Legs area

        # Sizing Guess
        # Boosters are more efficient. 0.08 = 8% inert mass.
        initial_delta=0.0675
    )

    return TwoStageConfig(s1_engine, s1_params, s2_engine, s2_params)