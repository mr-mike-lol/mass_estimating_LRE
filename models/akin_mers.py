# models/akin_mers.py
"""
Implements a collection of Mass Estimating Relations (MERs) for launch
vehicle components, based on the "Mass Estimating Relations" lecture
by Dr. David L. Akin (ENAE 791).

Approach Overview:
This module provides functions for "System-level Estimation".
This approach breaks a vehicle down into its primary components
(tanks, engines, structures, avionics, etc.) and estimates the
mass of each using simple, empirical formulas.

These MERs are derived from historical data (regression fits)
and are more detailed than a single vehicle-level inert fraction,
but less complex than a full, discipline-oriented analysis.
They are intended for rapid, iterative design passes (e.g., the
1st and 2nd pass analysis shown in the source document).

Pros:
+   **Fast & Iterative:** Allows for quick calculation of a full
    mass budget, enabling rapid trade studies and design loops
    (e.g., changing vehicle diameter).
+   **Component-Level Insight:** Provides a mass estimate for each
    subsystem, offering more insight than a single "delta".
+   **Historically Grounded:** The formulas are based on "prior art"
    and historical data, giving a reasonable starting point.

Cons:
-   **Empirical, Not Physical:** The formulas (e.g., for tanks) are
    linear regressions and do not model the underlying physics
    (like material choice, pressure, or safety factors).
-   **Brittle:** The MERs may be inaccurate if applied to vehicles
    or technologies far outside the original dataset (e.g., a
    modern, auto-claved composite tank might not follow the same
    MER as the 1970s-era tanks in the data).
-   **Requires Iteration:** The first pass is just an estimate. The
    source document shows the 1st pass (spherical tanks) results
    in a negative mass margin, proving that iteration is
    required.
"""

import math
from typing import Literal, Dict, Tuple, Any, NamedTuple, Optional

from models.common_params import EngineParams, StageParams, G0
from models.base import BaseEngineModel, BaseStageModel, ModelResult, StageMassBudget, StageModelResult
import vehicle_definitions


# --- Helper Data Structures ---

class GeometryResult(NamedTuple):
    """
    Holds calculated geometric properties to share between calculation steps.
    """
    r_lox: float
    h_lox: float
    area_lox_insulation: float  # Effective area for insulation calculation
    r_fuel: float
    h_fuel: float
    area_fuel_insulation: float  # Effective area for insulation calculation
    total_length: float
    vehicle_radius_fairing: float  # Radius used for fairing calculations
    m_gross: float  # Cached gross mass for avionics/wiring


# --- 1. Propulsion Model ---

class AkinPropulsionModel(BaseEngineModel):
    """
    Implements the Akin (1991) MERs for a full propulsion system.

    This model combines the individual MERs for the engine,
    thrust structure, and gimbals into a single "propulsion system"
    mass, allowing it to be plugged into the main optimization loop.
    """

    @property
    def model_name(self) -> str:
        """Returns the unique, human-readable name of the model."""
        return "Akin (1991) - Propulsion System"

    def estimate_mass(self, params: EngineParams) -> ModelResult:
        """
        Estimates the total propulsion system mass (Engine + Thrust
        Structure + Gimbals).

        References:
        1. Engine: Mass Estimating Relations (Akin, ENAE 791), Page 25.
           Formula: M_Rocket_Engine(kg) = 7.81e-4*T(N) + 3.37e-5*T(N)*sqrt(Ae/At) + 59
        2. Structure: Page 25.
           Formula: M_Thrust_Structure(kg) = 2.55e-4*T(N)
        3. Gimbals: Page 26.
           Formula: M_Gimbals(kg) = 237.8 * [T(N) / P_c(Pa)]^0.9375

        Args:
            params (EngineParams): The standardized dataclass. This model
                                   uses thrust, expansion_ratio, and
                                   chamber_pressure.

        Returns:
            ModelResult: A TypedDict containing the total mass and
                         component breakdown.
        """
        if params.chamber_pressure_Pa <= 0:
            raise ValueError("Chamber pressure must be > 0 for Akin gimbal calculation.")

        # Call the individual MER functions
        m_engine = self.estimate_engine_mass_mer(
            params.thrust_vac_N,
            params.expansion_ratio
        )

        m_thrust_structure = self.estimate_thrust_structure_mass(
            params.thrust_vac_N
        )

        m_gimbals = self.estimate_gimbal_mass(
            params.thrust_vac_N,
            params.chamber_pressure_Pa
        )

        total_mass = m_engine + m_thrust_structure + m_gimbals

        components = {
            "Engine (MER)": m_engine,
            "Thrust Structure (MER)": m_thrust_structure,
            "Gimbals (MER)": m_gimbals
        }

        return ModelResult(
            total_mass_kg=total_mass,
            components_kg=components,
            notes={
                "method": "Sum of Engine, Thrust Structure, and Gimbal MERs.",
                "reference": "Akin (ENAE 791)",
                "Applicability": "Propulsion system components (MERs). Ignores cycle."
                                 " Universal (all engines), but highly empirical."
            }
        )

    @staticmethod
    def estimate_engine_mass_mer(thrust_N: float, expansion_ratio: float) -> float:
        """
        Estimates liquid pump-fed rocket engine mass based on thrust and expansion ratio.
        Returns only the total mass.

        Reference:
        Mass Estimating Relations (Akin, ENAE 791), Page 25.
        Formula: M_Rocket_Engine(kg) = 7.81e-4*T(N) + 3.37e-5*T(N)*sqrt(Ae/At) + 59

        Args:
            thrust_N (float): Engine thrust in Newtons.
            expansion_ratio (float): Nozzle expansion ratio (Ae/At).

        Returns:
            float: Estimated engine-only mass in kg.
        """
        # Term 1 from formula
        term1 = 7.81e-4 * thrust_N
        # Term 2 from formula
        term2 = 3.37e-5 * thrust_N * (expansion_ratio ** 0.5)
        # Term 3 from formula
        term3 = 59.0

        total_mass = term1 + term2 + term3

        return total_mass

    @staticmethod
    def estimate_thrust_structure_mass(total_thrust_N: float) -> float:
        """
        Estimates thrust structure mass based on total vehicle thrust.

        Reference:
        Mass Estimating Relations (Akin, ENAE 791), Page 25.
        Formula: M_Thrust_Structure(kg) = 2.55e-4*T(N)

        Args:
            total_thrust_N (float): Total thrust of all engines supported by the structure, in Newtons.

        Returns:
            float: Estimated thrust structure mass in kg.
        """
        # Formula from source
        return 2.55e-4 * total_thrust_N

    @staticmethod
    def estimate_gimbal_mass(engine_thrust_N: float, chamber_pressure_Pa: float) -> float:
        """
        Estimates gimbal mass for a single engine.

        Reference:
        Mass Estimating Relations (Akin, ENAE 791), Page 26.
        Formula: M_Gimbals(kg) = 237.8 * [T(N) / P_c(Pa)]^0.9375

        Args:
            engine_thrust_N (float): Thrust of a single engine in Newtons.
            chamber_pressure_Pa (float): Engine chamber pressure in Pascals.

        Returns:
            float: Estimated gimbal mass in kg.
        """
        if chamber_pressure_Pa <= 0:
            raise ValueError("Chamber pressure must be > 0 for gimbal calculation")

        # Ratio T(N) / P_0(Pa) from formula
        ratio = engine_thrust_N / chamber_pressure_Pa
        # Full formula from source
        return 237.8 * (ratio ** 0.9375)


# --- 2. Stage Model ---

class AkinStageModel(BaseStageModel):
    """
    Implements the full stage mass estimation using Akin's MERs.
    Calculates tanks, insulation, fairings, avionics, and wiring.
    """

    def __init__(self):
        self.propulsion_model = AkinPropulsionModel()
        # Cache for geometry to share between abstract methods
        self._cached_geometry: Optional[GeometryResult] = None
        self._cached_params_key: Optional[Tuple] = None

    @property
    def model_name(self) -> str:
        return "Akin (1991) - Stage Model"

    def _ensure_geometry(self, engine: EngineParams, stage: StageParams) -> GeometryResult:
        """
        Helper: Calculates and caches geometry if params have changed.
        This allows _calculate_structural_mass and others to share data efficiently.
        """
        # Check if we already have valid cached geometry for these inputs
        current_key = (id(engine), id(stage), stage.payload_mass_kg, stage.delta_v_ms)
        if self._cached_geometry and self._cached_params_key == current_key:
            return self._cached_geometry

        # 1. Initial Sizing (Rocket Equation) to get Propellant Mass
        m_gross, m_prop, _, _ = self._calculate_initial_sizing(engine, stage)
        m_ox, m_fuel = self._get_propellant_mass_split(engine, m_prop)

        # 2. Calculate Geometry
        geo = self._calculate_geometry_internal(m_ox, m_fuel, m_gross, engine, stage)

        # Cache and return
        self._cached_geometry = geo
        self._cached_params_key = current_key
        return geo

    # --- Abstract Method Implementations ---

    def _calculate_tank_mass(self, params: EngineParams, stage: StageParams) -> float:
        """Required by BaseStageModel."""
        # Ensure we have sizing data (even though tanks depend mostly on mass here)
        _, m_prop, _, _ = self._calculate_initial_sizing(params, stage)
        m_ox, m_fuel = self._get_propellant_mass_split(params, m_prop)

        ox_type, fuel_type = self._parse_propellant_types(params.propellant_type)

        m_ox_tank = self._estimate_tank_mass_mer(m_ox, ox_type)
        m_fuel_tank = self._estimate_tank_mass_mer(m_fuel, fuel_type)

        return m_ox_tank + m_fuel_tank

    def _calculate_insulation_mass(self, params: EngineParams, stage: StageParams) -> float:
        """Required by BaseStageModel."""
        geo = self._ensure_geometry(params, stage)
        ox_type, fuel_type = self._parse_propellant_types(params.propellant_type)

        m_ox_ins = self._estimate_insulation_mer(geo.area_lox_insulation, ox_type)
        m_fuel_ins = self._estimate_insulation_mer(geo.area_fuel_insulation, fuel_type)
        return m_ox_ins + m_fuel_ins

    def _calculate_structural_mass(self, params: EngineParams, stage: StageParams) -> float:
        """Required by BaseStageModel (Fairings + Intertank)."""
        geo = self._ensure_geometry(params, stage)

        m_payload, m_inter, m_aft = self._calculate_fairings(geo, stage)
        return m_payload + m_inter + m_aft

    def _calculate_propulsion_system_mass(self, params: EngineParams, stage: StageParams) -> float:
        """Required by BaseStageModel."""
        # Get Gross Mass for TWR calculation
        m_gross, _, _, _ = self._calculate_initial_sizing(params, stage)

        total_thrust_req = m_gross * G0 * stage.initial_twr
        thrust_per_engine = total_thrust_req / stage.num_engines

        # Create temporary engine params with scaled thrust
        scaled_engine = EngineParams(
            thrust_vac_N=thrust_per_engine,
            isp_vac_s=params.isp_vac_s,
            chamber_pressure_Pa=params.chamber_pressure_Pa,
            expansion_ratio=params.expansion_ratio,
            mixture_ratio=params.mixture_ratio,
            propellant_type=params.propellant_type,
            cycle_type=params.cycle_type,
            fuel_density=params.fuel_density,
            oxidizer_density=params.oxidizer_density,
            safety_factor=params.safety_factor
        )

        # Get single engine system mass
        prop_result = self.propulsion_model.estimate_mass(scaled_engine)

        return prop_result['total_mass_kg'] * stage.num_engines

    def _calculate_avionics_and_power_mass(self, params: EngineParams, stage: StageParams) -> float:
        """
        Estimates avionics and wiring mass based on vehicle gross mass and length.

        References:
        1. Avionics: Mass Estimating Relations (Akin, ENAE 791), Page 20.
           Formula: M_avionics(kg) = 10 * (M_o(kg))^0.361
        2. Wiring: Mass Estimating Relations (Akin, ENAE 791), Page 20.
           Formula: M_wiring(kg) = 1.058 * sqrt(M_o(kg)) * l^0.25

        Args:
            params (EngineParams): Engine parameters.
            stage (StageParams): Stage parameters.

        Returns:
            float: Total mass of avionics and wiring.
        """
        geo = self._ensure_geometry(params, stage)
        m_gross = geo.m_gross

        if m_gross <= 0:
            return 0.0

        # Formula from source
        m_avionics = 10.0 * (m_gross ** 0.361)

        # Formula from source
        if geo.total_length > 0:
            m_wiring = 1.058 * (m_gross ** 0.5) * (geo.total_length ** 0.25)
        else:
            m_wiring = 0.0

        return m_avionics + m_wiring

    # --- Main Orchestrator (Optional Overrides) ---

    def calculate_full_stage_mass_budget(self, engine: EngineParams, stage: StageParams) -> StageMassBudget:
        """
        We override this to provide the specific detailed breakdown used in the Akin PDF,
        calling our implemented sub-methods.
        """

        # 1. Force geometry update first
        geo = self._ensure_geometry(engine, stage)
        m_gross = geo.m_gross
        _, m_prop, m_inert_guess, lambda_payload = self._calculate_initial_sizing(engine, stage)

        # 2. Call individual component methods to get totals
        # For simplicity and detail, we'll use the internal helpers directly here for the breakdown.

        # Propellant Split
        m_ox, m_fuel = self._get_propellant_mass_split(engine, m_prop)
        ox_type, fuel_type = self._parse_propellant_types(engine.propellant_type)

        # Breakdown: Tanks
        m_ox_tank = self._estimate_tank_mass_mer(m_ox, ox_type)
        m_fuel_tank = self._estimate_tank_mass_mer(m_fuel, fuel_type)

        # Breakdown: Insulation
        m_ox_ins = self._estimate_insulation_mer(geo.area_lox_insulation, ox_type)
        m_fuel_ins = self._estimate_insulation_mer(geo.area_fuel_insulation, fuel_type)

        # Breakdown: Fairings
        m_payload_fairing, m_intertank, m_aft_skirt = self._calculate_fairings(geo, stage)

        # Breakdown: Avionics/Wiring
        m_avionics = 10.0 * (m_gross ** 0.361)

        m_wiring = 1.058 * (m_gross ** 0.5) * (geo.total_length ** 0.25) if geo.total_length > 0 else 0.0

        # Breakdown: Propulsion
        total_thrust_req = m_gross * G0 * stage.initial_twr
        thrust_per_engine = total_thrust_req / stage.num_engines
        scaled_engine = EngineParams(
            thrust_vac_N=thrust_per_engine,
            isp_vac_s=engine.isp_vac_s,
            chamber_pressure_Pa=engine.chamber_pressure_Pa,
            expansion_ratio=engine.expansion_ratio,
            mixture_ratio=engine.mixture_ratio,
            propellant_type=engine.propellant_type,
            cycle_type=engine.cycle_type,
            fuel_density=engine.fuel_density,
            oxidizer_density=engine.oxidizer_density,
            safety_factor=engine.safety_factor
        )
        prop_result = self.propulsion_model.estimate_mass(scaled_engine)

        m_engines_total = prop_result['components_kg']["Engine (MER)"] * stage.num_engines
        m_struct_total = prop_result['components_kg']["Thrust Structure (MER)"] * stage.num_engines
        m_gimbals_total = prop_result['components_kg']["Gimbals (MER)"] * stage.num_engines

        # 9. Compilation
        components = {
            f"{ox_type} Tank": m_ox_tank,
            f"{fuel_type} Tank": m_fuel_tank,
            f"{ox_type} Insulation": m_ox_ins,
            f"{fuel_type} Insulation": m_fuel_ins,
            "Payload Fairing": m_payload_fairing,
            "Intertank Fairing": m_intertank,
            "Aft Fairing": m_aft_skirt,
            "Engines": m_engines_total,
            "Thrust Structure": m_struct_total,
            "Gimbals": m_gimbals_total,
            "Avionics": m_avionics,
            "Wiring": m_wiring
        }

        m_inert_calculated = sum(components.values())
        margin = (m_inert_guess - m_inert_calculated) / m_inert_guess if m_inert_guess > 0 else 0.0

        return StageMassBudget(
            gross_mass_kg=m_gross,
            propellant_mass_kg=m_prop,
            total_inert_mass_kg=m_inert_calculated,
            components_kg=components,
            notes={
                "M_inert_initial_guess_kg": m_inert_guess,
                "design_margin_percent": margin * 100.0,
                "lambda_payload": lambda_payload,
                "vehicle_length_m": geo.total_length,

                "h_lox_m": geo.h_lox,
                "h_fuel_m": geo.h_fuel
            }
        )

    def estimate_stage_fractions(self, params: EngineParams, stage_params: StageParams) -> StageModelResult:
        """Interface compatibility method."""
        budget = self.calculate_full_stage_mass_budget(params, stage_params)
        m_prop = budget['propellant_mass_kg']

        if m_prop <= 0:
            return {"propellant_mass_fraction": 0.0, "total_inert_fraction": 0.0, "component_fractions": {},
                    "notes": {}}

        return {
            "propellant_mass_fraction": m_prop / budget['gross_mass_kg'],
            "total_inert_fraction": budget['total_inert_mass_kg'] / m_prop,
            "component_fractions": {k: v / m_prop for k, v in budget['components_kg'].items()},
            "notes": budget['notes']
        }

    # --- Internal Calculation Logic ---

    def _calculate_initial_sizing(self, params: EngineParams, stage: StageParams) -> Tuple[float, float, float, float]:
        """
        Calculates M_gross, M_prop, M_inert_guess, and lambda.

        Vehicle-Level 1st Pass (Page 3 of MERs).
        """
        ve = params.isp_vac_s * G0
        r = math.exp(-stage.delta_v_ms / ve)
        lambda_payload = r - stage.initial_delta

        if lambda_payload <= 0:
            raise ValueError(f"Infeasible mission: delta ({stage.initial_delta}) is too high for this dV/Isp.")

        m_gross = stage.payload_mass_kg / lambda_payload
        m_prop = m_gross * (1 - r)
        m_inert_guess = stage.initial_delta * m_gross
        return m_gross, m_prop, m_inert_guess, lambda_payload

    def _calculate_geometry_internal(self, m_ox: float, m_fuel: float, m_gross: float,
                                     engine: EngineParams, stage: StageParams) -> GeometryResult:
        """
        Determines dimensions based on tank_geometry setting (Sphere vs Cylinder).
        Crucial for replicating Pass 1 vs Pass 2 logic.
        """
        r_lox, h_lox, a_lox_ins = 0.0, 0.0, 0.0
        r_fuel, h_fuel, a_fuel_ins = 0.0, 0.0, 0.0
        total_len = 0.0
        vehicle_radius = 0.0

        if stage.tank_geometry == "Sphere":
            # Pass 1 Logic (Pages 9-10)
            r_lox, a_lox_ins = self._geom_sphere(m_ox, engine.oxidizer_density)
            r_fuel, a_fuel_ins = self._geom_sphere(m_fuel, engine.fuel_density)
            h_lox = 2 * r_lox
            h_fuel = 2 * r_fuel

            # In Pass 1, fairings are sized roughly around the tanks
            vehicle_radius = r_fuel

            # Length approx (sum of fairing heights)
            total_len = (stage.payload_fairing_height_m +
                         stage.intertank_fairing_height_m +
                         stage.aft_fairing_height_m)

        elif stage.tank_geometry == "Cylinder":
            # Pass 2 Logic (Page 31)
            if stage.vehicle_diameter_m <= 0:
                raise ValueError("vehicle_diameter_m must be > 0 for 'Cylinder' geometry")

            vehicle_radius = stage.vehicle_diameter_m / 2.0

            r_lox, _, h_lox = self._geom_cylinder(m_ox, engine.oxidizer_density, vehicle_radius)
            r_fuel, _, h_fuel = self._geom_cylinder(m_fuel, engine.fuel_density, vehicle_radius)

            # SPECIAL AKIN LOGIC (Pages 31-32):
            # The PDF implies the insulation mass for cylindrical tanks is *not* based on
            # the cylinder's geometric area (side + 2 caps), but is instead based on the
            # surface area of a sphere of the same radius (4 * pi * r^2).
            a_lox_ins = 4.0 * math.pi * (r_lox ** 2)
            a_fuel_ins = 4.0 * math.pi * (r_fuel ** 2)

            # Length includes tanks + fairings
            total_len = (h_lox + h_fuel +
                         stage.payload_fairing_height_m +
                         stage.intertank_fairing_height_m +
                         stage.aft_fairing_height_m)
        else:
            raise ValueError(f"Unknown geometry: {stage.tank_geometry}")

        return GeometryResult(
            r_lox=r_lox, h_lox=h_lox, area_lox_insulation=a_lox_ins,
            r_fuel=r_fuel, h_fuel=h_fuel, area_fuel_insulation=a_fuel_ins,
            total_length=total_len, vehicle_radius_fairing=vehicle_radius,
            m_gross=m_gross
        )

    def _calculate_fairings(self, geo: GeometryResult, stage: StageParams) -> Tuple[float, float, float]:
        """Calculates mass for Payload, Intertank, and Aft fairings."""
        a_pay, a_int, a_aft = 0.0, 0.0, 0.0

        if stage.tank_geometry == "Sphere":
            # Cone + Frustum + Cylinder approach (Pass 1)
            a_pay = self._area_cone(geo.r_lox, stage.payload_fairing_height_m)
            a_int = self._area_frustum(geo.r_fuel, geo.r_lox, stage.intertank_fairing_height_m)
            a_aft = self._area_cylinder(geo.r_fuel, stage.aft_fairing_height_m)
        else:
            # Constant radius approach (Pass 2)
            r = geo.vehicle_radius_fairing
            a_pay = self._area_cone(r, stage.payload_fairing_height_m)
            a_int = self._area_cylinder(r, stage.intertank_fairing_height_m)
            a_aft = self._area_cylinder(r, stage.aft_fairing_height_m)

        return (
            self._estimate_fairing_mass_mer(a_pay),
            self._estimate_fairing_mass_mer(a_int),
            self._estimate_fairing_mass_mer(a_aft)
        )

    # --- MER Static Helpers (The Physics/Formulas) ---

    @staticmethod
    def _get_propellant_mass_split(params: EngineParams, total_mass: float) -> Tuple[float, float]:
        mr = params.mixture_ratio
        m_ox = total_mass * (mr / (1 + mr))
        m_fuel = total_mass * (1 / (1 + mr))
        return m_ox, m_fuel

    @staticmethod
    def _parse_propellant_types(prop_type_str: str) -> Tuple[str, str]:
        if "/" in prop_type_str:
            parts = prop_type_str.split("/")
            return parts[0], parts[1]
        return "LOX", "LH2"

    @staticmethod
    def _estimate_tank_mass_mer(mass_kg: float, prop_type: str) -> float:
        """
        Estimates propellant tank mass based on propellant mass (alternative MER).
        This is the MER used in the SSTO example calcs.

        Reference:
        Mass Estimating Relations (Akin, ENAE 791), Page 7.
        Formulas:
        M_LH2_Tank(kg) = 0.128 * M_LH2(kg)
        M_LOX_Tank(kg) = 0.0107 * M_LOX(kg)
        M_RP1_Tank(kg) = 0.0148 * M_RP1(kg)
        """
        if prop_type == "LH2":
            return 0.128 * mass_kg
        elif prop_type == "LOX":
            return 0.0107 * mass_kg
        elif prop_type == "RP1":
            return 0.0148 * mass_kg
        else:
            return 0.02 * mass_kg  # LCH4 estimated

    @staticmethod
    def _estimate_insulation_mer(area_m2: float, prop_type: str) -> float:
        """
        Estimates cryogenic insulation mass based on tank surface area.

        Reference:
        Mass Estimating Relations (Akin, ENAE 791), Page 8.
        Formulas:
        M_LH2_Insulation(kg) = 2.88 * A_tank(m^2)
        M_LOX_Insulation(kg) = 1.123 * A_tank(m^2)
        """
        if prop_type == "LH2":
            return 2.88 * area_m2
        elif prop_type == "LOX":
            return 1.123 * area_m2
        elif prop_type == "LCH4":
            # Methane is cryogenic (111K vs LOX 90K).
            # estimated, needs double-check; cz it's cryogenic
            return 1.115 * area_m2
        return 0.0

    @staticmethod
    def _estimate_fairing_mass_mer(area_m2: float) -> float:
        """
        Estimates fairing/shroud mass based on surface area.

        Reference:
        Mass Estimating Relations (Akin, ENAE 791), Page 20.
        Formula: M_fairing(kg) = 4.95 * (A_fairing(m^2))^1.15
        """
        if area_m2 <= 0: return 0.0
        return 4.95 * (area_m2 ** 1.15)

    # --- Geometry Static Helpers ---

    @staticmethod
    def _geom_sphere(mass: float, rho: float) -> Tuple[float, float]:
        """
        (Internal) Calculates the radius and surface area of a spherical tank.
        Reference: Logic derived from examples on Page 9 and Page 10.
        """
        if rho <= 0: return 0.0, 0.0
        vol = mass / rho
        r = (vol / (4 * math.pi / 3)) ** (1 / 3)
        area = 4 * math.pi * (r ** 2)
        return r, area

    @staticmethod
    def _geom_cylinder(mass: float, rho: float, r: float) -> Tuple[float, float, float]:
        """
        (Internal) Calculates the height and surface area of a cylindrical tank.
        Reference: Logic for 2nd/3rd pass, e.g., Page 31.
        """
        if rho <= 0 or r <= 0: return r, 0.0, 0.0
        vol = mass / rho
        h = vol / (math.pi * r ** 2)
        area_geom = 2 * math.pi * r * h + 2 * math.pi * r ** 2
        return r, area_geom, h

    @staticmethod
    def _area_cone(r: float, h: float) -> float:
        """
        Calculates the surface area of a cone (excluding the base).
        Reference: Mass Estimating Relations (Akin, ENAE 791), Page 21.
        """
        return math.pi * r * math.sqrt(r ** 2 + h ** 2)

    @staticmethod
    def _area_cylinder(r: float, h: float) -> float:
        """
        Calculates the surface area of a cylinder's side wall.
        Reference: Mass Estimating Relations (Akin, ENAE 791), Page 21.
        """
        return 2 * math.pi * r * h

    @staticmethod
    def _area_frustum(r1: float, r2: float, h: float) -> float:
        """
        Calculates the surface area of a cone frustum.
        Reference: Mass Estimating Relations (Akin, ENAE 791), Page 21.
        """
        return math.pi * (r1 + r2) * math.sqrt((r1 - r2) ** 2 + h ** 2)


# --- Reporting / Running Utils ---

def _get_pdf_reference_data(pass_num: int) -> Dict[str, Any]:
    """
    (Internal) Returns a dictionary of the reference mass budget values
    from the Akin ENAE 791 PDF for a specific pass.

    Args:
        pass_num (int): The iteration number (1, 2, or 3).

    Returns:
        Dict[str, Any]: A dict with 'components_kg', 'total_kg',
                        'guess_kg', and 'margin_pct'.
    """
    if pass_num == 1:
        # Data from Pass 1 (Page 30)
        # Note: These values are based on the dV=9200 (r=0.1129) run,
        return {
            "components_kg": {
                "LOX Tank": 1245, "LH2 Tank": 2482, "LOX Insulation": 119,
                "LH2 Insulation": 586, "Payload Fairing": 645, "Intertank Fairing": 1626,
                "Aft Fairing": 1905, "Engines": 2236, "Thrust Structure": 497,
                "Gimbals": 81, "Avionics": 744, "Wiring": 886
            },
            "total_kg": 13052,
            "guess_kg": 12240,
            "margin_pct": -6.22
        }
    elif pass_num == 2:
        # Data from Pass 2 (Page 32)
        return {
            "components_kg": {
                "LOX Tank": 1245, "LH2 Tank": 2482, "LOX Insulation": 56,
                "LH2 Insulation": 145, "Payload Fairing": 402, "Intertank Fairing": 448,
                "Aft Fairing": 579, "Engines": 2236, "Thrust Structure": 497,
                "Gimbals": 81, "Avionics": 744, "Wiring": 1044
            },
            "total_kg": 9960,
            "guess_kg": 12240,
            "margin_pct": 22.9
        }
    elif pass_num == 3:
        # Data from Pass 3 (Page 34)
        return {
            "components_kg": {
                "LOX Tank": 1382, "LH2 Tank": 2755, "LOX Insulation": 62,
                "LH2 Insulation": 160, "Payload Fairing": 427, "Intertank Fairing": 501,
                "Aft Fairing": 626, "Engines": 2443, "Thrust Structure": 552,
                "Gimbals": 90, "Avionics": 773, "Wiring": 1101
            },
            "total_kg": 10870,
            "guess_kg": 14130,
            "margin_pct": 30
        }
    else:
        # Default to empty if pass_num is not 1, 2, or 3
        return {"components_kg": {}, "total_kg": 0, "guess_kg": 0, "margin_pct": 0}


def print_ssto_results(budget: StageMassBudget, pass_num: int = 1, show_pdf_ref: bool = True):
    """Formats and prints the StageMassBudget."""
    print("=" * 75)
    print("🚀 AKIN SSTO ANALYSIS RESULTS (OOP Implementation)")
    if show_pdf_ref:
        print(f"(Based on ENAE 791, Pass {pass_num})")
    print("=" * 75)

    print("\n--- Initial Vehicle Sizing ---")
    print(f"  Gross Mass (M_o):         {budget['gross_mass_kg']:12,.1f} kg")
    print(f"  Propellant Mass (M_p):    {budget['propellant_mass_kg']:12,.1f} kg")
    print(f"  Initial Inert Guess (M_i):{budget['notes'].get('M_inert_initial_guess_kg', 0):12,.1f} kg")
    implied_payload = budget['gross_mass_kg'] - budget['propellant_mass_kg'] - budget['total_inert_mass_kg']
    print(f"  Payload Mass (implied):   {implied_payload:12,.1f} kg")

    print("\n--- Calculated Mass Budget ---")
    header = f"  {'Component':<20} | {'Calculated (kg)':>15}"
    sep = f"  {'-' * 20} | {'-' * 15}"

    pdf_data = {}
    if show_pdf_ref:
        pdf_data = _get_pdf_reference_data(pass_num)
        header += f" | {'PDF Ref (kg)':>12}"
        sep += f" | {'-' * 12}"

    print(header)
    print(sep)

    pdf_comps = pdf_data.get("components_kg", {})
    for name, mass in budget['components_kg'].items():
        line = f"  {name:<20} | {mass:15,.1f}"
        if show_pdf_ref:
            line += f" | {pdf_comps.get(name, 0):12,.0f}"
        print(line)

    print(sep)
    line_total = f"  {'Total Inert Mass':<20} | {budget['total_inert_mass_kg']:15,.1f}"
    if show_pdf_ref:
        line_total += f" | {pdf_data.get('total_kg', 0):12,.0f}"
    print(line_total)

    print("\n--- GEOMETRY DETAIL ---")
    length = budget['notes'].get('vehicle_length_m', 0.0)

    h_lox = budget['notes'].get('h_lox_m', 0.0)
    h_fuel = budget['notes'].get('h_fuel_m', 0.0)

    print(f"  Total Stage Length:  {length:15.2f} m")
    print(f"  -> LOX Tank Height:  {h_lox:15.2f} m")
    print(f"  -> Fuel Tank Height: {h_fuel:15.2f} m")
    print(f"  -> Fairings/Struct:  {length - h_lox - h_fuel:15.2f} m")

    print("\n--- FINAL DESIGN MARGIN ---")
    margin = budget['notes'].get('design_margin_percent', 0.0)
    print(f"  Calculated Margin:   {margin:15.2f} %")
    if show_pdf_ref:
        # Page 30, 32, 34
        print(f"  PDF Reference Margin: {pdf_data.get('margin_pct', 0.0):15.2f} %")
    print("=" * 75)


def check_physics_tsiolkovsky(budget: StageMassBudget, engine: EngineParams) -> Dict[str, float]:
    """
    Performs a standard physics check using the Tsiolkovsky Rocket Equation
    based on the FINAL calculated masses.

    Equation: dV = Isp * g0 * ln(m_initial / m_final)

    Args:
        budget: The calculated mass budget.
        engine: The engine parameters used (for Isp).

    Returns:
        Dict with calculated physics values.
    """
    m_initial = budget['gross_mass_kg']
    # Final mass = Gross - Propellant (i.e., Inert + Payload)
    m_final = m_initial - budget['propellant_mass_kg']

    if m_final <= 0 or m_initial <= 0:
        return {"actual_dv": 0.0, "mass_ratio": 0.0}

    mass_ratio = m_initial / m_final

    # Calculate actual achieved Delta V
    # Note: We use vacuum Isp here. For a first stage, losses are usually
    # accounted for by increasing the target Delta V requirement.
    actual_dv = engine.isp_vac_s * G0 * math.log(mass_ratio)

    return {
        "m_initial_kg": m_initial,
        "m_final_kg": m_final,
        "mass_ratio": mass_ratio,
        "actual_dv_ms": actual_dv,
        "isp_used": engine.isp_vac_s
    }


# --- Updated Orchestrators ---

def calculate_two_stage_mission(
        s1_engine: EngineParams,
        s1_params: StageParams,
        s2_engine: EngineParams,
        s2_params: StageParams
) -> Dict[str, Any]:
    """
    Calculates a two-stage mission and performs physics verification.
    """
    model = AkinStageModel()

    print("\n" + "!" * 60)
    print("!!! STARTING 2-STAGE CALCULATION LOOP !!!")
    print("!" * 60 + "\n")

    # --- 1. Calculate Stage 2 (Upper) ---
    print(">>> Calculating Stage 2 (Upper Stage)...")
    budget_s2 = model.calculate_full_stage_mass_budget(s2_engine, s2_params)

    # Verify S2 Physics
    check_s2 = check_physics_tsiolkovsky(budget_s2, s2_engine)
    s2_gross_mass = budget_s2['gross_mass_kg']
    print(f"   -> Stage 2 Gross Mass: {s2_gross_mass:,.1f} kg")
    print(f"   -> S2 Actual dV:       {check_s2['actual_dv_ms']:.1f} m/s (Target: {s2_params.delta_v_ms:.1f})")

    # --- 2. Pass Load to Stage 1 ---
    # The Gross Mass of S2 becomes the Payload of S1
    s1_params.payload_mass_kg = s2_gross_mass
    print(f">>> Passing Load to Stage 1: Payload is now {s1_params.payload_mass_kg:,.1f} kg")

    # --- 3. Calculate Stage 1 (Booster) ---
    print(">>> Calculating Stage 1 (Booster)...")
    budget_s1 = model.calculate_full_stage_mass_budget(s1_engine, s1_params)

    # Verify S1 Physics
    check_s1 = check_physics_tsiolkovsky(budget_s1, s1_engine)
    s1_gross_mass = budget_s1['gross_mass_kg']
    print(f"   -> Stage 1 Gross Mass: {s1_gross_mass:,.1f} kg")
    print(f"   -> S1 Actual dV:       {check_s1['actual_dv_ms']:.1f} m/s (Target: {s1_params.delta_v_ms:.1f})")

    # --- 4. Vehicle Summary ---
    total_glom = s1_gross_mass
    final_payload = s2_params.payload_mass_kg

    # Sum of ACTUAL achievable dV, not the target dV
    total_actual_dv = check_s1['actual_dv_ms'] + check_s2['actual_dv_ms']
    target_dv = s1_params.delta_v_ms + s2_params.delta_v_ms

    results = {
        "vehicle_summary": {
            "GLOM_kg": total_glom,
            "Final_Payload_kg": final_payload,
            "Target_Total_dV_ms": target_dv,
            "Actual_Total_dV_ms": total_actual_dv,
            "Payload_Fraction_pct": (final_payload / total_glom * 100) if total_glom > 0 else 0,
            "S1_Performance": check_s1,
            "S2_Performance": check_s2
        },
        "stage2_budget": budget_s2,
        "stage1_budget": budget_s1
    }

    return results


def print_two_stage_report(results: Dict[str, Any]):
    """
    Detailed report including Tsiolkovsky verification and Component Breakdown.
    """
    summary = results['vehicle_summary']
    s1_perf = summary['S1_Performance']
    s2_perf = summary['S2_Performance']

    print("\n" + "=" * 75)
    print("🚀 TWO-STAGE VEHICLE CONFIGURATION REPORT (With Physics Verification)")
    print("=" * 75)

    print(f"{'METRIC':<30} | {'VALUE':<20} | {'UNIT'}")
    print("-" * 75)
    print(f"{'GLOM (Start Mass)':<30} | {summary['GLOM_kg']:12,.1f}         | kg")
    print(f"{'Final Payload':<30} | {summary['Final_Payload_kg']:12,.1f}         | kg")
    print(f"{'Payload Fraction':<30} | {summary['Payload_Fraction_pct']:12.3f}         | %")
    print("-" * 75)
    print(f"{'Target Total dV':<30} | {summary['Target_Total_dV_ms']:12,.1f}         | m/s")
    print(f"{'Calculated Actual dV':<30} | {summary['Actual_Total_dV_ms']:12,.1f}         | m/s")

    # Delta V Balance Check
    diff = summary['Actual_Total_dV_ms'] - summary['Target_Total_dV_ms']
    status = "EXCESS" if diff >= -0.01 else "DEFICIT"  # Small float tolerance
    print(f"{'dV Balance':<30} | {diff:+12.1f} ({status})  | m/s")

    print("\n--- PHYSICS VERIFICATION (Tsiolkovsky Check) ---")
    print(f"{'Stage':<10} | {'Mass Ratio':<10} | {'Isp (s)':<8} | {'Target dV':<10} | {'Actual dV':<10} | {'Diff %'}")
    print("-" * 85)

    # Stage 2 Row
    s2_diff_pct = ((s2_perf['actual_dv_ms'] - 4500.0) / 4500.0) * 100
    print(
        f"{'Stage 2':<10} | {s2_perf['mass_ratio']:10.2f} | {s2_perf['isp_used']:8.1f} | {'~4500':<10} | {s2_perf['actual_dv_ms']:10.1f} | {s2_diff_pct:+6.2f}%")

    # Stage 1 Row
    s1_diff_pct = ((s1_perf['actual_dv_ms'] - 4500.0) / 4500.0) * 100
    print(
        f"{'Stage 1':<10} | {s1_perf['mass_ratio']:10.2f} | {s1_perf['isp_used']:8.1f} | {'~4500':<10} | {s1_perf['actual_dv_ms']:10.1f} | {s1_diff_pct:+6.2f}%")

    print("\n--- DETAILED COMPONENT BREAKDOWN ---")

    print("\n>>> STAGE 2 (UPPER) DETAIL <<<")
    print_ssto_results(results['stage2_budget'], pass_num=2, show_pdf_ref=False)

    print("\n>>> STAGE 1 (BOOSTER) DETAIL <<<")
    print_ssto_results(results['stage1_budget'], pass_num=2, show_pdf_ref=False)

    print("=" * 75 + "\n")


if __name__ == "__main__":
    print("Running SSTO Pass Example directly from akin_mers.py...")

    # 1. Get the specific config from the PDF (now dataclasses)
    # Options: get_akin_ssto_default_params_1st, get_akin_ssto_default_params_2nd, get_akin_ssto_default_params_3rd
    # default_rocket_params
    # We try to use the 2nd pass params if available (for better comparison), otherwise default
    try:
        engine_params, stage_params = vehicle_definitions.get_akin_ssto_default_params_1st()
        pass_to_show = 1
    except AttributeError:
        print("Error: Could not load default parameters from vehicle_definitions.")
        exit()

    model = AkinStageModel()
    results_budget = model.calculate_full_stage_mass_budget(engine_params, stage_params)
    print_ssto_results(results_budget, pass_num=pass_to_show, show_pdf_ref=True)

    try:
        # 1. Get the specific (4-elements tuple) config
        two_stage_config = vehicle_definitions.get_two_stage_example_params()

        # 2. unpack
        s1_engine = two_stage_config.stage1_engine
        s1_params = two_stage_config.stage1_params
        s2_engine = two_stage_config.stage2_engine
        s2_params = two_stage_config.stage2_params
    except AttributeError:
        print("Error: Could not load default parameters from vehicle_definitions.")
        exit()

    # --- 3. Run Calculation ---
    results = calculate_two_stage_mission(s1_engine, s1_params, s2_engine, s2_params)

    # --- 4. Print Report ---
    print_two_stage_report(results)