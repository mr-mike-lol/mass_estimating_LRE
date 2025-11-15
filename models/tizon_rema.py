# models/tizon_rema.py

"""
Implements the dimensionless engine mass model from Tizón & Román, 2017.

This model estimates engine mass by scaling the mass of each component
relative to a reference engine (SSME, LE-5, or RL10A) based on a set of
dimensionless parameter ratios.

Core Equation (Eq. 8 in the paper):
m_engine / m_engine_0 = SUM( alpha_i * (m_i / m_i_0) )

where:
- m_i / m_i_0 = PROD( (P_j / P_j_0) ^ a_ij )  (Eq. 7)
- alpha_i: Mass fraction coefficients (Tables 5 & 6)
- a_ij: Exponents (Tables 2 & 3)
- P_j: Engine parameters (Pc, m_dot, rt, etc.)

Reference:
Tizón, J. M., & Román, A. (2017). "A Mass Model for Liquid Propellant
Rocket Engines." 53rd AIAA/SAE/ASEE Joint Propulsion Conference.

--- Model Approach: Pros and Cons ---

Pros:
+ Component-Level Breakdown: Estimates mass for individual components
  (turbopumps, chambers, valves, etc.), enabling detailed trade studies.
+ Physics-Informed Scaling: Uses exponents (a_ij) derived from physical
  principles (e.g., hoop stress, flow rates), making it more robust
  than pure statistical regression.
+ Material Sensitivity: Explicitly accounts for changes in material
  density (rho_mat) and yield strength (sigma_y), which is critical
  for evaluating new alloys or manufacturing methods.
+ Cycle-Specific: Provides distinct coefficients (alphas) and reference
  engines for Staged Combustion (SC), Gas Generator (GG), and
  Expander (EX) cycles.

Cons:
- Reference-Dependent: The model's accuracy is tied to the similarity
  between the new engine and the reference engine (SSME, LE-5, RL10A).
  Extrapolation to vastly different sizes or architectures can be unreliable.
- Limited Propellant Baseline: All reference engines (SSME, LE-5, RL10A)
  are LOX/LH2. While propellant density is a parameter, the mass
  fractions (alphas) are derived from LH2 engines, which may not
  accurately represent Methane or Kerosene engines.
- Data Intensive: Requires a significant number of input parameters
  (Pc, Isp, O/F, expansion ratio, material properties) compared to
  simpler thrust-based models.
- `rt` Proxy: This implementation must use a proxy for the throat radius (rt)
  ratio, as `rt` is not a direct input. This introduces a minor
  assumption based on Thrust and Pc.
"""

import math
import inspect
from typing import Literal, Dict, Any, Optional

# Import from root (for constants and test functions)
import vehicle_definitions
from vehicle_definitions import (
    DENSITY_RP1, DENSITY_LH2, DENSITY_LOX, DENSITY_LCH4, G0
)

# Import from models/ (for data classes)
from models.common_params import EngineParams, CycleType, StageParams
from models.base import BaseEngineModel, ModelResult

# --- 1. Reference Engine Data ---
# Source: Tizón & Román, 2017, Table 7
REFERENCE_ENGINES: Dict[str, Dict[str, Any]] = {
    "SC": {
        "name": "SSME",
        "mass_kg": 3177.0,
        "thrust_vac_N": 2_280_000.0,
        "pc_Pa": 2.04e7,
        "rt_m": 0.138,  # Throat radius
        "m_dot_kg_s": 512.6,
        "of_ratio": 6.0,
        "expansion_ratio": 77.5,
        "prop_type": "LOX/LH2",
        "fuel_density": DENSITY_LH2,
        "oxidizer_density": DENSITY_LOX,
        "fs": 1.2  # Safety Factor, p. 15 (crewed)
    },
    "GG": {
        "name": "LE-5",
        "mass_kg": 255.0,
        "thrust_vac_N": 103_000.0,
        "pc_Pa": 3.65e6,
        "rt_m": 0.068,
        "m_dot_kg_s": 23.33,
        "of_ratio": 5.5,
        "expansion_ratio": 140.0,
        "prop_type": "LOX/LH2",
        "fuel_density": DENSITY_LH2,
        "oxidizer_density": DENSITY_LOX,
        "fs": 1.1  # Safety Factor, p. 15 (non-crewed)
    },
    "EX": {
        "name": "RL10A-3-A",
        "mass_kg": 138.0,
        "thrust_vac_N": 73_400.0,
        "pc_Pa": 3.20e6,
        "rt_m": 0.076,
        "m_dot_kg_s": 16.85,
        "of_ratio": 5.0,
        "expansion_ratio": 61.1,
        "prop_type": "LOX/LH2",
        "fuel_density": DENSITY_LH2,
        "oxidizer_density": DENSITY_LOX,
        "fs": 1.1  # Safety Factor, p. 15 (non-crewed)
    }
}

# Add bulk density and m_dot (if not specified) to reference data
for _cycle, engine in REFERENCE_ENGINES.items():
    engine['bulk_density'] = (1.0 + engine['of_ratio']) / \
                             ((1.0 / engine['fuel_density']) + (engine['of_ratio'] / engine['oxidizer_density']))
    if 'm_dot_kg_s' not in engine and 'isp_vac_s' in engine:
        engine['m_dot_kg_s'] = engine['thrust_vac_N'] / (engine['isp_vac_s'] * G0)

# --- 2. Mass Fraction Coefficients (Alpha_i) ---
# Source: Tizón & Román, 2017, Table 5 (Design)
ALPHAS_DESIGN: Dict[str, Dict[str, float]] = {
    "SC": {
        "tubes": 0.0640, "manifold": 0.1717, "jacket": 0.0470, "radiative_nozzle": 0.1566,
        "combustion_chamber": 0.1137, "gas_generator": 0.0551,
        "ox_turbopump": 0.1089, "fuel_turbopump": 0.1407,
        "ox_valve": 0.0698, "fuel_valve": 0.0690
    },
    "GG": {
        "tubes": 0.0987, "manifold": 0.0146, "jacket": 0.0589, "radiative_nozzle": 0.1969,
        "combustion_chamber": 0.1934, "gas_generator": 0.0438,
        "ox_turbopump": 0.1254, "fuel_turbopump": 0.1354,
        "ox_valve": 0.0352, "fuel_valve": 0.0352
    },
    "EX": {
        "tubes": 0.0629, "manifold": 0.2301, "jacket": 0.0210, "radiative_nozzle": 0.0,
        "combustion_chamber": 0.2007, "gas_generator": 0.0,  # Expander has no GG
        "ox_turbopump": 0.1411, "fuel_turbopump": 0.1411,
        "ox_valve": 0.0202, "fuel_valve": 0.0185
    }
}

# Source: Tizón & Román, 2017, Table 6 (Historical)
ALPHAS_HISTORICAL: Dict[str, Dict[str, float]] = {
    "SC": {
        "combustion_chamber": 0.0670, "ox_turbopump": 0.1044, "fuel_turbopump": 0.1346,
        "ox_valve": 0.0513, "fuel_valve": 0.0510, "structure": 0.1220, "auxiliary": 0.0380
    },
    "GG": {
        "combustion_chamber": 0.3545, "ox_turbopump": 0.1445, "fuel_turbopump": 0.1793,
        "ox_valve": 0.0579, "fuel_valve": 0.0579, "structure": 0.0871, "auxiliary": 0.1722
    },
    "EX": {
        "combustion_chamber": 0.2318, "ox_turbopump": 0.1411, "fuel_turbopump": 0.1411,
        "ox_valve": 0.0520, "fuel_valve": 0.0472, "structure": 0.0766, "auxiliary": 0.0879
    }
}

# --- 3. Parameter Exponents (a_ij) ---
# Source: Tizón & Román, 2017, Table 2 (Design)
# Keys: pc, expansion_ratio, rt, m_dot_prop, prop_density, rho_mat, sigma_y, fs exponents from Table 2.
EXPONENTS_DESIGN: Dict[str, Dict[str, float]] = {
    "tubes":              {"pc": 1.0, "expansion_ratio": 1.0, "rt": 2.0, "m_dot_prop": 0.0, "prop_density": 0.0, "rho_mat": 1.0, "sigma_y": -1.0, "fs": 1.0},
    "manifold":           {"pc": 1.0, "expansion_ratio": 0.0, "rt": 1.0, "m_dot_prop": 1.0, "prop_density": -1.0, "rho_mat": 1.0, "sigma_y": -1.0, "fs": 1.0},
    "jacket":             {"pc": 1.0, "expansion_ratio": 1.5, "rt": 3.0, "m_dot_prop": 0.0, "prop_density": 0.0, "rho_mat": 1.0, "sigma_y": -1.0, "fs": 1.0},
    "combustion_chamber": {"pc": 1.0, "expansion_ratio": 0.0, "rt": 2.0, "m_dot_prop": 0.0, "prop_density": 0.0, "rho_mat": 1.0, "sigma_y": -1.0, "fs": 1.0},
    "gas_generator":      {"pc": 1.0, "expansion_ratio": 0.0, "rt": 2.0, "m_dot_prop": 0.0, "prop_density": 0.0, "rho_mat": 1.0, "sigma_y": -1.0, "fs": 1.0},
    "ox_turbopump":       {"pc": 0.15, "expansion_ratio": 0.0, "rt": 0.0, "m_dot_prop": 0.9, "prop_density": -0.45, "rho_mat": 1.0, "sigma_y": -1.0, "fs": 1.0},
    "fuel_turbopump":     {"pc": 0.15, "expansion_ratio": 0.0, "rt": 0.0, "m_dot_prop": 0.9, "prop_density": -0.45, "rho_mat": 1.0, "sigma_y": -1.0, "fs": 1.0},
    "ox_valve":           {"pc": 1.0, "expansion_ratio": 0.0, "rt": 0.0, "m_dot_prop": 1.0, "prop_density": -1.0, "rho_mat": 1.0, "sigma_y": -1.0, "fs": 1.0},
    "fuel_valve":         {"pc": 1.0, "expansion_ratio": 0.0, "rt": 0.0, "m_dot_prop": 1.0, "prop_density": -1.0, "rho_mat": 1.0, "sigma_y": -1.0, "fs": 1.0},
}

# Source: Tizón & Román, 2017, Table 3 (Historical)
# Note: The paper lists "Turbo-pump" and "Valve" generically, but the
# alphas in Table 6 are split. We apply the generic exponents to both
# ox/fuel components for the historical method., rho_mat, sigma_y, fs exponents from Table 3.
EXPONENTS_HISTORICAL: Dict[str, Dict[str, float]] = {
    "combustion_chamber": {"pc": 1.0, "expansion_ratio": 0.0, "rt": 1.4, "m_dot_prop": 0.0, "prop_density": 0.0, "rho_mat": 1.0, "sigma_y": -1.0, "fs": 1.0},
    "ox_turbopump":       {"pc": 0.53, "expansion_ratio": 0.0, "rt": 0.0, "m_dot_prop": 0.53, "prop_density": -0.53, "rho_mat": 1.0, "sigma_y": -1.0, "fs": 1.0},
    "fuel_turbopump":     {"pc": 0.53, "expansion_ratio": 0.0, "rt": 0.0, "m_dot_prop": 0.53, "prop_density": -0.53, "rho_mat": 1.0, "sigma_y": -1.0, "fs": 1.0},
    "ox_valve":           {"pc": 0.3, "expansion_ratio": 0.0, "rt": 0.0, "m_dot_prop": 0.625, "prop_density": -0.625, "rho_mat": 1.0, "sigma_y": -1.0, "fs": 1.0},
    "fuel_valve":         {"pc": 0.3, "expansion_ratio": 0.0, "rt": 0.0, "m_dot_prop": 0.625, "prop_density": -0.625, "rho_mat": 1.0, "sigma_y": -1.0, "fs": 1.0},
    "structure":          {"pc": 0.92, "expansion_ratio": 0.0, "rt": 1.84, "m_dot_prop": 0.0, "prop_density": 0.0, "rho_mat": 1.0, "sigma_y": -1.0, "fs": 1.0},
    "auxiliary":          {"pc": 0.0, "expansion_ratio": 0.0, "rt": 1.0, "m_dot_prop": 0.0, "prop_density": 0.0, "rho_mat": 1.0, "sigma_y": -1.0, "fs": 1.0},
}

# Source: Tizón & Román, 2017, Table 8
# TODO: User needs to fill in the density (kg/m^3) and yield_stress (Pa)
# The values below are placeholders and for demonstration only.
MATERIAL_PROPERTIES: Dict[str, Dict[str, float]] = {
    "Inconel 718": {
        "rho_kg_m3": 8190.0,  # Placeholder, please verify
        "sigma_y_Pa": 1034e6, # Placeholder, please verify
    },
    "D6AC steel": {
        "rho_kg_m3": 7870.0,  # Placeholder, please verify
        "sigma_y_Pa": 1500e6, # Placeholder, please verify
    },
    "Al 7075 T6": {
        "rho_kg_m3": 2810.0,  # Placeholder, please verify
        "sigma_y_Pa": 503e6,  # Placeholder, please verify
    },
    # Example of a new material you might add for trade studies
    "MyNewAlloy": {
        "rho_kg_m3": 4500.0,  # e.g., A Titanium alloy
        "sigma_y_Pa": 1100e6, # e.g., A Titanium alloy
    }
}

# Maps reference engine (by cycle) components from Table 8
REFERENCE_MATERIALS: Dict[str, Dict[str, str]] = {
    "SC": { # SSME
        "nozzle": "Inconel 718",
        "combustion_chamber": "Inconel 718",
        "gas_generator": "Inconel 718",
        "turbopumps": "D6AC steel",
        "valves": "D6AC steel"
    },
    "GG": { # LE-5
        "nozzle": "Inconel 718",
        "combustion_chamber": "Inconel 718",
        "gas_generator": "Inconel 718",
        "turbopumps": "Al 7075 T6",
        "valves": "Al 7075 T6"
    },
    "EX": { # RL10A
        "nozzle": "Inconel 718",
        "combustion_chamber": "Inconel 718",
        "gas_generator": None, # No GG
        "turbopumps": "D6AC steel",
        "valves": "D6AC steel"
    }
}

# Maps the broad categories from Table 8 to the specific
# components in the Alpha tables (Tables 5 & 6)
COMPONENT_MATERIAL_MAP: Dict[str, str] = {
    # Nozzle components
    "tubes": "nozzle",
    "manifold": "nozzle",
    "jacket": "nozzle",
    "radiative_nozzle": "nozzle",
    # CC
    "combustion_chamber": "combustion_chamber",
    # GG
    "gas_generator": "gas_generator",
    # Turbopumps
    "ox_turbopump": "turbopumps",
    "fuel_turbopump": "turbopumps",
    # Valves
    "ox_valve": "valves",
    "fuel_valve": "valves",
    # Historical-only components (Not in Table 8, assume no change)
    "structure": "structure",
    "auxiliary": "auxiliary"
}


# --- TIZON ENGINE MASS MODEL (SECTION II-VI) ---

class TizonEngineModel(BaseEngineModel):
    """
    Implements the dimensionless mass model from Tizón & Román, 2017.

    This class estimates engine mass by scaling component masses relative
    to a cycle-specific reference engine.
    """

    def __init__(self, cycle_type: CycleType, method: Literal["design", "historical"] = "historical"):
        """
        Initializes the model for a specific cycle type and calculation method.

        Args:
            cycle_type (CycleType): The engine cycle ('GG', 'SC', or 'EX').
            method (Literal["design", "historical"]):
                - 'design': Uses exponents from Table 2 (design-based).
                - 'historical': Uses exponents from Table 3 (historical data regression).
        """
        if cycle_type not in REFERENCE_ENGINES:
            raise ValueError(f"Cycle type {cycle_type} not supported by Tizon model.")

        self.cycle_type = cycle_type
        self.method = method

        # 1. Select the reference engine
        self.ref_engine = REFERENCE_ENGINES[cycle_type]
        self.ref_mass_kg = self.ref_engine["mass_kg"]
        self.ref_materials = REFERENCE_MATERIALS[cycle_type]

        # 2. Select the alpha coefficients
        if method == "design":
            self.alphas = ALPHAS_DESIGN[cycle_type]
            self.exponents = EXPONENTS_DESIGN
        else:  # historical
            self.alphas = ALPHAS_HISTORICAL[cycle_type]
            self.exponents = EXPONENTS_HISTORICAL

    @property
    def model_name(self) -> str:
        """Returns the unique, human-readable name of the model."""
        return f"Tizon/RemA (2017) - {self.method.capitalize()}"

    def _calculate_param_ratios(self, params: EngineParams) -> Dict[str, Any]:
        """
        Calculates the ratios (P_j / P_j_0) for all parameters
        used in the model.
        """
        new_materials = params.materials

        # Calculate derived parameters for the 'new' engine
        m_dot_new = params.thrust_vac_N / (params.isp_vac_s * G0)
        bulk_density_new = params.bulk_density

        # --- Throat Radius (rt) Proxy Ratio ---
        # The model requires (rt / rt_0). Since EngineParams does not provide
        # rt, we must estimate it.
        # From F ~ Pc * At, we get At ~ F / Pc.
        # Since At = pi * rt^2, we have (rt/rt0)^2 = (At/At0) ~ (F/Pc) / (F0/Pc0).
        # This proxy is used in place of a direct rt measurement.
        # An alternative from Eq (12) would be (m*c* / Pc), but c* is unknown.
        rt_ratio_sq_proxy = (params.thrust_vac_N / params.chamber_pressure_Pa) / \
                            (self.ref_engine["thrust_vac_N"] / self.ref_engine["pc_Pa"])
        rt_ratio_proxy = rt_ratio_sq_proxy ** 0.5

        # Ratios for 'design' method (split turbopumps/valves)
        of_new = params.mixture_ratio
        of_ref = self.ref_engine['of_ratio']

        m_dot_ox_new = m_dot_new * (of_new / (1.0 + of_new))
        m_dot_fuel_new = m_dot_new * (1.0 / (1.0 + of_new))

        m_dot_ox_ref = self.ref_engine['m_dot_kg_s'] * (of_ref / (1.0 + of_ref))
        m_dot_fuel_ref = self.ref_engine['m_dot_kg_s'] * (1.0 / (1.0 + of_ref))

        # Prevent division by zero if ref flow is zero (e.g., placeholder)
        m_dot_ox_ref = m_dot_ox_ref if m_dot_ox_ref != 0 else 1.0
        m_dot_fuel_ref = m_dot_fuel_ref if m_dot_fuel_ref != 0 else 1.0

        # --- Material Property Ratios ---
        # Calculate rho_mat and sigma_y ratios for each component
        rho_mat_ratios: Dict[str, float] = {}
        sigma_y_ratios: Dict[str, float] = {}

        for component_name in self.alphas.keys():
            # Find the material category (e.g., "nozzle") for this component (e.g., "tubes")
            mat_category = COMPONENT_MATERIAL_MAP.get(component_name)

            if mat_category is None or mat_category not in self.ref_materials:
                # e.g., 'structure', 'auxiliary', or GG in EX cycle
                rho_mat_ratios[component_name] = 1.0
                sigma_y_ratios[component_name] = 1.0
                continue

            # 1. Get reference material properties
            ref_mat_name = self.ref_materials.get(mat_category)
            if ref_mat_name is None:  # Handle EX cycle's null GG
                rho_mat_ratios[component_name] = 1.0
                sigma_y_ratios[component_name] = 1.0
                continue

            ref_props = MATERIAL_PROPERTIES.get(ref_mat_name)
            if ref_props is None:
                raise ValueError(f"Reference material {ref_mat_name} not in MATERIAL_PROPERTIES")

            # 2. Get new material properties
            new_mat_name = None
            if new_materials:
                new_mat_name = new_materials.get(mat_category)

            # If no new material is specified for this category, use the reference material
            if new_mat_name is None:
                new_mat_name = ref_mat_name

            new_props = MATERIAL_PROPERTIES.get(new_mat_name)
            if new_props is None:
                raise ValueError(f"New material {new_mat_name} not in MATERIAL_PROPERTIES")

            # 3. Calculate and store the ratio
            rho_ref = ref_props.get("rho_kg_m3", 1.0)
            sig_ref = ref_props.get("sigma_y_Pa", 1.0)
            rho_new = new_props.get("rho_kg_m3", rho_ref)
            sig_new = new_props.get("sigma_y_Pa", sig_ref)

            rho_mat_ratios[component_name] = (rho_new / rho_ref) if rho_ref != 0 else 1.0
            sigma_y_ratios[component_name] = (sig_new / sig_ref) if sig_ref != 0 else 1.0

        return {
            "pc": params.chamber_pressure_Pa / self.ref_engine["pc_Pa"],
            "expansion_ratio": params.expansion_ratio / self.ref_engine["expansion_ratio"],
            "rt": rt_ratio_proxy,
            "fs": params.safety_factor / self.ref_engine["fs"],

            # --- Material Property Ratios ---
            # These are now dictionaries, calculated above
            "rho_mat_ratios": rho_mat_ratios,
            "sigma_y_ratios": sigma_y_ratios,

            # For 'design' model with split components
            "m_dot_ox": m_dot_ox_new / m_dot_ox_ref,
            "m_dot_fuel": m_dot_fuel_new / m_dot_fuel_ref,
            "density_ox": params.oxidizer_density / self.ref_engine["oxidizer_density"],
            "density_fuel": params.fuel_density / self.ref_engine["fuel_density"],

            # For 'historical' model with bulk components
            "m_dot_prop": m_dot_new / self.ref_engine["m_dot_kg_s"],
            "prop_density": bulk_density_new / self.ref_engine["bulk_density"],
        }

    def _calculate_component_ratio(self, component: str, param_ratios: Dict[str, Any]) -> float:
        """
        Calculates the (m_i / m_i_0) ratio for a single component using Eq. (7).

        Args:
            component (str): The name of the component (e.g., 'combustion_chamber').
            param_ratios (Dict[str, float]): Pre-calculated P_j/P_j_0 ratios.

        Returns:
            float: The (m_i / m_i_0) ratio.
        """

        if component not in self.exponents:
            # Component is in Alphas but not Exponents (e.g., 'radiative_nozzle')
            return 1.0

        exps = self.exponents[component]

        # Handle split vs. bulk properties for turbopumps and valves
        if self.method == "design":
            if component == "ox_turbopump" or component == "ox_valve":
                # Use oxidizer-specific m_dot and density
                m_dot_ratio = param_ratios["m_dot_ox"]
                density_ratio = param_ratios["density_ox"]
            elif component == "fuel_turbopump" or component == "fuel_valve":
                # Use fuel-specific m_dot and density
                m_dot_ratio = param_ratios["m_dot_fuel"]
                density_ratio = param_ratios["density_fuel"]
            else:
                # Default to bulk properties (though not used by design TPs/Valves)
                m_dot_ratio = param_ratios["m_dot_prop"]
                density_ratio = param_ratios["prop_density"]
        else:  # historical
            # Historical method uses bulk properties for all components
            m_dot_ratio = param_ratios["m_dot_prop"]
            density_ratio = param_ratios["prop_density"]

        # Calculate the product: PROD( (P_j / P_j_0) ^ a_ij )
        # Uses all exponents from Tizon 2017, Tables 2 & 3
        ratio = 1.0
        ratio *= (param_ratios["pc"] ** exps.get("pc", 0.0))
        ratio *= (param_ratios["expansion_ratio"] ** exps.get("expansion_ratio", 0.0))
        ratio *= (param_ratios["rt"] ** exps.get("rt", 0.0))
        ratio *= (m_dot_ratio ** exps.get("m_dot_prop", 0.0))
        ratio *= (density_ratio ** exps.get("prop_density", 0.0))

        # Add material property ratios
        rho_ratio = param_ratios["rho_mat_ratios"].get(component, 1.0)
        sigma_ratio = param_ratios["sigma_y_ratios"].get(component, 1.0)

        ratio *= (rho_ratio ** exps.get("rho_mat", 0.0))
        ratio *= (sigma_ratio ** exps.get("sigma_y", 0.0))

        ratio *= (param_ratios["fs"] ** exps.get("fs", 0.0))

        return ratio

    def estimate_mass(self, params: EngineParams) -> ModelResult:
        """
        Estimates the total engine mass and its component breakdown.

        Implements Eq. (8): m_engine = m_engine_0 * SUM( alpha_i * (m_i / m_i_0) )

        Args:
            params (EngineParams): The parameters of the new engine to estimate.

        Returns:
            ModelResult: A TypedDict containing 'total_mass_kg' and 'components_kg'.
        """

        # 1. Calculate all P_j / P_j_0 ratios once
        param_ratios = self._calculate_param_ratios(params)

        total_ratio_sum = 0.0
        component_masses: Dict[str, float] = {}

        # 2. Iterate over all components defined in the Alpha table
        for component, alpha_i in self.alphas.items():
            if alpha_i == 0.0:
                component_masses[component] = 0.0
                continue

            # 3. Calculate (m_i / m_i_0) for this component
            m_i_ratio = self._calculate_component_ratio(component, param_ratios)

            # 4. Add to the total sum (Eq. 8)
            total_ratio_sum += alpha_i * m_i_ratio

            # 5. Store component mass for breakdown
            # m_i = m_engine_0 * alpha_i * (m_i / m_i_0)
            component_masses[component] = self.ref_mass_kg * alpha_i * m_i_ratio

        # 6. Calculate final total mass
        total_mass = self.ref_mass_kg * total_ratio_sum

        # 7. Return the standardized ModelResult
        return ModelResult(
            total_mass_kg=total_mass,
            components_kg=component_masses,
            notes={
                "method": self.method,
                "reference_engine": self.ref_engine["name"]
            }
        )


# --- TIZON STAGE MASS MODELS (SECTION VII & VIII) ---

def calculate_propellant_mass(stage: StageParams) -> float:
    """
    Calculates the propellant mass required for a given mission,
    based on Tizón & Román, 2017, Section VII. [cite: 647]

    This implementation uses the rocket equation variant based on
    payload mass (m_pay) and inert mass fraction (delta = M_inert / M_gross),
    which is common in iterative sizing passes (like Akin's).

    M_prop = (M_pay * (MR - 1)) / (1 - MR * delta)
    where MR = exp(dV / (Isp * g0))

    Args:
        stage (StageParams): A StageParams object containing inputs:
            - delta_v_ms
            - payload_mass_kg
            - initial_delta (M_inert / M_gross)
            - engine.isp_vac_s

    Returns:
        float: The calculated propellant mass (m_prop) in kg.
    """
    if not stage.engine:
        raise ValueError("StageParams.engine must be set to calculate propellant mass.")

    ve = stage.engine.isp_vac_s * G0
    if ve == 0:
        return 0.0

    # Calculate Mass Ratio (MR) from the rocket equation [cite: 659]
    mass_ratio = math.exp(stage.delta_v_ms / ve)

    delta = stage.initial_delta
    m_pay = stage.payload_mass_kg

    # Re-arranged rocket equation using delta = M_inert / M_gross
    # This is a variation of Eq. 70 [cite: 663]
    numerator = m_pay * (mass_ratio - 1.0)
    denominator = 1.0 - (mass_ratio * delta)

    if denominator <= 0:
        raise ValueError(
            f"Rocket equation cannot be solved. "
            f"The combination of dV ({stage.delta_v_ms} m/s) and "
            f"inert fraction ({delta}) is physically impossible."
        )

    m_prop = numerator / denominator
    return m_prop


def calculate_tank_mass(stage: StageParams) -> float:
    """
    Calculates the mass of propellant tanks based on Tizón & Román, 2017,
    Section VIII. [cite: 678]

    This function calculates the mass for two separate tanks (fuel and oxidizer)
    and applies the geometry rules (Sphere or Cylinder) [cite: 684] and
    correction factor K. [cite: 696]

    Args:
        stage (StageParams): A StageParams object containing inputs:
            - propellant_mass_kg (must be calculated first)
            - engine (for densities and O/F ratio)
            - tank_geometry ("Sphere" or "Cylinder")
            - vehicle_diameter_m (if "Cylinder")
            - tank_ullage_factor [cite: 683]
            - tank_material_density_kg_m3
            - tank_material_allowable_stress_Pa
            - tank_correction_factor_K [cite: 696]

    Returns:
        float: The total calculated mass of both tanks (fuel + oxidizer) in kg.
    """
    if not stage.engine:
        raise ValueError("StageParams.engine must be set to calculate tank mass.")
    if stage.propellant_mass_kg <= 0:
        raise ValueError("StageParams.propellant_mass_kg must be calculated first.")

    # --- 1. Get parameters from Stage object ---
    m_prop = stage.propellant_mass_kg
    o_f_ratio = stage.engine.mixture_ratio
    fuel_dens = stage.engine.fuel_density
    ox_dens = stage.engine.oxidizer_density

    rho_tank = stage.tank_material_density_kg_m3
    sigma_zul = stage.tank_material_allowable_stress_Pa
    K = stage.tank_correction_factor_K
    ullage_factor = stage.tank_ullage_factor

    # --- 2. Calculate propellant volumes ---
    m_fuel = m_prop / (1.0 + o_f_ratio)  # [cite: 675]
    m_ox = m_prop * o_f_ratio / (1.0 + o_f_ratio)  # [cite: 676]

    v_fuel_net = m_fuel / fuel_dens
    v_ox_net = m_ox / ox_dens

    # Total tank volume (propellant + ullage + boil-off, etc.) [cite: 681]
    v_fuel_total = v_fuel_net * ullage_factor
    v_ox_total = v_ox_net * ullage_factor

    # --- 3. Calculate tank pressure (Eq. 76) [cite: 688] ---
    # p_tank = (10^(-0.1068 * (log10(V_tank) - 0.2588))) * 1e6
    try:
        p_tank_fuel = (10 ** (-0.1068 * (math.log10(v_fuel_total) - 0.2588))) * 1e6
        p_tank_ox = (10 ** (-0.1068 * (math.log10(v_ox_total) - 0.2588))) * 1e6
    except ValueError:
        raise ValueError("Failed to calculate tank pressure (Eq. 76). Check propellant volumes.")

    # --- 4. Calculate tank mass based on geometry ---
    total_tank_mass = 0.0

    if stage.tank_geometry == "Sphere":
        # --- Spherical Tank Calculation (Eq. 79, 80, 81) [cite: 702, 704] ---

        # Fuel Tank
        r_fuel = (3.0 * v_fuel_total / (4.0 * math.pi)) ** (1.0 / 3.0)  # [cite: 702]
        a_fuel = 4.0 * math.pi * r_fuel ** 2  # [cite: 702]
        t_fuel = (p_tank_fuel * r_fuel) / (2.0 * sigma_zul)  # [cite: 704]
        total_tank_mass += (a_fuel * t_fuel * rho_tank)

        # Oxidizer Tank
        r_ox = (3.0 * v_ox_total / (4.0 * math.pi)) ** (1.0 / 3.0)  # [cite: 702]
        a_ox = 4.0 * math.pi * r_ox ** 2  # [cite: 702]
        t_ox = (p_tank_ox * r_ox) / (2.0 * sigma_zul)  # [cite: 704]
        total_tank_mass += (a_ox * t_ox * rho_tank)

    elif stage.tank_geometry == "Cylinder":
        # --- Cylindrical Tank Calculation (Eq. 82, 83, 84) [cite: 706, 707, 709] ---
        # Assumes two separate tanks, each with spherical end caps.

        if stage.vehicle_diameter_m <= 0:
            raise ValueError("vehicle_diameter_m must be > 0 for Cylinder tanks.")

        # Fuel Tank
        r_fuel = stage.vehicle_diameter_m / 2.0
        v_caps_fuel = (4.0 / 3.0) * math.pi * r_fuel ** 3  # Volume of 2 end caps
        l_cyl_fuel = max(0, (v_fuel_total - v_caps_fuel) / (math.pi * r_fuel ** 2))  # [cite: 706]

        a_cyl_fuel = 2.0 * math.pi * r_fuel * l_cyl_fuel  # [cite: 707]
        a_sph_fuel = 4.0 * math.pi * r_fuel ** 2  # Area of 2 end caps

        t_cyl_fuel = (p_tank_fuel * r_fuel) / sigma_zul  # [cite: 709]
        t_sph_fuel = (p_tank_fuel * r_fuel) / (2.0 * sigma_zul)  # [cite: 704]

        total_tank_mass += (a_cyl_fuel * t_cyl_fuel * rho_tank) + (a_sph_fuel * t_sph_fuel * rho_tank)

        # Oxidizer Tank
        r_ox = stage.vehicle_diameter_m / 2.0
        v_caps_ox = (4.0 / 3.0) * math.pi * r_ox ** 3
        l_cyl_ox = max(0, (v_ox_total - v_caps_ox) / (math.pi * r_ox ** 2))  # [cite: 706]

        a_cyl_ox = 2.0 * math.pi * r_ox * l_cyl_ox  # [cite: 707]
        a_sph_ox = 4.0 * math.pi * r_ox ** 2

        t_cyl_ox = (p_tank_ox * r_ox) / sigma_zul  # [cite: 709]
        t_sph_ox = (p_tank_ox * r_ox) / (2.0 * sigma_zul)  # [cite: 704]

        total_tank_mass += (a_cyl_ox * t_cyl_ox * rho_tank) + (a_sph_ox * t_sph_ox * rho_tank)

    else:
        raise ValueError(f"Unknown tank_geometry: {stage.tank_geometry}")

    # --- 5. Apply Correction Factor K (Eq. 78) [cite: 696, 700] ---
    return total_tank_mass * K


# --- STANDALONE TEST RUNNERS ---


def run_full_stage_analysis(engine_params: EngineParams, stage_params: StageParams):
    """
    Runs a full analysis using functions from this file:
    1. Calculates Engine Mass (Tizon)
    2. Calculates Propellant Mass (Tizon Sec VII)
    3. Calculates Tank Mass (Tizon Sec VIII)
    """
    print("=" * 70)
    print(f"Running Full Stage Analysis (Tizon/RemA 2017)")
    print(f"  Prop:   {engine_params.propellant_type}")
    print(f"  Cycle:  {engine_params.cycle_type}")
    print(f"  dV:     {stage_params.delta_v_ms} m/s")
    print(f"  Payload: {stage_params.payload_mass_kg} kg")
    print(f"  Delta_0: {stage_params.initial_delta}")
    print("=" * 70)

    try:
        # --- 1. Calculate Engine Mass ---
        # We use 'historical' as it's generally preferred
        model = TizonEngineModel(engine_params.cycle_type, method="historical")
        engine_result = model.estimate_mass(engine_params)
        engine_mass_kg = engine_result['total_mass_kg']
        print(f"  --- 1. Engine Mass (Historical) ---")
        print(f"    Total Engine Mass (1x): {engine_mass_kg:,.1f} kg")

        # --- 2. Calculate Propellant Mass ---
        # We need the engine object inside the stage object
        stage_params.engine = engine_params
        m_prop = calculate_propellant_mass(stage_params)
        stage_params.propellant_mass_kg = m_prop  # Store for next step
        print(f"\n  --- 2. Propellant Mass (Sec VII) ---")
        print(f"    Total Propellant Mass: {m_prop:,.1f} kg")

        # --- 3. Calculate Tank Mass ---
        m_tanks = calculate_tank_mass(stage_params)
        print(f"\n  --- 3. Tank Mass (Sec VIII) ---")
        print(f"    Total Tank Mass: {m_tanks:,.1f} kg (K={stage_params.tank_correction_factor_K})")

        # --- 4. Summary ---
        m_inert_stage_only = m_tanks  # (In a real model, this would include structures, avionics, etc.)
        total_engine_mass = engine_mass_kg * stage_params.num_engines
        m_inert_total = m_inert_stage_only + total_engine_mass + stage_params.payload_mass_kg
        m_gross = m_inert_total + m_prop

        # This is the 'delta' calculated from our model
        final_delta = (m_inert_total - stage_params.payload_mass_kg) / m_gross

        print("-" * 70)
        print("  --- PRELIMINARY MASS SUMMARY ---")
        print(f"    Propellant Mass (M_p):   {m_prop:15,.1f} kg")
        print(f"    Engine Mass (N={stage_params.num_engines}):      {total_engine_mass:15,.1f} kg")
        print(f"    Tank Mass:               {m_tanks:15,.1f} kg")
        print(f"    Payload Mass:            {stage_params.payload_mass_kg:15,.1f} kg")
        print(f"    -------------------------------------")
        print(f"    GROSS MASS (M_o):        {m_gross:15,.1f} kg")
        print(f"    Inert / Gross (delta):   {final_delta:15.4f} (vs initial {stage_params.initial_delta:.4f})")


    except Exception as e:
        print(f"  ERROR during full stage analysis: {e}")
    print("=" * 70)


if __name__ == "__main__":
    """
    Standalone runner for the Tizon/RemA model.
    """
    # --- 1. CONFIGURE TEST RUN ---
    # This runs the ENGINE-ONLY analysis
    SPECIFIC_ENGINE_KEY_TO_TEST = "default"  # e.g., "ssme", "rd120", "default", or None for all

    # This runs the FULL-STAGE analysis (Engine + Propellant + Tanks)
    # Set to 'True' to run the full stage demo
    RUN_FULL_STAGE_DEMO = True

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
            print(f"[__main__] ERROR: Specific engine key '{SPECIFIC_ENGINE_KEY_TO_TEST}' not found.")
            print(f"[__main__] Available keys are: {list(available_engines.keys())}")
    else:
        # --- ELSE (All Engines) ---
        print(f"[__main__] Mode: Testing ALL {len(available_engines)} engines.")
        keys_to_test = list(available_engines.keys())

    # --- 4. Loop and run ENGINE-ONLY analysis ---
    if not keys_to_test:
        if not SPECIFIC_ENGINE_KEY_TO_TEST:
            print("[__main__] No engine getter functions found in vehicle_definitions.py.")
        else:
            print("[__main__] No analysis run due to configuration error.")

    for engine_key in keys_to_test:
        print(f"\n[__main__] Loading parameters for '{engine_key}'...")
        engine_getter = available_engines[engine_key]

        try:
            engine_params = engine_getter()

            # Instantiate and run for 'historical' method
            model_hist = TizonEngineModel(
                cycle_type=engine_params.cycle_type,
                method="historical"
            )
            # Call the INHERITED run_single_engine_analysis method
            model_hist.run_single_engine_analysis(engine_params)

            # Instantiate and run for 'design' method
            model_design = TizonEngineModel(
                cycle_type=engine_params.cycle_type,
                method="design"
            )
            # Call the INHERITED run_single_engine_analysis method
            model_design.run_single_engine_analysis(engine_params)

        except Exception as e:
            print(f"\n[__main__] ERROR processing engine '{engine_key}': {e}")

    print("-" * 70)

    # --- 5. Run FULL STAGE analysis (using default_rocket_params) ---
    if RUN_FULL_STAGE_DEMO:
        print("\n[__main__] Loading parameters for 'default_rocket_params'...")
        try:
            # Get the (engine, stage) tuple
            e_params, s_params = vehicle_definitions.default_rocket_params()
            # Run the new full stage analysis
            run_full_stage_analysis(e_params, s_params)
        except Exception as e:
            print(f"\n[__main__] ERROR processing 'default_rocket_params': {e}")

    print("-" * 70)
    print("[__main__] Standalone analysis complete.")