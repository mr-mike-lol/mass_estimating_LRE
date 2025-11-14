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
from models.common_params import EngineParams, CycleType
from typing import Literal, Dict, Any, Optional
from vehicle_definitions import G0, DENSITY_RP1, DENSITY_LH2, DENSITY_LOX, DENSITY_LCH4
import inspect
import vehicle_definitions

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
# Source: Tizón & Román, 2017, Table 5 (Design) [cite: 618]
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

# Source: Tizón & Román, 2017, Table 6 (Historical) [cite: 620]
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


class TizonEngineModel:
    """
    Implements the dimensionless mass model from Tizón & Román, 2017. [cite: 239, 266]

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
        # An alternative from Eq (12) [cite: 361] would be (m*c* / Pc), but c* is unknown.
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
        # Uses all exponents from Tizon 2017, Tables 2 & 3 [cite: 603, 605]
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

    def estimate_total_mass(self, params: EngineParams) -> Dict[str, Any]:
        """
        Estimates the total engine mass and its component breakdown.

        Implements Eq. (8): m_engine = m_engine_0 * SUM( alpha_i * (m_i / m_i_0) )

        Args:
            params (EngineParams): The parameters of the new engine to estimate.

        Returns:
            dict: A dictionary containing 'total_mass_kg' and 'components_kg'.
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

        return {
            "total_mass_kg": total_mass,
            "components_kg": component_masses,
            "notes": {
                "model": "Tizon/RemA (2017)",
                "method": self.method,
                "reference_engine": self.ref_engine["name"]
            }
        }


def run_single_engine_analysis(params: EngineParams):
    """
    Runs a full analysis for a single engine using both "design"
    and "historical" methods from Tizon/RemA (2017) and prints the results.

    Args:
        params (EngineParams): The engine parameters to analyze.
    """
    print("=" * 70)
    print(f"Running Tizon/RemA (2017) Analysis for:")
    print(f"  Thrust: {params.thrust_vac_N / 1_000_000:.3f} MN")
    print(f"  Prop:   {params.propellant_type}")
    print(f"  Cycle:  {params.cycle_type}")
    print(f"  Pc:     {params.chamber_pressure_Pa / 1e6:.2f} MPa")
    print(f"  Ae/At:  {params.expansion_ratio:.1f}")
    print(f"  FS:     {params.safety_factor}")
    print(f"  Mats:   {params.materials}")
    print("=" * 70)

    methods: list[Literal["design", "historical"]] = ["design", "historical"]

    try:
        for method in methods:
            # 1. Instantiate the model for this method
            model = TizonEngineModel(params.cycle_type, method=method)

            # 2. Estimate the mass
            result = model.estimate_total_mass(params)
            total = result['total_mass_kg']

            # 3. Print formatted results
            print(f"  --- Method: '{method}' (Reference: {model.ref_engine['name']}) ---")

            # Print component breakdown
            if total > 0:
                for component, mass in result['components_kg'].items():
                    percent = (mass / total * 100)
                    print(f"    {component:<30}: {mass:10,.1f} kg ({percent:5.1f}%)")

            # Print total
            print(f"    {'-' * 30}: {'-' * 10} -- {'-' * 7}")
            print(f"    {'TOTAL MASS':<30}: {total:10,.1f} kg (100.0%)")
            print("-" * 70)

    except ValueError as e:
        # Catch errors (e.g., "Cycle type X not supported")
        print(f"  ERROR processing {params.cycle_type} engine: {e}")
        print("=" * 70)


if __name__ == "__main__":
    """
    Standalone runner for the Tizon/RemA model.
    (Следует принципу zandbergen_engine.py)

    NOTE: This is intended for structured testing via `main.py`.
    Running this file directly will fail due to relative imports
    (e.g., `from .common_params ...`).
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
        print(f"\n[__main__] Loading parameters for '{engine_key}'...")
        engine_getter = available_engines[engine_key]

        try:
            engine_params = engine_getter()
            # Run the full analysis for this engine
            run_single_engine_analysis(engine_params)
        except Exception as e:
            print(f"\n[__main__] ERROR processing engine '{engine_key}': {e}")

    print("-" * 70)
    if keys_to_test:
        print(f"[__main__] Standalone analysis complete for {len(keys_to_test)} engine(s).")
    else:
        print("[__main__] No analysis was run.")