# models/kibbey_stage.py
"""
Provides an implementation of the "Next-Order Mass Model" for launch
vehicle stage mass estimation, as described by T. Kibbey (2015).

Approach Overview:
This module implements a semi-empirical scaling model, not a full
bottom-up physics model. It estimates the inert mass fraction of a
new stage by scaling from a well-defined reference stage (the Atlas V).

The model breaks the stage inert mass (relative to propellant mass)
[cite_start]into three key components[cite: 365]:
1.  f_i,E&S: Engine & Structure (thrust-dependent)
2.  f_i,ot: Oxidizer Tank (volume-dependent)
3.  f_i,ft: Fuel Tank (volume- and load-dependent)

Each component is scaled from the reference stage's values based on the
new stage's propellant properties (density, O/F ratio) and
[cite_start]mission parameters (T/W ratio, mass ratio)[cite: 333].

Pros:
+   Fast: Purely analytical, making it ideal for rapid trade studies
    [cite_start]and exploring large design spaces[cite: 628].
+   More Accurate (than Rho-Isp): Specifically designed to correct the
    flaws of simpler models, which fail to capture the structural
    [cite_start]penalties of low-density propellants like LH2[cite: 326, 437].
+   Captures Key Effects: Models how tank mass scales with propellant
    [cite_start]density and how the fuel tank mass is affected by launch loads[cite: 400, 404].

Cons:
-   Reference-Dependent: As a scaling model, its accuracy diminishes as
    the new design diverges from the reference (Atlas V).
-   [cite_start]Limited Scope: Relies on a simple correlation for engine T/W [cite: 373]
    [cite_start]and does not inherently model cost or size-scaling effects[cite: 637, 638].
"""

from models.common_params import EngineParams, StageParams
from models.common_params import (
    DENSITY_RP1, DENSITY_LH2, DENSITY_LOX, DENSITY_LCH4, G0
)
from models.base import BaseStageModel, StageModelResult
from typing import Dict, Any

# For standalone running
import inspect
import vehicle_definitions


class KibbeyStageModel(BaseStageModel):
    """
    [cite_start]Implements the 'Next-Order Mass Model' from Kibbey (2015) [cite: 331-332].

    This model estimates stage inert mass by scaling three components
    (Engine/Structure, Ox Tank, Fuel Tank) from a calibrated
    [cite_start]reference stage (Atlas V, LOX/RP1)[cite: 327, 342].
    """

    def __init__(self):
        """
        Initializes the model with reference constants from the Atlas V,
        [cite_start]as defined in Section III of the paper [cite: 342-360].
        """

        # --- Reference Constants (Atlas V) ---

        # [cite_start]Controlling constants [cite: 342-343]
        # f_i,tank,0: Inert fraction of *both* tanks / total propellant mass
        self.F_I_TANK_0 = 0.02  # [cite: 342]
        # rF: Multiplier on engine mass for engine-and-structure mass
        self.R_F = 2.735  # [cite: 343]

        # [cite_start]Independent variables (Mission) [cite: 350-353]
        self.R_0 = 6.62  # R: Initial to final mass ratio [cite: 352]
        self.PSI_0 = 1.28  # Psi: Initial vehicle vacuum thrust-to-weight [cite: 353]
        # [cite_start]R_mp: Propellant-to-Gross-Mass Ratio, derived from R [cite: 155]
        self.R_MP_0 = 1.0 - (1.0 / self.R_0)

        # [cite_start]Independent variables (Propellant) [cite: 354-360]
        self.OF_0 = 2.7  # O/F: Oxidizer to fuel mass ratio [cite: 355]
        self.FW_ENG_0 = 78.0  # Engine thrust-to-weight ratio [cite: 360]

        # Reference propellant densities (LOX/RP1)
        self.RHO_OX_0 = DENSITY_LOX  # kg/m^3
        self.RHO_FU_0 = DENSITY_RP1  # kg/m^3

        # Reference bulk density, calculated
        self.RHO_BULK_0 = self._calculate_bulk_density(self.RHO_OX_0, self.RHO_FU_0, self.OF_0)

        # Derived reference tank "specific inert fractions"
        # These are the f'i constants (mass_tank / mass_prop_in_tank)
        # derived from F_I_TANK_0 using Eq. [cite_start]15-16 [cite: 382-387].
        self.F_PRIME_I_OT_0 = 0.0180  # [cite: 388]
        self.F_PRIME_I_FT_0 = 0.0255  # [cite: 388]

    @property
    def model_name(self) -> str:
        """Returns the unique, human-readable name of the model."""
        return "Kibbey (2015) Stage Model"

    def _calculate_bulk_density(self, rho_ox: float, rho_fu: float, of_ratio: float) -> float:
        """Helper function to calculate propellant bulk density."""
        total_parts = 1.0 + of_ratio
        vol_fuel = 1.0 / rho_fu
        vol_ox = of_ratio / rho_ox
        if (vol_fuel + vol_ox) == 0:
            return 0.0
        return total_parts / (vol_fuel + vol_ox)

    def _calculate_bulk_density_ratio(self, params: EngineParams) -> float:
        """
        Calculates the new-to-reference bulk density ratio (r_rho).
        [cite_start]Implements Eq. 18[cite: 397].
        """
        rho_bulk_new = params.bulk_density
        if self.RHO_BULK_0 == 0:
            return 1.0
        return rho_bulk_new / self.RHO_BULK_0

    def _calculate_fi_E_S(self, params: EngineParams, r_new: float, psi_new: float) -> float:
        """
        Calculates the engine and structure inert fraction (f_i,E&S).
        Implements Eq. [cite_start]13 and 14[cite: 373, 376].

        Args:
            params (EngineParams): Parameters of the *new* engine.
            r_new (float): Mass ratio (R) of the *new* stage.
            psi_new (float): T/W (Psi) of the *new* stage.

        Returns:
            float: f_i,E&S (Engine/Structure mass / Total Propellant mass)
        """

        # [cite_start]1. Correlate new engine T/W (FW_Eng) with bulk density per Eq. 13 [cite: 373]
        rho_ratio = params.bulk_density / self.RHO_BULK_0
        if rho_ratio < 0:
            rho_ratio = 0  # Avoid issues with negative densities if inputs are bad
        fw_eng_new = 39.83 * rho_ratio + 38.17

        # 2. Calculate new propellant-to-gross-mass ratio (R_mp)
        if r_new <= 0:
            raise ValueError("Mission mass ratio (R) must be > 0.")
        if r_new == 1.0:
            # Avoid division by zero if R=1.0 (no propellant)
            r_mp_new = 0.0
        else:
            r_mp_new = 1.0 - (1.0 / r_new)

        if r_mp_new == 0 or fw_eng_new == 0:
            return 0.0

        # 3. Calculate f_i,E&S using Eq. [cite_start]14 [cite: 376]
        f_i_E_S = self.R_F * (psi_new / fw_eng_new) * (1.0 / r_mp_new)
        return f_i_E_S

    def _calculate_fi_ot(self, params: EngineParams) -> float:
        """
        Calculates the oxidizer tank inert fraction (f_i,ot).
        [cite_start]Implements Eq. 20[cite: 400].

        Args:
            params (EngineParams): Parameters of the *new* engine.

        Returns:
            float: f_i,ot (Oxidizer Tank mass / Total Propellant mass)
        """

        # [cite_start]1. Find oxidizer density ratio (r_rho_ox) per Eq. 17b [cite: 391]
        if self.RHO_OX_0 <= 0:
            raise ValueError("Reference oxidizer density must be > 0.")
        r_rho_ox = params.oxidizer_density / self.RHO_OX_0
        if r_rho_ox <= 0:
            raise ValueError("New oxidizer density must be > 0.")

        # 2. Calculate f_i,ot using Eq. [cite_start]20 (scales with 1/density) [cite: 400]
        of_new = params.mixture_ratio
        if (of_new + 1.0) == 0:
            return 0.0  # Avoid division by zero

        f_i_ot = self.F_PRIME_I_OT_0 * (1.0 / r_rho_ox) * (of_new / (of_new + 1.0))
        return f_i_ot

    def _calculate_fi_ft(self, params: EngineParams, r_new: float, psi_new: float) -> float:
        """
        Calculates the fuel tank inert fraction (f_i,ft),
        [cite_start]which depends on volume AND launch loads. Implements Eq. 25[cite: 425].

        Args:
            params (EngineParams): Parameters of the *new* engine.
            r_new (float): Mass ratio (R) of the *new* stage.
            psi_new (float): T/W (Psi) of the *new* stage.

        Returns:
            float: f_i,ft (Fuel Tank mass / Total Propellant mass)
        """

        # [cite_start]1. Calculate new bulk density ratio (r_p) per Eq. 18 [cite: 397]
        r_p = self._calculate_bulk_density_ratio(params)

        # 2. Calculate new fuel density ratio (r_rho_fu) per Eq. [cite_start]17a [cite: 390]
        if self.RHO_FU_0 <= 0:
            raise ValueError("Reference fuel density must be > 0.")
        r_rho_fu = params.fuel_density / self.RHO_FU_0
        if r_rho_fu <= 0:
            raise ValueError("New fuel density must be > 0.")

        # 3. Calculate new propellant-to-gross-mass ratio (R_mp)
        if r_new <= 0:
            raise ValueError("Mission mass ratio (R) must be > 0.")
        if r_new == 1.0:
            r_mp_new = 0.0
        else:
            r_mp_new = 1.0 - (1.0 / r_new)

        of_new = params.mixture_ratio
        if (of_new + 1.0) == 0:
            return 0.0  # Avoid division by zero

        # 4. Calculate the numerator and denominator of the load ratio (L/L0)
        numerator_L_ratio = (1.0 / r_mp_new if r_mp_new > 0 else 0) - (1.0 / (of_new + 1.0))
        denominator_L_ratio = 1.0 / self.R_MP_0 - (1.0 / (self.OF_0 + 1.0))

        if denominator_L_ratio == 0:
            return 0.0

        # This term represents (L/L0)
        load_ratio_term = (psi_new / self.PSI_0) * (numerator_L_ratio / denominator_L_ratio)

        # 5. Calculate f_i,ft using Eq. [cite_start]25 (with the corrected load term) [cite: 425]
        f_i_ft = self.F_PRIME_I_FT_0 * (r_p / r_rho_fu) * (1.0 / (of_new + 1.0)) * load_ratio_term
        return f_i_ft

    def estimate_stage_fractions(self,
                                 params: EngineParams,
                                 stage_params: StageParams) -> StageModelResult:
        """
        Estimates the stage propellant mass fraction (lambda) and inert fractions.

        This function assumes the input mission parameters are stored in
        the `stage_params` object.

        Args:
            params (EngineParams): Parameters of the *new* engine.
            stage_params (StageParams): Parameters of the *new* stage,
                                        including target dV, TWR, and payload.
                                        This model specifically uses:
                                        - `initial_twr` (as Psi)
                                        - `delta_v_ms` and `params.isp_vac_s` (to calculate R)

        Returns:
            StageModelResult: A dictionary containing the following key-value pairs:
                - 'propellant_mass_fraction': (lambda) Stage propellant mass / Stage total mass
                - 'total_inert_fraction': (f_i_total) Stage inert mass / Stage propellant mass
                - 'component_fractions': A nested dict with the 3 inert components
        """

        # Kibbey's model requires R (mass ratio) and Psi (T/W)
        # We get Psi directly from the stage parameters
        mission_T_W = stage_params.initial_twr

        # We must calculate the target mission_mass_ratio (R) from dV and Isp
        # R = exp(dV / (Isp * g0))
        if params.isp_vac_s <= 0:
            raise ValueError("Engine Isp must be > 0.")

        ve = params.isp_vac_s * vehicle_definitions.G0
        if ve <= 0:
            raise ValueError("Exhaust velocity (Isp) must be > 0.")

        mission_mass_ratio = (stage_params.delta_v_ms / ve)
        if mission_mass_ratio < 0:
            mission_mass_ratio = 0  # dV = 0 case

        # Calculate each of the three inert components
        f_i_E_S = self._calculate_fi_E_S(params, mission_mass_ratio, mission_T_W)
        f_i_ot = self._calculate_fi_ot(params)
        f_i_ft = self._calculate_fi_ft(params, mission_mass_ratio, mission_T_W)

        # Sum them to get f_i_total, as in Eq. [cite_start]11 [cite: 365]
        # (f_i_total = total inert mass / total propellant mass)
        f_i_total = f_i_E_S + f_i_ot + f_i_ft

        # Calculate lambda (propellant mass fraction), as in Eq. [cite_start]12 [cite: 366]
        # (lambda = total propellant mass / total stage mass)
        if (1.0 + f_i_total) == 0:
            lambda_fraction = 0.0
        else:
            lambda_fraction = 1.0 / (1.0 + f_i_total)

        # Return a structured dictionary with all results
        component_fractions: Dict[str, float] = {
            "f_i_E_S": f_i_E_S,  # the engine-and-structure inert mass per total propellant mass
            "f_i_ot": f_i_ot,  # the oxidizer tank inert mass per total propellant mass
            "f_i_ft": f_i_ft  # the fuel tank inert mass per total propellant mass
        }

        notes: Dict[str, Any] = {
            "reference_stage": "Atlas V (LOX/RP1)",
            "input_R": mission_mass_ratio,
            "input_Psi": mission_T_W
        }

        return StageModelResult(
            propellant_mass_fraction=lambda_fraction,
            total_inert_fraction=f_i_total,
            component_fractions=component_fractions,
            notes=notes
        )


if __name__ == "__main__":
    """
    Standalone runner for the Kibbey (2015) Stage Model.
    """

    print("[__main__] Running Kibbey Stage Model Standalone Analysis...")
    print("=" * 70)

    # --- 1. Get Engine & Stage Definitions ---
    # We will test this model with two different engines to see
    # the extrapolation effect, but using the *same* stage/mission.

    try:
        # Get the default custom rocket parameters
        engine_params, stage_params = vehicle_definitions.default_rocket_params()

        # Get a LOX/LH2 engine for comparison
        lh2_engine_params = vehicle_definitions.get_ssme_engine()

        # --- 2. Instantiate and Run Model ---
        model = KibbeyStageModel()

        # --- Test 1: Default (LOX/LCH4) Engine ---
        print("\n--- Test 1: Default (LOX/LCH4) Engine ---")
        # Use the inherited run method
        model.run_single_stage_analysis(engine_params, stage_params)

        # --- Test 2: SSME (LOX/LH2) Engine ---
        print("\n--- Test 2: SSME (LOX/LH2) Engine (Extrapolation) ---")
        # Use the same stage params, but the LH2 engine
        model.run_single_stage_analysis(lh2_engine_params, stage_params)

    except Exception as e:
        print(f"\n[__main__] ERROR during standalone analysis: {e}")

    print("=" * 70)
    print("[__main__] Standalone analysis complete.")