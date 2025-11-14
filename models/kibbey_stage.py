# models/kibbey_stage.py
"""

"""

from .common_params import EngineParams
from vehicle_definitions import DENSITY_RP1, DENSITY_LH2, DENSITY_LOX


class KibbeyStageModel:
    """
    Implements the 'Next-Order Mass Model' for estimating stage inert mass,
    based on scaling from a reference stage (Atlas V).

    Reference:
    Rho-Isp Revisited and Basic Stage Mass Estimating... (Kibbey, 2015), Section III.


    The model is calibrated to the Atlas V (LOX/RP1, SC) as the reference.
    """

    def __init__(self):
        """
        Initializes the model with reference constants from the Atlas V.
        """

        # --- Reference Constants (Atlas V) ---

        # Controlling constants
        self.F_I_TANK_0 = 0.02  # Inert fraction of *both* tanks
        self.R_F = 2.735  # Engine mass multiplier

        # Independent variables (Mission)
        self.R_0 = 6.62  # Mass ratio (mi/mf)
        self.PSI_0 = 1.28  # Overall T/W
        self.R_MP_0 = 1.0 - (1.0 / self.R_0)  # Propellant fraction (1 - 1/R)

        # Independent variables (Propellant)
        self.OF_0 = 2.7  # O/F ratio
        self.FW_ENG_0 = 78.0  # Engine T/W

        # These are standard values used for the LOX/RP1 reference.
        self.RHO_OX_0 = DENSITY_LOX
        self.RHO_FU_0 = DENSITY_RP1

        # Reference bulk density
        self.RHO_BULK_0 = self._calculate_bulk_density(self.RHO_OX_0, self.RHO_FU_0, self.OF_0)

        # Derived reference tank constants (f'i,ot,0 and f'i,ft,0)
        # Calculated based on Eq. 15-16, as in the paper
        self.F_PRIME_I_OT_0 = 0.0180  #
        self.F_PRIME_I_FT_0 = 0.0255  #

    def _calculate_bulk_density(self, rho_ox: float, rho_fu: float, of_ratio: float) -> float:
        """Helper function to calculate bulk density."""
        total_parts = 1.0 + of_ratio
        vol_fuel = 1.0 / rho_fu
        vol_ox = of_ratio / rho_ox
        if (vol_fuel + vol_ox) == 0:
            return 0.0
        return total_parts / (vol_fuel + vol_ox)

    def _calculate_bulk_density_ratio(self, params: EngineParams) -> float:
        """
        Calculates the bulk density ratio (r_rho), also used as r_p.
        Implements Eq. 18 .
        """
        rho_bulk_new = params.bulk_density
        if self.RHO_BULK_0 == 0:
            return 1.0
        return rho_bulk_new / self.RHO_BULK_0

    def _calculate_fi_E_S(self, params: EngineParams, r_new: float, psi_new: float) -> float:
        """
        Calculates the engine and structure inert fraction (f_i,E&S).
        Implements Eq. 13 and 14 .

        Args:
            params (EngineParams): Parameters of the *new* engine.
            r_new (float): Mass ratio (R) of the *new* stage.
            psi_new (float): T/W (Psi) of the *new* stage.

        Returns:
            float: f_i,E&S
        """

        # 1. Find the new engine T/W (FW_Eng) using Eq. 13
        rho_ratio = params.bulk_density / self.RHO_BULK_0
        fw_eng_new = 39.83 * rho_ratio + 38.17

        # 2. Find the R_mp of the new stage
        r_mp_new = 1.0 - (1.0 / r_new)
        if r_mp_new == 0 or fw_eng_new == 0:
            return 0.0

        # 3. Calculate f_i,E&S using Eq. 14
        f_i_E_S = self.R_F * (psi_new / fw_eng_new) * (1.0 / r_mp_new)
        return f_i_E_S

    def _calculate_fi_ot(self, params: EngineParams) -> float:
        """
        Calculates the oxidizer tank inert fraction (f_i,ot).
        Implements Eq. 20 .

        Args:
            params (EngineParams): Parameters of the *new* engine.

        Returns:
            float: f_i,ot
        """

        # 1. Find the oxidizer density ratio (r_rho_ox)
        r_rho_ox = params.oxidizer_density / self.RHO_OX_0
        if r_rho_ox == 0:
            return 0.0

        # 2. Calculate f_i,ot using Eq. 20
        of_new = params.mixture_ratio
        f_i_ot = self.F_PRIME_I_OT_0 * (1.0 / r_rho_ox) * (of_new / (of_new + 1.0))
        return f_i_ot

    def _calculate_fi_ft(self, params: EngineParams, r_new: float, psi_new: float) -> float:
        """
        Calculates the fuel tank inert fraction (f_i,ft),
        which depends on loads. Implements Eq. 25 .

        Args:
            params (EngineParams): Parameters of the *new* engine.
            r_new (float): Mass ratio (R) of the *new* stage.
            psi_new (float): T/W (Psi) of the *new* stage.

        Returns:
            float: f_i,ft
        """

        # 1. Calculate r_p (bulk density ratio, Eq. 18)
        r_p = self._calculate_bulk_density_ratio(params)

        # 2. Calculate r_rho_fu (fuel density ratio, Eq. 17a)
        r_rho_fu = params.fuel_density / self.RHO_FU_0
        if r_rho_fu == 0:
            return 0.0

        # 3. Find R_mp of the new stage
        r_mp_new = 1.0 - (1.0 / r_new)
        of_new = params.mixture_ratio

        # 4. Calculate the numerator and denominator of the load ratio (L/L0)
        numerator_L_ratio = 1.0 / r_mp_new - (1.0 / (of_new + 1.0))
        denominator_L_ratio = 1.0 / self.R_MP_0 - (1.0 / (self.OF_0 + 1.0))

        if denominator_L_ratio == 0:
            return 0.0

        # This term represents (L/L0) from Eq. 22, assuming F/F0 scales with (Psi/Psi_0)
        # per Eq. 24
        load_ratio_term = (psi_new / self.PSI_0) * (numerator_L_ratio / denominator_L_ratio)

        # 5. Calculate f_i,ft using Eq. 25
        f_i_ft = self.F_PRIME_I_FT_0 * (r_p / r_rho_fu) * (1.0 / (of_new + 1.0)) * load_ratio_term
        return f_i_ft

    def estimate_propellant_mass_fraction(self, params: EngineParams,
                                          mission_mass_ratio: float,
                                          mission_T_W: float) -> dict:
        """
        Estimates the stage propellant mass fraction (lambda) and inert mass.

        This function assumes that the input `mission_mass_ratio` (R) and
        `mission_T_W` (Psi) are the target parameters for the *new*
        stage being designed.

        Args:
            params (EngineParams): Parameters of the *new* engine.
            mission_mass_ratio (float): Target mass ratio (R = mi/mf)
                                        for the *new* stage.
            mission_T_W (float): Target thrust-to-weight (Psi = T/W)
                                 for the *new* stage.

        Returns:
            dict: A dictionary containing 'propellant_mass_fraction' (lambda)
                  and 'total_inert_fraction' (f_i_total).
        """

        # Calculate each of the three inert components
        f_i_E_S = self._calculate_fi_E_S(params, mission_mass_ratio, mission_T_W)
        f_i_ot = self._calculate_fi_ot(params)
        f_i_ft = self._calculate_fi_ft(params, mission_mass_ratio, mission_T_W)

        # Sum them, as in Eq. 11
        f_i_total = f_i_E_S + f_i_ot + f_i_ft

        # Calculate lambda (propellant mass fraction), as in Eq. 12
        if (1.0 + f_i_total) == 0:
            lambda_fraction = 0.0
        else:
            lambda_fraction = 1.0 / (1.0 + f_i_total)

        return {
            "propellant_mass_fraction": lambda_fraction,
            "total_inert_fraction: the stage propellant mass fraction": f_i_total,
            "component_fractions": {
                "f_i_E_S: the engine-and-structure inert mass per total propellant mass": f_i_E_S,
                "f_i_ot: the oxidizer tank inert mass per total propellant mass": f_i_ot,
                "f_i_ft: the fuel tank inert mass per total propellant mass": f_i_ft
            }
        }