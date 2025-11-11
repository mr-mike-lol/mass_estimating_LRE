# models/common_params.py

from dataclasses import dataclass
from typing import Literal, Optional

# Defines specific, allowed string values for propellant and cycle types
PropellantType = Literal["LOX/LH2", "LOX/RP1", "LOX/LCH4", "Storable"]
CycleType = Literal["GG", "SC", "EX"]  # Gas Generator, Staged Combustion, Expander


@dataclass
class EngineParams:
    """
    Stores common input parameters for various engine mass models.

    Args:
        thrust_vac_N (float): Thrust in vacuum, in Newtons.
        isp_vac_s (float): Specific impulse in vacuum, in seconds.
        chamber_pressure_Pa (float): Chamber pressure, in Pascals.
        propellant_type (PropellantType): The propellant pair.
        cycle_type (CycleType): The engine cycle type.
        mixture_ratio (float): Mixture ratio O/F (Oxidizer/Fuel).
        expansion_ratio (float): Nozzle expansion ratio (Ae/At).
        fuel_density (float): Density of the fuel, in kg/m^3.
        oxidizer_density (float): Density of the oxidizer, in kg/m^3.
        num_chambers (int): Number of chambers (for Zandbergen/Tizon models).
    """
    thrust_vac_N: float
    isp_vac_s: float
    chamber_pressure_Pa: float
    propellant_type: PropellantType
    cycle_type: CycleType
    mixture_ratio: float
    expansion_ratio: float

    # Densities in kg/m^3
    fuel_density: float
    oxidizer_density: float

    # Optional parameters for more detailed models
    num_chambers: int = 1

    @property
    def bulk_density(self) -> float:
        """
        Calculates the average (bulk) density of the propellant mixture.
        """
        total_parts = 1.0 + self.mixture_ratio
        vol_fuel = 1.0 / self.fuel_density
        vol_ox = self.mixture_ratio / self.oxidizer_density

        if (vol_fuel + vol_ox) == 0:
            return 0.0

        return total_parts / (vol_fuel + vol_ox)


@dataclass
class StageParams:
    """
    Stores parameters for stage-level models (e.g., Kibbey, Akin).

    For Akin's SSTO 1st-3rd Pass analysis, this dataclass holds
    both mission inputs (payload, dV, etc.) and is used to
    store calculated values (M_o, M_p, etc.) which are set to 0.0
    or None initially.

    Args:
        engine (EngineParams): The engine object for this stage.
        propellant_mass_kg (float): Propellant load mass for the stage, in kg.
        vehicle_gross_mass_kg (float): Gross mass of the entire launch vehicle (M_o), in kg.
        vehicle_length_m (float): Total length of the launch vehicle, in meters.
        stage_inert_mass_kg (float): Inert mass of the stage (excluding engine and payload), in kg.
        payload_mass_kg (float): Mass of the payload (and any upper stages), in kg.
        vehicle_diameter_m (float, optional): Main diameter of the vehicle/stage, in meters.
            Used for fairing and structural estimations. Defaults to 0.0.
        total_fairing_area_m2 (float, optional): Total wetted surface area of
            fairings (payload, intertank, aft). Used by `akin_mers.estimate_fairing_mass`.
            If 0, this MER is skipped. Defaults to 0.0.
    """
    # --- Core Mission & Payload Inputs ---
    payload_mass_kg: float
    delta_v_ms: float = 9200.0

    # --- 1st Pass Sizing Inputs (Akin) ---
    initial_delta: float = 0.08  # M_i / M_o guess
    initial_twr: float = 1.3
    num_engines: int = 6

    # --- Geometry Inputs (Akin) ---
    tank_geometry: Literal["Sphere", "Cylinder"] = "Sphere"
    vehicle_diameter_m: float = 0.0  # Required for 'Cylinder'
    payload_fairing_height_m: float = 7.0
    intertank_fairing_height_m: float = 7.0
    aft_fairing_height_m: float = 7.0

    # --- Calculated/Placeholder Values ---
    # These are calculated by the Akin model, not set as inputs.
    engine: Optional[EngineParams] = None  # Set during analysis
    propellant_mass_kg: float = 0.0
    vehicle_gross_mass_kg: float = 0.0
    vehicle_length_m: float = 0.0
    stage_inert_mass_kg: float = 0.0  # M_i (calculated)

    @property
    def engine_mass_kg(self) -> float:
        """
        Note: Engine mass must be calculated separately
        using one of the engine models.
        This property is here for data completeness.
        """
        # In a real application, this would be 0 and would be
        # calculated and added during analysis.
        return 0.0

    @property
    def total_stage_inert_mass_kg(self) -> float:
        """Total inert mass = Dry stage mass + Engine mass."""
        # We assume engine_mass_kg will be calculated later.
        return self.stage_inert_mass_kg + self.engine_mass_kg

@dataclass
class SolidRocketMotorParams:
    """
    Stores parameters for a Solid Rocket Motor (SRM).
    Used as input for SRM-specific MERs.

    Reference:
    [cite_start]MER for casing mass found in Akin (ENAE 791), Page 25[cite: 359, 361].

    Args:
        total_mass_kg (float): Total mass of the motor (casing + propellant).
        propellant_mass_kg (float): Mass of the propellant.
        avg_thrust_N (float): Average thrust in Newtons.
        burn_time_s (float): Burn duration in seconds.
    """
    total_mass_kg: float
    propellant_mass_kg: float
    avg_thrust_N: float
    burn_time_s: float

    @property
    def inert_mass_kg(self) -> float:
        """
        Returns the non-propellant mass (casing, nozzle, etc.).
        This is the 'inert mass' of the motor itself.
        """
        return self.total_mass_kg - self.propellant_mass_kg

    @property
    def total_impulse_Ns(self) -> float:
        """Calculates the total impulse."""
        return self.avg_thrust_N * self.burn_time_s

    @property
    def specific_impulse_s(self) -> float:
        """Calculates the effective specific impulse."""
        G0 = 9.80665  # Standard gravity
        if self.propellant_mass_kg == 0:
            return 0.0

        # Total Impulse = Isp * m_prop * g0
        return self.total_impulse_Ns / (self.propellant_mass_kg * G0)