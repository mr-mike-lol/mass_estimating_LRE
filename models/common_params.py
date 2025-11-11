# models/common_params.py

from dataclasses import dataclass
from typing import Literal

# Defines specific, allowed string values for propellant and cycle types
# Это улучшает проверку типов и автодополнение в IDE
PropellantType = Literal["LOX/LH2", "LOX/RP1", "LOX/LCH4", "Storable"]
CycleType = Literal["GG", "SC", "EX"]  # Gas Generator, Staged Combustion, Expander


@dataclass
class EngineParams:
    """
    Хранит общие входные параметры для различных моделей массы двигателя.

    Args:
        thrust_vac_N (float): Тяга в вакууме, в Ньютонах.
        isp_vac_s (float): Удельный импульс в вакууме, в секундах.
        chamber_pressure_Pa (float): Давление в камере сгорания, в Паскалях.
        propellant_type (PropellantType): Тип топливной пары.
        cycle_type (CycleType): Тип цикла двигателя.
        mixture_ratio (float): Соотношение компонентов O/F (Окислитель/Горючее).
        expansion_ratio (float): Степень расширения сопла (Ae/At).
        fuel_density (float): Плотность горючего, в кг/м^3.
        oxidizer_density (float): Плотность окислителя, в кг/м^3.
        num_chambers (int): Количество камер сгорания (для моделей Zandbergen/Tizon).
    """
    thrust_vac_N: float
    isp_vac_s: float
    chamber_pressure_Pa: float
    propellant_type: PropellantType
    cycle_type: CycleType
    mixture_ratio: float
    expansion_ratio: float

    # Плотности в кг/м^3
    fuel_density: float
    oxidizer_density: float

    # Опциональные параметры для более детальных моделей
    num_chambers: int = 1

    @property
    def bulk_density(self) -> float:
        """
        Рассчитывает среднюю (объемную) плотность компонентов топлива.
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
    Хранит параметры для моделей уровня ступени (например, Kibbey).

    Args:
        engine (EngineParams): Объект, описывающий двигатель этой ступени.
        propellant_mass_kg (float): Загрузка массы топлива в ступени, в кг.
        vehicle_gross_mass_kg (float): Полная стартовая масса ракетоносителя (M_o), в кг.
        vehicle_length_m (float): Общая длина ракетоносителя, в метрах.
        stage_inert_mass_kg (float): Инертная масса ступени (без двигателя и топлива), в кг.
        payload_mass_kg (float): Масса полезной нагрузки (и верхних ступеней), в кг.
    """
    engine: EngineParams
    propellant_mass_kg: float
    vehicle_gross_mass_kg: float
    vehicle_length_m: float
    stage_inert_mass_kg: float  # Инертная масса без двигателя
    payload_mass_kg: float  # Все, что над этой ступенью

    @property
    def engine_mass_kg(self) -> float:
        """
        Примечание: Масса двигателя должна рассчитываться отдельно
        с использованием одной из моделей.
        Это свойство здесь для полноты данных.
        """
        # В реальном приложении здесь будет 0, и она будет
        # рассчитана и добавлена во время анализа.
        return 0.0

    @property
    def total_stage_inert_mass_kg(self) -> float:
        """Общая инертная масса = Сухая масса ступени + Масса двигателя."""
        # Мы предполагаем, что engine_mass_kg будет рассчитана позже.
        return self.stage_inert_mass_kg + self.engine_mass_kg