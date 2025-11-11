# models/tizon_rema.py

import math
from .common_params import EngineParams, CycleType
from typing import Literal, Dict

# --- Глобальная константа ---
G0 = 9.80665  # Стандартное ускорение свободного падения, м/с^2

# --- 1. Данные эталонных двигателей ---
# [cite_start]Источник: Tizón & Román, 2017, Table 7 [cite: 2378]
REFERENCE_ENGINES = {
    "SC": {
        "name": "SSME",
        "mass_kg": 3177.0,
        "thrust_vac_N": 2_280_000.0,
        "pc_Pa": 2.04e7,
        "rt_m": 0.138,  # Используется для расчета (rt/rt0)
        "m_dot_kg_s": 512.6,
        "of_ratio": 6.0,
        "expansion_ratio": 77.5,
        "prop_type": "LOX/LH2",
        "fuel_density": 71.0,
        "oxidizer_density": 1140.0
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
        "fuel_density": 71.0,
        "oxidizer_density": 1140.0
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
        "fuel_density": 71.0,
        "oxidizer_density": 1140.0
    }
}
# Добавим объемную плотность и m_dot (если не задан) в эталонные данные
for engine in REFERENCE_ENGINES.values():
    engine['bulk_density'] = (1.0 + engine['of_ratio']) / \
                             ((1.0 / engine['fuel_density']) + (engine['of_ratio'] / engine['oxidizer_density']))
    if 'm_dot_kg_s' not in engine:
        engine['m_dot_kg_s'] = engine['thrust_vac_N'] / (engine['isp_vac_s'] * G0)

# --- 2. Коэффициенты долей массы (Alpha_i) ---
# [cite_start]Источник: Tizón & Román, 2017, Table 5 (Design) [cite: 2355]
ALPHAS_DESIGN = {
    "SC": {
        "tubes": 0.0640, "manifold": 0.1717, "jacket": 0.0470,
        "combustion_chamber": 0.1137, "gas_generator": 0.0551,
        "ox_turbopump": 0.1089, "fuel_turbopump": 0.1407,
        "ox_valve": 0.0698, "fuel_valve": 0.0690
        # "radiative_nozzle" не используется, т.к. SSME регенеративный
    },
    "GG": {
        "tubes": 0.0987, "manifold": 0.0146, "jacket": 0.0589,
        "combustion_chamber": 0.1934, "gas_generator": 0.0438,
        "ox_turbopump": 0.1254, "fuel_turbopump": 0.1354,
        "ox_valve": 0.0352, "fuel_valve": 0.0352
    },
    "EX": {
        "tubes": 0.0629, "manifold": 0.2301, "jacket": 0.0210,
        "combustion_chamber": 0.2007, "gas_generator": 0.0,  # У экспандера нет ГГ
        "ox_turbopump": 0.1411, "fuel_turbopump": 0.1411,
        "ox_valve": 0.0202, "fuel_valve": 0.0185
    }
}

# [cite_start]Источник: Tizón & Román, 2017, Table 6 (Historical) [cite: 2357]
ALPHAS_HISTORICAL = {
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

# --- 3. Степенные показатели (a_ij) ---
# [cite_start]Источник: Tizón & Román, 2017, Table 2 (Design) [cite: 2340]
# Ключи: pc, expansion_ratio, rt, m_dot_prop, prop_density
EXPONENTS_DESIGN = {
    "tubes": {"pc": 1.0, "expansion_ratio": 1.0, "rt": 2.0, "m_dot_prop": 0.0, "prop_density": 0.0},
    "manifold": {"pc": 1.0, "expansion_ratio": 0.0, "rt": 1.0, "m_dot_prop": 1.0, "prop_density": -1.0},
    "jacket": {"pc": 1.0, "expansion_ratio": 1.5, "rt": 3.0, "m_dot_prop": 0.0, "prop_density": 0.0},
    "combustion_chamber": {"pc": 1.0, "expansion_ratio": 0.0, "rt": 2.0, "m_dot_prop": 0.0, "prop_density": 0.0},
    "gas_generator": {"pc": 1.0, "expansion_ratio": 0.0, "rt": 2.0, "m_dot_prop": 0.0, "prop_density": 0.0},
    "ox_turbopump": {"pc": 0.15, "expansion_ratio": 0.0, "rt": 0.0, "m_dot_prop": 0.9, "prop_density": -0.45},
    "fuel_turbopump": {"pc": 0.15, "expansion_ratio": 0.0, "rt": 0.0, "m_dot_prop": 0.9, "prop_density": -0.45},
    "ox_valve": {"pc": 1.0, "expansion_ratio": 0.0, "rt": 0.0, "m_dot_prop": 1.0, "prop_density": -1.0},
    "fuel_valve": {"pc": 1.0, "expansion_ratio": 0.0, "rt": 0.0, "m_dot_prop": 1.0, "prop_density": -1.0},
}

# [cite_start]Источник: Tizón & Román, 2017, Table 3 (Historical) [cite: 2341]
EXPONENTS_HISTORICAL = {
    "combustion_chamber": {"pc": 1.0, "expansion_ratio": 0.0, "rt": 1.4, "m_dot_prop": 0.0, "prop_density": 0.0},
    "ox_turbopump": {"pc": 0.53, "expansion_ratio": 0.0, "rt": 0.0, "m_dot_prop": 0.53, "prop_density": -0.53},
    "fuel_turbopump": {"pc": 0.53, "expansion_ratio": 0.0, "rt": 0.0, "m_dot_prop": 0.53, "prop_density": -0.53},
    "ox_valve": {"pc": 0.3, "expansion_ratio": 0.0, "rt": 0.0, "m_dot_prop": 0.625, "prop_density": -0.625},
    "fuel_valve": {"pc": 0.3, "expansion_ratio": 0.0, "rt": 0.0, "m_dot_prop": 0.625, "prop_density": -0.625},
    "structure": {"pc": 0.92, "expansion_ratio": 0.0, "rt": 1.84, "m_dot_prop": 0.0, "prop_density": 0.0},
    "auxiliary": {"pc": 0.0, "expansion_ratio": 0.0, "rt": 1.0, "m_dot_prop": 0.0, "prop_density": 0.0},
}


class TizonEngineModel:
    """
    Реализует безразмерную модель массы из статьи Tizón & Román, 2017.

    Эта модель оценивает массу двигателя путем масштабирования массы
    каждого компонента относительно эталонного двигателя (SSME, LE-5,
    или RL10A) на основе набора безразмерных соотношений параметров.

    Reference:
    A Mass Model for Liquid Propellant Rocket Engines (Tizón & Román, 2017).
    [cite_start][cite: 1737-2648]
    """

    def __init__(self, cycle_type: CycleType, method: Literal["design", "historical"] = "historical"):
        """
        Инициализирует модель с заданным типом цикла и методом расчета.

        Args:
            cycle_type (CycleType): 'GG', 'SC', или 'EX'.
            method (Literal["design", "historical"]):
                - 'design': использует показатели из Таблицы 2.
                - 'historical': использует показатели из Таблицы 3.
        """
        if cycle_type not in REFERENCE_ENGINES:
            raise ValueError(f"Cycle type {cycle_type} not supported by Tizon model.")

        self.cycle_type = cycle_type
        self.method = method

        # 1. Выбор эталонного двигателя
        self.ref_engine = REFERENCE_ENGINES[cycle_type]
        self.ref_mass_kg = self.ref_engine["mass_kg"]

        # 2. Выбор коэффициентов Alpha
        if method == "design":
            self.alphas = ALPHAS_DESIGN[cycle_type]
            self.exponents = EXPONENTS_DESIGN
        else:  # historical
            self.alphas = ALPHAS_HISTORICAL[cycle_type]
            self.exponents = EXPONENTS_HISTORICAL

    def _calculate_param_ratios(self, params: EngineParams) -> Dict[str, float]:
        """
        Рассчитывает соотношения (P_j / P_j_0) для всех
        параметров, используемых в модели.
        """

        # Рассчитываем производные параметры для 'params'
        m_dot_new = params.thrust_vac_N / (params.isp_vac_s * G0)
        bulk_density_new = params.bulk_density

        # Прокси-соотношение для радиуса горловины (rt)
        # (rt/rt0)^2 = (At/At0) ~ (F/pc) / (F0/pc0)
        rt_ratio_sq_proxy = (params.thrust_vac_N / params.chamber_pressure_Pa) / \
                            (self.ref_engine["thrust_vac_N"] / self.ref_engine["pc_Pa"])

        # Соотношения для ТНА и клапанов зависят от m_dot и плотности
        # каждого компонента (окислителя или горючего)
        of_new = params.mixture_ratio
        of_ref = self.ref_engine['of_ratio']

        m_dot_ox_new = m_dot_new * (of_new / (1.0 + of_new))
        m_dot_fuel_new = m_dot_new * (1.0 / (1.0 + of_new))

        m_dot_ox_ref = self.ref_engine['m_dot_kg_s'] * (of_ref / (1.0 + of_ref))
        m_dot_fuel_ref = self.ref_engine['m_dot_kg_s'] * (1.0 / (1.0 + of_ref))

        return {
            "pc": params.chamber_pressure_Pa / self.ref_engine["pc_Pa"],
            "expansion_ratio": params.expansion_ratio / self.ref_engine["expansion_ratio"],
            "rt": rt_ratio_sq_proxy ** 0.5,

            # Для "design" модели с раздельными ТНА/клапанами
            "m_dot_ox": m_dot_ox_new / m_dot_ox_ref,
            "m_dot_fuel": m_dot_fuel_new / m_dot_fuel_ref,
            "density_ox": params.oxidizer_density / self.ref_engine["oxidizer_density"],
            "density_fuel": params.fuel_density / self.ref_engine["fuel_density"],

            # Для "historical" модели с общими ТНА/клапанами
            "m_dot_prop": m_dot_new / self.ref_engine["m_dot_kg_s"],
            "prop_density": bulk_density_new / self.ref_engine["bulk_density"],
        }

    def _calculate_component_ratio(self, component: str, param_ratios: Dict[str, float]) -> float:
        """
        Рассчитывает (m_i / m_i_0) для одного компонента, используя Eq. (7).

        Args:
            component (str): Название компонента (например, 'combustion_chamber').
            param_ratios (Dict[str, float]): Предварительно рассчитанные соотношения P_j/P_j_0.

        Returns:
            float: (m_i / m_i_0)
        """

        if component not in self.exponents:
            # Компонент присутствует в Alphas, но не в Exponents
            # (например, 'radiative_nozzle' в 'design')
            return 1.0

        exps = self.exponents[component]

        # Обработка особых случаев для ТНА и клапанов в "design"
        if self.method == "design":
            if component == "ox_turbopump" or component == "ox_valve":
                # Использовать m_dot_ox и density_ox
                m_dot_ratio = param_ratios["m_dot_ox"]
                density_ratio = param_ratios["density_ox"]
            elif component == "fuel_turbopump" or component == "fuel_valve":
                # Использовать m_dot_fuel и density_fuel
                m_dot_ratio = param_ratios["m_dot_fuel"]
                density_ratio = param_ratios["density_fuel"]
            else:
                # Стандартный m_dot_prop и prop_density
                m_dot_ratio = param_ratios["m_dot_prop"]
                density_ratio = param_ratios["prop_density"]
        else:  # historical
            m_dot_ratio = param_ratios["m_dot_prop"]
            density_ratio = param_ratios["prop_density"]

        # Рассчитываем произведение, Eq. (7)
        # (P1/P1_0)^a1 * (P2/P2_0)^a2 * ...
        ratio = 1.0
        ratio *= (param_ratios["pc"] ** exps.get("pc", 0.0))
        ratio *= (param_ratios["expansion_ratio"] ** exps.get("expansion_ratio", 0.0))
        ratio *= (param_ratios["rt"] ** exps.get("rt", 0.0))
        ratio *= (m_dot_ratio ** exps.get("m_dot_prop", 0.0))
        ratio *= (density_ratio ** exps.get("prop_density", 0.0))

        return ratio

    def estimate_total_mass(self, params: EngineParams) -> dict:
        """
        Оценивает общую массу двигателя и ее разбивку по компонентам.

        Реализует Eq. (8): m_engine = m_engine_0 * SUM(alpha_i * (m_i / m_i_0))

        Args:
            params (EngineParams): Параметры нового двигателя для оценки.

        Returns:
            dict: Словарь, содержащий 'total_mass_kg' и 'components_kg'.
        """

        # 1. Рассчитать все соотношения P_j / P_j_0 один раз
        param_ratios = self._calculate_param_ratios(params)

        total_ratio_sum = 0.0
        component_masses = {}

        # 2. Итерация по всем компонентам из таблицы Alphas
        for component, alpha_i in self.alphas.items():
            if alpha_i == 0.0:
                component_masses[component] = 0.0
                continue

            # 3. Рассчитать m_i / m_i_0 для этого компонента
            m_i_ratio = self._calculate_component_ratio(component, param_ratios)

            # 4. Добавить в общую сумму (Eq. 8)
            total_ratio_sum += alpha_i * m_i_ratio

            # 5. Сохранить массу компонента для детализации
            # m_i = m_engine_0 * alpha_i * (m_i / m_i_0)
            component_masses[component] = self.ref_mass_kg * alpha_i * m_i_ratio

        # 6. Рассчитать итоговую массу
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