# models/kibbey_stage.py

from .common_params import EngineParams


class KibbeyStageModel:
    """
    Реализует "Модель массы следующего порядка" (Next-Order Mass Model) для
    оценки инертной массы ступени, основанную на масштабировании от эталонной
    ступени (Atlas V).

    Reference:
    Rho-Isp Revisited and Basic Stage Mass Estimating... (Kibbey, 2015), Section III.
    [cite: 880-977]

    Модель откалибрована по Atlas V (LOX/RP1, SC) в качестве эталона.
    """

    def __init__(self):
        """
        Инициализирует модель с эталонными константами Atlas V.
        """

        # --- Эталонные константы (Atlas V) ---

        # Контролирующие константы [cite: 890-891]
        self.F_I_TANK_0 = 0.02  # Инертная доля массы *обоих* баков [cite: 890]
        self.R_F = 2.735  # Множитель массы двигателя [cite: 891]

        # Независимые переменные (Миссия) [cite: 898]
        self.R_0 = 6.62  # Отношение масс (mi/mf) [cite: 900]
        self.PSI_0 = 1.28  # Общая T/W [cite: 901]
        self.R_MP_0 = 1.0 - (1.0 / self.R_0)  # Доля топлива (1 - 1/R) [cite: 703]

        # Независимые переменные (Топливо) [cite: 902]
        self.OF_0 = 2.7  # O/F соотношение [cite: 903]
        self.FW_ENG_0 = 78.0  # T/W двигателя [cite: 908]

        # Эталонные плотности (LOX/RP1, из Akin, ENAE 791)
        self.RHO_OX_0 = 1140.0  # LOX
        self.RHO_FU_0 = 820.0  # RP1

        # Эталонная объемная плотность
        self.RHO_BULK_0 = self._calculate_bulk_density(self.RHO_OX_0, self.RHO_FU_0, self.OF_0)

        # Производные эталонные константы баков (f'i,ot,0 и f'i,ft,0)
        # Рассчитано на основе Eq. 15-16, как в статье [cite: 931-936]
        self.F_PRIME_I_OT_0 = 0.0180  # [cite: 936]
        self.F_PRIME_I_FT_0 = 0.0255  # [cite: 936]

    def _calculate_bulk_density(self, rho_ox: float, rho_fu: float, of_ratio: float) -> float:
        """Вспомогательная функция для расчета объемной плотности."""
        total_parts = 1.0 + of_ratio
        vol_fuel = 1.0 / rho_fu
        vol_ox = of_ratio / rho_ox
        return total_parts / (vol_fuel + vol_ox)

    def _calculate_bulk_density_ratio(self, params: EngineParams) -> float:
        """
        Рассчитывает отношение объемной плотности (r_rho), также используемое как r_p.
        Реализует Eq. 18[cite: 944].
        """
        rho_bulk_new = params.bulk_density
        if self.RHO_BULK_0 == 0:
            return 1.0
        return rho_bulk_new / self.RHO_BULK_0

    def _calculate_fi_E_S(self, params: EngineParams, r_new: float, psi_new: float) -> float:
        """
        Рассчитывает инертную долю двигателя и конструкции (f_i,E&S).
        Реализует Eq. 13 и 14[cite: 921, 924].

        Args:
            params (EngineParams): Параметры *нового* двигателя.
            r_new (float): Отношение масс (R) *новой* ступени.
            psi_new (float): T/W (Psi) *новой* ступени.

        Returns:
            float: f_i,E&S
        """

        # 1. Находим T/W нового двигателя (FW_Eng) используя Eq. 13 [cite: 921]
        rho_ratio = params.bulk_density / self.RHO_BULK_0
        fw_eng_new = 39.83 * rho_ratio + 38.17

        # 2. Находим R_mp новой ступени
        r_mp_new = 1.0 - (1.0 / r_new)
        if r_mp_new == 0 or fw_eng_new == 0:
            return 0.0

        # 3. Рассчитываем f_i,E&S используя Eq. 14 [cite: 924]
        f_i_E_S = self.R_F * (psi_new / fw_eng_new) * (1.0 / r_mp_new)
        return f_i_E_S

    def _calculate_fi_ot(self, params: EngineParams) -> float:
        """
        Рассчитывает инертную долю бака окислителя (f_i,ot).
        Реализует Eq. 20[cite: 948].

        Args:
            params (EngineParams): Параметры *нового* двигателя.

        Returns:
            float: f_i,ot
        """

        # 1. Находим отношение плотности окислителя (r_rho_ox) [cite: 939]
        r_rho_ox = params.oxidizer_density / self.RHO_OX_0
        if r_rho_ox == 0:
            return 0.0

        # 2. Рассчитываем f_i,ot используя Eq. 20 [cite: 948]
        of_new = params.mixture_ratio
        f_i_ot = self.F_PRIME_I_OT_0 * (1.0 / r_rho_ox) * (of_new / (of_new + 1.0))
        return f_i_ot

    def _calculate_fi_ft(self, params: EngineParams, r_new: float, psi_new: float) -> float:
        """
        Рассчитывает инертную долю бака горючего (f_i,ft),
        которая зависит от нагрузок. Реализует Eq. 25[cite: 973].

        Args:
            params (EngineParams): Параметры *нового* двигателя.
            r_new (float): Отношение масс (R) *новой* ступени.
            psi_new (float): T/W (Psi) *новой* ступени.

        Returns:
            float: f_i,ft
        """

        # 1. Рассчитываем r_p (отношение объемных плотностей, Eq. 18)
        r_p = self._calculate_bulk_density_ratio(params)

        # 2. Рассчитываем r_rho_fu (отношение плотностей горючего, Eq. 17b) [cite: 938]
        r_rho_fu = params.fuel_density / self.RHO_FU_0
        if r_rho_fu == 0:
            return 0.0

        # 3. Находим R_mp новой ступени
        r_mp_new = 1.0 - (1.0 / r_new)
        of_new = params.mixture_ratio

        # 4. Рассчитываем числитель и знаменатель отношения нагрузок (L/L0) [cite: 961]
        numerator_L_ratio = 1.0 - r_mp_new * (1.0 / (of_new + 1.0))
        denominator_L_ratio = 1.0 - self.R_MP_0 * (1.0 / (self.OF_0 + 1.0))

        if denominator_L_ratio == 0:
            return 0.0

        load_ratio_term = (psi_new / self.PSI_0) * (numerator_L_ratio / denominator_L_ratio)

        # 5. Рассчитываем f_i,ft используя Eq. 25 [cite: 973]
        f_i_ft = self.F_PRIME_I_FT_0 * (r_p / r_rho_fu) * (1.0 / (of_new + 1.0)) * load_ratio_term
        return f_i_ft

    def estimate_propellant_mass_fraction(self, params: EngineParams,
                                          mission_mass_ratio: float,
                                          mission_T_W: float) -> dict:
        """
        Оценивает долю массы топлива ступени (lambda) и инертную массу.

        Эта функция предполагает, что входные `mission_mass_ratio` (R) и
        `mission_T_W` (Psi) являются целевыми параметрами для *новой*
        ступени, которую мы проектируем.

        Args:
            params (EngineParams): Параметры *нового* двигателя.
            mission_mass_ratio (float): Целевое отношение масс (R = mi/mf)
                                        для *новой* ступени.
            mission_T_W (float): Целевая тяговооруженность (Psi = T/W)
                                 для *новой* ступени.

        Returns:
            dict: Словарь, содержащий 'propellant_mass_fraction' (lambda)
                  и 'total_inert_fraction' (f_i_total).
        """

        # Рассчитываем каждую из трех инертных компонент
        f_i_E_S = self._calculate_fi_E_S(params, mission_mass_ratio, mission_T_W)
        f_i_ot = self._calculate_fi_ot(params)
        f_i_ft = self._calculate_fi_ft(params, mission_mass_ratio, mission_T_W)

        # Суммируем их, как в Eq. 11 [cite: 913]
        f_i_total = f_i_E_S + f_i_ot + f_i_ft

        # Рассчитываем lambda (доля массы топлива), как в Eq. 12 [cite: 914]
        if (1.0 + f_i_total) == 0:
            lambda_fraction = 0.0
        else:
            lambda_fraction = 1.0 / (1.0 + f_i_total)

        return {
            "propellant_mass_fraction": lambda_fraction,
            "total_inert_fraction": f_i_total,
            "component_fractions": {
                "f_i_E_S": f_i_E_S,
                "f_i_ot": f_i_ot,
                "f_i_ft": f_i_ft
            }
        }