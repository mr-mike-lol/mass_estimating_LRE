# main.py

from models.common_params import EngineParams
# Мы предполагаем, что у нас есть эти модули
from models import akin_mers, mota_schlingloff
from models import zandbergen_engine, tizon_rema, kibbey_stage
from vehicle_definitions import get_le5_engine, get_ssme_engine, get_rl10a_engine, get_rd120_engine


def run_engine_mass_comparison(engine_name: str, engine_params: EngineParams):
    """
    Запускает все применимые модели ОЦЕНКИ МАССЫ ДВИГАТЕЛЯ
    для заданной конфигурации двигателя и выводит результаты,
    включая покомпонентную разбивку.
    """
    print("=" * 70)
    print(f"🚀 ЗАПУСК СРАВНЕНИЯ ДЛЯ: {engine_name}")
    print(f"   (Параметры: {engine_params.thrust_vac_N / 1_000_000:.3f} МН, "
          f"P_c: {engine_params.chamber_pressure_Pa / 1e6:.2f} МПа, "
          f"Ae/At: {engine_params.expansion_ratio:.1f})")
    print("=" * 70)

    # Словарь для хранения полных результатов
    results = {}

    # --- 1. Модель Akin (1991) ---
    try:
        results["Akin (1991)"] = akin_mers.estimate_engine_mass(
            engine_params.thrust_vac_N,
            engine_params.expansion_ratio
        )
        results["Akin (1991)"]["note"] = "Примечание: Упрощенная MER, не учитывает тип цикла или топлива."
    except Exception as e:
        results["Akin (1991)"] = {"error": str(e)}

    # --- 2. Модель Zandbergen (2015) ---
    try:
        results["Zandbergen (2015)"] = zandbergen_engine.estimate_engine_mass(engine_params)
        results["Zandbergen (2015)"]["note"] = "Примечание: Простая регрессия (Тяга + класс топлива). Монолитная."
    except Exception as e:
        results["Zandbergen (2015)"] = {"error": str(e)}

    # --- 3. Модель Tizón/RemA (2017) ---
    try:
        tizon_model = tizon_rema.TizonEngineModel(
            engine_params.cycle_type,
            method="historical"
        )
        tizon_result = tizon_model.estimate_total_mass(engine_params)
        ref_name = tizon_model.ref_engine['name']
        tizon_result["note"] = f"Примечание: *Относительная* модель. Эталон: ({ref_name})."
        results["Tizón/RemA (2017)"] = tizon_result
    except ValueError as e:
        results["Tizón/RemA (2017)"] = {"error": f"Неприменимо ({e})"}
    except Exception as e:
        results["Tizón/RemA (2017)"] = {"error": str(e)}

    # --- 4. Модель Mota/Schlingloff (2005) ---
    try:
        has_boost_pumps = (engine_params.cycle_type == "SC")
        mota_result = mota_schlingloff.estimate_total_engine_mass(
            engine_params,
            has_boost_pumps=has_boost_pumps
        )
        mota_result["note"] = "Примечание: Аналитико-статистическая модель[cite: 2876], учитывает Pc, O/F и др."
        results["Mota/Schlingloff (2005)"] = mota_result
    except Exception as e:
        results["Mota/Schlingloff (2005)"] = {"error": str(e)}

    # --- Вывод результатов ---
    print("--- 📊 Результаты оценки массы двигателя ---")

    # Сначала выводим сводную таблицу общих масс
    print("\n--- Сводка по общей массе ---")
    for model_name, res in results.items():
        if "error" in res:
            mass_str = f"ОШИБКА: {res['error']}"
            note = ""
        else:
            mass_str = f"{res['total_mass_kg']:,.1f} kg"
            note = res.get('note', '')

        print(f"{model_name:<28} | {mass_str:<20} | {note}")

    # Затем выводим покомпонентную разбивку
    print("\n--- Покомпонентная разбивка (где доступно) ---")
    for model_name, res in results.items():
        if "components_kg" in res:
            print(f"  --- {model_name} ---")
            total = res['total_mass_kg']
            for component, mass in res['components_kg'].items():
                percent = (mass / total * 100) if total > 0 else 0
                print(f"    {component:<30}: {mass:10,.1f} kg ({percent:5.1f}%)")
            print(f"    {'ИТОГО':<30}: {total:10,.1f} kg (100.0%)")
            print()  # Пустая строка для разделения

    print("\n")


def run_stage_model_example(engine_name: str, engine_params: EngineParams):
    """
    Запускает пример модели ОЦЕНКИ МАССЫ СТУПЕНИ (Kibbey).
    """
    print("=" * 70)
    print(f"🛰️  ЗАПУСК ПРИМЕРА МОДЕЛИ СТУПЕНИ (Kibbey, 2015)")
    print(f"   На базе двигателя: {engine_name}")
    print(f"   Эталонная ступень: Atlas V (LOX/RP1)")
    print("=" * 70)

    # Модель Kibbey требует миссионных параметров R (отношение масс) и Psi (T/W)
    # Мы используем эталонные значения Atlas V для демонстрации

    # Возьмем целевые параметры Atlas V как пример [cite: 899-901]
    example_mission_R = 6.62
    example_mission_Psi = 1.28

    try:
        kibbey_model = kibbey_stage.KibbeyStageModel()

        stage_fractions = kibbey_model.estimate_propellant_mass_fraction(
            engine_params,
            mission_mass_ratio=example_mission_R,
            mission_T_W=example_mission_Psi
        )

        lambda_frac = stage_fractions['propellant_mass_fraction']
        inert_frac = 1.0 - lambda_frac

        print(f"--- 📈 Результаты модели ступени Kibbey ---")

        # Добавляем примечание о релевантности
        if engine_params.propellant_type == "LOX/RP1":
            print("(!) ПРИМЕЧАНИЕ: Это наиболее релевантный тест, т.к. и двигатель (LOX/RP1),")
            print("    и эталон (Atlas V, LOX/RP1) используют одинаковый класс топлива.")
        else:
            print(f"(!) ПРИМЕЧАНИЕ: Это экстраполяция. Модель масштабирует ступень LOX/RP1")
            print(f"    для работы с параметрами двигателя {engine_params.propellant_type}.")

        print(f"\nДля ступени с R={example_mission_R} и T/W={example_mission_Psi}:")
        print(f"  Расчетная доля топлива (lambda): {lambda_frac:.4f}")
        print(f"  Расчетная инертная доля (1-lambda): {inert_frac:.4f}")
        print(f"  (f_i_total, инертная доля от массы топлива): {stage_fractions['total_inert_fraction']:.4f}")

        print("\n  Разбивка инертной доли (f_i):")
        total_f_i = stage_fractions['total_inert_fraction']
        for comp, val in stage_fractions['component_fractions'].items():
            percent = (val / total_f_i * 100) if total_f_i > 0 else 0
            print(f"    {comp:<10}: {val:8.5f} ({percent:5.1f}%)")
        print(f"    {'ИТОГО':<10}: {total_f_i:8.5f} (100.0%)")


    except Exception as e:
        print(f"ОШИБКА при выполнении модели Kibbey: {e}")
    print("\n")


if __name__ == "__main__":
    # --- 1. Получаем наши тестовые двигатели ---
    engine_le5 = get_le5_engine()  # LOX/LH2 GG
    engine_ssme = get_ssme_engine()  # LOX/LH2 SC
    engine_rl10 = get_rl10a_engine()  # LOX/LH2 EX
    engine_rd120 = get_rd120_engine()  # LOX/RP1 SC

    # --- 2. Запускаем сравнение МАССЫ ДВИГАТЕЛЯ ---
    print("||||" + "=" * 70 + "||||")
    print("||||   ЧАСТЬ 1: СРАВНЕНИЕ МОДЕЛЕЙ МАССЫ ДВИГАТЕЛЯ")
    print("||||" + "=" * 70 + "||||")

    run_engine_mass_comparison("LE-5 (LOX/LH2 GG)", engine_le5)
    run_engine_mass_comparison("SSME (LOX/LH2 SC)", engine_ssme)
    run_engine_mass_comparison("RL10A (LOX/LH2 EX)", engine_rl10)
    run_engine_mass_comparison("RD-120 (LOX/RP1 SC)", engine_rd120)

    # --- 3. Запускаем пример МОДЕЛИ МАССЫ СТУПЕНИ ---
    print("||||" + "=" * 70 + "||||")
    print("||||   ЧАСТЬ 2: ПРИМЕР МОДЕЛИ МАССЫ СТУПЕНИ (Kibbey)")
    print("||||" + "=" * 70 + "||||")

    # Модель Kibbey откалибрована по Atlas V (LOX/RP1) .
    # Запуск ее с двигателем RD-120 (также LOX/RP1)
    # является наиболее релевантным примером.
    run_stage_model_example("RD-120 (LOX/RP1 SC)", engine_rd120)

    # Также запустим для SSME, чтобы увидеть экстраполяцию
    run_stage_model_example("SSME (LOX/LH2 SC)", engine_ssme)