# isomer_analyzer.py — Модуль оценки стабильности ядерных изомеров
# Основан на Правиле 1 языка Crys: цветовая нейтральность
# Версия 1.0, Валентин Лебедкин, 2026

import math

# Фундаментальные константы модели СТПК
SIGMA = 0.75          # натяжение глюонной струны, ГэВ/фм
L_QCD = 1.0           # характерный масштаб конфайнмента, фм
K_B = 1.0             # тепловая энергия ядра, МэВ (kT ≈ 1 МэВ)
MEV_TO_GEV = 0.001    # перевод МэВ в ГэВ

def estimate_string_break_length(A, Z, spin_initial, parity_initial, spin_final, parity_final):
    """
    Оценка суммарной длины разрываемых глюонных трубок L_разрыва (в фм)
    при гипотетическом переходе изомера.
    Параметры:
        A – массовое число
        Z – заряд
        spin_initial, parity_initial – спин и чётность изомера
        spin_final, parity_final – спин и чётность основного состояния
    Возвращает:
        L_break – длина разрыва в фм
    """
    # Базовое число нуклонов, участвующих в переходе
    N = A - Z
    # Число нуклонов, которые должны перестроиться:
    # - изменение спина на ΔI требует переворота спина у ~ΔI нуклонов
    delta_I = abs(spin_initial - spin_final)
    # - изменение чётности (если есть) требует дополнительного перестроения
    parity_change = 1 if parity_initial != parity_final else 0
    
    # Каждый затронутый нуклон рвёт трубки с соседями (в среднем ~6 соседей)
    nucleons_affected = max(1, int(delta_I / 2) + parity_change * 3)
    # Средняя длина трубки между соседними нуклонами ~1 фм
    L_break = nucleons_affected * 6 * 1.0  # фм
    return L_break

def estimate_break_energy(L_break):
    """
    Энергия нарушения цветовой нейтральности (ГэВ).
    E_наруш = σ * L_разрыва
    """
    return SIGMA * L_break

def estimate_half_life(E_break, delta_E):
    """
    Оценка периода полураспада изомера (в секундах).
    Использует правило Аррениуса с подавлением exp(-E_break / kT_ядра).
    Параметры:
        E_break – энергия нарушения цветовой нейтральности (ГэВ)
        delta_E – энергия перехода (разность масс, ГэВ)
    Возвращает:
        период полураспада в секундах
    """
    # Если энергия нарушения больше доступной – переход сильно подавлен
    if E_break > delta_E:
        # Фундаментальная частота попыток: частота нуклонных колебаний ~10^22 Гц
        attempt_frequency = 1e22  # 1/с
        # Вероятность туннелирования через барьер
        tunneling_prob = math.exp(-E_break / (K_B * MEV_TO_GEV))
        # Среднее время жизни
        lifetime = 1.0 / (attempt_frequency * tunneling_prob)
        return lifetime
    else:
        # Если энергия нарушения меньше доступной – переход разрешён, время мало
        return 1e-12  # пикосекунды

def analyze_isomer(A, Z, spin_i, parity_i, spin_f, parity_f, delta_E):
    """
    Главная функция анализа изомера.
    Возвращает словарь с результатами.
    """
    L_break = estimate_string_break_length(A, Z, spin_i, parity_i, spin_f, parity_f)
    E_break = estimate_break_energy(L_break)
    half_life = estimate_half_life(E_break, delta_E)
    
    # Качественная оценка
    if half_life > 1e-6:  # микросекунды и более
        stability_class = "Изомер"
    elif half_life > 1e-9:  # наносекунды
        stability_class = "Метастабильное состояние"
    else:
        stability_class = "Короткоживущее возбуждение"
    
    return {
        "L_break_fm": L_break,
        "E_break_GeV": E_break,
        "E_break_MeV": E_break / MEV_TO_GEV,
        "delta_E_MeV": delta_E / MEV_TO_GEV,
        "half_life_s": half_life,
        "stability_class": stability_class
    }

# ========== ТЕСТОВЫЙ ПРИМЕР: Tc-99m ==========
if __name__ == "__main__":
    # Tc-99m: изомер со спином 1/2- -> основное состояние 9/2+
    # Энергия перехода 142 кэВ = 0.000142 ГэВ
    result = analyze_isomer(
        A=99, Z=43,
        spin_i=0.5, parity_i='-',   # изомер 1/2-
        spin_f=4.5, parity_f='+',   # основное 9/2+
        delta_E=0.000142             # ГэВ
    )
    
    print("=== Тест Tc-99m ===")
    print(f"Длина разрыва: {result['L_break_fm']:.1f} фм")
    print(f"Энергия нарушения: {result['E_break_MeV']:.1f} МэВ")
    print(f"Энергия перехода: {result['delta_E_MeV']:.3f} МэВ")
    print(f"Период полураспада: {result['half_life_s']:.1f} с")
    print(f"Класс: {result['stability_class']}")
    # Ожидаемый период Tc-99m ~ 6 часов = 21600 с