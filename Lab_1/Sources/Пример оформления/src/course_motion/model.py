"""Одномерное движение тележки с постоянным ускорением."""


def position(
    t_s: float,
    x0_m: float = 1.0,
    v0_m_s: float = 2.0,
    a_m_s2: float = 0.5,
) -> float:
    """Вернуть координату в метрах при t_s >= 0. Все аргументы — конечные числа."""
    if t_s < 0:
        raise ValueError("t_s must be non-negative")
    return x0_m + v0_m_s * t_s + 0.5 * a_m_s2 * t_s**2
