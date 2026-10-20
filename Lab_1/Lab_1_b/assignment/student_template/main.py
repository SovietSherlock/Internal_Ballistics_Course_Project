"""Точка входа лабораторной работы № 1.

Запуск:
    uv run python main.py

После выполнения в figures/ должны находиться все рисунки report.md.
Шаблонный рисунок требуется удалить после реализации собственной задачи.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
FIGURES_DIR = ROOT / "figures"
RANDOM_SEED = 20260921


@dataclass(frozen=True)
class Variant:
    """Данные одной строки таблицы вариантов.

    Замените поля и при необходимости расширьте структуру.
    """

    number: int
    d_m: float
    q_kg: float
    v_target_m_s: float
    gun_type: str
    model: str
    criterion: str
    p_max_pa: float
    l_max_calibers: float


@dataclass(frozen=True)
class Candidate:
    """Один проверочный набор свободных параметров."""

    name: str
    powder: str
    omega_over_q: float
    loading_density_kg_m3: float


def configure_plots() -> None:
    """Единый читаемый стиль всех авторских графиков."""

    plt.rcParams.update(
        {
            "figure.figsize": (7.2, 4.6),
            "figure.dpi": 140,
            "savefig.dpi": 220,
            "font.size": 11,
            "axes.grid": True,
            "grid.alpha": 0.28,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def evaluate_candidate(
    candidate: Candidate, variant: Variant, *, keep_trajectories: bool = True
) -> dict[str, Any]:
    """Выполнить полный прямой расчёт и вернуть диагностический результат.

    В собственной реализации рекомендуется вернуть:
    - входные параметры и статусы температурных расчётов;
    - значения всех показателей и ограничений;
    - исходный критерий, отдельные штрафы и штрафованный критерий;
    - траектории, необходимые для рисунков.
    """

    raise NotImplementedError("Реализуйте расчёт для своего варианта")


def save_template_check() -> Path:
    """Создать нейтральный рисунок для проверки конвейера файлов.

    Этот рисунок не является результатом лабораторной работы. После того как
    собственные графики готовы, удалите функцию и ссылку на рисунок из отчёта.
    """

    x = np.linspace(0.0, 1.0, 80)
    y = x
    fig, ax = plt.subplots()
    ax.plot(x, y, color="#2f5d8a", linewidth=2, label="тестовая линия")
    ax.set_xlabel("Безразмерная координата")
    ax.set_ylabel("Безразмерная величина")
    ax.set_title("Проверка воспроизводимого построения рисунка")
    ax.legend()
    fig.tight_layout()
    target = FIGURES_DIR / "template_check.png"
    fig.savefig(target)
    plt.close(fig)
    return target


def build_all_figures() -> list[Path]:
    """Пересоздать все рисунки отчёта и вернуть их пути."""

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # TODO: создайте Variant, набор Candidate, выполните расчёты и сохраните
    # каждый рисунок report.md. Для воспроизводимости фиксируйте seed:
    # rng = np.random.default_rng(RANDOM_SEED)

    return [save_template_check()]


def main() -> None:
    configure_plots()
    print(f"pyballistics={version('pyballistics')}")
    print(f"random_seed={RANDOM_SEED}")
    created = build_all_figures()
    for path in created:
        print(f"created={path.relative_to(ROOT)}")
    print("Шаблон запущен. Замените template_check.png результатами своего расчёта.")


if __name__ == "__main__":
    main()

