"""Лабораторная работа № 1: прямой расчёт внутренней баллистики."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from math import pi
from numbers import Number
from typing import Any

import numpy as np
from pyballistics import ozvb_lagrange


# Исходные данные варианта № 17 в единицах СИ.
CALIBER_M = 0.1524
PROJECTILE_MASS_KG = 43.56
TARGET_VELOCITY_M_S = 700.0
BORE_AREA_M2 = 0.26 * pi * CALIBER_M**2
FORCING_PRESSURE_PA = 30.0e6
FRICTION_FACTOR = 1.02
IGNITION_PRESSURE_PA = 5.0e6
NORMAL_TEMPERATURE_K = 293.15
MAX_BARREL_LENGTH_M = 50.0 * CALIBER_M

DEFAULT_CELL_COUNT = 150
DEFAULT_CFL = 0.9
MAX_MODEL_PRESSURE_PA = 700.0e6
MAX_MODEL_STEPS = 100_000


@dataclass(frozen=True)
class PowderCharge:
    """Одна навеска табличного пороха."""

    database_name: str
    mass_kg: float


@dataclass(frozen=True)
class BallisticInput:
    """Явное описание входных данных одного прямого расчёта."""

    name: str
    chamber_volume_m3: float
    charges: tuple[PowderCharge, ...]
    temperature_k: float = NORMAL_TEMPERATURE_K
    stop_position_m: float = MAX_BARREL_LENGTH_M
    cell_count: int = DEFAULT_CELL_COUNT
    cfl: float = DEFAULT_CFL


class RunStatus(StrEnum):
    """Численный статус, не зависящий от выполнения ограничений варианта."""

    SUCCESS = "success"
    SOLVER_ERROR = "solver_error"
    INVALID_RESULT = "invalid_result"


@dataclass(frozen=True)
class ModelRun:
    """Контролируемый результат вызова ``ozvb_lagrange``."""

    input_data: BallisticInput
    status: RunStatus
    stop_reason: str
    execution_time_s: float
    layers: tuple[dict[str, Any], ...]
    model_options: dict[str, Any]
    message: str

    @property
    def is_successful(self) -> bool:
        return self.status is RunStatus.SUCCESS


def build_model_options(input_data: BallisticInput) -> dict[str, Any]:
    """Преобразовать явное описание решения во входной словарь pyballistics."""

    return {
        "powders": [
            {"omega": charge.mass_kg, "dbname": charge.database_name}
            for charge in input_data.charges
        ],
        "init_conditions": {
            "q": PROJECTILE_MASS_KG,
            "d": CALIBER_M,
            "W_0": input_data.chamber_volume_m3,
            "T_0": input_data.temperature_k,
            "phi_1": FRICTION_FACTOR,
            "p_0": FORCING_PRESSURE_PA,
            "S": BORE_AREA_M2,
        },
        "igniter": {"p_ign_0": IGNITION_PRESSURE_PA},
        "meta_lagrange": {
            "n_cells": input_data.cell_count,
            "CFL": input_data.cfl,
        },
        "stop_conditions": {
            "x_p": input_data.stop_position_m,
            "p_max": MAX_MODEL_PRESSURE_PA,
            "steps_max": MAX_MODEL_STEPS,
        },
    }


def _validate_layers(layers: object) -> tuple[dict[str, Any], ...]:
    """Проверить структуру слоёв и конечность всех численных данных."""

    if not isinstance(layers, list) or not layers:
        raise ValueError("модель не вернула расчётные слои")

    required_keys = {"t", "step_count", "x", "u", "T", "rho", "p", "T_w", "k"}
    checked_layers: list[dict[str, Any]] = []
    for layer_index, layer in enumerate(layers):
        if not isinstance(layer, dict):
            raise TypeError(f"слой {layer_index} не является словарём")
        missing = required_keys.difference(layer)
        if missing:
            raise ValueError(
                f"в слое {layer_index} отсутствуют поля: {', '.join(sorted(missing))}"
            )

        for key, value in layer.items():
            if isinstance(value, np.ndarray):
                if value.size == 0:
                    raise ValueError(f"пустой массив {key!r} в слое {layer_index}")
                if not np.issubdtype(value.dtype, np.number):
                    raise TypeError(f"массив {key!r} в слое {layer_index} нечисловой")
                if not np.all(np.isfinite(value)):
                    raise ValueError(
                        f"NaN или бесконечность в поле {key!r} слоя {layer_index}"
                    )
            elif isinstance(value, Number) and not np.isfinite(float(value)):
                raise ValueError(
                    f"NaN или бесконечность в поле {key!r} слоя {layer_index}"
                )
        checked_layers.append(layer)

    return tuple(checked_layers)


def run_model_safely(input_data: BallisticInput) -> ModelRun:
    """Безопасно вызвать модель и преобразовать любую ошибку в ``ModelRun``."""

    options = build_model_options(input_data)
    try:
        raw_result = ozvb_lagrange(options)
    except Exception as error:  # защитная граница для ошибок вне pyballistics
        return ModelRun(
            input_data=input_data,
            status=RunStatus.SOLVER_ERROR,
            stop_reason="exception",
            execution_time_s=0.0,
            layers=(),
            model_options=options,
            message=f"ozvb_lagrange выбросила {type(error).__name__}: {error}",
        )

    if not isinstance(raw_result, dict):
        return ModelRun(
            input_data=input_data,
            status=RunStatus.INVALID_RESULT,
            stop_reason="invalid_result",
            execution_time_s=0.0,
            layers=(),
            model_options=options,
            message="ozvb_lagrange вернула результат, не являющийся словарём",
        )

    stop_reason = str(raw_result.get("stop_reason", "missing"))
    execution_time = float(raw_result.get("execution_time", 0.0))
    if stop_reason == "error":
        return ModelRun(
            input_data=input_data,
            status=RunStatus.SOLVER_ERROR,
            stop_reason=stop_reason,
            execution_time_s=execution_time,
            layers=(),
            model_options=options,
            message=str(raw_result.get("error_message", "неизвестная ошибка модели")),
        )

    try:
        layers = _validate_layers(raw_result.get("layers"))
    except (TypeError, ValueError, OverflowError) as error:
        return ModelRun(
            input_data=input_data,
            status=RunStatus.INVALID_RESULT,
            stop_reason=stop_reason,
            execution_time_s=execution_time,
            layers=(),
            model_options=options,
            message=f"некорректные расчётные данные: {error}",
        )

    return ModelRun(
        input_data=input_data,
        status=RunStatus.SUCCESS,
        stop_reason=stop_reason,
        execution_time_s=execution_time,
        layers=layers,
        model_options=options,
        message="расчёт завершён; ограничения на этом этапе не классифицируются",
    )


def main() -> None:
    """Выполнить пробный безопасный запуск квазиодномерной модели."""

    trial_input = BallisticInput(
        name="Пробный расчёт варианта 17",
        chamber_volume_m3=0.025,
        charges=(PowderCharge(database_name="11/7", mass_kg=11.5),),
    )
    run = run_model_safely(trial_input)
    print(f"{run.input_data.name}: status={run.status}, stop_reason={run.stop_reason}")
    print(f"Слоёв: {len(run.layers)}; время модели: {run.execution_time_s:.4f} с")
    print(f"Диагностика: {run.message}")


if __name__ == "__main__":
    main()
