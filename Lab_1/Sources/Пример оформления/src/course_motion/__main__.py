"""Запуск: uv run python -m course_motion из корня проекта."""

from course_motion.model import position


def main() -> None:
    for t_s in range(5):
        x_m = position(t_s)
        print(f"t = {t_s:2} s; x = {x_m:5.2f} m")


if __name__ == "__main__":
    main()
