"""Проверки контракта модели; достаточно стандартной библиотеки Python."""

import unittest

from course_motion.model import position


class TestPosition(unittest.TestCase):
    def test_initial_position(self):
        self.assertAlmostEqual(position(0.0, x0_m=7.0), 7.0)

    def test_uniform_motion(self):
        self.assertAlmostEqual(position(4.0, x0_m=1.0, v0_m_s=2.0, a_m_s2=0.0), 9.0)

    def test_accelerated_motion(self):
        # Независимый ручной расчет: 1 + 2*4 + 0.5*0.5*4**2 = 13.
        self.assertAlmostEqual(position(4.0), 13.0)

    def test_negative_acceleration(self):
        self.assertAlmostEqual(position(2.0, x0_m=0.0, v0_m_s=3.0, a_m_s2=-1.0), 4.0)

    def test_negative_time_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "non-negative"):
            position(-0.1)


if __name__ == "__main__":
    unittest.main()
