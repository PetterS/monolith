#!/usr/bin/python3
import unittest

from minimum.linear.ip import FakeIP
from minimum.linear.sum import Sum


class TestVariable(unittest.TestCase):
	def setUp(self) -> None:
		self.ip = FakeIP()
		self.v0 = self.ip.add_fixed_variable(0.0)
		self.v1 = self.ip.add_fixed_variable(1.0)
		self.v2 = self.ip.add_fixed_variable(2.0)
		self.v5 = self.ip.add_fixed_variable(5.0)
		self.v6 = self.ip.add_fixed_variable(6.0)

	def test_value(self) -> None:
		self.assertEqual(self.v5.value(), 5)

	def test_add(self) -> None:
		self.assertEqual((1 + self.v0).value(), 1.0)
		self.assertEqual((self.v0 + 1).value(), 1.0)
		self.assertEqual((self.v0 + self.v1).value(), 1.0)
		self.assertEqual((self.v2 + self.v5 + self.v2).value(), 9.0)
		s8 = self.v2 + self.v6
		self.assertEqual((s8 + s8).value(), 16.0)
		self.assertEqual(s8.value(), 8.0)
		self.assertEqual((s8 + self.v1).value(), 9.0)
		self.assertEqual((self.v1 + s8).value(), 9.0)
		self.assertEqual((s8 + 5).value(), 13.0)
		self.assertEqual((5 + s8).value(), 13.0)

	def test_subtract(self) -> None:
		self.assertEqual((self.v0 - self.v1).value(), -1.0)
		self.assertEqual((self.v2 - self.v5 - self.v2).value(), -5.0)
		s8 = self.v2 + self.v1 + 5
		self.assertEqual((s8 - s8).value(), 0.0)
		self.assertEqual(s8.value(), 8.0)
		self.assertEqual((s8 - self.v1).value(), 7.0)
		self.assertEqual((self.v1 - s8).value(), -7.0)
		self.assertEqual((s8 - 5).value(), 3.0)
		self.assertEqual((5 - s8).value(), -3.0)

	def test_negate(self) -> None:
		s3 = self.v1 + self.v2
		self.assertEqual(+self.v6, 6)
		self.assertEqual(+s3, 3)
		self.assertEqual(-self.v6, -6)
		self.assertEqual(-s3, -3)

	def test_multiply(self) -> None:
		self.assertEqual((4 * self.v5).value(), 20.0)
		self.assertEqual((self.v5 * 4).value(), 20.0)
		s8 = self.v2 + self.v6
		self.assertEqual((s8 * 4).value(), 32.0)
		self.assertEqual((4.0 * s8).value(), 32.0)

		# The tests below are caught by mypy if the comments
		# are removed.
		with self.assertRaises(AssertionError):
			self.v1 * self.v0  # type: ignore
		with self.assertRaises(AssertionError):
			self.v1 * s8  # type: ignore
		with self.assertRaises(AssertionError):
			s8 * self.v0  # type: ignore
		with self.assertRaises(AssertionError):
			s8 * s8  # type: ignore

	def test_str(self) -> None:
		self.assertEqual(
		    str(self.v1 + self.v1 + 5 + self.v0), "5 + 1*x0 + 2*x1")

	def test_constraint(self) -> None:
		self.assertEqual(str(self.v0 <= self.v1 + 5), "1*x0 + -1*x1 ≤ 5")
		self.assertEqual(str(self.v0 >= self.v1 + 5), "5 ≤ 1*x0 + -1*x1")
		self.assertEqual(str(self.v0 == self.v1 + 5), "5 ≤ 1*x0 + -1*x1 ≤ 5")
		s = self.v0 + self.v1
		self.assertEqual(str(s <= 1), "1*x0 + 1*x1 ≤ 1")
		self.assertEqual(str(s >= 1), "1 ≤ 1*x0 + 1*x1")
		self.assertEqual(str(s == 1), "1 ≤ 1*x0 + 1*x1 ≤ 1")

		with self.assertRaises(AssertionError):
			Sum(0) >= 1
		Sum(0) <= 1


if __name__ == '__main__':
	unittest.main()
