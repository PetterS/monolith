import numbers
from typing import no_type_check

import minimum.linear.constraint
import minimum.linear.ip
import minimum.linear.sum


class Variable:
	def __init__(self, index: int, ip: "minimum.linear.ip.IP" = None) -> None:
		"""Normally not used directly. Instead, Variables
		are created from the IP class."""
		self.index = index
		self.ip = ip

	def value(self) -> float:
		assert self.ip is not None
		return self.ip._value(self.index)

	def __add__(self, other) -> "minimum.linear.sum.Sum":
		s = minimum.linear.sum.Sum(0)
		s += self
		s += other
		return s

	def __radd__(self, other) -> "minimum.linear.sum.Sum":
		return self + other

	def __pos__(self) -> "Variable":
		return self

	def __sub__(self, other) -> "minimum.linear.sum.Sum":
		s = minimum.linear.sum.Sum(0)
		s += self
		s -= other
		return s

	def __neg__(self) -> "minimum.linear.sum.Sum":
		s = minimum.linear.sum.Sum(0)
		s -= self
		return s

	def __mul__(self, other: float) -> "minimum.linear.sum.Sum":
		assert isinstance(other, numbers.Number)
		s = minimum.linear.sum.Sum(0)
		s.coefficients[self.index] = s.coefficients.get(self.index, 0) + other
		s.ip = self.ip
		return s

	def __rmul__(self, other: float) -> "minimum.linear.sum.Sum":
		return self * other

	def __le__(self, other) -> "minimum.linear.constraint.Constraint":
		s = minimum.linear.sum.Sum(0)
		s += self
		s -= other
		return minimum.linear.constraint.Constraint(None, s, 0)

	@no_type_check
	def __eq__(self, other) -> "minimum.linear.constraint.Constraint":
		s = minimum.linear.sum.Sum(0)
		s += self
		s -= other
		return minimum.linear.constraint.Constraint(0, s, 0)

	def __ge__(self, other) -> "minimum.linear.constraint.Constraint":
		s = minimum.linear.sum.Sum(0)
		s += self
		s -= other
		return minimum.linear.constraint.Constraint(0, s, None)

	def __repr__(self) -> str:
		return "x" + str(self.index)
