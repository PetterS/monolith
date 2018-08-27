import numbers
from typing import no_type_check, Dict, Optional, List

import minimum.linear.ip
import minimum.linear.variable


class Sum:
	def __init__(self, constant=0) -> None:
		self.constant = constant
		self.coefficients: Dict[int, float] = {}
		self.ip: Optional[minimum.linear.ip.IP] = None

	def copy(self) -> "Sum":
		new = Sum(self.constant)
		new.coefficients = self.coefficients.copy()
		new.ip = self.ip
		new._clean()
		return new

	def value(self) -> float:
		assert self.ip
		result = self.constant
		for index, coefficient in self.coefficients.items():
			result += coefficient * self.ip._value(index)
		return result

	def __str__(self) -> str:
		self._clean()
		s: List[str] = []
		if self.constant != 0:
			s += [str(self.constant)]
		indices = [k for k in self.coefficients.keys()]
		indices.sort()
		for index in indices:
			s += [str(self.coefficients[index]) + "*x" + str(index)]
		return " + ".join(s)

	def __repr__(self) -> str:
		return "Sum(" + self.__str__() + ")"

	def __iadd__(self, other) -> "Sum":
		if isinstance(other, minimum.linear.variable.Variable):
			self.coefficients[other.index] = self.coefficients.get(
			    other.index, 0) + 1
			self._combine_ip(other.ip)
		elif isinstance(other, Sum):
			self.constant += other.constant
			for index, coefficient in other.coefficients.items():
				self.coefficients[index] = self.coefficients.get(
				    index, 0) + coefficient
			self._combine_ip(other.ip)
		else:
			self.constant += other
		return self

	def __isub__(self, other) -> "Sum":
		if isinstance(other, minimum.linear.variable.Variable):
			self.coefficients[other.index] = self.coefficients.get(
			    other.index, 0) - 1
			self._combine_ip(other.ip)
		elif isinstance(other, Sum):
			self.constant -= other.constant
			for index, coefficient in other.coefficients.items():
				self.coefficients[index] = self.coefficients.get(
				    index, 0) - coefficient
			self._combine_ip(other.ip)
		else:
			self.constant -= other
		return self

	def __add__(self, other) -> "Sum":
		result = self.copy()
		result += other
		return result

	def __radd__(self, other) -> "Sum":
		return self + other

	def __pos__(self) -> "Sum":
		return self.copy()

	def __sub__(self, other) -> "Sum":
		result = self.copy()
		result -= other
		return result

	def __rsub__(self, other) -> "Sum":
		result = Sum(0)
		result += other
		result -= self
		return result

	def __neg__(self) -> "Sum":
		result = Sum(0)
		result -= self
		return result

	def __imul__(self, other: float) -> "Sum":
		assert isinstance(other, numbers.Number)
		self.constant *= other
		for index, _ in self.coefficients.items():
			self.coefficients[index] *= other
		return self

	def __mul__(self, other: float) -> "Sum":
		result = self.copy()
		result *= other
		return result

	def __rmul__(self, other: float) -> "Sum":
		return self * other

	def __le__(self, other) -> "minimum.linear.constraint.Constraint":
		s = self.copy()
		s -= other
		return minimum.linear.constraint.Constraint(None, s, 0)

	@no_type_check
	def __eq__(self, other) -> "minimum.linear.constraint.Constraint":
		s = self.copy()
		s -= other
		return minimum.linear.constraint.Constraint(0, s, 0)

	def __ge__(self, other) -> "minimum.linear.constraint.Constraint":
		s = self.copy()
		s -= other
		return minimum.linear.constraint.Constraint(0, s, None)

	def _combine_ip(self, ip) -> None:
		if self.ip is None:
			self.ip = ip
		assert ip is None or ip == self.ip, (
		    "Can not combine variables from different IPs.")

	def _clean(self) -> None:
		new_coefficients = {}
		for index, coefficient in self.coefficients.items():
			if coefficient != 0:
				new_coefficients[index] = coefficient
		self.coefficients = new_coefficients
