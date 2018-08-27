from typing import Dict, Optional, List, Union

import numpy as np  # type: ignore

import minimum.bibliotek.c_api
from minimum.bibliotek.c_api import minimum_linear_ip, minimum_linear_ip_destroy, Ip_p
import minimum.linear.constraint
import minimum.linear.proto
import minimum.linear.solver
import minimum.linear.sum
import minimum.linear.variable


class RawIP:
	def __init__(self, ip_proto: minimum.linear.proto.IP) -> None:
		# TODO: Why is this type not deduced by mypy?
		self.ip = minimum_linear_ip(ip_proto.SerializeToString(
		))  # tasype: Optional[minimum.bibliotek.c_api.Ip_p]

	def clear(self) -> None:
		if self.ip != Ip_p(0):
			minimum_linear_ip_destroy(self.ip)
			self.ip = Ip_p(0)

	def __del__(self) -> None:
		self.clear()


class IP:
	def __init__(self) -> None:
		self.proto = minimum.linear.proto.IP()  # type: minimum.linear.proto.IP
		# We keep the solution outside the proto for now.
		self.solution = None  # type: np.ndarray

	def add_variable(
	    self, is_integer: bool = False) -> minimum.linear.variable.Variable:
		index = len(self.proto.variable)
		var = self.proto.variable.add()
		var.bound.lower = -np.inf
		var.bound.upper = np.inf
		var.cost = 0
		if is_integer:
			var.type = minimum.linear.proto.Variable.INTEGER

		return minimum.linear.variable.Variable(index, self)

	def add_boolean(self) -> minimum.linear.variable.Variable:
		x = self.add_variable(True)
		self.add_bounds(0, x, 1)
		return x

	def add_boolean_vector(self,
	                       n: int) -> List[minimum.linear.variable.Variable]:
		v = []
		for i in range(n):
			v.append(self.add_boolean())
		return np.array(v, dtype=object)

	def add_boolean_grid(
	    self, m: int, n: int) -> List[List[minimum.linear.variable.Variable]]:
		v = []
		for i in range(m):
			v.append(self.add_boolean_vector(n))
		return np.array(v, dtype=object)

	def add_boolean_cube(
	    self, m: int, n: int,
	    o: int) -> List[List[List[minimum.linear.variable.Variable]]]:
		v = []
		for i in range(m):
			v.append(self.add_boolean_grid(n, o))
		return np.array(v, dtype=object)

	def add_bounds(self, lb: Optional[float],
	               x: minimum.linear.variable.Variable,
	               ub: Optional[float]) -> None:
		var = self.proto.variable[x.index]
		if lb is not None:
			var.bound.lower = max(lb, var.bound.lower)
		if ub is not None:
			var.bound.upper = min(ub, var.bound.upper)
		assert var.bound.lower <= var.bound.upper, "Can not add impossible contraint."

	def add_constraint(self, c: minimum.linear.constraint.Constraint) -> None:
		assert isinstance(c, minimum.linear.constraint.Constraint)
		assert len(c.sum.coefficients) > 0

		if len(c.sum.coefficients) == 1:
			index, coefficient = list(c.sum.coefficients.items())[0]
			assert coefficient != 0
			lb, ub = None, None
			if c.lb is not None:
				lb = c.lb / coefficient
			if c.ub is not None:
				ub = c.ub / coefficient
			var = minimum.linear.variable.Variable(index, self)
			self.add_bounds(lb, var, ub)
			return

		constraint = self.proto.constraint.add()
		for index, coefficient in c.sum.coefficients.items():
			entry = constraint.sum.add()
			entry.coefficient = coefficient
			entry.variable = index

		if c.lb is not None:
			constraint.bound.lower = c.lb
		else:
			constraint.bound.lower = -np.inf

		if c.ub is not None:
			constraint.bound.upper = c.ub
		else:
			constraint.bound.upper = np.inf

	def add_objective(self, s: Union[minimum.linear.variable.Variable,
	                                 minimum.linear.sum.Sum, float]) -> None:
		if isinstance(s, minimum.linear.variable.Variable):
			assert s.ip == self

			self.proto.variable[s.index].cost += 1
		elif isinstance(s, minimum.linear.sum.Sum):
			assert s.ip == self

			self.proto.objective_constant += s.constant
			for index, coefficient in s.coefficients.items():
				self.proto.variable[index].cost = coefficient
		else:
			self.proto.objective_constant += s

	def _value(self, index: int) -> float:
		assert self.solution is not None, "There is no solution."
		return self.solution[index]


class FakeIP(IP):
	"Fake IP class used for testing."

	def __init__(self) -> None:
		self.variables: Dict[int, float] = {}

	def add_fixed_variable(self,
	                       value: float) -> minimum.linear.variable.Variable:
		index = len(self.variables)
		self.variables[index] = value
		return minimum.linear.variable.Variable(index, self)

	def _value(self, index) -> float:
		return self.variables[index]
