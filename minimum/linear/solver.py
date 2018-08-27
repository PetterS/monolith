from enum import Enum

import numpy as np  # type: ignore

from minimum.bibliotek.c_api import error_string, minimum_linear_solver_default, minimum_linear_solver_destroy, minimum_linear_solver_minisat, minimum_linear_solutions, minimum_linear_solutions_get, minimum_linear_solutions_destroy
import minimum.linear.ip


class Type(Enum):
	default = 1
	minisat = 2


class Solutions:
	def __init__(self, solver, ip):
		self.ip = ip  # type: minimum.linear.ip.IP
		self.rip = minimum.linear.ip.RawIP(ip.proto)
		self.solver = solver
		self.solutions = minimum_linear_solutions(solver.solver, self.rip.ip)

	def __del__(self):
		if self.solutions is not None:
			minimum_linear_solutions_destroy(self.solutions)
			self.solutions = None

	def get(self):
		if self.ip.solution is None:
			self.ip.solution = np.zeros(
			    len(self.ip.proto.variable), dtype=np.float64)
		res = minimum_linear_solutions_get(self.solutions, self.rip.ip,
		                                   self.ip.solution)
		self._check_result(res)
		return res == 0

	def _check_result(self, res: int) -> None:
		if res == 2:
			assert False, error_string()
		assert res == 0 or res == 1, "res=" + str(res)


class Solver:
	def __init__(self, type: Type = Type.default) -> None:

		self.solver = None
		if type == Type.default:
			self.solver = minimum_linear_solver_default()
		elif type == Type.minisat:
			self.solver = minimum_linear_solver_minisat()
		else:
			assert False, "Unknown solver type."

	def __del__(self):
		if self.solver is not None:
			minimum_linear_solver_destroy(self.solver)
			self.solver = None

	def solutions(self, ip) -> Solutions:
		return Solutions(self, ip)
