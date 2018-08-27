#!/usr/bin/python3
from typing import no_type_check
import unittest

import numpy as np  # type: ignore

from minimum.linear.ip import IP
from minimum.linear.solver import Solver, Type
from minimum.linear.sum import Sum


class TestIP(unittest.TestCase):
	def test_basic(self) -> None:
		ip = IP()
		x = ip.add_variable()
		y = ip.add_variable()
		ip.add_constraint(x + y <= 1)
		ip.add_constraint(0 <= x)
		ip.add_constraint(x <= 1)
		ip.add_constraint(0 <= y)
		ip.add_constraint(y <= 1)
		ip.add_objective(-x)

		solver = Solver()
		self.assertTrue(solver.solutions(ip).get())

		self.assertEqual(x.value(), 1)
		self.assertEqual(y.value(), 0)

	def test_basic_bounds(self) -> None:
		ip = IP()
		x = ip.add_variable()
		y = ip.add_variable()
		ip.add_bounds(0, x, 1)
		ip.add_bounds(0, y, 1)
		ip.add_constraint(x + y <= 1)
		ip.add_objective(-x)
		self.assertTrue(Solver().solutions(ip).get())

		self.assertEqual(x.value(), 1)
		self.assertEqual(y.value(), 0)

	def test_basic_variable_type(self) -> None:
		ip = IP()
		x = ip.add_variable(True)
		y = ip.add_boolean()
		ip.add_bounds(0, x, 0.5)
		ip.add_objective(-x - y)
		self.assertTrue(Solver().solutions(ip).get())

		self.assertEqual(x.value(), 0)
		self.assertEqual(y.value(), 1)

	def test_sudoku(self) -> None:
		n = 3
		ip = IP()
		x = ip.add_boolean_cube(n * n, n * n, n * n)

		for i in range(n * n):
			for j in range(n * n):
				s = Sum(0)
				for k in range(n * n):
					s += x[i][j][k]
				ip.add_constraint(s == 1)

		for i in range(n * n):
			for k in range(n * n):
				s = Sum(0)
				for j in range(n * n):
					s += x[i][j][k]
				ip.add_constraint(s == 1)

		for j in range(n * n):
			for k in range(n * n):
				s = Sum(0)
				for i in range(n * n):
					s += x[i][j][k]
				ip.add_constraint(s == 1)

		for i in range(n):
			for j in range(n):
				for k in range(n * n):
					s = Sum(0)
					for i2 in range(n):
						for j2 in range(n):
							s += x[n * i + i2][n * j + j2][k]
					ip.add_constraint(s == 1)

		solver = Solver()
		solutions = solver.solutions(ip)
		self.assertTrue(solutions.get())

		d = np.zeros((n * n, n * n))
		for i in range(n * n):
			for j in range(n * n):
				for k in range(n * n):
					d[i, j] += (k + 1) * x[i][j][k].value()

		for i in range(n * n):
			self.assertEqual(len(set(d[:, i])), n * n)
			self.assertEqual(len(set(d[i, :])), n * n)
		for i in range(n):
			for j in range(n):
				self.assertEqual(
				    len(set(d[n * i:n * i + n, n * j:n * j + n].flatten())),
				    n * n)

	def test_minisat(self) -> None:
		ip = IP()
		x = ip.add_boolean()
		y = ip.add_boolean()
		ip.add_objective(-x - y)
		solver = Solver(Type.minisat)
		self.assertTrue(solver.solutions(ip).get())
		self.assertAlmostEqual(x.value(), 1)
		self.assertAlmostEqual(y.value(), 1)

	def test_next(self) -> None:
		ip = IP()
		x = ip.add_boolean()
		y = ip.add_boolean()
		z = ip.add_boolean()
		ip.add_constraint(x + y + z >= 1)
		solver = Solver()
		solutions = solver.solutions(ip)
		n = 0
		while solutions.get():
			n += 1
		self.assertEqual(n, 7)

	def test_impossible_constraint(self) -> None:
		ip = IP()
		x = ip.add_boolean()
		with self.assertRaises(AssertionError):
			ip.add_constraint(x >= 2)
		with self.assertRaises(AssertionError):
			ip.add_constraint(0 * x == 0)

	def test_objective_constant(self) -> None:
		ip = IP()
		ip.add_objective(3)

	@no_type_check
	def test_sum(self) -> None:
		ip = IP()
		x = ip.add_boolean_grid(5, 5)
		self.assertIsInstance(sum(x[:][1]), Sum)
		self.assertIsInstance(sum(x[2][:]), Sum)

		x3 = ip.add_boolean_cube(5, 5, 6)
		self.assertIsInstance(sum(sum(x3[:][2][:])), Sum)

	def test_infeasible(self) -> None:
		ip = IP()
		x = ip.add_variable()
		y = ip.add_variable()
		ip.add_constraint(x + y <= 1)
		ip.add_constraint(x + y >= 2)
		ip.add_constraint(0 <= x)
		ip.add_constraint(x <= 3)
		ip.add_constraint(0 <= y)
		ip.add_constraint(y <= 3)
		ip.add_objective(-x)

		solver = Solver()
		self.assertFalse(solver.solutions(ip).get())


if __name__ == '__main__':
	unittest.main()
