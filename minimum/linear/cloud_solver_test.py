#!/usr/bin/python3
from threading import Thread
import unittest
from xmlrpc.client import Binary
from xmlrpc.server import SimpleXMLRPCServer

import numpy as np

from minimum.linear.cloud_solver import CloudSolver
from minimum.linear.ip import IP
from minimum.linear.solver import Solver
from minimum.linear.sum import Sum


class FakeNeosService:
	def __init__(self, solution):
		self.solution = solution

	def ping(self):
		return "NeosServer is alive\n"

	def listSolversInCategory(self, category):
		return [category + ".CPLEX.MPS"]

	def printQueue(self):
		return "<empty queue>"

	def submitJob(self, xml):
		return 1, "abcd"

	def getIntermediateResults(self, jobNumber, password, offset):
		return Binary(b"progress"), 0

	def getJobStatus(self, jobNumber, password):
		if jobNumber == 1 and password == "abcd":
			return "Done"
		else:
			return "Error"

	def getFinalResults(self, jobNumber, password):
		msg = b"MIP - Integer optimal solution"
		for i in range(len(self.solution)):
			if self.solution[i] != 0:
				msg += b'C00' + str(i).encode() + b'      ' + str(
				    self.solution[i]).encode() + b'\n'
		return Binary(msg)


class FakeNeosServer:
	def __init__(self, solution):
		self.service = FakeNeosService(solution)

		# Create server
		self.server = SimpleXMLRPCServer(
		    ("localhost", 8000), logRequests=False)
		self.server.register_introspection_functions()
		self.server.register_instance(self.service)
		self.server.register_multicall_functions()

		self.thread = Thread(target=self.run)
		self.thread.start()

	def run(self):
		# Run the server's main loop
		self.server.serve_forever()

	def stop(self):
		if self.server is not None:
			self.server.shutdown()
			self.server.server_close()
			self.server = None


class TestCloudSolver(unittest.TestCase):
	def setUp(self) -> None:
		self.n = 3
		n = self.n
		self.ip = IP()
		x = self.ip.add_boolean_cube(n * n, n * n, n * n)
		self.x = x

		for i in range(n * n):
			for j in range(n * n):
				s = Sum(0)
				for k in range(n * n):
					s += x[i][j][k]
				self.ip.add_constraint(s == 1)

		for i in range(n * n):
			for k in range(n * n):
				s = Sum(0)
				for j in range(n * n):
					s += x[i][j][k]
				self.ip.add_constraint(s == 1)

		for j in range(n * n):
			for k in range(n * n):
				s = Sum(0)
				for i in range(n * n):
					s += x[i][j][k]
				self.ip.add_constraint(s == 1)

		for i in range(n):
			for j in range(n):
				for k in range(n * n):
					s = Sum(0)
					for i2 in range(n):
						for j2 in range(n):
							s += x[n * i + i2][n * j + j2][k]
					self.ip.add_constraint(s == 1)

	def verify_solution(self):
		n = self.n
		x = self.x
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

	@unittest.skip("Not a good unit test. Sends email.")
	def test_cloud(self) -> None:
		solver = CloudSolver()
		solutions = solver.solutions(self.ip)
		self.assertTrue(solutions.get())
		self.verify_solution()

	@unittest.skip("Slow.")
	def test_fake_cloud(self) -> None:
		solver = Solver()
		solutions = solver.solutions(self.ip)
		self.assertTrue(solutions.get())
		self.verify_solution()

		fake_server = FakeNeosServer(self.ip.solution)
		self.ip.solution = None

		cloud_solver = CloudSolver("http://localhost:8000")
		solutions = cloud_solver.solutions(self.ip)
		self.assertTrue(solutions.get())
		self.verify_solution()

		fake_server.stop()


if __name__ == '__main__':
	unittest.main()
