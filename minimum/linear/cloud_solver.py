import os
import re
import tempfile
import xmlrpc.client

import numpy as np

from minimum.bibliotek.c_api import minimum_linear_ip_save_mps
from minimum.linear.ip import RawIP


class CloudSolutions:
	def __init__(self, ip, address):
		self.ip = ip
		self.address = address

	def get(self):
		(handle, file_name) = tempfile.mkstemp()
		os.close(handle)
		os.unlink(file_name)

		rip = RawIP(self.ip.proto)
		minimum_linear_ip_save_mps(rip.ip, file_name)
		file_name += ".mps"
		with open(file_name, "r") as f:
			mps = f.read()
		os.unlink(file_name)

		CONTACT_EMAIL = "petter.strandmark@gmail.com"
		xml = """<document>
<category>milp</category>
<solver>CPLEX</solver>
<inputMethod>MPS</inputMethod>
<email>%s</email>
<MPS>%s</MPS>
<options></options>
<post>display solution variables -</post>
<comments></comments>
</document>""" % (CONTACT_EMAIL, mps)

		# Connect to NEOS
		neos = xmlrpc.client.ServerProxy(self.address)
		# Verify the connection was successful
		test = neos.ping()
		if test != "NeosServer is alive\n":
			raise RuntimeError("Could not make connection to the NEOS server")

		print(neos.system.listMethods())
		print(neos.listSolversInCategory("milp"))
		print("NEOS status:")
		print(neos.printQueue())
		# print(neos.getSolverTemplate("milp", "Gurobi", "MPS"))
		# print(neos.system.methodSignature("submitJob"))

		(jobNumber, password) = neos.submitJob(xml)
		# print("Submitted jobNumber =", jobNumber)
		if jobNumber == 0:
			raise RuntimeError("NEOS error: %s" % password)

		offset = 0
		status = ""
		# Print out partial job output while job is running
		while status != "Done":
			(msg, offset) = neos.getIntermediateResults(
			    jobNumber, password, offset)
			# print(msg.data.decode())
			status = neos.getJobStatus(jobNumber, password)

		print("Job done.")
		# Print out the final result
		msg = neos.getFinalResults(jobNumber, password).data.decode()
		# print(msg)

		if "MIP - Integer optimal solution" not in msg:
			return False

		self.ip.solution = np.zeros(
		    len(self.ip.proto.variable), dtype=np.float64)
		for match in re.finditer(r"^C(\d+)\s+(-?\d+.?\d*)$", msg,
		                         re.MULTILINE):
			column = int(match.group(1))
			value = float(match.group(2))
			self.ip.solution[column] = value

		return True


class CloudSolver:
	def __init__(self, address="https://neos-server.org:3333"):
		self.address = address

	def solutions(self, ip):
		return CloudSolutions(ip, self.address)
