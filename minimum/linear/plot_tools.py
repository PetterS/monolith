import re
from typing import List, Tuple


def extract_primal_dual_results(
        s: str) -> Tuple[List[float], List[float], List[float], List[float]]:
	"""Parses a string containing logs from a primal-dual solver run and returns
	the iteration count, objective function values, the infeasibility, and
	time stamps."""

	float_group = r'(-?\d\.\d+e[+-]\d\d)'
	pattern = r'^\s*(\d+)\s+' + float_group + r'\s+' + float_group + r'\s+' + float_group + r'\s+' + float_group + r'\s+' + float_group + r'\s*$'
	iteration = []
	objective = []
	infeasibility = []
	time = []
	for match in re.finditer(pattern, s, re.MULTILINE):
		iteration.append(float(match.group(1)))
		objective.append(float(match.group(2)))
		infeasibility.append(float(match.group(5)))
		time.append(float(match.group(6)))
	return iteration, objective, infeasibility, time


def extract_scs_results(
        s: str) -> Tuple[List[float], List[float], List[float], List[float]]:
	"""Parses a string containing logs from an SCS solver run and returns
	the iteration count, objective function values, the duality gap, and
	time stamps."""

	float_group = r'(-?\d\.\d+e[+-]\d\d)'

	try:
		setup_time_pattern = r'Setup time: ' + float_group + 's'
		match = re.search(setup_time_pattern, s)
		assert match, "Did not find setup time."
		setup_time = float(match.group(1))
	except AttributeError:
		return [], [], [], []

	pattern = r'^\s*(\d+)\|'
	for i in range(7):
		pattern += r'\s+' + float_group
	pattern += r'\s*$'

	iteration = []
	objective = []
	gap = []
	time = []
	for match in re.finditer(pattern, s, re.MULTILINE):
		iteration.append(float(match.group(1)))
		objective.append(float(match.group(5)))
		gap.append(float(match.group(4)))
		time.append(float(match.group(8)) + setup_time)
	return iteration, objective, gap, time
