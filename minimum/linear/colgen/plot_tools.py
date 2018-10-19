import collections
import io
import re
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from minimum.core import record_stream
from minimum.linear.colgen.log_pb2 import LogEntry


class Event:
	def __init__(self):
		self.time = 0
		self.iteration = 0
		self.objective = 0
		self.plot_token = 'ko'
		self.markersize = 4

	def plot(self, ax):
		ax.plot(
		    self.time,
		    self.objective,
		    self.plot_token,
		    markersize=self.markersize)


class ColgenResults:
	def __init__(self):
		self.iteration = []  # type: List[float]
		self.time = []
		self.objective = []  # type: List[float]
		self.int_objective = []  # type: Optional[List[float]]
		self.events = []  # type: List[Event]
		self.problem_size = []  # type: List[int]

	def plot_generated_columns_bars(self, ax) -> None:
		if not self.problem_size:
			return
		bar_width = np.diff(np.array(self.time + [self.time[-1]]))
		ax.bar(
		    self.time,
		    self.problem_size,
		    width=bar_width,
		    color="#dddddd",
		    align='edge')
		ax.set_ylabel("Number of roster columns in LP")

	def data_frame(self):
		return pd.DataFrame({
		    "Iteration": self.iteration,
		    "Objective": self.objective,
		    "IntegerObjective": self.int_objective,
		    "ProblemSize": self.problem_size
		},
		                    index=pd.Series(self.time, name="Time"))


def extract_colgen_results(s: str) -> ColgenResults:
	"""Parses a string containing logs from a colgen solver run and returns
	the iteration count and objective function values."""

	sci_float_group = r'(-?\d\.\d+e[+-]\d\d|nan|n/a)'
	float_group = r'(-?\d+\.?\d*)'
	int_group = r'(-?\d+)'
	sep = r'\s+'
	pattern = (
	    r'^\s*' + int_group + sep + float_group + sep + int_group + '?' + sep +
	    int_group + sep + int_group + sep + int_group + sep + int_group + sep +
	    int_group + sep + sci_float_group + sep + sci_float_group + sep +
	    sci_float_group + sep + sci_float_group + sep + r'\d+\%' + r'\s*$')
	fix_shift_pattern = r'-- Fixed for '
	fix_column_pattern = r'-- Fixing column'
	done_pattern = 'Colgen done'
	results = ColgenResults()

	previous_iterations = 0
	previous_times = 0.
	last_iter = 0
	last_objective = 0.
	last_time = 0.
	for line in s.splitlines():
		match = re.search(pattern, line)
		if match:
			last_iter = int(match.group(1)) + previous_iterations
			last_objective = float(match.group(2))
			time_str = match.group(12)
			if time_str != "n/a":
				last_time = float(time_str) + previous_times

			results.iteration.append(last_iter)
			results.objective.append(last_objective)
			results.time.append(last_time)

			if (match.group(3) is not None
			    and results.int_objective is not None):
				results.int_objective.append(float(match.group(3)))
			else:
				results.int_objective = None
			results.problem_size.append(int(match.group(4)))
		if re.search(fix_shift_pattern, line):
			event = Event()
			event.iteration = last_iter
			event.time = last_time
			event.objective = last_objective
			event.plot_token = 'ko'
			results.events.append(event)
		if re.search(fix_column_pattern, line):
			event = Event()
			event.iteration = last_iter
			event.time = last_time
			event.objective = last_objective
			event.plot_token = 'kx'
			results.events.append(event)
		if re.search(done_pattern, line):
			event = Event()
			event.iteration = last_iter
			event.time = last_time
			event.objective = last_objective
			event.plot_token = 'ro'
			event.markersize = 8
			results.events.append(event)
			previous_iterations = last_iter
			previous_times = last_time

	return results


def process_proto_log(f: io.BufferedReader) \
 -> Tuple[ColgenResults, Dict[str, np.ndarray]]:
	"""Parses a file containing logs from a colgen solver run and returns
	a matrix for every worker."""

	results = ColgenResults()
	all_data = collections.defaultdict(list)  # type: Dict[str, list]
	num_iterations = 0
	time_of_last_solution = 0
	while f.peek(1):
		entry = LogEntry()
		entry.ParseFromString(record_stream.read(f))
		for member in entry.set_partitioning.member:
			all_data[member.id].append(list(member.constraint))

		num_shifts = len(member.constraint)
		num_iterations += 1

		t = time_of_last_solution + entry.cumulative_time
		if results.time and t < results.time[-1]:
			time_of_last_solution = results.time[-1]
			t = time_of_last_solution + entry.cumulative_time
		results.time.append(t)
		results.iteration.append(num_iterations)
		results.objective.append(entry.fractional_objective)
		results.int_objective.append(entry.integer_objective)

		for proto_event in entry.event:
			event = Event()
			event.iteration = num_iterations
			event.time = t
			event.objective = entry.fractional_objective
			event_type = proto_event.WhichOneof("event_type")
			if event_type == "fix_shift":
				event.plot_token = 'ko'
				results.events.append(event)
			elif event_type == "fix_column":
				event.plot_token = 'kx'
				event.markersize = 6
				results.events.append(event)
			elif event_type == "integer_solution":
				event.plot_token = 'ro'
				event.markersize = 8
				results.events.append(event)

	matrices = {}
	All = np.zeros((num_shifts, num_iterations))
	for name, all_values in all_data.items():
		assert (len(all_values[0]) == num_shifts)
		assert (len(all_values) == num_iterations)

		A = np.zeros((num_shifts, num_iterations, 3), dtype=np.uint8)
		for t in range(num_iterations):
			for i in range(num_shifts):
				d = all_values[t][i]
				if d.fixed and d.value == 1:
					A[i, t, :] = [255, 0, 0]
				elif d.fixed and d.value == 0:
					A[i, t, :] = [0, 0, 255]
				else:
					value = float(d.value) * 255
					A[i, t, :] = [value, value, value]

				All[i, t] += d.value
		matrices[name] = A

	matrices["all"] = All
	return results, matrices


def plot_fractional_solution(Aorg: np.ndarray, ax) -> None:
	A = Aorg.copy()
	for i in range(A.shape[0]):
		if np.all(A[i, -1, :] == 0):
			A[i, -1, :] = [0, 0, 255]
	lastA = A[:, -1, :].reshape((A.shape[0], 1, A.shape[2]))
	num_end_pixels = int(0.03 * A.shape[1])
	A = np.concatenate((A, lastA[:, num_end_pixels * [0], :]), axis=1)

	ax.imshow(A, interpolation='none', aspect='auto')
	ax.set_xlabel("Iteration number")
	ax.set_ylabel("Working on shift")
	ax.set_xlim([0, ax.get_xlim()[1]])


def plot_data(results: ColgenResults, ax1, name="Colgen") -> None:
	ax2 = ax1.twinx()

	time = np.array(results.time)
	objective = np.array(results.objective)

	ax1.plot(time, objective, "-")
	ax1.set_xlabel("Time (s)")
	ax1.set_ylabel("Objective value")
	ax1.set_title(name)

	ax1.set_xlim([0, 1.05 * np.max(time)])
	ax2.set_xlim([0, 1.05 * np.max(time)])
	ax1.set_ylim([0.9 * np.min(objective), 2.0 * np.min(objective)])

	if results.int_objective is not None:
		ax1.plot(time, results.int_objective, "o", color="#aaaaaa")
		m = np.min(results.int_objective)
		ax1.plot([time[0], time[-1]], [m, m], 'k:')

	for event in results.events:
		event.plot(ax1)

	results.plot_generated_columns_bars(ax2)

	ax1.set_zorder(ax2.get_zorder() + 1)  # put ax in front of ax2
	ax1.patch.set_visible(False)  # hide the 'canvas'


def set_ylim(results_list: List[ColgenResults], ax) -> None:
	top = -1e100
	bottom = 1e100
	for results in results_list:
		bottom = min(bottom, 0.9 * min(results.objective))
		top = max(top, 2.0 * min(results.objective))
	ax.set_ylim([bottom, top])
