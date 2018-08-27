#!/usr/bin/env python3

from glob import glob
import os

# Set this to the location of NSPLib.
nsplib = r"C:\Users\Petter\Dropbox\Datasets\NSPLib"


def run_solver(data_set, case):
	case_file = os.path.join(nsplib, "Cases", str(case) + ".gen")
	log_file = data_set + "." + str(case) + ".output.log"

	files = glob(os.path.join(nsplib, data_set, "*.nsp"))
	names = [f.split(".")[0] for f in files]
	names = [n.split(os.path.sep)[-1] for n in names]
	nums = sorted([int(n) for n in names])
	files = [os.path.join(nsplib, data_set, str(n) + ".nsp") for n in nums]

	try:
		os.unlink(log_file)
	except FileNotFoundError:
		pass

	for f in files:
		print(case_file)
		print(f)
		print(log_file)

		# This may need to change depending on shell.
		os.system("nurses " + f + " " + case_file + " >> " + log_file)


for data_set in ["N25", "N50", "N75", "N100"]:
	for case in [1, 2, 3, 4, 5, 6, 7, 8]:
		run_solver(data_set, case)

for data_set in ["N30", "N60"]:
	for case in [9, 10, 11, 12, 13, 14, 15, 16]:
		run_solver(data_set, case)
