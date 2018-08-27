import argparse
import os
import re
import subprocess
from typing import List, Union

import minimum.linear.data.util

parser = argparse.ArgumentParser(description='Run colgen scheduling tests.')
parser.add_argument(
    '--more',
    default=False,
    action='store_const',
    const=True,
    help='Run more tests')
args = parser.parse_args()

gent = "Gent"
tests_to_run: List[Union[str, int]] = [i for i in range(1, 11)]
if args.more:
	tests_to_run = [i for i in range(1, 20)]
tests_to_run += [gent]

dir = minimum.linear.data.util.get_directory()

for i in tests_to_run:
	file_name = "Instance" + str(i) + ".txt"
	print("Running", file_name, "...")

	name = "Instance" + str(i)
	full_name = os.path.join(dir, "shift_scheduling", name + ".txt")
	output_file_name = os.path.join(
	    os.path.dirname(__file__), file_name + ".out")
	error_file_name = os.path.join(
	    os.path.dirname(__file__), file_name + ".err")
	org_error_file_name = os.path.join(
	    os.path.dirname(__file__), file_name + "-org.err")

	start_args = [
	    "shift_scheduling_colgen", "--nosave_solution", "--num_solutions=2",
	    f"--proto_log_file={name}.record", full_name
	]
	if i == gent:
		start_args += ["--read_gent_instance"]

	completed = subprocess.run(
	    start_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout = completed.stdout
	stderr = completed.stderr

	with open(org_error_file_name, "wb") as f:
		f.write(stderr)

	stderr = re.sub(
	    br"\[ WAIT \].*  OK  \].*$", b"[  OK  ]", stderr, flags=re.MULTILINE)

	stderr = re.sub(
	    br"^( *\d+ *\d+\.?\d* *\d+ *\d+ *-?\d+ *\d+ *\d+ *\d+).+( *\d+%)",
	    br"\1        n/a       n/a     n/a     n/a          \2",
	    stderr,
	    flags=re.MULTILINE)

	stderr = re.sub(
	    br"^Elapsed time : \d+\.?\d*s\.",
	    b"Elapsed time : n/a.",
	    stderr,
	    flags=re.MULTILINE)

	stderr = stderr.replace(b"\r\n", b"\n")
	stdout = stdout.replace(b"\r\n", b"\n")

	with open(output_file_name, "wb") as f:
		f.write(stdout)
	with open(error_file_name, "wb") as f:
		f.write(stderr)
