#!/usr/bin/env python3
"""
	Pre-commit hook for Python, C++, etc.

If run like

	python3 precommit.py install

in Linux or Windows it will install itself as a git hook.

When run, it check staged files for style, lint errors, type errors
etc. It will also run Python unit tests if any Python tests are
modified.

It requires the following Python packages:
	mypy pyflakes pytest yapf
If requires the following LLVM programs:
	clang-format
It requires the following Node packages:
	eslint prettier
"""
import os
import subprocess
import stat
import sys
import time
from typing import List

import colorama
from colorama import Fore, Style

SOURCE_TO_CHECK = ["minimum", "misc"]

colorama.init()


class PrecommitFailure(Exception):
	"""Means that a precommit failed and the program will terminate."""

	def __init__(self, message: str) -> None:
		self.message = message


class PrecommitWarning(Exception):
	"""Issues a warning and the program will continue."""

	def __init__(self, message: str) -> None:
		self.message = message


class Precommit:
	"""Base class that describes a precommit."""

	def run(self) -> None:
		"""Runs the precommit and raises PrecommitFailure on failure."""
		pass

	def check_for_modifications(self, message: str) -> None:
		"""Fails the test with an error if there are local modifications."""
		run = subprocess.run(["git", "ls-files", "-m"], stdout=subprocess.PIPE)
		assert run.returncode == 0
		if len(run.stdout) > 0:
			raise PrecommitFailure(message)

	def run_with_check(self, args: List[str], errormessage: str) -> None:
		"""Runs the list of arguments and exits the current process on failure."""
		run = subprocess.run(
		    args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		if run.returncode != 0:
			raise PrecommitFailure(
			    run.stdout.decode("utf-8") + "\n\n" + errormessage)


class NoModifiedFiles(Precommit):
	"""Checks that no files are modified.

	The presubmit script assumes that files are either staged or unmodified.
	"""

	def run(self):
		self.check_for_modifications(
		    "There are local modifications. Please stage/stash them.")


class LineEndings(Precommit):
	"""Checks that all files do not contain CRLF line endings."""

	def __init__(self, files):
		self.files = files

	def run(self):
		failure_message = ""
		for name in self.files:
			with open(name, "rb") as f:
				data = f.read()
			if b"\r\n" in data:
				data = data.replace(b'\r\n', b'\n')
				with open(name, "wb") as f:
					f.write(data)
				failure_message += f"{name} contains CRLF line endings.\n"
		if failure_message:
			raise PrecommitFailure(failure_message +
			                       "\nThe files have been fixed." +
			                       " Please stage them if OK.")


class ValidUtf8(Precommit):
	"""Checks that all files decode as UTF-8."""

	def __init__(self, files):
		self.files = files

	def run(self):
		for name in self.files:
			try:
				with open(name, "rb") as f:
					f.read().decode("utf-8")
			except UnicodeDecodeError:
				raise PrecommitFailure(f"{name} is not valid UTF-8.")


class UnixLineEndings(Precommit):
	"""Checks that all files decode as UTF-8."""

	def __init__(self, files):
		self.files = files

	def run(self):
		for name in self.files:
			with open(name, "rb") as f:
				data = f.read()
			unix_data = data.replace(b"\r\n", b"\n")
			if unix_data != data:
				with open(name, "rb") as f:
					f.write(unix_data)

		self.check_for_modifications("Line endings have been changed to UNIX. "
		                             + "Please stage them if OK.")


class Pyflakes(Precommit):
	"""Runs Pyflakes.

	Pyflakes is very fast and has essentially zero false positives.
	"""

	def __init__(self, files):
		self.files = files

	def run(self):
		self.run_with_check([sys.executable, "-m", "pyflakes"] + self.files,
		                    "Pyflakes failed.")


class MyPy(Precommit):
	"""Runs the type-checker Mypy on the provided Python files."""

	def __init__(self, files):
		self.files = files

	def run(self):
		self.run_with_check(
		    [sys.executable, "-m", "mypy", "--ignore-missing-imports"
		     ] + self.files, "Mypy failed.")


class PyTest(Precommit):
	"""Runs the Python unit tests in the code base."""

	def run(self):
		self.run_with_check(
		    [sys.executable, "-m", "pytest", "-v"] + SOURCE_TO_CHECK,
		    "Pytest failed.")


class Yapf(Precommit):
	"""Checks that all provided Python files are unmodifed by the formatter."""

	def __init__(self, files):
		self.files = files

	def run(self):
		self.run_with_check([
		    sys.executable, "-m", "yapf", "--recursive", "--in-place",
		    "--parallel", "-vv"
		] + self.files, "Yapf failed.")
		# If Yapf made any changes, we fail the precommit to let the user review
		# them.
		self.check_for_modifications("Python formatter made modifications. " +
		                             "Please stage them if OK.")


class Prettier(Precommit):
	def __init__(self, files):
		self.files = files

	def run(self):
		subprocess.run(
		    ["prettier", "--write", "--loglevel", "log"] + self.files,
		    check=True,
		    shell=True)
		self.check_for_modifications(
		    "Javascript formatter made modifications.")


class ClangFormat(Precommit):
	"""Checks that all provided C++ files are unmodifed by the formatter."""

	def __init__(self, files):
		self.files = files

	def run(self):
		try:
			self.run_with_check(
			    ["clang-format", "-i", "-style=file"] + self.files,
			    "Clang-format failed.")
		except FileNotFoundError:
			raise PrecommitWarning("clang-format not found.")
		# If Clang made any changes, we fail the precommit to let the user review
		# them.
		self.check_for_modifications("C++ formatter made modifications. " +
		                             "Please stage them if OK.")


def get_staged_files(pattern: str) -> List[str]:
	"""Returns all staged files about to be committed."""
	run = subprocess.run(
	    [
	        "git", "diff", "--cached", "--name-only", "--diff-filter=ACM",
	        pattern
	    ],
	    stdout=subprocess.PIPE)
	assert run.returncode == 0
	return [
	    file for file in run.stdout.decode("utf-8").split("\n")
	    if len(file) > 0 and not "third-party/" in file
	]


def print_files(name: str, files: List[str]) -> None:
	if not files:
		return
	print(Style.BRIGHT, "The following", Fore.YELLOW, name, Fore.WHITE,
	      "files will be checked:")
	for f in files:
		print(Fore.YELLOW, "   ", f)
	print(Style.RESET_ALL)


def install_if_needed() -> None:
	precommit_hook = os.path.join(
	    os.path.dirname(__file__), ".git", "hooks", "pre-commit")
	if os.path.isfile(precommit_hook):
		return

	interpreter = "python3"
	if os.name == "nt":
		interpreter = "py -3"
	with open(precommit_hook, "w") as f:
		print("#!/bin/bash", file=f)
		print("set -e", file=f)
		print("{} precommit.py".format(interpreter), file=f)

	# Make the script executable. Required on Linux.
	st = os.stat(precommit_hook)
	os.chmod(precommit_hook, st.st_mode | stat.S_IEXEC)

	print("Installed pre-commit hook in", precommit_hook)


def main():
	if not os.path.isdir(".git"):
		print("Need to run in the git root directory.")
		sys.exit(1)
	install_if_needed()

	all_files = get_staged_files("*")
	if not all_files:
		print("No files to commit. OK.")
		return

	python_files = get_staged_files("*.py")
	cpp_files = get_staged_files("*.cpp") + get_staged_files("*.h")
	js_files = get_staged_files("*.js")
	txt_files = get_staged_files("*.txt")
	cmake_files = get_staged_files("*.cmake")
	all_text_files = (
	    python_files + cpp_files + js_files + txt_files + cmake_files)

	precommits: List[Precommit] = [
	    NoModifiedFiles(),
	    LineEndings(all_text_files),
	    ValidUtf8(all_text_files),
	    UnixLineEndings(all_text_files),
	]
	if python_files:
		precommits += [
		    Pyflakes(python_files),
		    MyPy(python_files),
		    PyTest(),
		    Yapf(python_files)
		]
	if cpp_files:
		precommits += [ClangFormat(cpp_files)]
	if js_files:
		precommits += [Prettier(js_files)]

	print(Style.BRIGHT, "Pre-commit run starting.")
	print_files("Python", python_files)
	print_files("C++", cpp_files)
	print_files("Javascript", js_files)

	for precommit in precommits:
		try:
			print(precommit.__class__.__name__, "...", end="")
			sys.stdout.flush()
			start_time = time.time()
			precommit.run()
			elapsed = time.time() - start_time
			print(Fore.GREEN, Style.BRIGHT, "OK", Style.RESET_ALL,
			      "{:.1f}s.".format(elapsed))
		except PrecommitWarning as warning:
			print(Fore.YELLOW, Style.BRIGHT, "WARNING", Style.RESET_ALL,
			      warning.message)
		except PrecommitFailure as failure:
			print(Fore.RED, Style.BRIGHT, "FAIL")
			print(failure.message, Style.RESET_ALL)
			sys.exit(1)


if __name__ == "__main__":
	main()
