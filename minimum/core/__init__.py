import os
from importlib.machinery import SourceFileLoader


def load_module(module_name, filename):
	for prefix in [""] + os.environ["PATH"].split(os.pathsep):
		for prefix2 in ["", ".."]:
			try:
				full_filename = os.path.join(prefix, prefix2, filename)
				return SourceFileLoader(module_name,
				                        full_filename).load_module()
			except FileNotFoundError:
				pass
	raise FileNotFoundError(filename + " not found.")
