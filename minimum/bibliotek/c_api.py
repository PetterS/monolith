import ctypes  # type: ignore
import numpy as np  # type: ignore
import os


def load_library(name):
	for prefix in [""] + os.environ["PATH"].split(os.pathsep):
		try:
			# print(os.path.join(prefix, name))
			return ctypes.cdll.LoadLibrary(os.path.join(prefix, name))
		except OSError:
			# print("    ", err)
			pass

		try:
			# print(os.path.join(prefix, "lib" + name + ".so"))
			return ctypes.cdll.LoadLibrary(
			    os.path.join(prefix, "lib" + name + ".so"))
		except OSError:
			# print("    ", err)
			pass
	raise OSError(name + " not found.")


minimum_bibliotek = load_library("minimum_bibliotek")

minimum_bibliotek_error_string_ptr = ctypes.c_char_p.in_dll(
    minimum_bibliotek, "error_string")


def error_string() -> str:
	val = minimum_bibliotek_error_string_ptr.value
	if val:
		return val.decode("utf-8")
	else:
		return ""


ctypes_minimum_linear_test = minimum_bibliotek.minimum_linear_test
ctypes_minimum_linear_test.restype = ctypes.c_int


def minimum_linear_test(data):
	assert data.dtype == np.int32
	return ctypes_minimum_linear_test(
	    data.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
	    ctypes.c_int(data.size))


class Ip_p(ctypes.c_void_p):
	# http://stackoverflow.com/questions/17840144/why-does-setting-ctypes-dll-function-restype-c-void-p-return-long
	pass


class Solver_p(ctypes.c_void_p):
	pass


class Solutions_p(ctypes.c_void_p):
	pass


ctypes_minimum_linear_ip = minimum_bibliotek.minimum_linear_ip
ctypes_minimum_linear_ip.restype = Ip_p


def minimum_linear_ip(proto_byte_string: bytes) -> Ip_p:
	return ctypes_minimum_linear_ip(
	    ctypes.c_char_p(proto_byte_string),
	    ctypes.c_int32(len(proto_byte_string)))


ctypes_minimum_linear_ip_destroy = minimum_bibliotek.minimum_linear_ip_destroy
ctypes_minimum_linear_ip_destroy.restype = None


def minimum_linear_ip_destroy(ip: Ip_p):
	assert type(ip) == Ip_p
	return ctypes_minimum_linear_ip_destroy(ip)


ctypes_minimum_linear_ip_save_mps = minimum_bibliotek.minimum_linear_ip_save_mps
ctypes_minimum_linear_ip_save_mps.restype = None


def minimum_linear_ip_save_mps(ip: Ip_p, file_name: str):
	assert type(ip) == Ip_p
	assert type(file_name) == str
	ctypes_minimum_linear_ip_save_mps(
	    ip, ctypes.c_char_p(file_name.encode("utf-8")))


ctypes_minimum_linear_solver_default = minimum_bibliotek.minimum_linear_solver_default
ctypes_minimum_linear_solver_default.restype = Solver_p


def minimum_linear_solver_default():
	return ctypes_minimum_linear_solver_default()


ctypes_minimum_linear_solver_minisat = minimum_bibliotek.minimum_linear_solver_minisat
ctypes_minimum_linear_solver_minisat.restype = Solver_p


def minimum_linear_solver_minisat():
	return ctypes_minimum_linear_solver_minisat()


ctypes_minimum_linear_solver_destroy = minimum_bibliotek.minimum_linear_solver_destroy
ctypes_minimum_linear_solver_destroy.restype = None


def minimum_linear_solver_destroy(solver: Solver_p):
	assert type(solver) == Solver_p
	return ctypes_minimum_linear_solver_destroy(solver)


ctypes_minimum_linear_solutions = minimum_bibliotek.minimum_linear_solutions
ctypes_minimum_linear_solutions.restype = Solutions_p


def minimum_linear_solutions(solver: Solver_p, ip: Ip_p):
	assert type(ip) == Ip_p
	assert type(solver) == Solver_p
	solutions = ctypes_minimum_linear_solutions(solver, ip)
	assert solutions, error_string()
	return solutions


ctypes_minimum_linear_solutions_get = minimum_bibliotek.minimum_linear_solutions_get
ctypes_minimum_linear_solutions_get.restype = ctypes.c_int


def minimum_linear_solutions_get(solutions: Solutions_p, ip: Ip_p, solution):
	assert type(solutions) == Solutions_p
	assert type(ip) == Ip_p
	assert solution.dtype == np.float64
	return ctypes_minimum_linear_solutions_get(
	    solutions, ip, solution.ctypes.data_as(
	        ctypes.POINTER(ctypes.c_double)), ctypes.c_int32(solution.size))


ctypes_minimum_linear_solutions_destroy = minimum_bibliotek.minimum_linear_solutions_destroy
ctypes_minimum_linear_solutions_destroy.restype = None


def minimum_linear_solutions_destroy(solutions: Solutions_p):
	assert type(solutions) == Solutions_p
	return ctypes_minimum_linear_solutions_destroy(solutions)


ctypes_sum_function = minimum_bibliotek.sum_function


def sum_function(data):
	return ctypes_sum_function(
	    data.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
	    ctypes.c_int32(data.size))


ctypes_add_function = minimum_bibliotek.add_function


def add_function(a, b):
	return ctypes_add_function(ctypes.c_int32(a), ctypes.c_int32(b))


# Run a quick test to make sure everything is OK.
_a = np.array([1, 2, 3, 4, 5, 6], dtype=np.int32)
assert minimum_linear_test(_a) == 21
