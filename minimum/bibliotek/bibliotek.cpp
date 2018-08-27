
#include <algorithm>
#include <iostream>

#include <minimum/bibliotek/bibliotek.h>
#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/core/time.h>
#include <minimum/linear/constraint_solver.h>
#include <minimum/linear/glpk.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/sat_solver.h>
#include <minimum/linear/solver.h>
using minimum::core::check;
using minimum::core::range;

static_assert(sizeof(int) == 4, "int needs to be 32 bits.");
static_assert(sizeof(double) == 8, "double needs to be 32 bits.");

const char* error_string = "";
std::string error_string_data;

namespace {
void set_error_string(const std::exception& error) {
	error_string_data = error.what();
	error_string = error_string_data.c_str();
}
}  // namespace

int minimum_linear_test(int* data, int n) {
	int sum = 0;
	for (int i : range(n)) {
		sum += data[i];
	}
	return sum;
}

void* minimum_linear_ip(const char* data, int size) { return new minimum::linear::IP(data, size); }

void minimum_linear_ip_destroy(void* ip_) {
	auto ip = reinterpret_cast<minimum::linear::IP*>(ip_);
	delete ip;
}

void minimum_linear_ip_save_mps(void* ip_, const char* file_name) {
	minimum_core_assert(ip_ != nullptr);
	auto ip = reinterpret_cast<minimum::linear::IP*>(ip_);
	minimum::linear::IPSolver solver;
	solver.save_MPS(*ip, file_name);
}

void* minimum_linear_solver_default() {
	auto solver = new minimum::linear::IPSolver;
	solver->set_silent(true);
	return solver;
}

void* minimum_linear_solver_glpk() {
	auto solver = new minimum::linear::GlpkSolver();
	solver->set_silent(true);
	return solver;
}

void* minimum_linear_solver_minisat() {
	auto solver =
	    new minimum::linear::IpToSatSolver(std::bind(minimum::linear::minisat_solver, true));
	solver->set_silent(true);
	return solver;
}

void* minimum_linear_solver_constraint() {
	auto solver = new minimum::linear::ConstraintSolver();
	solver->set_silent(true);
	return solver;
}

void minimum_linear_solver_destroy(void* solver_) {
	auto solver = reinterpret_cast<minimum::linear::Solver*>(solver_);
	delete solver;
}

void* minimum_linear_solutions(void* solver_, void* ip_) {
	try {
		auto ip = reinterpret_cast<minimum::linear::IP*>(ip_);
		auto solver = reinterpret_cast<minimum::linear::Solver*>(solver_);
		auto solutions = new minimum::linear::SolutionsPointer(solver->solutions(ip));
		return solutions;
	} catch (std::exception& err) {
		set_error_string(err);
		return nullptr;
	}
}

int minimum_linear_solutions_get(void* solutions_, void* ip_, double* solution, int len) {
	auto solutions = reinterpret_cast<minimum::linear::SolutionsPointer*>(solutions_);
	auto ip = reinterpret_cast<minimum::linear::IP*>(ip_);
	try {
		if (!(*solutions)->get()) {
			return 1;
		}
	} catch (std::exception& err) {
		set_error_string(err);
		return 2;
	} catch (...) {
		return 20;
	}

	if (len < ip->get().variable_size()) {
		std::cerr << len << "<" << ip->get().variable_size() << "\n";
		return 3;
	}

	for (auto j : range(ip->get().variable_size())) {
		solution[j] = ip->get_solution().primal(j);
	}
	return 0;
}

void minimum_linear_solutions_destroy(void* solutions_) {
	auto solutions = reinterpret_cast<minimum::linear::SolutionsPointer*>(solutions_);
	delete solutions;
}

int add_function(int a, int b) { return a + b + 2; }

int sum_function(int* d, int n) {
	std::cout << "n=" << n << "\n";
	int sum = 0;
	for (int i : range(n)) {
		std::cout << "d[" << i << "]=" << d[i] << "\n";
		sum += d[i];
	}
	std::cout << "Sum is " << sum << "\n";
	return sum;
}
