#pragma once

#include <memory>
#include <vector>

#include <gflags/gflags.h>

#include <minimum/core/string.h>
#include <minimum/linear/first_order_solver.h>
#include <minimum/linear/glpk.h>
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/scs.h>
#include <minimum/linear/solver.h>

std::vector<std::string> allowed = {"ip", "glpk", "minisat", "scs", "primal-dual"};

static bool validate_solver(const char* flagname, const std::string& value) {
	for (auto& allow : allowed) {
		if (value == allow) {
			return true;
		}
	}
	if (value == "") {
		return false;
	}
	std::cerr << value << " is an invalid value for \'solver.\' It needs to be one of "
	          << minimum::core::to_string(allowed) << "." << std::endl;
	return false;
}

DEFINE_string(solver, "ip", "The solver to use.");
DEFINE_validator(solver, validate_solver);

namespace minimum {
namespace linear {

// Creates a solver for a given IP. The solver used is determined via
// the command line option --solver.
std::unique_ptr<Solver> create_solver_from_command_line() {
	if (FLAGS_solver == "ip") {
		return std::make_unique<IPSolver>();
	} else if (FLAGS_solver == "glpk") {
		return std::make_unique<GlpkSolver>();
	} else if (FLAGS_solver == "minisat") {
		return std::make_unique<IpToSatSolver>(bind(minisat_solver, false));
	} else if (FLAGS_solver == "scs") {
		return std::make_unique<ScsSolver>();
	} else if (FLAGS_solver == "primal-dual") {
		return std::make_unique<PrimalDualSolver>();
	} else {
		minimum::core::check(false, "Unknown solver.");
		return nullptr;
	}
}
}  // namespace linear
}  // namespace minimum
