#include <iterator>
#include <string>

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/core/scope_guard.h>
#include <minimum/core/string.h>
#include <minimum/linear/easy-ip-internal.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/minisatp.h>
#include <minimum/linear/sat_solver.h>

// Include this header last, since it contains simple names without
// namespaces.
#include <third-party/minisatp/PbSolver.h>

using minimum::core::check;
using minimum::core::range;

namespace minimum {
namespace linear {

class MinisatpSolutions : public Solutions {
   public:
	MinisatpSolutions(IP* ip_) : ip(ip_) {}

	virtual bool get() override {
		check(first_time, "Minisatp can only be called once.");
		first_time = false;
		return solve_minisatp();
	}

   private:
	bool solve_minisatp();
	bool solve_minisatp_single();

	IP* ip;
	bool first_time = true;
};

SolutionsPointer MinisatpSolver::solutions(IP* ip) const {
	return {std::make_unique<MinisatpSolutions>(ip)};
}

bool MinisatpSolutions::solve_minisatp() {
	bool no_cost = true;
	for (size_t j = 0; j < ip->get().variable_size(); ++j) {
		if (ip->get().variable(j).cost() != 0) {
			no_cost = false;
		}
	}
	check(no_cost, "Cost function currently not supported.");
	return solve_minisatp_single();

#if 0
	int upper = 0;
	int lower = 0;

	// Add constraint objective <= 0.
	// 0 here is arbitrary. It will be changed below.
	auto new_constraint = ip->get().add_constraint();
	new_constraint->mutable_bound()->set_lower(-1e100);
	new_constraint->mutable_bound()->set_upper(0);
	at_scope_exit({ ip->get().mutable_constraint()->RemoveLast(); });

	for (size_t j = 0; j < ip->get().variable_size(); ++j) {
		if (ip->get().variable(j).cost() != 0 /* && !allow_ignoring_cost_function */) {
			int icost = ip->get().variable(j).cost();
			check(icost == ip->get().variable(j).cost(), "SAT requires integer costs.");
			auto entry = new_constraint->add_sum();
			entry->set_variable(j);
			entry->set_coefficient(icost);
			if (icost < 0) {
				lower += icost;
			} else {
				upper += icost;
			}
			no_cost = false;
		}
	}

	if (no_cost) {
		return solve_minisatp_single();
	}

	// Minisatp does a linear search by default. This code
	// performs a binary search for the optimum.

	std::vector<double> best_solution;
	do {
		int current = lower + upper;
		if (current < 0) {
			current = current / 2 - 1;
		} else {
			current = current / 2;
		}
		std::clog << "Objective value in [" << lower + ip->get().objective_constant() << ", "
		          << upper + ip->get().objective_constant() << "]." << std::endl;
		if (lower >= upper) {
			break;
		}
		std::clog << "-- Trying " << current + ip->get().objective_constant() << "... ";

		new_constraint->mutable_bound()->set_upper(current);
		bool ok = solve_minisatp_single();

		if (ok) {
			upper = current;
			best_solution.resize(ip->get().variable_size());
			for (auto j : range(ip->get().variable_size())) {
				best_solution[j] = ip->get().variable(j).solution();
			}
			std::clog << "SAT." << std::endl;
		} else {
			lower = current + 1;
			std::clog << "UNSAT." << std::endl;
		}
	} while (true);

	if (!best_solution.empty()) {
		for (auto j : range(ip->get().variable_size())) {
			ip->set_solution(j, best_solution[j]);
		}
	}
	return !best_solution.empty();
#endif
}

bool MinisatpSolutions::solve_minisatp_single() {
	reset_minisatp_globals();
	opt_verbosity = 0;  // Silent mode.
	PbSolver minisatp_solver(true);
	vector<Minisat::Lit> minisatp_literals;

	bool cost_function_added = false;
	for (size_t j = 0; j < ip->get().variable_size(); ++j) {
		auto lb = ip->get().variable(j).bound().lower();
		auto ub = ip->get().variable(j).bound().upper();
		check((lb == 0 || lb == 1) && (ub == 0 || ub == 1),
		      "SAT solver requires boolean variables.");
		minimum_core_assert(lb <= ub);

		std::string varname = "x" + minimum::core::to_string(j);
		minisatp_literals.push_back(Minisat::mkLit(minisatp_solver.getVar(varname.c_str())));
		vec<Minisat::Lit> lvec;
		lvec.push(minisatp_literals.back());
		if (lb == 1) {
			vec<Int> cvec;
			cvec.push(1);
			minisatp_solver.addConstr(lvec, cvec, 1, 1);
		}
		if (ub == 0) {
			vec<Int> cvec;
			cvec.push(1);
			minisatp_solver.addConstr(lvec, cvec, 0, -1);
		}
	}

	PbSolver::solve_Command solve_command = PbSolver::sc_Minimize;
	if (cost_function_added) {
		solve_command = PbSolver::sc_FirstSolution;
	}
	minimum_core_assert(minisatp_literals.size() == ip->get().variable_size());

	const auto num_constraints = ip->get().constraint_size();
	vector<vec<Minisat::Lit>> lit_rows(num_constraints);
	vector<vec<Int>> coeff_rows(num_constraints);
	for (size_t i = 0; i < num_constraints; ++i) {
		auto& constraint = ip->get().constraint(i);
		for (auto& entry : constraint.sum()) {
			auto var = minisatp_literals.at(entry.variable());
			int coeff = entry.coefficient();
			check(coeff == entry.coefficient(),
			      "SAT solver requires integer coefficients in constraints.");
			lit_rows.at(i).push(var);
			coeff_rows.at(i).push(coeff);
		}
	}

	for (size_t i = 0; i < num_constraints; ++i) {
		const int limit = 1000 * 1000 * 1000;
		auto to_int = [limit](double rhs) {
			if (rhs > limit) {
				return limit;
			} else if (rhs < -limit) {
				return -limit;
			} else {
				int irhs = rhs;
				minimum_core_assert(rhs == irhs);
				return irhs;
			}
		};

		auto lower = to_int(ip->get().constraint(i).bound().lower());
		auto upper = to_int(ip->get().constraint(i).bound().upper());

		minimum_core_assert(lit_rows.at(i).size() == coeff_rows.at(i).size());
		/*
		std::clog << "Adding: " << lower << " <= ";
		for (int x = 0; x < lit_rows.at(i).size(); ++x) {
		std::clog << coeff_rows.at(i)[x].data << "*x" << lit_rows.at(i)[x].x << " + ";
		}
		std::clog << " <= " << upper << "\n";
		*/

		if (lower == upper) {
			minisatp_solver.addConstr(lit_rows.at(i), coeff_rows.at(i), lower, 0);
		} else {
			if (lower > -limit) {
				minisatp_solver.addConstr(lit_rows.at(i), coeff_rows.at(i), lower, 1);
			}
			if (upper < limit) {
				minisatp_solver.addConstr(lit_rows.at(i), coeff_rows.at(i), upper, -1);
			}
		}
	}

	minisatp_solver.solve(solve_command);

	// std::cerr << "Petter: " << minisatp_solver->best_goalvalue.data << "\n";
	if (minisatp_solver.best_goalvalue == Int_MAX) {
		return false;
	}

	minimum_core_assert(minisatp_solver.best_model.size() == ip->get().variable_size());
	for (size_t j = 0; j < ip->get().variable_size(); ++j) {
		ip->set_solution(j, minisatp_solver.best_model[j]);
	}

	return true;
}
}  // namespace linear
}  // namespace minimum
