//
// Converts an integer program to a SAT problem.
//

#include <iomanip>
#include <iterator>
#include <string>

#include <minimum/core/check.h>
#include <minimum/core/scope_guard.h>
#include <minimum/core/string.h>
#include <minimum/core/time.h>
#include <minimum/linear/easy-ip-internal.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/sat_solver.h>

using minimum::core::check;

namespace minimum {
namespace linear {

class SatSolutions : public Solutions {
   public:
	SatSolutions(const IpToSatSolver* creator_, IP* ip_, std::unique_ptr<SatSolver> sat_solver_)
	    : ip(ip_), solver(std::move(sat_solver_)), creator(creator_) {
		convert_to_sat();
	}

	virtual bool get() {
		auto result = next(first_time);
		first_time = false;
		return result;
	}

	void convert_to_sat();
	bool next(bool start_over);

	IP* ip = nullptr;
	std::unique_ptr<SatSolver> solver;

	vector<SatSolver::Var> literals;
	vector<SatSolver::Literal> objective_function_literals;
	vector<SatSolver::Literal> objective_function_slack_literals;
	int sat_objective_offset = 0;

	const IpToSatSolver* creator = nullptr;
	bool first_time = true;
};

IpToSatSolver::IpToSatSolver(std::function<std::unique_ptr<SatSolver>()> sat_solver_factory_)
    : sat_solver_factory(std::move(sat_solver_factory_)) {}

IpToSatSolver::~IpToSatSolver() {}

SolutionsPointer IpToSatSolver::solutions(IP* ip) const {
	ip->linearize_pseudoboolean_terms();
	return {std::make_unique<SatSolutions>(this, ip, sat_solver_factory())};
}

std::unique_ptr<SatSolver> IpToSatSolver::convert(IP* ip) {
	SatSolutions solutions(this, ip, sat_solver_factory());
	auto sat_solver = move(solutions.solver);
	return sat_solver;
}

std::pair<SatSolver::Var, bool> negate(std::pair<SatSolver::Var, bool> literal) {
	return {literal.first, !literal.second};
}

// The binomial approach is only effective for small k.
void add_at_most_k_constraint_binomial(SatSolver* solver,
                                       const vector<SatSolver::Literal>& literals,
                                       int k) {
	vector<int> index_set;
	for (size_t ix = 0; ix < literals.size(); ++ix) {
		index_set.emplace_back(ix);
	}

	vector<vector<int>> subsets;
	generate_subsets(index_set, k + 1, &subsets);

	vector<SatSolver::Literal> clause(k + 1);
	for (auto& subset : subsets) {
		int clause_index = 0;
		for (int ix : subset) {
			clause[clause_index++] = negate(literals[ix]);
		}
		solver->add_clause(clause);
	}
}

void add_at_most_k_constraint(SatSolver* solver, const vector<SatSolver::Literal>& X, int k) {
	minimum_core_assert(k >= 0);
	auto n = X.size();

	// Perhaps the binomial is better when e.g. n <= 7?
	if (n <= 1 || k == 0) {
		add_at_most_k_constraint_binomial(solver, X, k);
		return;
	}

	// This implementation follows
	// Carsten Sinz,
	// “Towards an Optimal CNF Encoding of Boolean Cardinality Constraints,”
	// Principles and Practice of Constraint Programming, 2005.
	// §2, page 2.

	vector<vector<SatSolver::Literal>> s;
	for (size_t i = 0; i < n; ++i) {
		s.emplace_back();
		for (size_t j = 0; j < k; ++j) {
			s.back().emplace_back(solver->add_helper_variable(), true);
		}
	}

	// (¬x1 ∨ s1,1)
	solver->add_clause(negate(X[0]), s[0][0]);

	// for 1 < j ≤ k
	for (size_t j = 1; j < k; ++j) {
		// (¬s1,j)
		solver->add_clause(negate(s[0][j]));
	}

	// for 1 < i < n
	for (size_t i = 1; i < n - 1; ++i) {
		// (¬xi ∨ si,1)
		solver->add_clause(negate(X[i]), s[i][0]);
		// (¬si−1,1 ∨ si,1)
		solver->add_clause(negate(s[i - 1][0]), s[i][0]);

		// for 1 < j ≤ k
		for (size_t j = 1; j < k; ++j) {
			// (¬xi ∨ ¬si−1,j−1 ∨ si,j)
			solver->add_clause({negate(X[i]), negate(s[i - 1][j - 1]), s[i][j]});
			// (¬si−1,j ∨ si,j)
			solver->add_clause(negate(s[i - 1][j]), s[i][j]);
		}
		// (¬xi ∨ ¬si−1,k)
		solver->add_clause(negate(X[i]), negate(s[i - 1][k - 1]));
	}
	// (¬xn ∨ ¬sn−1,k)
	solver->add_clause(negate(X[n - 1]), negate(s[n - 2][k - 1]));
}

void add_1_or_fewer_constraint(SatSolver* solver, const vector<SatSolver::Literal>& X) {
	for (int j1 = 0; j1 < X.size(); ++j1) {
		for (int j2 = j1 + 1; j2 < X.size(); ++j2) {
			solver->add_clause(negate(X[j1]), negate(X[j2]));
		}
	}
}

void add_equal_to_1_constraint(SatSolver* solver, const vector<SatSolver::Literal>& X) {
	add_1_or_fewer_constraint(solver, X);
	solver->add_clause(X);
}

void SatSolutions::convert_to_sat() {
	using namespace std;

	sat_objective_offset = 0;
	literals.clear();
	objective_function_literals.clear();
	for (auto& var : ip->get().variable()) {
		auto lb = var.bound().lower();
		auto ub = var.bound().upper();
		check((lb == 0 || lb == 1) && (ub == 0 || ub == 1),
		      "SAT solver requires boolean variables. Bound is [",
		      lb,
		      ", ",
		      ub,
		      "].");
		minimum_core_assert(lb <= ub);

		literals.push_back(solver->add_variable());
		if (lb == 1) {
			solver->add_clause(make_pair(literals.back(), true));
		}
		if (ub == 0) {
			solver->add_clause(make_pair(literals.back(), false));
		}

		if (!creator->allow_ignoring_cost_function) {
			if (var.cost() != 0) {
				int icost = var.cost();
				check(icost == var.cost(), "SAT requires integer costs.");

				// Add new literals equivalent to the variable and add them to
				// the vector of cost literals.
				for (int count = 1; count <= std::abs(icost); ++count) {
					auto lit = solver->add_variable();
					solver->add_clause(make_pair(literals.back(), true), make_pair(lit, false));
					solver->add_clause(make_pair(literals.back(), false), make_pair(lit, true));
					bool is_positive = true;
					if (icost < 0) {
						is_positive = false;
						sat_objective_offset -= 1;
					}
					objective_function_literals.emplace_back(lit, is_positive);
				}
			}
		}
	}

	if (objective_function_literals.size() > 0) {
		// An objective function of
		//   3x + y
		// is modelled as
		//   x1 + x2 + x3 + y1
		// where
		//   x1 ⇔ x
		//   x2 ⇔ x
		//   x3 ⇔ x
		//   y1 ⇔ y.
		// Then slack variables are added with an upper bound:
		//   x1 + x2 + x3 + y1 + s1 + s2 + s3 + s4 ≤ 4.
		// By assuming a different number of slack variables = 1, different
		// objective functions value can be tested for satisfiability.

		// First, we need to add the slack literals to the objective function.
		objective_function_slack_literals.clear();
		for (size_t i = 0; i < objective_function_literals.size(); ++i) {
			auto lit = solver->add_variable();
			objective_function_slack_literals.emplace_back(lit, true);
		}

		std::vector<SatSolver::Literal> objective_clause;
		for (auto& lit : objective_function_literals) {
			objective_clause.emplace_back(lit);
		}
		for (auto& lit : objective_function_slack_literals) {
			objective_clause.emplace_back(lit);
		}
		add_at_most_k_constraint(
		    solver.get(), objective_clause, objective_function_literals.size());
	}

	auto num_constraints = ip->get().constraint_size();
	std::vector<int> lower(num_constraints);
	std::vector<int> upper(num_constraints);
	for (size_t i = 0; i < num_constraints; ++i) {
		auto to_int = [](double rhs) {
			const int limit = 1000 * 1000 * 1000;
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

		lower[i] = to_int(ip->get().constraint(i).bound().lower());
		upper[i] = to_int(ip->get().constraint(i).bound().upper());
	}

	vector<vector<SatSolver::Literal>> lit_rows(num_constraints);
	for (size_t i = 0; i < num_constraints; ++i) {
		auto& constraint = ip->get().constraint(i);
		for (auto& entry : constraint.sum()) {
			auto var = literals.at(entry.variable());
			int coeff = entry.coefficient();
			check(coeff == entry.coefficient(),
			      "SAT solver requires integer coefficients in constraints.");

			if (coeff == 1) {
				lit_rows.at(i).emplace_back(var, true);
			} else if (coeff > 1) {
				// Add new literals equivalent to the variable and add them to
				// the vector of cost literals.
				for (int count = 1; count <= coeff; ++count) {
					auto lit = solver->add_helper_variable();
					solver->add_clause(make_pair(var, true), make_pair(lit, false));
					solver->add_clause(make_pair(var, false), make_pair(lit, true));
					lit_rows.at(i).emplace_back(lit, true);
				}
			} else if (coeff == -1) {
				lower[i] += 1;
				upper[i] += 1;
				lit_rows.at(i).emplace_back(var, false);
			} else if (coeff < -1) {
				// Add new literals equivalent to the variable and add them to
				// the vector of cost literals.
				for (int count = 1; count <= -coeff; ++count) {
					auto lit = solver->add_helper_variable();
					solver->add_clause(make_pair(var, true), make_pair(lit, false));
					solver->add_clause(make_pair(var, false), make_pair(lit, true));
					lit_rows.at(i).emplace_back(lit, false);

					lower[i] += 1;
					upper[i] += 1;
				}
			}
		}
	}

	vector<vector<int>> subsets;
	// std::map<tuple<int, int, int>, int> count;
	for (size_t i = 0; i < num_constraints; ++i) {
		int num_literals = lit_rows.at(i).size();

		// count[make_tuple(std::max(-1, lower[i]), min<int>(upper[i], num_literals), num_literals)]
		// += 1;

		// Some common cases are treated specially here for
		// increased performance.
		if (lower[i] == 1 && upper[i] == 1) {
			// x1 + ... + xn == 1.
			add_equal_to_1_constraint(solver.get(), lit_rows[i]);
		} else if (lower[i] == 1 && upper[i] >= num_literals) {
			// x1 + ... + xn >= 1.
			solver->add_clause(lit_rows[i]);
		} else if (lower[i] <= 0 && upper[i] == 1) {
			// x1 + ... + xn <= 1.
			add_1_or_fewer_constraint(solver.get(), lit_rows[i]);
		} else if (lower[i] <= 0 && upper[i] == num_literals - 1) {
			// x1 + ... + xn <= n - 1.
			auto neg_lit_row = lit_rows[i];
			for (auto& lit : neg_lit_row) {
				lit = negate(lit);
			}
			solver->add_clause(neg_lit_row);
		} else {
			if (lower[i] > 0) {
				auto neg_lit_row = lit_rows[i];
				for (auto& lit : neg_lit_row) {
					lit = negate(lit);
				}
				int bound = int(neg_lit_row.size()) - lower[i];
				if (bound < 0) {
					// This constraint is always infeasible.
					auto x = solver->add_helper_variable();
					solver->add_clause(make_pair(x, false));
					solver->add_clause(make_pair(x, true));
				} else {
					add_at_most_k_constraint(solver.get(), neg_lit_row, bound);
				}
			}

			if (upper[i] < num_literals) {
				if (upper[i] < 0) {
					// This constraint is always infeasible.
					auto x = solver->add_helper_variable();
					solver->add_clause(make_pair(x, false));
					solver->add_clause(make_pair(x, true));
				} else {
					add_at_most_k_constraint(solver.get(), lit_rows[i], upper[i]);
				}
			}
		}
	}

	// for (auto& elem : count) {
	//	std::cerr << core::to_string(elem.first) << ": " << elem.second << "\n";
	//}
}

bool SatSolutions::next(bool start_over) {
	minimum_core_assert(solver);

	std::vector<int> sat_solution;
	if (!start_over) {
		// Forbid previous solution.
		vector<SatSolver::Literal> negated_solution;
		for (size_t j = 0; j < ip->get().variable_size(); ++j) {
			if (ip->get().variable(j).is_helper()) {
				continue;
			}

			if (ip->get_solution().primal(j) == 1) {
				negated_solution.push_back(std::make_pair(literals[j], false));
			} else {
				negated_solution.push_back(std::make_pair(literals[j], true));
			}
		}
		solver->add_clause(negated_solution);
	} else {
		// If we have a solution, use it for SAT heuristics.
		if (ip->is_feasible()) {
			for (int i = 0; i < ip->get().variable_size(); ++i) {
				if (ip->get_solution().primal(i) >= 0.5) {
					sat_solution.push_back(true);
				} else {
					sat_solution.push_back(false);
				}
			}
			if (!creator->silent) {
				std::vector<double> solution;
				solution.reserve(ip->get().variable_size());
				for (size_t j = 0; j < ip->get().variable_size(); ++j) {
					solution.push_back(ip->get_solution().primal(j));
				}
				print_histogram_01(solution.begin(), solution.end());
			}
		}
	}

	double start_time = core::wall_time();
	bool result = solver->solve(&sat_solution);
	if (!creator->silent && literals.size() > 30) {
		std::clog << "Time solving SAT: " << core::wall_time() - start_time << "\n";
	}

	if (!result) {
		return false;
	}

	if (start_over && objective_function_literals.size() > 0) {
		int upper = objective_function_slack_literals.size();
		int lower = 0;
		bool ok;

		do {
			int current = (lower + upper) / 2;
			if (!creator->silent) {
				std::clog << "Objective value in ["
				          << lower + sat_objective_offset + ip->get_objective_constant() << ", "
				          << upper + sat_objective_offset + ip->get_objective_constant() << "]."
				          << std::endl;
			}
			if (lower >= upper) {
				break;
			}
			if (!creator->silent) {
				std::clog << "-- Trying "
				          << current + sat_objective_offset + ip->get_objective_constant()
				          << "... ";
			}

			std::vector<SatSolver::Literal> assumptions;
			for (int i = 0; i < objective_function_literals.size() - current; ++i) {
				assumptions.push_back(objective_function_slack_literals.at(i));
			}
			ok = solver->solve(assumptions, &sat_solution);

			if (ok) {
				upper = current;
				if (!creator->silent) {
					std::clog << "SAT.\n";
				}
			} else {
				lower = current + 1;
				if (!creator->silent) {
					std::clog << "UNSAT.\n";
				}
			}
		} while (true);

		// Add this objective as an assumption when resolving.
		std::vector<SatSolver::Literal> neg_objective_function_literals =
		    objective_function_literals;
		for (auto& lit : neg_objective_function_literals) {
			lit = negate(lit);
		}
		add_at_most_k_constraint(solver.get(), objective_function_literals, upper);
		add_at_most_k_constraint(solver.get(),
		                         neg_objective_function_literals,
		                         objective_function_literals.size() - upper);
	}

	for (size_t j = 0; j < ip->get().variable_size(); ++j) {
		auto value = sat_solution.at(literals.at(j));
		minimum_core_assert(value == 1 || value == 0);
		ip->set_solution(j, value);
	}

	// Check feasibility just to make sure everything is alright.
	minimum_core_assert(ip->is_feasible());
	return true;
}

void internal_subset(const std::vector<int>& set,
                     int left,
                     int index,
                     std::vector<int>* scratch_space,
                     std::vector<std::vector<int>>* all_subsets) {
	if (left == 0) {
		all_subsets->push_back(*scratch_space);
		return;
	}
	if (left > set.size() - index) {
		// We don’t have enough elements left to create a subset.
		return;
	}
	for (std::size_t i = index; i < set.size(); i++) {
		scratch_space->push_back(set[i]);
		internal_subset(set, left - 1, i + 1, scratch_space, all_subsets);
		scratch_space->pop_back();
	}
}

size_t choose(size_t n, size_t k) {
	if (k == 0) {
		return 1;
	}
	return (n * choose(n - 1, k - 1)) / k;
}

void generate_subsets(const std::vector<int>& set,
                      int subset_size,
                      std::vector<std::vector<int>>* output) {
	size_t num_subsets = choose(set.size(), subset_size);
	if (num_subsets > 50000000) {
		// Maybe change this limit in the future.
		throw std::runtime_error("Too many subsets. Choose a better algorithm.");
	}

	output->clear();
	output->reserve(num_subsets);
	std::vector<int> scratch_space;
	scratch_space.reserve(subset_size);
	internal_subset(set, subset_size, 0, &scratch_space, output);
}

std::unique_ptr<IP> read_CNF(std::istream& in) {
	using namespace std;
	auto ip = make_unique<IP>();

	string line;
	int nvars = 0;
	int nclauses = 0;
	while (getline(in, line)) {
		if (line.empty() || line[0] == 'c') {
			continue;
		}

		stringstream lin(line);
		if (line[0] == 'p') {
			char _;
			string cnf;
			lin >> _ >> cnf >> nvars >> nclauses;
			check(bool(lin), "Could not parse \"", line, "\".");
			check(cnf == "cnf", "Not a CNF file.");
			break;
		}
	}

	auto vars = ip->add_boolean_vector(nvars);
	for (int i = 0; i < nclauses; ++i) {
		Sum clause;
		while (true) {
			int var;
			in >> var;
			check(bool(in), "Could not read a number from the CNF file.");
			if (var > 0) {
				clause += vars.at(var - 1);
			} else if (var < 0) {
				clause += 1 - vars.at(-var - 1);
			} else {
				break;
			}
		}
		ip->add_constraint(clause >= 1);
	}

	return ip;
}
}  // namespace linear
}  // namespace minimum
