#include <algorithm>
#include <iomanip>
#include <random>
#include <utility>
#include <vector>
using std::find_if;
using std::pair;
using std::size_t;
using std::vector;

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/core/scope_guard.h>
#include <minimum/core/string.h>
#include <minimum/core/time.h>
#include <minimum/linear/constraint_solver.h>
using namespace minimum::core;

namespace minimum {
namespace linear {
namespace constraint_solver {

using Index = int;
using Value = int;

// A variable in the optimization problem. Keeps track of bounds
// and assigned value.
struct Variable {
	Variable() : lower(0), upper(1), org_lower(0), org_upper(1) {}
	Variable(Value lower_, Value upper_)
	    : lower(lower_), upper(upper_), org_lower(lower_), org_upper(upper_) {}
	Value lower;
	Value upper;
	Value value = 0;
	Value org_lower;
	Value org_upper;
};

// A constraint in the optimization problem.
//
// The data structures in this file are a bit peculiar (for speed). A Constraint
// does not know its variables and coefficients – only the lhs, rhs and the minimum
// and maximum possible values.
class Constraint {
   public:
	Constraint(const vector<pair<const Variable*, Value>>& coefficients, Value lhs, Value rhs);
	void assume(const Variable* var, Value value, Value coef);
	void unassume(const Variable* var, Value value, Value coef);
	bool possible() const;
	void print() const;

	Value get_min() const { return min; }
	Value get_max() const { return max; }
	Value get_rhs() const { return rhs; }
	void set_rhs(Value rhs_) { rhs = rhs_; }

   private:
	Value min;
	Value max;
	Value lhs;
	Value rhs;
};

class Searcher {
   public:
	Searcher(int num_variables);
	void add_constraint(vector<pair<Index, Value>> coefficients, Value lhs, Value rhs);
	void set_objective(vector<pair<Index, Value>> coefficients);
	void set_bounds(Index var, Value lower, Value upper);
	bool search();
	const vector<Variable>& get_variables() const { return variables; }

	long long num_assumptions = 0;
	long long num_backtracks = 0;

   private:
	bool has_preprocessed = false;
	void preprocess();

	bool search_variable(Index v);
	bool assume(Index var, Value value);
	void unassume();

	void print();

	// Constant after search is started
	// ================================

	vector<Variable> variables;
	vector<Constraint> constraints;
	vector<Index> variable_order;
	struct ConstraintRef {
		ConstraintRef(Index c, Value c2) : constraint(c), coefficient_of_variable(c2) {}
		Index constraint;
		Value coefficient_of_variable;
	};
	// Holds all constraints indexed by variables.
	vector<vector<ConstraintRef>> variable_to_constraints;
	Index objective = -1;
	Value objective_value = 0;

	// Search state
	// ============

	vector<pair<Index, Value>> assumptions;
	bool first_time_with_objective = false;
	vector<int> variable_failures;
};

Constraint::Constraint(const vector<pair<const Variable*, Value>>& coefficients,
                       Value lhs_,
                       Value rhs_)
    : lhs(lhs_), rhs(rhs_) {
	minimum_core_assert(lhs <= rhs);
	min = 0;
	max = 0;
	for (auto& coef : coefficients) {
		if (coef.second >= 0) {
			min += coef.second * coef.first->lower;
			max += coef.second * coef.first->upper;
		} else {
			max += coef.second * coef.first->lower;
			min += coef.second * coef.first->upper;
		}
	}
}

void Constraint::assume(const Variable* var, Value value, Value coef) {
	if (coef >= 0) {
		min -= coef * (var->org_lower - value);
		max -= coef * (var->org_upper - value);
	} else {
		max -= coef * (var->org_lower - value);
		min -= coef * (var->org_upper - value);
	}
}

void Constraint::unassume(const Variable* var, Value value, Value coef) {
	if (coef >= 0) {
		min += coef * (var->org_lower - value);
		max += coef * (var->org_upper - value);
	} else {
		max += coef * (var->org_lower - value);
		min += coef * (var->org_upper - value);
	}
}

bool Constraint::possible() const { return lhs <= max && min <= rhs; }

void Constraint::print() const {
	std::cout << lhs << " <= ";
	// for (auto& coef : coefficients) {
	//	std::cout << coef.second << "*x" << coef.first << " ";
	//}
	std::cout << " <= " << rhs << ". ";
	std::cout << "min=" << min << ", ";
	std::cout << "max=" << max << "\n";
}

Searcher::Searcher(int num_variables)
    : variables(num_variables),
      variable_to_constraints(num_variables),
      variable_failures(num_variables, 0) {}

void Searcher::add_constraint(vector<pair<Index, Value>> coefficients, Value lhs, Value rhs) {
	has_preprocessed = false;

	vector<pair<const Variable*, Value>> internal_coefficients;
	for (auto& coef : coefficients) {
		internal_coefficients.emplace_back(&variables[coef.first], coef.second);
		variable_to_constraints[coef.first].emplace_back(constraints.size(), coef.second);
	}
	constraints.emplace_back(std::move(internal_coefficients), lhs, rhs);
}

void Searcher::set_objective(vector<pair<Index, Value>> coefficients) {
	check(objective < 0, "Can only set objective once.");
	add_constraint(std::move(coefficients),
	               std::numeric_limits<Value>::min(),
	               std::numeric_limits<Value>::max());
	objective = constraints.size() - 1;
}

void Searcher::set_bounds(Index var, Value lower, Value upper) {
	has_preprocessed = false;

	variables[var].lower = lower;
	variables[var].org_lower = lower;
	variables[var].upper = upper;
	variables[var].org_upper = upper;
}

void Searcher::preprocess() {
	has_preprocessed = true;

	variable_order.clear();
	for (auto i : range(variables.size())) {
		variable_order.push_back(i);
	}

	// Sort most used to least.
	std::stable_sort(variable_order.begin(), variable_order.end(), [this](Index a, Index b) {
		return variable_to_constraints[a].size() > variable_to_constraints[b].size();
	});

	// Remove fixed variables.
	auto itr = variable_order.begin();
	while (itr != variable_order.end()) {
		auto& var = variables[*itr];
		if (var.org_lower == var.org_upper) {
			itr = variable_order.erase(itr);
			var.value = var.org_lower;
		} else {
			++itr;
		}
	}
}

bool Searcher::search() {
	num_assumptions = 0;
	if (!has_preprocessed) {
		preprocess();
	}
	first_time_with_objective =
	    objective >= 0 && constraints[objective].get_rhs() == std::numeric_limits<Value>::max();

	bool result = search_variable(0);
	minimum_core_assert(assumptions.empty());

	if (first_time_with_objective) {
		// A search optimizing the objective always returns false since it
		// always tries to improve the objective until it fails.
		// If it found any feasible solution, the search is considered a success.
		result = constraints[objective].get_rhs() < std::numeric_limits<Value>::max();
		constraints[objective].set_rhs(objective_value);
	}

	for (auto v : variable_order) {
		// Change the lower bound of the variable so that we do not visit
		// this solution again in subsequent searches.
		variables[v].lower = variables[v].value;
	}
	// Set the bound of the last variable to be one past the current
	// value, which we have already searched.
	variables[variable_order.back()].lower = variables[variable_order.back()].value + 1;

	return result;
}

bool Searcher::search_variable(Index k) {
	auto v = variable_order[k];

	bool all_failed = true;
	for (auto value : range(variables[v].lower, variables[v].upper + 1)) {
		bool possible = assume(v, value);
		if (possible) {
			at_scope_exit(unassume(););
			all_failed = false;

			if (k + 1 < variable_order.size()) {
				bool result = search_variable(k + 1);
				if (result) {
					return true;
				}
			} else {
				// Solution found.
				for (auto& assumption : assumptions) {
					variables[assumption.first].value = assumption.second;
				}

				if (first_time_with_objective) {
					// When all variables are assigned, min and max for a
					// constraint are the same.
					objective_value = constraints[objective].get_max();
					// Try to look for a better solution from now on.
					constraints[objective].set_rhs(objective_value - 1);
				} else {
					return true;
				}
			}
		}
	}

	if (all_failed) {
		variable_failures[v]++;
	}
	num_backtracks++;

	// Restore the bounds to the original values. They may have been changed
	// in order to fast-forward the searcher to a certain state. Once we have
	// used the modified bounds once, we should restore the original ones.
	variables[v].lower = variables[v].org_lower;
	variables[v].upper = variables[v].org_upper;
	return false;
}

bool Searcher::assume(Index var, Value value) {
	num_assumptions++;

	auto variable = &variables[var];
	for (auto& c : variable_to_constraints[var]) {
		// Nice assertion to have when developing, but this function is a
		// huge hotspot.
		// minimum_core_assert(constraints[c.first].possible());

		constraints[c.constraint].assume(variable, value, c.coefficient_of_variable);
		if (!constraints[c.constraint].possible()) {
			// Roll back the assumptions already made.
			for (auto c2 = variable_to_constraints[var].begin();; ++c2) {
				constraints[c2->constraint].unassume(variable, value, c2->coefficient_of_variable);
				if (c2->constraint == c.constraint) {
					break;
				}
			}
			return false;
		}
	}
	assumptions.emplace_back(var, value);
	return true;
}

void Searcher::unassume() {
	minimum_core_assert(!assumptions.empty());
	auto& assumption = assumptions.back();
	auto variable = &variables[assumption.first];
	for (auto& c : variable_to_constraints[assumption.first]) {
		constraints[c.constraint].unassume(variable, assumption.second, c.coefficient_of_variable);
	}
	assumptions.pop_back();
}

void Searcher::print() {
	for (auto& constraint : constraints) {
		constraint.print();
	}
}
}  // namespace constraint_solver

class ConstraintSolutions : public Solutions {
   public:
	ConstraintSolutions(IP* ip_, bool silent_)
	    : ip(ip_), searcher(ip_->get_number_of_variables()), silent(silent_) {
		using namespace constraint_solver;

		auto num_rows = ip->get().constraint_size();

		for (int j = 0; j < ip->get_number_of_variables(); ++j) {
			auto& bound = ip->get().variable(j).bound();
			Value lower = bound.lower();
			check(lower == bound.lower(),
			      "ConstraintSolver::solve: Need bounded integer variables.");
			Value upper = bound.upper();
			check(upper == bound.upper(),
			      "ConstraintSolver::solve: Need bounded integer variables.");
			searcher.set_bounds(j, lower, upper);
		}

		vector<vector<pair<Index, Value>>> constraints(num_rows);
		for (auto i : range(num_rows)) {
			auto& constraint = ip->get().constraint(i);
			for (auto& entry : constraint.sum()) {
				Value v = entry.coefficient();
				minimum_core_assert(v == entry.coefficient());
				constraints[i].emplace_back(entry.variable(), v);
			}
		}
		for (auto i : range(num_rows)) {
			auto& constraint = ip->get().constraint(i);
			Value lhs = std::numeric_limits<Value>::min();
			Value rhs = std::numeric_limits<Value>::max();
			if (lhs < constraint.bound().lower()) {
				lhs = constraint.bound().lower();
				minimum_core_assert(lhs == constraint.bound().lower());
			}
			if (rhs > constraint.bound().upper()) {
				rhs = constraint.bound().upper();
				minimum_core_assert(rhs == constraint.bound().upper());
			}
			searcher.add_constraint(std::move(constraints[i]), lhs, rhs);
		}

		vector<pair<Index, Value>> objective_coefs;
		for (int j = 0; j < ip->get_number_of_variables(); ++j) {
			if (ip->get().variable(j).cost() != 0) {
				Value cost = ip->get().variable(j).cost();
				check(cost == ip->get().variable(j).cost(),
				      "ConstraintSolver::solve: Need integer costs.");
				objective_coefs.emplace_back(j, cost);
			}
		}
		if (!objective_coefs.empty()) {
			searcher.set_objective(std::move(objective_coefs));
		}
	}

	// Inherited via Solutions
	virtual bool get() override {
		auto start_time = wall_time();
		bool result = searcher.search();
		auto elapsed = wall_time() - start_time;

		if (!silent) {
			std::cerr << ip->get_number_of_variables() << " vars, " << ip->get().constraint_size()
			          << " rows. Time: " << elapsed << "\n";
			std::cerr << "-- " << searcher.num_backtracks / 1000 << "k backtracks, "
			          << searcher.num_assumptions / 1000 << "k assumptions ("
			          << to_string(int(1e9 * elapsed / searcher.num_assumptions))
			          << " ns/assump).\n";
		}

		if (result) {
			auto& variables = searcher.get_variables();
			for (auto j : range(ip->get_number_of_variables())) {
				ip->set_solution(j, variables[j].value);
			}
			minimum_core_assert(ip->is_feasible(), "Solution from searcher not feasible.");
		}
		return result;
	}

   private:
	IP* ip;
	constraint_solver::Searcher searcher;
	bool silent;
};

ConstraintSolver::~ConstraintSolver() {}

SolutionsPointer ConstraintSolver::solutions(IP* ip) const {
	ip->linearize_pseudoboolean_terms();
	return {std::make_unique<ConstraintSolutions>(ip, silent)};
}
}  // namespace linear
}  // namespace minimum
