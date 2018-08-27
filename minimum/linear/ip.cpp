
#include <cmath>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <minimum/core/check.h>
#include <minimum/core/numeric.h>
#include <minimum/core/range.h>
#include <minimum/linear/constraint.h>
#include <minimum/linear/easy-ip-internal.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/pseudoboolean.h>
#include <minimum/linear/pseudoboolean_constraint.h>
#include <minimum/linear/sum.h>
#include <minimum/linear/variable.h>

using minimum::core::check;
using minimum::core::range;
using minimum::core::to_string;

namespace minimum {
namespace linear {

IP::IP() : impl(new Implementation(this)) {}

IP::IP(const void* proto_data, int length) : impl(new Implementation(this)) {
	impl->ip.ParseFromArray(proto_data, length);
	impl->solution.mutable_primal()->Resize(impl->ip.variable_size(),
	                                        std::numeric_limits<double>::quiet_NaN());
	impl->solution.mutable_dual()->Resize(impl->ip.constraint_size(),
	                                      std::numeric_limits<double>::quiet_NaN());
}

IP::~IP() {
	if (impl) {
		delete impl;
		impl = nullptr;
	}
}

Variable IP::add_variable(VariableType type, double this_cost) {
	auto var = impl->ip.add_variable();
	impl->solution.add_primal(nan);
	auto bound = var->mutable_bound();

	if (type == Boolean) {
		bound->set_lower(0.0);
		bound->set_upper(1.0);
	} else {
		bound->set_lower(-1e100);
		bound->set_upper(1e100);
	}
	impl->variable_types.push_back(type);
	var->set_cost(this_cost);

	auto index = impl->ip.variable_size() - 1;
	if (type == Boolean || type == Integer) {
		var->set_type(proto::Variable_Type_INTEGER);
	}
	return Variable(index, this);
}

BooleanVariable IP::add_boolean(double this_cost) {
	auto variable = add_variable(Boolean, this_cost);
	return BooleanVariable(variable);
}

Sum IP::add_variable_as_booleans(int lower_bound, int upper_bound) {
	Sum sum = 0;
	Sum constraint = 0;
	for (int i = lower_bound; i <= upper_bound; ++i) {
		auto var = add_boolean();
		sum += i * var;
		constraint += var;
	}
	add_constraint(std::move(constraint) == 1);
	return sum;
}

Sum IP::add_variable_as_booleans(const std::initializer_list<int>& values) {
	Sum sum = 0;
	Sum constraint = 0;
	for (int i : values) {
		auto var = add_boolean();
		sum += i * var;
		constraint += var;
	}
	add_constraint(std::move(constraint) == 1);
	return sum;
}

vector<Variable> IP::add_vector(int n, VariableType type, double this_cost) {
	vector<Variable> v;
	for (int i = 0; i < n; ++i) {
		v.push_back(add_variable(type, this_cost));
	}
	return v;
}

vector<BooleanVariable> IP::add_boolean_vector(int n, double this_cost) {
	vector<BooleanVariable> v;
	for (int i = 0; i < n; ++i) {
		v.push_back(add_boolean(this_cost));
	}
	return v;
}

template <typename Var>
vector<vector<Var>> grid_creator(int m, int n, std::function<Var()> var_creator) {
	minimum_core_assert(m >= 0 && n >= 0);
	vector<vector<Var>> grid;

	for (int i = 0; i < m; ++i) {
		grid.push_back(vector<Var>());
		for (int j = 0; j < n; ++j) {
			grid.back().push_back(var_creator());
		}
	}

	return grid;
}

vector<vector<Variable>> IP::add_grid(int m, int n, VariableType type, double this_cost) {
	auto var_creator = [&]() { return add_variable(type, this_cost); };
	return grid_creator<Variable>(m, n, var_creator);
}

vector<vector<BooleanVariable>> IP::add_boolean_grid(int m, int n, double this_cost) {
	auto var_creator = [&]() { return add_boolean(this_cost); };
	return grid_creator<BooleanVariable>(m, n, var_creator);
}

template <typename Var>
vector<vector<vector<Var>>> cube_creator(int m, int n, int o, std::function<Var()> var_creator) {
	minimum_core_assert(m >= 0 && n >= 0 && o >= 0);
	vector<vector<vector<Var>>> grid;

	for (int i = 0; i < m; ++i) {
		grid.push_back(vector<vector<Var>>());
		for (int j = 0; j < n; ++j) {
			grid.back().push_back(vector<Var>());
			for (int k = 0; k < o; ++k) {
				grid.back().back().push_back(var_creator());
			}
		}
	}

	return grid;
}

vector<vector<vector<Variable>>> IP::add_cube(
    int m, int n, int o, VariableType type, double this_cost) {
	auto var_creator = [&]() { return add_variable(type, this_cost); };
	return cube_creator<Variable>(m, n, o, var_creator);
}

vector<vector<vector<BooleanVariable>>> IP::add_boolean_cube(int m,
                                                             int n,
                                                             int o,
                                                             double this_cost) {
	auto var_creator = [&]() { return add_boolean(this_cost); };
	return cube_creator<BooleanVariable>(m, n, o, var_creator);
}

IP::VariableType IP::get_type(Variable x) const {
	impl->check_creator(x);
	return impl->variable_types.at(x.index);
}

void IP::mark_variable_as_helper(Variable x) {
	impl->ip.mutable_variable(x.index)->set_is_helper(true);
}

// Adds the constraint
//    L ≤ constraint ≤ U
DualVariable IP::add_constraint(double L, const Sum& sum, double U) {
	impl->check_creator(sum);

	// Special case for a constant sum.
	if (sum.size() == 0) {
		check(L <= sum.constant() && sum.constant() <= U,
		      "A constraint that is always false may not be added.");
		return {};
	}

	auto indices = sum.indices();
	auto values = sum.values();
	vector<std::pair<int, double>> entries_with_duplicates;
	entries_with_duplicates.reserve(sum.size());
	for (size_t i = 0; i < sum.size(); ++i) {
		entries_with_duplicates.emplace_back(indices[i], values[i]);
	}

	// Remove duplicates and add the entries.
	std::sort(entries_with_duplicates.begin(), entries_with_duplicates.end());
	vector<std::pair<int, double>> entries_to_add;
	entries_to_add.reserve(sum.size());
	for (int ind = 0; ind < entries_with_duplicates.size(); ++ind) {
		int col = entries_with_duplicates[ind].first;
		double this_value = entries_with_duplicates[ind].second;

		ind++;
		for (; ind < entries_with_duplicates.size() && col == entries_with_duplicates[ind].first;
		     ++ind) {
			this_value += entries_with_duplicates[ind].second;
		}
		ind--;

		if (this_value != 0) {
			entries_to_add.emplace_back(col, this_value);
		}
	}

	// Special case for a single-variable sum.
	if (entries_to_add.size() == 0) {
		// The entries can now be empty if, e.g., the constraint x <= x is added.
		check(L <= sum.constant() && sum.constant() <= U,
		      "A constraint that is always false may not be added.");
		return {};
	} else if (entries_to_add.size() == 1 && !impl->is_in_exists_mode) {
		auto col = entries_to_add[0].first;
		double value = entries_to_add[0].second;
		if (value == 0) {
			check(L <= sum.constant() && sum.constant() <= U,
			      "A constraint that is always false may not be added.");
		}

		if (impl->ip.variable(col).is_convex()) {
			check(value >= 0, "Can not make constraint convex.");
			check(L <= -1e100, "Can not make constraint convex.");
		}

		auto bound = impl->ip.mutable_variable(col)->mutable_bound();
		if (value > 0) {
			// add_bounds((L - sum.constant()) / value, Variable(col, this), (U - sum.constant()) /
			// value);
			bound->set_lower(std::max((L - sum.constant()) / value, bound->lower()));
			bound->set_upper(std::min((U - sum.constant()) / value, bound->upper()));
		} else {
			// add_bounds((U - sum.constant()) / value, Variable(col, this), (L - sum.constant()) /
			// value);
			bound->set_lower(std::max((U - sum.constant()) / value, bound->lower()));
			bound->set_upper(std::min((L - sum.constant()) / value, bound->upper()));
		}
		return {};
	}

	auto constraint = impl->ip.add_constraint();
	impl->solution.add_dual(nan);
	auto bound = constraint->mutable_bound();
	constraint->mutable_sum()->Reserve(entries_to_add.size());
	bound->set_lower(L - sum.constant());
	bound->set_upper(U - sum.constant());

	int row_index = static_cast<int>(impl->ip.constraint_size() - 1);
	bool all_variables_have_bounds = true;
	double constraint_min = 0;
	double constraint_max = 0;

	for (auto& entry : entries_to_add) {
		int col = entry.first;
		double value = entry.second;

		if (impl->ip.variable(col).is_convex()) {
			check(value >= 0, "Can not make constraint convex.");
			check(L <= -1e100, "Can not make constraint convex.");
		}

		auto sum_entry = constraint->add_sum();
		sum_entry->set_coefficient(value);
		sum_entry->set_variable(col);

		const double lb = impl->ip.variable(col).bound().lower();
		const double ub = impl->ip.variable(col).bound().upper();
		if (lb < -1e10 || ub > 1e10) {
			all_variables_have_bounds = false;
		}
		if (value > 0) {
			constraint_min += value * lb;
			constraint_max += value * ub;
		} else {
			constraint_min += value * ub;
			constraint_max += value * lb;
		}
	}

	if (impl->is_in_exists_mode) {
		auto slack = add_variable(Real);
		auto indicator = impl->exists_variables.back();
		auto entry = constraint->add_sum();
		entry->set_coefficient(1.0);
		entry->set_variable(slack.index);

		double M = 10000;
		if (all_variables_have_bounds) {
			M = 2 * (constraint_max - constraint_min);
		}
		minimum_core_assert(M >= 0);
		impl->is_in_exists_mode = false;
		// If indicator is zero, the slack must be zero. If it is one,
		// the slack can be anything.
		add_constraint(slack <= M * indicator);
		add_constraint(slack >= -M * indicator);
		impl->is_in_exists_mode = true;
	}

	return {std::size_t(row_index), this};
}

DualVariable IP::add_constraint(const Constraint& constraint) {
	return add_constraint(constraint.lower_bound, constraint.sum, constraint.upper_bound);
}

std::vector<DualVariable> IP::add_constraint(const ConstraintList& list) {
	std::vector<DualVariable> duals;
	for (size_t i = 0; i < list.size(); ++i) {
		duals.emplace_back(add_constraint(list.constraints()[i]));
	}
	return duals;
}

void IP::add_constraint(const BooleanVariable& variable) { add_bounds(1, variable, 1); }

void IP::add_pseudoboolean_constraint(PseudoBooleanConstraint&& constraint) {
	// std::clog << "Adding " << constraint.lower_bound << " <= " <<  constraint.sum << " <= " <<
	// constraint.upper_bound << "\n";
	impl->pseudoboolean_constraints.emplace_back(std::move(constraint));
}

// Adds the constraint
//    L ≤ variable ≤ U
void IP::add_bounds(double L, const Variable& variable, double U) {
	impl->check_creator(variable);
	if (get_type(variable) == Boolean) {
		check(L == 0 || L == 1, "Lower bound of a boolean variable needs to be 0 or 1.");
		check(U == 0 || U == 1, "Lower bound of a boolean variable needs to be 0 or 1.");
	}
	auto bound = impl->ip.mutable_variable(variable.index)->mutable_bound();
	L = std::max(L, bound->lower());
	U = std::min(U, bound->upper());
	check(L <= U, "Lower bound can not be higher than the upper bound: ", L, " > ", U, ".");
	bound->set_lower(L);
	bound->set_upper(U);
}

void IP::add_objective(const Sum& sum) {
	impl->check_creator(sum);

	auto indices = sum.size() > 0 ? sum.indices() : nullptr;
	auto values = sum.size() > 0 ? sum.values() : nullptr;
	for (size_t i = 0; i < sum.size(); ++i) {
		auto index = indices[i];
		auto value = values[i];
		if (impl->ip.variable(index).is_convex()) {
			check(value >= 0, "Can not add a convex term with a negative coefficient.");
		}
		auto var = impl->ip.mutable_variable(index);
		var->set_cost(var->cost() + value);
	}

	impl->ip.set_objective_constant(impl->ip.objective_constant() + sum.constant());
}

void IP::clear_objective() {
	for (auto j : range(impl->ip.variable_size())) {
		impl->ip.mutable_variable(j)->clear_cost();
	}
}

void IP::add_pseudoboolean_objective(PseudoBoolean pb) {
	impl->pseudoboolean_objective.emplace_back(std::move(pb));
}

void IP::Implementation::check_creator(const Variable& t) const {
	check(t.creator == nullptr || t.creator == creator, "Variable comes from a different solver.");
}

void IP::Implementation::check_creator(const Sum& t) const {
	check(t.creator() == nullptr || t.creator() == creator, "Sum comes from a different solver.");
}

double IP::get_solution(const Variable& variable) const {
	impl->check_creator(variable);
	minimum_core_assert(variable.index < impl->solution.primal_size());
	return impl->solution.primal(variable.index);
}

bool IP::get_solution(const BooleanVariable& variable) const {
	impl->check_creator(variable);
	minimum_core_assert(variable.index < impl->solution.primal_size());
	return impl->solution.primal(variable.index) > 0.5;
}

double IP::get_solution(const Sum& sum) const {
	impl->check_creator(sum);

	double value = sum.constant();
	for (size_t i = 0; i < sum.size(); ++i) {
		auto ind = sum.indices()[i];
		minimum_core_assert(0 <= ind && ind < impl->solution.primal_size());
		value += sum.values()[i] * impl->solution.primal(ind);
	}
	return value;
}

bool IP::get_solution(const LogicalExpression& expression) const {
	return get_solution(expression.sum) > 0.5;
}

double IP::get_entire_objective() const {
	minimum_core_assert(impl->ip.variable_size() == impl->solution.primal_size());
	double value = impl->ip.objective_constant();
	for (auto i : range(impl->ip.variable_size())) {
		value += impl->ip.variable(i).cost() * impl->solution.primal(i);
	}
	return value;
}

double IP::get_solution(const DualVariable& variable) const {
	minimum_core_assert(variable.get_index() < impl->solution.dual_size());
	return impl->solution.dual(variable.get_index());
}

void IP::clear_solution() {
	for (int j = 0; j < impl->solution.primal_size(); ++j) {
		impl->solution.set_primal(j, nan);
	}
	for (int i = 0; i < impl->solution.dual_size(); ++i) {
		impl->solution.set_dual(i, nan);
	}
}

size_t IP::get_number_of_variables() const { return impl->ip.variable_size(); }

void IP::add_max_consecutive_constraints(int N, const std::vector<Sum>& variables) {
	if (N >= variables.size()) {
		return;
	}

	for (int d = 0; d < variables.size() - N; ++d) {
		Sum active_in_window = 0;
		for (int d2 = d; d2 < d + N + 1; ++d2) {
			active_in_window += variables.at(d2);
		}
		add_constraint(std::move(active_in_window) <= N);
	}
}

void IP::add_min_consecutive_constraints(int N,
                                         const std::vector<Sum>& variables,
                                         bool OK_at_the_border) {
	if (N <= 1) {
		return;
	}
	minimum_core_assert(N <= variables.size());
	if (N == variables.size()) {
		for (auto& var : variables) {
			add_constraint(var == 1);
		}
		return;
	}

	for (int window_size = 1; window_size <= N - 1; ++window_size) {
		// Look for windows of size minimum - 1 along with the
		// surrounding slots.
		//
		// […] [x1] [y1] [x2] […]
		//
		// x1 = 0 ∧ x2 = 0 ⇒ y1 = 0
		// ⇔
		// x1 + x2 - y1 ≥ 0
		//
		// Then add windows with more y variables. E.g.
		//
		// x1 + x2 - y1 - y2 - y3 ≥ -2.

		for (int window_start = 0; window_start < variables.size() - window_size + 1;
		     ++window_start) {
			Sum constraint = 0;

			if (window_start - 1 >= 0) {
				constraint += variables.at(window_start - 1);
			} else if (OK_at_the_border) {
				continue;
			}

			for (int i = window_start; i < window_start + window_size; ++i) {
				constraint -= variables.at(i);
			}

			if (window_start + window_size < variables.size()) {
				constraint += variables.at(window_start + window_size);
			} else if (OK_at_the_border) {
				continue;
			}

			add_constraint(std::move(constraint) >= -window_size + 1);
		}
	}
}

std::size_t IP::matrix_size() const {
	std::size_t size = 0;
	for (auto& constraint : impl->ip.constraint()) {
		size += constraint.sum_size();
	}
	return size;
}

double IP::get_objective_constant() const { return impl->ip.objective_constant(); }

void IP::set_solution(std::size_t index, double value) {
	minimum_core_assert(index < impl->solution.primal_size());
	impl->solution.set_primal(index, value);
}

void IP::set_dual_solution(std::size_t index, double value) {
	minimum_core_assert(
	    index < impl->solution.dual_size(), index, ">=", impl->solution.dual_size());
	impl->solution.set_dual(index, value);
}

bool IP::check_invariants() const {
	minimum_core_assert(impl->ip.variable_size() == impl->solution.primal_size());
	minimum_core_assert(impl->ip.constraint_size() == impl->solution.dual_size());

	std::unordered_set<int> vars;
	for (auto& constraint : impl->ip.constraint()) {
		if (constraint.bound().lower() > constraint.bound().upper()) {
			return false;
		}

		vars.clear();
		for (auto& entry : constraint.sum()) {
			minimum_core_assert(0 <= entry.variable()
			                    && entry.variable() < impl->ip.variable_size());
			minimum_core_assert(vars.find(entry.variable()) == vars.end());
			vars.insert(entry.variable());
		}
	}
	for (auto& variable : impl->ip.variable()) {
		if (variable.bound().lower() > variable.bound().upper()) {
			return false;
		}
	}
	return true;
}

bool IP::is_feasible(double eps) const {
	check_invariants();

	for (auto j : range(impl->ip.variable_size())) {
		auto& variable = impl->ip.variable(j);
		auto solution = impl->solution.primal(j);
		if (solution != solution) {
			return false;
		}
		if (!minimum::core::is_feasible(
		        variable.bound().lower(), solution, variable.bound().upper(), eps)) {
			return false;
		}
	}

	for (auto& constraint : impl->ip.constraint()) {
		double row_sum = 0;
		for (auto& entry : constraint.sum()) {
			auto var = impl->solution.primal(entry.variable());
			row_sum += entry.coefficient() * var;
		}
		if (!minimum::core::is_feasible(
		        constraint.bound().lower(), row_sum, constraint.bound().upper(), eps)) {
			return false;
		}
	}
	return true;
}

bool IP::is_feasible_and_integral(double feasible_eps, double integral_eps) const {
	if (!is_feasible(feasible_eps)) {
		return false;
	}

	for (auto j : range(impl->ip.variable_size())) {
		if (impl->ip.variable(j).type() == proto::Variable_Type_INTEGER) {
			auto sol = impl->solution.primal(j);
			int int_sol = int(sol + 0.5);
			if (std::abs(sol - int_sol) > integral_eps) {
				return false;
			}
		}
	}
	return true;
}

bool IP::is_dual_feasible(double eps) const {
	minimum_core_assert(impl->ip.variable_size() == impl->solution.primal_size());
	minimum_core_assert(impl->ip.constraint_size() == impl->solution.dual_size());

	vector<double> ATy(impl->ip.variable_size(), 0);
	for (auto i : range(impl->ip.constraint_size())) {
		auto& constraint = impl->ip.constraint(i);

		for (auto& entry : constraint.sum()) {
			ATy.at(entry.variable()) += entry.coefficient() * impl->solution.dual(i);
		}
		if (constraint.bound().upper() >= 1e100
		    && !minimum::core::is_feasible(0., impl->solution.dual(i), 1e100, eps)) {
			// std::clog << "Constraint upper " << constraint.bound().DebugString() << " but " <<
			// constraint.dual_solution() << "\n";
			return false;
		}
		if (constraint.bound().lower() <= -1e100
		    && !minimum::core::is_feasible(-1e100, impl->solution.dual(i), 0., eps)) {
			// std::clog << "Constraint lower " << constraint.bound().DebugString() << " but " <<
			// constraint.dual_solution() << "\n";
			return false;
		}
	}

	for (auto j : range(impl->ip.variable_size())) {
		auto primal = impl->solution.primal(j);
		auto& bound = impl->ip.variable(j).bound();
		bool is_boundary = bound.lower() + eps >= primal || primal >= bound.upper() - eps;
		if (is_boundary) {
			// We don’t have duals for box constraints, so we can not check this row.
			continue;
		}
		// But for interior points, we know that dual is zero and can proceed.
		if (!minimum::core::is_equal(ATy[j], impl->ip.variable(j).cost(), eps)) {
			return false;
		}
	}

	return true;
}

void IP::start_exists() {
	if (impl->is_in_exists_mode) {
		throw std::runtime_error("Nested IP::exists are not allowed.");
	}
	impl->is_in_exists_mode = true;
	next_exists();
}

void IP::end_exists() {
	minimum_core_assert(impl->is_in_exists_mode);
	impl->is_in_exists_mode = false;

	// We have a number of indicator variables, indicating if
	// a particular set of constraints are broken. They can
	// not all be zero.
	Sum sum = 0;
	for (auto& var : impl->exists_variables) {
		sum += var;
	}
	add_constraint(std::move(sum) <= impl->exists_variables.size() - 1);
	impl->exists_variables.clear();
}

void IP::next_exists() {
	minimum_core_assert(impl->is_in_exists_mode);
	impl->exists_variables.emplace_back(add_variable(Boolean));
}

const Sum& IP::linearize_pseudoboolean_term(const std::vector<int> indices) {
	auto itr = impl->monomial_to_sum.find(indices);
	if (itr != end(impl->monomial_to_sum)) {
		return itr->second;
	}

	auto check_bounds = [this](int index) {
		auto& bound = impl->ip.variable(index).bound();
		check((bound.lower() == 0 || bound.lower() == 1)
		          && (bound.upper() == 0 || bound.upper() == 1),
		      "Variables in a PseudoBoolean expression need to be in {0, 1}.");
	};

	Sum result;
	if (indices.empty()) {
		result = 1.0;
	} else if (indices.size() == 1) {
		check_bounds(indices.back());
		result = Variable(indices.back(), this);
	} else {
		// Add a variable y equivalent to the monomial.
		auto y = add_boolean();
		Sum x_sum = 0;
		for (auto index : indices) {
			check_bounds(index);
			Variable x(index, this);
			add_constraint(y <= x);
			x_sum += x;
		}
		add_constraint(y >= x_sum - indices.size() + 1);
		result = y;
	}

	// Put the result in the map and return a reference to the sum.
	auto emplace_result = impl->monomial_to_sum.emplace(make_pair(indices, std::move(result)));
	minimum_core_assert(emplace_result.second);
	return emplace_result.first->second;
}

void IP::linearize_pseudoboolean_terms() {
	Sum objective_sum = 0;
	for (const auto& objective : impl->pseudoboolean_objective) {
		for (auto& term : objective.indices()) {
			auto& monomial = term.first;
			objective_sum += term.second * linearize_pseudoboolean_term(monomial);
		}
	}
	add_objective(objective_sum);
	impl->pseudoboolean_objective.clear();

	for (const auto& constraint : impl->pseudoboolean_constraints) {
		Sum constraint_sum = 0;
		for (auto& term : constraint.sum.indices()) {
			auto& monomial = term.first;
			constraint_sum += term.second * linearize_pseudoboolean_term(monomial);
		}
		add_constraint(constraint.lower_bound, constraint_sum, constraint.upper_bound);
	}
	impl->pseudoboolean_constraints.clear();
}

Sum IP::max_impl(std::initializer_list<Sum> terms) {
	auto y = add_variable(Real);
	for (auto&& sum : terms) {
		add_constraint(y >= sum);
	}
	impl->ip.mutable_variable(impl->ip.variable_size() - 1)->set_is_convex(true);
	return {y};
}

Sum IP::abs_impl(const Sum& sum) {
	auto y = add_variable(Real);
	add_constraint(y >= sum);
	add_constraint(y >= -sum);
	impl->ip.mutable_variable(impl->ip.variable_size() - 1)->set_is_convex(true);
	return {y};
}

const proto::IP& IP::get() const { return impl->ip; }
const proto::Solution& IP::get_solution() const { return impl->solution; }
}  // namespace linear
}  // namespace minimum
