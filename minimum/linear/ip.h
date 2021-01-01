// Petter Strandmark
// petter.strandmark@gmail.com
//
/// This header provides a simple modelling language in C++.
///
/// Example:
///
///		IP ip;
///		auto x = ip.add_variable();
///		auto y = ip.add_variable();
///
///		ip.add_objective(3.0*x + y);
///		ip.add_constraint(x + y <= 1);
///
///		solve(&ip);
///		cout << x << " " << y << endl;
///

#pragma once

#include <functional>
#include <initializer_list>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

#include <minimum/linear/constraint.h>
#include <minimum/linear/export.h>
#include <minimum/linear/proto.h>
#include <minimum/linear/pseudoboolean.h>
#include <minimum/linear/pseudoboolean_constraint.h>
#include <minimum/linear/sum.h>
#include <minimum/linear/variable.h>

namespace minimum {
namespace linear {

using std::size_t;
using std::vector;

constexpr double nan = std::numeric_limits<double>::quiet_NaN();
constexpr double infinity = std::numeric_limits<double>::infinity();

// Generates all subsets of a given size.
// Expose this externally for easier testing.
MINIMUM_LINEAR_API void generate_subsets(const std::vector<int>& set,
                                         int subset_size,
                                         std::vector<std::vector<int>>* output);

namespace internal {
template <typename Iterator>
class ExistsIterator {
   public:
	ExistsIterator(Iterator iterator_, Iterator end_, IP& ip_)
	    : iterator(iterator_), end(end_), ip(ip_) {}

	ExistsIterator& operator++();

	bool operator!=(const ExistsIterator& rhs) { return iterator != rhs.iterator; }

	auto operator*() { return *iterator; }

   private:
	Iterator iterator, end;
	IP& ip;
};

template <typename T>
class ExistsHelper {
   public:
	template <typename Collection>
	ExistsHelper(Collection&& collection_, IP& ip_);
	~ExistsHelper();

	auto begin_iteration() const {
		typedef decltype(begin(collection)) Iterator;
		return ExistsIterator<Iterator>(begin(collection), end(collection), ip);
	}

	auto end_iteration() const {
		typedef decltype(end(collection)) Iterator;
		return ExistsIterator<Iterator>(end(collection), end(collection), ip);
	}

   private:
	std::vector<T> collection;
	IP& ip;
};

template <typename T>
auto begin(const minimum::linear::internal::ExistsHelper<T>& helper) {
	return helper.begin_iteration();
}

template <typename T>
auto end(const minimum::linear::internal::ExistsHelper<T>& helper) {
	return helper.end_iteration();
}

}  // namespace internal.

/// Represents an integer (linear) program.
///
///
class MINIMUM_LINEAR_API IP {
   public:
	template <typename T>
	friend class internal::ExistsHelper;
	template <typename T>
	friend class internal::ExistsIterator;

	enum VariableType {
		Boolean,
		Binary = Boolean,  // The variable can only be 0 or 1.
		Integer,           // The variable can only be an integer.
		Real,              // The variable can be any real number.
		Convex             // Created internally.
		                   //
		                   // The variable can be any real number, but represents
		                   // the value of a convex function and can not, e.g.,
		                   // appear with a negative coefficient in the objective
		                   // function.
	};

	IP();
	IP(const void* proto_data, int length);
	IP(const IP&) = delete;
	~IP();

	/// Adds a variable to the optimization problems. Variables must not
	/// have their values queried after the creating IP class has been
	/// destroyed.
	Variable add_variable(VariableType type = Boolean, double this_cost = 0.0);

	/// Adds a boolean variable to the optimization problems.
	BooleanVariable add_boolean(double this_cost = 0.0);

	/// Adds an integer variable which is implemented as multiple boolean
	/// variables under the hood. The returned Sum can be used just as a
	/// variable.
	Sum add_variable_as_booleans(int lower_bound, int upper_bound);
	Sum add_variable_as_booleans(const std::initializer_list<int>& values);

	/// Creates a vector of variables.
	vector<Variable> add_vector(int n, VariableType type = Boolean, double this_cost = 0.0);
	/// Creates a vector of logical variables.
	vector<BooleanVariable> add_boolean_vector(int n, double this_cost = 0.0);

	/// Creates a grid of variables.
	vector<vector<Variable>> add_grid(int m,
	                                  int n,
	                                  VariableType type = Boolean,
	                                  double this_cost = 0.0);
	/// Creates a grid of logical variables.
	vector<vector<BooleanVariable>> add_boolean_grid(int m, int n, double this_cost = 0.0);

	/// Creates a 3D grid of variables.
	vector<vector<vector<Variable>>> add_cube(
	    int m, int n, int o, VariableType type = Boolean, double this_cost = 0.0);
	/// Creates a 3D grid of logical variables.
	vector<vector<vector<BooleanVariable>>> add_boolean_cube(int m,
	                                                         int n,
	                                                         int o,
	                                                         double this_cost = 0.0);

	// Marks the variable as a helper variable. This only means that it will
	// not be taken into account when computing the next solution.
	void mark_variable_as_helper(Variable x);

	/// Adds the constraint
	///    L <= constraint <= U.
	///
	/// Returns the number of constraints added. This number
	/// may be 0 if the constraint is trivial, e.g. 0 = 0.
	DualVariable add_constraint(double L, const Sum& sum, double U);
	DualVariable add_constraint(const Constraint& constraint);
	std::vector<DualVariable> add_constraint(const ConstraintList& list);
	void add_constraint(const BooleanVariable& variable);

	/// Adds the constraint
	///    L <= variable <= U
	void add_bounds(double L, const Variable& variable, double U);

	/// Adds a linear function (plus a constant) to the
	/// objective function.
	void add_objective(const Sum& sum);
	void clear_objective();

	/// Retrieves the value of a variable in the solution.
	double get_solution(const Variable& variable) const;
	/// Retrieves the value of a variable in the solution.
	bool get_solution(const BooleanVariable& variable) const;
	/// Evaluates a sum of variables in the solution.
	double get_solution(const Sum& sum) const;
	/// Evaluates an expression with the solution.
	bool get_solution(const LogicalExpression& expression) const;
	/// Evaluates an expression with the solution.
	double get_entire_objective() const;

	/// Retrieves the value of a dual variable in the solution.
	double get_solution(const DualVariable& variable) const;

	/// Clears the primal and dual solution.
	void clear_solution();

	size_t get_number_of_variables() const;

	//
	// Special modelling member functions.
	//

	// Adds the constraint that at most N of the 0/1-variables can be 1 in
	// a rows. This will add several linear constraints.
	void add_max_consecutive_constraints(int N, const std::vector<Sum>& variables);
	// Adds the constraint that at least N of the 0/1-variables must be consecutive
	// whenever they are 1.
	// If OK_at_the_border is true, the imaginary variables surrounding the vector
	// are treated as 1.
	void add_min_consecutive_constraints(int N,
	                                     const std::vector<Sum>& variables,
	                                     bool OK_at_the_border);

	/// Allows to constraints of the type ∃x: constaint(x, y).
	///
	/// Example:
	///
	///   for (auto x: ip.exists({0.0, 1.0, 2.0, 3.0})) {
	///     ip.add_constraint(y <= sin(x));
	///   }
	///
	/// Now there must exist an x such that y ≤ sin(x).
	///
	/// Note: Currently implemented in such a way that the constraints have to
	/// be in the range ±10000, unless all variables involved have explicit
	/// bounds.
	template <typename T>
	internal::ExistsHelper<T> exists(const std::vector<T>& values) {
		return {values, *this};
	}
	template <typename T>
	internal::ExistsHelper<T> exists(const std::initializer_list<T>& values) {
		return {values, *this};
	}

	// Models the maximum function and returns it as a convex sum.
	template <typename... Args>
	Sum max(Args&&... convertible_to_sums) {
		std::initializer_list<Sum> sums = {Sum(std::forward<Args>(convertible_to_sums))...};
		return max_impl(sums);
	}

	template <typename Arg>
	Sum abs(Arg&& convertible_to_sum) {
		Sum sum(std::forward<Arg>(convertible_to_sum));
		return abs_impl(sum);
	}

	//
	// More advanced functionality
	//

	// Adds a pseudo-boolean function to the objective. This function has a different
	// name to emphasize that pseudo-boolean objectives are more expensive.
	void add_pseudoboolean_objective(PseudoBoolean pb);

	// Adds a pseudo-boolean constraint. This function has a different name to
	// emphasize that pseudo-boolean constraints are more expensive.
	void add_pseudoboolean_constraint(PseudoBooleanConstraint&& constraint);

	// Replaces all pseudoboolean terms with linear constraints and new variables.
	//
	// E.g., x1*x2*x3 is replaced by the new variable y and the constraints
	//
	//     y ≤ x1,    y ≤ x2,    y ≤ x3,    y ≥ x1 + x2 + x3 - 2.
	//
	// If a solver can not handle pseudo-boolean constraints in any other way,
	// it will call this method (or perhaps fail). This method does not normally
	// need to be called.
	void linearize_pseudoboolean_terms();

	/// Converts the problem to SAT (no objective) and saves as CNF.
	void save_CNF(const std::string& file_name);

	// Obtain a read-only copy of the internal representation.
	const proto::IP& get() const;
	const proto::Solution& get_solution() const;

	// Get the variable with the specific index. Normally not needed.
	Variable get_variable(std::size_t index) const { return Variable(index, this); }

	// Number of entries in the system matrix for all constraints.
	std::size_t matrix_size() const;

	// Low-level function for solver to set the solution.
	void set_solution(std::size_t index, double value);
	void set_dual_solution(std::size_t index, double value);

	double get_objective_constant() const;

	// Internal constistency check. Does nothing.
	// Returns false if the problem is guaranteed infeasible, otherwise true.
	// Throws if there are severe issues like mismatched sizes.
	bool check_invariants() const;
	// Whether the solution is set to a feasible one.
	bool is_feasible(double eps = 1e-9) const;
	// Whether the solution is set to a feasible one and integrality constraints
	// are fulfilled.
	bool is_feasible_and_integral(double feasible_eps = 1e-9, double integral_eps = 1e-6) const;
	// Whether the dual solution is set to a feasible one.
	bool is_dual_feasible(double eps = 1e-9) const;

   protected:
	size_t get_variable_index(const Variable& x) const { return x.index; }

	// Called when entering a new exists() block.
	void start_exists();
	// Called for each entry in an exists() block.
	void next_exists();
	// Called when finishing an exists() block.
	void end_exists();

   private:
	Sum max_impl(std::initializer_list<Sum> terms);
	// Models the absolute value function and returns it as a sum.
	Sum abs_impl(const Sum& sum);

	VariableType get_type(Variable x) const;

	const Sum& linearize_pseudoboolean_term(const std::vector<int> indices);

	class Implementation;
	Implementation* impl;
};

template <typename T>
template <typename Collection>
internal::ExistsHelper<T>::ExistsHelper(Collection&& collection_, IP& ip_)
    : collection(std::forward<Collection>(collection_)), ip(ip_) {
	ip.start_exists();
}

template <typename T>
internal::ExistsHelper<T>::~ExistsHelper() {
	ip.end_exists();
}

template <typename Iterator>
internal::ExistsIterator<Iterator>& internal::ExistsIterator<Iterator>::operator++() {
	iterator++;
	if (iterator != end) {
		ip.next_exists();
	}
	return *this;
}

// Reads an MPS file and returns an IP instance.
MINIMUM_LINEAR_API std::unique_ptr<IP> read_MPS(std::istream& in);
// Reads an CNF file and returns an IP instance.
MINIMUM_LINEAR_API std::unique_ptr<IP> read_CNF(std::istream& in);
}  // namespace linear
}  // namespace minimum
