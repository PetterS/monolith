// Petter Strandmark
// petter.strandmark@gmail.com
//
// [1] Chambolle, A., & Pock, T. (2011). A first-order primal-dual algorithm for convex
//     problems with applications to imaging. Journal of Mathematical Imaging and Vision,
//     40(1), 120-145.
//
// [2] Pock, T., & Chambolle, A. (2011, November). Diagonal preconditioning for first
//     order primal-dual algorithms in convex optimization. In Computer Vision (ICCV),
//     2011 IEEE International Conference on (pp. 1762-1769). IEEE.
//

#include <cstdio>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <random>
#include <string>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/core/string.h>
#include <minimum/core/time.h>
#include <minimum/linear/first_order_solver.h>
using minimum::core::range;
using minimum::core::to_string;
using minimum::core::to_string_with_separator;
using minimum::core::wall_time;

namespace minimum {
namespace linear {

bool should_print(std::size_t iteration, const FirstOrderOptions& options) {
	if (iteration <= 1) {
		return true;
	}
	std::size_t print_interval = options.print_interval;
	if (print_interval == 0) {
		print_interval = pow(10.0, floor(log10(iteration))) + 0.5;
	}
	return iteration % print_interval == 0;
}

double get_feasibility_error(const Eigen::VectorXd& x,
                             Eigen::VectorXd* temp_storage,
                             const Eigen::VectorXd& lb,
                             const Eigen::VectorXd& ub,
                             const Matrix& A,
                             const Eigen::VectorXd& b,
                             const std::vector<LinearConstraintType>& constraint_types) {
	using namespace std;
	const auto m = A.rows();
	const auto n = A.cols();
	minimum_core_assert(lb.rows() == n);
	minimum_core_assert(ub.rows() == n);
	minimum_core_assert(x.rows() == n);
	minimum_core_assert(b.rows() == m);
	minimum_core_assert(constraint_types.size() == m);

	double feasibility_error = 0;
	*temp_storage = A * x - b;

	for (ptrdiff_t i = 0; i < m; ++i) {
		if (constraint_types[i] == LinearConstraintType::Equality) {
			feasibility_error = max(feasibility_error, abs((*temp_storage)(i)));
		} else if (constraint_types[i] == LinearConstraintType::LessThan) {
			feasibility_error = max(feasibility_error, (*temp_storage)(i));
		} else if (constraint_types[i] == LinearConstraintType::GreaterThan) {
			feasibility_error = max(feasibility_error, -(*temp_storage)(i));
		}
	}

	for (ptrdiff_t j = 0; j < n; ++j) {
		feasibility_error = max(feasibility_error, lb(j) - x(j));
		feasibility_error = max(feasibility_error, x(j) - ub(j));
	}

	double denominator = max(x.maxCoeff(), -x.minCoeff());
	if (denominator <= 1e-10) {
		denominator = 1;
	}

	return feasibility_error / denominator;
}

// Will use x_prev as a temporary storage after
// examining it.
bool check_convergence_and_log(std::ptrdiff_t iteration,
                               const Eigen::VectorXd& x,
                               Eigen::VectorXd* x_prev,
                               const Eigen::VectorXd& y,
                               const Eigen::VectorXd& y_prev,
                               const Eigen::VectorXd& c,
                               const Eigen::VectorXd& lb,
                               const Eigen::VectorXd& ub,
                               const Matrix& A,
                               const Eigen::VectorXd& b,
                               const std::vector<LinearConstraintType>& constraint_types,
                               const FirstOrderOptions& options,
                               double start_time) {
	using namespace Eigen;
	using namespace std;

	bool is_converged = false;

	// cerr << "x      = " << (*x).transpose() << endl;
	// cerr << "Ax     = " << (A*(*x)).transpose() << endl;
	// cerr << "Ax - b = " << (A*(*x) - b).transpose() << endl;

	auto get_relative_change = [iteration](const Eigen::VectorXd& x,
	                                       const Eigen::VectorXd& x_prev) -> double {
		double relative_change;
		if (iteration <= 0) {
			relative_change = std::numeric_limits<double>::quiet_NaN();
		} else {
			relative_change = (x - x_prev).norm() / (x.norm() + x_prev.norm());
			if (relative_change != relative_change) {
				// Both x and x_prev were null vectors.
				relative_change = 0;
			}
		}
		return relative_change;
	};

	double relative_change_x = get_relative_change(x, *x_prev);
	double relative_change_y = get_relative_change(y, y_prev);
	double feasibility_error = get_feasibility_error(x, x_prev, lb, ub, A, b, constraint_types);

	if (c.isZero(1e-30) && feasibility_error < 100 * options.tolerance) {
		// If there is no objective function, only feasibility is important.
		is_converged = true;
	} else if (relative_change_x < options.tolerance && relative_change_y < options.tolerance) {
		is_converged = true;
	}

	if (options.log_function) {
		double elapsed_time = wall_time() - start_time;

		ostringstream message;
		message << setw(9);
		if (iteration == -1) {
			message << "end";
		} else {
			message << iteration;
		}
		message << "   " << setw(15) << setprecision(6) << scientific << c.dot(x) << " " << setw(12)
		        << setprecision(3) << scientific << relative_change_x << " " << setw(12)
		        << setprecision(3) << scientific << relative_change_y << " " << setw(15)
		        << setprecision(6) << scientific << feasibility_error << " " << setw(8)
		        << setprecision(1) << scientific << elapsed_time;
		options.log_function(message.str());
	}

	return is_converged;
}

bool first_order_primal_dual_solve(Eigen::VectorXd* x_ptr,        /// Primal variables (in/out).
                                   Eigen::VectorXd* y_ptr,        /// Dual variables (in/out).
                                   const Eigen::VectorXd& c_org,  /// Objective function.
                                   const Eigen::VectorXd& lb,     /// Lower bound on x.
                                   const Eigen::VectorXd& ub,     /// Upper bound on x.
                                   const Matrix& A,               /// Equality constraint matrix.
                                   const Eigen::VectorXd& b,  /// Right-hand side of constraints.
                                   const std::vector<LinearConstraintType>& constraint_types,
                                   const FirstOrderOptions& options) {
	using namespace Eigen;
	using namespace std;

	double start_time = wall_time();

	VectorXd& x = *x_ptr;
	VectorXd& y = *y_ptr;

	const auto n = x.size();
	const auto m = y.size();
	minimum_core_assert(c_org.size() == n);
	minimum_core_assert(b.size() == m);
	minimum_core_assert(A.rows() == m);
	minimum_core_assert(A.cols() == n);
	minimum_core_assert(constraint_types.size() == m);

	VectorXd x_prev(n);
	VectorXd y_prev(m);

	const Matrix AT = A.transpose();

	Eigen::VectorXd c = c_org;
	double c_change_factor = 1.0;
	if (options.rescale_c) {
		auto bnorm = b.norm();
		auto cnorm = c.norm();
		if (bnorm > 1e-6 && cnorm > 1e-6) {
			c_change_factor = b.norm() / c.norm();
			c *= c_change_factor;
			y *= c_change_factor;
		}
	}

	double cnorm_1 = c.cwiseAbs().sum();
	double cnorm_2 = c.norm();
	double cnorm_inf = c.cwiseAbs().maxCoeff();

	double bnorm_1 = b.cwiseAbs().sum();
	double bnorm_2 = b.norm();
	double bnorm_inf = b.cwiseAbs().maxCoeff();

	if (options.log_function) {
		options.log_function("Problem size: " + to_string_with_separator(A.rows()) + " x "
		                     + to_string_with_separator(A.cols()) + " ("
		                     + to_string_with_separator(A.nonZeros()) + " non-zeros)");
		options.log_function(
		    to_string("|c|_1 = ", cnorm_1, ", |c|_2 = ", cnorm_2, ", |c|_∞ = ", cnorm_inf));
		options.log_function(
		    to_string("|b|_1 = ", bnorm_1, ", |b|_2 = ", bnorm_2, ", |b|_∞ = ", bnorm_inf));
		options.log_function(
		    "   Iter         Objective     Rel. ch. x   Rel. ch. y   Infeasibility    Time ");
		options.log_function(
		    "------------------------------------------------------------------------------");
		x_prev.setConstant(std::numeric_limits<double>::quiet_NaN());
		y_prev.setConstant(std::numeric_limits<double>::quiet_NaN());
		check_convergence_and_log(
		    0, x, &x_prev, y, y_prev, c_org, lb, ub, A, b, constraint_types, options, start_time);
	}

	// Compute preconditioners as in eq. (10) from [2], with alpha = 1.
	VectorXd Tvec(n);
	VectorXd Svec(m);
	Tvec.setZero();
	Svec.setZero();
	for (int k = 0; k < A.outerSize(); ++k) {
		for (Matrix::InnerIterator it(A, k); it; ++it) {
			auto i = it.row();
			auto j = it.col();
			auto value = abs(it.value());
			Tvec(j) += value;
			Svec(i) += value;
		}
	}

	for (ptrdiff_t j = 0; j < n; ++j) {
		Tvec(j) = 1.0 / Tvec(j);
	}

	for (ptrdiff_t i = 0; i < m; ++i) {
		Svec(i) = 1.0 / Svec(i);
	}

	size_t iteration;
	for (iteration = 1; iteration <= options.maximum_iterations; ++iteration) {
		bool should_check_convergence = should_print(iteration, options);

		x_prev = x;
		if (should_check_convergence) {
			y_prev = y;
		}

		// See eq. (18) from [2].

		x = x - Tvec.asDiagonal() * (AT * y + c);

		for (ptrdiff_t j = 0; j < n; ++j) {
			x(j) = max(lb(j), min(ub(j), x(j)));
		}

		y = y + Svec.asDiagonal() * (A * (2 * x - x_prev) - b);

		for (size_t i = 0; i < m; ++i) {
			if (constraint_types[i] == LinearConstraintType::LessThan) {
				y(i) = max(y(i), 0.0);
			} else if (constraint_types[i] == LinearConstraintType::GreaterThan) {
				y(i) = min(y(i), 0.0);
			}
		}

		if (should_check_convergence) {
			if (check_convergence_and_log(iteration,
			                              x,
			                              &x_prev,
			                              y,
			                              y_prev,
			                              c_org,
			                              lb,
			                              ub,
			                              A,
			                              b,
			                              constraint_types,
			                              options,
			                              start_time)) {
				break;
			}
		}
	}

	check_convergence_and_log(
	    -1, x, &x_prev, y, y_prev, c_org, lb, ub, A, b, constraint_types, options, start_time);

	double feasibility_error = get_feasibility_error(x, &x_prev, lb, ub, A, b, constraint_types);

	if (options.rescale_c) {
		y /= c_change_factor;
	}

	return feasibility_error < 100 * options.tolerance;
}

bool MINIMUM_LINEAR_API first_order_admm_solve(Eigen::VectorXd* x_ptr,
                                               const Eigen::VectorXd& c,
                                               const Eigen::VectorXd& lb,
                                               const Eigen::VectorXd& ub,
                                               const Matrix& A,
                                               const Eigen::VectorXd& b,
                                               const FirstOrderOptions& options) {
	using namespace Eigen;
	using namespace std;

	double start_time = wall_time();

	VectorXd& x = *x_ptr;

	const auto n = x.size();
	const auto m = A.rows();
	minimum_core_assert(c.size() == n);
	minimum_core_assert(A.cols() == n);

	VectorXd z(n);
	z.setZero();
	VectorXd u(n);
	u.setZero();

	const double rho = options.rho;

	std::vector<LinearConstraintType> constraint_types(m, LinearConstraintType::Equality);

	vector<Triplet<double>> triplets;

	for (auto j = decltype(n)(0); j < n; ++j) {
		triplets.emplace_back(j, j, rho);
	}

	for (int k = 0; k < A.outerSize(); ++k) {
		for (Matrix::InnerIterator it(A, k); it; ++it) {
			auto i = it.row();
			auto j = it.col();
			auto value = it.value();

			triplets.emplace_back(j, i + n, value);
			triplets.emplace_back(n + i, j, value);
		}
	}

	typedef SparseMatrix<double>::Index index;
	SparseMatrix<double, ColMajor> System(index(n + m), index(n + m));
	System.setFromTriplets(triplets.begin(), triplets.end());
	System.makeCompressed();

	/*cerr << "System = " << endl << System.toDense() << endl << endl;*/

	SparseLU<SparseMatrix<double, ColMajor>> solver;
	solver.analyzePattern(System);
	solver.factorize(System);
	auto computation_info = solver.info();
	if (computation_info != Success) {
		if (computation_info == NumericalIssue) {
			throw runtime_error("Eigen::NumericalIssue");
		} else if (computation_info == NoConvergence) {
			throw runtime_error("Eigen::NoConvergence ");
		} else if (computation_info == InvalidInput) {
			throw runtime_error("Eigen::InvalidInput ");
		} else {
			throw runtime_error("Unknown Eigen error.");
		}
	}

	VectorXd x_prev(n);
	VectorXd z_prev(m);

	VectorXd lhs(m + n);
	VectorXd xv(m + n);

	if (options.log_function) {
		options.log_function(
		    "   Iter         Objective     Rel. ch. x   Rel. ch. z   ||Ax - b||_inf   Time ");
		options.log_function(
		    "------------------------------------------------------------------------------");
		x_prev.setConstant(std::numeric_limits<double>::quiet_NaN());
		z_prev.setConstant(std::numeric_limits<double>::quiet_NaN());
		check_convergence_and_log(
		    0, x, &x_prev, z, z_prev, c, lb, ub, A, b, constraint_types, options, start_time);
	}

	size_t iteration;
	for (iteration = 1; iteration <= options.maximum_iterations; ++iteration) {
		bool should_check_convergence = should_print(iteration, options);

		if (should_check_convergence) {
			x_prev = x;
			z_prev = z;
		}

		lhs.block(0, 0, n, 1) = -c + rho * (z - u);
		lhs.block(n, 0, m, 1) = b;
		xv = solver.solve(lhs);
		x = xv.block(0, 0, n, 1);

		z = x + u;
		for (ptrdiff_t i = 0; i < n; ++i) {
			z(i) = max(lb(i), min(ub(i), z(i)));
		}

		u = u + x - z;

		if (should_check_convergence) {
			if (check_convergence_and_log(iteration,
			                              x,
			                              &x_prev,
			                              z,
			                              z_prev,
			                              c,
			                              lb,
			                              ub,
			                              A,
			                              b,
			                              constraint_types,
			                              options,
			                              start_time)) {
				// TODO
				break;
			}
		}
	}

	check_convergence_and_log(
	    -1, x, &x_prev, z, z_prev, c, lb, ub, A, b, constraint_types, options, start_time);

	double feasibility_error = get_feasibility_error(x, &x_prev, lb, ub, A, b, constraint_types);
	return feasibility_error < 100 * options.tolerance;
}

class PrimalDualSolutions : public Solutions {
   public:
	PrimalDualSolutions(IP* ip_, FirstOrderOptions opt_) : ip(ip_), opt(std::move(opt_)) {}

	void get_system_matrix(Matrix* A);
	void check_invariants() const;

	virtual bool get() override;

   private:
	IP* ip;
	FirstOrderOptions opt;
};

void PrimalDualSolutions::check_invariants() const {
	ip->check_invariants();
	for (auto& variable : ip->get().variable()) {
		minimum_core_assert(variable.type() == proto::Variable_Type_CONTINUOUS);
	}
}

SolutionsPointer PrimalDualSolver::solutions(IP* ip) const {
	return {std::make_unique<PrimalDualSolutions>(ip, opt)};
}

void PrimalDualSolutions::get_system_matrix(Matrix* A) {
	using namespace Eigen;
	using namespace std;

	auto n = ip->get().variable_size();
	auto m = ip->get().constraint_size();
	check_invariants();

	typedef unsigned int index;
	vector<Triplet<double, index>> sparse_indices;
	sparse_indices.reserve(ip->matrix_size());
	for (auto i : range(m)) {
		auto& constraint = ip->get().constraint(i);
		for (auto& entry : constraint.sum()) {
			sparse_indices.emplace_back(index(i), index(entry.variable()), entry.coefficient());
		}
	}

	A->resize(static_cast<int>(m), static_cast<int>(n));
	A->setFromTriplets(sparse_indices.begin(), sparse_indices.end());

	check_invariants();
}

bool PrimalDualSolutions::get() {
	using namespace Eigen;
	using namespace std;

	Matrix A;
	get_system_matrix(&A);

	auto n = ip->get().variable_size();
	auto m = ip->get().constraint_size();

	VectorXd lb(n);
	VectorXd ub(n);
	VectorXd b(m);
	VectorXd c(n);
	VectorXd x(n);
	VectorXd y(m);

	for (auto j : range(n)) {
		auto& variable = ip->get().variable(j);
		lb[j] = variable.bound().lower();
		ub[j] = variable.bound().upper();
		c[j] = variable.cost();
	}

	std::vector<LinearConstraintType> constraint_types;
	for (size_t i = 0; i < m; ++i) {
		auto lb = ip->get().constraint(i).bound().lower();
		auto ub = ip->get().constraint(i).bound().upper();
		if (lb == ub) {
			b[i] = lb;
			constraint_types.push_back(LinearConstraintType::Equality);
		} else if (ub < 1e100) {
			minimum_core_assert(lb <= -1e100);  // Cannot convert if both lb and ub exists.
			b[i] = ub;
			constraint_types.push_back(LinearConstraintType::LessThan);
		} else {
			minimum_core_assert(lb > -1e100);
			b[i] = lb;
			constraint_types.push_back(LinearConstraintType::GreaterThan);
		}
	}

	bool can_warmstart = true;
	for (auto j : range(n)) {
		auto sol = ip->get_solution().primal(j);
		if (sol != sol) {
			can_warmstart = false;
			break;
		}
		x[j] = sol;
	}
	for (auto i : range(m)) {
		auto sol = ip->get_solution().dual(i);
		if (sol != sol) {
			can_warmstart = false;
			break;
		}
		// Negative sign for the duals in order to match the linear programming
		// duals.
		y[i] = -sol;
	}
	if (!can_warmstart) {
		x.setZero();
		y.setZero();
	}

	bool feasible = first_order_primal_dual_solve(&x, &y, c, lb, ub, A, b, constraint_types, opt);

	for (size_t j = 0; j < n; ++j) {
		ip->set_solution(j, x(j));
	}
	for (size_t i = 0; i < m; ++i) {
		// Negative sign for the duals in order to match the linear programming
		// duals.
		ip->set_dual_solution(i, -y(i));
	}

	return feasible;
}

SolutionsPointer AdmmSolver::solutions(IP*) const {
	// Convert to equality-constrained problem.

	// for (size_t i = 0; i < m; ++i) {
	//	auto& lb = get_rhs_lower()[i];
	//	auto& ub = get_rhs_upper()[i];
	//	if (lb != ub) {
	//		auto slack = add_variable(IP::Real);
	//		add_bounds(0.0, slack, std::numeric_limits<double>::infinity());
	//		get_rows().push_back(int(i));
	//		get_cols().push_back(int(get_variable_index(slack)));

	//		if (ub < 1e100) {
	//			minimum_core_assert(lb <= -1e100);  // Cannot convert if both lb and ub exists.
	//			lb = ub;
	//			get_values().push_back(1.0);
	//		} else {
	//			minimum_core_assert(lb > -1e100);
	//			ub = lb;
	//			get_values().push_back(-1.0);
	//		}
	//	}
	//}

	return {nullptr};
}
}  // namespace linear
}  // namespace minimum
