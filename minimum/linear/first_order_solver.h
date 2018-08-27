#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <minimum/linear/ip.h>
#include <minimum/linear/solver.h>

namespace minimum {
namespace linear {

struct FirstOrderOptions;
enum class LinearConstraintType : char { Equality, LessThan, GreaterThan };

using Matrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;

// Solves the linear program
//
//   minimize c·x
//   such that Ax “op” b,  where “op” can be ≥, ≤, or = for each element.
//             l ≤ x ≤ u.
//
// Implements the first-order solver described in [1], with the preconditioner of [2].
//
// [1] Chambolle, A., & Pock, T. (2011). A first-order primal-dual algorithm for convex
//     problems with applications to imaging. Journal of Mathematical Imaging and Vision,
//     40(1), 120-145.
//
// [2] Pock, T., & Chambolle, A. (2011, November). Diagonal preconditioning for first
//     order primal-dual algorithms in convex optimization. In Computer Vision (ICCV),
//     2011 IEEE International Conference on (pp. 1762-1769). IEEE.
//
bool MINIMUM_LINEAR_API
first_order_primal_dual_solve(Eigen::VectorXd* x,         /// Primal variables (in/out).
                              Eigen::VectorXd* y,         /// Dual variables (in/out).
                              const Eigen::VectorXd& c,   /// Objective function.
                              const Eigen::VectorXd& lb,  /// Lower bound on x.
                              const Eigen::VectorXd& ub,  /// Upper bound on x.
                              const Matrix& A,            /// Constraint matrix.
                              const Eigen::VectorXd& b,   /// Right-hand side of constraints.
                              const std::vector<LinearConstraintType>& constraint_types,
                              const FirstOrderOptions& options);

// Solves the linear program
//
//   minimize c·x
//   such that Ax = b,
//             l ≤ x ≤ u.
bool MINIMUM_LINEAR_API
first_order_admm_solve(Eigen::VectorXd* x,         /// Primal variables (in/out).
                       const Eigen::VectorXd& c,   /// Objective function.
                       const Eigen::VectorXd& lb,  /// Lower bound on x.
                       const Eigen::VectorXd& ub,  /// Upper bound on x.
                       const Matrix& A,            /// Equality constraint matrix.
                       const Eigen::VectorXd& b,   /// Right-hand side of constraints.
                       const FirstOrderOptions& options);

struct MINIMUM_LINEAR_API FirstOrderOptions {
	std::size_t maximum_iterations = 5000;
	// Print interval determines how often information should be logged
	// and how often convergence should be checked.
	//
	// If this value is 0, the print interval will be gradually increased
	// as time goes by.
	std::size_t print_interval = 0;
	std::function<void(const std::string&)> log_function = nullptr;

	// Stops the solver if the relative changes in the primal and dual
	// variables are both lower than this value.
	// Also the solution is considered feasible if the maximum error
	// is less than 100 times this value.
	double tolerance = 1e-9;

	// Only used by the ADMM algorithm.
	double rho = 10;

	// Will rescale the norm of the objective function vector to match the
	// right-hand side of the constraints. This is done internally and the
	// reported objective function values are not affected.
	bool rescale_c = true;
};

class MINIMUM_LINEAR_API PrimalDualSolver : public Solver {
   public:
	PrimalDualSolver() : Solver() {}

	virtual SolutionsPointer solutions(IP* ip) const override;

	FirstOrderOptions& options() { return opt; }

   private:
	FirstOrderOptions opt;
};

class MINIMUM_LINEAR_API AdmmSolver : public Solver {
   public:
	AdmmSolver() : Solver() {}

	virtual SolutionsPointer solutions(IP* ip) const override;

	FirstOrderOptions& options() { return opt; }

   private:
	FirstOrderOptions opt;
};
}  // namespace linear
}  // namespace minimum
