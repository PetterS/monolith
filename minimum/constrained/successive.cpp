#include <iomanip>
#include <iostream>
#include <sstream>
using namespace std;

#include <absl/strings/substitute.h>

#include <minimum/constrained/successive.h>
#include <minimum/core/numeric.h>
#include <minimum/core/range.h>
#include <minimum/core/time.h>
#include <minimum/linear/glpk.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/solver.h>
using namespace minimum::core;

using minimum::linear::IP;
using minimum::linear::IPSolver;
using minimum::linear::Sum;
using minimum::nonlinear::ConstrainedFunction;

namespace minimum::constrained {
bool successive_linear_programming_step(minimum::linear::Solver* solver,
                                        std::optional<minimum::linear::SolutionsPointer>* solutions,
                                        ConstrainedFunction* constrained_function,
                                        double step_size) {
	const int n = constrained_function->all().get_number_of_scalars();
	Eigen::VectorXd x(n);
	constrained_function->all().copy_user_to_global(&x);

	IP ip;
	auto d = ip.add_vector(n, IP::VariableType::Real);
	for (auto i : range(n)) {
		ip.add_bounds(-step_size, d[i], step_size);
	}

	// Do some trickery in order to obtain the the global indices for the variables
	// in the objective function.
	//   (i) Prepare a vector with all indices in order.
	//   (ii) Load it into user space.
	//   (iii) Load it back into a vector via the objective function.
	// That vector now holds the global indices of the objective function.
	auto n_objective = constrained_function->objective().get_number_of_scalars();
	Eigen::VectorXd indices(n);
	for (auto i : range(n)) {
		indices[i] = i;
	}
	constrained_function->all().copy_global_to_user(indices);
	Eigen::VectorXd objective_indices(n_objective);
	constrained_function->objective().copy_user_to_global(&objective_indices);

	// Evaluate the objective function and its gradient.
	Eigen::VectorXd x_objective(n_objective);
	Eigen::VectorXd g_objective(n_objective);
	constrained_function->all().copy_global_to_user(x);
	constrained_function->objective().copy_user_to_global(&x_objective);
	auto obj = constrained_function->objective().evaluate(x_objective, &g_objective);
	ip.add_objective(obj);
	for (auto i : range(n_objective)) {
		int global_index = static_cast<int>(objective_indices[i]);
		ip.add_objective(g_objective[i] * d.at(global_index));
	}

	// Linearize every constraint and add it to the LP.
	for (auto& entry : constrained_function->constraints()) {
		auto& constraint = entry.second;
		auto m = constraint.get_number_of_scalars();
		Eigen::VectorXd x_local(m);
		Eigen::VectorXd g_local(m);
		for (auto i : range(m)) {
			x_local[i] = x[constraint.global_indices.at(i)];
		}
		double value = constraint.evaluate(x_local, &g_local);

		Sum linearization = value;
		for (auto i : range(m)) {
			linearization += g_local[i] * d[constraint.global_indices.at(i)];
		}

		if (constraint.type == minimum::nonlinear::Constraint::Type::LESS_THAN_OR_EQUAL) {
			ip.add_constraint(linearization <= 0);
		} else if (constraint.type == minimum::nonlinear::Constraint::Type::EQUAL) {
			ip.add_constraint(linearization == 0);
		}
	}

#if 0
	// If we don't have a single constraint, add one because GLPK does not
	// like a matrix with zero rows.
	if (ip.get().constraint_size() == 0) {
		ip.add_constraint(sum(d) <= d.size() * step_size);
	}

	if (glp->num_cols() == ip.get_number_of_variables()
	    && glp->num_rows() == ip.get().constraint_size()) {
		// We can warm-start this problem.
		glp->set_presolve(false);
	} else {
		// This is either the first time we solve or we need to reset the problem
		// because it is incompatible.
		// IP removes coefficients that are exactly zero so for some problems this
		// can cause a constraint to get removed. There could be edge cases where
		// one constraint is removed at the same time another one reappears. This
		// is currently not handled.
		if (glp->num_cols() != 0 || glp->num_rows() != 0) {
			std::cerr << absl::Substitute("-- Existing problem is ($0, $1); need ($2, $3). Recreating.\n",
										  glp->num_cols(),
										  glp->num_rows(),
										  ip.get_number_of_variables(),
										  ip.get().constraint_size());
		}
		*glp = minimum::linear::GlpProblem();
		glp->add_columns(ip.get_number_of_variables());
		glp->add_rows(ip.get().constraint_size());
		glp->set_presolve(true);
	}
	glp->set_silent();
	glp->set_cost(ip);
	glp->set_constant(ip.get_objective_constant());
	glp->set_col_bounds(ip);
	glp->set_row_bounds(ip);
	glp->load_matrix(ip);
	if (!glp->solve()) {
		return false;
	}
	glp->extract_solution(&ip);
#else

	if (!solutions->has_value()) {
		*solutions = solver->solutions(&ip);
	} else {
		solutions->value()->warm_start(&ip);
	}
	if (!solutions->value()->get()) {
		return false;
	}
#endif

	for (auto i : range(n)) {
		x[i] += d[i].value();
	}
	constrained_function->all().copy_global_to_user(x);
	return true;
}

namespace {
struct Info {
	double objective = 0;
	double max_violation = 0;
	double solve_time = 0;

	double sum() const { return objective + max_violation; }
};

Info compute_info(ConstrainedFunction& constrained_function) {
	Info i;
	i.objective = constrained_function.objective().evaluate();
	i.max_violation = 0;
	for (auto& itr : constrained_function.constraints()) {
		itr.second.evaluate_and_cache();
		i.max_violation = max(i.max_violation, itr.second.violation());
	}
	return i;
}

void log_info(std::function<void(std::string_view)> log_function,
              string_view iteration,
              const Info& info,
              double step,
              string_view status) {
	stringstream sout;
	sout << setw(6) << right << iteration << " " << setw(12) << right << info.objective << " "
	     << setw(12) << right << info.max_violation << " " << setw(12) << right << step << " "
	     << setw(10) << status << " " << right << setw(12) << info.solve_time;
	log_function(sout.str());
}
}  // namespace

SuccessiveLinearProgrammingSolver::Result SuccessiveLinearProgrammingSolver::solve(
    minimum::nonlinear::ConstrainedFunction* function) {
	std::optional<minimum::linear::SolutionsPointer> solutions;
	IPSolver solver;
	solver.silent = true;

	double step = initial_step_size;
	Info info = compute_info(*function);
	if (log_function) {
		log_function("| Iter  | Objective | Infeas.  |  Step size  | Action    |   LP time  |");
		log_function("+-------+-----------+----------+-------------+-----------+------------+");
		log_info(log_function, "start", info, step, "START");
	}
	for (int iteration = 1;; ++iteration) {
		double start_time = wall_time();
		while (!minimum::constrained::successive_linear_programming_step(
		    &solver, &solutions, function, step)) {
			step *= 2;
			if (step >= 1e9) {
				return Result::INFEASIBLE;
			}
			log_function("-- Infeasible problem; increasing step size.");
		}
		double time_taken = wall_time() - start_time;

		Info new_info = compute_info(*function);
		new_info.solve_time = time_taken;
		string_view status = "N/A";

		if (relative_error(new_info.objective, info.objective) < function_improvement_tolerance
		    && function->is_feasible()) {
			status = "CONVERGED";
		} else if (new_info.sum() >= info.sum()) {
			step /= 2;
			status = "DECREASE";
		} else {
			step *= 1.2;
			status = "INCREASE";
			if (step >= 1e9) {
				return Result::INFEASIBLE;
			}
		}
		if (log_function) {
			log_info(log_function, to_string(iteration), new_info, step, status);
		}
		if (status == "CONVERGED") {
			break;
		}
		info = new_info;

		if (iteration >= max_iterations) {
			return Result::NO_CONVERGENCE;
		}
	}
	return Result::OK;
}
}  // namespace minimum::constrained
