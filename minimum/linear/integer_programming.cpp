// Petter Strandmark 2014
// petter.strandmark@gmail.com

#include <atomic>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
//#include <typeid>
#include <vector>

// Using Clp as the solver
// clang-format off
#include <coin/OsiClpSolverInterface.hpp>

#include <coin/CbcModel.hpp>
#include <coin/CbcEventHandler.hpp>
#include <coin/CbcCutGenerator.hpp>
#include <coin/CbcStrategy.hpp>
#include <coin/CbcHeuristicLocal.hpp>

#include <coin/CglGomory.hpp>
#include <coin/CglProbing.hpp>
#include <coin/CglKnapsackCover.hpp>
#include <coin/CglRedSplit.hpp>
#include <coin/CglClique.hpp>
#include <coin/CglOddHole.hpp>
#include <coin/CglFlowCover.hpp>
#include <coin/CglMixedIntegerRounding2.hpp>
#include <coin/CglPreProcess.hpp>
// clang-format on

#ifdef HAS_CPLEX
#include <coin/OsiCpxSolverInterface.hpp>
#endif

#ifdef HAS_MOSEK
#include <coin/OsiMskSolverInterface.hpp>
#endif

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/core/scope_guard.h>
#include <minimum/linear/glpk.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/easy-ip-internal.h>
#include <minimum/linear/sat_solver.h>
#include <minimum/linear/solver.h>

using minimum::core::check;
using minimum::core::range;
using minimum::core::to_string;

namespace minimum {
namespace linear {

bool solve(IP* ip) { return IPSolver().solutions(ip)->get(); }

bool solve_minisat(IP* ip) {
	IpToSatSolver solver(std::bind(minisat_solver, true));
	return solver.solutions(ip)->get();
}

class IPSolutions : public Solutions {
   public:
	IPSolutions(IP* ip_, const IPSolver& solver_) : ip(ip_), solver(solver_) { warm_start(ip); }

	// For empty problems, we return true the first time.
	virtual bool get() override {
		if (is_guaranteed_infeasible) {
			return false;
		}

		minimum_core_assert(!solve_started.load(), "Clp/Cbc is not reentrant.");
		solve_started.store(true);
		at_scope_exit(solve_started.store(false));

		if (has_solution) {
			if (!problem) {
				return false;
			}
			return next();
		} else {
			has_solution = true;
			if (!problem) {
				return true;
			}
			if (should_resolve) {
				return resolve();
			} else {
				return solve();
			}
		}
	}

	void warm_start(IP* ip_) override {
		ip = ip_;
		is_guaranteed_infeasible = false;
		has_integer_variables = false;
		should_resolve = has_solution;
		has_solution = false;

		if (ip->get().variable_size() == 0) {
			// Allow empty problem without failing.
			return;
		}

		for (auto& var : ip->get().variable()) {
			if (var.type() == proto::Variable_Type_INTEGER) {
				has_integer_variables = true;
				break;
			}
		}

		if (!ip->check_invariants()) {
			is_guaranteed_infeasible = true;
			return;
		}

		solver.get_problem(*ip, problem);
	}

	std::unique_ptr<OsiSolverInterface> problem;
	std::unique_ptr<CbcModel> model;
	bool has_integer_variables = false;

   private:
	// Cbc/Clp is unfortunately not reentrant. :-(
	static std::atomic<bool> solve_started;

	bool parse_solution();
	bool solve();
	bool resolve();
	bool next();

	IP* ip;
	const IPSolver& solver;
	bool is_guaranteed_infeasible = false;
	bool has_solution = false;
	bool should_resolve = false;
};
std::atomic<bool> IPSolutions::solve_started = ATOMIC_FLAG_INIT;

class MyEventHandler : public CbcEventHandler {
   public:
	virtual CbcAction event(CbcEvent whichEvent) {
		// If in sub tree carry on
		if (model_->parentModel()) {
			return noAction;
		}

		if (!(whichEvent == CbcEventHandler::solution
		      || whichEvent == CbcEventHandler::heuristicSolution)) {
			return noAction;
		}

		if (callback_function) {
			auto n = model_->getNumCols();
			auto best_solution = model_->bestSolution();
			const int* org_columns = model_->originalColumns();

			if (org_columns) {
				for (int i = 0; i < ip->get_number_of_variables(); ++i) {
					ip->set_solution(i, 0);
				}
				for (int i = 0; i < n; ++i) {
					ip->set_solution(org_columns[i], best_solution[i]);
				}
			} else {
				for (int i = 0; i < n; ++i) {
					ip->set_solution(i, best_solution[i]);
				}
			}

			minimum_core_assert(ip->is_feasible_and_integral());
			callback_function();
		}

		return noAction;  // carry on
	}

	MyEventHandler(const IPSolver::CallBack& callback_function_, IP* ip_, const CbcModel* model_)
	    : callback_function(callback_function_), ip(ip_), model(model_) {}

	MyEventHandler(const MyEventHandler& rhs)
	    : CbcEventHandler(rhs),
	      callback_function(rhs.callback_function),
	      ip(rhs.ip),
	      model(rhs.model) {}

	virtual ~MyEventHandler() {}

	MyEventHandler& operator=(const MyEventHandler& rhs) {
		if (this != &rhs) {
			callback_function = rhs.callback_function;
			ip = rhs.ip;
			model = rhs.model;
		}
		return *this;
	}

	virtual CbcEventHandler* clone() const { return new MyEventHandler(*this); }

   protected:
	IPSolver::CallBack callback_function;
	IP* ip;
	const CbcModel* model;
};

void IPSolver::get_problem(const IP& ip_to_solve,
                           std::unique_ptr<OsiSolverInterface>& problem) const {
	auto& ip_data = ip_to_solve.get();
	auto m = ip_data.constraint_size();
	auto n = ip_data.variable_size();

	bool has_integer_variables = false;
	for (auto& var : ip_data.variable()) {
		if (var.type() == proto::Variable_Type_INTEGER) {
			has_integer_variables = true;
			break;
		}
	}

	if (!problem) {
		auto clp_problem = new OsiClpSolverInterface;
		problem.reset(clp_problem);
		// Turn off information from the LP solver if we are
		// using branch and cut/bound since this means a lot
		// of LPs.
		if (has_integer_variables) {
			clp_problem->setLogLevel(0);
		}
	}

	vector<int> rows;
	vector<int> cols;
	vector<double> values;
	vector<double> rhs_lower;
	vector<double> rhs_upper;
	vector<double> var_lb;
	vector<double> var_ub;
	vector<double> cost;

	auto matrix_size = ip_to_solve.matrix_size();
	rows.reserve(matrix_size);
	cols.reserve(matrix_size);
	values.reserve(matrix_size);
	rhs_lower.reserve(m);
	rhs_upper.reserve(m);
	var_lb.reserve(n);
	var_ub.reserve(n);
	cost.reserve(n);

	bool last_index_present = false;
	for (auto i : range(m)) {
		auto& constraint = ip_data.constraint(i);
		for (auto& entry : constraint.sum()) {
			values.push_back(entry.coefficient());
			cols.push_back(entry.variable());
			rows.push_back(i);

			// Check if last_index is present.
			if (entry.variable() == n - 1) {
				last_index_present = true;
			}
		}
		rhs_lower.push_back(constraint.bound().lower());
		rhs_upper.push_back(constraint.bound().upper());
	}

	for (auto j : range(n)) {
		auto& var = ip_data.variable(j);
		var_lb.push_back(var.bound().lower());
		var_ub.push_back(var.bound().upper());
		cost.push_back(var.cost());
	}

	// Maybe add a dummy constraint to satisfy Clp.
	if (!last_index_present) {
		int row_index = static_cast<int>(rhs_upper.size());
		rhs_lower.push_back(-1e100);
		rhs_upper.push_back(1e100);
		rows.push_back(row_index);
		cols.push_back(n - 1);
		values.push_back(100);
	}

	bool is_warmstart =
	    problem->getNumCols() == cost.size() && problem->getNumRows() == rhs_lower.size();
	CoinWarmStart* warmstart = nullptr;
	if (is_warmstart) {
		warmstart = problem->getWarmStart();
	}

	CoinBigIndex coin_matrix_size(values.size());
	minimum_core_assert(coin_matrix_size == matrix_size || coin_matrix_size == matrix_size + 1);
	CoinPackedMatrix coinMatrix(false, rows.data(), cols.data(), values.data(), coin_matrix_size);

	problem->loadProblem(
	    coinMatrix, var_lb.data(), var_ub.data(), cost.data(), rhs_lower.data(), rhs_upper.data());
	problem->setDblParam(OsiObjOffset, -ip_to_solve.get_objective_constant());

	for (auto j : range(n)) {
		auto& var = ip_data.variable(j);
		if (var.type() == proto::Variable_Type_INTEGER) {
			problem->setInteger(static_cast<int>(j));
		}
	}

	if (is_warmstart) {
		// TODO: Just solve the problem from scratch if warm-starting failed.
		check(problem->setWarmStart(warmstart), "Warm-starting solver failed.");
	}
}

SolutionsPointer IPSolver::solutions(IP* ip_to_solve) const {
	ip_to_solve->linearize_pseudoboolean_terms();
	return {std::make_unique<IPSolutions>(ip_to_solve, *this)};
}

bool IPSolver::solve_relaxation(IP* ip_to_solve) const {
	auto solutions = std::make_unique<IPSolutions>(ip_to_solve, *this);
	solutions->has_integer_variables = false;
	return solutions->get();
}

bool IPSolutions::parse_solution() {
	OsiSolverInterface* solved_problem = nullptr;

	if (!has_integer_variables) {
		solved_problem = problem.get();

		if (solved_problem->isAbandoned()) {
			std::cerr << "-- Abandoned." << std::endl;
			return false;
		} else if (solved_problem->isProvenPrimalInfeasible()) {
			std::cerr << "-- Infeasible." << std::endl;
			return false;
		} else if (solved_problem->isProvenDualInfeasible()) {
			std::cerr << "-- Unbounded." << std::endl;
			return false;
		} else if (solved_problem->isPrimalObjectiveLimitReached()) {
			std::cerr << "-- Primal objective limit." << std::endl;
			return false;
		} else if (solved_problem->isDualObjectiveLimitReached()) {
			std::cerr << "-- Dual objective limit." << std::endl;
			return false;
		} else if (solved_problem->isIterationLimitReached()) {
			std::cerr << "-- Iteration limit." << std::endl;
			return false;
		} else if (has_integer_variables && !solved_problem->isProvenOptimal()) {
			std::cerr << "-- Not optimal." << std::endl;
			return false;
		}
	} else {
		if (model->isInitialSolveAbandoned()) {
			return false;
		}

		if (!model->isInitialSolveProvenOptimal()) {
			return false;
		}

		if (model->isProvenInfeasible()) {
			// throw std::runtime_error("Problem infeasible.");
			return false;
		}

		if (model->isContinuousUnbounded()) {
			return false;
		}

		if (model->isProvenDualInfeasible()) {
			// throw std::runtime_error("Problem unbounded.");
			return false;
		}

		if (!model->isProvenOptimal()) {
			// throw std::runtime_error("Time limit reached.");
			return false;
		}

		solved_problem = model->solver();
	}

	int numberColumns = solved_problem->getNumCols();
	const double* raw_solution = solved_problem->getColSolution();

	for (size_t i = 0; i < numberColumns; ++i) {
		ip->set_solution(i, raw_solution[i]);
	}

	int numberRows = solved_problem->getNumRows();
	const double* raw_dual_solution = solved_problem->getRowPrice();
	for (size_t i = 0; i < numberRows && i < ip->get().constraint_size(); ++i) {
		ip->set_dual_solution(i, raw_dual_solution[i]);
	}

	return true;
}

bool IPSolutions::solve() {
	if (!has_integer_variables) {
		if (solver.silent) {
			problem->messageHandler()->setLogLevel(0);
		}

		problem->initialSolve();
	} else {
		// Pass the solver with the problem to be solved to CbcModel
		model.reset(new CbcModel(*problem.get()));

		CbcMain0(*model);

		// Only the most important log messages.
		model->setLogLevel(1);
		if (solver.silent) {
			model->setLogLevel(0);
		}

		if (solver.time_limit_in_seconds > 0) {
			model->setDblParam(CbcModel::CbcMaximumSeconds, solver.time_limit_in_seconds);
		}

		if (solver.callback_function) {
			MyEventHandler my_event_handler(solver.callback_function, ip, model.get());
			model->passInEventHandler(&my_event_handler);
		}
		const char* argv2[] = {"minimum_linear", "-solve", "-quit"};
		CbcMain1(3, argv2, *model);
	}

	return parse_solution();
}

bool IPSolutions::resolve() {
	check(!has_integer_variables, "Can not warm-start integer programs.");
	if (solver.silent) {
		problem->messageHandler()->setLogLevel(0);
	}
	problem->resolve();
	return parse_solution();
}

bool IPSolutions::next() {
	if (!has_integer_variables) {
		return false;
	}

	OsiSolverInterface* refSolver = nullptr;
	OsiSolverInterface* solver = nullptr;
	const double* objective = nullptr;

	minimum_core_assert(model);
	refSolver = model->referenceSolver();
	solver = model->solver();

	objective = refSolver->getObjCoefficients();

	//
	// We add two new rows to the problem in order to get the
	// next solution. If the current solution is x = (1, 0, 1),
	//
	//    (a)  (1 - x1) + x2 + (1 - x3) >= 1
	//    (b)  objective(x) == *optimal*.
	//
	CoinPackedVector solution_cut, objective_cut;
	double best_objective = 0;
	double solution_rhs = 1.0;
	for (int iColumn = 0; iColumn < ip->get().variable_size(); iColumn++) {
		auto& variable = ip->get().variable(iColumn);
		double value = ip->get_solution().primal(iColumn);
		if (solver->isInteger(iColumn) && !variable.is_helper()) {
			// only works for 0-1 variables
			auto& bound = ip->get().variable(iColumn).bound();
			minimum_core_assert(bound.lower() == 0.0 || bound.lower() == 1.0);
			minimum_core_assert(bound.upper() == 0.0 || bound.upper() == 1.0);
			// double check integer
			minimum_core_assert(fabs(floor(value + 0.5) - value) < 1.0e-5);
			if (value > 0.5) {
				// at 1.0
				solution_cut.insert(iColumn, -1.0);
				solution_rhs -= 1.0;
			} else {
				// at 0.0
				solution_cut.insert(iColumn, 1.0);
			}
		}

		best_objective += value * objective[iColumn];
		objective_cut.insert(iColumn, objective[iColumn]);
		refSolver->setObjCoeff(iColumn, 0.0);
	}

	// now add cut
	refSolver->addRow(solution_cut, solution_rhs, COIN_DBL_MAX);
	refSolver->addRow(objective_cut, best_objective, best_objective);

	if (!has_integer_variables) {
		refSolver->branchAndBound();
	} else {
		model->resetToReferenceSolver();
		// model->setHotstartSolution(bestSolution, nullptr);

		// Do complete search
		model->branchAndBound();
	}

	return parse_solution();
}

void IPSolver::set_callback(const CallBack& callback_function_) {
	callback_function = callback_function_;
}

void IPSolver::save_MPS(const IP& ip_to_save, const std::string& file_name) const {
	std::unique_ptr<OsiSolverInterface> new_problem;
	get_problem(ip_to_save, new_problem);
	new_problem->writeMps(file_name.c_str());
}
}  // namespace linear
}  // namespace minimum
