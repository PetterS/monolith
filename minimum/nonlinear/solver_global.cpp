// Petter Strandmark.
//
// [1] Stig Skelboe, Computation of Rational Interval Functions, BIT 14, 1974.

#include <algorithm>
#include <cstdio>
#include <queue>
#include <set>
#include <sstream>

#include <minimum/core/check.h>
#include <minimum/core/time.h>
#include <minimum/nonlinear/solver.h>
using minimum::core::check;
using minimum::core::wall_time;

namespace minimum {
namespace nonlinear {

void midpoint(const IntervalVector& x, Eigen::VectorXd* x_mid) {
	x_mid->resize(x.size());
	for (int i = 0; i < x.size(); ++i) {
		(*x_mid)[i] = (x[i].get_upper() + x[i].get_lower()) / 2.0;
	}
}

double volume(const IntervalVector& x) {
	if (x.empty()) {
		return 0.0;
	}

	double vol = 1.0;
	for (auto itr = x.begin(); itr != x.end(); ++itr) {
		const auto& interval = *itr;
		vol *= interval.get_upper() - interval.get_lower();
	}
	return vol;
}

struct GlobalQueueEntry {
	IntervalVector box;
	Interval<double> bounds;

	bool operator<(const GlobalQueueEntry& rhs) const {
		return volume(this->box) < volume(rhs.box);

		// From [1].
		// return this->bounds.get_lower() > rhs.bounds.get_lower();
	}
};

typedef std::vector<GlobalQueueEntry> IntervalQueue;

std::ostream& operator<<(std::ostream& out, IntervalQueue queue) {
	while (!queue.empty()) {
		const auto box = queue.front().box;
		std::pop_heap(begin(queue), end(queue));
		queue.pop_back();
		out << box[0] << "; ";
	}
	out << std::endl;
	return out;
}

IntervalVector get_bounding_box(const IntervalQueue& queue_in,
                                double* function_lower_bound,
                                double* sum_of_volumes) {
	*sum_of_volumes = 0;

	if (queue_in.empty()) {
		return IntervalVector();
	}

	auto n = queue_in.front().box.size();
	std::vector<double> upper_bound(n, -1e100);
	std::vector<double> lower_bound(n, 1e100);

	*function_lower_bound = std::numeric_limits<double>::infinity();

	for (const auto& elem : queue_in) {
		const auto& box = elem.box;
		for (int i = 0; i < n; ++i) {
			lower_bound[i] = std::min(lower_bound[i], box[i].get_lower());
			upper_bound[i] = std::max(upper_bound[i], box[i].get_upper());
		}
		*sum_of_volumes += volume(box);
		*function_lower_bound = std::min(*function_lower_bound, elem.bounds.get_lower());
	}

	IntervalVector out(n);
	for (int i = 0; i < n; ++i) {
		out[i] = Interval<double>(lower_bound[i], upper_bound[i]);
	}

	return out;
}

//
// Splits an interval into 2^n subintervals and adds them all to the
// queue.
//
int split_interval(const Function& function,
                   const IntervalVector& x,
                   IntervalQueue* queue,
                   double upper_bound) {
	auto n = x.size();
	std::vector<int> split(n, 0);

	Eigen::VectorXd mid;
	midpoint(x, &mid);

	int evaluations = 0;

	while (true) {
		queue->emplace_back();
		GlobalQueueEntry& entry = queue->back();
		IntervalVector& x_split = entry.box;
		x_split.resize(n);

		for (int i = 0; i < n; ++i) {
			double a, b;
			if (split[i] == 0) {
				a = x[i].get_lower();
				b = mid[i];
			} else {
				a = mid[i];
				b = x[i].get_upper();
			}
			x_split[i] = Interval<double>(a, b);
		}

		entry.bounds = function.evaluate(entry.box);
		evaluations++;

		if (entry.bounds.get_lower() <= upper_bound) {
			std::push_heap(begin(*queue), end(*queue));
		} else {
			queue->pop_back();
		}

		// Move to the next binary vector
		//  000001
		//  000010
		//  000011
		//   ...
		//  111111
		int i = 0;
		split[i]++;
		while (split[i] > 1) {
			split[i] = 0;
			i++;
			if (i < n) {
				split[i]++;
			} else {
				break;
			}
		}
		if (i == n) {
			break;
		}
	}

	return evaluations;
}

void GlobalSolver::solve(const Function& function, SolverResults* results) const {
	solve_global(function, results);
}

IntervalVector GlobalSolver::solve_global(const Function& function, SolverResults* results) const {
	auto x_interval = function.get_all_variable_bounds();
	for (const auto& interval : x_interval) {
		check(interval.get_lower() > -Interval<double>::infinity
		          && interval.get_upper() < Interval<double>::infinity,
		      "GlobalSolver::solve: Need bounds for all variables.");
	}
	return solve_global(function, x_interval, results);
}

IntervalVector GlobalSolver::solve_global(const Function& function,
                                          const IntervalVector& x_interval,
                                          SolverResults* results) const {
	return private_solve(function, {}, 0, x_interval, results);
}

IntervalVector GlobalSolver::solve_global(const ConstrainedFunction& function,
                                          const IntervalVector& start_box,
                                          SolverResults* results) const {
	return private_solve(function.objective(),
	                     function.constraints(),
	                     function.feasibility_tolerance,
	                     start_box,
	                     results);
}

IntervalVector GlobalSolver::private_solve(const Function& function,
                                           const std::map<std::string, Constraint> constraints,
                                           const double feasibility_tolerance,
                                           const IntervalVector& x_interval,
                                           SolverResults* results) const {
	using namespace std;
	double global_start_time = wall_time();

	check(x_interval.size() == function.get_number_of_scalars(),
	      "solve_global: input vector does not match the function's number of scalars");
	auto n = x_interval.size();

	IntervalQueue queue;
	queue.reserve(2 * this->maximum_iterations);

	GlobalQueueEntry entry;
	entry.bounds = function.evaluate(x_interval);
	entry.box = x_interval;
	queue.push_back(entry);

	double upper_bound = Interval<double>::infinity;
	if (constraints.empty()) {
		upper_bound = entry.bounds.get_upper();
	}

	Eigen::VectorXd best_x(n);
	IntervalVector best_interval;

	int number_of_function_evaluations = 0;
	int iterations = 0;
	results->exit_condition = SolverResults::INTERNAL_ERROR;

	while (!queue.empty()) {
		double start_time = wall_time();
		iterations++;

		const auto box = queue.front().box;
		const auto bounds = queue.front().bounds;
		// Remove current element from queue.
		pop_heap(begin(queue), end(queue));
		queue.pop_back();

		// cerr << "-- Processing " << box << " with bounds " << bounds << ". Global upper bound is
		// " << upper_bound << endl;

		bool maybe_feasible = true;
		for (auto itr : constraints) {
			auto& constraint = itr.second;
			auto range = constraint.evaluate(box);
			if (constraint.type == Constraint::Type::EQUAL) {
				maybe_feasible = range.contains(0.0);
			} else {
				maybe_feasible = range.get_lower() <= 0;
			}
			if (!maybe_feasible) {
				break;
			}
		}

		if (maybe_feasible) {
			best_interval = box;
		}

		// Could this interval containt a feasible point and a point
		// with a better objective value?
		if (maybe_feasible && bounds.get_lower() < upper_bound) {
			// Evaluate middle point.
			Eigen::VectorXd x(box.size());
			midpoint(box, &x);
			double value = function.evaluate(x);

			bool is_feasible = true;
			for (auto itr : constraints) {
				double value = itr.second.evaluate(x);
				if (itr.second.type == Constraint::Type::EQUAL) {
					is_feasible = std::abs(value) <= feasibility_tolerance;
				} else {
					is_feasible = value <= feasibility_tolerance;
				}
				if (!is_feasible) {
					break;
				}
			}

			number_of_function_evaluations++;

			if (is_feasible && value < upper_bound) {
				upper_bound = value;
				best_x = x;
			}

			// Add new elements to queue.
			number_of_function_evaluations += split_interval(function, box, &queue, upper_bound);
		}
		results->function_evaluation_time += wall_time() - start_time;

		if (iterations % 100 == 0) {
			// Remove intervals from queue which cannot contain the optimum.
			bool removed = false;
			for (ptrdiff_t i = ptrdiff_t(queue.size()) - 1; i >= 0; --i) {
				if (queue[i].bounds.get_lower() > upper_bound) {
					queue.erase(queue.begin() + i);
					removed = true;
				}
			}
			if (removed) {
				make_heap(begin(queue), end(queue));
			}
		}

		start_time = wall_time();
		int log_interval = 1;
		if (iterations > 20) {
			log_interval = 10;
		}
		if (iterations > 200) {
			log_interval = 100;
		}
		if (iterations > 2000) {
			log_interval = 1000;
		}
		if (iterations > 20000) {
			log_interval = 10000;
		}
		if (iterations > 200000) {
			log_interval = 100000;
		}
		if (number_of_function_evaluations >= this->maximum_iterations) {
			log_interval = 1;
		}

		if (iterations % log_interval == 0) {
			double volumes_sum;
			double lower_bound =
			    upper_bound;  // Default lower bound if queue is empty (problem is solved).
			auto bounding_box = get_bounding_box(queue, &lower_bound, &volumes_sum);
			double vol_bounding = volume(bounding_box);

			double avg_magnitude = (std::abs(lower_bound) + std::abs(upper_bound)) / 2.0;
			double relative_gap = (upper_bound - lower_bound) / avg_magnitude;

			results->stopping_criteria_time += wall_time() - start_time;
			start_time = wall_time();

			if (this->log_function) {
				if (iterations == 1) {
					this->log_function(
					    "Iteration Q-size   l.bound   u.bound     rel.gap    bounding   volume");
					this->log_function(
					    "----------------------------------------------------------------------");
				}

				char tmp[1024];
				sprintf(tmp,
				        "%9d %6d %+10.3e %+10.3e %10.2e %10.3e %10.3e",
				        iterations,
				        int(queue.size()),
				        lower_bound,
				        upper_bound,
				        relative_gap,
				        vol_bounding,
				        volumes_sum);
				this->log_function(tmp);
			}

			results->log_time += wall_time() - start_time;

			if (relative_gap <= this->function_improvement_tolerance) {
				results->exit_condition = SolverResults::FUNCTION_TOLERANCE;
				break;
			}

			if (vol_bounding <= this->argument_improvement_tolerance) {
				results->exit_condition = SolverResults::ARGUMENT_TOLERANCE;
				break;
			}
		}

		if (number_of_function_evaluations >= this->maximum_iterations) {
			results->exit_condition = SolverResults::NO_CONVERGENCE;
			break;
		}
	}

	double tmp = 0;
	double lower_bound = upper_bound;
	auto bounding_box = get_bounding_box(queue, &lower_bound, &tmp);
	results->optimum_lower = lower_bound;
	results->optimum_upper = upper_bound;

	if (bounding_box.empty()) {
		// Problem was solved exactly (queue empty)
		if (log_function) {
			log_function("Problem solved exactly.");
		}
		minimum_core_assert(queue.empty());
		results->exit_condition = SolverResults::FUNCTION_TOLERANCE;
		bounding_box = best_interval;
	}

	function.copy_global_to_user(best_x);

	results->total_time = wall_time() - global_start_time;
	return bounding_box;
}
}  // namespace nonlinear
}  // namespace minimum
