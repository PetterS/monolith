#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <streambuf>
#include <string>
#include <tuple>
#include <vector>
using namespace std;

#include <gflags/gflags.h>
#include <google/protobuf/text_format.h>

#include <minimum/core/check.h>
#include <minimum/core/color.h>
#include <minimum/core/flamegraph.h>
#include <minimum/core/grid.h>
#include <minimum/core/main.h>
#include <minimum/core/numeric.h>
#include <minimum/core/openmp.h>
#include <minimum/core/range.h>
#include <minimum/core/time.h>
#include <minimum/linear/colgen/problem.h>
#include <minimum/linear/colgen/set_partitioning_problem.h>
#include <minimum/linear/data/util.h>
#include <minimum/linear/glpk.h>
#include <minimum/linear/proto.h>
using namespace minimum::core;
using namespace minimum::linear;
using namespace minimum::linear::colgen;

DEFINE_int32(num_solutions, 1, "Number of solutions to compute. Default: 1.");
DEFINE_double(verify_solution,
              -1,
              "If set, checks that the solution is at least as good as this value.");

class ShiftGenerationColgenProblem : public Problem {
   public:
	ShiftGenerationColgenProblem(const minimum::linear::proto::ShiftGenerationProblem& problem_)
	    : Problem(problem_.period_size()), problem(problem_) {
		FLAMEGRAPH_LOG_FUNCTION;
		for (int row : range(problem.period_size())) {
			Column under_slack(problem.under_cost(), 0, 1e100);
			under_slack.add_coefficient(row, 1.0);
			under_slack.set_integer(false);
			pool.add(move(under_slack));

			Column over_slack(problem.over_cost(), 0, 1e100);
			over_slack.add_coefficient(row, -1.0);
			over_slack.set_integer(false);
			pool.add(move(over_slack));

			set_row_lower_bound(row, problem.period(row).demand());
			set_row_upper_bound(row, problem.period(row).demand());
		}

		// Generate a few random shifts to start with.
		add_shift(0, problem.shift_length_bound().min());
		add_shift(0, problem.shift_length_bound().max());
	}

	virtual void generate(const std::vector<double>& dual_variables) override {
		FLAMEGRAPH_LOG_FUNCTION;

		int best_start = -1;
		int best_end = -1;
		double best_cost = -1e100;
		for (int start : range(problem.period_size())) {
			for (int length : range(problem.shift_length_bound().min(),
			                        problem.shift_length_bound().max() + 1)) {
				if (start + length > problem.period_size()) {
					continue;
				}
				double cost = 0;
				for (int i : range(start, start + length)) {
					cost += dual_variables[i];
				}
				if (cost > best_cost) {
					best_cost = cost;
					best_start = start;
					best_end = start + length;
				}
			}
		}
		add_shift(best_start, best_end);
	}

	virtual int fix(const FixInformation& information) override {
		FLAMEGRAPH_LOG_FUNCTION;

		if (information.iteration <= 4
		    || information.objective_change < -FLAGS_objective_change_before_fixing
		    || information.objective_change > 0) {
			return 0;
		}

		int best_index_to_fix = -1;
		double best_score = 1;
		for (int i : shift_column_indices) {
			const auto& column = pool.at(i);

			double residual = integer_residual(column.solution_value);
			// Try to fix columns that are close to integral. But not exactly, since
			// those are not too interesting.
			double score = abs(residual - 0.1);
			if (score < best_score) {
				best_score = score;
				best_index_to_fix = i;
			}
		}
		int fix_to = lround(pool.at(best_index_to_fix).solution_value);
		if (abs(pool.at(best_index_to_fix).solution_value - fix_to) > 1e-6) {
			cerr << "-- Fixing column with value " << pool.at(best_index_to_fix).solution_value
			     << " to " << fix_to << ".\n";
		}
		pool.at(best_index_to_fix).fix(fix_to);
		return 1;
	}

	virtual double integral_solution_value() override {
		FLAMEGRAPH_LOG_FUNCTION;

		auto solution = get_solution();
		double cost = 0;
		vector<int> fulfilled(problem.period_size(), 0);
		for (auto& shift : solution) {
			for (int i : range(shift.first, shift.second)) {
				fulfilled[i] += 1;
			}
			cost += problem.shift_cost();
		}
		for (auto i : range(fulfilled.size())) {
			auto diff = fulfilled[i] - problem.period(i).demand();
			if (diff > 0) {
				cost += diff * problem.over_cost();
			} else if (diff < 0) {
				cost -= diff * problem.under_cost();
			}
		}
		return cost;
	}

	void add_shift(int start, int end) {
		FLAMEGRAPH_LOG_FUNCTION;
		minimum_core_assert(0 <= start && start < problem.period_size());
		minimum_core_assert(start < end && end <= problem.period_size());
		minimum_core_assert(problem.shift_length_bound().min() <= end - start);
		minimum_core_assert(end - start <= problem.shift_length_bound().max());

		Column column(problem.shift_cost(), 0, 1000'000);
		for (int row : range(start, end)) {
			column.add_coefficient(row, 1.0);
		}

		auto s1 = pool.size();
		pool.add(std::move(column));
		auto s2 = pool.size();
		if (s1 < s2) {
			shift_column_indices.push_back(s1);
		}
	}

	vector<pair<int, int>> get_solution() const {
		FLAMEGRAPH_LOG_FUNCTION;
		vector<pair<int, int>> result;
		for (int i : shift_column_indices) {
			const auto& column = pool.at(i);
			int min = problem.period_size();
			int max = 0;
			for (auto& entry : column) {
				min = std::min(min, entry.row);
				max = std::max(max, entry.row);
			}

			for (int i = 0; i < lround(column.solution_value); ++i) {
				result.emplace_back(min, max + 1);
			}
		}
		sort(result.begin(), result.end(), [](const auto& p1, const auto& p2) {
			if (p1.first != p2.first) {
				return p1.first < p2.first;
			}
			return p1.second < p2.second;
		});
		return result;
	}

	double fraction_used() const {
		double used = 0;
		for (int i : shift_column_indices) {
			if (lround(pool.at(i).solution_value) > 0) {
				used += 1;
			}
		}
		return used / shift_column_indices.size();
	}

   private:
	const minimum::linear::proto::ShiftGenerationProblem& problem;
	vector<int> shift_column_indices;
};

void cout_row(const minimum::linear::proto::ShiftGenerationProblem& problem,
              function<void(int)> cell_content) {
	for (int i = 0; i < problem.period_size(); ++i) {
		cout << "|" << setw(3);
		cell_content(i);
	}
	cout << "|";
}

void cout_row(const minimum::linear::proto::ShiftGenerationProblem& problem) {
	for (int i = 0; i < problem.period_size(); ++i) {
		cout << "+---";
	}
	cout << "+\n";
}

int main_program(int num_args, char* args[]) {
	using namespace std;

	auto start_time = wall_time();
	minimum::linear::proto::ShiftGenerationProblem problem;

	if (num_args <= 1) {
		problem.add_period()->set_demand(1);
		problem.add_period()->set_demand(2);
		problem.add_period()->set_demand(4);
		problem.add_period()->set_demand(7);
		problem.add_period()->set_demand(8);
		problem.add_period()->set_demand(6);
		problem.add_period()->set_demand(6);
		problem.add_period()->set_demand(1);

		problem.set_over_cost(1.0);
		problem.set_under_cost(1.0);

		problem.mutable_shift_length_bound()->set_min(3);
		problem.mutable_shift_length_bound()->set_max(5);
	} else {
		Timer t(to_string("Reading file ", args[1]));

		ifstream fin(args[1]);
		check(!!fin, "Could not open ", args[1]);
		string file_content(istreambuf_iterator<char>{fin}, istreambuf_iterator<char>{});
		if (!google::protobuf::TextFormat::ParseFromString(file_content, &problem)) {
			t.fail();
			return 1;
		}
		t.OK();
	}

	ShiftGenerationColgenProblem colgen_problem(problem);
	vector<double> solution_objectives;

	for (int i = 1; i <= FLAGS_num_solutions; ++i) {
		colgen_problem.solve();
		auto elapsed_time = wall_time() - start_time;
		cerr << colgen_problem.fraction_used() * 100.0 << "% of generated columns used.\n";
		cerr << "Colgen done in " << elapsed_time << "s.\n\n\n";

		auto result = colgen_problem.get_solution();
		cout_row(problem);
		vector<int> fulfilled(problem.period_size(), 0);
		for (auto& shift : result) {
			cout_row(problem, [&shift](int i) {
				if (shift.first <= i && i < shift.second) {
					cout << YELLOW << "XXX" << NORMAL;
				} else {
					cout << "";
				}
			});
			cout << endl;
			for (int i : range(shift.first, shift.second)) {
				fulfilled[i] += 1;
			}
		}
		cout_row(problem);
		cout_row(problem, [&problem](int i) { cout << problem.period(i).demand(); });
		cout << " demand.\n";
		cout_row(problem);
		cout_row(problem, [&fulfilled](int i) { cout << fulfilled.at(i); });
		cout << " fullfilled.\n";
		cout_row(problem);
		cout_row(problem, [&problem, &fulfilled](int i) {
			int diff = problem.period(i).demand() - fulfilled.at(i);
			if (diff == 0) {
				cout << GREEN;
			} else {
				cout << RED;
			}
			cout << setw(3) << diff;
			cout << NORMAL;
		});
		cout << " abs. difference.\n";
		cout_row(problem);
		auto solution_value = colgen_problem.integral_solution_value();
		cout << "Solution value: " << solution_value << endl;
		if (FLAGS_verify_solution != -1) {
			check(solution_value <= FLAGS_verify_solution, "Solution is not good enough.");
		}
	}
	return 0;
}

int main(int num_args, char* args[]) {
	if (num_args >= 2 && ends_with(to_string(args[1]), "-help")) {
		gflags::SetUsageMessage(string("Runs colgen method for shift generation.\n   ")
		                        + string(args[0]) + string(" <input_file>"));
		gflags::ShowUsageWithFlagsRestrict(args[0], "colgen");
		return 0;
	}
	return main_runner(main_program, num_args, args);
}
