
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>
using namespace std;

#include <minimum/algorithms/knapsack.h>
#include <minimum/core/check.h>
#include <minimum/core/numeric.h>
#include <minimum/core/string.h>
#include <minimum/linear/colgen/problem.h>
#include <minimum/linear/data/util.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/solver.h>
using namespace minimum::core;
using namespace minimum::linear;
using namespace minimum::linear::colgen;

int roll_width = 5600;

struct Requirement {
	Requirement(int w, int a) : width(w), amount(a) {}
	int width = 0;
	int amount = 0;
};
std::vector<Requirement> requirements = {{1380, 22},
                                         {1520, 25},
                                         {1560, 12},
                                         {1710, 14},
                                         {1820, 18},
                                         {1880, 18},
                                         {1930, 20},
                                         {2000, 10},
                                         {2050, 12},
                                         {2100, 14},
                                         {2140, 16},
                                         {2150, 18},
                                         {2200, 20}};

class CuttingStockProblem : public Problem {
   public:
	CuttingStockProblem(int roll_width_, vector<Requirement> requirements_)
	    : Problem(requirements_.size()),
	      roll_width(roll_width_),
	      requirements(requirements_),
	      rng(1u) {
		// Add trivial columns. One for each roll.
		for (int r = 0; r < requirements.size(); ++r) {
			Column column(1.0, 0, 1e100);
			column.add_coefficient(r, 1.0);
			pool.add(move(column));
		}

		for (int i = 0; i < requirements.size(); ++i) {
			set_row_lower_bound(i, requirements[i].amount);
		}

		for (int i = 0; i < requirements.size(); ++i) {
			widths.push_back(requirements[i].width);
		}
	}

	virtual void generate(const std::vector<double>& dual_variables_in) override {
		auto dual_variables = dual_variables_in;
		minimum_core_assert(dual_variables.size() == requirements.size());

		// Solve knapsack problem.
		int added = 0;
		while (added < 10) {
			vector<int> solution;
			minimum::algorithms::solve_knapsack(roll_width, widths, dual_variables, &solution);

			Column column(1.0, 0, 1e100);
			for (int i = 0; i < requirements.size(); ++i) {
				if (solution[i] > 0) {
					column.add_coefficient(i, solution[i]);
					dual_variables[i] = -1e10;
				}
			}

			if (column.reduced_cost(dual_variables_in) < -1e-6) {
				pool.add(move(column));
				added++;
			} else {
				break;
			}
		}
	}

	virtual int fix(const FixInformation& information) override {
		int fixed = 0;
		if (information.objective_change >= -0.05 && information.objective_change <= 0) {
			vector<int> indices(active_columns().size());
			iota(indices.begin(), indices.end(), 0);

			vector<int> widths;
			for (auto i : active_columns()) {
				auto& column = pool.at(i);
				int width = 0;
				for (auto& entry : column) {
					width += entry.coef * requirements[entry.row].width;
				}
				widths.push_back(width);
			}

			sort(indices.begin(), indices.end(), [&widths](int i1, int i2) {
				return widths[i1] > widths[i2];
			});

			for (int i : indices) {
				auto& column = pool.at(active_columns().at(i));

				if (integer_residual(column.solution_value) > 1e-5) {
					column.fix(int(column.solution_value + 0.5));

					if (false) {
						cerr << "-- Fixing [";
						bool first = true;
						for (auto& entry : column) {
							if (!first) {
								cerr << " ";
							}
							first = false;
							if (entry.coef != 1) {
								cerr << entry.coef << "x";
							}
							cerr << requirements[entry.row].width;
						}
						cerr << "] = " << int(column.solution_value + 0.5) << "\n";
					}

					fixed++;
					break;
				}
			}
		}
		return fixed;
	}

   private:
	const int roll_width;
	const vector<Requirement> requirements;
	vector<int> widths;

	mt19937_64 rng;
};

double run_colgen(int roll_width, const vector<Requirement>& requirements) {
	CuttingStockProblem problem(roll_width, requirements);

	double best_num_rolls = numeric_limits<double>::max();
	for (int iter = 1; iter <= 2; ++iter) {
		check(problem.solve(), "Colgen failed.");

		double num_rolls = 0;
		for (auto i : problem.active_columns()) {
			auto& column = problem.column_at(i);

			if (column.solution_value > 0) {
				vector<int> widths;
				for (auto& entry : column) {
					for (int i = 0; i < int(entry.coef + 0.5); ++i) {
						widths.push_back(requirements[entry.row].width);
					}
				}
				cout << to_string(widths) << ": " << column.solution_value << "\n";
				num_rolls += column.solution_value;
			}
		}
		cout << "Total of " << to_string(num_rolls) << " rolls produced.\n";
		best_num_rolls = min(best_num_rolls, num_rolls);
	}

	return best_num_rolls;
}

double run_direct() {
	IP ip;
	vector<Sum> roll_constraints(requirements.size());
	vector<vector<int>> patterns;
	vector<Variable> produced;

	// Add all patterns with three rolls.
	for (int i = 0; i < requirements.size(); ++i) {
		for (int j = i; j < requirements.size(); ++j) {
			for (int k = j; k < requirements.size(); ++k) {
				if (requirements[i].width + requirements[j].width + requirements[k].width
				    <= roll_width) {
					patterns.emplace_back(vector<int>{i, j, k});
					produced.emplace_back(ip.add_variable(IP::Integer));
					ip.add_constraint(produced.back() >= 0);
					roll_constraints[i] += produced.back();
					roll_constraints[j] += produced.back();
					roll_constraints[k] += produced.back();
					ip.add_objective(produced.back());
				}
			}
		}
	}

	// Add one pattern with four rolls.
	if (4 * requirements[0].width <= roll_width) {
		patterns.emplace_back(vector<int>{0, 0, 0, 0});
		produced.emplace_back(ip.add_variable(IP::Integer));
		ip.add_constraint(produced.back() >= 0);
		roll_constraints[0] += 4 * produced.back();
		ip.add_objective(produced.back());
	}

	cout << patterns.size() << " patterns.\n";

	for (int r = 0; r < requirements.size(); ++r) {
		ip.add_constraint(roll_constraints[r] >= requirements[r].amount);
	}

	IPSolver solver;
	solver.set_silent(true);

	auto print_solution = [&]() {
		double num_rolls = 0;
		for (int p = 0; p < patterns.size(); ++p) {
			if (produced[p].value() > 0) {
				vector<int> widths;
				for (int i : patterns[p]) {
					widths.push_back(requirements[i].width);
				}
				cout << to_string(widths) << ": " << produced[p].value() << "\n";
				num_rolls += produced[p].value();
			}
		}
		cout << "Total of " << to_string(num_rolls) << " rolls produced.\n";
	};

	cout << "LP solution:\n";
	minimum_core_assert(solver.solve_relaxation(&ip));
	print_solution();
	cout << "\n";

	cout << "IP solution:\n";
	minimum_core_assert(solver.solutions(&ip)->get());
	print_solution();

	return ip.get_entire_objective();
}

int main(int argc, char* argv[]) {
	try {
		bool run_default = false;
		if (argc >= 2 && string(argv[1]) == "default") {
			run_default = true;
		}

		if (run_default) {
			auto direct = run_direct();
			cout << "\nColgen\n";
			auto colgen = run_colgen(roll_width, requirements);
			// We can probably expect optimality for this simple problem.
			check(abs(direct - colgen) < 1e-6, "Colgen solution not optimal.");
		} else {
			ifstream fin(data::get_directory() + "/cutting_stock/hard28.txt");
			vector<pair<string, double>> results;
			while (true) {
				string problem_name;
				getline(fin, problem_name);
				cout << "Running " << problem_name << "\n";
				int n;
				fin >> n >> roll_width;
				requirements.clear();
				for (int i = 0; i < n; ++i) {
					Requirement req(0, 0);
					fin >> req.width >> req.amount;
					requirements.push_back(req);
				}
				string s;
				getline(fin, s);
				if (!fin) {
					break;
				}

				cout << "\nColgen\n";
				auto result = run_colgen(roll_width, requirements);
				results.emplace_back(problem_name, result);

				if (results.size() > 10000) {
					break;
				}
			}

			for (auto& entry : results) {
				cout << setw(20) << left << entry.first << ": " << entry.second << "\n";
			}
		}
	} catch (std::exception& err) {
		std::cerr << "ERROR: " << err.what() << std::endl;
		return 1;
	}
}
