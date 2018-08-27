// Suggested parameters:
//
// --min_rmp_iterations=30 or larger.
//
// See
// https://github.com/google/or-tools/blob/master/examples/cpp/strawberry_fields_with_column_generation.cc
// for a description of the problem.
#include <fstream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <tuple>
#include <vector>
using namespace std;

#include <absl/strings/substitute.h>
#include <gflags/gflags.h>

#include <minimum/core/check.h>
#include <minimum/core/grid.h>
#include <minimum/core/main.h>
#include <minimum/core/range.h>
#include <minimum/linear/colgen/problem.h>
using namespace minimum::core;
using namespace minimum::linear;

DEFINE_int32(colgen_instance, -1, "Which instance to solve (0 - 9)");

// Code in operations_research is adapted from or-tools, file
// strawberry_fields_with_column_generation.cc
namespace operations_research {

struct Instance {
	int max_boxes;
	int width;
	int height;
	const char* grid;

	char get(int x, int y) const { return grid[pos(x, y)]; }

	constexpr int size() const { return width * height; }

	int pos(int x, int y) const {
		minimum_core_assert(0 <= x && x < width);
		minimum_core_assert(0 <= y && y < height);
		return x + width * y;
	}
	int x(int pos) const {
		minimum_core_assert(0 <= pos && pos < size());
		return pos % width;
	}
	int y(int pos) const {
		minimum_core_assert(0 <= pos && pos < size());
		return pos / width;
	}

	int number_of_strawberries() const {
		return std::count_if(grid, grid + size(), [](char ch) { return ch == '@'; });
	}

	int strawberry_index(int x, int y) const {
		int p = pos(x, y);
		minimum_core_assert(grid[p] == '@');
		return std::count_if(grid, grid + p, [](char ch) { return ch == '@'; });
	}

	void print() const {
		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				std::cout << get(x, y);
			}
			std::cout << "\n";
		}
	}

	void print(const std::vector<int>& solution) const {
		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				int p = pos(x, y);
				if (solution[p] < 0) {
					std::cout << ".";
				} else {
					std::cout << char('A' + solution[p]);
				}
			}
			std::cout << "\n";
		}
	}
};

Instance kInstances[] = {{4,
                          22,
                          6,
                          "..@@@@@..............."
                          "..@@@@@@........@@@..."
                          ".....@@@@@......@@@..."
                          ".......@@@@@@@@@@@@..."
                          ".........@@@@@........"
                          ".........@@@@@........"},
                         {3,
                          13,
                          10,
                          "............."
                          "............."
                          "............."
                          "...@@@@......"
                          "...@@@@......"
                          "...@@@@......"
                          ".......@@@..."
                          ".......@@@..."
                          ".......@@@..."
                          "............."},
                         {4,
                          13,
                          9,
                          "............."
                          "..@.@.@......"
                          "...@.@.@....."
                          "..@.@.@......"
                          "..@.@.@......"
                          "...@.@.@....."
                          "....@.@......"
                          "..........@@@"
                          "..........@@@"},
                         {4,
                          13,
                          9,
                          ".........@..."
                          ".........@..."
                          "@@@@@@@@@@..."
                          "..@......@..."
                          "..@......@..."
                          "..@......@..."
                          "..@@@@@@@@@@@"
                          "..@.........."
                          "..@.........."},
                         {7,
                          25,
                          14,
                          "........................."
                          "..@@@@@@@@@@@@@@@@@@@@..."
                          "..@@@@@@@@@@@@@@@@@@@@..."
                          "..@@.................@..."
                          "..@@.................@..."
                          "..@@.......@@@.......@.@."
                          "..@@.......@@@.......@..."
                          "..@@...@@@@@@@@@@@@@@@..."
                          "..@@...@@@@@@@@@@@@@@@..."
                          "..@@.......@@@.......@..."
                          "..@@.......@@@.......@..."
                          "..@@.................@..."
                          "..@@.................@..."
                          "........................."},
                         {6,
                          25,
                          16,
                          "........................."
                          "......@@@@@@@@@@@@@......"
                          "........................."
                          ".....@..........@........"
                          ".....@..........@........"
                          ".....@......@............"
                          ".....@......@.@@@@@@@...."
                          ".....@......@............"
                          ".....@......@.@@@@@@@...."
                          ".....@......@............"
                          "....@@@@....@............"
                          "....@@@@....@............"
                          "..@@@@@@....@............"
                          "..@@@.......@............"
                          "..@@@...................."
                          "..@@@@@@@@@@@@@@@@@@@@@@@"},
                         {5,
                          40,
                          18,
                          "........................................"
                          "........................................"
                          "...@@@@@@..............................."
                          "...@@@@@@..............................."
                          "...@@@@@@..............................."
                          "...@@@@@@.........@@@@@@@@@@............"
                          "...@@@@@@.........@@@@@@@@@@............"
                          "..................@@@@@@@@@@............"
                          "..................@@@@@@@@@@............"
                          ".............@@@@@@@@@@@@@@@............"
                          ".............@@@@@@@@@@@@@@@............"
                          "........@@@@@@@@@@@@...................."
                          "........@@@@@@@@@@@@...................."
                          "........@@@@@@.........................."
                          "........@@@@@@.........................."
                          "........................................"
                          "........................................"
                          "........................................"},
                         {8,
                          40,
                          18,
                          "........................................"
                          "..@@.@.@.@.............................."
                          "..@@.@.@.@...............@.............."
                          "..@@.@.@.@............@................."
                          "..@@.@.@.@.............................."
                          "..@@.@.@.@.................@............"
                          "..@@.@..................@..............."
                          "..@@.@.................................."
                          "..@@.@.................................."
                          "..@@.@................@@@@.............."
                          "..@@.@..............@@@@@@@@............"
                          "..@@.@.................................."
                          "..@@.@..............@@@@@@@@............"
                          "..@@.@.................................."
                          "..@@.@................@@@@.............."
                          "..@@.@.................................."
                          "..@@.@.................................."
                          "........................................"},
                         {10,
                          40,
                          19,
                          "@@@@@..................................."
                          "@@@@@..................................."
                          "@@@@@..................................."
                          "@@@@@..................................."
                          "@@@@@..................................."
                          "@@@@@...........@@@@@@@@@@@............."
                          "@@@@@...........@@@@@@@@@@@............."
                          "....................@@@@................"
                          "....................@@@@................"
                          "....................@@@@................"
                          "....................@@@@................"
                          "....................@@@@................"
                          "...............@@@@@@@@@@@@@@..........."
                          "...............@@@@@@@@@@@@@@..........."
                          ".......@@@@@@@@@@@@@@@@@@@@@@..........."
                          ".......@@@@@@@@@........................"
                          "........................................"
                          "........................................"
                          "........................................"},
                         {10,
                          40,
                          25,
                          "...................@...................."
                          "...............@@@@@@@@@................"
                          "............@@@.........@@@............."
                          "...........@...............@............"
                          "..........@.................@..........."
                          ".........@...................@.........."
                          ".........@...................@.........."
                          ".........@.....@@......@@....@.........."
                          "........@.....@@@@....@@@@....@........."
                          "........@.....................@........."
                          "........@.....................@........."
                          "........@..........@@.........@........."
                          ".......@@..........@@.........@@........"
                          "........@.....................@........."
                          "........@.....................@........."
                          "........@......@@@@@@@@@......@........."
                          "........@......@@@@@@@@@......@........."
                          ".........@...................@.........."
                          ".........@...................@.........."
                          ".........@...................@.........."
                          "..........@.................@..........."
                          "...........@...............@............"
                          "............@@@.........@@@............."
                          "...............@@@@@@@@@................"
                          "...................@...................."}};

const int kInstanceCount = 10;

class Box {
   public:
	static const int kAreaCost = 1;
	static const int kFixedCost = 10;

	Box() {}
	Box(int x_min, int x_max, int y_min, int y_max)
	    : x_min_(x_min), x_max_(x_max), y_min_(y_min), y_max_(y_max) {
		minimum_core_assert(x_max >= x_min);
		minimum_core_assert(y_max >= y_min);
	}

	int x_min() const { return x_min_; }
	int x_max() const { return x_max_; }
	int y_min() const { return y_min_; }
	int y_max() const { return y_max_; }

	// Lexicographic order
	int Compare(const Box& box) const {
		int c;
		if ((c = (x_min() - box.x_min())) != 0) return c;
		if ((c = (x_max() - box.x_max())) != 0) return c;
		if ((c = (y_min() - box.y_min())) != 0) return c;
		return y_max() - box.y_max();
	}

	bool Contains(int x, int y) const {
		return x >= x_min() && x <= x_max() && y >= y_min() && y <= y_max();
	}

	int Cost() const {
		return kAreaCost * (x_max() - x_min() + 1) * (y_max() - y_min() + 1) + kFixedCost;
	}

	std::string DebugString() const {
		return absl::Substitute("[$0,$1x$2,$3]c$4", x_min(), y_min(), x_max(), y_max(), Cost());
	}

   private:
	int x_min_;
	int x_max_;
	int y_min_;
	int y_max_;
};

struct BoxLessThan {
	bool operator()(const Box& b1, const Box& b2) const { return b1.Compare(b2) < 0; }
};
}  // namespace operations_research

class CoveringProblem : public colgen::Problem {
   public:
	CoveringProblem(const operations_research::Instance& instance_)
	    : Problem(number_of_rows(instance_)), instance(instance_) {
		// Add a huge column covering everything.
		colgen::Column large =
		    make_column(operations_research::Box{0, instance.width - 1, 0, instance.height - 1});
		pool.add(move(large));

		for (int i = 0; i < instance.size(); ++i) {
			if (instance.get(instance.x(i), instance.y(i)) == '@') {
				set_row_lower_bound(i, 1);
			} else {
				set_row_lower_bound(i, 0);
			}
			set_row_upper_bound(i, 1);
		}
		set_row_upper_bound(instance.size(), instance.max_boxes);
	}

	static int number_of_rows(const operations_research::Instance& instance) {
		return instance.size() + 1;  // + 1 is for the limit on the number of houses.
	}

	void generate(const std::vector<double>& dual_variables) override {
		int added = 0;
		const double max_boxes_dual = dual_variables.at(instance.size());

		auto integral_image =
		    make_grid<double>(instance.width + 1, instance.height + 1, []() { return 0; });
		for (int x = 0; x < instance.width; ++x) {
			for (int y = 0; y < instance.height; ++y) {
				integral_image[x + 1][y + 1] = dual_variables.at(instance.pos(x, y))
				                               + integral_image[x][y + 1] + integral_image[x + 1][y]
				                               - integral_image[x][y];
			}
		}

		for (int y_min = 0; y_min < instance.height; ++y_min) {
			for (int x_min = 0; x_min < instance.width; ++x_min) {
				for (int y_max = y_min; y_max < instance.height; ++y_max) {
					for (int x_max = x_min; x_max < instance.width; ++x_max) {
						operations_research::Box box(x_min, x_max, y_min, y_max);

						double cell_coverage_dual =
						    integral_image[x_max + 1][y_max + 1] - integral_image[x_min][y_max + 1]
						    - integral_image[x_max + 1][y_min] + integral_image[x_min][y_min];

						// All coefficients for new column are 1, so no need to
						// multiply constraint duals by any coefficients when
						// computing the reduced cost.
						const double reduced_cost =
						    box.Cost() - (cell_coverage_dual + max_boxes_dual);

						if (reduced_cost < 0) {
							auto col = make_column(box);
							minimum_core_assert(
							    abs(reduced_cost - col.reduced_cost(dual_variables)) < 1e-6,
							    reduced_cost,
							    " =/= ",
							    col.reduced_cost(dual_variables));
							pool.add(move(col));
							added++;
							if (added > 100) {
								return;
							}
						}
					}
				}
			}
		}
	}

	auto compute_solution() {
		std::vector<int> result(instance.size(), -1);
		int piece_index = 0;
		for (auto i : active_columns()) {
			auto& column = pool.at(i);
			if (column.solution_value > 0.5) {
				for (auto& entry : column) {
					if (entry.row < instance.size()) {
						result[entry.row] = piece_index;
					}
				}
				piece_index++;
			}
		}
		return result;
	}

	int fix(const FixInformation& information) override { return 0; }

   private:
	colgen::Column make_column(const operations_research::Box& box) const {
		colgen::Column column(box.Cost(), 0, 1);
		for (int x = box.x_min(); x <= box.x_max(); ++x) {
			for (int y = box.y_min(); y <= box.y_max(); ++y) {
				column.add_coefficient(instance.pos(x, y), 1);
			}
		}
		column.add_coefficient(instance.size(), 1);
		return column;
	}

	const operations_research::Instance& instance;
};

void solve_instance(const operations_research::Instance& instance) {
	CoveringProblem problem(instance);
	double objective_value = problem.solve();
	auto solution = problem.compute_solution();
	instance.print();
	instance.print(solution);
	double worst_solution = operations_research::Box::kFixedCost
	                        + operations_research::Box::kAreaCost * instance.size();
	std::cout << "Worst solution is " << worst_solution << ", found " << objective_value << "\n";
	check(objective_value < worst_solution - 1, "Solution not good enough.");
}

int main_program(int num_args, char* args[]) {
	using namespace std;

	if (FLAGS_colgen_instance == -1) {
		for (int i = 0; i < operations_research::kInstanceCount; ++i) {
			const operations_research::Instance& instance = operations_research::kInstances[i];
			solve_instance(instance);
		}
	} else {
		check(FLAGS_colgen_instance >= 0
		          && FLAGS_colgen_instance < operations_research::kInstanceCount,
		      "Invalid colgen instance.");
		const operations_research::Instance& instance =
		    operations_research::kInstances[FLAGS_colgen_instance];
		solve_instance(instance);
	}

	return 0;
}

int main(int num_args, char* args[]) { return main_runner(main_program, num_args, args); }
