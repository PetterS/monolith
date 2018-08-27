//
// A solver for Shikaku puzzles.
// http://en.wikipedia.org/wiki/Shikaku
//

#include <iomanip>
#include <iostream>
#include <vector>

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/solver.h>
using namespace minimum::core;
using namespace minimum::linear;

int main_program() {
	using namespace std;

	const int _ = 0;
	const vector<vector<int>> given = {{9, _, _, _, 12, _, _, 5, _, _},
	                                   {_, _, _, _, _, _, _, _, _, _},
	                                   {_, _, _, _, _, _, _, _, _, 6},
	                                   {8, _, 6, _, 8, _, _, _, _, _},
	                                   {_, _, _, _, _, _, _, _, _, _},
	                                   {_, _, _, _, _, _, _, _, _, _},
	                                   {_, _, _, _, _, 6, _, 8, _, 12},
	                                   {4, _, _, _, _, _, _, _, _, _},
	                                   {_, _, _, _, _, _, _, _, _, _},
	                                   {_, _, 3, _, _, 9, _, _, _, 4}};

	int n = given.size();

	struct Rect {
		int i1;
		int i2;
		int j1;
		int j2;

		int area() { return (i2 - i1) * (j2 - j1); }
		bool contains(int i, int j) { return i1 <= i && i < i2 && j1 <= j && j < j2; }
	};

	IP ip;
	vector<Rect> rectangles;
	vector<BooleanVariable> variables;

	// Generate all possible rectangles. A brute-force search
	// is fine here and the obvious optimizations not needed.
	int num_rects = 0;
	for (int i : range(n)) {
		for (int j : range(n)) {
			if (given[i][j] != 0) {
				auto area = given[i][j];
				// Generate all rectangles here.
				for (int i1 : range(n)) {
					for (int j1 : range(n)) {
						for (int i2 : range(i1 + 1, n + 1)) {
							for (int j2 : range(j1 + 1, n + 1)) {
								Rect rect{i1, i2, j1, j2};
								if (rect.contains(i, j) && rect.area() == area) {
									num_rects++;
									rectangles.emplace_back(rect);
									variables.emplace_back(ip.add_boolean());
								}
							}
						}
					}
				}
			}
		}
	}
	cout << "Found " << num_rects << " rectangles.\n";

	// The only constraints we have is that each grid square is contained
	// in exactly one rectangle. We already took care of the area
	// requirement when we generated them.
	for (int i : range(n)) {
		for (int j : range(n)) {
			Sum ij_sum = 0;
			for (auto r : range(rectangles.size())) {
				if (rectangles[r].contains(i, j)) {
					ij_sum += variables.at(r);
				}
			}
			ip.add_constraint(ij_sum == 1);
		}
	}

	IPSolver solver;
	solver.set_silent(true);
	auto solutions = solver.solutions(&ip);
	while (solutions->get()) {
		// Print the solution (a bit messy to generate an ascii table).

		vector<vector<char>> vertical_border(n, vector<char>(n + 1, ' '));
		vector<string> horizontal_border(n + 1, string(n, ' '));
		horizontal_border[0] = string(n, '-');
		horizontal_border[n] = string(n, '-');
		for (int i : range(n)) {
			vertical_border[i][0] = '|';
			vertical_border[i][n] = '|';
		}

		// For each rectangle in the solution, add its borders to the
		// vectors of borders to print.
		for (auto r : range(rectangles.size())) {
			if (variables.at(r).bool_value()) {
				auto& rect = rectangles.at(r);
				for (int i : range(rect.i1, rect.i2)) {
					vertical_border[i][rect.j1] = '|';
					vertical_border[i][rect.j2] = '|';
				}
				for (int j : range(rect.j1, rect.j2)) {
					horizontal_border[rect.i1][j] = '-';
					horizontal_border[rect.i2][j] = '-';
				}
			}
		}

		// Print all given numbers separated with borders. In each empty grid
		// space, print a dot.
		cout << "\nSolution:\n";
		cout << setw(n * 3 + 1) << setfill('-') << "" << setfill(' ') << "\n";
		for (int i : range(n)) {
			cout << vertical_border[i][0];
			for (int j : range(n)) {
				cout << setw(2) << (given[i][j] != 0 ? to_string(given[i][j]) : ".")
				     << vertical_border[i][j + 1];
			}
			cout << "\n";
			cout << vertical_border[i][0];
			for (int j : range(n)) {
				auto cpy = horizontal_border[i + 1][j];
				cout << setw(2) << setfill(cpy) << "" << setfill(' ');
				if (cpy == '-') {
					cout << cpy;
				} else {
					cout << vertical_border[i][j + 1];
				}
			}
			cout << "\n";
		}
	}

	return 0;
}

int main() {
	try {
		return main_program();
	} catch (std::exception& err) {
		std::cerr << "ERROR: " << err.what() << std::endl;
		return 1;
	}
}
