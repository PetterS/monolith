#include <iostream>
#include <random>
using namespace std;

#include <minimum/core/check.h>
#include <minimum/core/string.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/solver.h>
using namespace minimum::linear;

int main() {
	auto I = {0, 1, 2, 3};
	auto J = {0, 1, 2, 3};

	vector<vector<int>> s0 = {{1, 0, 0, 0},   //
	                          {0, 0, 0, 0},   //
	                          {0, 1, 0, 0},   //
	                          {0, 0, 1, 0}};  //

	IP ip;

	auto s = ip.add_boolean_grid(4, 4);
	auto y = ip.add_grid(4, 4, IP::Integer);
	Sum n = 0;

	for (auto i : I) {
		for (auto j : J) {
			Sum should_be_even = s0[i][j];
			for (auto i2 : I) {
				should_be_even += s[i2][j];
			}
			for (auto j2 : J) {
				should_be_even += s[i][j2];
			}
			should_be_even -= s[i][j];

			ip.add_constraint(2 * y[i][j] == should_be_even);

			n += s[i][j];
		}
	}

	ip.add_objective(n);
	solve(&ip);

	cout << "Solution in " << n << " moves:\n";
	for (auto i : I) {
		cout << minimum::core::to_string(s[i]) << '\n';
	}
}
