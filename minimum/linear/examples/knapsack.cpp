// Petter Strandmark 2013
// petter.strandmark@gmail.com

#include <iostream>
using namespace std;

#include <minimum/linear/ip.h>
#include <minimum/linear/solver.h>
using namespace minimum::linear;

int main() {
	IP ip;

	double weights_[12] = {332, 980, 158, 422, 869, 765, 529, 786, 40, 268, 910, 269};
	vector<double> weight(weights_, weights_ + 12);

	double values_[12] = {30, 31, 20, 79, 96, 69, 88, 31, 84, 2, 24, 14};
	vector<double> value(values_, values_ + 12);

	auto amount = ip.add_boolean_vector(12);

	Sum total_weight = 0;
	for (int i = 0; i < 12; ++i) {
		ip.add_objective(-value[i] * amount[i]);  // Negative because we want to maximize.
		total_weight += weight[i] * amount[i];
	}

	ip.add_constraint(total_weight <= 3000);

	solve(&ip);

	for (int i = 0; i < 12; ++i) {
		if (amount[i].bool_value()) {
			cout << "Item #" << i + 1 << endl;
		}
	}
	cout << "Total weight: " << total_weight << endl;
}
