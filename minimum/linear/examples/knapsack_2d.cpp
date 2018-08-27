// Two-dimensional knapsack
// http://pyeasyga.readthedocs.org/en/latest/examples.html

#include <iostream>
using namespace std;

#include <minimum/core/check.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/solver.h>
using namespace minimum::core;
using namespace minimum::linear;

int main_program() {
	struct Data {
		double weight;
		double volume;
		double price;
	};
	vector<Data> data = {{821, 0.8, 118}, {1144, 1, 322},   {634, 0.7, 166},  {701, 0.9, 195},
	                     {291, 0.9, 100}, {1702, 0.8, 142}, {1633, 0.7, 100}, {1086, 0.6, 145},
	                     {124, 0.6, 100}, {718, 0.9, 208},  {976, 0.6, 100},  {1438, 0.7, 312},
	                     {910, 1, 198},   {148, 0.7, 171},  {1636, 0.9, 117}, {237, 0.6, 100},
	                     {771, 0.9, 329}, {604, 0.6, 391},  {1078, 0.6, 100}, {640, 0.8, 120},
	                     {1510, 1, 188},  {741, 0.6, 271},  {1358, 0.9, 334}, {1682, 0.7, 153},
	                     {993, 0.7, 130}, {99, 0.7, 100},   {1068, 0.8, 154}, {1669, 1, 289}};

	IP ip;
	auto x = ip.add_boolean_vector(data.size());
	Sum weight = 0;
	Sum volume = 0;
	Sum price = 0;
	for (int i = 0; i < data.size(); ++i) {
		weight += x[i] * data[i].weight;
		volume += x[i] * data[i].volume;
		price += x[i] * data[i].price;
	}
	ip.add_constraint(weight <= 12210);
	ip.add_constraint(volume <= 12);
	ip.add_objective(-price);

	check(solve(&ip), "Could not solve problem.");

	cout << "Weight: " << weight.value() << "\n";
	cout << "Volume: " << volume.value() << "\n";
	cout << "Price : " << price.value() << "\n";
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
