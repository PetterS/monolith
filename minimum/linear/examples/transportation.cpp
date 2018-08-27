// Petter Strandmark 2013
// petter.strandmark@gmail.com
//
// Dantzig, G B, chapter 3.3 in Linear Programming and Extensions,
// Princeton University Press, Princeton, New Jersey, 1963.
// http://en.wikipedia.org/wiki/AMPL

#include <iostream>
#include <map>
#include <set>
#include <string>
using namespace std;

#include <minimum/core/check.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/solver.h>
using namespace minimum::linear;

// In order to have two-dimensional "arrays" (maps).
pair<string, string> operator,(string p1, string p2) { return make_pair(p1, p2); }

void main_program() {
	IP lp;

	typedef string Place;
	Place seattle = "Seattle";
	Place san_diego = "San Diego";
	Place new_york = "New York";
	Place chicago = "Chicago";
	Place topeka = "Topeka";

	set<Place> all_places;
	all_places.insert(seattle);
	all_places.insert(san_diego);
	all_places.insert(new_york);
	all_places.insert(chicago);
	all_places.insert(topeka);

	map<Place, double> market_demand;
	market_demand[new_york] = 325;
	market_demand[chicago] = 300;
	market_demand[topeka] = 275;

	map<Place, double> plant_capacity;
	plant_capacity[seattle] = 350;
	plant_capacity[san_diego] = 600;

	map<pair<Place, Place>, double> distance;
	distance[seattle, new_york] = 2.5;
	distance[san_diego, new_york] = 2.5;

	distance[seattle, chicago] = 1.7;
	distance[san_diego, chicago] = 1.8;

	distance[seattle, topeka] = 2.5;
	distance[san_diego, topeka] = 2.5;

	map<pair<Place, Place>, Variable> shipment;

	Sum cost = 0;
	for (auto p : all_places) {
		for (auto m : all_places) {
			shipment[p, m] = lp.add_variable(IP::Real);
			lp.add_bounds(0, shipment[p, m], 1e10);

			cost += shipment[p, m] * 90 * distance[p, m] / 1000;
		}
	}

	// Observe supply limit at p.
	for (auto p : all_places) {
		Sum leaving_plant = 0;
		for (auto m : all_places) {
			leaving_plant += shipment[p, m];
		}
		lp.add_constraint(leaving_plant <= plant_capacity[p]);
	}

	// Satisfy demand at market m.
	for (auto m : all_places) {
		Sum arriving_to_market = 0;
		for (auto p : all_places) {
			arriving_to_market += shipment[p, m];
		}
		lp.add_constraint(arriving_to_market >= market_demand[m]);
	}

	// Solve linear program.
	solve(&lp);

	cout << endl << endl;
	for (auto p : all_places) {
		for (auto m : all_places) {
			if (shipment[p, m].value() > 0.1) {
				std::cout << p << " ships " << shipment[p, m] << " units to " << m << std::endl;
			}
		}
	}
	cout << "Cost is " << cost << endl;
}

int main() {
	try {
		main_program();
	} catch (std::exception& err) {
		std::cerr << "ERROR: " << err.what() << std::endl;
		return 1;
	}
}
