// Petter Strandmark 2013
// petter.strandmark@gmail.com

#include <iostream>
#include <random>
using namespace std;

#include <minimum/core/check.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/solver.h>
using namespace minimum::linear;

int main() {
	const int K = 10;  // Number of passengers.
	const int S = 5;   // Number of seats.
	const int T = 10;  // Number of parts into which the distance is divided.

	// Generate random bookings.
	mt19937_64 rng(unsigned(0));
	uniform_int_distribution<int> booking(0, T - 1);
	typedef pair<int, int> Booking;
	vector<Booking> all_bookings(K);
	for (int k = 0; k < K; ++k) {
		auto& Tk = all_bookings[k];
		auto start = booking(rng);
		auto end = booking(rng);
		if (start > end) {
			swap(start, end);
		}
		Tk = make_pair(start, end);
		cout << "Passenger " << k << " booked stations " << Tk.first << "--" << Tk.second << endl;
	}
	cout << endl;

	IP ip;
	auto x = ip.add_boolean_cube(K, S, T);

	// Respect bookings.
	for (int k = 0; k < K; ++k) {
		auto Tk = all_bookings[k];

		for (int s = 0; s < S; ++s) {
			for (int t = 0; t < Tk.first; ++t) {
				ip.add_bounds(0, x[k][s][t], 0);
			}
			for (int t = Tk.second + 1; t < T; ++t) {
				ip.add_bounds(0, x[k][s][t], 0);
			}
		}
	}

	// Passengers only take up one seat.
	for (int k = 0; k < K; ++k) {
		for (int t = 0; t < T; ++t) {
			Sum sum = 0;
			for (int s = 0; s < S; ++s) {
				sum += x[k][s][t];
			}
			ip.add_constraint(sum <= 1);
		}
	}

	// One passenger per seat.
	for (int s = 0; s < S; ++s) {
		for (int t = 0; t < T; ++t) {
			Sum sum = 0;
			for (int k = 0; k < K; ++k) {
				sum += x[k][s][t];
			}
			ip.add_constraint(sum <= 1);
		}
	}

	// Same seat throughout the journey.
	for (int k = 0; k < K; ++k) {
		auto Tk = all_bookings[k];

		for (int s = 0; s < S; ++s) {
			for (int t = Tk.first + 1; t <= Tk.second; ++t) {
				ip.add_constraint(x[k][s][Tk.first] == x[k][s][t]);
			}
		}
	}

	// The objective is to maximize bookings.
	Sum objective = 0;
	for (int k = 0; k < K; ++k) {
		for (int t = 0; t < T; ++t) {
			for (int s = 0; s < S; ++s) {
				objective += x[k][s][t];
			}
		}
	}
	// Negative because we want to maximize.
	ip.add_objective(-objective);

	cout << "Solving program..." << endl;
	solve(&ip);

	for (int s = 0; s < S; ++s) {
		cout << "Seat " << s << ": [";
		for (int t = 0; t < T; ++t) {
			bool found = false;
			for (int k = 0; k < K; ++k) {
				if (x[k][s][t].bool_value()) {
					minimum::core::check(!found);
					cout << k % 10;  // Will always print 0--9.
					found = true;
				}
			}
			if (!found) {
				cout << "-";
			}
		}
		cout << "]" << endl;
	}
}
