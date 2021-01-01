// Petter Strandmark
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
	auto x = ip.add_boolean_grid(K, S);

	// Passengers only take up one seat.
	for (int k = 0; k < K; ++k) {
		Sum sum = 0;
		for (int s = 0; s < S; ++s) {
			sum += x[k][s];
		}
		ip.add_constraint(sum <= 1);
	}

	// One passenger per seat.
	for (int k1 = 0; k1 < K; ++k1) {
		auto Tk1 = all_bookings[k1];
		for (int k2 = k1 + 1; k2 < K; ++k2) {
			auto Tk2 = all_bookings[k2];
			if ((Tk1.first <= Tk2.first && Tk1.second >= Tk2.first)
			    || (Tk1.first >= Tk2.first && Tk1.first <= Tk2.second)) {
				// Passengers overlap. They can not be assigned the
				// same seat.
				for (int s = 0; s < S; ++s) {
					ip.add_constraint(x[k1][s] + x[k2][s] <= 1);
				}
			}
		}
	}

	// The objective is to maximize bookings.
	Sum objective = 0;
	for (int k = 0; k < K; ++k) {
		auto Tk = all_bookings[k];
		for (int s = 0; s < S; ++s) {
			objective += (Tk.second - Tk.first + 1) * x[k][s];
		}
	}
	// Negative because we want to maximize.
	ip.add_objective(-objective);

	cout << "Solving program..." << endl;
	solve(&ip);

	for (int s = 0; s < S; ++s) {
		vector<char> bookings(T, '-');
		for (int k = 0; k < K; ++k) {
			auto Tk = all_bookings[k];
			if (x[k][s].bool_value()) {
				for (int t = Tk.first; t <= Tk.second; ++t) {
					minimum::core::check(bookings[t] == '-');
					bookings[t] = '0' + k % 10;
				}
			}
		}

		cout << "Seat " << s << ": [";
		for (int t = 0; t < T; ++t) {
			cout << bookings[t];
		}
		cout << "]" << endl;
	}
}
