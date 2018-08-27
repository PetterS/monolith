// Written by Petter Strandmark a long time ago (way before he knew
// about integer or constraint programming).

#include <ctime>
#include <iostream>
#include <vector>
using namespace std;

#include <minimum/core/time.h>
using namespace minimum::core;

//
// Assumptions:   1 <= L <= x_i <= U
//
class IntegerFeasibility {
   public:
	IntegerFeasibility(int _n)
	    : n(_n),
	      minval(1),
	      maxval(_n),
	      point_to_constraint(n),
	      all_different(false),
	      used(n, false),
	      V(n, 0) {}

	void set_L(int L) { minval = L; }

	void set_U(int U) {
		maxval = U;
		used.resize(U, false);
	}

	void add_equality_constraint(int rhs, int i) {
		int ind[] = {i};
		add_equality_constraint(rhs, vector<int>(ind, ind + 1));
	}
	void add_equality_constraint(int rhs, int i, int j) {
		int ind[] = {i, j};
		add_equality_constraint(rhs, vector<int>(ind, ind + 2));
	}
	void add_equality_constraint(int rhs, int i, int j, int k) {
		int ind[] = {i, j, k};
		add_equality_constraint(rhs, vector<int>(ind, ind + 3));
	}
	void add_equality_constraint(int rhs, int i, int j, int k, int l) {
		int ind[] = {i, j, k, l};
		add_equality_constraint(rhs, vector<int>(ind, ind + 4));
	}
	void add_equality_constraint(int rhs, int i, int j, int k, int l, int m) {
		int ind[] = {i, j, k, l, m};
		add_equality_constraint(rhs, vector<int>(ind, ind + 5));
	}
	void add_equality_constraint(int rhs, const vector<int>& ind) {
		vector<int> coef(ind.size(), 1);
		add_equality_constraint(rhs, ind, coef);
	}
	void add_equality_constraint(int rhs, const vector<int>& ind, const vector<int>& coef) {
		// TODO: use coef
		LHS.push_back(0);
		LHS_size.push_back(ind.size());
		LHS_used.push_back(0);
		RHS.push_back(rhs);
		int con = RHS.size() - 1;
		for (int i = 0; i < ind.size(); ++i) {
			point_to_constraint.at(ind[i]).push_back(con);
		}
	}

	void set_all_different() { all_different = true; }

	void print() {
		cout << "(";
		for (int i = 0; i < n - 1; ++i) {
			cout << V[i] << ",";
		}
		cout << V[n - 1] << ")\n";
	}

	void search(int k) {
		const auto& constraints = point_to_constraint[k];
		for (int i = 0; i < constraints.size(); ++i) {
			int c = constraints[i];  // index of constraint
			LHS_used[c]++;
		}

		// Go though all values of this variable
		for (int v = minval; v <= maxval; v++) {
			// If this value has not already been used
			if (!all_different || !used[v]) {
				used[v] = true;
				V[k] = v;

				// Go through all constraints with this point
				bool ok = true;
				bool nomore = false;
				for (int i = 0; i < constraints.size(); ++i) {
					int c = constraints[i];  // index of constraint
					LHS[c] += v;

					// Has this constraint overflowed?
					if (LHS[c] > RHS[c]) {
						ok = false;
						nomore = true;
					}

					// Have we assigned all variable to this constraint and it is
					// incorrect?
					if (LHS_used[c] == LHS_size[c] && LHS[c] != RHS[c]) {
						ok = false;
					}
				}

				// Only continue if no constraints have been violated
				if (ok) {
					// Have we reached the end of the tree?
					if (k == n - 1) {
						print();
						count++;
					} else {
						// Search one step further down
						search(k + 1);
					}
				}

				// Restore all constraints with this point
				for (int i = 0; i < constraints.size(); ++i) {
					int c = constraints[i];  // index of constraint
					LHS[c] -= v;
				}

				used[v] = false;

				if (nomore) {
					break;
				}
			}
		}

		for (int i = 0; i < constraints.size(); ++i) {
			int c = constraints[i];  // index of constraint
			LHS_used[c]--;
		}
	}

	long long solve() {
		count = 0;
		search(0);
		return count;
	}

   private:
	int n;
	int minval;
	int maxval;
	long long count;
	vector<vector<int>> point_to_constraint;
	bool all_different;

	vector<bool> used;
	vector<int> V;

	vector<int> RHS;
	vector<int> LHS;
	vector<int> LHS_size;
	vector<int> LHS_used;
};

int main() {
	// 0 1 2
	// 3 4 5
	// 6 7 8
	IntegerFeasibility squareprob(9);
	int rhs = 15;
	squareprob.set_L(1);
	squareprob.set_U(9);

	squareprob.add_equality_constraint(rhs, 0, 1, 2);
	squareprob.add_equality_constraint(rhs, 3, 4, 5);
	squareprob.add_equality_constraint(rhs, 6, 7, 8);

	squareprob.add_equality_constraint(rhs, 0, 3, 6);
	squareprob.add_equality_constraint(rhs, 1, 4, 7);
	squareprob.add_equality_constraint(rhs, 2, 5, 8);

	squareprob.add_equality_constraint(rhs, 0, 4, 8);
	squareprob.add_equality_constraint(rhs, 2, 4, 6);

	squareprob.set_all_different();

	auto start = wall_time();
	long long count = squareprob.solve();
	double time = wall_time() - start;
	cout << "count = " << count << "  Time = " << time << endl;

	//       0   1   2
	//     11  12  13  3
	//   10  17  18  14  4
	//     9  16  15   5
	//       8   7   6
	IntegerFeasibility hexprob2(19);
	rhs = 38;
	hexprob2.add_equality_constraint(rhs, 0, 1, 2);
	hexprob2.add_equality_constraint(rhs, 11, 12, 13, 3);
	hexprob2.add_equality_constraint(rhs, 10, 17, 18, 14, 4);
	hexprob2.add_equality_constraint(rhs, 9, 16, 15, 5);
	hexprob2.add_equality_constraint(rhs, 8, 7, 6);

	hexprob2.add_equality_constraint(rhs, 0, 11, 10);
	hexprob2.add_equality_constraint(rhs, 1, 12, 17, 9);
	hexprob2.add_equality_constraint(rhs, 2, 13, 18, 16, 8);
	hexprob2.add_equality_constraint(rhs, 3, 14, 15, 7);
	hexprob2.add_equality_constraint(rhs, 4, 5, 6);

	hexprob2.add_equality_constraint(rhs, 2, 3, 4);
	hexprob2.add_equality_constraint(rhs, 1, 13, 14, 5);
	hexprob2.add_equality_constraint(rhs, 0, 12, 18, 15, 6);
	hexprob2.add_equality_constraint(rhs, 11, 17, 16, 7);
	hexprob2.add_equality_constraint(rhs, 10, 9, 8);

	hexprob2.set_all_different();

	start = wall_time();
	count = hexprob2.solve();
	time = wall_time() - start;
	cout << "count = " << count << "  Time = " << time << endl;

	//       0   1   2
	//     3   4   5   6
	//   7   8   9  10  11
	//    12  13  14  15
	//      16  17  18
	//   IntegerFeasibility hexprob(19);
	//   rhs = 38;
	//   hexprob.add_equality_constraint(rhs,0, 1, 2);
	//   hexprob.add_equality_constraint(rhs,3, 4, 5, 6);
	//   hexprob.add_equality_constraint(rhs,7, 8, 9, 10, 11);
	//   hexprob.add_equality_constraint(rhs,12, 13, 14, 15);
	//   hexprob.add_equality_constraint(rhs,16, 17, 18);
	//
	//   hexprob.add_equality_constraint(rhs,7, 12, 16);
	//   hexprob.add_equality_constraint(rhs,3, 8, 13, 17);
	//   hexprob.add_equality_constraint(rhs,0, 4, 9, 14, 18);
	//   hexprob.add_equality_constraint(rhs,1, 5, 10, 15);
	//   hexprob.add_equality_constraint(rhs,2, 6, 11);
	//
	//   hexprob.add_equality_constraint(rhs,0, 3, 7);
	//   hexprob.add_equality_constraint(rhs,1, 4, 8, 12);
	//   hexprob.add_equality_constraint(rhs,2, 5, 9, 13, 16);
	//   hexprob.add_equality_constraint(rhs,6, 10, 14, 17);
	//   hexprob.add_equality_constraint(rhs,11, 15, 18);
	//   hexprob.set_all_different();
	//   cout << "count = " << hexprob.solve() << endl;
}
