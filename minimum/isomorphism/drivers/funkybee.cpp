
// Tough graphs from
// http://funkybee.narod.ru/graphs.htm

#include <iostream>
#include <random>

#include <minimum/core/time.h>

#include <minimum/isomorphism/graph-isomorphism.h>
using namespace minimum::isomorphism;

std::vector<bool> ground_truth = {
    true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, true,  false, true,  true,  false, false, false};

int main_function() {
	using namespace std;

	ifstream funkybee("funkybee.txt");

	int iter = 0;
	while (true) {
		int n = -1;
		funkybee >> n;
		if (!funkybee || n <= 0) {
			break;
		}

		Eigen::SparseMatrix<size_t> graph1(n, n);
		for (int i = 0; i < n; ++i) {
			string row;
			funkybee >> row;
			if (row.size() != n) {
				throw runtime_error("Invalid file.");
			}
			for (int j = 0; j < n; ++j) {
				if (row[j] == '1') {
					add_edge(graph1, i, j);
				}
			}
		}

		Eigen::SparseMatrix<size_t> graph2(n, n);
		for (int i = 0; i < n; ++i) {
			string row;
			funkybee >> row;
			if (row.size() != n) {
				throw runtime_error("Invalid file.");
			}
			for (int j = 0; j < n; ++j) {
				if (row[j] == '1') {
					add_edge(graph2, i, j);
				}
			}
		}

		double start_time = minimum::core::wall_time();
		vector<int> isomorphism;
		int num_isomorphisms = is_isomorphic(graph1, graph2, &isomorphism);
		double elapsed_time = minimum::core::wall_time() - start_time;
		cout << setw(3) << right << iter + 1 << ") n = " << setw(3) << right << n << "  ";
		if (num_isomorphisms == 0) {
			cout << "NON-isomorphic ";
			if (ground_truth.at(iter)) {
				throw logic_error("Invalid result.");
			}
		} else if (!ground_truth.at(iter)) {
			throw logic_error("Invalid result.");
		}

		cout << setprecision(3) << elapsed_time << endl;

		++iter;
	}

	return 0;
}

int main() {
	using namespace std;
	try {
		return main_function();
	} catch (exception& err) {
		cerr << endl << "ERROR: " << err.what() << endl;
	}
	return 1;
}
