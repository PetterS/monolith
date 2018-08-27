#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#ifdef _X
#undef _X
#endif
#include <Eigen/Sparse>

#include <bliss/graph.hh>
#include <bliss/utils.hh>

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/core/scope_guard.h>
#include <minimum/isomorphism/graph-isomorphism.h>
using minimum::core::make_scope_guard;
using minimum::core::range;

namespace minimum {
namespace isomorphism {

void read_graph(const std::string& filename, Eigen::SparseMatrix<size_t>* graph) {
	using namespace std;
	ifstream fin(filename, ios::binary);
	if (!fin) {
		throw runtime_error("File does not exist.");
	}

	auto read_word = [&fin]() {
		unsigned char c1 = fin.get();
		unsigned char c2 = fin.get();
		if (!fin) {
			throw runtime_error("Could not read file.");
		}
		unsigned short word = c2;
		word <<= 8;
		word += c1;
		// clog << "Read " << word << endl;
		return word;
	};

	unsigned short num_nodes = read_word();

	vector<Eigen::Triplet<size_t>> triplets;
	for (int i = 0; i < num_nodes; ++i) {
		unsigned short num_edges = read_word();
		for (int e = 0; e < num_edges; ++e) {
			unsigned short j = read_word();
			triplets.emplace_back(i, j, 1);
		}
	}

	graph->resize(num_nodes, num_nodes);
	graph->setFromTriplets(triplets.begin(), triplets.end());
}

void graph_to_unit_test(const Eigen::SparseMatrix<size_t>& graph, std::string varname) {
	using namespace std;
	auto n = graph.rows();
	cout << "Eigen::SparseMatrix<size_t> " << varname << "(" << n << ", " << n << ");" << endl;
	for (int k = 0; k < graph.outerSize(); ++k) {
		for (Eigen::SparseMatrix<size_t>::InnerIterator it(graph, k); it; ++it) {
			auto i = it.row();
			auto j = it.col();
			cout << "add_edge(" << varname << ", " << i << ", " << j << ");" << endl;
		}
	}
	cout << endl;
}

void check_sizes(const Eigen::SparseMatrix<size_t>& graph1,
                 const Eigen::SparseMatrix<size_t>& graph2,
                 const std::vector<int>& isomorphism) {
	if ((graph1.rows() != graph1.cols()) || (graph2.rows() != graph2.cols())
	    || (graph1.rows() != graph2.rows()))
		throw std::runtime_error("Graphs do not have the same size.");

	auto n = graph1.rows();
	if (isomorphism.size() != n) {
		throw std::runtime_error("Isomorphism vector has incorrect size.");
	}
}

void check_graphs_and_isomorphism(const Eigen::SparseMatrix<size_t>& graph1,
                                  const Eigen::SparseMatrix<size_t>& graph2,
                                  const std::vector<int>& isomorphism) {
	using namespace std;

	check_sizes(graph1, graph2, isomorphism);
	auto n = graph1.rows();
	// Count the referenced elements to determine whether
	// we have a bijection.
	vector<int> element_count(n, 0);
	for (int e : isomorphism) {
		if (e < 0 || e >= n) {
			throw std::logic_error("Isomorphism vector has an incorrect element.");
		}
		element_count[e]++;
	}
	for (int e : element_count) {
		if (e != 1) {
			throw std::runtime_error("Invalid isomorphism.");
		}
	}
}

bool check_isomorphism(const Eigen::SparseMatrix<size_t>& graph1,
                       const Eigen::SparseMatrix<size_t>& graph2,
                       const std::vector<int>& isomorphism) {
	if (graph1.rows() != graph1.cols()) return false;
	if (graph2.rows() != graph2.cols()) return false;
	if (graph1.rows() != graph2.rows()) return false;
	if (graph1.nonZeros() != graph2.nonZeros()) return false;
	check_graphs_and_isomorphism(graph1, graph2, isomorphism);

	for (int k = 0; k < graph1.outerSize(); ++k) {
		for (Eigen::SparseMatrix<size_t>::InnerIterator it(graph1, k); it; ++it) {
			auto i = it.row();
			auto j = it.col();
			auto val = it.value();

			auto i2 = isomorphism[i];
			auto j2 = isomorphism[j];
			auto val2 = graph2.coeff(i2, j2);

			if (val != val2) {
				return false;
			}
		}
	}

	return true;
}

inline void hash_combine(std::size_t& seed, size_t v) {
	// http://stackoverflow.com/questions/2590677/how-do-i-combine-hash-values-in-c0x
	seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// returns possible, where possible[i] is a list of values that
// can be assigned to isomorphism[i].
std::vector<std::vector<int>> get_possible_assignments(const Eigen::SparseMatrix<size_t>& graph1,
                                                       const Eigen::SparseMatrix<size_t>& graph2,
                                                       const std::vector<int>& isomorphism) {
	using namespace std;

	check_sizes(graph1, graph2, isomorphism);
	auto n = graph1.rows();

	// Find all indices that have not been used yet.
	vector<int> non_assigned_indices;
	for (int i = 0; i < n; ++i) {
		non_assigned_indices.emplace_back(i);
	}
	for (int e : isomorphism) {
		if (e >= 0) {
			auto itr = find(non_assigned_indices.begin(), non_assigned_indices.end(), e);
			non_assigned_indices.erase(itr);
		}
	}

	// Set up neighbors vector and add preassignments.
	typedef Eigen::Matrix<size_t, Eigen::Dynamic, 1> Vector;
	Vector ones(n);
	ones.setConstant(1);
	Vector neighbors1 = ones;
	Vector neighbors2 = ones;
	hash<int> hasher;
	for (int i = 0; i < n; ++i) {
		// Has this node already been assigned?
		if (isomorphism[i] >= 0) {
			auto h = hasher(i);
			neighbors1[i] = h;
			neighbors2[isomorphism[i]] = h;
		}
	}

	vector<vector<int>> possible(n);
	vector<size_t> hash_copy_1;
	vector<size_t> hash_copy_2;

	vector<size_t> hashes1(n, 0);
	vector<size_t> hashes2(n, 0);
	for (int iter = 1; iter <= n; ++iter) {
		neighbors1 = graph1 * neighbors1;
		neighbors2 = graph2 * neighbors2;

		for (int i = 0; i < n; ++i) {
			hash_combine(hashes1[i], neighbors1[i]);
			hash_combine(hashes2[i], neighbors2[i]);
		}

		for (auto& p : possible) {
			p.clear();
		}

		// The sets of hashes have to be identical for an
		// isomorphism to be possible.
		hash_copy_1 = hashes1;
		hash_copy_2 = hashes2;
		sort(hash_copy_1.begin(), hash_copy_1.end());
		sort(hash_copy_2.begin(), hash_copy_2.end());
		if (hash_copy_1 != hash_copy_2) {
			// Return empty possibilities. Will result in early exit.
			return possible;
		}

		hash_copy_2.clear();
		for (int j : non_assigned_indices) {
			hash_copy_2.push_back(hashes2[j]);
		}

		bool no_choices = true;
		for (int i = 0; i < n; ++i) {
			// Do we need to assign this isomorphism index?
			if (isomorphism[i] < 0) {
				// See what indices we can assign here.
				for (size_t ind = 0; ind < hash_copy_2.size(); ++ind) {
					if (hashes1[i] == hash_copy_2[ind]) {
						possible[i].emplace_back(non_assigned_indices[ind]);
					}
				}
				if (possible[i].size() > 1) {
					no_choices = false;
				}
			}
		}

		if (no_choices) {
			// We don't need to continue and compute a more
			// advanced hash. We already have a uniquely determined
			// assignment. (early exit)
			return possible;
		}
	}

	return possible;
}

int internal_is_isomorphic(const Eigen::SparseMatrix<size_t>& graph1,
                           const Eigen::SparseMatrix<size_t>& graph2,
                           std::vector<int>* isomorphism,
                           std::vector<std::vector<int>>* all_isomorphisms,
                           bool count_all_isomorphisms,
                           int indentation) {
	using namespace std;

	const int verbosity = 0;
	auto n = graph1.rows();

	if (verbosity >= 3) {
		cout << setw(indentation) << setfill('-') << "";
		cout << setw(n - indentation) << setfill(' ') << "";
		cout << "[";
		for (int e : *isomorphism) {
			cout << setw(2) << setfill(' ') << right;
			if (e >= 0)
				cout << e;
			else
				cout << '*';
			cout << ' ';
		}
		cout << "]" << endl;
	}

	auto possible = get_possible_assignments(graph1, graph2, *isomorphism);

	if (verbosity >= 3) {
		cout << setw(indentation) << setfill('-') << "";
		cout << setw(n - indentation) << setfill(' ') << "";
		cout << " Possibilities: [";
		for (int i = 0; i < n; ++i) {
			cout << setw(2) << setfill(' ') << right;
			if ((*isomorphism)[i] < 0)
				cout << possible[i].size();
			else
				cout << '-';
			cout << ' ';
		}
		cout << "]" << endl;
	}

	for (int i = 0; i < n; ++i) {
		if ((*isomorphism)[i] < 0 && possible[i].empty()) {
			// Assignment is needed for i, but not possible (early exit).
			return 0;
		}
	}

	if (indentation == 0 && verbosity >= 1) {
		cout << setw(indentation) << setfill('=') << "" << setfill(' ');
		cout << "-- Root node stats: ";
		int indices_to_choose = 0;
		double sum = 0;
		for (const auto& p : possible) {
			if (p.size() > 1) {
				indices_to_choose++;
			}
			sum += p.size();
		}
		cout << indices_to_choose << " indetermined.";
		if (indices_to_choose > 0) {
			cout << " avg. poss. = " << setprecision(2) << sum / double(indices_to_choose) << "."
			     << endl;
		}
	}

	// Make all assignments that only have one option.
	for (size_t i = 0; i < isomorphism->size(); ++i) {
		if (possible[i].size() == 1) {
			(*isomorphism)[i] = possible[i].front();
		}
	}
	// Undo this when we leave the function.
	auto restore_isomorphism = make_scope_guard([&]() {
		for (size_t i = 0; i < isomorphism->size(); ++i) {
			if (possible[i].size() == 1) {
				(*isomorphism)[i] = -1;
			}
		}
	});

	// Check that the isomorphism is still sane.
	vector<int> element_count(n, 0);
	for (int e : *isomorphism) {
		if (e >= 0) {
			element_count[e]++;
		}
	}
	for (int e : element_count) {
		if (e > 1) {
			// Invalid isomorphism after fixing mandatory assignments (early exit).
			return 0;
		}
	}

	// Find a position in the vector that still needs to
	// be assigned.
	int assignable_index = -1;
	for (size_t i = 0; i < isomorphism->size(); ++i) {
		if (possible[i].size() > 1) {
			assignable_index = i;
			break;
		}
	}

	// Do we have nothing left to assign?
	if (assignable_index < 0) {
		// Is this an invalid isomorphism? (early exit)
		bool incomplete_isomorphism =
		    find_if(isomorphism->begin(), isomorphism->end(), [](int i) { return i < 0; })
		    != isomorphism->end();
		if (incomplete_isomorphism) {
			return 0;
		}

		if (verbosity >= 2) {
			cout << setw(indentation) << setfill('-') << "";
			cout << setw(n - indentation) << setfill(' ') << "";
			cout << "Checking: ";
		}
		bool isomorphic = check_isomorphism(graph1, graph2, *isomorphism);
		if (verbosity >= 2) {
			cout << (isomorphic ? "yes" : "no") << endl;
		}
		if (!count_all_isomorphisms && isomorphic) {
			// Keep the modifications to the isomorphism, because
			// we want to return the isomorphism to the user.
			restore_isomorphism.dismiss();
		}
		if (isomorphic && all_isomorphisms != nullptr) {
			all_isomorphisms->emplace_back(*isomorphism);
		}
		return isomorphic ? 1 : 0;
	}

	int num_isomorphisms = 0;
	for (int e : possible[assignable_index]) {
		(*isomorphism)[assignable_index] = e;
		int this_num_isomorphisms = internal_is_isomorphic(
		    graph1, graph2, isomorphism, all_isomorphisms, count_all_isomorphisms, indentation + 1);
		if (indentation == 0 && verbosity >= 1) {
			cout << "-- p(" << assignable_index << ") = " << e << " gave " << this_num_isomorphisms
			     << " isomorphisms." << endl;
		}
		if (!count_all_isomorphisms && this_num_isomorphisms > 0) {
			// Keep the modifications to the isomorphism, because
			// we want to return the isomorphism to the user.
			restore_isomorphism.dismiss();
			return 1;
		}
		num_isomorphisms += this_num_isomorphisms;
		(*isomorphism)[assignable_index] = -1;
	}

	return num_isomorphisms;
}

bool is_possibly_isomorphic(const Eigen::SparseMatrix<size_t>& graph1,
                            const Eigen::SparseMatrix<size_t>& graph2) {
	// size_t matrix multiplication does not compile.
	Eigen::SparseMatrix<ptrdiff_t> G1 = graph1.cast<ptrdiff_t>();
	Eigen::SparseMatrix<ptrdiff_t> G2 = graph2.cast<ptrdiff_t>();
	for (int iter = 1; iter <= graph1.rows(); ++iter) {
		if (G1.diagonal().sum() != G2.diagonal().sum()) {
			return false;
		}
		if (G1.nonZeros() != G2.nonZeros()) {
			return false;
		}
		if (G1.sum() != G2.sum()) {
			return false;
		}
		if (G1.nonZeros() == 0) {
			break;
		}
		G1 = G1 * graph1.cast<ptrdiff_t>();
		G2 = G2 * graph2.cast<ptrdiff_t>();
	}
	return true;

	// Use this code when Eigen >= 3.3 works
	//
	// Eigen::SparseMatrix<size_t> G1 = graph1;
	// Eigen::SparseMatrix<size_t> G2 = graph2;
	// for (int iter = 1; iter <= graph1.rows(); ++iter) {
	//	if (G1.nonZeros() != G2.nonZeros()) {
	//		return false;
	//	}
	//	if (G1.nonZeros() == 0) {
	//		break;
	//	}
	//	if (G1.sum() != G2.sum()) {
	//		return false;
	//	}

	//	// TOOD: Compare trace.
	//	//if (G1.diagonal().sum() != G2.diagonal().sum()) {
	//	//	return false;
	//	//}
	//
	//	G1 = G1 * graph1;
	//	G2 = G2 * graph2;
	//}
	// return true;
}

int is_isomorphic(const Eigen::SparseMatrix<size_t>& graph1,
                  const Eigen::SparseMatrix<size_t>& graph2,
                  std::vector<int>* isomorphism) {
	using namespace std;

	if (!is_possibly_isomorphic(graph1, graph2)) {
		// cout << "-- Can not possibly be isomorphic." << endl;
		return 0;
	}

	bool count_all_isomorphisms = isomorphism == nullptr;
	vector<int> local_isomorphism;
	if (count_all_isomorphisms) {
		isomorphism = &local_isomorphism;
	}

	auto n = graph1.rows();
	isomorphism->clear();
	isomorphism->resize(n, -1);

	return internal_is_isomorphic(graph1, graph2, isomorphism, nullptr, count_all_isomorphisms, 0);
}

std::vector<std::vector<int>> all_isomorphisms(const Eigen::SparseMatrix<size_t>& graph1,
                                               const Eigen::SparseMatrix<size_t>& graph2) {
	using namespace std;

	if (!is_possibly_isomorphic(graph1, graph2)) {
		return {};
	}

	auto n = graph1.rows();
	vector<int> isomorphism(n, -1);
	std::vector<std::vector<int>> all;

	int num_isomorphisms = internal_is_isomorphic(graph1, graph2, &isomorphism, &all, true, 0);
	minimum_core_assert(num_isomorphisms == all.size());

	return all;
}

std::vector<Permutation> bliss_automorphism(const Eigen::SparseMatrix<size_t>& graph) {
	bliss::Digraph bliss_graph(graph.rows());

	for (int k = 0; k < graph.outerSize(); ++k) {
		for (Eigen::SparseMatrix<size_t>::InnerIterator it(graph, k); it; ++it) {
			auto i = it.row();
			auto j = it.col();
			auto val = it.value();

			if (i == j) {
				// Loops are converted to colors.
				bliss_graph.change_color(i, val);
			} else {
				auto k = bliss_graph.add_vertex(val);
				bliss_graph.add_edge(i, k);
				bliss_graph.add_edge(k, j);
			}
		}
	}

	struct CallbackData {
		int n;
		std::vector<Permutation> generators;
	};

	CallbackData callback_data;
	callback_data.n = graph.rows();

	auto callback = [](void* data, unsigned int n, const unsigned int* aut) -> void {
		CallbackData* callback_data = (CallbackData*)data;
		std::vector<int> generator(n);
		std::copy(aut, aut + n, generator.begin());
		callback_data->generators.emplace_back(std::move(generator));

		// bliss::print_permutation(stderr, n, aut);
		// std::cerr << "\n";
	};

	bliss::Stats stats;
	bliss_graph.find_automorphisms(stats, callback, &callback_data);
	std::cerr << "-- Bliss: group size = " << stats.get_group_size_approx() << "\n";

	return std::move(callback_data.generators);
}
}  // namespace isomorphism
}  // namespace minimum
