#pragma once

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <set>
#include <string>
#include <vector>

#ifdef _X
#undef _X
#endif
#include <Eigen/Sparse>

#include <minimum/isomorphism/export.h>

namespace minimum {
namespace isomorphism {

MINIMUM_ISOMORPHISM_API void read_graph(const std::string& filename,
                                        Eigen::SparseMatrix<size_t>* graph);

namespace {
void add_edge(Eigen::SparseMatrix<size_t>& matrix,
              Eigen::SparseMatrix<size_t>::Index i,
              Eigen::SparseMatrix<size_t>::Index j) {
	matrix.coeffRef(i, j) += 1;
}

void add_bidirectional_edge(Eigen::SparseMatrix<size_t>& matrix,
                            Eigen::SparseMatrix<size_t>::Index i,
                            Eigen::SparseMatrix<size_t>::Index j) {
	matrix.coeffRef(i, j) += 1;
	matrix.coeffRef(j, i) += 1;
}
}  // namespace

template <typename RandomEngine>
Eigen::SparseMatrix<size_t> graph_from_string(const std::string& graph_string,
                                              RandomEngine& engine) {
	using namespace std;
	size_t num_nodes = 0;
	{
		stringstream sin(graph_string);
		string row;
		sin >> row;
		num_nodes = row.size();
	}
	if (num_nodes <= 0) {
		return {};
	}

	stringstream sin(graph_string);
	Eigen::SparseMatrix<size_t> graph(num_nodes, num_nodes);

	vector<int> translation;
	for (int i = 0; i < num_nodes; ++i) {
		translation.push_back(i);
	}
	shuffle(translation.begin(), translation.end(), engine);

	for (int i = 0; i < num_nodes; ++i) {
		string row;
		sin >> row;
		if (row.size() != num_nodes) {
			throw runtime_error("Invalid string.");
		}
		for (int j = 0; j < num_nodes; ++j) {
			if (row[j] == '1') {
				add_edge(graph, translation[i], translation[j]);
			}
		}
	}
	return graph;
}

MINIMUM_ISOMORPHISM_API void graph_to_unit_test(const Eigen::SparseMatrix<size_t>& graph,
                                                std::string varname);

MINIMUM_ISOMORPHISM_API void check_sizes(const Eigen::SparseMatrix<size_t>& graph1,
                                         const Eigen::SparseMatrix<size_t>& graph2,
                                         const std::vector<int>& isomorphism);

MINIMUM_ISOMORPHISM_API void check_graphs_and_isomorphism(const Eigen::SparseMatrix<size_t>& graph1,
                                                          const Eigen::SparseMatrix<size_t>& graph2,
                                                          const std::vector<int>& isomorphism);

MINIMUM_ISOMORPHISM_API bool check_isomorphism(const Eigen::SparseMatrix<size_t>& graph1,
                                               const Eigen::SparseMatrix<size_t>& graph2,
                                               const std::vector<int>& isomorphism);

// returns possible, where possible[i] is a list of values that
// can be assigned to isomorphism[i].
MINIMUM_ISOMORPHISM_API std::vector<std::vector<int>> get_possible_assignments(
    const Eigen::SparseMatrix<size_t>& graph1,
    const Eigen::SparseMatrix<size_t>& graph2,
    const std::vector<int>& isomorphism);

MINIMUM_ISOMORPHISM_API int is_isomorphic(const Eigen::SparseMatrix<size_t>& graph1,
                                          const Eigen::SparseMatrix<size_t>& graph2,
                                          std::vector<int>* isomorphism);

MINIMUM_ISOMORPHISM_API std::vector<std::vector<int>> all_isomorphisms(
    const Eigen::SparseMatrix<size_t>& graph1, const Eigen::SparseMatrix<size_t>& graph2);

class MINIMUM_ISOMORPHISM_API Permutation {
   public:
	Permutation(std::vector<int> perm_) : perm(std::move(perm_)) {}
	int permute(int i) const { return perm.at(i); }

   private:
	std::vector<int> perm;
};

MINIMUM_ISOMORPHISM_API std::vector<Permutation> bliss_automorphism(
    const Eigen::SparseMatrix<size_t>& graph);
}  // namespace isomorphism
}  // namespace minimum
