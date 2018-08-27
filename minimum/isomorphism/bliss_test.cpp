
#include <catch.hpp>

#include <minimum/isomorphism/graph-isomorphism.h>
using namespace minimum::isomorphism;

TEST_CASE("three_nodes") {
	Eigen::SparseMatrix<size_t> graph1(3, 3);
	add_bidirectional_edge(graph1, 0, 1);
	add_bidirectional_edge(graph1, 0, 2);

	bliss_automorphism(graph1);
}
