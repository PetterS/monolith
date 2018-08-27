
#include <random>
using namespace std;

#include <catch.hpp>

#include <minimum/isomorphism/graph-isomorphism.h>
using namespace minimum::isomorphism;

// Crashes in GCC with latest Eigen (3.3.1). Works with Eigen 3.2.
// Disable the tests for now as this library is not terribly
// important.
#if !defined(__GNUC__)

void verify(const Eigen::SparseMatrix<size_t>& graph1,
            const Eigen::SparseMatrix<size_t>& graph2,
            int expected_count) {
	auto isomorphisms = all_isomorphisms(graph1, graph2);
	CHECK(isomorphisms.size() == expected_count);
	for (auto& isomorphism : isomorphisms) {
		CHECK(check_isomorphism(graph1, graph2, isomorphism));
	}
}

TEST_CASE("three_nodes") {
	Eigen::SparseMatrix<size_t> graph1(3, 3);
	add_bidirectional_edge(graph1, 0, 1);
	add_bidirectional_edge(graph1, 0, 2);

	Eigen::SparseMatrix<size_t> graph2(3, 3);
	add_bidirectional_edge(graph2, 1, 2);
	add_bidirectional_edge(graph2, 1, 0);

	verify(graph1, graph2, 2);
}

TEST_CASE("fixing_mandatory") {
	Eigen::SparseMatrix<size_t> graph1(16, 16);
	add_edge(graph1, 2, 0);
	add_edge(graph1, 7, 0);
	add_edge(graph1, 6, 1);
	add_edge(graph1, 13, 2);
	add_edge(graph1, 1, 3);
	add_edge(graph1, 12, 3);
	add_edge(graph1, 6, 4);
	add_edge(graph1, 2, 5);
	add_edge(graph1, 14, 5);
	add_edge(graph1, 13, 7);
	add_edge(graph1, 4, 8);
	add_edge(graph1, 12, 8);
	add_edge(graph1, 7, 9);
	add_edge(graph1, 14, 9);
	add_edge(graph1, 3, 10);
	add_edge(graph1, 8, 10);
	add_edge(graph1, 15, 10);
	add_edge(graph1, 0, 11);
	add_edge(graph1, 5, 11);
	add_edge(graph1, 9, 11);
	add_edge(graph1, 6, 12);
	add_edge(graph1, 13, 14);
	add_edge(graph1, 1, 15);
	add_edge(graph1, 4, 15);

	Eigen::SparseMatrix<size_t> graph2(16, 16);
	add_edge(graph2, 2, 0);
	add_edge(graph2, 3, 0);
	add_edge(graph2, 5, 2);
	add_edge(graph2, 5, 3);
	add_edge(graph2, 5, 4);
	add_edge(graph2, 1, 6);
	add_edge(graph2, 2, 7);
	add_edge(graph2, 4, 7);
	add_edge(graph2, 3, 8);
	add_edge(graph2, 4, 8);
	add_edge(graph2, 0, 9);
	add_edge(graph2, 7, 9);
	add_edge(graph2, 8, 9);
	add_edge(graph2, 14, 10);
	add_edge(graph2, 15, 10);
	add_edge(graph2, 6, 11);
	add_edge(graph2, 15, 11);
	add_edge(graph2, 6, 12);
	add_edge(graph2, 14, 12);
	add_edge(graph2, 10, 13);
	add_edge(graph2, 11, 13);
	add_edge(graph2, 12, 13);
	add_edge(graph2, 1, 14);
	add_edge(graph2, 1, 15);

	verify(graph1, graph2, 72);
}

#endif
