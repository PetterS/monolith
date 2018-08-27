
#include <catch.hpp>

#include <minimum/isomorphism/graph-isomorphism.h>
using namespace minimum::isomorphism;

TEST_CASE("three_nodes") {
	Eigen::SparseMatrix<size_t> graph1(3, 3);
	add_bidirectional_edge(graph1, 0, 1);
	add_bidirectional_edge(graph1, 0, 2);

	Eigen::SparseMatrix<size_t> graph2(3, 3);
	add_bidirectional_edge(graph2, 1, 2);
	add_bidirectional_edge(graph2, 1, 0);

	CHECK(check_isomorphism(graph1, graph2, {1, 2, 0}));
	CHECK(check_isomorphism(graph1, graph2, {1, 0, 2}));

	CHECK_FALSE(check_isomorphism(graph1, graph2, {0, 1, 2}));
}

TEST_CASE("three_nodes_ii") {
	Eigen::SparseMatrix<size_t> graph1(3, 3);
	add_bidirectional_edge(graph1, 0, 1);
	add_bidirectional_edge(graph1, 0, 2);

	Eigen::SparseMatrix<size_t> graph2(3, 3);
	add_bidirectional_edge(graph2, 0, 1);
	add_bidirectional_edge(graph2, 0, 2);
	add_bidirectional_edge(graph2, 1, 2);

	CHECK_FALSE(check_isomorphism(graph1, graph2, {0, 1, 2}));
}

TEST_CASE("three_nodes_false") {
	Eigen::SparseMatrix<size_t> graph1(3, 3);
	add_bidirectional_edge(graph1, 0, 1);
	add_bidirectional_edge(graph1, 0, 2);

	Eigen::SparseMatrix<size_t> graph2(3, 3);
	add_bidirectional_edge(graph2, 0, 1);
	add_bidirectional_edge(graph2, 0, 2);
	add_bidirectional_edge(graph2, 1, 2);

	CHECK_FALSE(check_isomorphism(graph1, graph2, {0, 1, 2}));
	CHECK_FALSE(check_isomorphism(graph2, graph1, {0, 1, 2}));
}

TEST_CASE("four_nodes") {
	Eigen::SparseMatrix<size_t> graph1(4, 4);
	add_bidirectional_edge(graph1, 0, 1);
	add_bidirectional_edge(graph1, 1, 2);
	add_bidirectional_edge(graph1, 2, 3);
	add_bidirectional_edge(graph1, 3, 0);

	Eigen::SparseMatrix<size_t> graph2(4, 4);
	add_bidirectional_edge(graph2, 0, 2);
	add_bidirectional_edge(graph2, 2, 3);
	add_bidirectional_edge(graph2, 3, 1);
	add_bidirectional_edge(graph2, 1, 0);

	CHECK(check_isomorphism(graph1, graph2, {0, 2, 3, 1}));
	CHECK(check_isomorphism(graph1, graph2, {1, 3, 2, 0}));

	CHECK_FALSE(check_isomorphism(graph1, graph2, {0, 1, 2, 3}));
	CHECK_FALSE(check_isomorphism(graph1, graph2, {3, 2, 1, 0}));
}

TEST_CASE("four_nodes_double_edge") {
	Eigen::SparseMatrix<size_t> graph1(4, 4);
	add_bidirectional_edge(graph1, 0, 1);
	add_bidirectional_edge(graph1, 0, 1);
	add_bidirectional_edge(graph1, 1, 2);
	add_bidirectional_edge(graph1, 2, 3);
	add_bidirectional_edge(graph1, 3, 0);

	Eigen::SparseMatrix<size_t> graph2(4, 4);
	add_bidirectional_edge(graph2, 0, 2);
	add_bidirectional_edge(graph2, 2, 3);
	add_bidirectional_edge(graph2, 3, 1);
	add_bidirectional_edge(graph2, 3, 1);
	add_bidirectional_edge(graph2, 1, 0);

	CHECK(check_isomorphism(graph1, graph2, {3, 1, 0, 2}));
	CHECK_FALSE(check_isomorphism(graph1, graph2, {0, 1, 2, 3}));
}
