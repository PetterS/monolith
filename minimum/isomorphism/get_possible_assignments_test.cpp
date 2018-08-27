
#include <random>

#include "catch.hpp"

#include <minimum/isomorphism/graph-isomorphism.h>
using namespace minimum::isomorphism;

TEST_CASE("three_nodes") {
	Eigen::SparseMatrix<size_t> graph1(3, 3);
	add_bidirectional_edge(graph1, 0, 1);
	add_bidirectional_edge(graph1, 0, 2);

	Eigen::SparseMatrix<size_t> graph2(3, 3);
	add_bidirectional_edge(graph2, 1, 2);
	add_bidirectional_edge(graph2, 1, 0);

	std::vector<int> isomorphism;
	CHECK_THROWS(get_possible_assignments(graph1, graph2, isomorphism));
	isomorphism.resize(3, -1);
	auto possible = get_possible_assignments(graph1, graph2, isomorphism);
	REQUIRE(possible.size() == 3);
	CHECK(possible[0].size() == 1);
	CHECK(possible[1].size() == 2);
	CHECK(possible[2].size() == 2);
}

TEST_CASE("ring") {
	int n = 6;
	std::vector<int> indices;
	for (int i = 0; i < n; ++i) {
		indices.emplace_back(i);
	}

	Eigen::SparseMatrix<size_t> graph1(n, n);
	for (int i = 0; i < n - 1; ++i) {
		add_bidirectional_edge(graph1, indices[i], indices[i + 1]);
	}
	add_bidirectional_edge(graph1, indices[n - 1], indices[0]);

	std::mt19937 engine(0u);
	std::shuffle(indices.begin(), indices.end(), engine);

	Eigen::SparseMatrix<size_t> graph2(n, n);
	for (int i = 0; i < n - 1; ++i) {
		add_bidirectional_edge(graph2, indices[i], indices[i + 1]);
	}
	add_bidirectional_edge(graph2, indices[n - 1], indices[0]);

	std::vector<int> isomorphism(n, -1);
	auto possible = get_possible_assignments(graph1, graph2, isomorphism);
	REQUIRE(possible.size() == n);
	for (int i = 0; i < n; ++i) {
		CHECK(possible[i].size() == n);
	}
}

TEST_CASE("single_possibility") {
	// This test checks that the hashing function can establish long-range
	// relationships between nodes. Locally, several nodes may seem similar,
	// but globally there is only one option.

	int n = 100;
	std::vector<int> indices;
	for (int i = 0; i < n; ++i) {
		indices.emplace_back(i);
	}

	auto create_graph = [&]() {
		Eigen::SparseMatrix<size_t> graph(n, n);
		// Create ring.
		for (int i = 0; i < n - 4; ++i) {
			add_bidirectional_edge(graph, indices[i], indices[i + 1]);
		}
		add_bidirectional_edge(graph, indices[n - 4], indices[0]);
		// Add three arms to the ring.
		add_bidirectional_edge(graph, indices[n - 3], indices[0]);
		add_bidirectional_edge(graph, indices[n - 2], indices[3]);
		add_bidirectional_edge(graph, indices[n - 1], indices[18]);
		return graph;
	};

	auto graph1 = create_graph();
	std::mt19937 engine(0u);
	std::shuffle(indices.begin(), indices.end(), engine);
	auto graph2 = create_graph();

	std::vector<int> isomorphism(n, -1);
	auto possible = get_possible_assignments(graph1, graph2, isomorphism);
	REQUIRE(possible.size() == n);
	for (int i = 0; i < n; ++i) {
		CHECK(possible[i].size() == 1);
	}
}
