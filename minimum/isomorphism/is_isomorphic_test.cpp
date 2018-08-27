
#include <random>

#include "catch.hpp"

#include <minimum/isomorphism/graph-isomorphism.h>
using namespace minimum::isomorphism;

// Crashes in GCC with latest Eigen (3.3.1). Works with Eigen 3.2.
// Disable the tests for now as this library is not terribly
// important.
#if !defined(__GNUC__)

TEST_CASE("three_nodes") {
	Eigen::SparseMatrix<size_t> graph1(3, 3);
	add_bidirectional_edge(graph1, 0, 1);
	add_bidirectional_edge(graph1, 0, 2);

	Eigen::SparseMatrix<size_t> graph2(3, 3);
	add_bidirectional_edge(graph2, 1, 2);
	add_bidirectional_edge(graph2, 1, 0);

	std::vector<int> isomorphism;
	REQUIRE(is_isomorphic(graph1, graph2, &isomorphism) > 0);
	CHECK(check_isomorphism(graph1, graph2, isomorphism));

	CHECK(is_isomorphic(graph1, graph2, nullptr) == 2);
}

TEST_CASE("three_nodes_false") {
	Eigen::SparseMatrix<size_t> graph1(3, 3);
	add_bidirectional_edge(graph1, 0, 1);
	add_bidirectional_edge(graph1, 0, 2);

	Eigen::SparseMatrix<size_t> graph2(3, 3);
	add_bidirectional_edge(graph2, 1, 2);

	std::vector<int> isomorphism;
	CHECK(is_isomorphic(graph1, graph2, &isomorphism) == 0);
	CHECK(is_isomorphic(graph1, graph2, nullptr) == 0);

	Eigen::SparseMatrix<size_t> graph3(4, 4);
	add_bidirectional_edge(graph3, 0, 1);
	add_bidirectional_edge(graph3, 0, 2);
	CHECK(is_isomorphic(graph1, graph2, &isomorphism) == 0);
	CHECK(is_isomorphic(graph1, graph2, nullptr) == 0);
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

	std::vector<int> isomorphism;
	REQUIRE(is_isomorphic(graph1, graph2, &isomorphism) > 0);
	CHECK(check_isomorphism(graph1, graph2, isomorphism));

	CHECK(is_isomorphic(graph1, graph2, nullptr) == 8);
}

TEST_CASE("four_nodes_false") {
	Eigen::SparseMatrix<size_t> graph1(4, 4);
	add_bidirectional_edge(graph1, 0, 1);
	add_bidirectional_edge(graph1, 1, 2);
	add_bidirectional_edge(graph1, 2, 3);
	add_bidirectional_edge(graph1, 3, 0);

	Eigen::SparseMatrix<size_t> graph2(4, 4);
	add_bidirectional_edge(graph2, 0, 1);
	add_bidirectional_edge(graph2, 0, 2);
	add_bidirectional_edge(graph2, 0, 3);
	add_bidirectional_edge(graph2, 1, 2);

	std::vector<int> isomorphism;
	CHECK(is_isomorphic(graph1, graph2, &isomorphism) == 0);
	CHECK(is_isomorphic(graph1, graph2, nullptr) == 0);
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

	std::vector<int> isomorphism;
	REQUIRE(is_isomorphic(graph1, graph2, &isomorphism) > 0);
	CHECK(check_isomorphism(graph1, graph2, isomorphism));

	CHECK(is_isomorphic(graph1, graph2, nullptr) == 2 * n);
}

TEST_CASE("ring_plus_extra") {
	int n = 6;
	std::vector<int> indices;
	for (int i = 0; i < n + 1; ++i) {
		indices.emplace_back(i);
	}

	Eigen::SparseMatrix<size_t> graph1(n + 1, n + 1);
	for (int i = 1; i < n; ++i) {
		add_bidirectional_edge(graph1, indices[i], indices[i + 1]);
	}
	add_bidirectional_edge(graph1, indices[n], indices[1]);

	std::mt19937 engine(0u);
	std::shuffle(indices.begin(), indices.end(), engine);

	Eigen::SparseMatrix<size_t> graph2(n + 1, n + 1);
	for (int i = 0; i < n - 1; ++i) {
		add_bidirectional_edge(graph2, indices[i], indices[i + 1]);
	}
	add_bidirectional_edge(graph2, indices[n - 1], indices[0]);

	std::vector<int> isomorphism;
	REQUIRE(is_isomorphic(graph1, graph2, &isomorphism) > 0);
	CHECK(check_isomorphism(graph1, graph2, isomorphism));
}

TEST_CASE("large_ring") {
	int n = 20;
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

	std::vector<int> isomorphism;
	REQUIRE(is_isomorphic(graph1, graph2, &isomorphism) > 0);
	CHECK(check_isomorphism(graph1, graph2, isomorphism));

	CHECK(is_isomorphic(graph1, graph2, nullptr) == 2 * n);
}

// General example where hashing produces false matches due to
// global topology.
TEST_CASE("ring_vs_mobius") {
	int n = 7;
	std::vector<int> indices;
	for (int i = 0; i < 2 * n; ++i) {
		indices.emplace_back(i);
	}

	auto get_ribbon = [&]() {
		Eigen::SparseMatrix<size_t> graph(2 * n, 2 * n);
		// o
		// |
		// o
		add_bidirectional_edge(graph, indices[0], indices[n]);
		for (int i = 0; i < n - 1; ++i) {
			//     o--o
			// ...    |
			//     o--o
			add_bidirectional_edge(graph, indices[i], indices[i + 1]);
			add_bidirectional_edge(graph, indices[n + i], indices[n + i + 1]);
			add_bidirectional_edge(graph, indices[i + 1], indices[n + i + 1]);
		}
		return graph;
	};

	// Make ring strip.
	auto graph1 = get_ribbon();
	add_bidirectional_edge(graph1, indices[0], indices[n - 1]);
	add_bidirectional_edge(graph1, indices[n], indices[2 * n - 1]);

	std::mt19937 engine(0u);
	std::shuffle(indices.begin(), indices.end(), engine);

	// Make Möbius strip.
	auto graph2 = get_ribbon();
	add_bidirectional_edge(graph2, indices[0], indices[2 * n - 1]);
	add_bidirectional_edge(graph2, indices[n], indices[n - 1]);

	CHECK(is_isomorphic(graph1, graph2, nullptr) == 0);
}

// Small graph where the hashing produces false matches.
TEST_CASE("Funkybee-11") {
	using namespace std;
	string graph_data = R"(
		01111110000
		10111101000
		11011010100
		11100001110
		11100001101
		11000010111
		10100101011
		01011010011
		00111100011
		00010111101
		00001111110
		)";
	mt19937 engine(5u);
	auto graph1 = graph_from_string(graph_data, engine);
	auto graph2 = graph_from_string(graph_data, engine);

	std::vector<int> isomorphism;
	// REQUIRE(is_isomorphic(graph1, graph2, &isomorphism) > 0);
	// CHECK(check_isomorphism(graph1, graph2, isomorphism));
	CHECK(is_isomorphic(graph1, graph2, nullptr) > 0);
}

// Strongly regular graph.
// http://www.maths.gla.ac.uk/~es/srgraphs.php
TEST_CASE("SRG-10-3-0-1") {
	using namespace std;
	string graph_data = R"(
		0111000000
		1000110000
		1000001100
		1000000011
		0100001010
		0100000101
		0010100001
		0010010010
		0001100100
		0001011000
		)";
	mt19937 engine(0u);
	auto graph1 = graph_from_string(graph_data, engine);
	auto graph2 = graph_from_string(graph_data, engine);

	std::vector<int> isomorphism;
	REQUIRE(is_isomorphic(graph1, graph2, &isomorphism) > 0);
	CHECK(check_isomorphism(graph1, graph2, isomorphism));
	CHECK(is_isomorphic(graph1, graph2, nullptr) > 0);
}

// Strongly regular graph.
// http://www.maths.gla.ac.uk/~es/srgraphs.php
//
// http://en.wikipedia.org/wiki/File:Paley13.svg
//
TEST_CASE("SRG-13-6-2-3") {
	using namespace std;
	string graph_data = R"(
		0111111000000
		1011000111000
		1100100100110
		1100010010101
		1010001010011
		1001001001110
		1000110101001
		0110001001101
		0101100001011
		0100011110010
		0011010100011
		0010110011100
		0001101110100
		)";
	mt19937 engine(0u);
	auto graph1 = graph_from_string(graph_data, engine);
	auto graph2 = graph_from_string(graph_data, engine);

	std::vector<int> isomorphism;
	REQUIRE(is_isomorphic(graph1, graph2, &isomorphism) > 0);
	CHECK(check_isomorphism(graph1, graph2, isomorphism));
	CHECK(is_isomorphic(graph1, graph2, nullptr) > 0);
}

// Strongly regular graph.
// http://www.maths.gla.ac.uk/~es/srgraphs.php
TEST_CASE("SRG-25-12-5-6") {
	using namespace std;
	string graph_data = R"(
		0111111111111000000000000
		1011111000000111111000000
		1101111000000000000111111
		1110000111000111000111000
		1110000100110100110100110
		1110000010101010101010101
		1110000001011001011001011
		1001100011100101001000111
		1001010101010010110010011
		1001001110001000111101100
		1000110100011110001101001
		1000101010101011010110010
		1000011001110101100011100
		0101100100101001110011001
		0101010010110001101101010
		0101001100011110001010110
		0100110011001110010001110
		0100101011010100101110001
		0100011101100011010100101
		0011100001110010011011100
		0011010010011101010100101
		0011001001101110100100011
		0010110101001001101110010
		0010101110010011100001101
		0010011110100100011011010
		)";
	mt19937 engine(0u);
	auto graph1 = graph_from_string(graph_data, engine);
	auto graph2 = graph_from_string(graph_data, engine);

	std::vector<int> isomorphism;
	// REQUIRE(is_isomorphic(graph1, graph2, &isomorphism) > 0);
	// CHECK(check_isomorphism(graph1, graph2, isomorphism));
	// CHECK(is_isomorphic(graph1, graph2, nullptr) > 0);
}

// Strongly regular graph.
// http://www.maths.gla.ac.uk/~es/srgraphs.php
TEST_CASE("SRG-64-18-2-6") {
	using namespace std;
	string graph_data = R"(
		0111111111111111111000000000000000000000000000000000000000000000
		1011000000000000000111111111111111000000000000000000000000000000
		1101000000000000000000000000000000111111111111111000000000000000
		1110000000000000000000000000000000000000000000000111111111111111
		1000011000000000000111110000000000111110000000000111110000000000
		1000101000000000000000001111100000000001111100000000001111100000
		1000110000000000000000000000011111000000000011111000000000011111
		1000000011000000000111001000010000100001110010000100001000011100
		1000000101000000000000100111001000010000001001110011100100000010
		1000000110000000000000010000100111001110000100001000010011100001
		1000000000011000000111000100000100001001001100001010000100010011
		1000000000101000000000100010011010110100100010000001000011101000
		1000000000110000000000011001100001000010010001110100111000000100
		1000000000000011000100110010001000000011110001000000010010010011
		1000000000000101000010001001110000001000000110101011100001001000
		1000000000000110000001000100000111110100001000010100001100100100
		1000000000000000011100110000100001010001001100010001000001011100
		1000000000000000101010001000011010001110010001000010001100100010
		1000000000000000110001000111000100100000100010101100110010000001
		0100100100100100100000000001000010000100000011011000101111000000
		0100100100100010010000000010000001010000101001100000010001100101
		0100100100100001001000000000101000000010010100110001000010101010
		0100100010010100100000001000000100001000010110100100000100100101
		0100100001001100100000000100010000100000101000101010001000101010
		0100010100001010010000100000000100100100001000011001010010010010
		0100010010100001001000010000010000001100010011000001010001010100
		0100010010010100001010000000000001000110000100011110001000011000
		0100010010001010001100000000000010111011000000000000000000101111
		0100010001001010100001000000001000010100100011000110000100010001
		0100001100010010010000010100000000010011000100010100100110000001
		0100001010010100010001000000100000101001001000001000111001000100
		0100001001100001001000101000000000010011100001000010101001001000
		0100001001010001010100000001000000000000111100100111010000010000
		0100001001001001100010000010000000101001010010000001100110000010
		0010100100010001001000011001001001000000000101000010000101010001
		0010100010010001100010000001110100000000010000001000011010010010
		0010100001100010010000100101001001000000100000010100001010011000
		0010100001010001010100001110100000000001000000100000100000001111
		0010100001001100010001000011010100000000001010000001000101010100
		0010010100100100100000000001011101000100000000100111010000100000
		0010010100010100001010010000100110001000000000010001100100000110
		0010010100001100010001100100000011010000000000001010100001001001
		0010010010100001100010011000001010000010000010000100100010001001
		0010010001100010100001100010010010100000000001000000111000000110
		0010001100010010001100100100100001000010001000000010011000100010
		0010001010001100010110000100100100100000000100000101000010101000
		0010001010001010001011110000000010000101000000000000001111010000
		0010001010001001100101001010010000001000100000000010010001100001
		0010001001100010001100011010001000010000010000000101000100100100
		0001100100001001001000100010110010001001001001001000000001000010
		0001100010100010010000010010100110100001010010010000000010000100
		0001100010010010100001001100000011000011100001001000001000000001
		0001100010001010001100000000011101000100111100000000000000110000
		0001100001001100001010001100001010010001000110010000000100001000
		0001010100001001010100010010001100011000000110100001000000000001
		0001010010100001010100100000110001100010100000101000010000001000
		0001010001010100001101001000010001011000001001100010000000000100
		0001010001010010100110000100001100100010010000110100000000000010
		0001010001010001010011110001000000000001000011011000100000010000
		0001001100100100100000001110100010111010000000100000100000100000
		0001001100010010100001010011000100001100011001000000010100000000
		0001001100001001100010100101001000000110100100001010000010000000
		0001001010100100010001011001000001010100100110000100000001000000
		0001001001100100001010100001110000100100011000010001001000000000
		)";
	mt19937 engine(0u);
	auto graph1 = graph_from_string(graph_data, engine);
	auto graph2 = graph_from_string(graph_data, engine);

	std::vector<int> isomorphism;
	REQUIRE(is_isomorphic(graph1, graph2, &isomorphism) > 0);
	CHECK(check_isomorphism(graph1, graph2, isomorphism));
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

	CHECK(is_isomorphic(graph1, graph2, nullptr) == 72);
}

TEST_CASE("iso_r001_s20") {
	// Auto-generated.

	Eigen::SparseMatrix<size_t> graph1(20, 20);
	add_edge(graph1, 1, 0);
	add_edge(graph1, 3, 0);
	add_edge(graph1, 0, 1);
	add_edge(graph1, 2, 1);
	add_edge(graph1, 4, 1);
	add_edge(graph1, 1, 2);
	add_edge(graph1, 5, 2);
	add_edge(graph1, 15, 2);
	add_edge(graph1, 0, 3);
	add_edge(graph1, 18, 3);
	add_edge(graph1, 1, 4);
	add_edge(graph1, 6, 4);
	add_edge(graph1, 7, 4);
	add_edge(graph1, 11, 4);
	add_edge(graph1, 17, 4);
	add_edge(graph1, 2, 5);
	add_edge(graph1, 8, 5);
	add_edge(graph1, 16, 5);
	add_edge(graph1, 4, 6);
	add_edge(graph1, 4, 7);
	add_edge(graph1, 9, 7);
	add_edge(graph1, 5, 8);
	add_edge(graph1, 10, 8);
	add_edge(graph1, 7, 9);
	add_edge(graph1, 8, 10);
	add_edge(graph1, 12, 10);
	add_edge(graph1, 14, 10);
	add_edge(graph1, 4, 11);
	add_edge(graph1, 10, 12);
	add_edge(graph1, 13, 12);
	add_edge(graph1, 12, 13);
	add_edge(graph1, 10, 14);
	add_edge(graph1, 2, 15);
	add_edge(graph1, 5, 16);
	add_edge(graph1, 19, 16);
	add_edge(graph1, 4, 17);
	add_edge(graph1, 3, 18);
	add_edge(graph1, 16, 19);

	Eigen::SparseMatrix<size_t> graph2(20, 20);
	add_edge(graph2, 4, 0);
	add_edge(graph2, 12, 0);
	add_edge(graph2, 14, 0);
	add_edge(graph2, 7, 1);
	add_edge(graph2, 14, 2);
	add_edge(graph2, 14, 3);
	add_edge(graph2, 0, 4);
	add_edge(graph2, 11, 4);
	add_edge(graph2, 17, 4);
	add_edge(graph2, 11, 5);
	add_edge(graph2, 19, 5);
	add_edge(graph2, 14, 6);
	add_edge(graph2, 1, 7);
	add_edge(graph2, 10, 7);
	add_edge(graph2, 12, 8);
	add_edge(graph2, 16, 8);
	add_edge(graph2, 10, 9);
	add_edge(graph2, 7, 10);
	add_edge(graph2, 9, 10);
	add_edge(graph2, 18, 10);
	add_edge(graph2, 4, 11);
	add_edge(graph2, 5, 11);
	add_edge(graph2, 18, 11);
	add_edge(graph2, 0, 12);
	add_edge(graph2, 8, 12);
	add_edge(graph2, 14, 13);
	add_edge(graph2, 15, 13);
	add_edge(graph2, 0, 14);
	add_edge(graph2, 2, 14);
	add_edge(graph2, 3, 14);
	add_edge(graph2, 6, 14);
	add_edge(graph2, 13, 14);
	add_edge(graph2, 13, 15);
	add_edge(graph2, 8, 16);
	add_edge(graph2, 4, 17);
	add_edge(graph2, 10, 18);
	add_edge(graph2, 11, 18);
	add_edge(graph2, 5, 19);

	std::vector<int> isomorphism;
	REQUIRE(is_isomorphic(graph1, graph2, &isomorphism) > 0);
	CHECK(check_isomorphism(graph1, graph2, isomorphism));

	CHECK(is_isomorphic(graph1, graph2, nullptr) == 6);
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

	CHECK(is_isomorphic(graph1, graph2, nullptr) == 2);
}

TEST_CASE("four_nodes_double_edge_self_loop") {
	Eigen::SparseMatrix<size_t> graph1(4, 4);
	add_bidirectional_edge(graph1, 0, 1);
	add_bidirectional_edge(graph1, 0, 1);
	add_bidirectional_edge(graph1, 1, 2);
	add_bidirectional_edge(graph1, 2, 3);
	add_bidirectional_edge(graph1, 3, 0);
	add_edge(graph1, 0, 0);

	Eigen::SparseMatrix<size_t> graph2(4, 4);
	add_bidirectional_edge(graph2, 0, 2);
	add_bidirectional_edge(graph2, 2, 3);
	add_bidirectional_edge(graph2, 3, 1);
	add_bidirectional_edge(graph2, 3, 1);
	add_bidirectional_edge(graph2, 1, 0);
	add_edge(graph2, 3, 3);

	CHECK(is_isomorphic(graph1, graph2, nullptr) == 1);
}

#endif
