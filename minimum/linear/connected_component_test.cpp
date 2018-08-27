#include <iostream>

#include <catch.hpp>

#include <minimum/linear/connected_component.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/solver.h>
using namespace minimum::linear;

void line_with_edges_always_active(IP* ip, ConnectedComponent<>* connected) {
	const int n = 6;
	auto x = ip->add_boolean_vector(n);

	for (int i = 0; i < n; ++i) {
		connected->add_node(x[i]);
		if (i > 0) {
			connected->add_edge(x[i - 1], x[i]);
		}
	}
	connected->set_first_node(x[1]);
	connected->finish();

	ip->add_objective(100 * x[0] - 10 * x[1] + x[2] + 2 * x[3] - 10 * x[4] + 50 * x[5]);
	REQUIRE(solve(ip));
	connected->debug();

	CHECK_FALSE(x[0].bool_value());
	CHECK(x[1].bool_value());
	CHECK(x[2].bool_value());
	CHECK(x[3].bool_value());
	CHECK(x[4].bool_value());
	CHECK_FALSE(x[5].bool_value());
}

TEST_CASE("flow_line_with_edges_always_active") {
	IP ip;
	ConnectedComponentFlow<> connected(&ip);
	line_with_edges_always_active(&ip, &connected);
}

TEST_CASE("noflow_line_with_edges_always_active") {
	IP ip;
	ConnectedComponentNoFlow<> connected(&ip, 6);
	line_with_edges_always_active(&ip, &connected);
}

void connected_with_edge_variables(IP* ip, ConnectedComponent<>* connected) {
	//
	//       0
	//   0------1
	//   |  4 / |
	// 2 |  /   | 1
	//   |/     |
	//   2------3
	//       3
	//
	auto e = ip->add_boolean_vector(5);
	connected->add_node(0);
	connected->add_node(1);
	connected->add_node(2);
	connected->add_node(3);
	connected->add_edge(0, 1, e[0]);
	connected->add_edge(1, 3, e[1]);
	connected->add_edge(0, 2, e[2]);
	connected->add_edge(2, 3, e[3]);
	connected->add_edge(1, 2, e[4]);

	connected->finish();
	ip->add_objective(-e[0] - e[1] + e[4] + 10 * e[2] + 10 * e[3]);
	REQUIRE(solve(ip));
	connected->debug();

	CHECK(e[0].bool_value());
	CHECK(e[1].bool_value());
	CHECK_FALSE(e[2].bool_value());
	CHECK_FALSE(e[3].bool_value());
	CHECK(e[4].bool_value());
}

TEST_CASE("flow_connected_with_edge_variables") {
	IP ip;
	ConnectedComponentFlow<> connected(&ip);
	connected_with_edge_variables(&ip, &connected);
}

TEST_CASE("noflow_connected_with_edge_variables") {
	IP ip;
	ConnectedComponentNoFlow<> connected(&ip, 4);
	connected_with_edge_variables(&ip, &connected);
}
