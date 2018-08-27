#pragma once

#include <map>
#include <vector>

#include <minimum/core/check.h>
#include <minimum/linear/ip.h>

namespace minimum {
namespace linear {

// Enforces that a network of nodes should be connected.
//
// A node can be added in two ways:
//  1) A node that is always part of the connected component.
//  2) A node that is only part of the connected component if
//     a variable is true.
//
// Edges can also be added in two ways:
//  1) Edges that are always allowed.
//  2) Edges that are only allowed when a variable is true.
//
// Apart from the unit tests, see parking.cpp, hasiwokakero.cpp,
// hitori.cpp, and nurikabe.cpp for different ways to use this
// class.
template <typename NodeIndex = std::size_t>
class ConnectedComponent {
   public:
	// Adds a node that should always be part of the connected component.
	virtual void add_node(NodeIndex node) = 0;
	// Adds a node that should be part of the connected component whenever
	// the boolean variable is 1.
	virtual void add_node(BooleanVariable part_of_connected_component) = 0;
	// Specifies one node that is guaranteed to be part of the connected
	// component. This is needed if all nodes have been added with
	// variables.
	virtual void set_first_node(BooleanVariable part_of_connected_component) = 0;
	// Specifies that one of the variables has to be part of the
	// connected component (could be all of the variables).
	virtual void set_first_node(std::vector<BooleanVariable> part_of_connected_component) = 0;

	// Adds an edge that always is allowed to connect two nodes.
	virtual void add_edge(NodeIndex node1, NodeIndex node2) = 0;
	// Adds an edge that always is allowed to connect two nodes.
	virtual void add_edge(BooleanVariable node1, BooleanVariable node2) = 0;
	// Adds an edge that is allowed to connect two nodes iff the variable
	// is equal to 1.
	virtual void add_edge(NodeIndex node1, NodeIndex node2, LogicalExpression edge_allowed) = 0;
	// Adds an edge that is allowed to connect two nodes iff the variable
	// is equal to 1.
	virtual void add_edge(BooleanVariable node1,
	                      BooleanVariable node2,
	                      LogicalExpression edge_allowed) = 0;

	// Adds the neccessary constraints to the problem.
	virtual void finish(int number_of_segments = 1) = 0;

	// Prints debug info to std::clog.
	virtual void debug(){};
};

// Implements ConnectedComponent based on flow in the graph.
//
// The connected component can not be completely empty.
template <typename NodeIndex = std::size_t>
class ConnectedComponentFlow : public ConnectedComponent<NodeIndex> {
   public:
	ConnectedComponentFlow(IP* ip);

	void add_node(NodeIndex node) override;
	void add_node(BooleanVariable part_of_connected_component) override;
	void set_first_node(BooleanVariable part_of_connected_component) override;
	void set_first_node(std::vector<BooleanVariable> part_of_connected_component) override;

	void add_edge(NodeIndex node1, NodeIndex node2) override;
	void add_edge(BooleanVariable node1, BooleanVariable node2) override;
	void add_edge(NodeIndex node1, NodeIndex node2, LogicalExpression edge_allowed) override;
	void add_edge(BooleanVariable node1,
	              BooleanVariable node2,
	              LogicalExpression edge_allowed) override;

	void finish(int number_of_segments = 1) override;

   private:
	void add_edge_internal(NodeIndex node1, NodeIndex node2, Sum edge_allowed);
	void add_initial_flow();

	// Whether an index is part of the connected component.
	std::map<NodeIndex, Sum> node_variables;
	std::map<NodeIndex, Sum> node_balance;
	std::vector<NodeIndex> start_nodes;
	std::vector<BooleanVariable> start_nodes_variables;

	// Whether an edge is allowed to be used.
	std::map<std::pair<NodeIndex, NodeIndex>, Sum> edge_variables;
	// The flow along an edge.
	std::map<std::pair<NodeIndex, NodeIndex>, Variable> all_flow;
	// The flow along an edge in the opposite direction.
	std::map<std::pair<NodeIndex, NodeIndex>, Variable> all_reverse_flow;

	bool is_finished = false;
	IP* ip;

	ConnectedComponentFlow(const ConnectedComponentFlow&) = delete;
	ConnectedComponentFlow& operator=(const ConnectedComponentFlow&) = delete;
};

// Implements ConnectedComponent with many layers of boolean variables.
//
// Forces connectivity over n nodes by adding n×n new boolean variables y.
//   • If x[i]=1, then one of y[i][k] must be 1. x[i] = ∑k y[i][k].
//   •⁠ If y[i][k]=1, y[i'][k-1] must be 1 for at least one neighbor i'.
//   • Only one of y[i][0]=1.
//
// This allows using connectivity constraints with a SAT solver. If you
// are using an IP solver, flow-based connectivity constraints are
// generally better.
template <typename NodeIndex = std::size_t>
class ConnectedComponentNoFlow : public ConnectedComponent<NodeIndex> {
   public:
	// max_with is the maximum possible width of a connected component
	// in the solution. It is good to keep this down.
	ConnectedComponentNoFlow(IP* ip, int max_width);

	void add_node(NodeIndex node) override;
	void add_node(BooleanVariable part_of_connected_component) override;
	void set_first_node(BooleanVariable part_of_connected_component) override;
	void set_first_node(std::vector<BooleanVariable> part_of_connected_component) override;

	void add_edge(NodeIndex node1, NodeIndex node2) override;
	void add_edge(BooleanVariable node1, BooleanVariable node2) override;
	void add_edge(NodeIndex node1, NodeIndex node2, LogicalExpression edge_allowed) override;
	void add_edge(BooleanVariable node1,
	              BooleanVariable node2,
	              LogicalExpression edge_allowed) override;

	void finish(int number_of_segments = 1) override;
	virtual void debug() override;

   private:
	void add_edge_internal(const std::vector<BooleanVariable>& ys1,
	                       const std::vector<BooleanVariable>& ys2);
	void add_edge_internal(const std::vector<BooleanVariable>& ys1,
	                       const std::vector<BooleanVariable>& ys2,
	                       LogicalExpression edge_allowed);

	std::map<Key<BooleanVariable>, std::vector<BooleanVariable>> helpers;
	std::map<NodeIndex, std::vector<BooleanVariable>> node_index_helpers;
	std::map<Key<BooleanVariable>, Sum> sum_of_helper_neighbors;
	IP* const ip;
	const int max_width;
	Sum level0_sum;
	bool has_start_variables = false;

	ConnectedComponentNoFlow(const ConnectedComponentNoFlow&) = delete;
	ConnectedComponentNoFlow& operator=(const ConnectedComponentNoFlow&) = delete;
};

namespace internal {
// Helper function to always get a NodeIndex from a variable. Because we
// are implementing an abstract class, also unused methods need to be
// available.
template <typename NodeIndex>
NodeIndex node_index(const BooleanVariable& x) {
	minimum_core_assert(
	    false, "NodeIndex must be size_t (default) when referencing nodes with BooleanVariables.");
	return NodeIndex();
}

template <>
inline std::size_t node_index(const BooleanVariable& x) {
	return x.get_index();
}
}  // namespace internal

template <typename NodeIndex>
ConnectedComponentFlow<NodeIndex>::ConnectedComponentFlow(IP* ip_) : ip(ip_) {}

template <typename NodeIndex>
void ConnectedComponentFlow<NodeIndex>::add_node(NodeIndex node) {
	minimum_core_assert(!is_finished);

	if (node_balance.find(node) == end(node_balance)) {
		node_balance.emplace(node, -1);
	}
	if (start_nodes.empty() == 1) {
		start_nodes.push_back(node);
	}

	// This node should always be part of the connected component, so
	// represent it with a sum equal to 1.
	node_variables[node] = 1;
}

template <typename NodeIndex>
void ConnectedComponentFlow<NodeIndex>::add_node(BooleanVariable part_of_connected_component) {
	minimum_core_assert(!is_finished);

	auto index = internal::node_index<NodeIndex>(part_of_connected_component);
	if (node_balance.find(index) == end(node_balance)) {
		node_balance.emplace(index, -part_of_connected_component);
	}
	node_variables[index] = part_of_connected_component;
}

template <typename NodeIndex>
void ConnectedComponentFlow<NodeIndex>::set_first_node(
    BooleanVariable part_of_connected_component) {
	minimum::core::check(start_nodes.empty(), "Can not set first node -- already have one.");
	start_nodes.push_back(internal::node_index<NodeIndex>(part_of_connected_component));
	ip->add_bounds(1, part_of_connected_component, 1);
}

template <typename NodeIndex>
void ConnectedComponentFlow<NodeIndex>::set_first_node(
    std::vector<BooleanVariable> input_start_nodes) {
	minimum::core::check(start_nodes.empty(), "Can not set first node -- already have one.");
	for (auto& var : input_start_nodes) {
		start_nodes.push_back(internal::node_index<NodeIndex>(var));
	}
	start_nodes_variables = std::move(input_start_nodes);
}

template <typename NodeIndex>
void ConnectedComponentFlow<NodeIndex>::add_edge(NodeIndex node1, NodeIndex node2) {
	// This edge is always allowed to be used, so represent it with
	// a sum equal to 1.
	add_edge_internal(node1, node2, 1);
}

template <typename NodeIndex>
void ConnectedComponentFlow<NodeIndex>::add_edge(BooleanVariable node1, BooleanVariable node2) {
	add_edge(internal::node_index<NodeIndex>(node1), internal::node_index<NodeIndex>(node2));
}

template <typename NodeIndex>
void ConnectedComponentFlow<NodeIndex>::add_edge(NodeIndex node1,
                                                 NodeIndex node2,
                                                 LogicalExpression edge_allowed) {
	add_edge_internal(node1, node2, edge_allowed.get_sum());
}

template <typename NodeIndex>
void ConnectedComponentFlow<NodeIndex>::add_edge(BooleanVariable node1,
                                                 BooleanVariable node2,
                                                 LogicalExpression edge_allowed) {
	minimum_core_assert(
	    (std::is_same<NodeIndex, std::size_t>::value),
	    "NodeIndex must be size_t (default) when referencing nodes with BooleanVariables.");
	add_edge(internal::node_index<NodeIndex>(node1),
	         internal::node_index<NodeIndex>(node2),
	         edge_allowed);
}

template <typename NodeIndex>
void ConnectedComponentFlow<NodeIndex>::add_edge_internal(NodeIndex node1,
                                                          NodeIndex node2,
                                                          Sum edge_allowed) {
	minimum_core_assert(!is_finished);

	core::check(node_balance.find(node1) != end(node_balance), "Add nodes before adding the edge.");
	core::check(node_balance.find(node2) != end(node_balance), "Add nodes before adding the edge.");

	edge_variables[std::make_pair(node1, node2)] = edge_allowed;

	Sum& node1_balance = node_balance[node1];
	Sum& node2_balance = node_balance[node2];

	auto flow = ip->add_variable(IP::Real);
	ip->add_bounds(0, flow, 1e100);
	auto reverse_flow = ip->add_variable(IP::Real);
	ip->add_bounds(0, reverse_flow, 1e100);
	node1_balance -= flow;
	node1_balance += reverse_flow;
	node2_balance += flow;
	node2_balance -= reverse_flow;

	all_flow[std::make_pair(node1, node2)] = flow;
	all_reverse_flow[std::make_pair(node1, node2)] = reverse_flow;
}

template <typename NodeIndex>
void ConnectedComponentFlow<NodeIndex>::finish(int number_of_segments) {
	core::check(number_of_segments == 1,
	            "ConnectedComponentFlow does not yet support number of segments > 1.");

	is_finished = true;

	add_initial_flow();

	// The balance at every node should be 0 or more.
	// Excess is allowed because then we don’t have to
	// the exact right amount at the source.
	for (const auto& balance : node_balance) {
		ip->add_constraint(balance.second >= 0);
	}

	// Add three constraints to each flow variable.
	//  1. The flow is zero if the start node is not used.
	//  2. The flow is zero if the end node is not used.
	//  3. The flow is zero if the edge is not used.
	auto max_possible_flow = node_variables.size();
	for (const auto& elem : all_flow) {
		auto node1 = node_variables[elem.first.first];
		auto node2 = node_variables[elem.first.second];
		auto flow = elem.second;
		auto edge_var = edge_variables[elem.first];
		ip->add_constraint(flow <= max_possible_flow * node1);
		ip->add_constraint(flow <= max_possible_flow * node2);
		ip->add_constraint(flow <= max_possible_flow * edge_var);
	}
	for (const auto& elem : all_reverse_flow) {
		auto node1 = node_variables[elem.first.first];
		auto node2 = node_variables[elem.first.second];
		auto reverse_flow = elem.second;
		auto edge_var = edge_variables[elem.first];
		ip->add_constraint(reverse_flow <= max_possible_flow * node1);
		ip->add_constraint(reverse_flow <= max_possible_flow * node2);
		ip->add_constraint(reverse_flow <= max_possible_flow * edge_var);
	}
}

template <>
inline void ConnectedComponentFlow<std::size_t>::add_initial_flow() {
	core::check(!start_nodes.empty(), "Need to call set_first_node.");
	double max_possible_flow = double(node_variables.size());

	// Add some flow to the first node (that we know is part of the
	// connected component).
	if (start_nodes.size() == 1) {
		node_balance[start_nodes.back()] += max_possible_flow;
	} else {
		minimum_core_assert(start_nodes.size() == start_nodes_variables.size());

		Sum sos1 = 0;
		for (auto& x : start_nodes_variables) {
			auto y = ip->add_boolean();
			sos1 += y;
			// If x is not part of the connected component, then it
			// can not get any flow.
			ip->add_constraint(x + (1 - y) >= 0);
			node_balance[x.get_index()] += max_possible_flow * y;
		}

		// If this constraint is changed, the number of components could be
		// changed to more than two.
		ip->add_constraint(sos1 == 1);
	}
}

template <typename NodeIndex>
void ConnectedComponentFlow<NodeIndex>::add_initial_flow() {
	// We know we have a single start node because we are not
	// using size_t.
	minimum_core_assert(start_nodes.size() == 1, "Need a single start node if size_t is not used.");

	// Add some flow to the first node (that we know is part of the
	// connected component).
	auto max_possible_flow = node_variables.size();
	node_balance[start_nodes.back()] += max_possible_flow;
}

template <typename NodeIndex>
ConnectedComponentNoFlow<NodeIndex>::ConnectedComponentNoFlow(IP* ip_, int max_width_)
    : ip(ip_), max_width(max_width_), level0_sum(0) {
	core::check(max_width >= 2, "Max width need to be >= 2");
}

template <typename NodeIndex>
void ConnectedComponentNoFlow<NodeIndex>::add_node(BooleanVariable x) {
	auto ys = ip->add_boolean_vector(max_width);
	Sum y_sum = 0;
	for (auto& y : ys) {
		ip->mark_variable_as_helper(y);
		sum_of_helper_neighbors[y] = 0;
		y_sum += y;
	}
	ip->add_constraint(y_sum == x);
	level0_sum += ys[0];
	helpers[x] = std::move(ys);
}

template <typename NodeIndex>
void ConnectedComponentNoFlow<NodeIndex>::add_node(NodeIndex node) {
	auto ys = ip->add_boolean_vector(max_width);
	Sum y_sum = 0;
	for (auto& y : ys) {
		ip->mark_variable_as_helper(y);
		sum_of_helper_neighbors[y] = 0;
		y_sum += y;
	}
	ip->add_constraint(y_sum == 1);
	level0_sum += ys[0];
	node_index_helpers[node] = std::move(ys);
}

template <typename NodeIndex>
void ConnectedComponentNoFlow<NodeIndex>::set_first_node(
    BooleanVariable part_of_connected_component) {
	minimum::core::check(!has_start_variables, "Can not set first node -- already have one.");
	auto& ys = helpers[part_of_connected_component];
	core::check(ys.size() == max_width, "Need to add node first.");
	ip->add_bounds(1, ys[0], 1);
	has_start_variables = true;
}

template <typename NodeIndex>
void ConnectedComponentNoFlow<NodeIndex>::set_first_node(
    std::vector<BooleanVariable> input_start_nodes) {
	minimum::core::check(!has_start_variables, "Can not set first nodes -- already have one.");
	Sum sum = 0;
	for (auto& var : input_start_nodes) {
		auto& ys = helpers[var];
		core::check(ys.size() == max_width, "Need to add node first.");
		sum += ys[0];
	}
	ip->add_constraint(sum == 1);
	has_start_variables = true;
}

template <typename NodeIndex>
void ConnectedComponentNoFlow<NodeIndex>::add_edge(NodeIndex node1, NodeIndex node2) {
	auto& ys1 = node_index_helpers[node1];
	auto& ys2 = node_index_helpers[node2];
	add_edge_internal(ys1, ys2);
}

template <typename NodeIndex>
void ConnectedComponentNoFlow<NodeIndex>::add_edge(BooleanVariable x1, BooleanVariable x2) {
	auto& ys1 = helpers[x1];
	auto& ys2 = helpers[x2];
	add_edge_internal(ys1, ys2);
}

template <typename NodeIndex>
void ConnectedComponentNoFlow<NodeIndex>::add_edge(NodeIndex node1,
                                                   NodeIndex node2,
                                                   LogicalExpression edge_allowed) {
	auto& ys1 = node_index_helpers[node1];
	auto& ys2 = node_index_helpers[node2];
	add_edge_internal(ys1, ys2, edge_allowed);
}

template <typename NodeIndex>
void ConnectedComponentNoFlow<NodeIndex>::add_edge(BooleanVariable x1,
                                                   BooleanVariable x2,
                                                   LogicalExpression edge_allowed) {
	auto& ys1 = helpers[x1];
	auto& ys2 = helpers[x2];
	add_edge_internal(ys1, ys2, edge_allowed);
}

template <typename NodeIndex>
void ConnectedComponentNoFlow<NodeIndex>::add_edge_internal(
    const std::vector<BooleanVariable>& ys1, const std::vector<BooleanVariable>& ys2) {
	core::check(ys1.size() == max_width && ys2.size() == max_width, "Need to add nodes first.");
	for (int i = 1; i < max_width; ++i) {
		sum_of_helper_neighbors[ys1[i]] += ys2[i - 1];
		sum_of_helper_neighbors[ys2[i]] += ys1[i - 1];
	}
}

template <typename NodeIndex>
void ConnectedComponentNoFlow<NodeIndex>::add_edge_internal(const std::vector<BooleanVariable>& ys1,
                                                            const std::vector<BooleanVariable>& ys2,
                                                            LogicalExpression edge_allowed) {
	core::check(ys1.size() == max_width && ys2.size() == max_width, "Need to add nodes first.");
	for (int i = 1; i < max_width; ++i) {
		auto z1 = ip->add_boolean();
		ip->add_constraint(z1 - edge_allowed.get_sum() <= 0);
		ip->add_constraint(z1 - ys2[i - 1] <= 0);
		ip->add_constraint(edge_allowed.get_sum() + ys2[i - 1] - z1 <= 1);
		sum_of_helper_neighbors[ys1[i]] += z1;

		auto z2 = ip->add_boolean();
		ip->add_constraint(z2 - edge_allowed.get_sum() <= 0);
		ip->add_constraint(z2 - ys1[i - 1] <= 0);
		ip->add_constraint(edge_allowed.get_sum() + ys1[i - 1] - z2 <= 1);
		sum_of_helper_neighbors[ys2[i]] += z2;
	}
}

template <typename NodeIndex>
void ConnectedComponentNoFlow<NodeIndex>::finish(int number_of_segments) {
	ip->add_constraint(level0_sum == number_of_segments);

	for (auto& elem : helpers) {
		auto& ys = elem.second;
		for (int i = 1; i < max_width; ++i) {
			// If y[i]=1, then at least one of its neigbors (at level i-1)
			// has to be 1.
			ip->add_constraint(sum_of_helper_neighbors[ys[i]] - ys[i] >= 0);
		}
	}

	for (auto& elem : node_index_helpers) {
		auto& ys = elem.second;
		for (int i = 1; i < max_width; ++i) {
			// If y[i]=1, then at least one of its neigbors (at level i-1)
			// has to be 1.
			ip->add_constraint(sum_of_helper_neighbors[ys[i]] - ys[i] >= 0);
		}
	}
}

template <typename NodeIndex>
void ConnectedComponentNoFlow<NodeIndex>::debug() {
	for (auto& elem : helpers) {
		auto& x = elem.first.get_var();
		std::clog << "Variable: " << x.bool_value() << ": ";
		auto& ys = elem.second;
		for (auto& y : ys) {
			std::clog << y.bool_value() << " ";
		}
		std::clog << "\n";
	}
	for (auto& elem : node_index_helpers) {
		auto& x = elem.first;
		std::clog << "Variable " << core::to_string(x) << ": ";
		auto& ys = elem.second;
		for (auto& y : ys) {
			std::clog << y.bool_value() << " ";
		}
		std::clog << "\n";
	}
}
}  // namespace linear
}  // namespace minimum
