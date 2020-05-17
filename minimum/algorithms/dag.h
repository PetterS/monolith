#pragma once
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <utility>
#include <vector>

#include <minimum/core/check.h>
#include <minimum/core/grid.h>
#include <minimum/core/range.h>

namespace minimum {
namespace algorithms {

template <typename T, int num_weights>
class DAGEdge {
	static_assert(num_weights > 0);

   public:
	DAGEdge(int to_, T cost_) : to(to_), cost(cost_) {
		for (int i = 0; i < num_weights; ++i) {
			weights[i] = 0;
		}
	}
	int to;
	T cost;
	std::array<int, num_weights> weights;
	bool operator<(const DAGEdge& other) const { return to < other.to; }
	bool operator==(const DAGEdge& other) const { return to == other.to; }
};

template <typename T>
class DAGEdge<T, 0> {
   public:
	DAGEdge(int to_, T cost_) : to(to_), cost(cost_) {}
	int to;
	T cost;
	bool operator<(const DAGEdge& other) const { return to < other.to; }
	bool operator==(const DAGEdge& other) const { return to == other.to; }
};

// A DAG that is already topologically sorted, i.e. only has edges
// from i → j, where i < j.
template <typename T = double, int num_weights = 1, int num_edge_weights = 0>
class SortedDAG {
	static_assert(num_weights >= 0, "num_weights can not be negative.");

   public:
	using Edge = DAGEdge<T, num_edge_weights>;
	struct Node {
		T cost = 0;
		std::array<int, num_weights> weights = {};
		std::vector<Edge> edges;
	};

	SortedDAG(int num_nodes) : nodes(num_nodes) {}

	// A copy of master_dag with only a subset of the nodes.
	SortedDAG(const SortedDAG<T, num_weights>& master_dag, int begin, int end) {
		for (int i = begin; i < end; ++i) {
			nodes.push_back(master_dag.nodes[i]);
		}
		nodes.emplace_back();
		for (auto& node : nodes) {
			for (auto& edge : node.edges) {
				edge.to -= begin;
				if (edge.to >= size() - 1) {
					edge.to = size() - 1;
				}
			}
		}
	}

	int size() const { return nodes.size(); }

	const Node& get_node(int i) const { return nodes.at(i); }

	void set_node_cost(int i, T new_cost) { nodes.at(i).cost = new_cost; }

	void change_node_cost(int i, T delta) { nodes.at(i).cost += delta; }

	void set_node_weight(int i, int w, int new_weight) { nodes.at(i).weights[w] = new_weight; }

	void disconnect_node(int i) {
		minimum_core_assert(i >= 0 && i < nodes.size());
		for (int j = 0; j < i; ++j) {
			auto& edges = nodes[j].edges;
			edges.erase(
			    std::remove_if(edges.begin(), edges.end(), [i](auto& e) { return e.to == i; }),
			    edges.end());
		}
		nodes[i].edges.clear();
	}

	Edge& add_edge(int i, int j, T cost = 0) {
		minimum_core_assert(i >= 0 && j >= 0 && i < nodes.size() && j < nodes.size());
		minimum::core::check(i < j, "Graph needs to be topologically sorted.");
		return nodes[i].edges.emplace_back(j, cost);
	}

	std::string str() const {
		std::string out;
		for (int i : minimum::core::range(nodes.size())) {
			std::vector<int> dest;
			for (auto& e : nodes[i].edges) {
				dest.push_back(e.to);
			}
			out += minimum::core::to_string(i, ": ", dest, "\n");
		}
		return out;
	}

	// Mapping from and to a new set of indices.
	struct Translator {
		bool made_changes = false;
		std::vector<int> old_to_new;
		std::vector<std::vector<int>> new_to_old;
	};

	// Reduces the graph, i.e. removes all nodes that can not be part of a path
	// from the first to the last node.
	//
	// Returns a translator from and to the new set of indices. If no changes
	// were mode, made_changes will be false in the translator.
	Translator reduce_graph() {
		if (size() <= 2) {
			return {};
		}

		// Find the forward reachable nodes.
		std::vector<int> forward_reachable(nodes.size(), 0);
		std::vector<int> stack;
		stack.push_back(0);
		int forward_nodes_found = 0;
		while (!stack.empty()) {
			int i = stack.back();
			stack.pop_back();
			if (forward_reachable[i] == 1) {
				continue;
			}
			forward_nodes_found++;
			forward_reachable[i] = 1;
			for (auto& edge : nodes[i].edges) {
				stack.push_back(edge.to);
			}
		}
		minimum_core_assert(forward_reachable[nodes.size() - 1], "There is no path.");

		// Find the backwards reachable nodes.
		std::vector<int> backward_reachable(nodes.size(), 0);
		int backwards_nodes_found = 1;
		backward_reachable[nodes.size() - 1] = 1;
		for (int i = nodes.size() - 2; i >= 0; --i) {
			for (auto& edge : nodes[i].edges) {
				if (backward_reachable[edge.to]) {
					backward_reachable[i] = 1;
					backwards_nodes_found++;
					break;
				}
			}
		}

		// Return if all nodes are reachable.
		if (forward_nodes_found == nodes.size() && backwards_nodes_found == nodes.size()) {
			return {};
		}

		// Create the new set of nodes from those that are forward and
		// backwards reachable.
		Translator translator;
		translator.made_changes = true;
		translator.old_to_new.resize(nodes.size(), -1);
		std::vector<Node> new_nodes;
		for (int i : minimum::core::range(nodes.size())) {
			if (forward_reachable[i] == 1 && backward_reachable[i] == 1) {
				translator.old_to_new[i] = new_nodes.size();
				new_nodes.emplace_back(std::move(nodes[i]));
			}
		}

		// Now switch to the new set of nodes and translate all the edges.
		nodes = std::move(new_nodes);
		for (auto& node : nodes) {
			for (auto& edge : node.edges) {
				edge.to = translator.old_to_new[edge.to];
			}
			node.edges.erase(
			    std::remove_if(
			        node.edges.begin(), node.edges.end(), [](auto& e) { return e.to < 0; }),
			    node.edges.end());
		}

		translator.new_to_old.resize(nodes.size());
		for (int i : minimum::core::range(translator.old_to_new.size())) {
			if (translator.old_to_new[i] >= 0) {
				minimum_core_assert(translator.new_to_old[translator.old_to_new[i]].empty());
				translator.new_to_old[translator.old_to_new[i]].push_back(i);
			}
		}

		return translator;
	}

	// Merges nodes that have to be traversed anyway. Simply adds the weights, so is
	// not always appropriate.
	// Edges are merged without taking their cost into account.
	Translator merge_graph(std::function<bool(int)> can_merge) {
		Translator translator;

		std::vector<int> mergeable(nodes.size(), 1);
		for (int i : minimum::core::range(size())) {
			for (auto& edge : nodes[i].edges) {
				for (int j = i + 1; j < edge.to; ++j) {
					mergeable[j] = 0;
				}
			}
		}

		std::vector<Node> new_nodes;
		translator.old_to_new.resize(size(), -1);
		for (int i = 0; i < size(); ++i) {
			if (mergeable[i] == 0 || !can_merge(i)) {
				// Simply add this node.
				new_nodes.emplace_back(std::move(nodes[i]));
				translator.old_to_new[i] = new_nodes.size() - 1;
			} else {
				translator.made_changes = true;

				new_nodes.emplace_back();
				for (; i < size() && mergeable[i] == 1 && can_merge(i); ++i) {
					translator.old_to_new[i] = new_nodes.size() - 1;
					// Add edges
					for (auto& edge : nodes[i].edges) {
						new_nodes.back().edges.push_back(edge);
					}
					for (int w = 0; w < num_weights; ++w) {
						new_nodes.back().weights[w] += nodes[i].weights[w];
					}
				}
				// Decrement i, because it is being incremented in
				// both for loops.
				--i;
			}
		}

		// Now switch to the new set of nodes and translate all the edges.
		nodes = std::move(new_nodes);
		for (int i = 0; i < size(); ++i) {
			auto& node = nodes[i];
			for (auto& edge : node.edges) {
				edge.to = translator.old_to_new[edge.to];
			}
			std::sort(node.edges.begin(), node.edges.end());
			node.edges.erase(std::unique(node.edges.begin(), node.edges.end()), node.edges.end());
			node.edges.erase(
			    std::remove_if(
			        node.edges.begin(), node.edges.end(), [i](auto& edge) { return edge.to == i; }),
			    node.edges.end());
		}

		translator.new_to_old.resize(nodes.size());
		for (int i : minimum::core::range(translator.old_to_new.size())) {
			if (translator.old_to_new[i] >= 0) {
				translator.new_to_old[translator.old_to_new[i]].push_back(i);
			}
		}

		return translator;
	}

	void write_dot(std::ostream& out,
	               std::function<std::string(int)> node_name,
	               std::function<std::string(int)> node_color,
	               std::function<int(int)> rank = nullptr) const {
		out << "digraph dag {\n";
		for (auto i : minimum::core::range(nodes.size())) {
			for (auto& edge : nodes[i].edges) {
				out << i << " -> " << edge.to << ";\n";
			}
		}

		for (auto i : minimum::core::range(nodes.size())) {
			out << i << "[label=\"" << node_name(i)
			    << "\", style=filled, fillcolor=" << node_color(i) << "];\n";
		}

		if (rank) {
			std::map<int, std::vector<int>> ranks;
			for (auto i : minimum::core::range(nodes.size())) {
				ranks[rank(i)].push_back(i);
			}
			out << "edge[style=invis,dir=none];\n";
			for (auto& r : ranks) {
				out << "{ rank=same; ";
				bool first = true;
				for (auto i : r.second) {
					if (!first) {
						out << " -> ";
					}
					first = false;
					out << i;
				}
				out << "}\n";
			}
		}

		out << "}\n";
	}

   private:
	std::vector<Node> nodes;
};

template <typename T = double>
struct SolutionEntry {
	int prev = -1;
	T cost = std::numeric_limits<T>::max() / 2;
	bool operator==(const SolutionEntry& rhs) const { return prev == rhs.prev; }
};

template <typename T, int num_weights, int num_edge_weights>
T solution_cost(const SortedDAG<T, num_weights, num_edge_weights>& dag,
                const std::vector<int>& solution) {
	if (solution.empty()) {
		return 0;
	}
	auto cost = dag.get_node(solution[0]).cost;
	for (auto i : minimum::core::range(std::size_t(1), solution.size())) {
		auto from = solution[i - 1];
		auto to = solution[i];
		for (auto& edge : dag.get_node(from).edges) {
			if (edge.to == to) {
				cost += dag.get_node(to).cost + edge.cost;
				break;
			}
		}
	}
	return cost;
}

template <typename T, int num_weights, int num_edge_weights>
T shortest_path(const SortedDAG<T, num_weights, num_edge_weights>& dag,
                std::vector<SolutionEntry<T>>* solution_ptr) {
	using minimum::core::range;

	auto& solution = *solution_ptr;
	solution.clear();
	if (dag.size() == 0) {
		return 0;
	}

	solution.resize(dag.size());
	solution[0].cost = dag.get_node(0).cost;
	for (auto i : range(dag.size())) {
		for (auto& edge : dag.get_node(i).edges) {
			auto node_cost = dag.get_node(edge.to).cost;
			auto cost = solution[i].cost + node_cost + edge.cost;
			if (cost < solution[edge.to].cost) {
				solution[edge.to].cost = cost;
				solution[edge.to].prev = i;
			}
		}
	}
	return solution.back().cost;
}

template <typename T, int num_weights>
T resource_constrained_shortest_path(SortedDAG<T, num_weights, 0> dag,
                                     std::array<T, num_weights> upper_bounds,
                                     std::vector<std::vector<SolutionEntry<T>>>* solutions) {
	using minimum::core::range;

	constexpr bool verbose = false;
	if (verbose) {
		std::cerr << "\n\n";
	}

	static_assert(num_weights >= 1, "Need at least one resource constraint.");
	solutions->clear();
	solutions->emplace_back();

	// Save the original node costs so we can restore them later.
	std::vector<T> org_costs(dag.size());
	for (int i : range(dag.size())) {
		org_costs[i] = dag.get_node(i).cost;
	}

	// Uses the same relaxation as in
	// http://faculty.nps.edu/kwood/docs/ORSNZ03CarlyleWood.pdf
	//
	// We want to solve
	//  min c'x
	//      x forms a path through the DAG.
	//      r_i'x <= b_i,  for i=1..num_weights.
	//
	// The dual function is
	//
	// d(λ) = inf_{x forms a path}( c'x  + ∑ λ_i·(r_i'x - b_i) )
	//
	// The subgradient is g_i = r_i'x*.

	// The dual variables.
	std::array<T, num_weights> lambda{};

	for (int iter = 1; iter <= 20; ++iter) {
		T lambda_sum = 0;  // ∑ λ_i·b_i
		for (int w : range(num_weights)) {
			lambda_sum += lambda[w] * upper_bounds[w];
		}

		// Solve the inf_x problem. (The term ∑ λ_i·b_i does not depend
		// on x and is left out here.)
		auto& solution = solutions->back();
		auto solution_value = shortest_path(dag, &solution);

		// Form the LHS of all constraints.
		std::array<T, num_weights> lhs{};
		int i = dag.size() - 1;
		while (i >= 0) {
			for (int w : range(num_weights)) {
				lhs[w] += dag.get_node(i).weights[w];
			}
			i = solution[i].prev;
		}

		// Step size for updating λ.
		auto tau = 0.5 / (iter + 1);

		// Update λ and check if the solution is feasible.
		bool feasible = true;
		for (int w = 0; w < num_weights; ++w) {
			if (lhs[w] > upper_bounds[w]) {
				feasible = false;
			}

			auto super_gradient = lhs[w] - upper_bounds[w];
			lambda[w] += tau * super_gradient;
			lambda[w] = std::max(T(0), lambda[w]);
		}

		// Set the costs of the shortest path problem to c'x + ∑ λ_i·(r_i'x).
		for (int i : range(dag.size())) {
			T node_cost = org_costs[i];
			for (int w : range(num_weights)) {
				node_cost += dag.get_node(i).weights[w] * lambda[w];
			}
			dag.set_node_cost(i, node_cost);
		}

		if (verbose) {
			std::cerr << "dual=" << solution_value - lambda_sum << ", "
			          << minimum::core::to_string(lhs)
			          << " <= " << minimum::core::to_string(upper_bounds)
			          << ", updated lambda=" << minimum::core::to_string(lambda);
			if (feasible) {
				T cost = 0;
				int i = dag.size() - 1;
				while (i >= 0) {
					cost += org_costs[i];
					i = solutions->back()[i].prev;
				}
				std::cerr << ", feasible=" << cost;
			}
			std::cerr << "\n";
		}

		if (feasible) {
			// Only save this solution if it is different from the previous one.
			if (solutions->size() >= 2) {
				if (solutions->at(solutions->size() - 1) != solutions->at(solutions->size() - 2)) {
					solutions->emplace_back();
				} else {
					solutions->back().clear();
				}
			} else {
				solutions->emplace_back();
			}
		}
	}

	// Restore original costs.
	for (int i : range(dag.size())) {
		dag.set_node_cost(i, org_costs[i]);
	}
	// Remove temp solution space.
	solutions->pop_back();
	if (solutions->empty()) {
		return 0;
	}

	// Compute cost of final feasible solution with the original costs.
	int i = dag.size() - 1;
	T cost = 0;
	while (i >= 0) {
		cost += dag.get_node(i).cost;
		i = solutions->back()[i].prev;
	}

	return cost;
}

template <typename T, int num_weights, int num_edge_weights>
T resource_constrained_shortest_path(const SortedDAG<T, num_weights, num_edge_weights>& dag,
                                     int lower_bound,
                                     int upper_bound,
                                     std::vector<int>* solution_ptr) {
	static_assert(num_weights >= 1, "Need weights for resource constraints.");
	static_assert(num_edge_weights <= 1,
	              "Edge weights for consecutive constraint is not supported.");
	using minimum::core::check;
	using minimum::core::range;

	check(lower_bound <= upper_bound, "Invalid bounds: ", lower_bound, " > ", upper_bound);
	auto& solution = *solution_ptr;
	solution.clear();
	if (dag.size() == 0) {
		return 0;
	} else if (dag.size() == 1) {
		solution.emplace_back(0);
		auto w = dag.get_node(0).weights[0];
		check(lower_bound <= w && w <= upper_bound, "Infeasible.");
		return dag.get_node(0).cost;
	}

	// partial[i][c] is the smallest cost to reach node i consuming c
	// resource.
	auto partial = minimum::core::make_grid<SolutionEntry<T>>(dag.size(), upper_bound + 1);

	partial[0][dag.get_node(0).weights[0]].cost = dag.get_node(0).cost;

	for (auto i : range(dag.size())) {
		for (auto c : range(upper_bound + 1)) {
			if (i > 0 && partial[i][c].prev < 0) {
				// We cannot reach node i with c weight.
				continue;
			}
			for (auto& edge : dag.get_node(i).edges) {
				// Weight of the path up to and including the edge’s destination.
				auto weight = c + dag.get_node(edge.to).weights[0];
				if constexpr (num_edge_weights > 0) {
					weight += edge.weights[0];
				}
				check(weight >= 0, "Negative resource encountered along path.");

				if (weight > upper_bound) {
					continue;
				}
				// Cost of the path up to and including the edge’s destination.
				auto cost = partial[i][c].cost + dag.get_node(edge.to).cost + edge.cost;
				if (cost < partial[edge.to][weight].cost) {
					partial[edge.to][weight].cost = cost;
					partial[edge.to][weight].prev = i;
				}
			}
		}
	}

	T best = std::numeric_limits<T>::max();
	int best_c = -1;
	for (auto c : range(lower_bound, upper_bound + 1)) {
		if (partial.back()[c].prev >= 0) {
			// std::cerr << "-- There is a path with cost " << partial.back()[c].cost
			//           << " and " << c << " weight.\n";
			auto cost = partial.back()[c].cost;
			if (cost < best) {
				best = cost;
				best_c = c;
			}
		}
	}
	check(best_c >= 0, "Could not find a feasible path.");

	int i = dag.size() - 1;
	solution.push_back(i);
	int current_c = best_c;
	while (true) {
		if (partial[i][current_c].prev < 0) {
			break;
		}

		int i_prev = i;
		i = partial[i][current_c].prev;
		current_c -= dag.get_node(i_prev).weights[0];
		if constexpr (num_edge_weights > 0) {
			for (auto& edge : dag.get_node(i).edges) {
				if (edge.to == i_prev) {
					current_c -= edge.weights[0];
					break;
				}
			}
		}
		solution.push_back(i);
	}
	reverse(solution.begin(), solution.end());
	return partial.back()[best_c].cost;
}

namespace internal {
template <typename T>
struct Entry {
	std::tuple<int, int, int> prev = std::make_tuple(-1, -1, -1);
	T cost = std::numeric_limits<T>::max();
	bool valid() { return get<0>(prev) >= 0; }
};
}  // namespace internal

template <typename T, int num_weights, int num_edge_weights>
T resource_constrained_shortest_path(const SortedDAG<T, num_weights, num_edge_weights>& dag,
                                     int lower_bound,
                                     int upper_bound,
                                     int min_consecutive,
                                     int max_consecutive,
                                     std::vector<int>* solution_ptr) {
	static_assert(num_weights >= 2, "Need weights for resource and consecutive constraints.");
	static_assert(num_edge_weights <= 1,
	              "Edge weights for consecutive constraint is not supported.");
	using minimum::core::check;
	using minimum::core::range;

	auto& solution = *solution_ptr;
	solution.clear();
	if (dag.size() == 0) {
		return 0;
	} else if (dag.size() == 1) {
		solution.emplace_back(0);
		auto w = dag.get_node(0).weights[0];
		check(lower_bound <= w && w <= upper_bound, "Infeasible.");
		return dag.get_node(0).cost;
	}
	// We allow a negative number as lower bound, but set it to zero since
	// weights are not allowed to go negative anyway.
	lower_bound = std::max(lower_bound, 0);
	minimum_core_assert(max_consecutive >= 1);

	// partial(i, c, d) is the smallest cost to reach node i consuming c
	// resource ending in d consecutive active nodes.
	minimum::core::Grid3D<internal::Entry<T>, int> partial(
	    dag.size(), upper_bound + 1, max_consecutive + 1);

	partial(0, dag.get_node(0).weights[0], 0).cost = dag.get_node(0).cost;

	for (auto i : range(dag.size())) {
		for (auto c : range(upper_bound + 1)) {
			for (auto d : range(max_consecutive + 1)) {
				if (i > 0 && !partial(i, c, d).valid()) {
					// We cannot reach node i with c weight and d consecutive.
					continue;
				}
				for (auto& edge : dag.get_node(i).edges) {
					auto& to_node = dag.get_node(edge.to);

					// Weight of the path up to and including the edge’s destination.
					auto weight = c + to_node.weights[0];
					if constexpr (num_edge_weights > 0) {
						weight += edge.weights[0];
					}
					if (weight < 0) {
						throw std::runtime_error("Negative resource encountered along path.");
					}
					if (weight > upper_bound) {
						// This path is not allowed.
						continue;
					}

					auto this_consecutive = to_node.weights[1];
					// Number of consecutive costs up to and including the edge’s
					// destination.
					auto consecutive = d + this_consecutive;
					if (this_consecutive == 0) {
						// If this node stops a segment with consecutive cost, we
						// should make sure the minimum is fulfilled.
						if (consecutive > 0 && consecutive < min_consecutive) {
							// This path is not allowed.
							continue;
						}
						consecutive = 0;
					} else if (this_consecutive > 0) {
						if (consecutive > max_consecutive) {
							// This path is not allowed.
							continue;
						}
					} else {
						// Reset the consecutive counter and always allow
						// this path.
						consecutive = 0;
					}

					// Cost of the path up to and including the edge’s destination.
					auto cost = partial(i, c, d).cost + to_node.cost + edge.cost;
					if (cost < partial(edge.to, weight, consecutive).cost) {
						partial(edge.to, weight, consecutive).cost = cost;
						partial(edge.to, weight, consecutive).prev = std::make_tuple(i, c, d);
					}
				}
			}
		}
	}

	// Find the best feasible path.
	T best = std::numeric_limits<T>::max();
	int best_c = -1;
	int best_d = -1;
	for (auto c : range(lower_bound, upper_bound + 1)) {
		for (auto d : range(max_consecutive + 1)) {
			if (partial(dag.size() - 1, c, d).valid()) {
				auto cost = partial(dag.size() - 1, c, d).cost;
				if (cost < best) {
					best = cost;
					best_c = c;
					best_d = d;
				}
			}
		}
	}
	check(best_c >= 0, "Could not find a feasible path.");

	// Recover solution.
	int i = dag.size() - 1;
	solution.push_back(i);
	int c = best_c;
	int d = best_d;
	while (i >= 0) {
		const auto& prev = partial(i, c, d).prev;
		if (get<0>(prev) < 0) {
			break;
		}
		solution.push_back(get<0>(prev));
		std::tie(i, c, d) = prev;
	}
	reverse(solution.begin(), solution.end());
	return best;
}

template <typename T, int num_weights, int num_edge_weights>
void resource_constrained_shortest_path_partial(
    const SortedDAG<T, num_weights, num_edge_weights>& dag,
    int lower_bound,
    int upper_bound,
    int min_consecutive,
    int max_consecutive,
    int num_splits,
    std::vector<int>* solution_ptr) {
	minimum::core::check(
	    1 <= num_splits && 5 * num_splits <= dag.size(), "Can not split problem into ", num_splits);
	if (num_splits == 1) {
		resource_constrained_shortest_path(
		    dag, lower_bound, upper_bound, min_consecutive, max_consecutive, solution_ptr);
		return;
	}
	solution_ptr->clear();

	int s = 0;
	int t = dag.size() / num_splits;
	std::vector<int> subsolution;
	for (int part = 0; part < num_splits; ++part) {
		// Create and solve subproblem.
		SortedDAG<T, num_weights, num_edge_weights> subdag(dag, s, t);
		int lb = ceil(double(lower_bound) / num_splits);
		int ub = floor(double(upper_bound) / num_splits);

		if (part == num_splits - 1) {
			int resource = 0;
			for (auto i : *solution_ptr) {
				resource += dag.get_node(i).weights[0];
			}
			lb = lower_bound - resource;
			ub = upper_bound - resource;
		}

		resource_constrained_shortest_path(
		    subdag, lb, ub, min_consecutive, max_consecutive, &subsolution);
		subsolution.pop_back();

		// Translate solution.
		for (auto i : subsolution) {
			int global_index = i + s;
			if (!solution_ptr->empty() && global_index == solution_ptr->back()) {
				continue;
			}
			solution_ptr->push_back(global_index);
		}

		if (part < num_splits - 1) {
			// Move window to next problem.
			// Pick the last node with zero consecutive penalty to start the next
			// segment.
			s = solution_ptr->back();
			while (dag.get_node(s).weights[1] > 0) {
				solution_ptr->pop_back();
				s = solution_ptr->back();
			}
			t += dag.size() / num_splits;
			if (part == num_splits - 2) {
				t = dag.size() - 1;
			}
		}
	}

	// std::cerr << "F: " << to_string(*solution_ptr) << ".\n";
	int resource = 0;
	int consecutive = 0;
	for (auto i : *solution_ptr) {
		resource += dag.get_node(i).weights[0];
		if (dag.get_node(i).weights[1] < 0) {
			consecutive = 0;
		} else if (dag.get_node(i).weights[1] == 0) {
			minimum_core_assert(consecutive == 0 || min_consecutive <= consecutive,
			                    "Solution min consecutive not feasible: ",
			                    consecutive);
			consecutive = 0;
		} else {
			consecutive += dag.get_node(i).weights[1];
			minimum_core_assert(consecutive <= max_consecutive,
			                    "Solution max consecutive not feasible: ",
			                    consecutive);
		}
	}
	minimum_core_assert(lower_bound <= resource && resource <= upper_bound,
	                    "Solution resource not feasible: ",
	                    lower_bound,
	                    " <= ",
	                    resource,
	                    " <= ",
	                    upper_bound);
}
}  // namespace algorithms
}  // namespace minimum
