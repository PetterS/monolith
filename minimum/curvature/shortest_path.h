// Petter Strandmark.
//
// Computes the shortest path between two sets of
// nodes in a dynamically generated graph. Uses
// Dijkstra's algorithm or A*.
//
#ifndef CURVE_EXTRACTION_SHORTEST_PATH_H
#define CURVE_EXTRACTION_SHORTEST_PATH_H

#include <functional>
#include <set>
#include <vector>

#include <minimum/curvature/export.h>

namespace minimum {
namespace curvature {

struct Neighbor {
	Neighbor() {}

	Neighbor(int dest, double d) {
		this->destination = dest;
		this->distance = d;
	}

	Neighbor& operator=(const Neighbor& right) {
		this->destination = right.destination;
		this->distance = right.distance;
		return *this;
	}

	int destination;
	double distance;
};

struct ShortestPathOptions {
	ShortestPathOptions()
	    : print_progress(false),
	      maximum_queue_size(0),
	      compute_all_distances(false),
	      store_visited(false),
	      store_parents(false) {}
	// Prints progress to stderr about the number of
	// nodes visited.
	bool print_progress;
	// If the queue grows beyond this size, the algorithm
	// terminates with an error.
	std::size_t maximum_queue_size;
	// The distance to all nodes from the start set is
	// calculated. The end set can then be empty.
	// The distances are stored in the vector.
	bool compute_all_distances;
	mutable std::vector<float> distance;
	// Store the time at which all nodes in the graph
	// were visited.
	bool store_visited;
	mutable std::vector<int> visit_time;
	// Store the shortest path graph, i.e. the parent
	// of every node.
	// This options is probably most useful in
	// combination with compute_all_distances.
	bool store_parents;
	mutable std::vector<int> parents;
};

// Computes the shortest path between two sets of nodes in a graph.
// The graph does not have to be explicitly known; it can be computed
// on the fly. Uses Dijkstra's algorithm, or A* if a heuristic is
// provided.
RESEARCH_CURVATURE_API double shortest_path(  // The number of nodes in the graph.
    int n,
    // The set of nodes from which to compute the
    // path.
    const std::set<int>& start_set,
    // The set of nodes where the path can end.
    const std::set<int>& end_set,
    // Oracle which returns the set of neighbors of
    // a given node, along with their distances.
    const std::function<void(int, std::vector<Neighbor>* neighbors)>& get_neighbors,
    // (output) Recieves the shortest path.
    std::vector<int>* path,
    // (optional) Heuristic function for A*
    const std::function<double(int)>* get_lower_bound = 0,
    // (optional) Options to the solver.
    const ShortestPathOptions& options = ShortestPathOptions());

// Identical function except that the last argument is
// a const reference instead of a pointer.
RESEARCH_CURVATURE_API double shortest_path(
    int n,
    const std::set<int>& start_set,
    const std::set<int>& end_set,
    const std::function<void(int, std::vector<Neighbor>* neighbors)>& get_neighbors,
    std::vector<int>* path,
    const std::function<double(int)>& get_lower_bound,
    const ShortestPathOptions& options = ShortestPathOptions());

RESEARCH_CURVATURE_API double bidirectional_shortest_path(
    // The number of nodes in the graph.
    int n,
    // The set of nodes from which to compute the
    // path.
    const std::set<int>& start_set,
    // The set of nodes where the path can end.
    const std::set<int>& end_set,
    // Oracle which returns the set of neighbors of
    // a given node, along with their distances.
    const std::function<void(int, std::vector<Neighbor>* neighbors)>& get_neighbors,
    // (output) Recieves the shortest path.
    std::vector<int>* path);

//
// Same functionality as above, but will use less memory per node added
// to the priority queue. On the other hand, it will always allocate
// 12*n bytes of memory.
//
// This function will allocate no more than 16*n + O(1) bytes.
//
// Does not support lower bound heuristics.
//
// ** NOT THREAD-SAFE **
//
RESEARCH_CURVATURE_API double shortest_path_memory_efficient(
    int n,
    const std::set<int>& start_set,
    const std::set<int>& end_set,
    const std::function<void(int, std::vector<Neighbor>* neighbors)>& neighbors,
    std::vector<int>* path,
    const ShortestPathOptions& options = ShortestPathOptions());
}  // namespace curvature
}  // namespace minimum

#endif
