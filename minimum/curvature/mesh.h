// Petter Strandmark 2013.
#ifndef CURVE_EXTRACTION_MESH_H
#define CURVE_EXTRACTION_MESH_H

#include <functional>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <minimum/curvature/export.h>

namespace minimum {
namespace curvature {

class RESEARCH_CURVATURE_API Mesh {
   public:
	// Represents a 3D point in the mesh. Contains
	// an adjacency list, representing the undirected
	// added to the mesh.
	struct Point {
		Point(float x, float y, float z);
		float x, y, z;
		std::vector<int> adjacent_points;

		bool operator<(const Point& rhs) const;
	};

	typedef std::pair<int, int> Edge;
	typedef std::tuple<int, int, int> EdgePair;

	Mesh() {}
	~Mesh();

	Mesh(const Mesh&) = delete;
	Mesh& operator=(const Mesh&) = delete;

	// Adds a point to the mesh. If the same point is
	// added twice, only one point will be stored.
	void add_point(float x, float y, float z);

	// Returns the index of a point or -1 if it is
	// not present in the mesh.
	// The non-const version sorts the points and forbids the addition of any
	// subsequent points. The const version throws if the points are not
	// sorted.
	int find_point(float x, float y, float z) const;
	int find_point(float x, float y, float z);

	int number_of_points() const;
	// Returns a point in the mesh. Throws an exception
	// if the index is invalid.
	const Point& get_point(int p) const;

	// Transforms all points in the mesh according to
	//   x = tx + sx * x
	//   x = ty + sy * y
	//   x = tz + sz * z
	void transform_points(float tx, float ty, float tz, float sx, float sy, float sz);

	// Adds an edge to the mesh. The points must have been
	// previously added.
	//
	// NOTE: After adding edges, no more points can be
	//       added.
	void add_edge(float x1, float y1, float z1, float x2, float y2, float z2);

	// Adds edges between all points that are no more
	// than dmax apart.
	void add_edges(float dmax);

	// Same as above, but ignores edges according to
	// the ignore function.
	void add_edges(float dmax,
	               const std::function<bool(float, float, float, float, float, float)>& ignore);

	// Returns the number of edges in the mesh.
	int number_of_edges() const;

	// Returns a specific edge. Throws an exception
	// if the index is invalid.
	const Edge& get_edge(int e) const;

	// Finds an edge. Throws an exception
	// if the edge is not found.
	int find_edge(int p1, int p2) const;

	// Computes all edges adjacent to a given edge.
	void get_adjacent_edges(int e, std::vector<int>* adjacent) const;

	// Finished the mesh, i.e. builds adjacency lists
	// for edges and creates edge pairs. After finishing,
	// the mesh can not be modified.
	void finish(bool create_pairs = false);

	// Returns the maximum node degree.
	int get_connectivity() const;

	// The number of edge pairs in the mesh.
	int number_of_edge_pairs() const;

	const EdgePair& get_edge_pair(int ep) const;
	// Finds an edge pair. Throws an exception
	// if the pair is not found.
	int find_edge_pair(int p1, int p2, int p3) const;

	// Computes all adjacent pairs to a given pair.
	void get_adjacent_pairs(int ep, std::vector<int>* adjacent) const;

	// Length of an edge.
	float edge_length(int) const;

	// Starts an SVG file, by drawing the entire mesh.
	// After starting, additional curves may be drawn
	// in the same file by using draw_path.
	void start_SVG(const std::string& filename);

	// Same as above, but the function edge_style
	// specifies the SVG style of each edge.
	void start_SVG(const std::string& filename, const std::function<std::string(int)>& edge_style);

	// Draws a path with a specific color in the started SVG file.
	void draw_path(const std::vector<int>& path, const std::string& color = "ff0000");
	// Finished the SVG file. Will be called automatically
	// when the mesh is destructed.
	void end_SVG();

	// Input vector contains edge indices.
	std::vector<Point> edgepath_to_points(const std::vector<int>&) const;

	// Input vector contains edge pair indices.
	std::vector<Point> pairpath_to_points(const std::vector<int>&) const;

   private:
	void finalize_points();

	bool finished = false;
	bool no_more_points = false;
	std::vector<Point> points;
	std::vector<Edge> edges;
	std::vector<EdgePair> edge_pairs;

	int connectivity = -1;

	std::ofstream* fout = nullptr;
};
}  // namespace curvature
}  // namespace minimum

#endif
