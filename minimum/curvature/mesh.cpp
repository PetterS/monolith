// Petter Strandmark 2013.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <tuple>

using std::ignore;
using std::tie;

#include <minimum/curvature/mesh.h>

namespace minimum {
namespace curvature {

Mesh::Point::Point(float x, float y, float z) {
	this->x = x;
	this->y = y;
	this->z = z;
}

bool Mesh::Point::operator<(const Point& rhs) const {
	const float tol = 1e-4f;

	if (x < rhs.x - tol) {
		return true;
	}
	if (std::abs(x - rhs.x) <= tol && y < rhs.y - tol) {
		return true;
	}
	if (std::abs(x - rhs.x) <= tol && std::abs(y - rhs.y) <= tol && z < rhs.z - tol) {
		return true;
	}
	return false;
}

Mesh::~Mesh() { this->end_SVG(); }

void Mesh::add_point(float x, float y, float z) {
	if (this->no_more_points) {
		throw std::runtime_error("Mesh::add_point: Can not add more points.");
	}

	this->points.push_back(Point(x, y, z));
}

int Mesh::find_point(float x, float y, float z) {
	if (!this->no_more_points) {
		this->finalize_points();
	}

	Point point(x, y, z);
	auto itr = std::lower_bound(points.begin(), points.end(), point);
	if (itr != points.end() && !(point < *itr)) {
		// Found the point; return its index.
		return int(itr - points.begin());
	} else {
		return -1;
	}
}

int Mesh::find_point(float x, float y, float z) const {
	if (!this->no_more_points) {
		throw std::runtime_error("find_point: mesh not finalized.");
	}

	Point point(x, y, z);
	auto itr = std::lower_bound(points.begin(), points.end(), point);
	if (itr != points.end() && !(point < *itr)) {
		// Found the point; return its index.
		return int(itr - points.begin());
	} else {
		return -1;
	}
}

int Mesh::number_of_points() const { return int(this->points.size()); }

const Mesh::Point& Mesh::get_point(int p) const { return this->points.at(p); }

void Mesh::finalize_points() {
	std::sort(points.begin(), points.end());
	this->no_more_points = true;
}

void Mesh::transform_points(float tx, float ty, float tz, float sx, float sy, float sz) {
	for (auto& point : points) {
		point.x = tx + sx * point.x;
		point.y = ty + sy * point.y;
		point.z = tz + sz * point.z;
	}
}

void Mesh::add_edge(float x1, float y1, float z1, float x2, float y2, float z2) {
	if (this->finished) {
		throw std::runtime_error("Mesh::add_edge: Can not add more edges.");
	}

	if (!this->no_more_points) {
		this->finalize_points();
	}

	int p1 = this->find_point(x1, y1, z1);
	int p2 = this->find_point(x2, y2, z2);

	if (p1 < 0 || p2 < 0) {
		throw std::runtime_error("Mesh:add_edge: point not found.");
	}
	if (p1 == p2) {
		throw std::runtime_error("Mesh:add_edge: p1 == p2.");
	}

	auto prev = points[p1].adjacent_points.size();

	auto insert = [](int i, std::vector<int>& v) {
		if (std::find(v.begin(), v.end(), i) == v.end()) {
			v.push_back(i);
		}
	};
	insert(p2, points[p1].adjacent_points);
	insert(p1, points[p2].adjacent_points);
	if (prev != points[p1].adjacent_points.size()) {
		// This edge was not previously present.
		edges.push_back(Edge(p1, p2));
		edges.push_back(Edge(p2, p1));
	}
}

void Mesh::add_edges(float dmax) {
	this->add_edges(dmax, [](float, float, float, float, float, float) { return false; });
}

struct UnitVector {
	float coordinates[3];

	UnitVector(float x, float y, float z) {
		float length = std::sqrt(x * x + y * y + z * z);
		coordinates[0] = x / length;
		coordinates[1] = y / length;
		coordinates[2] = z / length;
	}

	// Approximate comparison.
	bool operator<(const UnitVector& rhs) const {
		const float tol = 1e-4f;
		for (int i = 0; i < 3; ++i) {
			if (coordinates[i] < rhs.coordinates[i] - tol) {
				return true;
			}
			if (rhs.coordinates[i] < coordinates[i] - tol) {
				return false;
			}
		}
		return false;
	}
};

void Mesh::add_edges(float dmax,
                     const std::function<bool(float, float, float, float, float, float)>& ignore) {
	// Brute-force O(n^2) algorithm, should be
	// improved.
	for (int p1 = 0; p1 < points.size(); ++p1) {
		std::map<UnitVector, float> smallest_vector;
		std::map<UnitVector, std::vector<float>> edge_to_add;

		for (int p2 = p1 + 1; p2 < points.size(); ++p2) {
			float x1 = points[p1].x;
			float y1 = points[p1].y;
			float z1 = points[p1].z;

			float x2 = points[p2].x;
			float y2 = points[p2].y;
			float z2 = points[p2].z;

			float dx = x1 - x2;
			float dy = y1 - y2;
			float dz = z1 - z2;

			float d = dx * dx + dy * dy + dz * dz;
			if (d <= dmax * dmax) {
				UnitVector unit_vector(dx, dy, dz);
				if (!ignore(x1, y1, z1, x2, y2, z2)) {
					// Is this direction already added?
					auto itr = smallest_vector.find(unit_vector);
					if (itr == smallest_vector.end()) {
						smallest_vector[unit_vector] = d;
						auto& vec = edge_to_add[unit_vector];
						vec.resize(6);
						vec[0] = x1;
						vec[1] = y1;
						vec[2] = z1;
						vec[3] = x2;
						vec[4] = y2;
						vec[5] = z2;
					} else {
						if (d < itr->second) {
							itr->second = d;
							auto& vec = edge_to_add[unit_vector];
							vec.resize(6);
							vec[0] = x1;
							vec[1] = y1;
							vec[2] = z1;
							vec[3] = x2;
							vec[4] = y2;
							vec[5] = z2;
						}
					}
				}
			}
		}

		for (const auto& itr : edge_to_add) {
			const auto& vec = itr.second;
			this->add_edge(vec[0], vec[1], vec[2], vec[3], vec[4], vec[5]);
		}
	}
}

int Mesh::number_of_edges() const { return int(this->edges.size()); }

const Mesh::Edge& Mesh::get_edge(int e) const { return edges.at(e); }

int Mesh::find_edge(int p1, int p2) const {
	if (!this->finished) {
		throw std::runtime_error("find_edge: Need to finish mesh first.");
	}
	Edge edge(p1, p2);
	auto itr = std::lower_bound(edges.begin(), edges.end(), edge);
	if (itr == edges.end() || edge < *itr) {
		throw std::runtime_error("find_edge: edge not found.");
	}
	return int(itr - edges.begin());
}

void Mesh::get_adjacent_edges(int e, std::vector<int>* adjacent) const {
	adjacent->clear();
	int p1 = edges[e].first;
	int p2 = edges[e].second;

	// Find reverse edge.
	int e_rev = this->find_edge(p2, p1);

	auto& adjacency = points[p2].adjacent_points;
	for (auto itr = adjacency.begin(); itr != adjacency.end(); ++itr) {
		int e2 = this->find_edge(p2, *itr);
		if (e2 != e && e2 != e_rev) {
			adjacent->push_back(e2);
		}
	}
}

void Mesh::finish(bool create_pairs) {
	if (this->finished) {
		return;
	}

	if (!this->no_more_points) {
		this->finalize_points();
	}

	this->finished = true;

	// Compute the maximum node degree. This is fast compared
	// to the operations below.
	this->connectivity = 0;
	for (const auto& p : points) {
		this->connectivity = std::max(this->connectivity, int(p.adjacent_points.size()));
	}

	// Sort the array of edges. This defines which index each edge has.
	std::sort(edges.begin(), edges.end());

	if (create_pairs) {
		std::vector<int> adjacency;
		adjacency.reserve(100);

		// Count edge pairs.
		std::size_t num_edge_pairs = 0;
		for (int e = 0; e < edges.size(); ++e) {
			this->get_adjacent_edges(e, &adjacency);
			num_edge_pairs += adjacency.size();
		}

		this->edge_pairs.reserve(num_edge_pairs);

		// Create edge pairs.
		for (int e = 0; e < edges.size(); ++e) {
			int p1 = edges[e].first;
			int p2 = edges[e].second;

			this->get_adjacent_edges(e, &adjacency);
			for (auto e : adjacency) {
				int p2b = edges[e].first;
				int p3 = edges[e].second;

				// Debug check
				if (p2 != p2b) {
					throw std::runtime_error("Mesh::finish: internal mesh error");
				}

				this->edge_pairs.push_back(EdgePair(p1, p2, p3));
			}
		}

		// Sort the array of edge pairs. This defines which index each pair has.
		std::sort(edge_pairs.begin(), edge_pairs.end());
	}
}

int Mesh::get_connectivity() const {
	if (!this->finished) {
		throw std::runtime_error("get_connectivity: Need to finish mesh first.");
	}
	return this->connectivity;
}

int Mesh::find_edge_pair(int p1, int p2, int p3) const {
	if (!this->finished) {
		throw std::runtime_error("find_edge_pair: Need to finish mesh first.");
	}
	EdgePair edge_pair(p1, p2, p3);
	auto itr = std::lower_bound(edge_pairs.begin(), edge_pairs.end(), edge_pair);
	if (itr == edge_pairs.end() || edge_pair < *itr) {
		throw std::runtime_error("find_edge_pair: edge not found.");
	}
	return int(itr - edge_pairs.begin());
}

void Mesh::get_adjacent_pairs(int ep, std::vector<int>* adjacent) const {
	adjacent->clear();
	int p1, p2, p3;
	tie(p1, p2, p3) = edge_pairs[ep];

	int e2 = this->find_edge(p2, p3);

	// Static scratch space to avoid allocations.
	// NOT thread-safe.
	static std::vector<int> adjacent_edges;

	this->get_adjacent_edges(e2, &adjacent_edges);
	for (auto e : adjacent_edges) {
		int p3b = edges[e].first;
		int p4 = edges[e].second;
		// Debug check
		if (p3 != p3b) {
			throw std::runtime_error("Mesh::get_adjacent_pairs: internal mesh error");
		}

		int ep = find_edge_pair(p2, p3, p4);

		if (p4 != p1 && p4 != p2 && p4 != p3) {
			adjacent->push_back(ep);
		}
	}
}

int Mesh::number_of_edge_pairs() const { return int(this->edge_pairs.size()); }

const Mesh::EdgePair& Mesh::get_edge_pair(int ep) const { return this->edge_pairs.at(ep); }

float Mesh::edge_length(int e) const {
	const Edge& edge = this->edges.at(e);
	float dx = points[edge.first].x - points[edge.second].x;
	float dy = points[edge.first].y - points[edge.second].y;
	float dz = points[edge.first].z - points[edge.second].z;
	return std::sqrt(dx * dx + dy * dy + dz * dz);
}

void Mesh::start_SVG(const std::string& filename) {
	auto style = [](int) -> std::string { return "stroke-width:0.05;stroke:#4f4f4f;"; };
	start_SVG(filename, style);
}

void Mesh::start_SVG(const std::string& filename,
                     const std::function<std::string(int)>& edge_style) {
	if (fout) {
		end_SVG();
	}

	// Compute maximum dimensions.
	float max_x = 0;
	float max_y = 0;
	for (const auto& point : points) {
		max_x = std::max(max_x, point.x);
		max_y = std::max(max_y, point.y);
	}

	// Output file.
	this->fout = new std::ofstream(filename.c_str());

	// Draw SVG header.
	*fout << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
	*fout << "<svg width=\"" << max_x << "\" height=\"" << max_y
	      << "\" id=\"svg2\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n";

	// Draw background mesh.
	//<line style="stroke-width:0.05;stroke:#4f4f4f;"  x1="0" y1="48.15" x2="12.025" y2="56.175" />
	for (int e = 0; e < edges.size(); ++e) {
		const auto& edge = edges[e];
		// There are zero or two edges for each pair of points.
		// Only one needs to be drawn.
		if (edge.first < edge.second) {
			std::string style = edge_style(e);
			if (style.length() > 0) {
				*fout << "<line style=\"" << style << "\" ";
				*fout << "x1=\"" << points[edge.first].x << "\" ";
				*fout << "y1=\"" << points[edge.first].y << "\" ";
				*fout << "x2=\"" << points[edge.second].x << "\" ";
				*fout << "y2=\"" << points[edge.second].y << "\" ";
				*fout << "/>\n";
			}
		}
	}
}

auto Mesh::edgepath_to_points(const std::vector<int>& path) const -> std::vector<Point> {
	std::vector<Mesh::Point> point_vector;

	// Start point
	if (path.size() > 0) {
		const auto& edge = edges[path[0]];
		point_vector.push_back(points[edge.first]);
	}

	for (auto edge_index : path) {
		const auto& edge = edges[edge_index];
		point_vector.push_back(points[edge.second]);
	}

	return point_vector;
}

auto Mesh::pairpath_to_points(const std::vector<int>& path) const -> std::vector<Point> {
	std::vector<Mesh::Point> point_vector;

	// Start points
	if (path.size() > 0) {
		int p1, p2;
		tie(p1, p2, ignore) = this->edge_pairs.at(path[0]);
		point_vector.push_back(points[p1]);
		point_vector.push_back(points[p2]);
	}

	for (auto pair_index : path) {
		int p3;
		tie(ignore, ignore, p3) = this->edge_pairs.at(pair_index);
		point_vector.push_back(points[p3]);
	}

	return point_vector;
}

void Mesh::draw_path(const std::vector<int>& path, const std::string& color) {
	if (!fout) {
		throw std::runtime_error("draw_path: call start_SVG first.");
	}

	*fout << "<g>\n";
	for (auto e : path) {
		const auto& edge = edges.at(e);
		*fout << "<line style=\"stroke-width:0.2;stroke:#" << color << ";\" ";
		*fout << "x1=\"" << points[edge.first].x << "\" ";
		*fout << "y1=\"" << points[edge.first].y << "\" ";
		*fout << "x2=\"" << points[edge.second].x << "\" ";
		*fout << "y2=\"" << points[edge.second].y << "\" ";
		*fout << "/>\n";
	}
	*fout << "</g>\n";
}

void Mesh::end_SVG() {
	if (fout) {
		*fout << "</svg>";
		delete fout;
		fout = 0;
	}
}
}  // namespace curvature
}  // namespace minimum