// Petter Strandmark 2013.
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

using std::max;
using std::min;

#include <minimum/curvature/grid_mesh.h>

namespace minimum {
namespace curvature {

struct UnitVector {
	float coordinates[3];

	UnitVector(float x, float y, float z) {
		float length = std::sqrt(x * x + y * y + z * z);
		if (length < 1e-4) {
			throw std::runtime_error("UnitVector: too short.");
		}
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

void GridMesh::initialize(
    int x_dim,
    int y_dim,
    int z_dim,
    double dmax,
    const std::function<bool(float, float, float, float, float, float)>* ignore,
    bool create_pairs) {
	for (int x = 0; x < x_dim; ++x) {
		for (int y = 0; y < y_dim; ++y) {
			for (int z = 0; z < z_dim; ++z) {
				this->add_point(x, y, z);
			}
		}
	}

	for (int x1 = 0; x1 < x_dim; ++x1) {
		for (int y1 = 0; y1 < y_dim; ++y1) {
			for (int z1 = 0; z1 < z_dim; ++z1) {
				std::map<UnitVector, float> smallest_vector;
				std::map<UnitVector, std::vector<float>> edge_to_add;

				for (int x2 = x1; x2 <= min(x_dim - 1.0, x1 + dmax); ++x2) {
					for (int y2 = max(0.0, y1 - dmax); y2 <= min(y_dim - 1.0, y1 + dmax); ++y2) {
						for (int z2 = max(0.0, z1 - dmax); z2 <= min(z_dim - 1.0, z1 + dmax);
						     ++z2) {
							double dx = x1 - x2;
							double dy = y1 - y2;
							double dz = z1 - z2;
							double d = dx * dx + dy * dy + dz * dz;
							if (d > 0.1 && d <= dmax * dmax + 1e-6) {
								if (!ignore || !(*ignore)(x1, y1, z1, x2, y2, z2)) {
									UnitVector unit_vector(dx, dy, dz);

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
					}
				}

				for (auto itr = edge_to_add.begin(); itr != edge_to_add.end(); ++itr) {
					const auto& vec = itr->second;
					this->add_edge(vec[0], vec[1], vec[2], vec[3], vec[4], vec[5]);
				}
			}
		}
	}

	this->finish(create_pairs);
}
}  // namespace curvature
}  // namespace minimum
