// Petter Strandmark.

#include <algorithm>

#include <minimum/nonlinear/auto_diff_term.h>
using minimum::nonlinear::to_double;

#include <minimum/curvature/data_term.h>

namespace minimum {
namespace curvature {

TriLinear::TriLinear(
    const double* unary_, int M_, int N_, int O_, const std::vector<double>& voxeldimensions_)
    : M(M_), N(N_), O(O_), unary(unary_), voxeldimensions(voxeldimensions_) {}

//
//  pixel ix = 5 has real coordinates [4.5, 5.5).
//
int TriLinear::xyz_to_ind(double x, double y, double z) const {
	int ix = std::max(std::min(int(x + 0.5), M - 1), 0);
	int iy = std::max(std::min(int(y + 0.5), N - 1), 0);
	int iz = std::max(std::min(int(z + 0.5), O - 1), 0);
	return ix + M * iy + M * N * iz;
}

template <typename R>
R TriLinear::evaluate(R x, R y, R z) const {
	return 0;
}

template <typename R>
R TriLinear::evaluate_line_integral(R sx, R sy, R sz, R ex, R ey, R ez) const {
	using std::floor;
	using std::sqrt;

	int interpolation_points = 25;

	R dx = ((ex - sx) * voxeldimensions[0]);
	R dy = ((ey - sy) * voxeldimensions[1]);
	R dz = ((ez - sz) * voxeldimensions[2]);

	R segment_length = sqrt(dx * dx + dy * dy + dz * dz) / interpolation_points;

	dx /= (interpolation_points - 1);
	dy /= (interpolation_points - 1);
	dz /= (interpolation_points - 1);

	R cost = 0;

	// Goes along the line connecting the start and end point
	// in each point we perform trilinear interpolation
	for (int i = 0; i < interpolation_points; i++) {
		// Current location along the line connecting start and stop
		// +0.5 since we define (0,0,0) to be center of voxel 0 0 0.
		R cx = sx + dx * R(i) + R(0.5);
		R cy = sy + dy * R(i) + R(0.5);
		R cz = sz + dz * R(i) + R(0.5);

		if (dx == 0) cx = 0;

		if (dy == 0) cy = 0;

		if (dz == 0) cz = 0;

		// What voxel are we currently inside
		int vx = floor(to_double(cx));
		int vy = floor(to_double(cy));
		int vz = floor(to_double(cz));

		// Relative to center of voxel how far away is
		// each dimension?
		R rx = cx - R(vx);
		R ry = cy - R(vy);
		R rz = cz - R(vz);

		// Set values
		auto get_intensity = [&](int x, int y, int z) -> double {
			return unary[xyz_to_ind(x, y, z)];
		};

		// Intensity at neighboring voxels
		double V000 = get_intensity(vx + 0, vy + 0, vz + 0);
		double V001 = get_intensity(vx + 0, vy + 0, vz + 1);
		double V010 = get_intensity(vx + 0, vy + 1, vz + 0);
		double V011 = get_intensity(vx + 0, vy + 1, vz + 1);
		double V100 = get_intensity(vx + 1, vy + 0, vz + 0);
		double V101 = get_intensity(vx + 1, vy + 0, vz + 1);
		double V110 = get_intensity(vx + 1, vy + 1, vz + 0);
		double V111 = get_intensity(vx + 1, vy + 1, vz + 1);

		// Trilinear interpolation
		cost += V000 * (1 - rx) * (1 - ry) * (1 - rz) + V100 * (0 + rx) * (1 - ry) * (1 - rz)
		        + V010 * (1 - rx) * (0 + ry) * (1 - rz) + V001 * (1 - rx) * (1 - ry) * (0 + rz)
		        + V101 * (0 + rx) * (1 - ry) * (0 + rz) + V011 * (1 - rx) * (0 + ry) * (0 + rz)
		        + V110 * (0 + rx) * (0 + ry) * (1 - rz) + V111 * (0 + rx) * (0 + ry) * (1 - rz);
	}

	cost *= segment_length;
	return cost;
}
}  // namespace curvature
}  // namespace minimum

#define INSTANTIATE(classname, R)                                                             \
	template R minimum::curvature::classname::evaluate_line_integral(R, R, R, R, R, R) const; \
	template R minimum::curvature::classname::evaluate(R, R, R) const;

typedef fadbad::F<double, 3> F3;
typedef fadbad::F<double, 6> F6;
typedef fadbad::F<fadbad::F<double, 3>, 3> FF3;
typedef fadbad::F<fadbad::F<double, 6>, 6> FF6;
typedef fadbad::F<fadbad::F<double, 9>, 9> FF9;
INSTANTIATE(TriLinear, double)
INSTANTIATE(TriLinear, F3)
INSTANTIATE(TriLinear, F6)
INSTANTIATE(TriLinear, FF3)
INSTANTIATE(TriLinear, FF6)
INSTANTIATE(TriLinear, FF9)
#undef INSTANTIATE
