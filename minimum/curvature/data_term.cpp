// Petter Strandmark 2013.
#include <algorithm>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <minimum/nonlinear/auto_diff_term.h>
using minimum::nonlinear::to_double;

#include <minimum/curvature/data_term.h>

namespace fadbad {

template <unsigned n>
F<double, n> abs(F<double, n> x) {
	if (to_double(x) < 0) {
		return -x;
	} else {
		return x;
	}
}

template <unsigned n>
F<F<double, n>, n> abs(F<F<double, n>, n> x) {
	if (to_double(x) < 0) {
		return -x;
	} else {
		return x;
	}
}
}  // namespace fadbad

namespace minimum {
namespace curvature {

PieceWiseConstant::PieceWiseConstant(
    const double* unary_, int M_, int N_, int O_, const std::vector<double>& voxeldimensions_)

    : M(M_), N(N_), O(O_), unary(unary_), voxeldimensions(voxeldimensions_) {}

//
//  pixel ix = 5 has real coordinates [4.5, 5.5).
//
int PieceWiseConstant::xyz_to_ind(double x, double y, double z) const {
	int ix = std::max(std::min(int(x + 0.5), M - 1), 0);
	int iy = std::max(std::min(int(y + 0.5), N - 1), 0);
	int iz = std::max(std::min(int(z + 0.5), O - 1), 0);
	return ix + M * iy + M * N * iz;
}

template <typename R>
R PieceWiseConstant::evaluate(R x, R y, R z) const {
	return unary[xyz_to_ind(to_double(x), to_double(y), to_double(z))];
}

template <typename R>
R fractional_part(R x) {
	double integer_part;
	modf(minimum::nonlinear::to_double(x), &integer_part);
	return x - integer_part;
}

template <typename R>
bool PieceWiseConstant::inside_volume(R x, R y, R z) const {
	if (to_double(x) < -0.5) return false;
	if (to_double(y) < -0.5) return false;
	if (to_double(z) < -0.5) return false;

	if (to_double(x) > M - 0.5) return false;
	if (to_double(y) > N - 0.5) return false;
	if (to_double(z) > O - 0.5) return false;

	return true;
}

template <typename R>
R PieceWiseConstant::evaluate_line_integral(R sx, R sy, R sz, R ex, R ey, R ez) const {
	using std::abs;
	using std::sqrt;

	if (!inside_volume(sx, sy, sz)) return R(std::numeric_limits<double>::infinity());

	if (!inside_volume(ex, ey, ez)) return R(std::numeric_limits<double>::infinity());

	R dx = (ex - sx) * voxeldimensions[0];
	R dy = (ey - sy) * voxeldimensions[1];
	R dz = (ez - sz) * voxeldimensions[2];

	R line_length = sqrt(dx * dx + dy * dy + dz * dz);

	typedef std::pair<R, int> crossing;

	// Integral invariant of direction
	auto set_crossings =
	    [](int index_change, R k, R dx, R start, std::vector<crossing>& crossings) -> void {
		if (abs(to_double(k)) <= 1e-4f) return;

		auto shifted_start = start + 0.5;
		R current = 1.0 - fractional_part(shifted_start);

		if (k < 0) {
			current = 1.0 - current;
			index_change = -index_change;
		}

		for (; to_double(current) <= abs(to_double(dx)); current = current + 1.0) {
			crossings.push_back(crossing(abs(current / k), index_change));
		}
	};

	// Have a static vector for temporary storage.
	// Avoids memory allocations for almost all calls.
	std::vector<crossing>* scratch_space;
#ifdef USE_OPENMP
	// Maximum 100 threads.
	// Support for C++11 thread_local is not very good yet.
	static std::vector<crossing> local_space[100];
	scratch_space = &local_space[omp_get_thread_num()];
#else
	static std::vector<crossing> local_space;
	scratch_space = &local_space;
#endif

	scratch_space->clear();
	scratch_space->push_back(crossing(0.0f, 0));  // start

	// x
	R kx = dx / line_length;
	set_crossings(1, kx, dx, sx, *scratch_space);

	// y
	R ky = dy / line_length;
	set_crossings(M, ky, dy, sy, *scratch_space);

	// z
	R kz = dz / line_length;
	set_crossings(M * N, kz, dz, sz, *scratch_space);

	// last stretch
	scratch_space->push_back(crossing(line_length, 0));

	// Sort vector.
	std::sort(scratch_space->begin(), scratch_space->end());
	auto crossings_end = scratch_space->end();

	auto prev = scratch_space->begin();
	auto next = scratch_space->begin();

	int source_id = xyz_to_ind(to_double(sx), to_double(sy), to_double(sz));

	R cost = 0;
	for (next++; next != crossings_end; prev++, next++) {
		R distance = (next->first - prev->first);

		// Removes small distances without this small numerical
		// error might lead the data cost to sample in the wrong region.
		// If the cost is \inf this is a problem no matter how tiny the distance.
		if (distance > 1e-6) cost += distance * unary[source_id];

		// Which dimension did we cross?
		source_id += next->second;
	}

	return cost;
}
}  // namespace curvature
}  // namespace minimum

#ifndef WIN32
#define INSTANTIATE(classname, R)                                                             \
	template R minimum::curvature::classname::evaluate_line_integral(R, R, R, R, R, R) const; \
	template R minimum::curvature::classname::evaluate(R, R, R) const;

typedef fadbad::F<double, 3> F3;
typedef fadbad::F<double, 6> F6;
typedef fadbad::F<fadbad::F<double, 3>, 3> FF3;
typedef fadbad::F<fadbad::F<double, 6>, 6> FF6;
typedef fadbad::F<fadbad::F<double, 9>, 9> FF9;
INSTANTIATE(PieceWiseConstant, double)
INSTANTIATE(PieceWiseConstant, F3)
INSTANTIATE(PieceWiseConstant, F6)
INSTANTIATE(PieceWiseConstant, FF3)
INSTANTIATE(PieceWiseConstant, FF6)
INSTANTIATE(PieceWiseConstant, FF9)
#endif
