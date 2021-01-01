// Petter Strandmark.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>

// to_double
#include <minimum/nonlinear/auto_diff_term.h>

#include <minimum/curvature/curvature.h>

namespace fadbad {

template <unsigned n>
F<double, n> abs(F<double, n> x) {
	if (minimum::nonlinear::to_double(x) < 0) {
		return -x;
	} else {
		return x;
	}
}

template <unsigned n>
F<F<double, n>, n> abs(F<F<double, n>, n> x) {
	if (minimum::nonlinear::to_double(x) < 0) {
		return -x;
	} else {
		return x;
	}
}
}  // namespace fadbad

// This struct is a tuple of floating point values.
// Two tuples are considered equal if no two elements
// differ by more than 1e-4.
template <typename R, int n>
struct FloatingPointCacheEntry {
	R data[n];

	bool operator<(const FloatingPointCacheEntry& rhs) const {
		const R tol = R(1e-4);
		for (int i = 0; i < n; ++i) {
			if (data[i] < rhs.data[i] - tol) {
				return true;
			}
			if (rhs.data[i] < data[i] - tol) {
				return false;
			}
		}
		return false;
	}
};

template <typename R>
struct use_cache {
	static const bool value = false;
};

template <>
struct use_cache<double> {
	static const bool value = true;
};

template <>
struct use_cache<float> {
	static const bool value = true;
};

template <typename R>
R compute_curvature_internal(R x1, R y1, R z1, R x2, R y2, R z2, R x3, R y3, R z3, R p, int n) {
	using std::abs;
	using std::atan2;
	using std::pow;
	using std::sqrt;

	const R a = x1 - 2 * x2 + x3;
	const R b = y1 - 2 * y2 + y3;
	const R c = z1 - 2 * z2 + z3;
	const R dx = x1 - x2;
	const R dy = y1 - y2;
	const R dz = z1 - z2;

	const R sq1 = c * dy - b * dz;
	const R sq2 = c * dx - a * dz;
	const R sq3 = b * dx - a * dy;

	const R numerator = pow(sq1 * sq1 + sq2 * sq2 + sq3 * sq3, p / R(2.0));

	// k(t)^p * ds
	auto kp_ds = [a, b, c, dx, dy, dz, numerator, p](R t) -> R {
		const R sq4 = a * t - dx;
		const R sq5 = b * t - dy;
		const R sq6 = c * t - dz;
		const R exponent = 3 * p / R(2.0) - R(0.5);
		return numerator / pow(sq4 * sq4 + sq5 * sq5 + sq6 * sq6, exponent);
	};

	// Trapezoidal rule.
	R sum = 0;
	sum += kp_ds(0) / R(2.0);
	sum += kp_ds(1) / R(2.0);
	for (int k = 1; k <= n - 1; ++k) {
		sum += kp_ds(R(k) / R(n));
	}
	sum /= R(n);

	return sum;
}

int minimum::curvature::curvature_cache_hits = 0;
int minimum::curvature::curvature_cache_misses = 0;

template <typename R>
R minimum::curvature::compute_curvature(
    R x1, R y1, R z1, R x2, R y2, R z2, R x3, R y3, R z3, R p, bool writable_cache, int n) {
	static std::map<FloatingPointCacheEntry<R, 8>, R> curvature_cache;
	const std::size_t max_cache_size = 1000000;
	FloatingPointCacheEntry<R, 8> entry;

	if (use_cache<R>::value) {
		// The static cache of curvature values allows
		// for much faster computations if the coordinates
		// come from a regular grid.
		// Create the entry by subtracting the mean from the coordinates.
		entry.data[0] = x2 - x1;
		entry.data[1] = x3 - x1;
		entry.data[2] = y2 - y1;
		entry.data[3] = y3 - y1;
		entry.data[4] = z2 - z1;
		entry.data[5] = z3 - z1;
		// p and n should also be in the cache.
		entry.data[6] = p;
		entry.data[7] = R(n);
		// If the value is in the cache, return it.
		auto itr = curvature_cache.find(entry);
		if (itr != curvature_cache.end()) {
			curvature_cache_hits++;
			return itr->second;
		} else {
			curvature_cache_misses++;
		}
	}

	R value = compute_curvature_internal(x1, y1, z1, x2, y2, z2, x3, y3, z3, p, n);
	if (use_cache<R>::value) {
		if (writable_cache && curvature_cache.size() < max_cache_size) {
			curvature_cache[entry] = value;
		}
	}

	return value;
}

int minimum::curvature::torsion_cache_hits = 0;
int minimum::curvature::torsion_cache_misses = 0;

template <typename R>
R minimum::curvature::compute_torsion(R x1,
                                      R y1,
                                      R z1,
                                      R x2,
                                      R y2,
                                      R z2,
                                      R x3,
                                      R y3,
                                      R z3,
                                      R x4,
                                      R y4,
                                      R z4,
                                      R p,
                                      bool writable_cache,
                                      int n) {
	using std::abs;
	using std::atan2;
	using std::pow;
	using std::sqrt;

	// The static cache of torsion values allows
	// for much faster computations if the coordinates
	// come from a regular grid.
	static std::map<FloatingPointCacheEntry<R, 14>, R> torsion_cache;
	const std::size_t max_cache_size = 10000000;
	FloatingPointCacheEntry<R, 14> entry;

	if (use_cache<R>::value) {
		// Create the entry by subtracting the mean from the coordinates.
		R mx = (x1 + x2 + x3 + x4) / R(4.0);
		R my = (y1 + y2 + y3 + y4) / R(4.0);
		R mz = (z1 + z2 + z3 + z4) / R(4.0);
		entry.data[0] = x1 - mx;
		entry.data[1] = x2 - mx;
		entry.data[2] = x3 - mx;
		entry.data[3] = x4 - mx;
		entry.data[4] = y1 - my;
		entry.data[5] = y2 - my;
		entry.data[6] = y3 - my;
		entry.data[7] = y4 - my;
		entry.data[8] = z1 - mz;
		entry.data[9] = z2 - mz;
		entry.data[10] = z3 - mz;
		entry.data[11] = z4 - mz;
		// p and n should also be in the cache.
		entry.data[12] = p;
		entry.data[13] = R(n);
		// If the value is in the cache, return it.
		auto itr = torsion_cache.find(entry);
		if (itr != torsion_cache.end()) {
			torsion_cache_hits++;
			return itr->second;
		} else {
			torsion_cache_misses++;
		}
	}

	R numerator =
	    (y3 * (x2 * z4 - x2 * z1 + x4 * (z1 - z2)) - y4 * (x2 * z3 - x2 * z1 + x1 * (z2 - z3))
	     - x3 * (y1 * (z2 - z4) + y4 * (z1 - z2))
	     - y2 * (x1 * (z3 - z4) - x3 * (z1 - z4) + x4 * (z1 - z3)) + y1 * (x2 * z3 - x2 * z4)
	     + x1 * y3 * (z2 - z4) + x4 * y1 * (z2 - z3));

	R sum = 0;
	if (std::abs(minimum::nonlinear::to_double(numerator)) > 1e-9) {
		auto tp_ds = [=](R t) -> R {
			R s1 = ((z1 - 3 * z2 + 3 * z3 - z4) / R(2)) * t * t + (2 * z2 - z1 - z3) * t
			       + (z1 - z3) / R(2);
			R s2 = ((y1 - 3 * y2 + 3 * y3 - y4) / R(2)) * t * t + (2 * y2 - y1 - y3) * t
			       + (y1 - y3) / R(2);
			R s3 = ((x1 - 3 * x2 + 3 * x3 - x4) / R(2)) * t * t + (2 * x2 - x1 - x3) * t
			       + (x1 - x3) / R(2);
			R s4 = z1 - 2 * z2 + z3 - t * (z1 - 3 * z2 + 3 * z3 - z4);
			R s5 = y1 - 2 * y2 + y3 - t * (y1 - 3 * y2 + 3 * y3 - y4);
			R s6 = x1 - 2 * x2 + x3 - t * (x1 - 3 * x2 + 3 * x3 - x4);
			R sq1 = (s5 * s3 - s2 * s6);
			R sq2 = (s4 * s3 - s1 * s6);
			R sq3 = (s2 * s4 - s1 * s5);
			R denominator = sq1 * sq1 + sq2 * sq2 + sq3 * sq3;

			return sqrt(s3 * s3 + s2 * s2 + s1 * s1)        // ds
			       * pow(abs(numerator / denominator), p);  // |tau|^p
		};

		// Trapezoidal rule.
		sum = 0;
		sum += tp_ds(0) / R(2.0);
		sum += tp_ds(1) / R(2.0);
		for (int k = 1; k <= n - 1; ++k) {
			sum += tp_ds(R(k) / R(n));
		}
		sum /= R(n);

		// Debug: Check for nan.
		if (minimum::nonlinear::to_double(sum) != minimum::nonlinear::to_double(sum)) {
			throw std::runtime_error("compute_torsion: nan");
		}
	}

	if (use_cache<R>::value) {
		// Set the cache and return.
		if (writable_cache && torsion_cache.size() < max_cache_size) {
			torsion_cache[entry] = sum;
		}
	}

	return sum;
}

#ifndef WIN32
#define INSTANTIATE_CURVATURE(R) \
	template R minimum::curvature::compute_curvature(R, R, R, R, R, R, R, R, R, R, bool, int);
INSTANTIATE_CURVATURE(double);
INSTANTIATE_CURVATURE(float);

#define INSTANTIATE_TORSION(R)                      \
	template R minimum::curvature::compute_torsion( \
	    R, R, R, R, R, R, R, R, R, R, R, R, R, bool, int);
INSTANTIATE_TORSION(double);
INSTANTIATE_TORSION(float);

// Create instantiations for auto-diff.
#define INSTANTIATE_FN(n)                             \
	typedef fadbad::F<double, n> F##n;                \
	INSTANTIATE_CURVATURE(F##n);                      \
	INSTANTIATE_TORSION(F##n);                        \
	typedef fadbad::F<fadbad::F<double, n>, n> FF##n; \
	INSTANTIATE_CURVATURE(FF##n);                     \
	INSTANTIATE_TORSION(FF##n);
INSTANTIATE_FN(1);
INSTANTIATE_FN(2);
INSTANTIATE_FN(3);
INSTANTIATE_FN(6);
INSTANTIATE_FN(9);
INSTANTIATE_FN(12);
#endif
