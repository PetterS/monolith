// Petter Strandmark.
#ifndef CURVE_EXTRACTION_CURVATURE_H
#define CURVE_EXTRACTION_CURVATURE_H

#include <minimum/curvature/export.h>

namespace minimum {
namespace curvature {

RESEARCH_CURVATURE_API extern int curvature_cache_hits;
RESEARCH_CURVATURE_API extern int curvature_cache_misses;
template <typename R>
R compute_curvature(R x1,
                    R y1,
                    R z1,
                    R x2,
                    R y2,
                    R z2,
                    R x3,
                    R y3,
                    R z3,
                    R power = 2.0,
                    bool writable_cache = true,
                    int n_approximation_points = 200);

RESEARCH_CURVATURE_API extern int torsion_cache_hits;
RESEARCH_CURVATURE_API extern int torsion_cache_misses;

template <typename R>
R compute_torsion(R x1,
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
                  R power = 2.0,
                  bool writable_cache = true,
                  int n_approximation_points = 200);
}  // namespace curvature
}  // namespace minimum

#ifdef WIN32
#define INSTANTIATE_CURVATURE(R)                                             \
	template R RESEARCH_CURVATURE_API minimum::curvature::compute_curvature( \
	    R, R, R, R, R, R, R, R, R, R, bool, int);
INSTANTIATE_CURVATURE(double);
INSTANTIATE_CURVATURE(float);

#define INSTANTIATE_TORSION(R)                                             \
	template R RESEARCH_CURVATURE_API minimum::curvature::compute_torsion( \
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

#endif
