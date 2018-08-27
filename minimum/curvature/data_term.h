// Petter Strandmark 2013.
#ifndef CURVE_EXTRACTION_DATA_TERM_H
#define CURVE_EXTRACTION_DATA_TERM_H

#include <vector>

#include <minimum/curvature/export.h>

namespace minimum {
namespace curvature {

class RESEARCH_CURVATURE_API PieceWiseConstant {
   public:
	PieceWiseConstant(
	    const double* unary, int M, int N, int O, const std::vector<double>& voxeldimensions);

	template <typename R>
	R evaluate(R x, R y, R z = 0.0) const;

	template <typename R>
	R evaluate_line_integral(R x1, R y1, R z1, R x2, R y2, R z2) const;

   private:
	int xyz_to_ind(double x, double y, double z) const;
	template <typename R>
	bool inside_volume(R x, R y, R s) const;
	int M, N, O;
	const double* unary;
	const std::vector<double> voxeldimensions;
};

class RESEARCH_CURVATURE_API TriLinear {
   public:
	TriLinear(const double* unary, int M, int N, int O, const std::vector<double>& voxeldimensions);

	template <typename R>
	R evaluate(R x, R y, R z = 0.0) const;

	template <typename R>
	R evaluate_line_integral(R x1, R y1, R z1, R x2, R y2, R z2) const;

   private:
	int xyz_to_ind(double x, double y, double z) const;
	int M, N, O;
	const double* unary;
	const std::vector<double> voxeldimensions;
};

#ifdef WIN32
#define INSTANTIATE(classname, R)                                                            \
	template R RESEARCH_CURVATURE_API minimum::curvature::classname::evaluate_line_integral( \
	    R, R, R, R, R, R) const;                                                             \
	template R RESEARCH_CURVATURE_API minimum::curvature::classname::evaluate(R, R, R) const;

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
}  // namespace curvature
}  // namespace minimum

#endif
