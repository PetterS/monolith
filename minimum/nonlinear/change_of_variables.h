// Petter Strandmark.
#ifndef MINIMUM_NONLINEAR_CHANGE_OF_VARIABLES_H
#define MINIMUM_NONLINEAR_CHANGE_OF_VARIABLES_H

#include <cstring>
using std::size_t;

namespace minimum {
namespace nonlinear {

// This class allows for a change of variables to be added.
// For example,
//
//    minimize f(x1, x2)
//    subject to x1 >= 0
//
// can be formulated as
//
//    minimize f(exp(t), x2)
//
class ChangeOfVariables {
   public:
	virtual ~ChangeOfVariables(){};
	virtual void t_to_x(double* x, const double* t) const = 0;
	virtual void x_to_t(double* t, const double* xt) const = 0;
	virtual int x_dimension() const = 0;
	virtual int t_dimension() const = 0;
	virtual void update_gradient(double* t_gradient,
	                             const double* t_input,
	                             const double* x_gradient) const = 0;
	void update_hessian() const;
};
}  // namespace nonlinear
}  // namespace minimum

#endif
