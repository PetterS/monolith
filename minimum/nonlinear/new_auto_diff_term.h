#ifndef MINIMUM_NONLINEAR_NEW_AUTO_DIFF_TERM_H
#define MINIMUM_NONLINEAR_NEW_AUTO_DIFF_TERM_H

#include <type_traits>
#include <typeinfo>

#include <badiff.h>
#include <fadiff.h>

#include <minimum/nonlinear/term.h>

namespace minimum {
namespace nonlinear {

template <typename Functor, typename R, int... D>
struct FunctorCaller;

template <typename Functor, typename R, int D0, int... DN>
struct FunctorCaller<Functor, R, D0, DN...> {
	template <typename... T>
	R call(const Functor& functor, R* const* const variables, T... rest) {
		R x[D0];
		for (int i = 0; i < D0; ++i) {
			x[i] = (*variables)[i];
		}

		FunctorCaller<Functor, R, DN...> next_caller;
		return next_caller.call(functor, variables + 1, rest..., x);
	}
};

template <typename Functor, typename R>
struct FunctorCaller<Functor, R> {
	template <typename... T>
	R call(const Functor& functor, R* const* const variables, T... rest) {
		minimum_core_assert(false, "Not implemented!");
		return functor(rest...);
	}
};

template <int... D>
struct Product;

template <int D0, int... DN>
struct Product<D0, DN...> {
	static const int value = D0 * Product<DN...>::value;
};

template <>
struct Product<> {
	static const int value = 1;
};

static_assert(Product<5>::value == 5, "Product test failed.");
static_assert(Product<5, 2>::value == 5 * 2, "Product test failed.");
static_assert(Product<5, 2, 3>::value == 5 * 2 * 3, "Product test failed.");
static_assert(Product<5, 2, 3, 5>::value == 5 * 2 * 3 * 5, "Product test failed.");

// Double specialization without allocation of temporary
// variables on the stack.
template <typename Functor, int D0, int... DN>
struct FunctorCaller<Functor, double, D0, DN...> {
	template <typename... T>
	double call(const Functor& functor, double* const* const variables, T... rest) {
		FunctorCaller<Functor, double, DN...> next_caller;
		return next_caller.call(functor, variables + 1, rest..., variables[0]);
	}
};

template <typename Functor>
struct FunctorCaller<Functor, double> {
	template <typename... T>
	double call(const Functor& functor, double* const* const variables, T... rest) {
		return functor(rest...);
	}
};

//
// Term which allows for automatic computation of derivatives. It is
// used in the following way:
//
//   auto term = make_shared<AutoDiffTerm<Functor, 1>>(arg1, arg2, ...)
//
// where arg1, arg2, etc. are arguments to the constructor of Functor.
//
template <typename Functor, int... D>
class NewAutoDiffTerm : public SizedTerm<D...> {
   public:
	template <typename... Args>
	NewAutoDiffTerm(Args&&... args) : functor(std::forward<Args>(args)...) {}

	virtual double evaluate(double* const* const variables) const override {
		FunctorCaller<Functor, double, D...> caller;
		return caller.call(functor, variables);
	}

	virtual double evaluate(double* const* const variables,
	                        std::vector<Eigen::VectorXd>* gradient) const override {
		check(false, "gradient not implemented.");
		return 0;
	}

	virtual double evaluate(double* const* const variables,
	                        std::vector<Eigen::VectorXd>* gradient,
	                        std::vector<std::vector<Eigen::MatrixXd>>* hessian) const override {
		check(false, "hessian not implemented.");
		return 0;
	}

   protected:
	Functor functor;
};
}  // namespace nonlinear
}  // namespace minimum
.

#endif
