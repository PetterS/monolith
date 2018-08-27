#ifndef MINIMUM_NONLINEAR_INTERVAL_TERM
#define MINIMUM_NONLINEAR_INTERVAL_TERM

#include <memory>

#include <minimum/nonlinear/term.h>

namespace minimum {
namespace nonlinear {

template <int... D>
class IntervalTermBase : public SizedTerm<D...> {
   public:
	virtual double evaluate(double* const* const variables,
	                        std::vector<Eigen::VectorXd>* gradient) const override {
		core::check(false, "evaluate with gradient: not implemented.");
		return 0;
	}

	virtual double evaluate(double* const* const variables,
	                        std::vector<Eigen::VectorXd>* gradient,
	                        std::vector<std::vector<Eigen::MatrixXd>>* hessian) const override {
		core::check(false, "evaluate with hessian: not implemented.");
		return 0;
	}
};

template <typename Functor, int... D>
class IntervalTerm : public IntervalTermBase<D...> {};

// Works the same way as make_differentiable.
//
// See global_test.cpp for examples.
template <int... arg_sizes, typename Functor>
std::shared_ptr<Term> make_interval_term(Functor&& lambda) {
	static_assert(sizeof...(arg_sizes) > 0,
	              "make_interval_term<n0, n1, ...> must at least provide n0.");
	typedef typename std::remove_reference<Functor>::type FunctorClass;
	typedef IntervalTerm<FunctorClass, arg_sizes...> TermType;
	return std::make_shared<TermType>(std::forward<Functor>(lambda));
}

//
// 1-variable specialization
//
template <typename Functor, int D0>
class IntervalTerm<Functor, D0> : public IntervalTermBase<D0> {
   public:
	template <typename... Args>
	IntervalTerm(Args&&... args) : functor(std::forward<Args>(args)...) {}

	virtual double evaluate(double* const* const variables) const override {
		return functor(variables[0]);
	}

	virtual Interval<double> evaluate_interval(
	    const Interval<double>* const* const variables) const override {
		return functor(variables[0]);
	};

   private:
	Functor functor;
};

//
// 2-variable specialization
//
template <typename Functor, int D0, int D1>
class IntervalTerm<Functor, D0, D1> : public IntervalTermBase<D0, D1> {
   public:
	template <typename... Args>
	IntervalTerm(Args&&... args) : functor(std::forward<Args>(args)...) {}

	virtual double evaluate(double* const* const variables) const override {
		return functor(variables[0], variables[1]);
	}

	virtual Interval<double> evaluate_interval(
	    const Interval<double>* const* const variables) const override {
		return functor(variables[0], variables[1]);
	};

   private:
	Functor functor;
};
}  // namespace nonlinear
}  // namespace minimum
#endif
