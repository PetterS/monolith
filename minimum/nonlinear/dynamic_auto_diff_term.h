#ifndef MINIMUM_NONLINEAR_LARGE_AUTO_DIFF_TERM_H
#define MINIMUM_NONLINEAR_LARGE_AUTO_DIFF_TERM_H
//
// This header defines DynamicAutoDiffTerm, in which both the number of variables and their sizes
// are known only at runtime.
//
//    class Functor {
//     public:
//      template <typename R>
//      R operator()(const std::vector<int>& dimensions, const R* const* const x) const {
//        // ...
//      }
//    };
//    vector<int> dimensions = {2, 3, 5};
//    DynamicAutoDiffTerm<Functor> my_term(dimensions);
//
// is equivalent to
//
//    class Functor {
//     public:
//      template <typename R>
//      R operator()(const R* const x, const R* const y, const R* const z) const {
//        // ...
//      }
//    };
//    AutoDiffTerm<Functor, 2, 3, 5> my_term();
//
// For a large number of scalars, DynamicAutoDiffTerm will efficiently compute gradients because a
// backwards differentiation is used. For a large number of variables, computing Hessians is not
// recommended due to the matrix necessarily becoming dense.
//
#include <vector>

#include <badiff.h>
#include <fadiff.h>

#include <minimum/core/check.h>
#include <minimum/nonlinear/term.h>
using minimum::core::check;

namespace minimum {
namespace nonlinear {

template <typename Functor>
class DynamicAutoDiffTerm final : public Term {
   public:
	template <typename... Args>
	DynamicAutoDiffTerm(std::vector<int> dimensions_, Args&&... args)
	    : functor(std::forward<Args>(args)...),
	      dimensions(std::move(dimensions_)),
	      total_size(create_total_size(dimensions)) {}

	static int create_total_size(const std::vector<int>& dimensions) {
		check(dimensions.size() >= 1, "Number of variables can not be 0.");
		int total_size = 0;
		for (auto d : dimensions) {
			check(d >= 1, "A variable dimension must be 1 or greater.");
			total_size += d;
		}
		return total_size;
	}

	int number_of_variables() const override { return dimensions.size(); }

	int variable_dimension(int var) const override { return dimensions[var]; }

	double evaluate(double* const* const variables) const override {
		return functor(dimensions, variables);
	}

	double evaluate(double* const* const variables,
	                std::vector<Eigen::VectorXd>* gradient) const override {
		using R = fadbad::B<double>;

		std::vector<R> x_data(total_size);
		std::vector<R*> x(dimensions.size());
		int pos = 0;
		for (int i = 0; i < dimensions.size(); ++i) {
			auto d = dimensions[i];
			x[i] = &x_data[pos];
			for (int j = 0; j < d; ++j) {
				x[i][j] = variables[i][j];
			}
			pos += d;
		}

		R f = functor(dimensions, x.data());
		f.diff(0, 1);

		for (int i = 0; i < dimensions.size(); ++i) {
			auto d = dimensions[i];
			for (int j = 0; j < d; ++j) {
				(*gradient)[i][j] = x[i][j].d(0);
			}
			pos += d;
		}

		return f.val();
	}

	double evaluate(double* const* const variables,
	                std::vector<Eigen::VectorXd>* gradient,
	                std::vector<std::vector<Eigen::MatrixXd>>* hessian) const override {
		using R = fadbad::B<fadbad::F<double>>;

		std::vector<R> x_data(total_size);
		std::vector<R*> x(dimensions.size());
		int pos = 0;
		for (int i = 0; i < dimensions.size(); ++i) {
			auto d = dimensions[i];
			x[i] = &x_data[pos];
			for (int j = 0; j < d; ++j) {
				x[i][j] = variables[i][j];
				x[i][j].x().diff(pos, total_size);
				pos++;
			}
		}

		R f = functor(dimensions, x.data());
		f.diff(0, 1);

		int pos1 = 0;
		for (int i1 = 0; i1 < dimensions.size(); ++i1) {
			auto d1 = dimensions[i1];
			for (int j1 = 0; j1 < d1; ++j1) {
				(*gradient)[i1][j1] = x[i1][j1].d(0).val();

				int pos2 = 0;
				for (int i2 = 0; i2 < dimensions.size(); ++i2) {
					auto d2 = dimensions[i2];
					for (int j2 = 0; j2 < d2; ++j2) {
						(*hessian)[i1][i2](j1, j2) = x[i1][j1].d(0).d(pos2);
						pos2++;
					}
				}
				pos1++;
			}
		}

		return f.val().val();
	}

   private:
	const Functor functor;
	const std::vector<int> dimensions;
	const int total_size;
};
}  // namespace nonlinear
}  // namespace minimum

#endif
