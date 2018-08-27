#pragma once

#include <minimum/core/parser.h>
#include <minimum/core/range.h>
using minimum::core::range;

namespace minimum::nonlinear {
// Functor that can be used with AutoDiffTerm or DynamicAutoDiffTerm.
//
// It parses a string at runtime to determine the expression to evaluate.
class StringFunctor1 {
   public:
	StringFunctor1(std::string_view expr, std::vector<std::string> names_)
	    : expression(minimum::core::Parser(expr).parse()), names(names_) {}
	StringFunctor1(minimum::core::Expression expression_, std::vector<std::string> names_)
	    : expression(std::move(expression_)), names(names_) {}

	template <typename R>
	R operator()(const R* const x) const {
		std::unordered_map<std::string, R> vars;
		for (auto i : range(names.size())) {
			vars[names[i]] = x[i];
		}
		return expression.evaluate<R>(vars);
	}

	template <typename R>
	R operator()(const std::vector<int>& dimensions, const R* const* const x) const {
		minimum_core_assert(dimensions.size() == 1);
		minimum_core_assert(dimensions[0] == names.size());

		std::unordered_map<std::string, R> vars;
		for (auto i : range(names.size())) {
			vars[names[i]] = x[0][i];
		}
		return expression.evaluate<R>(vars);
	}

   private:
	const minimum::core::Expression expression;
	const std::vector<std::string> names;
};

class StringFunctorN {
   public:
	StringFunctorN(std::string_view expr, std::vector<std::string> names_)
	    : expression(minimum::core::Parser(expr).parse()), names(names_) {}
	StringFunctorN(minimum::core::Expression expression_, std::vector<std::string> names_)
	    : expression(std::move(expression_)), names(names_) {}

	template <typename R>
	R operator()(const std::vector<int>& dimensions, const R* const* const x) const {
		minimum_core_assert(dimensions.size() == names.size());

		std::unordered_map<std::string, R> vars;
		for (auto i : range(names.size())) {
			minimum_core_assert(dimensions[i] == 1);
			vars[names[i]] = x[i][0];
		}
		return expression.evaluate<R>(vars);
	}

   private:
	const minimum::core::Expression expression;
	const std::vector<std::string> names;
};
}  // namespace minimum::nonlinear
