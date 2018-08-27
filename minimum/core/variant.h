#pragma once
#include <utility>

namespace minimum::core {

// Combines multiple lambdas into one lambda that has
// overloaded operator().
// Useful for std::visit.
template <typename... Lambdas>
struct overload : Lambdas... {
	using Lambdas::operator()...;
	constexpr overload(Lambdas&&... lambdas) : Lambdas(std::forward<Lambdas>(lambdas))... {}
};

// Switch for variants.
// switch_(opt,
//	[](int i) { ...},
//  [](string_view s) { ... });
template <typename Variant, typename... Lambdas>
constexpr auto switch_(Variant&& var, Lambdas&&... lambdas) {
	return std::visit(overload{std::forward<Lambdas>(lambdas)...}, std::forward<Variant>(var));
}
}  // namespace minimum::core
