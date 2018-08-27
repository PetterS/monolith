#pragma once
#include <cmath>

namespace minimum {
namespace core {

template <typename Float>
Float relative_error(Float value, Float truth) {
	return std::abs(value - truth) / (1 + std::abs(truth));
}

// Whether value lies roughly within the bounds.
template <typename Float>
bool is_feasible(Float lower_bound, Float value, Float upper_bound, Float eps = 1e-9) {
	if (value < lower_bound) {
		return relative_error(value, lower_bound) <= eps;
	} else if (value > upper_bound) {
		return relative_error(value, upper_bound) <= eps;
	} else {
		return true;
	}
}

// Whether value is roughly equal to truth.
template <typename Float>
bool is_equal(Float value, Float truth, Float eps = 1e-9) {
	return is_feasible(truth, value, truth, eps);
}

// The distance from the value to the nearest integer.
template <typename Float>
Float integer_residual(Float value) {
	auto integer_value = std::round(value);
	return std::abs(value - integer_value);
}
}  // namespace core
}  // namespace minimum
