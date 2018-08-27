#pragma once

#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <minimum/core/check.h>

namespace minimum {
namespace core {

namespace internal {
constexpr bool is_separator(char ch) {
	return ch == ' ' || ch == ';' || ch == ',' || ch == '\xA0' || ch == '\n' || ch == '\r'
	       || ch == '\t';
}

constexpr bool is_comment(char ch) { return ch == '#' || ch == '/'; }

constexpr bool is_floating_point_character(char ch) {
	return ('0' <= ch && ch <= '9') || ch == 'n' || ch == 'N' || ch == '+' || ch == '-' || ch == 'a'
	       || ch == 'A' || ch == 'e' || ch == 'E' || ch == '.' || ch == 'i' || ch == 'I'
	       || ch == 'f' || ch == 'F';
}
}  // namespace internal

// Reads a matrix from a stream where each row in the stream will be one row in the matrix.
// The numbers are separated by e.g. spaces, tabls, commas, or semicolons. A blank line can
// either finish the matrix or be ignored.
//
// This function is hand-rolled in a seemingly complicated way after careful profiling.
template <typename Real>
std::vector<std::vector<Real>> read_matrix(std::istream& fin, bool stop_after_blank_line = false) {
	using std::vector;

	minimum_core_assert(fin);
	size_t size = 10000;
	if (stop_after_blank_line) {
		size = 1;
	}
	vector<char> buffer(size);

	vector<vector<Real>> matrix;
	matrix.emplace_back();
	size_t N = 0;

	Real number = 0;
	bool negative = false;
	bool has_number = false;

	size_t read = 0;
	bool in_comment = false;
	do {
		fin.read(buffer.data(), size);
		read = fin.gcount();
		for (size_t i = 0; i < read; ++i) {
			auto ch = buffer[i];

			if (in_comment && ch == '\n') {
				in_comment = false;
			} else if (in_comment) {
				continue;
			}

			if (internal::is_comment(ch)) {
				in_comment = true;
			} else if ('0' <= ch && ch <= '9') {
				has_number = true;
				if (number >= std::numeric_limits<Real>::max() / 10 - 10) {
					throw std::runtime_error("Number too big.");
				}
				number *= Real(10);
				number += Real(ch - '0');
			} else if (ch == '-') {
				minimum_core_assert(!has_number && !negative);
				negative = true;
			} else if (ch == '+') {
				minimum_core_assert(!has_number);
			} else {
				minimum_core_assert(internal::is_separator(ch));

				if (has_number) {
					matrix.back().emplace_back(negative ? -number : number);
					number = 0;
					negative = false;
					has_number = false;
				}
				if (ch == '\n') {
					if (!matrix.back().empty()) {
						if (N == 0) {
							N = matrix.back().size();
						} else {
							minimum_core_assert(N == matrix.back().size());
						}
						matrix.emplace_back();
						matrix.back().reserve(N);
					} else if (stop_after_blank_line) {
						read = 0;  // Will break outer loop.
						break;
					}
				}
			}
		}
	} while (read > 0);

	if (has_number) {
		matrix.back().emplace_back(negative ? -number : number);
	}

	if (!matrix.back().empty()) {
		if (N > 0) {
			minimum_core_assert(N == matrix.back().size());
		}
	} else {
		matrix.pop_back();
	}

	return matrix;
}

namespace internal {
template <typename Real>
std::vector<std::vector<Real>> read_matrix_floating_point(std::istream& fin,
                                                          bool stop_after_blank_line) {
	using std::vector;

	minimum_core_assert(fin);
	size_t size = 10000;
	if (stop_after_blank_line) {
		size = 1;
	}
	vector<char> buffer(size);

	vector<vector<Real>> matrix;
	matrix.emplace_back();
	size_t N = 0;

	vector<char> number_data;
	number_data.reserve(100);

	size_t read = 0;
	bool in_comment = false;
	do {
		fin.read(buffer.data(), size);
		read = fin.gcount();
		for (size_t i = 0; i < read; ++i) {
			auto ch = buffer[i];

			if (in_comment && ch == '\n') {
				in_comment = false;
			} else if (in_comment) {
				continue;
			}

			if (is_comment(ch)) {
				in_comment = true;
			} else if (is_floating_point_character(ch)) {
				number_data.emplace_back(ch);
			} else {
				minimum_core_assert(is_separator(ch));

				if (!number_data.empty()) {
					number_data.emplace_back('\0');
					matrix.back().emplace_back((Real)std::atof(number_data.data()));
					number_data.clear();
				}
				if (ch == '\n') {
					if (!matrix.back().empty()) {
						if (N == 0) {
							N = matrix.back().size();
						} else {
							minimum_core_assert(N == matrix.back().size());
						}
						matrix.emplace_back();
						matrix.back().reserve(N);
					} else if (stop_after_blank_line) {
						read = 0;  // Will break outer loop.
						break;
					}
				}
			}
		}
	} while (read > 0);

	if (!number_data.empty()) {
		number_data.emplace_back('\0');
		matrix.back().emplace_back((Real)std::atof(number_data.data()));
	}

	if (!matrix.back().empty()) {
		if (N > 0) {
			minimum_core_assert(N == matrix.back().size());
		}
	} else {
		matrix.pop_back();
	}

	return matrix;
}
}  // namespace internal

template <>
std::vector<std::vector<double>> read_matrix(std::istream& fin, bool stop_after_blank_line) {
	return internal::read_matrix_floating_point<double>(fin, stop_after_blank_line);
}

template <>
std::vector<std::vector<float>> read_matrix(std::istream& fin, bool stop_after_blank_line) {
	return internal::read_matrix_floating_point<float>(fin, stop_after_blank_line);
}
}  // namespace core
}  // namespace minimum
