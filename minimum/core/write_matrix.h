#pragma once

#include <iostream>

namespace minimum {
namespace core {

// Writes a matrix to a stream. Not heavily optimized like read_matrix.
template <typename Container2D>
void write_matrix(std::ostream& fin, Container2D&& matrix) {
	for (auto&& row : matrix) {
		for (auto&& elem : row) {
			fin << elem << "\t";
		}
		fin << "\n";
	}
}
}  // namespace core
}  // namespace minimum
