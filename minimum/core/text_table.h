#pragma once
#include <iomanip>
#include <iostream>

#include <minimum/core/check.h>
#include <minimum/core/range.h>

template <typename T>
void print_table(std::ostream& out,
                 const std::vector<std::vector<T>> data,
                 const std::vector<std::vector<bool>>& vertical_border,
                 const std::vector<std::vector<bool>>& horizontal_border) {
	using minimum::core::check;
	using minimum::core::range;

	const auto m = data.size();
	check(m >= 1, "There is no data.");
	const auto n = data[0].size();
	for (auto& row : data) {
		check(row.size() == n, "Input table must be rectangular.");
	}
	check(vertical_border.size() == m, "Vertical border must be m x n+1.");
	for (auto& row : vertical_border) {
		check(row.size() == n + 1, "Vertical border must be m x n+1.");
	}
	check(horizontal_border.size() == m + 1, "Horizontal border border must be m+1 x n.");
	for (auto& row : horizontal_border) {
		check(row.size() == n, "Horizontal border must be m+1 x n.");
	}

	for (int j : range(n)) {
		auto horz = horizontal_border[0][j] ? '-' : ' ';
		auto vert = vertical_border[0][j] ? '|' : ' ';
		if (vert != ' ') {
			out << '+';
		} else {
			out << horz;
		}
		out << horz << horz;
	}
	out << "+\n";

	for (int i : range(m)) {
		out << (vertical_border[i][0] ? '|' : ' ');
		for (int j : range(n)) {
			out << std::setw(2) << data[i][j] << (vertical_border[i][j + 1] ? '|' : ' ');
		}
		out << "\n";
		if (!horizontal_border[i + 1][0]) {
			out << (vertical_border[i][0] ? '|' : ' ');
		} else {
			out << '+';
		}
		for (int j : range(n)) {
			char horz1 = horizontal_border[i + 1][j] ? '-' : ' ';
			char horz2 = ' ';
			if (j < n - 1) {
				horz2 = horizontal_border[i + 1][j + 1] ? '-' : ' ';
			}
			out << horz1 << horz1;
			char vert1 = vertical_border[i][j + 1] ? '|' : ' ';
			char vert2 = ' ';
			if (i < m - 1) {
				vert2 = vertical_border[i + 1][j + 1] ? '|' : ' ';
			}
			bool any_vert = vert1 != ' ' || vert2 != ' ';
			bool any_horz = horz1 != ' ' || horz2 != ' ';
			if (any_vert && any_horz) {
				out << '+';
			} else if (vert1 != ' ') {
				out << '|';
			} else if (any_horz) {
				out << '-';
			} else {
				out << ' ';
			}
		}
		out << "\n";
	}
}

template <typename T>
void print_table(std::ostream& out, const std::vector<std::vector<T>> data) {
	using minimum::core::check;
	using minimum::core::range;

	const auto m = data.size();
	check(m >= 1, "There is no data.");
	const auto n = data[0].size();
	for (auto& row : data) {
		check(row.size() == n, "Input table must be square.");
	}

	std::vector<std::vector<bool>> vertical_border(m, std::vector<bool>(n + 1, false));
	std::vector<std::vector<bool>> horizontal_border(m + 1, std::vector<bool>(n, false));
	horizontal_border[0] = std::vector<bool>(n, true);
	horizontal_border[m] = std::vector<bool>(n, true);
	for (int i : range(m)) {
		vertical_border[i][0] = true;
		vertical_border[i][n] = true;
	}
	for (int i : range(m)) {
		for (int j : range(n - 1)) {
			if (data[i][j] != data[i][j + 1]) {
				vertical_border[i][j + 1] = true;
			}
		}
	}
	for (int i : range(m - 1)) {
		for (int j : range(n)) {
			if (data[i][j] != data[i + 1][j]) {
				horizontal_border[i + 1][j] = true;
			}
		}
	}

	print_table(out, data, vertical_border, horizontal_border);
}
