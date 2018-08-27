#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>

#include <minimum/core/check.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/proto.h>
#include <minimum/linear/sat_solver.h>
#include <minimum/linear/solver.h>
#include <minimum/linear/sum.h>
#include <minimum/linear/variable.h>

namespace minimum {
namespace linear {

class IP::Implementation {
   public:
	Implementation(IP* creator_) : creator(creator_) {}

	proto::IP ip;
	proto::Solution solution;

	vector<PseudoBooleanConstraint> pseudoboolean_constraints;
	vector<PseudoBoolean> pseudoboolean_objective;
	// Holds the translation of monomials to extra added variables.
	// Returns a sum because the empty monomial is translated as 1.
	std::map<vector<int>, const Sum> monomial_to_sum;

	// The type of each variable.
	vector<VariableType> variable_types;

	bool is_in_exists_mode = false;
	std::vector<Variable> exists_variables;

	void check_creator(const Variable& t) const;
	void check_creator(const Sum& t) const;
	IP* creator;
};

template <typename Iterator>
void print_histogram_01(Iterator begin, Iterator end) {
	std::vector<int> bins(10, 0);
	int zeros = 0;
	int ones = 0;
	int n = 0;
	for (; begin != end; ++begin) {
		n++;
		auto value = *begin;
		if (value <= 1e-7) {
			zeros++;
		} else if (std::abs(value - 1) <= 1e-7) {
			ones++;
		}
		int bin = int(std::max(std::min(value, 0.999), 0.0) * 10);
		bins.at(bin)++;
	}

	std::clog
	    << "(ex. 0)  0    0.1    0.2    0.3    0.4    0.5    0.6    0.7    0.8    0.9    1.0   "
	       "(ex. 1)\n";
	auto percent = [n](double count) {
		return core::to_string(std::right,
		                       std::setw(3),
		                       std::fixed,
		                       std::setprecision(0),
		                       100.0 * double(count) / n,
		                       "%");
	};
	auto exact = [](int count) {
		char postfix = 'G';
		if (count <= 999) {
			postfix = ' ';
		} else if (count <= 999'999) {
			postfix = 'k';
		} else if (count <= 999'999'999) {
			postfix = 'M';
		}
		while (count >= 1000) {
			count /= 1000;
		}
		return core::to_string(std::right, std::setw(3), count, postfix);
	};

	std::clog << " " << percent(zeros) << "   |";
	for (int i = 0; i < 10; ++i) {
		std::clog << percent(bins.at(i)) << "  |";
	}
	std::clog << "    " << percent(ones) << "\n";
	std::clog << " " << exact(zeros) << "   |";
	for (int i = 0; i < 10; ++i) {
		std::clog << " " << exact(bins.at(i)) << " |";
	}
	std::clog << "    " << exact(ones) << "\n";
}
}  // namespace linear
}  // namespace minimum
