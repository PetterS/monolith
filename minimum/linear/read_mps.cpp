// Documentation of MPS:
// http://lpsolve.sourceforge.net/5.5/mps-format.htm

#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <minimum/core/check.h>
#include <minimum/core/string.h>
#include <minimum/linear/ip.h>

using minimum::core::check;
using minimum::core::from_string;
using minimum::core::to_string;

namespace minimum {
namespace linear {

namespace {
std::string trim(const std::string& s) {
	std::string out;
	auto first = s.find_first_not_of(' ');
	auto last = s.find_last_not_of(' ') + 1;
	if (first < last) {
		return s.substr(first, last - first);
	}
	return "";
}
}  // namespace

std::unique_ptr<IP> read_MPS(std::istream& in) {
	using namespace std;
	auto ip = make_unique<IP>();

	check(bool(in), "Not a valid stream.");
	struct RowInfo {
		char type = '-';
		Sum sum = 0;
		double rhs = 0;
		double range = std::numeric_limits<double>::quiet_NaN();
	};
	unordered_map<string, RowInfo> rows;
	map<string, Variable> vars;
	IP::VariableType variable_type = IP::Real;
	unordered_set<string> variables_not_mentioned;
	unordered_set<string> bound_types;

	string line;
	string section;
	while (getline(in, line)) {
		if (line.empty()) {
			continue;
		}
		if (line[0] == '*') {
			continue;
		}

		if (line[0] != ' ') {
			stringstream lin(line);
			lin >> ws >> section;
			continue;
		}

		if (line.length() < 61) {
			line.resize(61, ' ');
		}
		string field1 = trim(line.substr(1, 2));
		string field2 = trim(line.substr(4, 8));
		string field3 = trim(line.substr(14, 8));
		string field4 = trim(line.substr(24, 12));
		string field5 = trim(line.substr(39, 8));
		string field6 = trim(line.substr(49, 12));

		if (section == "NAME") {
			// Do not process this section.
		} else if (section == "ROWS") {
			auto& type = field1;
			auto& name = field2;
			check(type.length() == 1, "Invalid row type: ", type);
			check(!name.empty(), "Need a row name.");
			auto& row = rows[name];
			row.type = type[0];
		} else if (section == "COLUMNS") {
			auto& varname = field2;
			check(!varname.empty(), "Need a variable name.");

			if (field3 == "\'MARKER\'") {
				auto& marker_type = field5;
				if (marker_type == "\'INTORG\'") {
					variable_type = IP::Integer;
				} else if (marker_type == "\'INTEND\'") {
					variable_type = IP::Real;
				} else {
					check(false, "Unknown marker type: \"", marker_type, "\"");
				}
			} else {
				string rowname[2] = {field3, field5};
				string value_string[2] = {field4, field6};

				for (int a = 0; a < 2; ++a) {
					if (!rowname[a].empty()) {
						auto variter = vars.find(varname);
						if (variter == end(vars)) {
							variter =
							    vars.emplace(make_pair(varname, ip->add_variable(variable_type)))
							        .first;
							ip->add_bounds(-1e100, variter->second, 1e100);
							variables_not_mentioned.insert(varname);
						}
						auto row = rows.find(rowname[a]);
						check(row != end(rows), "Could not find row: ", rowname[a]);
						row->second.sum += from_string<double>(value_string[a]) * variter->second;
					}
				}
			}
		} else if (section == "RHS") {
			string rowname[2] = {field3, field5};
			string value_string[2] = {field4, field6};
			for (int a = 0; a < 2; ++a) {
				if (!rowname[a].empty()) {
					auto row = rows.find(rowname[a]);
					check(row != end(rows), "Could not find row: ", rowname[a]);
					row->second.rhs = from_string<double>(value_string[a]);
				}
			}
		} else if (section == "BOUNDS") {
			auto& type = field1;
			auto& varname = field3;
			auto& bound_string = field4;

			bound_types.insert(type);
			check(!varname.empty(), "Could not parse bound: \"", line, "\".");
			auto variter = vars.find(varname);
			check(variter != end(vars), "Could not find variable: ", varname);

			if (type == "PL") {
				// Default bound 0 <= x is added afterwards.
				continue;
			}

			if (type == "FR") {
				// Free variable.
				variables_not_mentioned.erase(varname);
				continue;
			} else if (type == "BV") {
				// Binary variable.
				variables_not_mentioned.erase(varname);
				ip->add_bounds(0, variter->second, 1);
				// TODO: Set binary.
				continue;
			}

			auto bound = from_string<double>(bound_string);
			if (type == "FX") {
				variables_not_mentioned.erase(varname);
				ip->add_bounds(bound, variter->second, bound);
			} else if (type == "UP") {
				ip->add_constraint(variter->second <= bound);
			} else if (type == "LO") {
				variables_not_mentioned.erase(varname);
				ip->add_constraint(bound <= variter->second);
			} else {
				check(false, "Unknown bound type: ", type);
			}
		} else if (section == "RANGES") {
			string rowname[2] = {field3, field5};
			string value_string[2] = {field4, field6};
			for (int a = 0; a < 2; ++a) {
				if (!rowname[a].empty()) {
					auto row = rows.find(rowname[a]);
					check(row != end(rows), "Could not find row: ", rowname[a]);
					row->second.range = from_string<double>(value_string[a]);
				}
			}
		} else {
			check(false, "Unknown section: ", section);
		}
	}

	for (auto& name : variables_not_mentioned) {
		ip->add_bounds(0, vars[name], 1e100);
	}

	for (auto& row : rows) {
		double lower_bound = -1e100;
		double upper_bound = 1e100;

		if (row.second.type == 'E') {
			lower_bound = row.second.rhs;
			upper_bound = row.second.rhs;
		} else if (row.second.type == 'L') {
			upper_bound = row.second.rhs;
		} else if (row.second.type == 'G') {
			lower_bound = row.second.rhs;
		} else if (row.second.type == 'N') {
			// E226 becomes incorrect if we add this.
			// row.second.sum += row.second.rhs;
		} else {
			check(false, "Unknown row type: ", row.second.type);
		}

		if (row.second.range == row.second.range) {
			if (row.second.type == 'E') {
				if (row.second.range >= 0) {
					upper_bound += row.second.range;
				} else {
					lower_bound += row.second.range;
				}
			} else if (row.second.type == 'L') {
				lower_bound = upper_bound - std::abs(row.second.range);
			} else if (row.second.type == 'G') {
				upper_bound = lower_bound + std::abs(row.second.range);
			} else {
				check(false, "Unknown RANGES type.");
			}
		}

		if (row.second.type == 'N') {
			ip->add_objective(row.second.sum);
		} else {
			// TODO: Make the first-order solver do this conversion instead.
			if (lower_bound > -1e100 && upper_bound < 1e100) {
				ip->add_constraint(row.second.sum <= upper_bound);
				ip->add_constraint(lower_bound <= row.second.sum);
			} else {
				ip->add_constraint(lower_bound, row.second.sum, upper_bound);
			}
		}
	}
	return ip;
}
}  // namespace linear
}  // namespace minimum
