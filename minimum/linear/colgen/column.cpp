#include <cmath>
#include <vector>
using namespace std;

#include <minimum/core/check.h>
#include <minimum/linear/colgen/column.h>
using namespace minimum::core;

namespace minimum {
namespace linear {
namespace colgen {

class Column::Implementation {
   public:
	vector<RowEntry> rows;
	double cost = 0;
	double lower_bound = 0;
	double upper_bound = 1e100;
	int fixed_value = -1;
	bool is_integer = true;
};

Column::Column() : impl(nullptr) {}

Column::Column(double cost, double lower_bound, double upper_bound) : impl(new Implementation) {
	impl->cost = cost;
	impl->lower_bound = lower_bound;
	impl->upper_bound = upper_bound;
}

Column::Column(Column&& rhs) noexcept {
	impl = rhs.impl;
	solution_value = rhs.solution_value;
	rhs.impl = nullptr;
}

Column::~Column() { delete impl; }

Column& Column::operator=(Column&& rhs) {
	delete impl;
	impl = rhs.impl;
	solution_value = rhs.solution_value;
	rhs.impl = nullptr;
	return *this;
}

double Column::lower_bound() const {
	if (is_fixed()) {
		return impl->fixed_value;
	}
	return impl->lower_bound;
}

double Column::upper_bound() const {
	if (is_fixed()) {
		return impl->fixed_value;
	}
	return impl->upper_bound;
}

double Column::cost() const { return impl->cost; }

void Column::add_coefficient(int row, double coef) { impl->rows.emplace_back(row, coef); }

void Column::set_integer(bool is_integer) {
	check(!(is_fixed() && !is_integer), "Can not set a fixed column to be real-valued.");
	impl->is_integer = is_integer;
}

bool Column::is_integer() const { return impl->is_integer; }

void Column::fix(int value) {
	check(value >= 0, "Can not fix column to negative value.");
	check(impl->is_integer, "Can not fix real-valued column.");
	solution_value = value;
	impl->fixed_value = value;
}

void Column::unfix() { impl->fixed_value = -1; }

bool Column::is_fixed() const { return impl->fixed_value >= 0; }

double Column::reduced_cost(const std::vector<double>& dual_variables) const {
	double value = impl->cost;
	for (auto& entry : impl->rows) {
		value -= entry.coef * dual_variables[entry.row];
	}
	return value;
}

RowEntry* Column::begin() { return impl->rows.data(); }

RowEntry* Column::end() { return impl->rows.data() + impl->rows.size(); }

const RowEntry* Column::begin() const { return impl->rows.data(); }

const RowEntry* Column::end() const { return impl->rows.data() + impl->rows.size(); }

bool Column::operator==(const Column& rhs) const {
	return abs(impl->cost - rhs.impl->cost) < 1e-9
	       && abs(impl->upper_bound - rhs.impl->upper_bound) < 1e-9
	       && abs(impl->lower_bound - rhs.impl->lower_bound) < 1e-9 && impl->rows == rhs.impl->rows;
}

proto::Column Column::to_proto() const {
	proto::Column proto_column;
	proto_column.set_cost(impl->cost);
	proto_column.set_lower_bound(impl->lower_bound);
	proto_column.set_upper_bound(impl->upper_bound);
	for (auto& entry : impl->rows) {
		auto proto_entry = proto_column.add_entry();
		proto_entry->set_row(entry.row);
		proto_entry->set_value(entry.coef);
	}
	return proto_column;
}

Column Column::from_proto(const proto::Column& proto_column) {
	Column column(proto_column.cost(), proto_column.lower_bound(), proto_column.upper_bound());
	for (auto& entry : proto_column.entry()) {
		column.add_coefficient(entry.row(), entry.value());
	}
	return column;
}
}  // namespace colgen
}  // namespace linear
}  // namespace minimum
