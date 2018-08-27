#include <algorithm>
#include <iostream>
using namespace std;

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/core/record_stream.h>
#include <minimum/linear/colgen/column_pool.h>
#include <minimum/linear/colgen/proto.h>
using namespace minimum::core;

namespace minimum {
namespace linear {
namespace colgen {

class ColumnPool::Implementation {
   public:
	vector<ColumnScore> column_scores;
	vector<Column> columns;
};

ColumnPool::ColumnPool() : impl(new Implementation) {}

ColumnPool::ColumnPool(std::istream* in) : impl(new Implementation) {
	string tmp;
	read_record(in, &tmp);
	check(tmp == "ColumnPool");
	while (true) {
		read_record(in, &tmp);
		if (tmp.empty()) {
			break;
		}
		proto::Column column;
		check(column.ParseFromArray(tmp.data(), tmp.size()), "Could not parse column.");
		impl->columns.emplace_back(Column::from_proto(column));
	}
}

ColumnPool::~ColumnPool() { delete impl; }

void ColumnPool::add(Column&& column) {
	// TODO: Check if the columns exists in a smarter way.
	if (find(impl->columns.begin(), impl->columns.end(), column) != impl->columns.end()) {
		return;
	}
	impl->columns.emplace_back(move(column));
}

size_t ColumnPool::size() const { return impl->columns.size(); }

const Column& ColumnPool::at(size_t index) const { return impl->columns.at(index); }

Column& ColumnPool::at(std::size_t index) { return impl->columns.at(index); }

size_t ColumnPool::allowed_size() const {
	size_t num = 0;
	for (auto i : range(size())) {
		// Do not include columns that are fixed to zero.
		if (impl->columns[i].upper_bound() > 1e-9) {
			num++;
		}
	}
	return num;
}

const vector<ColumnScore>& ColumnPool::get_sorted(const vector<double>& dual_variables,
                                                  std::size_t max_count) {
	impl->column_scores.clear();
	for (auto i : range(size())) {
		// Do not include columns that are fixed to zero.
		if (impl->columns[i].upper_bound() > 1e-9) {
			impl->column_scores.emplace_back();
			impl->column_scores.back().index = i;
			impl->column_scores.back().reduced_cost = impl->columns[i].reduced_cost(dual_variables);
		}
	}
	sort(impl->column_scores.begin(), impl->column_scores.end());

	if (max_count > 0 && impl->column_scores.size() > max_count) {
		impl->column_scores.resize(max_count);
	}

	for (auto i : range(impl->column_scores.size())) {
		if (impl->column_scores[i].reduced_cost >= 0) {
			impl->column_scores.resize(i);
			break;
		}
	}

	return impl->column_scores;
}

Column* ColumnPool::begin() { return impl->columns.data(); }

Column* ColumnPool::end() { return impl->columns.data() + impl->columns.size(); }

const Column* ColumnPool::begin() const { return impl->columns.data(); }

const Column* ColumnPool::end() const { return impl->columns.data() + impl->columns.size(); }

void ColumnPool::clear() {
	impl->columns.clear();
	impl->column_scores.clear();
}

void ColumnPool::save_to_stream(std::ostream* out) const {
	write_record(out, "ColumnPool");
	string data;
	for (auto& column : impl->columns) {
		column.to_proto().SerializeToString(&data);
		write_record(out, data);
	}
	write_record(out, "");
}
}  // namespace colgen
}  // namespace linear
}  // namespace minimum
