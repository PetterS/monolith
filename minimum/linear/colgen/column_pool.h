#pragma once

#include <iostream>
#include <vector>

#include <minimum/linear/colgen/column.h>
#include <minimum/linear/colgen/export.h>

namespace minimum {
namespace linear {
namespace colgen {

struct ColumnScore {
	std::size_t index = -1;
	double reduced_cost = 0;
	bool operator<(const ColumnScore& rhs) const { return reduced_cost < rhs.reduced_cost; }
	bool operator==(const ColumnScore& rhs) const {
		return reduced_cost == rhs.reduced_cost && index == rhs.index;
	}
};

// This container owns all columns for a complete run of column generation. No column is
// ever purged from the pool.
class MINIMUM_LINEAR_COLGEN_API ColumnPool {
   public:
	ColumnPool();
	ColumnPool(std::istream* in);
	~ColumnPool();
	ColumnPool(const ColumnPool&) = delete;

	std::size_t size() const;
	const Column& at(std::size_t index) const;
	Column& at(std::size_t index);

	// Number of columns that are not fixed to zero.
	std::size_t allowed_size() const;

	// Adds a column to the pool. If an equal column already exists, a duplicate one
	// will not be added.
	void add(Column&& column);

	// Returns the columns with negative reduced costs, sorted with the most negative
	// reduced cost first.
	// If max_count > 0, no more than that are returned.
	const std::vector<ColumnScore>& get_sorted(const std::vector<double>& dual_variables,
	                                           std::size_t max_count = 0);

	const Column* begin() const;
	const Column* end() const;
	Column* begin();
	Column* end();

	// Resets the pool.
	void clear();
	// Save all columns to a stream. Does not save fix state or the
	// solution value.
	void save_to_stream(std::ostream* out) const;

   private:
	class Implementation;
	Implementation* impl;
};
}  // namespace colgen
}  // namespace linear
}  // namespace minimum
