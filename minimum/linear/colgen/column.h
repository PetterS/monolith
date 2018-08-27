#pragma once
#include <vector>

#include <minimum/linear/colgen/export.h>
#include <minimum/linear/colgen/proto.h>

namespace minimum {
namespace linear {
namespace colgen {

struct RowEntry {
	RowEntry(int row_, double coef_) : row(row_), coef(coef_) {}
	int row;
	double coef;
	bool operator==(const RowEntry& rhs) const { return row == rhs.row && coef == rhs.coef; }
};

// Holds one column in a column generation problem.
class MINIMUM_LINEAR_COLGEN_API Column {
   public:
	double solution_value = 0;

	// Creates an empty placeholder column. It is undefined to perform any operations
	// on it except assigning to it.
	Column();
	Column(double cost, double lower_bound, double upper_bound);
	Column(Column&&) noexcept;
	Column(const Column&) = delete;
	~Column();
	Column& operator=(Column&&);
	Column& operator=(const Column&) = delete;

	double lower_bound() const;
	double upper_bound() const;
	double cost() const;

	void add_coefficient(int row, double coef);

	// Whether this column should be an integer in the final
	// solution. (This affects when the optimization stops.)
	//  Default: true.
	void set_integer(bool is_integer);
	bool is_integer() const;

	// Sets the lower and upper bounds to a non-negative integer. Can
	// be undone with unfix().
	void fix(int value);
	void unfix();
	bool is_fixed() const;

	double reduced_cost(const std::vector<double>& dual_variables) const;

	RowEntry* begin();
	RowEntry* end();
	const RowEntry* begin() const;
	const RowEntry* end() const;

	// Only the column entries in the system matrix will be considered for
	// comparison. If only the cost or bounds are different, the columns
	// are considered equal.
	bool operator==(const Column& rhs) const;

	// The proto does not store fix state or the solution value.
	proto::Column to_proto() const;
	static Column from_proto(const proto::Column& proto_column);

   private:
	class Implementation;
	Implementation* impl;
};
}  // namespace colgen
}  // namespace linear
}  // namespace minimum
