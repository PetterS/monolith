#include <sstream>
#include <vector>
using namespace std;

#include <catch.hpp>

#include <minimum/linear/colgen/column_pool.h>
using namespace minimum::linear::colgen;

static_assert(std::is_nothrow_move_constructible<Column>::value,
              "Column should have a noexcept move constructor to work with std::vector optimally.");

TEST_CASE("small") {
	Column column0(0, 0, 1);
	column0.add_coefficient(0, 1);
	Column column1(0, 0, 1);
	column1.add_coefficient(1, 1);
	Column column2(0, 0, 1);
	column2.add_coefficient(2, -1);

	ColumnPool org_pool;
	org_pool.add(move(column0));
	org_pool.add(move(column1));
	org_pool.add(move(column2));

	REQUIRE(org_pool.size() == 3);
	stringstream str;
	org_pool.save_to_stream(&str);
	org_pool.clear();
	REQUIRE(org_pool.size() == 0);

	ColumnPool pool(&str);
	REQUIRE(pool.size() == 3);

	vector<double> dual = {1, 10, 100};
	auto sorted = pool.get_sorted(dual);
	REQUIRE(sorted.size() == 2);
	CHECK(sorted[0].index == 1);
	CHECK(sorted[0].reduced_cost == -10);
	CHECK(sorted[1].index == 0);
	CHECK(sorted[1].reduced_cost == -1);

	auto sorted1 = pool.get_sorted(dual, 1);
	REQUIRE(sorted1.size() == 1);
	CHECK(sorted1[0] == sorted[0]);

	// Test begin/end.
	int c = 0;
	for (auto& column : pool) {
		c++;
	}
	CHECK(c == 3);

	// Do not return zero-fixed columns.
	pool.at(0).fix(0);
	pool.at(1).fix(0);
	CHECK(pool.get_sorted(dual).empty());
}

TEST_CASE("duplicates") {
	Column column0(0, 0, 1);
	column0.add_coefficient(0, 1);
	Column column1(0, 0, 1);
	column1.add_coefficient(1, 1);
	Column column00(0, 0, 1);
	column00.add_coefficient(0, 1);
	Column column_different_cost(5, 0, 1);
	column_different_cost.add_coefficient(0, 1);
	Column column_different_lower_bound(0, 0.5, 1);
	column_different_lower_bound.add_coefficient(0, 1);
	Column column_different_upper_bound(0, 0, 0.5);
	column_different_upper_bound.add_coefficient(0, 1);

	ColumnPool pool;
	pool.add(move(column0));
	pool.add(move(column1));
	pool.add(move(column00));
	pool.add(move(column_different_cost));
	pool.add(move(column_different_lower_bound));
	pool.add(move(column_different_upper_bound));

	REQUIRE(pool.size() == 5);
}

TEST_CASE("const_member_functions") {
	ColumnPool pool;
	pool.add(Column(0, 0, 1));
	const ColumnPool& cpool = pool;
	CHECK(*pool.begin() == *cpool.begin());
	CHECK(pool.at(0) == cpool.at(0));
	CHECK((pool.end() - pool.begin()) == (cpool.end() - cpool.begin()));
}

TEST_CASE("can_not_fix_real_valued_column") {
	Column column(0, 0, 1);
	column.set_integer(false);
	CHECK_THROWS(column.fix(1));
}

TEST_CASE("can_not_set_fixed_column_to_be_read_valued") {
	Column column(0, 0, 1);
	column.fix(1);
	CHECK_THROWS(column.set_integer(false));
}
