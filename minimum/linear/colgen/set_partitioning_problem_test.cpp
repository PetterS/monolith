#include <numeric>
#include <vector>
using namespace std;

#include <gflags/gflags_declare.h>
#include <catch.hpp>

#include <minimum/linear/colgen/set_partitioning_problem.h>
using namespace minimum::linear::colgen;

DECLARE_double(fix_threshold);
DECLARE_double(objective_change_before_fixing);

class TestSetPartitioningProblem : public SetPartitioningProblem {
   public:
	using SetPartitioningProblem::add_column;
	using SetPartitioningProblem::column_allowed;
	using SetPartitioningProblem::compute_fractional_solution;
	using SetPartitioningProblem::fix_using_columns;
	using SetPartitioningProblem::fixes_for_member;
	using SetPartitioningProblem::fractional_solution;
	using SetPartitioningProblem::integral_solution_value;
	using SetPartitioningProblem::number_of_rows;
	using SetPartitioningProblem::SetPartitioningProblem;

	using Problem::pool;

	virtual void generate(const std::vector<double>& dual_variables) override {}
};

TEST_CASE("number_of_rows") {
	TestSetPartitioningProblem problem(10, 20);
	CHECK(problem.number_of_rows() == 30);
}

TEST_CASE("fractional_solution") {
	int num_groups = 10;
	TestSetPartitioningProblem problem(num_groups, 20);

	auto add_column = [&](vector<int> constraints) {
		Column column(0, 0, 1);
		column.add_coefficient(0, 1.0);
		for (auto c : constraints) {
			column.add_coefficient(num_groups + c, 1.0);
		}
		problem.add_column(move(column));
	};
	add_column({11});
	add_column({11, 12, 13});
	add_column({11, 12, 13, 14});

	problem.pool.at(0).solution_value = 0.5;
	problem.pool.at(1).solution_value = 0.3;
	problem.pool.at(2).solution_value = 0.2;

	problem.compute_fractional_solution({0, 1, 2});

	CHECK(problem.fractional_solution(0, 11) == Approx(1));
	CHECK(problem.fractional_solution(0, 12) == Approx(0.5));
	CHECK(problem.fractional_solution(0, 13) == Approx(0.5));
	CHECK(problem.fractional_solution(0, 14) == Approx(0.2));
}

TEST_CASE("numerical_error") {
	int num_groups = 10;
	TestSetPartitioningProblem problem(num_groups, 20);

	auto add_column = [&](vector<int> constraints) {
		Column column(0, 0, 1);
		column.add_coefficient(0, 1.0);
		for (auto c : constraints) {
			column.add_coefficient(num_groups + c, 1.0);
		}
		problem.add_column(move(column));
	};
	add_column({11});
	add_column({11, 12});

	problem.pool.at(0).solution_value = 0;
	problem.pool.at(1).solution_value = 0;

	CHECK_THROWS_AS(problem.compute_fractional_solution(0, {0, 1}), std::logic_error);
}

// The normal case. If a constraint is assigned to a member a lot, that member
// will be fixed to that constraint.
TEST_CASE("fix") {
	FLAGS_objective_change_before_fixing = 1.0;
	FLAGS_fix_threshold = 0.75;

	int num_groups = 1;
	TestSetPartitioningProblem problem(num_groups, 20);

	auto add_column = [&](vector<int> constraints) {
		Column column(0, 0, 1);
		column.add_coefficient(0, 1.0);
		for (auto c : constraints) {
			column.add_coefficient(num_groups + c, 1.0);
		}
		problem.add_column(move(column));
	};
	add_column({11});
	add_column({11, 12, 13});
	add_column({12, 13, 14});

	problem.pool.at(0).solution_value = 0.5;
	problem.pool.at(1).solution_value = 0.3;
	problem.pool.at(2).solution_value = 0.2;

	vector<size_t> all_columns = {0, 1, 2};

	Problem::FixInformation information;
	information.iteration = 1;
	information.objective_change = -0.1;
	// Do not fix first iteration.
	CHECK(problem.fix_using_columns(information, all_columns) == 0);

	information.iteration = 100;
	information.objective_change = -10;
	// Do not fix if there is a large change.
	CHECK(problem.fix_using_columns(information, all_columns) == 0);

	information.iteration = 100;
	information.objective_change = -0.1;
	CHECK(problem.fix_using_columns(information, all_columns) == 1);

	CHECK(problem.column_allowed(0));
	CHECK_FALSE(problem.pool.at(0).is_fixed());
	CHECK(problem.pool.at(0).solution_value == 0.5 / (0.3 + 0.5));

	CHECK(problem.column_allowed(1));
	CHECK_FALSE(problem.pool.at(1).is_fixed());
	CHECK(problem.pool.at(1).solution_value == 0.3 / (0.3 + 0.5));

	CHECK_FALSE(problem.column_allowed(2));  // This column does not contain 11.
	CHECK(problem.pool.at(2).is_fixed());
	CHECK(problem.pool.at(2).upper_bound() == 0);
	CHECK(problem.pool.at(2).solution_value == 0);

	CHECK(!problem.member_fully_fixed(0));

	vector<int> expected_fixes(20, -1);
	expected_fixes[11] = 1;
	CHECK(problem.fixes_for_member(0) == expected_fixes);

	problem.unfix_all();
	CHECK(problem.fixes_for_member(0) == (vector<int>(20, -1)));
	CHECK(problem.column_allowed(0));
	CHECK(problem.column_allowed(1));
	CHECK(problem.column_allowed(2));
}

TEST_CASE("fix_to_zero") {
	FLAGS_objective_change_before_fixing = 1.0;
	FLAGS_fix_threshold = 0.75;

	int num_groups = 1;
	TestSetPartitioningProblem problem(num_groups, 20);

	auto add_column = [&](vector<int> constraints) {
		Column column(0, 0, 1);
		column.add_coefficient(0, 1.0);
		for (auto c : constraints) {
			column.add_coefficient(num_groups + c, 1.0);
		}
		problem.add_column(move(column));
	};
	add_column({11});
	add_column({11, 12, 13});
	add_column({12, 13, 14});
	add_column({1, 2, 3, 11});

	problem.pool.at(0).solution_value = 0.5;
	problem.pool.at(1).solution_value = 0.3;
	problem.pool.at(2).solution_value = 0.2;
	problem.pool.at(3).solution_value = 0.0;

	vector<size_t> all_columns = {0, 1, 2, 3};

	Problem::FixInformation information;
	information.iteration = 1;
	information.objective_change = -0.1;

	// Doing this 100 times will converge the average fractional
	// solution. Many columns will have an average close to 0.
	for (int i = 0; i < 100; ++i) {
		CHECK(problem.fix_using_columns(information, all_columns) == 0);
	}

	information.iteration = 100;
	information.objective_change = -0.1;
	CHECK(problem.fix_using_columns(information, all_columns) == 1);
	CHECK(problem.column_allowed(0));
	CHECK_FALSE(problem.column_allowed(2));  // Because does not contain 11.
	CHECK_FALSE(problem.column_allowed(3));  // Because of fixes to 0.

	// We expect most columns fixed to zero now.
	vector<int> expected_fixes(20, 0);
	expected_fixes[11] = 1;
	expected_fixes[12] = -1;
	expected_fixes[13] = -1;
	expected_fixes[14] = -1;
	CHECK(problem.fixes_for_member(0) == expected_fixes);

	problem.unfix_all();
	CHECK(problem.column_allowed(0));
	CHECK(problem.column_allowed(2));
	CHECK(problem.column_allowed(3));

	// This time there will be no fixes to 0, since the average fractional
	// solution is not converged.
	CHECK(problem.fix_using_columns(information, all_columns) == 1);
	CHECK_FALSE(problem.column_allowed(2));  // Because does not contain 11.
	CHECK(problem.column_allowed(3));        // No fixes to zero yet.
}

// If no member could be fixed to a constraint, a column with value >= 0.5 is
// fixed to one.
TEST_CASE("fix_column") {
	int num_groups = 1;
	TestSetPartitioningProblem problem(num_groups, 20);

	auto add_column = [&](vector<int> constraints) {
		Column column(0, 0, 1);
		column.add_coefficient(0, 1.0);
		for (auto c : constraints) {
			column.add_coefficient(num_groups + c, 1.0);
		}
		problem.add_column(move(column));
	};
	add_column({11});
	add_column({12, 13});
	add_column({12, 13, 14});

	problem.pool.at(0).solution_value = 0.5;
	problem.pool.at(1).solution_value = 0.25;
	problem.pool.at(2).solution_value = 0.25;

	vector<size_t> all_columns = {0, 1, 2};

	Problem::FixInformation information;
	information.iteration = 100;
	information.objective_change = -0.1;
	CHECK(problem.fix_using_columns(information, all_columns) == 1);

	CHECK(problem.pool.at(0).is_fixed());
	CHECK(problem.pool.at(0).lower_bound() == 1);
	CHECK(problem.pool.at(0).solution_value == 1);

	CHECK(problem.member_fully_fixed(0));

	// This is a column fix which is not reflected in the
	// fixes_for_member call.
	CHECK(problem.fixes_for_member(0) == (vector<int>(20, -1)));

	problem.unfix_all();
	CHECK(!problem.pool.at(0).is_fixed());
}

// The last option is to just fix a column to one.
TEST_CASE("fix_column_as_last_resort") {
	int num_groups = 1;
	TestSetPartitioningProblem problem(num_groups, 20);

	auto add_column = [&](vector<int> constraints) {
		Column column(0, 0, 1);
		column.add_coefficient(0, 1.0);
		for (auto c : constraints) {
			column.add_coefficient(num_groups + c, 1.0);
		}
		problem.add_column(move(column));
	};
	add_column({11});
	add_column({11, 12, 13});
	add_column({12, 13, 14});
	add_column({14, 15});

	problem.pool.at(0).solution_value = 0.25;
	problem.pool.at(1).solution_value = 0.25;
	problem.pool.at(2).solution_value = 0.25;
	problem.pool.at(3).solution_value = 0.25;

	vector<size_t> all_columns = {0, 1, 2, 3};

	Problem::FixInformation information;
	information.iteration = 100;
	information.objective_change = -0.1;
	CHECK(problem.fix_using_columns(information, all_columns) == 4);

	CHECK(problem.pool.at(0).is_fixed());
	CHECK(problem.pool.at(0).lower_bound() == 1);
	CHECK(problem.pool.at(0).solution_value == 1);

	CHECK(problem.pool.at(1).solution_value == 0);
	CHECK(problem.pool.at(2).solution_value == 0);
	CHECK(problem.pool.at(3).solution_value == 0);

	CHECK(problem.member_fully_fixed(0));

	// This is a column fix which is not reflected in the
	// fixes_for_member call.
	CHECK(problem.fixes_for_member(0) == (vector<int>(20, -1)));

	problem.unfix_all();
	CHECK(!problem.pool.at(0).is_fixed());
}

// Of nothing can be fixed, return -1.
TEST_CASE("fix_done") {
	TestSetPartitioningProblem problem(1, 20);

	Column column(0, 0, 1);
	column.add_coefficient(0, 1.0);
	column.solution_value = 1.0;
	problem.add_column(move(column));

	Problem::FixInformation information;
	information.iteration = 100;
	information.objective_change = -0.1;
	CHECK(problem.fix_using_columns(information, {0}) == -1);
}

TEST_CASE("invalid_columns") {
	TestSetPartitioningProblem problem(10, 20);

	CHECK_THROWS(problem.add_column(Column(0, 0, 1)));

	{
		Column column(0, 1, 2);
		column.add_coefficient(0, 1.0);
		CHECK_THROWS(problem.add_column(move(column)));
	}

	{
		Column column(0, 0, 1);
		column.add_coefficient(15, 1.0);
		CHECK_THROWS(problem.add_column(move(column)));
	}

	{
		Column column(0, 0, 1);
		column.add_coefficient(0, 1.0);
		CHECK_NOTHROW(problem.add_column(move(column)));
	}
}

TEST_CASE("make_solution_integer") {
	int num_groups = 2;
	TestSetPartitioningProblem problem(num_groups, 20);

	auto add_column = [&](int p, vector<int> constraints, double value) {
		Column column(0, 0, 1);
		column.add_coefficient(p, 1.0);
		for (auto c : constraints) {
			column.add_coefficient(num_groups + c, 1.0);
		}
		column.solution_value = value;
		problem.add_column(move(column));
	};
	add_column(0, {11}, 0.9);
	add_column(0, {12}, 0.1);
	add_column(1, {14, 15}, 0.3);
	add_column(1, {11, 12, 13}, 0.4);
	add_column(1, {14, 15, 16}, 0.3);

	vector<size_t> indices(problem.pool.size());
	iota(indices.begin(), indices.end(), 0);
	problem.integral_solution_value(indices);

	CHECK(problem.fractional_solution(0, 11) == Approx(1));
	CHECK(problem.fractional_solution(0, 12) == Approx(0));

	CHECK(problem.fractional_solution(1, 11) == Approx(1));
	CHECK(problem.fractional_solution(1, 12) == Approx(1));
	CHECK(problem.fractional_solution(1, 13) == Approx(1));
	CHECK(problem.fractional_solution(1, 14) == Approx(0));
	CHECK(problem.fractional_solution(1, 15) == Approx(0));
	CHECK(problem.fractional_solution(1, 16) == Approx(0));

	// Make sure the pool was not modified.
	CHECK(problem.pool.at(0).solution_value == 0.9);
	CHECK(problem.pool.at(1).solution_value == 0.1);
	CHECK(problem.pool.at(2).solution_value == 0.3);
	CHECK(problem.pool.at(3).solution_value == 0.4);
	CHECK(problem.pool.at(4).solution_value == 0.3);
}

TEST_CASE("integral_solution_value") {
	int num_groups = 2;
	TestSetPartitioningProblem problem(num_groups, 20);

	problem.initialize_constraint(11, 1, 1, 100.0, 1000.0);
	problem.initialize_constraint(12, 1, 1, 100.0, 1000.0);
	problem.initialize_constraint(13, 1, 1, 100.0, 1000.0);
	problem.initialize_constraint(14, 1, 1, 100.0, 1000.0);
	problem.initialize_constraint(15, 1, 1, 100.0, 1000.0);
	problem.initialize_constraint(16, 1, 1, 100.0, 1000.0);

	auto add_column = [&](int p, vector<int> constraints, double value) {
		Column column(1.0, 0, 1);
		column.add_coefficient(p, 1.0);
		for (auto c : constraints) {
			column.add_coefficient(num_groups + c, 1.0);
		}
		column.solution_value = value;
		problem.add_column(move(column));
	};
	add_column(0, {11}, 0.9);
	add_column(0, {12}, 0.1);
	add_column(1, {14, 15}, 0.3);
	add_column(1, {11, 12, 13}, 0.4);
	add_column(1, {14, 15, 16}, 0.3);

	vector<size_t> indices(problem.pool.size());
	iota(indices.begin(), indices.end(), 0);
	auto cost = problem.integral_solution_value(indices);

	CHECK(problem.fractional_solution(0, 11) == Approx(1));
	CHECK(problem.fractional_solution(0, 12) == Approx(0));

	CHECK(problem.fractional_solution(1, 11) == Approx(1));
	CHECK(problem.fractional_solution(1, 12) == Approx(1));
	CHECK(problem.fractional_solution(1, 13) == Approx(1));
	CHECK(problem.fractional_solution(1, 14) == Approx(0));
	CHECK(problem.fractional_solution(1, 15) == Approx(0));
	CHECK(problem.fractional_solution(1, 16) == Approx(0));

	double expected_cost = 1.0 + 1.0  // Column costs.
	                       + 100.0    // Constraint 11.
	                       + 0        // Constraint 12.
	                       + 0        // Constraint 13.
	                       + 1000.0   // Constraint 14.
	                       + 1000.0   // Constraint 15.
	                       + 1000.0;  // Constraint 16.
	CHECK(cost == Approx(expected_cost));
}
