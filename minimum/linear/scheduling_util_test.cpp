#include <fstream>
#include <sstream>
using namespace std;

#include <catch.hpp>

#include <minimum/core/grid.h>
#include <minimum/core/range.h>
#include <minimum/linear/data/util.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/scheduling_util.h>
#include <minimum/linear/solver.h>
using minimum::linear::IP;
using minimum::linear::IPSolver;
using namespace minimum::core;
using namespace minimum::linear;

TEST_CASE("instance1") {
	auto data = R"(#
# This is a comment. Comments start with #
SECTION_HORIZON
# All instances start on a Monday
# The horizon length in days:
14

SECTION_SHIFTS
# ShiftID, Length in mins, Shifts which cannot follow this shift | separated
D,480,

SECTION_STAFF
# ID, MaxShifts, MaxTotalMinutes, MinTotalMinutes, MaxConsecutiveShifts, MinConsecutiveShifts, MinConsecutiveDaysOff, MaxWeekends
A,D=14,4320,3360,5,2,2,1
B,D=14,4320,3360,5,2,2,1
C,D=14,4320,3360,5,2,2,1
D,D=14,4320,3360,5,2,2,1
E,D=14,4320,3360,5,2,2,1
F,D=14,4320,3360,5,2,2,1
G,D=14,4320,3360,5,2,2,1
H,D=14,4320,3360,5,2,2,1

SECTION_DAYS_OFF
# EmployeeID, DayIndexes (start at zero)
A,0
B,5
C,8
D,2
E,9
F,5
G,1
H,7

SECTION_SHIFT_ON_REQUESTS
# EmployeeID, Day, ShiftID, Weight
A,2,D,2
A,3,D,2
B,0,D,3
B,1,D,3
B,2,D,3
B,3,D,3
B,4,D,3
C,0,D,1
C,1,D,1
C,2,D,1
C,3,D,1
C,4,D,1
D,8,D,2
D,9,D,2
F,0,D,2
F,1,D,2
H,9,D,1
H,10,D,1
H,11,D,1
H,12,D,1
H,13,D,1

SECTION_SHIFT_OFF_REQUESTS
# EmployeeID, Day, ShiftID, Weight
C,12,D,1
C,13,D,1
F,8,D,3
H,2,D,3
H,3,D,3

SECTION_COVER
# Day, ShiftID, Requirement, Weight for under, Weight for over
0,D,5,100,1
1,D,7,100,1
2,D,6,100,1
3,D,4,100,1
4,D,5,100,1
5,D,5,100,1
6,D,5,100,1
7,D,6,100,1
8,D,7,100,1
9,D,4,100,1
10,D,2,100,1
11,D,5,100,1
12,D,6,100,1
13,D,4,100,1
)";
	stringstream sin(data);
	auto problem = read_nottingham_instance(sin);
	INFO(problem.DebugString());

	REQUIRE(problem.num_days() == 14);
	REQUIRE(problem.shift_size() == 1);
	REQUIRE(problem.worker_size() == 8);
	REQUIRE(problem.requirement_size() == 14);

	// Make sure that the time number is reduced with gcd.
	CHECK(problem.shift(0).duration() == 1);
	CHECK(problem.worker(0).time_limit().max() == 4320 / 480);
	CHECK(problem.worker(0).time_limit().min() == 3360 / 480);

	// Create and solve an IP.
	IP ip;
	auto x = ip.add_boolean_cube(problem.worker_size(), problem.num_days(), problem.shift_size());
	auto objective = create_ip(ip, problem, x);
	auto solution = make_grid<int>(problem.worker_size(), problem.num_days(), problem.shift_size());

	// Test disabled because it takes more than a second to run.
	if (false) {
		IPSolver solver;
		REQUIRE(solver.solutions(&ip)->get());

		for (int p : range(problem.worker_size())) {
			for (int d : range(problem.num_days())) {
				for (int s : range(problem.shift_size())) {
					solution[p][d][s] = x[p][d][s].value() > 0.5 ? 1 : 0;
				}
			}
		}

		stringstream sol;
		print_solution(sol, problem, solution);
		INFO(sol.str());

		CHECK(objective.value() == 607);
	}
}

TEST_CASE("gent") {
	auto dir = data::get_directory();
	ifstream problem_file(dir + "/NSPLib/N25/1.nsp");
	ifstream case_file(dir + "/NSPLib/Cases/6.gen");
	auto problem = read_gent_instance(problem_file, case_file);
	INFO(problem.DebugString());

	REQUIRE(problem.num_days() == 7);
	REQUIRE(problem.shift_size() == 3);
	REQUIRE(problem.worker_size() == 25);
	REQUIRE(problem.requirement_size() == 7 * 3);

	IP ip;
	auto x = ip.add_boolean_cube(problem.worker_size(), problem.num_days(), problem.shift_size());
	auto objective = create_ip(ip, problem, x);
	auto solution = make_grid<int>(problem.worker_size(), problem.num_days(), problem.shift_size());

	// Test disabled because it takes more than a second to run.
	IPSolver solver;
	REQUIRE(solver.solutions(&ip)->get());

	for (int p : range(problem.worker_size())) {
		for (int d : range(problem.num_days())) {
			for (int s : range(problem.shift_size())) {
				solution[p][d][s] = x[p][d][s].value() > 0.5 ? 1 : 0;
			}
		}
	}

	stringstream sol;
	print_solution(sol, problem, solution);
	INFO(sol.str());
	CHECK(objective.value() == 301);
}

TEST_CASE("verify_scheduling_problem") {
	minimum::linear::proto::SchedulingProblem empty_problem;
	CHECK_THROWS(verify_scheduling_problem(empty_problem));

	auto dir = data::get_directory();
	ifstream problem_file(dir + "/NSPLib/N25/1.nsp");
	ifstream case_file(dir + "/NSPLib/Cases/6.gen");
	auto problem = read_gent_instance(problem_file, case_file);

	CHECK_NOTHROW(verify_scheduling_problem(problem));

	{
		auto p = problem;
		p.clear_num_days();
		CHECK_THROWS(verify_scheduling_problem(p));
	}

	{
		auto p = problem;
		p.clear_shift();
		CHECK_THROWS(verify_scheduling_problem(p));
	}

	{
		auto p = problem;
		p.clear_worker();
		CHECK_THROWS(verify_scheduling_problem(p));
	}

	{
		auto p = problem;
		p.clear_requirement();
		CHECK_THROWS(verify_scheduling_problem(p));
	}

	{
		auto p = problem;
		p.add_shift();
		CHECK_THROWS(verify_scheduling_problem(p));
	}

	{
		auto p = problem;
		p.mutable_worker(1)->add_shift_limit();
		CHECK_THROWS(verify_scheduling_problem(p));
	}
}

TEST_CASE("save_load_solution") {
	auto dir = data::get_directory();
	ifstream problem_file(dir + "/NSPLib/N25/1.nsp");
	ifstream case_file(dir + "/NSPLib/Cases/6.gen");
	auto problem = read_gent_instance(problem_file, case_file);

	auto sol = make_grid<int>(problem.worker_size(), problem.num_days(), problem.shift_size());
	for (int p : range(problem.worker_size())) {
		for (int d : range(problem.num_days())) {
			if (d % 2 == 0) {
				sol[p][d][(p + d) % problem.shift_size()] = 1;
			}
		}
	}

	stringstream data;
	save_solution(data, problem, sol);
	auto read = load_solution(data, problem);

	CHECK(sol == read);
}
