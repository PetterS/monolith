#include <vector>
using namespace std;

#include <catch.hpp>

#include <minimum/core/grid.h>
#include <minimum/core/range.h>
#include <minimum/linear/colgen/shift_scheduling_pricing.h>
using namespace minimum::core;
using namespace minimum::linear::colgen;

// Creates a basic problem with one shift and one worker.
// The worker has time_limit and consecutive_shifts_limit
// set, but in a way so that they do not have an effect.
minimum::linear::proto::SchedulingProblem basic_problem() {
	minimum::linear::proto::SchedulingProblem problem;
	problem.set_num_days(14);

	auto shift = problem.add_shift();
	shift->set_duration(1);
	shift->set_id("S");
	shift->add_ok_after(0);

	auto worker = problem.add_worker();
	worker->set_id("A");
	auto shift_limit = worker->add_shift_limit();
	shift_limit->set_min(0);
	shift_limit->set_max(14);
	auto time_limit = worker->mutable_time_limit();
	time_limit->set_max(problem.num_days());
	auto consecutive_limit = worker->mutable_consecutive_shifts_limit();
	consecutive_limit->set_max(problem.num_days());

	for (int d : range(problem.num_days())) {
		auto req = problem.add_requirement();
		req->set_day(d);
		req->set_shift(0);
	}
	return problem;
}

std::mt19937_64 rng;

TEST_CASE("basic") {
	auto problem = basic_problem();
	auto fixes = make_grid<int>(problem.num_days(), 1, []() { return -1; });
	auto solution = make_grid<int>(problem.num_days(), 1);
	vector<double> duals(problem.worker_size() + problem.num_days(), 0);
	REQUIRE(create_roster_cspp(problem, duals, 0, fixes, &solution, &rng));
}

TEST_CASE("time_limit") {
	auto problem = basic_problem();
	auto fixes = make_grid<int>(problem.num_days(), 1, []() { return -1; });
	auto solution = make_grid<int>(problem.num_days(), 1);
	vector<double> duals(problem.worker_size() + problem.num_days(), 0);

	duals[1 + 3] = 1;
	duals[1 + 4] = 1;
	duals[1 + 5] = 2;
	duals[1 + 6] = 3;
	duals[1 + 7] = 2;
	duals[1 + 8] = 2;
	duals[1 + 9] = 1;
	auto limit = problem.mutable_worker(0)->mutable_time_limit();
	limit->set_min(4);
	limit->set_max(4);
	REQUIRE(create_roster_cspp(problem, duals, 0, fixes, &solution, &rng));
	int assigned = 0;
	for (int d : range(problem.num_days())) {
		assigned += solution[d][0];
	}
	CHECK(assigned == 4);
	CHECK(solution[5][0] == 1);
	CHECK(solution[6][0] == 1);
	CHECK(solution[7][0] == 1);
	CHECK(solution[8][0] == 1);
}

TEST_CASE("consequtive_limit") {
	auto problem = basic_problem();
	auto fixes = make_grid<int>(problem.num_days(), 1, []() { return -1; });
	auto solution = make_grid<int>(problem.num_days(), 1);
	vector<double> duals(problem.worker_size() + problem.num_days(), 0);

	duals[1 + 3] = -2;
	duals[1 + 4] = -1;
	duals[1 + 5] = 100;
	duals[1 + 6] = 100;
	duals[1 + 7] = -1;
	duals[1 + 8] = -2;

	REQUIRE(create_roster_cspp(problem, duals, 0, fixes, &solution, &rng));
	CHECK(solution[3][0] == 0);
	CHECK(solution[4][0] == 0);
	CHECK(solution[5][0] == 1);
	CHECK(solution[6][0] == 1);
	CHECK(solution[7][0] == 0);
	CHECK(solution[8][0] == 0);

	auto limit = problem.mutable_worker(0)->mutable_consecutive_shifts_limit();
	limit->set_min(4);
	limit->set_max(10000);
	REQUIRE(create_roster_cspp(problem, duals, 0, fixes, &solution, &rng));
	CHECK(solution[3][0] == 0);
	CHECK(solution[4][0] == 1);
	CHECK(solution[5][0] == 1);
	CHECK(solution[6][0] == 1);
	CHECK(solution[7][0] == 1);
	CHECK(solution[8][0] == 0);
}

TEST_CASE("consequtive_days_off") {
	auto problem = basic_problem();
	auto fixes = make_grid<int>(problem.num_days(), 1, []() { return -1; });
	auto solution = make_grid<int>(problem.num_days(), 1);
	vector<double> duals(problem.worker_size() + problem.num_days(), 0);

	duals[1 + 3] = 2;
	duals[1 + 4] = 2;
	duals[1 + 5] = 2;
	duals[1 + 6] = -1;
	duals[1 + 7] = 2;
	duals[1 + 8] = 2;
	duals[1 + 9] = 2;

	REQUIRE(create_roster_cspp(problem, duals, 0, fixes, &solution, &rng));
	CHECK(solution[6][0] == 0);

	auto limit = problem.mutable_worker(0)->mutable_consecutive_days_off_limit();
	limit->set_min(2);
	limit->set_max(10000);
	REQUIRE(create_roster_cspp(problem, duals, 0, fixes, &solution, &rng));
	CHECK(solution[6][0] == 1);
}

TEST_CASE("max_weekends") {
	auto problem = basic_problem();
	auto fixes = make_grid<int>(problem.num_days(), 1, []() { return -1; });
	auto solution = make_grid<int>(problem.num_days(), 1);
	vector<double> duals(problem.worker_size() + problem.num_days(), 1);

	REQUIRE(create_roster_cspp(problem, duals, 0, fixes, &solution, &rng));
	CHECK(solution[5][0] == 1);      // Saturday
	CHECK(solution[7 + 5][0] == 1);  // Saturday

	auto limit = problem.mutable_worker(0)->mutable_working_weekends_limit();
	limit->set_max(1);
	REQUIRE(create_roster_cspp(problem, duals, 0, fixes, &solution, &rng));
	CHECK((solution[5][0] == 0 || solution[7 + 5][0] == 0));  // Saturdays
	CHECK(solution[5][0] == solution[6][0]);
	CHECK(solution[7 + 5][0] == solution[7 + 6][0]);  // Sundays
}

TEST_CASE("shift_limit") {
	auto problem = basic_problem();
	auto fixes = make_grid<int>(problem.num_days(), 1, []() { return -1; });
	auto solution = make_grid<int>(problem.num_days(), 1);
	vector<double> duals(problem.worker_size() + problem.num_days(), 1);

	REQUIRE(create_roster_cspp(problem, duals, 0, fixes, &solution, &rng));
	int assigned = 0;
	for (int d : range(problem.num_days())) {
		assigned += solution[d][0];
	}
	CHECK(assigned == problem.num_days());

	auto limit = problem.mutable_worker(0)->mutable_shift_limit(0);
	limit->set_min(0);
	limit->set_max(5);
	REQUIRE(create_roster_cspp(problem, duals, 0, fixes, &solution, &rng));
	assigned = 0;
	for (int d : range(problem.num_days())) {
		assigned += solution[d][0];
	}
	CHECK(assigned == 5);
}
