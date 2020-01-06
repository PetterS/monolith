#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
using namespace std;

#include <minimum/core/check.h>
#include <minimum/core/grid.h>
#include <minimum/core/range.h>
#include <minimum/core/string.h>
#include <minimum/linear/glpk.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/scheduling_util.h>
#include <minimum/linear/solver.h>
using namespace minimum::core;
using namespace minimum::linear;

namespace {
string read_line(istream& file) {
	string line;
	getline(file, line);
	if (!file) {
		return "";
	}

	// Check for comment.
	for (auto ch : line) {
		if (ch == '#') {
			return read_line(file);
		} else if (ch != ' ') {
			break;
		}
	}

	if (line.size() > 0 && line[line.size() - 1] == '\r') {
		line = line.substr(0, line.size() - 1);
	}

	return line;
}

template <typename Int>
Int gcd(Int a, Int b) {
	minimum_core_assert(b >= 0);
	if (a < b) {
		return gcd(b, a);
	}
	if (a == 1 || b == 1) {
		return 1;
	} else if (b == 0) {
		return a;
	} else {
		return gcd(b, a - b);
	}
}

string remove_space(string in) {
	in.erase(remove_if(in.begin(), in.end(), ::isspace), in.end());
	return in;
}
}  // namespace

namespace minimum {
namespace linear {

MINIMUM_LINEAR_API void verify_scheduling_problem(
    const minimum::linear::proto::SchedulingProblem& problem) {
	minimum_core_assert(problem.worker_size() > 0);
	minimum_core_assert(problem.shift_size() > 0);
	minimum_core_assert(problem.requirement_size() > 0);
	minimum_core_assert(problem.num_days() > 0);

	auto verify_bound = [](const minimum::linear::proto::IntegerBound& bound) {
		minimum_core_assert(bound.min() <= bound.max());
	};
	auto verify_day = [&problem](int d) { minimum_core_assert(0 <= d && d < problem.num_days()); };
	auto verify_shift = [&problem](int s) {
		minimum_core_assert(0 <= s && s < problem.shift_size());
	};

	for (auto& worker : problem.worker()) {
		if (worker.has_time_limit()) {
			verify_bound(worker.time_limit());
		}
		if (worker.has_consecutive_shifts_limit()) {
			verify_bound(worker.consecutive_shifts_limit());
		}
		if (worker.has_consecutive_days_off_limit()) {
			verify_bound(worker.consecutive_days_off_limit());
		}
		if (worker.has_working_weekends_limit()) {
			verify_bound(worker.working_weekends_limit());
		}

		minimum_core_assert(worker.shift_limit_size() == problem.shift_size());

		for (int day : worker.required_day_off()) {
			verify_day(day);
		}

		for (auto& pref : worker.shift_preference()) {
			verify_day(pref.day());
			verify_shift(pref.shift());
		}

		for (auto& pref : worker.day_off_preference()) {
			verify_day(pref.day());
		}
	}

	for (auto& shift : problem.shift()) {
		minimum_core_assert(shift.duration() >= 0);
		set<int> other_shifts;
		for (int s : shift.ok_after()) {
			verify_shift(s);
			other_shifts.insert(s);
		}
		for (int s : shift.forbidden_after()) {
			verify_shift(s);
			other_shifts.insert(s);
		}
		minimum_core_assert(other_shifts.size() == problem.shift_size());
		if (shift.has_consecutive_limit()) {
			verify_bound(shift.consecutive_limit());
		}
	}

	for (auto& requirement : problem.requirement()) {
		verify_day(requirement.day());
		verify_shift(requirement.shift());
	}
}

minimum::linear::proto::SchedulingProblem read_nottingham_instance(std::istream& in, bool verbose) {
	check(bool(in), "File not found or invalid stream.");

	minimum::linear::proto::SchedulingProblem problem;
	map<string, int> id_to_shift;
	map<string, int> id_to_staff;

	minimum_core_assert(read_line(in) == "SECTION_HORIZON");
	problem.set_num_days(from_string<int>(read_line(in)));
	minimum_core_assert(problem.num_days() > 0, "Need days.");
	if (verbose) {
		clog << "Number of days: " << problem.num_days() << endl;
	}
	minimum_core_assert(read_line(in) == "");

	minimum_core_assert(read_line(in) == "SECTION_SHIFTS");
	auto shift_string = read_line(in);
	vector<pair<int, string>> forbidden_after_pairs;
	while (shift_string != "") {
		auto parts = split(shift_string, ',');
		auto shift = problem.add_shift();

		int e = 0;
		shift->set_id(parts.at(e++));
		shift->set_duration(from_string<int>(parts.at(e++)));
		if (e < parts.size()) {
			forbidden_after_pairs.emplace_back(problem.shift_size() - 1, parts.at(e++));
		}
		minimum_core_assert(e == parts.size());

		if (verbose) {
			clog << "-- Shift " << shift->id() << ", " << shift->duration() << " min.\n";
		}

		id_to_shift[shift->id()] = problem.shift_size() - 1;

		shift_string = read_line(in);
	}
	minimum_core_assert(problem.shift_size() > 0, "Need shifts.");
	if (verbose) {
		clog << "Number of shifts: " << problem.shift_size() << endl;
	}

	// We need to parse the forbidden combinations afterwards,
	// when we have seen all shift IDs.
	for (auto& ix_string : forbidden_after_pairs) {
		auto forbidden_after = split(ix_string.second, '|');
		for (auto& forbidden : forbidden_after) {
			auto s = id_to_shift.at(forbidden);
			problem.mutable_shift(ix_string.first)->add_forbidden_after(s);
		}
	}
	for (int s1 = 0; s1 < problem.shift_size(); ++s1) {
		for (int s2 = 0; s2 < problem.shift_size(); ++s2) {
			if (find(problem.shift(s1).forbidden_after().begin(),
			         problem.shift(s1).forbidden_after().end(),
			         s2)
			    == problem.shift(s1).forbidden_after().end()) {
				problem.mutable_shift(s1)->add_ok_after(s2);
			}
		}

		if (verbose) {
			clog << "-- Shift " << s1 << ". OK after = " << to_string(problem.shift(s1).ok_after())
			     << " forbidden_after = " << to_string(problem.shift(s1).forbidden_after()) << "\n";
		}
	}

	minimum_core_assert(read_line(in) == "SECTION_STAFF");
	auto staff_string = read_line(in);
	while (staff_string != "") {
		auto parts = split(staff_string, ',');
		auto worker = problem.add_worker();

		int e = 0;
		worker->set_id(parts.at(e++));

		for (int s : range(problem.shift_size())) {
			auto limit = worker->add_shift_limit();
			limit->set_min(0);
			limit->set_max(problem.num_days());
		}

		vector<int> print_max_shifts;
		auto max_shifts = split(parts.at(e++), '|');
		for (auto& max_shift : max_shifts) {
			auto id_value = split(max_shift, '=');
			auto s = id_to_shift[id_value.at(0)];
			auto shift_limit = worker->mutable_shift_limit(s);
			shift_limit->set_max(from_string<int>(id_value.at(1)));
			print_max_shifts.push_back(from_string<int>(id_value.at(1)));
		}

		worker->mutable_time_limit()->set_max(from_string<int>(parts.at(e++)));
		worker->mutable_time_limit()->set_min(from_string<int>(parts.at(e++)));
		worker->mutable_consecutive_shifts_limit()->set_max(from_string<int>(parts.at(e++)));
		worker->mutable_consecutive_shifts_limit()->set_min(from_string<int>(parts.at(e++)));
		worker->mutable_consecutive_days_off_limit()->set_max(problem.num_days());  // Not used.
		worker->mutable_consecutive_days_off_limit()->set_min(from_string<int>(parts.at(e++)));
		worker->mutable_working_weekends_limit()->set_max(from_string<int>(parts.at(e++)));
		worker->mutable_working_weekends_limit()->set_min(0);  // Not used.
		minimum_core_assert(e == parts.size());

		if (verbose) {
			clog << "-- Staff " << setw(2) << worker->id() << ": Consequtive limit ["
			     << worker->consecutive_shifts_limit().min() << ", "
			     << worker->consecutive_shifts_limit().max() << "]. "
			     << "Time limit [" << worker->time_limit().min() << ", "
			     << worker->time_limit().max()
			     << "]. Max weekends: " << worker->working_weekends_limit().max()
			     << ". Max shifts: " << to_string(print_max_shifts) << "." << endl;
		}

		id_to_staff[worker->id()] = problem.worker_size() - 1;

		staff_string = read_line(in);
	}
	minimum_core_assert(problem.worker_size() > 0, "Need workers.");
	if (verbose) {
		clog << "Number of staff: " << problem.worker_size() << endl;
	}

	minimum_core_assert(read_line(in) == "SECTION_DAYS_OFF");
	auto days_off_string = read_line(in);
	while (days_off_string != "") {
		auto days_off = split(days_off_string, ',');
		auto w = id_to_staff.at(days_off.at(0));
		for (auto itr = days_off.begin() + 1; itr != days_off.end(); ++itr) {
			problem.mutable_worker(w)->add_required_day_off(from_string<int>(*itr));
		}
		days_off_string = read_line(in);
	}

	minimum_core_assert(read_line(in) == "SECTION_SHIFT_ON_REQUESTS");
	int on_requests_read = 0;
	auto shift_on_request_string = read_line(in);
	while (shift_on_request_string != "") {
		auto shift_on_request = split(shift_on_request_string, ',');
		minimum_core_assert(shift_on_request.size() == 4);
		auto w = id_to_staff.at(shift_on_request.at(0));

		auto request = problem.mutable_worker(w)->add_shift_preference();
		request->set_day(from_string<int>(shift_on_request.at(1)));
		request->set_shift(id_to_shift.at(shift_on_request.at(2)));
		auto weight = from_string<double>(shift_on_request.at(3));
		problem.set_objective_constant(problem.objective_constant() + weight);
		request->set_cost(-weight);

		on_requests_read++;
		shift_on_request_string = read_line(in);
	}
	if (verbose) {
		clog << "Read " << on_requests_read << " shift on requests." << endl;
	}

	minimum_core_assert(read_line(in) == "SECTION_SHIFT_OFF_REQUESTS");
	int off_requests_read = 0;
	auto shift_off_request_string = read_line(in);
	while (shift_off_request_string != "") {
		auto shift_off_request = split(shift_off_request_string, ',');
		minimum_core_assert(shift_off_request.size() == 4);
		auto w = id_to_staff.at(shift_off_request.at(0));

		auto request = problem.mutable_worker(w)->add_shift_preference();
		request->set_day(from_string<int>(shift_off_request.at(1)));
		request->set_shift(id_to_shift.at(shift_off_request.at(2)));
		auto weight = from_string<double>(shift_off_request.at(3));
		request->set_cost(weight);

		off_requests_read++;
		shift_off_request_string = read_line(in);
	}
	if (verbose) {
		clog << "Read " << off_requests_read << " shift off requests." << endl;
	}

	minimum_core_assert(read_line(in) == "SECTION_COVER");
	auto found_constraint =
	    make_grid<bool>(problem.num_days(), problem.shift_size(), []() { return false; });
	auto cover_string = read_line(in);
	while (cover_string != "") {
		auto cover = split(cover_string, ',');
		auto requirement = problem.add_requirement();

		int e = 0;
		auto day = from_string<int>(cover.at(e++));
		auto shift = id_to_shift.at(cover.at(e++));
		requirement->set_day(day);
		requirement->set_shift(shift);
		requirement->set_wanted(from_string<int>(cover.at(e++)));
		requirement->set_under_cost(from_string<int>(cover.at(e++)));
		requirement->set_over_cost(from_string<int>(cover.at(e++)));
		minimum_core_assert(e == cover.size());

		minimum_core_assert(0 <= day && day < problem.num_days());
		minimum_core_assert(0 <= shift && shift < problem.shift_size());
		minimum_core_assert(!found_constraint[day][shift]);
		found_constraint[day][shift] = true;

		cover_string = read_line(in);
	}
	if (verbose) {
		clog << "Read " << problem.requirement_size() << " cover constraints." << endl;
	}

	// Add dummy constraints for all missing cover constraints.
	int dummy_added = 0;
	for (auto d : range(problem.num_days())) {
		for (auto s : range(problem.shift_size())) {
			if (!found_constraint[d][s]) {
				auto requirement = problem.add_requirement();
				requirement->set_day(d);
				requirement->set_shift(s);
				dummy_added++;
			}
		}
	}
	if (dummy_added > 0 && verbose) {
		clog << "Added " << dummy_added << " dummy (empty) cover constraints.\n";
	}

	minimum_core_assert(id_to_staff.size() == problem.worker_size());
	minimum_core_assert(id_to_shift.size() == problem.shift_size());

	// Simplify the problem by reducing the size of the coefficients.
	int divisor = 0;
	for (auto& shift : problem.shift()) {
		divisor = gcd(divisor, shift.duration());
	}
	if (divisor > 1) {
		if (verbose) {
			clog << "Dividing all minutes by " << divisor << endl;
		}
		for (auto& shift : *problem.mutable_shift()) {
			minimum_core_assert(shift.duration() % divisor == 0);
			shift.set_duration(shift.duration() / divisor);
		}
		for (auto& staff : *problem.mutable_worker()) {
			// Round max down.
			staff.mutable_time_limit()->set_max(
			    floor(double(staff.time_limit().max()) / double(divisor)));
			// Round min up.
			staff.mutable_time_limit()->set_min(
			    ceil(double(staff.time_limit().min()) / double(divisor)));
		}
	}

	verify_scheduling_problem(problem);
	return problem;
}

MINIMUM_LINEAR_API minimum::linear::proto::SchedulingProblem read_gent_instance(
    std::istream& problem_data, std::istream& case_data) {
	minimum::linear::proto::SchedulingProblem problem;
	minimum_core_assert(problem_data);

	int num_nurses;
	int num_days;
	int num_shifts;
	problem_data >> num_nurses >> num_days >> num_shifts;
	minimum_core_assert(problem_data);
	problem.set_num_days(num_days);

	// The last shift represents a free day, so do not add
	// that as a shift.
	for (int s = 0; s < num_shifts - 1; ++s) {
		auto shift = problem.add_shift();
		shift->set_id(to_string(s));
		shift->set_duration(1);

		for (int s2 = 0; s2 < num_shifts - 1; ++s2) {
			shift->add_ok_after(s2);
		}
	}

	for (int d = 0; d < num_days; ++d) {
		for (int s = 0; s < num_shifts - 1; ++s) {
			auto requirement = problem.add_requirement();
			int coverage = -1;
			problem_data >> coverage;
			minimum_core_assert(problem_data);

			requirement->set_day(d);
			requirement->set_shift(s);
			requirement->set_wanted(coverage);
			requirement->set_over_cost(0);
			requirement->set_under_cost(10'000);
		}

		int coverage = -1;
		problem_data >> coverage;
		check(coverage == 0, "Can not have coverage for free shift.");
		minimum_core_assert(problem_data);
	}

	for (int n = 0; n < num_nurses; ++n) {
		auto worker = problem.add_worker();
		worker->set_id(to_string('A' + char(n)));

		for (int d = 0; d < num_days; ++d) {
			for (int s = 0; s < num_shifts - 1; ++s) {
				auto pref = worker->add_shift_preference();
				int cost = -1;
				problem_data >> cost;
				minimum_core_assert(problem_data);
				pref->set_day(d);
				pref->set_shift(s);
				pref->set_cost(cost);
			}

			auto pref = worker->add_day_off_preference();
			int cost = -1;
			problem_data >> cost;
			minimum_core_assert(problem_data);
			pref->set_day(d);
			pref->set_cost(cost);
		}
	}

	// Check that there isn’t more information.
	int dummy = 0;
	problem_data >> dummy;
	minimum_core_assert(!problem_data);

	int num_days_tmp = -1;
	int num_shifts_tmp = -1;
	case_data >> num_days_tmp >> num_shifts_tmp;
	minimum_core_assert(case_data);
	check(num_days_tmp == num_days && num_shifts_tmp == num_shifts,
	      "cases_data file does not match problem_data file.");

	int min_assignments;
	int max_assignments;
	case_data >> min_assignments >> max_assignments;
	minimum_core_assert(case_data);
	for (int n = 0; n < num_nurses; ++n) {
		problem.mutable_worker(n)->mutable_time_limit()->set_max(max_assignments);
		problem.mutable_worker(n)->mutable_time_limit()->set_min(min_assignments);
	}

	int min_consecutive;
	int max_consecutive;
	case_data >> min_consecutive >> max_consecutive;
	minimum_core_assert(case_data);
	for (int n = 0; n < num_nurses; ++n) {
		problem.mutable_worker(n)->mutable_consecutive_shifts_limit()->set_max(max_consecutive);
		problem.mutable_worker(n)->mutable_consecutive_shifts_limit()->set_min(min_consecutive);
	}

	int min_consecutive_shifts;
	int max_consecutive_shifts;
	int min_assignments_per_shift;
	int max_assignments_per_shift;
	for (int s = 0; s < num_shifts - 1; ++s) {
		case_data >> min_consecutive_shifts;
		case_data >> max_consecutive_shifts;
		auto consecutive_limit = problem.mutable_shift(s)->mutable_consecutive_limit();
		consecutive_limit->set_min(min_consecutive_shifts);
		consecutive_limit->set_max(max_consecutive_shifts);

		case_data >> min_assignments_per_shift;
		case_data >> max_assignments_per_shift;
		for (int n = 0; n < num_nurses; ++n) {
			auto limit = problem.mutable_worker(n)->add_shift_limit();
			limit->set_max(max_assignments_per_shift);
			limit->set_min(min_assignments_per_shift);
		}
		minimum_core_assert(case_data);
	}
	case_data >> min_consecutive_shifts;
	case_data >> max_consecutive_shifts;
	for (int n = 0; n < num_nurses; ++n) {
		auto limit = problem.mutable_worker(n)->mutable_consecutive_days_off_limit();
		limit->set_min(min_consecutive_shifts);
		limit->set_max(max_consecutive_shifts);
	}
	int min_consecutive_days_off;
	int max_consecutive_days_off;
	case_data >> min_consecutive_days_off;
	case_data >> max_consecutive_days_off;
	for (int n = 0; n < num_nurses; ++n) {
		problem.mutable_worker(n)->mutable_consecutive_days_off_limit()->set_min(
		    min_consecutive_days_off);
		problem.mutable_worker(n)->mutable_consecutive_days_off_limit()->set_max(
		    max_consecutive_days_off);
	}
	minimum_core_assert(case_data);

	// Check that there isn’t more information.
	case_data >> dummy;
	minimum_core_assert(!case_data);

	verify_scheduling_problem(problem);
	return problem;
}

minimum::linear::proto::TaskSchedulingProblem read_ptask_instance(std::istream& data) {
	minimum_core_assert(data);
	auto line = read_line(data);
	auto tokens = split(line, '=');
	minimum_core_assert(remove_space(tokens.at(0)) == "Type");
	minimum_core_assert(remove_space(tokens.at(1)) == "1");
	line = read_line(data);
	tokens = split(line, '=');
	minimum_core_assert(remove_space(tokens.at(0)) == "Jobs");
	int num_tasks = from_string<int>(tokens.at(1));
	minimum_core_assert(num_tasks > 0);
	clog << "-- " << num_tasks << " tasks.\n";

	proto::TaskSchedulingProblem problem;
	for (int t = 0; t < num_tasks; ++t) {
		line = read_line(data);
		minimum_core_assert(data);
		stringstream sin(line);
		int start = -1;
		int end = -1;
		sin >> start >> end;
		minimum_core_assert(sin);

		auto task = problem.add_task();
		task->set_start(start);
		task->set_end(end);
	}
	line = read_line(data);
	tokens = split(line, '=');
	minimum_core_assert(remove_space(tokens.at(0)) == "Qualifications");
	int num_machines = from_string<int>(tokens.at(1));
	minimum_core_assert(num_machines > 0);
	clog << "-- " << num_machines << " machines.\n";

	for (int m = 0; m < num_machines; ++m) {
		line = read_line(data);
		minimum_core_assert(data);
		tokens = split(line, ':');
		int num_tasks_for_machine = from_string<int>(tokens.at(0));

		auto machine = problem.add_machine();

		stringstream sin(tokens.at(1));
		vector<int> tasks_for_machine;
		tasks_for_machine.reserve(num_tasks_for_machine);
		for (int t = 0; t < num_tasks_for_machine; ++t) {
			int task = -1;
			sin >> task;
			minimum_core_assert(sin);
			minimum_core_assert(0 <= task && task < num_tasks);
			tasks_for_machine.push_back(task);
		}

		sort(tasks_for_machine.begin(), tasks_for_machine.end(), [&problem](int t1, int t2) {
			return problem.task(t1).start() < problem.task(t2).start();
		});

		for (auto task : tasks_for_machine) {
			machine->add_qualification(task);
		}
	}

	return problem;
}

void save_solution(ostream& out,
                   const minimum::linear::proto::SchedulingProblem& problem,
                   const std::vector<std::vector<std::vector<int>>>& solution) {
	for (int p = 0; p < problem.worker_size(); ++p) {
		for (int d = 0; d < problem.num_days(); ++d) {
			for (int s = 0; s < problem.shift_size(); ++s) {
				if (solution[p][d][s] == 1) {
					out << problem.shift(s).id();
				}
			}
			out << '\t';
		}
		out << endl;
	}
}

vector<vector<vector<int>>> load_solution(
    istream& in, const minimum::linear::proto::SchedulingProblem& problem) {
	vector<vector<vector<int>>> solution(
	    problem.worker_size(),
	    vector<vector<int>>(problem.num_days(), vector<int>(problem.shift_size(), 0)));
	string line;
	for (int p : range(problem.worker_size())) {
		getline(in, line);
		check(bool(in), "Could not read line from file.");
		auto shifts = split(line, '\t');
		check(shifts.size() == problem.num_days(),
		      "Incorrect number of shifts for ",
		      problem.worker(p).id(),
		      ".");
		for (int d : range(problem.num_days())) {
			auto& shift = shifts[d];
			for (int s : range(problem.shift_size())) {
				if (shift == problem.shift(s).id()) {
					solution[p][d][s] = 1;
				}
			}
		}
	}
	return solution;
}

double print_solution(ostream& out,
                      const minimum::linear::proto::SchedulingProblem& problem,
                      const std::vector<std::vector<std::vector<int>>>& solution) {
	for (int p = 0; p < problem.worker_size(); ++p) {
		out << setw(3) << right << problem.worker(p).id() << ": ";
		for (int d = 0; d < problem.num_days(); ++d) {
			bool printed = false;
			for (int s = 0; s < problem.shift_size(); ++s) {
				if (solution[p][d][s] == 1) {
					minimum_core_assert(!printed,
					                    "Solution has more than one shift assigned on one day.");
					out << to_string(setw(2), setfill('='), left, problem.shift(s).id());
					printed = true;
				}
			}
			if (!printed) {
				out << setw(2) << ".";
			}
		}
		out << endl;
	}
	auto objective = objective_value(problem, solution);
	out << "Final solution (checked for feasibility): " << objective << "\n";
	out << "-- Cover costs  : " << cover_cost(problem, solution) << "\n";
	out << "-- Roster costs : " << objective - cover_cost(problem, solution) << "\n";
	return objective;
}

int cover_cost(const minimum::linear::proto::SchedulingProblem& problem,
               const std::vector<std::vector<std::vector<int>>>& solution) {
	int constraint_cost = 0;
	for (auto& requirement : problem.requirement()) {
		int constraint_sum = 0;
		for (auto p : range(problem.worker_size())) {
			constraint_sum += solution[p][requirement.day()][requirement.shift()];
		}
		if (constraint_sum < requirement.wanted()) {
			constraint_cost += (requirement.wanted() - constraint_sum) * requirement.under_cost();
		} else if (constraint_sum > requirement.wanted()) {
			constraint_cost += (constraint_sum - requirement.wanted()) * requirement.over_cost();
		}
	}
	return constraint_cost;
}

Sum create_ip(IP& ip,
              const minimum::linear::proto::SchedulingProblem& problem,
              const vector<vector<vector<Sum>>>& x,
              const std::vector<int>& staff_indices,
              const std::vector<int>& day_indices) {
	// Returns a variable (Sum) that is 1 iff a person is working
	// on a particular day.
	auto working_on_day = [&](int person, int day) {
		Sum working_on_day = 0;
		for (int s = 0; s < problem.shift_size(); ++s) {
			working_on_day += x[person][day][s];
		}
		return working_on_day;
	};

	try {
		// At most one shift per day per person.
		for (int p : staff_indices) {
			for (int d : day_indices) {
				Sum num_day_shifts = 0;
				for (int s = 0; s < problem.shift_size(); ++s) {
					num_day_shifts += x[p][d][s];
				}
				ip.add_constraint(num_day_shifts <= 1);
			}
		}
	} catch (std::exception& e) {
		std::clog << "-- Infeasible due to at most one shift per day per person. " << e.what()
		          << "\n";
		throw;
	}

	// Maximum shifts.
	for (int s = 0; s < problem.shift_size(); ++s) {
		for (int p : staff_indices) {
			Sum working_shifts = 0;
			for (int d : day_indices) {
				working_shifts += x[p][d][s];
			}
			try {
				ip.add_constraint(working_shifts >= problem.worker(p).shift_limit(s).min());
				ip.add_constraint(working_shifts <= problem.worker(p).shift_limit(s).max());
			} catch (std::exception& e) {
				std::clog << "-- Infeasible due to maximum shifts. " << e.what() << "\n";
				throw;
			}
			try {
				ip.add_constraint(working_shifts >= problem.worker(p).shift_limit(s).min());
			} catch (std::exception& e) {
				std::clog << "-- Infeasible due to minimum shifts. " << e.what() << "\n";
				throw;
			}
		}
	}

	// Max/min total minutes.
	try {
		for (int p : staff_indices) {
			if (!problem.worker(p).has_time_limit()) {
				continue;
			}

			Sum working_minutes = 0;
			for (int d : day_indices) {
				for (int s = 0; s < problem.shift_size(); ++s) {
					working_minutes += problem.shift(s).duration() * x[p][d][s];
				}
			}

			ip.add_constraint(problem.worker(p).time_limit().min(),
			                  working_minutes,
			                  problem.worker(p).time_limit().max());
		}
	} catch (std::exception& e) {
		std::clog << "-- Infeasible due to min/max total minutes. " << e.what() << "\n";
		throw;
	}

	// Shifts that cannot follow other shifts.
	try {
		for (int s1 = 0; s1 < problem.shift_size(); ++s1) {
			for (int s2 : problem.shift(s1).forbidden_after()) {
				for (int p : staff_indices) {
					for (int d = 0; d < problem.num_days() - 1; ++d) {
						ip.add_constraint(x[p][d][s1] + x[p][d + 1][s2] <= 1);
					}
				}
			}
		}
	} catch (std::exception& e) {
		std::clog << "-- Infeasible due to shift that can not follow another shift. " << e.what()
		          << "\n";
		throw;
	}

	// Min/max consecutive days working.
	bool reported_min_consecutive_failure = false;
	bool reported_max_consecutive_failure = false;
	for (int p : staff_indices) {
		if (!problem.worker(p).has_consecutive_shifts_limit()) {
			continue;
		}

		int min_consecutive = problem.worker(p).consecutive_shifts_limit().min();
		int max_consecutive = problem.worker(p).consecutive_shifts_limit().max();

		vector<Sum> working_on_days;
		for (int d : day_indices) {
			working_on_days.emplace_back(working_on_day(p, d));
		}
		try {
			ip.add_min_consecutive_constraints(min_consecutive, working_on_days, true);
		} catch (std::exception& e) {
			if (!reported_min_consecutive_failure) {
				std::clog << "-- Infeasible due to min consecutive days working. " << e.what()
				          << "\n";
				reported_min_consecutive_failure = true;
			}
		}
		try {
			ip.add_max_consecutive_constraints(max_consecutive, working_on_days);
		} catch (std::exception& e) {
			if (!reported_max_consecutive_failure) {
				std::clog << "-- Infeasible due to max consecutive days working. " << e.what()
				          << "\n";
				reported_max_consecutive_failure = true;
			}
		}
	}
	if (reported_min_consecutive_failure || reported_max_consecutive_failure) {
		check(false, "Problem infeasible.");
	}

	// Min/max consecutive days working a specific shift.
	for (int s = 0; s < problem.shift_size(); ++s) {
		if (!problem.shift(s).has_consecutive_limit()) {
			continue;
		}
		for (int p : staff_indices) {
			vector<Sum> shifts;
			for (int d : range(problem.num_days())) {
				shifts.emplace_back(x[p][d][s]);
			}
			try {
				ip.add_min_consecutive_constraints(
				    problem.shift(s).consecutive_limit().min(), shifts, true);
				ip.add_max_consecutive_constraints(problem.shift(s).consecutive_limit().max(),
				                                   shifts);
			} catch (std::exception& e) {
				std::clog << "-- Infeasible due to consecuitive individual shift. " << e.what()
				          << "\n";
				throw;
			}
		}
	}

	// Min consecutive days off.
	try {
		for (int p : staff_indices) {
			if (!problem.worker(p).has_consecutive_days_off_limit()) {
				continue;
			}

			int max_consequtuve = problem.worker(p).consecutive_days_off_limit().max();
			check(max_consequtuve >= problem.num_days(),
			      "Max consecutive days off constraint is not implemented.");

			int min_consecutive = problem.worker(p).consecutive_days_off_limit().min();
			if (min_consecutive > 1) {
				vector<Sum> day_off;
				for (int d : day_indices) {
					day_off.emplace_back(1 - working_on_day(p, d));
				}
				ip.add_min_consecutive_constraints(min_consecutive, day_off, true);
			}
		}
	} catch (std::exception& e) {
		std::clog << "-- Infeasible due to min consequitive days off. " << e.what() << "\n";
		throw;
	}

	// Max working weekends.
	for (int p : staff_indices) {
		if (!problem.worker(p).has_working_weekends_limit()) {
			continue;
		}

		Sum working_weekends = 0;
		for (int d = 5; d < problem.num_days(); d += 7) {
			auto saturday = working_on_day(p, d);
			minimum_core_assert(d + 1 < problem.num_days());
			auto sunday = working_on_day(p, d + 1);

			auto weekend = ip.add_boolean();
			ip.add_constraint(saturday <= weekend);
			ip.add_constraint(sunday <= weekend);
			working_weekends += weekend;
		}
		check(problem.worker(p).working_weekends_limit().min() <= 0,
		      "Min working weekends constraint is not implemented.");
		ip.add_constraint(working_weekends <= problem.worker(p).working_weekends_limit().max());
	}

	// Days off.
	try {
		for (int p : staff_indices) {
			for (int day_off : problem.worker(p).required_day_off()) {
				ip.add_constraint(working_on_day(p, day_off) == 0);
			}
		}
	} catch (std::exception& e) {
		std::clog << "-- Infeasible due to days off. " << e.what() << "\n";
		throw;
	}

	Sum objective = problem.objective_constant();
	for (int p : staff_indices) {
		for (const auto& pref : problem.worker(p).shift_preference()) {
			objective += pref.cost() * x[p][pref.day()][pref.shift()];
		}

		for (const auto& pref : problem.worker(p).day_off_preference()) {
			objective += pref.cost() * (1 - working_on_day(p, pref.day()));
		}
	}

	// Cover
	for (auto& requirement : problem.requirement()) {
		Sum num_working = 0;
		for (int p : staff_indices) {
			num_working += x[p][requirement.day()][requirement.shift()];
		}
		auto under_slack = ip.add_variable(IP::Real);
		objective += requirement.under_cost() * under_slack;
		auto over_slack = ip.add_variable(IP::Real);
		objective += requirement.over_cost() * over_slack;
		ip.add_constraint(over_slack >= 0);
		ip.add_constraint(under_slack >= 0);
		ip.add_constraint(num_working + under_slack - over_slack == requirement.wanted());
	}
	ip.add_objective(objective);

	return objective;
}

Sum create_ip(IP& ip,
              const minimum::linear::proto::SchedulingProblem& problem,
              const vector<vector<vector<Sum>>>& x) {
	std::vector<int> staff_indices, day_indices;
	for (int p = 0; p < problem.worker_size(); ++p) {
		staff_indices.push_back(p);
	}
	for (int d = 0; d < problem.num_days(); ++d) {
		day_indices.push_back(d);
	}
	return create_ip(ip, problem, x, staff_indices, day_indices);
}

Sum create_ip(IP& ip,
              const minimum::linear::proto::SchedulingProblem& problem,
              const vector<vector<vector<BooleanVariable>>>& x) {
	vector<vector<vector<Sum>>> sums;
	for (auto& layer : x) {
		sums.emplace_back();
		for (auto& row : layer) {
			sums.back().emplace_back();
			for (auto& var : row) {
				sums.back().back().emplace_back(var);
			}
		}
	}
	return create_ip(ip, problem, sums);
}

int objective_value(const minimum::linear::proto::SchedulingProblem& problem,
                    const vector<vector<vector<int>>>& current_solution) {
	vector<vector<vector<Sum>>> x;
	for (int p = 0; p < problem.worker_size(); ++p) {
		x.emplace_back();
		for (int d = 0; d < problem.num_days(); ++d) {
			x.back().emplace_back();
			for (int s = 0; s < problem.shift_size(); ++s) {
				x.back().back().emplace_back(current_solution[p][d][s]);
			}
		}
	}
	IP ip;
	auto objective = create_ip(ip, problem, x);  // Will throw if infeasible.
	GlpkSolver solver;
	solver.set_silent(true);
	check(solver.solutions(&ip)->get(), "Problem is infeasible.");

	return objective.value();
}

std::string quick_solution_improvement(const proto::SchedulingProblem& problem,
                                       vector<vector<vector<int>>>& current_solution) {
	vector<vector<int>> current_assignment;
	vector<int> zero_shifts(problem.shift_size(), 0);
	current_assignment.resize(problem.num_days(), zero_shifts);
	auto wanted_assignment = current_assignment;
	auto over_cost = current_assignment;
	auto under_cost = current_assignment;

	for (auto& requirement : problem.requirement()) {
		wanted_assignment[requirement.day()][requirement.shift()] = requirement.wanted();
		over_cost[requirement.day()][requirement.shift()] = requirement.over_cost();
		under_cost[requirement.day()][requirement.shift()] = requirement.under_cost();
	}

	vector<vector<int>> working_on_day;
	vector<int> zero_days(problem.num_days(), 0);
	working_on_day.resize(problem.worker_size(), zero_days);
	vector<int> working_weekends(problem.worker_size(), 0);
	vector<int> working_minutes(problem.worker_size(), 0);
	vector<vector<int>> working_shift(problem.worker_size(), zero_shifts);

	for (int p : range(problem.worker_size())) {
		for (int d : range(problem.num_days())) {
			for (int s : range(problem.shift_size())) {
				if (current_solution[p][d][s] > 0) {
					working_on_day[p][d] = 1;
					working_minutes[p] += problem.shift(s).duration();
					working_shift[p][s] += 1;
					current_assignment[d][s] += 1;
				}
			}
		}

		for (int d = 5; d < problem.num_days(); d += 7) {
			if (working_on_day[p][d] || working_on_day[p][d + 1]) {
				working_weekends[p] += 1;
			}
		}
	}

	for (const int d : range(problem.num_days())) {
		for (const int s : range(problem.shift_size())) {
			if (wanted_assignment[d][s] - current_assignment[d][s] > 0) {
				// We would like this shift assigned. Can we find someone for it?

				for (const int p : range(problem.worker_size())) {
					// Would it be feasible and good to assign this shift to this worker?
					// If so, do that and return.

					int prev_shift = -1;
					if (working_on_day[p][d]) {
						for (int s2 : range(problem.shift_size())) {
							if (current_solution[p][d][s2] > 0) {
								prev_shift = s2;
								break;
							}
						}
					}
					if (prev_shift == s) {
						continue;
					}

					if (prev_shift < 0 && problem.worker(p).has_consecutive_shifts_limit()) {
						int consecutive = 1;
						for (int d2 = d - 1; d2 >= 0 && working_on_day[p][d2]; --d2, consecutive++)
							;
						for (int d2 = d + 1; d2 < problem.num_days() && working_on_day[p][d2];
						     ++d2, consecutive++)
							;
						if (consecutive > problem.worker(p).consecutive_shifts_limit().max()) {
							continue;
						}
						if (consecutive < problem.worker(p).consecutive_shifts_limit().min()) {
							continue;
						}
					}

					if (prev_shift < 0 && problem.worker(p).has_consecutive_days_off_limit()) {
						int min_consecutive = problem.worker(p).consecutive_days_off_limit().min();
						int consecutive = 0;
						for (int d2 = d - 1; d2 >= 0 && !working_on_day[p][d2]; --d2, consecutive++)
							;
						if (consecutive > 0 && consecutive < min_consecutive) {
							continue;
						}
						consecutive = 0;
						for (int d2 = d + 1; d2 < problem.num_days() && !working_on_day[p][d2];
						     ++d2, consecutive++)
							;
						if (consecutive > 0 && consecutive < min_consecutive) {
							continue;
						}
					}

					if (prev_shift >= 0
					    && working_shift[p][prev_shift] - 1
					           < problem.worker(p).shift_limit(prev_shift).min()) {
						continue;
					}
					if (working_shift[p][s] + 1 > problem.worker(p).shift_limit(s).max()) {
						continue;
					}

					if (problem.shift(s).has_consecutive_limit()) {
						minimum_core_assert(false, "Not implemented.");
					}
					if (prev_shift >= 0 && problem.shift(prev_shift).has_consecutive_limit()) {
						minimum_core_assert(false, "Not implemented.");
					}

					if (prev_shift < 0 && problem.worker(p).has_working_weekends_limit()) {
						bool new_weekend = false;
						if (d % 7 == 5) {
							if (!working_on_day[p][d] && !working_on_day[p][d + 1]) {
								new_weekend = true;
							}
						}
						if (d % 7 == 6) {
							if (!working_on_day[p][d - 1] && !working_on_day[p][d]) {
								new_weekend = true;
							}
						}
						if (new_weekend
						    && working_weekends[p]
						           >= problem.worker(p).working_weekends_limit().max()) {
							continue;
						}
					}

					if (problem.worker(p).has_time_limit()) {
						int new_duration = working_minutes[p] + problem.shift(s).duration();
						if (prev_shift >= 0) {
							new_duration -= problem.shift(prev_shift).duration();
						}
						if (new_duration > problem.worker(p).time_limit().max()) {
							continue;
						} else if (new_duration < problem.worker(p).time_limit().min()) {
							continue;
						}
					}

					if (prev_shift < 0) {
						bool feasible_day_off = true;
						for (int day_off : problem.worker(p).required_day_off()) {
							if (day_off == d) {
								feasible_day_off = false;
								break;
							}
						}
						if (!feasible_day_off) {
							continue;
						}
					}

					if (d > 0 && working_on_day[p][d - 1]) {
						int s1 = 0;
						for (; s1 < problem.shift_size(); ++s1) {
							if (current_solution[p][d - 1][s1] > 0) {
								break;
							}
						}

						bool feasible_shift_combination = true;
						for (int s2 : problem.shift(s1).forbidden_after()) {
							if (s2 == s) {
								feasible_shift_combination = false;
								break;
							}
						}
						if (!feasible_shift_combination) {
							continue;
						}
					}

					if (d < problem.num_days() - 1 && working_on_day[p][d + 1]) {
						int s1 = 0;
						for (; s1 < problem.shift_size(); ++s1) {
							if (current_solution[p][d + 1][s1] > 0) {
								break;
							}
						}

						bool feasible_shift_combination = true;
						for (int s2 : problem.shift(s).forbidden_after()) {
							if (s2 == s1) {
								feasible_shift_combination = false;
								break;
							}
						}
						if (!feasible_shift_combination) {
							continue;
						}
					}

					double delta_objective = -under_cost[d][s];
					for (const auto& pref : problem.worker(p).shift_preference()) {
						if (pref.day() == d && pref.shift() == s) {
							delta_objective += pref.cost();
						}
					}

					std::string operation;
					if (prev_shift >= 0) {
						int diff =
						    current_assignment[d][prev_shift] - wanted_assignment[d][prev_shift];
						if (diff > 0) {
							delta_objective -= over_cost[d][prev_shift];
						} else {
							delta_objective += under_cost[d][prev_shift];
						}

						for (const auto& pref : problem.worker(p).shift_preference()) {
							if (pref.day() == d && pref.shift() == prev_shift) {
								delta_objective -= pref.cost();
							}
						}

						operation = to_string("Changing shift ",
						                      problem.shift(prev_shift).id(),
						                      " to ",
						                      problem.shift(s).id(),
						                      " for ",
						                      problem.worker(p).id(),
						                      " on day ",
						                      d,
						                      " (",
						                      delta_objective,
						                      ").");
					} else {
						for (const auto& pref : problem.worker(p).day_off_preference()) {
							if (pref.day() == d) {
								delta_objective += pref.cost();
							}
						}

						operation = to_string("Assigning shift ",
						                      problem.shift(s).id(),
						                      " for ",
						                      problem.worker(p).id(),
						                      " on day ",
						                      d,
						                      " (",
						                      delta_objective,
						                      ").");
					}

					if (delta_objective < 0) {
						// This change is feasible and will improve the objective.
						// Do it and return.
						current_solution[p][d][s] = 1;
						if (prev_shift >= 0) {
							current_solution[p][d][prev_shift] = 0;
						}
						return operation;
					}
				}
			}
		}
	}

	for (const int p : range(problem.worker_size())) {
		for (const auto& pref : problem.worker(p).day_off_preference()) {
			int d = pref.day();
			int s = 0;
			for (; current_solution[p][d][s] == 0; ++s)
				;
			minimum_core_assert(s < problem.shift_size());
			if (working_on_day[p][d] && wanted_assignment[d][s] - current_assignment[d][s] >= 0) {
				// Maybe we can deassign this?

				minimum_core_assert(false, "Constraints for deassigning are not implemented.");

				current_solution[p][d][s] = 0;
				return to_string("Deassigning shift ",
				                 problem.shift(s).id(),
				                 " for ",
				                 problem.worker(p).id(),
				                 " on day ",
				                 d,
				                 ".");
			}
		}
	}
	return "";
}

}  // namespace linear
}  // namespace minimum