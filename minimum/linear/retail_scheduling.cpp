#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>

#include <absl/strings/str_format.h>
#include <rapidxml.hpp>

#include <minimum/core/check.h>
#include <minimum/core/grid.h>
#include <minimum/core/range.h>
#include <minimum/core/string.h>
#include <minimum/linear/retail_scheduling.h>

using namespace std;
using namespace minimum::core;

namespace minimum {
namespace linear {

namespace {
std::string read_line(std::istream& file) {
	std::string line;
	std::getline(file, line);
	if (!file) {
		return "";
	}

	// Check for comment.
	for (auto ch : line) {
		if (ch == '#') {
			return read_line(file);
		} else if (::isspace(ch)) {
			break;
		}
	}

	if (!line.empty() && line[line.size() - 1] == '\r') {
		line = line.substr(0, line.size() - 1);
	}

	return line;
}

std::string remove_space(std::string in) {
	in.erase(std::remove_if(in.begin(), in.end(), ::isspace), in.end());
	return in;
}
}  // namespace

int RetailProblem::num_cover_constraints() const { return num_tasks * periods.size(); }

int RetailProblem::cover_constraint_index(int period, int task) const {
	minimum_core_assert(0 <= period && period < periods.size());
	minimum_core_assert(0 <= task && task < num_tasks);
	return period * num_tasks + task;
}
int RetailProblem::cover_constraint_to_period(int c) const { return c / num_tasks; }
int RetailProblem::cover_constraint_to_task(int c) const { return c % num_tasks; }

void RetailProblem::expect_line(std::istream& file, std::string_view expected) {
	std::string tag;
	if (!expected.empty()) {
		for (int i = 0; i < 1000 && tag.empty(); ++i) {
			tag = remove_space(read_line(file));
		}
	}
	minimum::core::check(tag == expected, "Expected ", expected, " but got \"", tag, "\".");
}

RetailProblem::RetailProblem(std::istream& file) {
	minimum::core::check(!!file, "Can not read from stream.");
	expect_line(file, "SECTION_HORIZON");
	num_days = minimum::core::from_string<int>(read_line(file));

	expect_line(file, "SECTION_TASKS");
	num_tasks = minimum::core::from_string<int>(read_line(file));

	expect_line(file, "SECTION_STAFF");
	while (true) {
		auto line = read_line(file);
		if (line.empty()) {
			break;
		}
		auto parts = minimum::core::split(line, ',');
		minimum::core::check(parts.size() == 3, "Invalid staff line: ", line);
		auto id = parts[0];
		Staff person;
		person.id = parts[0];
		person.min_minutes = minimum::core::from_string<int>(parts[1]);
		person.max_minutes = minimum::core::from_string<int>(parts[2]);

		if (person.min_minutes % 15 != 0) {
			// A minimum number of minutes of 29 is equivalent to 30 because
			// all shifts are divisible by 15.
			person.min_minutes += (15 - (person.min_minutes % 15));
		}
		if (person.max_minutes % 15 != 0) {
			// A maximum number of minutes of 32 is equivalent to 30 because
			// all shifts are divisible by 15.
			person.max_minutes -= person.max_minutes % 15;
		}

		staff.emplace_back(person);
	}

	int num_periods = num_days * 24 * 4;
	Period period;
	period.min_cover.resize(num_tasks, 0);
	period.max_cover.resize(num_tasks, 0);
	periods.resize(num_periods, period);

	expect_line(file, "SECTION_COVER");
	while (true) {
		auto line = read_line(file);
		if (line.empty()) {
			break;
		}
		auto parts = minimum::core::split(line, ',');
		minimum::core::check(parts.size() == 5, "Invalid cover line: ", line);
		int day = minimum::core::from_string<int>(parts[0]) - 1;
		auto time = parts[1];
		int task = minimum::core::from_string<int>(parts[2]) - 1;
		minimum_core_assert(0 <= task && task < num_tasks);

		int hour4 = time_to_hour4(time);
		int p = time_to_period(day, hour4);
		auto& period = periods.at(p);
		minimum_core_assert(period.human_readable_time == "" || period.human_readable_time == time);
		minimum_core_assert(period.day < 0 || period.day == day);
		period.human_readable_time = time;
		period.hour4_since_start = 24 * 4 * day + hour4;
		period.day = day;
		period.min_cover[task] = minimum::core::from_string<int>(parts[3]);
		period.max_cover[task] = minimum::core::from_string<int>(parts[4]);
		period.can_start = hour4 == 0 || (6 * 4 <= hour4 && hour4 <= 10 * 4)
		                   || (14 * 4 <= hour4 && hour4 <= 18 * 4) || (20 * 4 <= hour4);
	}
	// Sanity check that all periods have been present in the problem formulation.
	for (auto& period : periods) {
		minimum_core_assert(period.day >= 0);
	}

	for (int i = 0; i < 10; ++i) {
		expect_line(file, "");
	}
}

int RetailProblem::time_to_hour4(const std::string& s) const {
	auto parts = minimum::core::split(s, '-');
	minimum::core::check(parts.size() == 2, "Invalid period: ", s);
	auto hm_parts = minimum::core::split(parts[0], ':');
	minimum::core::check(hm_parts.size() == 2, "Invalid period: ", s);
	int hour = minimum::core::from_string<int>(hm_parts[0]);
	int minute = minimum::core::from_string<int>(hm_parts[1]);
	minimum::core::check(
	    minute == 0 || minute == 15 || minute == 30 || minute == 45, "Invalid period: ", s);
	return 4 * hour + minute / 15;
}

std::string RetailProblem::format_hour4(int hour4) const {
	int hour = hour4 / 4;
	int minutes = (hour4 % 4) * 15;
	return absl::StrFormat("%02d:%02d:00", hour, minutes);
}

int RetailProblem::time_to_period(int day, int hour4) const {
	// First day starts at 06:00.
	int period = (hour4 - 6 * 4) + (day * 24 * 4);
	minimum_core_assert(0 <= period && period < periods.size());
	return period;
}

// Returns the index of the first period of the day.
int RetailProblem::day_start(int day) const {
	if (day == 0) {
		return 0;
	}
	if (day == num_days) {
		return periods.size();
	}
	// First day starts at 06:00.
	int period = 18 * 4 + ((day - 1) * 24 * 4);
	minimum_core_assert(time_to_hour4(periods.at(period).human_readable_time) == 0);
	return period;
}

void RetailProblem::print_info() const {
	std::cerr << "Number of days: " << num_days << '\n';
	std::cerr << "Number of periods: " << periods.size() << '\n';
	std::cerr << "Number of tasks: " << num_tasks << '\n';
	std::cerr << "Number of staff: " << staff.size() << '\n';
	std::cerr << "Number of cover constaints: " << num_tasks << " * " << periods.size() << " = "
	          << num_cover_constraints() << '\n';
	for (auto& member : staff) {
		std::cerr << "    Staff " << member.id << ": " << member.min_minutes << ", "
		          << member.max_minutes << " minutes of work.\n";
	}
}

void RetailProblem::check_feasibility_for_staff(
    const int staff_index, const std::vector<std::vector<int>>& solution) const {
	auto& member = staff.at(staff_index);
	int minutes_worked = 0;
	// Shift length.
	int current_task = -1;
	int current_period = -1;
	for (int p = 0; p < periods.size(); ++p) {
		auto& period = periods.at(p);

		int task = -1;
		for (int t : range(num_tasks)) {
			if (solution[p][t] != 0) {
				check(task == -1, "Multiple tasks in the same period.");
				task = t;
				minutes_worked += 15;
			}
		}

		bool last_period = p == periods.size() - 1;
		if (current_task == -1 && task >= 0) {
			// Start of the first task of the shift.
			current_task = task;
			current_period = p;
		} else if (current_task != -1 && (current_task != task || last_period)) {
			int final_p = p - 1;
			if (last_period) {
				final_p = p;
			}
			int task_length = (final_p - current_period + 1) * 15;
			check(task_length >= 60, "Task length of ", task_length, " is too short.");
		}
	}
	check(member.min_minutes <= minutes_worked && minutes_worked <= member.max_minutes,
	      member.min_minutes,
	      " <= ",
	      minutes_worked,
	      " <= ",
	      member.max_minutes,
	      " is not fulfilled for staff ",
	      member.id);
}

// Returns the solution value and throws if it is infeasible.
int RetailProblem::check_feasibility(
    const std::vector<std::vector<std::vector<int>>>& solution) const {
	using minimum::core::make_grid;
	using minimum::core::range;
	std::string errors;
	auto cover = make_grid<int>(periods.size(), num_tasks);

	for (int staff_index : range(staff.size())) {
		for (int period_index : range(periods.size())) {
			for (int t : range(num_tasks)) {
				if (solution.at(staff_index).at(period_index).at(t) == 1) {
					cover.at(period_index).at(t) += 1;
				}
			}
		}
		check_feasibility_for_staff(staff_index, solution.at(staff_index));
	}

	int objective = 0;

	for (int period_index : range(periods.size())) {
		auto& period = periods.at(period_index);
		for (int t : range(num_tasks)) {
			auto current = cover.at(period_index).at(t);
			if (current < period.min_cover.at(t)) {
				errors += minimum::core::to_string("\nCover violation for task ",
				                                   t,
				                                   " on day ",
				                                   period.day,
				                                   " at ",
				                                   period.human_readable_time,
				                                   ". Wanted ",
				                                   period.min_cover.at(t),
				                                   "; have ",
				                                   current,
				                                   ".");
			} else if (current > period.max_cover.at(t)) {
				auto delta = current - period.max_cover.at(t);
				objective += delta * delta;
			}
		}
	}

	if (errors.length() > 0) {
		throw std::runtime_error("Feasibility errors:" + errors);
	}

	return objective;
}

std::string RetailProblem::solution_to_string(
    const std::vector<std::vector<std::vector<int>>>& solution) const {
	check_feasibility(solution);
	ostringstream out;
	for (int staff_index = 0; staff_index < staff.size(); ++staff_index) {
		for (int p = 0; p < periods.size(); ++p) {
			for (int t : range(num_tasks)) {
				out << solution.at(staff_index).at(p).at(t) << " ";
			}
		}
	}
	return out.str();
}

std::vector<std::vector<std::vector<int>>> RetailProblem::string_to_solution(
    const std::string input) const {
	auto solution = make_grid<int>(staff.size(), periods.size(), num_tasks);
	istringstream in(input);
	for (int staff_index = 0; staff_index < staff.size(); ++staff_index) {
		for (int p = 0; p < periods.size(); ++p) {
			for (int t : range(num_tasks)) {
				in >> solution.at(staff_index).at(p).at(t);
			}
		}
	}
	check_feasibility(solution);
	return solution;
}

std::string RetailProblem::save_solution(std::string problem_filename,
                                 const std::vector<std::vector<std::vector<int>>>& solution,
                                 double solution_time,
                                 std::string timestamp,
                                 std::string solution_method) const {
	auto objective_value = check_feasibility(solution);

	std::ostringstream file;
	file << "<Solution xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
	        "xsi:noNamespaceSchemaLocation=\"Solution.xsd\">\n";
	file << "<ProblemFile>" << problem_filename << "</ProblemFile>\n";
	file << absl::StrFormat("<Penalty>%d</Penalty>\n", objective_value);
	file << "<TimeStamp>" << timestamp << "</TimeStamp>\n";
	file << "<Author>Petter Strandmark</Author>\n";
	file << "<Reference>\n";
	file << "  " << solution_method << "\n";
	file << "</Reference>\n";
	file << "<Machine></Machine>\n";
	int days = solution_time / 60 / 60 / 24;
	int hours = solution_time / 60 / 60 - days * 24;
	int minutes = solution_time / 60 - hours * 60 - days * 24 * 60;
	int seconds = solution_time - minutes * 60 - hours * 60 * 60 - days * 24 * 60 * 60;
	file << absl::StrFormat(
	    "<SolveTime>%d.%02d:%02d:%02d</SolveTime>\n", days, hours, minutes, seconds);
	for (int staff_index = 0; staff_index < staff.size(); ++staff_index) {
		auto& current = solution[staff_index];
		file << "<Employee ID=\"" << staff[staff_index].id << "\">\n";
		int working_minutes = 0;
		int working_minutes_written = 0;
		int current_task = -1;
		int current_period = -1;
		for (int p = 0; p < periods.size(); ++p) {
			auto& period = periods.at(p);
			int start_on_day_hour4 = time_to_hour4(period.human_readable_time);

			int task = -1;
			for (int t : range(num_tasks)) {
				if (current[p][t] != 0) {
					minimum_core_assert(task == -1);
					task = t;
					working_minutes += 15;
				}
			}

			bool last_period = p == periods.size() - 1;
			if (current_task == -1 && task >= 0) {
				// Start of the first task of the shift.
				file << "    <Shift>\n";
				file << absl::StrFormat("        <Day>%d</Day>\n", period.day);
				file << "        <Start>" << format_hour4(start_on_day_hour4) << "</Start>\n";
				current_task = task;
				current_period = p;
			} else if (current_task != -1 && (current_task != task || last_period)) {
				int final_p = p - 1;
				if (current_task == task && last_period) {
					final_p = p;
				}
				int length = (final_p - current_period + 1) * 15;

				file << absl::StrFormat(
				    "        <Task "
				    "WorkLength=\"%d\"><ID>%d</ID><Length>%d</Length></Task>\n",
				    length,
				    current_task + 1,
				    length);
				working_minutes_written += length;
				current_task = task;
				current_period = p;
				if (task == -1 || last_period) {
					file << "    </Shift>\n";
				}
			}
		}
		file << "</Employee>\n";
		minimum_core_assert(working_minutes_written == working_minutes,
		                    "Solution error for employee ",
		                    staff[staff_index].id);
	}
	file << "</Solution>\n";
	return file.str();
}

std::vector<std::vector<std::vector<int>>> RetailProblem::load_solution(std::istream& in) const {
	auto solution = make_grid<int>(staff.size(), periods.size(), num_tasks);
	check(!!in, "Invalid stream");

	string contents(istreambuf_iterator<char>(in), {});
	contents += '\0';
	rapidxml::xml_document<> doc;
	doc.parse<0>(contents.data());
	auto root = doc.first_node("Solution");
	check(root != nullptr, "Did not find <Solution>.");
	for (auto employee = root->first_node("Employee"); employee;
	     employee = employee->next_sibling("Employee")) {
		auto id_attr = employee->first_attribute("ID");
		check(id_attr != nullptr, "<Employee> without ID");
		auto id = id_attr->value();
		int staff_index = -1;
		for (int i : range(staff.size())) {
			if (staff[i].id == id) {
				staff_index = i;
			}
		}
		check(staff_index >= 0, "Invalid staff ID");

		for (auto shift = employee->first_node("Shift"); shift;
		     shift = shift->next_sibling("Shift")) {
			auto day_node = shift->first_node("Day");
			check(day_node != nullptr, "No <Day> node.");
			int day = from_string<int>(day_node->value());
			check(0 <= day && day <= num_days, "Invalid day");

			auto start = shift->first_node("Start");
			check(start != nullptr, "No <Start> node.");

			auto hms_parts = minimum::core::split(start->value(), ':');
			check(hms_parts.size() == 3, "Invalid time: ", start->value());
			int hour = from_string<int>(hms_parts[0]);
			int minute = from_string<int>(hms_parts[1]);
			check(minute == 0 || minute == 15 || minute == 30 || minute == 45,
			      "Invalid period: ",
			      start->value());
			int hour4 = 4 * hour + minute / 15;
			int period = time_to_period(day, hour4);
			for (auto task = shift->first_node("Task"); task; task = task->next_sibling("Task")) {
				auto id_node = task->first_node("ID");
				check(id_node != nullptr, "No <ID> node.");
				int task_id = from_string<int>(id_node->value()) - 1;
				check(0 <= task_id && task_id < num_tasks, "Invalid task ID");

				auto length_node = task->first_node("Length");
				check(length_node != nullptr, "No <Length> node.");
				int length = from_string<int>(length_node->value());

				check(length % 15 == 0, "Invalid task length");
				int finish = period + length / 15;
				for (; period < finish; ++period) {
					solution[staff_index][period][task_id] = 1;
				}
			}
		}
	}
	return solution;
}

}  // namespace linear
}  // namespace minimum
