#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>

#include <absl/strings/str_format.h>

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/core/string.h>

namespace minimum {
namespace linear {

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

class RetailProblem {
   public:
	int num_days = 0;
	int num_tasks = 0;
	const int min_shift_length = 6 * 4;
	const int max_shift_length = 10 * 4;

	struct Staff {
		std::string id;
		int min_minutes = -1;
		int max_minutes = -1;
	};

	std::vector<Staff> staff;

	struct Period {
		std::string human_readable_time;
		int day = -1;
		int hour4_since_start = -1;
		std::vector<int> min_cover;
		std::vector<int> max_cover;
		bool can_start = false;
	};

	std::vector<Period> periods;

	int num_cover_constraints() const { return num_tasks * periods.size(); }

	int cover_constraint_index(int period, int task) const {
		minimum_core_assert(0 <= period && period < periods.size());
		minimum_core_assert(0 <= task && task < num_tasks);
		return period * num_tasks + task;
	}
	int cover_constraint_to_period(int c) const { return c / num_tasks; }
	int cover_constraint_to_task(int c) const { return c % num_tasks; }

	void expect_line(std::istream& file, std::string_view expected) {
		std::string tag;
		if (!expected.empty()) {
			for (int i = 0; i < 1000 && tag.empty(); ++i) {
				tag = remove_space(read_line(file));
			}
		}
		minimum::core::check(tag == expected, "Expected ", expected, " but got \"", tag, "\".");
	}

	RetailProblem(std::istream& file) {
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
			minimum_core_assert(period.human_readable_time == ""
			                    || period.human_readable_time == time);
			minimum_core_assert(period.day < 0 || period.day == day);
			period.human_readable_time = time;
			period.hour4_since_start = 24 * 4 * day + hour4;
			period.day = day;
			period.min_cover[task] = minimum::core::from_string<int>(parts[3]);
			period.max_cover[task] = minimum::core::from_string<int>(parts[4]);
			period.can_start = day < num_days
			                   && ((6 * 4 <= hour4 && hour4 <= 10 * 4)
			                       || (14 * 4 <= hour4 && hour4 <= 18 * 4) || (20 * 4 <= hour4));
		}
		// Sanity check that all periods have been present in the problem formulation.
		for (auto& period : periods) {
			minimum_core_assert(period.day >= 0);
		}

		for (int i = 0; i < 10; ++i) {
			expect_line(file, "");
		}
	}

	int time_to_hour4(const std::string& s) const {
		auto parts = minimum::core::split(s, '-');
		minimum::core::check(parts.size() == 2, "Invalid period: ", s);
		auto hm_parts = minimum::core::split(parts[0], ':');
		minimum::core::check(parts.size() == 2, "Invalid period: ", s);
		int hour = minimum::core::from_string<int>(hm_parts[0]);
		int minute = minimum::core::from_string<int>(hm_parts[1]);
		minimum::core::check(
		    minute == 0 || minute == 15 || minute == 30 || minute == 45, "Invalid period: ", s);
		return 4 * hour + minute / 15;
	}

	std::string format_hour4(int hour4) const {
		int hour = hour4 / 4;
		int minutes = (hour4 % 4) * 15;
		return absl::StrFormat("%02d:%02d:00", hour, minutes);
	}

	int time_to_period(int day, int hour4) {
		// First day starts at 06:00.
		int period = (hour4 - 6 * 4) + (day * 24 * 4);
		minimum_core_assert(0 <= period && period < periods.size());
		return period;
	}

	// Returns the index of the first period of the day.
	int day_start(int day) const {
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

	void print_info() const {
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

	// Returns the solution value and throws if it is infeasible.
	int check_feasibility(const std::vector<std::vector<std::vector<int>>>& solution) const {
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

				// TODO: Check individual roster feasibility.
			}
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

	int save_solution(std::string problem_filename,
	                  std::string filename,
	                  const std::vector<std::vector<std::vector<int>>>& solution) const {
		using namespace minimum::core;

		auto objective_value = check_feasibility(solution);

		std::ofstream file(filename);
		file << "<Solution xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
		        "xsi:noNamespaceSchemaLocation=\"Solution.xsd\">\n";
		file << "<ProblemFile>" << problem_filename << "</ProblemFile>\n";
		file << absl::StrFormat("<Penalty>%d</Penalty>\n", objective_value);
		file << "<TimeStamp>2019-10-30T12:43:35</TimeStamp>\n";
		file << "<Author>Petter Strandmark</Author>\n";
		file << "<Reference>\n";
		file << "</Reference>\n";
		file << "<Machine></Machine>\n";
		file << "<SolveTime>0.00:01:00</SolveTime>\n";
		for (int staff_index = 0; staff_index < staff.size(); ++staff_index) {
			auto& current = solution[staff_index];
			file << "<Employee ID = \"" << staff[staff_index].id << "\">\n";

			int current_task = -1;
			int current_period = -1;
			for (int p = 0; p < periods.size(); ++p) {
				auto& period = periods.at(p);
				int hour4 = period.hour4_since_start;
				int start_on_day_hour4 = time_to_hour4(period.human_readable_time);

				int task = -1;
				for (int t : range(num_tasks)) {
					if (current[p][t] != 0) {
						check(task == -1, "Multiple tasks in the same period.");
						task = t;
					}
				}

				// if (staff_index == 0) {
				//	if (p == 0) {
				//		for (int i = 0; i < 6 * 4; ++i) {
				//			std::cerr << " ";
				//		}
				//	}
				//	if (p == day_start(period.day)) {
				//		std::cerr << period.human_readable_time << " ";
				//	}
				//	std::cerr << (task == -1 ? "." : "1");
				//	if (period.day < num_days && p == day_start(period.day + 1) - 1) {
				//		std::cerr << " " << period.human_readable_time << "\n";
				//	}
				//}

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
					if (last_period) {
						final_p = p;
					}
					int length = (periods[final_p].hour4_since_start
					              - periods[current_period].hour4_since_start + 1)
					             * 15;

					file << absl::StrFormat(
					    "        <Task "
					    "WorkLength=\"%d\"><ID>%d</ID><Length>%d</Length></Task>\n",
					    length,
					    current_task + 1,
					    length);
					current_task = task;
					current_period = hour4;
					if (task == -1 || last_period) {
						file << "    </Shift>\n";
					}
				}
			}
			file << "</Employee>\n";
		}
		file << "</Solution>\n";
		return objective_value;
	}
};

}  // namespace linear
}  // namespace minimum
