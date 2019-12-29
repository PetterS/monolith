#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <string_view>

#include <gflags/gflags.h>

#include <minimum/core/check.h>
#include <minimum/core/main.h>
#include <minimum/core/string.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/scheduling_util.h>
#include <minimum/linear/scs.h>
#include <minimum/linear/solver.h>

using namespace minimum::core;

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
	int num_periods = -1;

	struct Staff {
		std::string id;
		int min_minutes = -1;
		int max_minutes = -1;
	};

	std::map<std::string, Staff> staff;

	struct Period {
		std::string human_readable_time;
		int day = -1;
		std::vector<int> min_cover;
		std::vector<int> max_cover;
	};

	std::vector<Period> periods;

	void expect_line(std::istream& file, std::string_view expected) {
		std::string tag;
		if (!expected.empty()) {
			for (int i = 0; i < 1000 && tag.empty(); ++i) {
				tag = remove_space(read_line(file));
			}
		}
		check(tag == expected, "Expected ", expected, " but got \"", tag, "\".");
	}

	RetailProblem(std::istream& file) {
		expect_line(file, "SECTION_HORIZON");
		num_days = from_string<int>(read_line(file));

		expect_line(file, "SECTION_TASKS");
		num_tasks = from_string<int>(read_line(file));

		expect_line(file, "SECTION_STAFF");
		while (true) {
			auto line = read_line(file);
			if (line.empty()) {
				break;
			}
			auto parts = split(line, ',');
			check(parts.size() == 3, "Invalid staff line: ", line);
			auto id = parts[0];
			auto& person = staff[id];
			person.id = parts[0];
			person.min_minutes = from_string<int>(parts[1]);
			person.max_minutes = from_string<int>(parts[2]);
		}

		num_periods = num_days * 24 * 4;
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
			auto parts = split(line, ',');
			check(parts.size() == 5, "Invalid cover line: ", line);
			int day = from_string<int>(parts[0]) - 1;
			auto time = parts[1];
			int task = from_string<int>(parts[2]) - 1;
			minimum_core_assert(0 <= task && task < num_tasks);

			int p = time_to_period(day, time);
			auto& period = periods.at(p);
			minimum_core_assert(period.human_readable_time == ""
			                    || period.human_readable_time == time);
			minimum_core_assert(period.day < 0 || period.day == day);
			period.human_readable_time = time;
			period.day = day;
			period.min_cover[task] = from_string<int>(parts[3]);
			period.max_cover[task] = from_string<int>(parts[4]);
		}

		for (int i = 0; i < 10; ++i) {
			expect_line(file, "");
		}
	}

	int time_to_period(int day, const std::string& s) {
		auto parts = split(s, '-');
		check(parts.size() == 2, "Invalid period: ", s);
		auto hm_parts = split(parts[0], ':');
		check(parts.size() == 2, "Invalid period: ", s);
		int hour = from_string<int>(hm_parts[0]);
		int minute = from_string<int>(hm_parts[1]);
		check(minute == 0 || minute == 15 || minute == 30 || minute == 45, "Invalid period: ", s);

		int hour4 = 4 * hour + minute / 15;

		// First day starts at 06:00.
		int period = (hour4 - 6 * 4) + (day * 24 * 4);
		minimum_core_assert(0 <= period && period < num_periods);
		return period;
	}

	void print_info() const {
		std::cerr << "Number of days: " << num_days << '\n';
		std::cerr << "Number of periods: " << num_periods << '\n';
		std::cerr << "Number of tasks: " << num_tasks << '\n';
		std::cerr << "Number of staff: " << staff.size() << '\n';
	}
};

int main_program(int num_args, char* args[]) {
	using namespace std;
	if (num_args <= 1) {
		cerr << "Need a file name.\n";
		return 1;
	}
	ifstream input(args[1]);
	RetailProblem problem(input);
	problem.print_info();

	return 0;
}

int main(int num_args, char* args[]) { return main_runner(main_program, num_args, args); }
