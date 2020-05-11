#pragma once
#include <iostream>
#include <string>
#include <string_view>

#include <minimum/core/check.h>
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
		bool can_start = false;
	};

	std::vector<Period> periods;

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
			auto& person = staff[id];
			person.id = parts[0];
			person.min_minutes = minimum::core::from_string<int>(parts[1]);
			person.max_minutes = minimum::core::from_string<int>(parts[2]);
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
			period.day = day;
			period.min_cover[task] = minimum::core::from_string<int>(parts[3]);
			period.max_cover[task] = minimum::core::from_string<int>(parts[4]);
			period.can_start = (6 * 4 <= hour4 && hour4 <= 10 * 4)
			                   || (14 * 4 <= hour4 && hour4 <= 18 * 4) || (20 * 4 <= hour4);
		}
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

	int time_to_period(int day, int hour4) {
		// First day starts at 06:00.
		int period = (hour4 - 6 * 4) + (day * 24 * 4);
		minimum_core_assert(0 <= period && period < periods.size());
		return period;
	}

	void print_info() const {
		std::cerr << "Number of days: " << num_days << '\n';
		std::cerr << "Number of periods: " << periods.size() << '\n';
		std::cerr << "Number of tasks: " << num_tasks << '\n';
		std::cerr << "Number of staff: " << staff.size() << '\n';
	}
};

}  // namespace linear
}  // namespace minimum
