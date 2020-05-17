#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>

#include <minimum/core/check.h>
#include <minimum/core/grid.h>
#include <minimum/core/range.h>
#include <minimum/core/string.h>
#include <minimum/linear/export.h>

namespace minimum {
namespace linear {

class MINIMUM_LINEAR_API RetailProblem {
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

	int num_cover_constraints() const;

	int cover_constraint_index(int period, int task) const;
	int cover_constraint_to_period(int c) const;
	int cover_constraint_to_task(int c) const;

	RetailProblem(std::istream& file);
	int time_to_hour4(const std::string& s) const;

	std::string format_hour4(int hour4) const;

	int time_to_period(int day, int hour4);

	// Returns the index of the first period of the day.
	int day_start(int day) const;

	void print_info() const;

	void check_feasibility_for_staff(const int staff_index,
	                                 const std::vector<std::vector<int>>& solution) const;

	// Returns the solution value and throws if it is infeasible.
	int check_feasibility(const std::vector<std::vector<std::vector<int>>>& solution) const;

	std::string solution_to_string(
	    const std::vector<std::vector<std::vector<int>>>& solution) const;
	std::vector<std::vector<std::vector<int>>> string_to_solution(const std::string input) const;

	int save_solution(std::string problem_filename,
	                  std::string filename,
	                  const std::vector<std::vector<std::vector<int>>>& solution) const;

   private:
	void expect_line(std::istream& file, std::string_view expected);
};

}  // namespace linear
}  // namespace minimum
