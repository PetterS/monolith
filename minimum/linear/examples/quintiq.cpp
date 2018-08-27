// Solves the Quintiq Puzzle.
// http://careers.quintiq.com/puzzle.html
//
// Press print screen before running the program and it
// will parse that image.
//
// Need to run as administrator in order to use the mouse
// to fill in the solution afterwards.

#include <cmath>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

#ifdef WIN32
#include <windows.h>
#endif

#include <minimum/core/check.h>
#include <minimum/core/color.h>
#include <minimum/core/range.h>
#include <minimum/core/scope_guard.h>
#include <minimum/core/time.h>
#include <minimum/linear/glpk.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/solver.h>

using namespace std;
using namespace minimum::core;
using namespace minimum::linear;

constexpr int num_drivers = 11;
constexpr int num_days = 14;
constexpr int num_lines = 3;

constexpr double shift_preference_granted = -3;
constexpr double driver_want_off = -4;
constexpr double long_rest_bonus = -5;
constexpr double unassigned_shift_cost = 20;
constexpr double late_shift_plus_early = 30;

constexpr int max_late_shifts = 3;
constexpr double four_consecutive_late_shifts = 10;

constexpr int target_late_shifts = 4;
constexpr double shift_assigned_not_four = 8;

typedef std::vector<std::vector<std::vector<float>>> Image;

//
// Captures the image in the clipboard.
//
Image get_clipboard() {
#ifdef WIN32
	using namespace std;
	minimum_core_assert(OpenClipboard(nullptr) != 0);
	at_scope_exit(CloseClipboard());

	HBITMAP hbitmap = (HBITMAP)GetClipboardData(CF_BITMAP);
	if (hbitmap != NULL) {
		BITMAP bitmap;
		GetObject(hbitmap, sizeof(BITMAP), &bitmap);

		auto bytes_per_pixel = bitmap.bmBitsPixel / 8;
		minimum_core_assert(bytes_per_pixel >= 3);

		vector<unsigned char> data(bitmap.bmHeight * bitmap.bmWidthBytes);
		minimum_core_assert(
		    GetBitmapBits(hbitmap, bitmap.bmHeight * bitmap.bmWidthBytes, data.data()) != 0);
		Image image;
		image.reserve(bitmap.bmHeight);
		for (int r = 0; r < bitmap.bmHeight; ++r) {
			image.emplace_back();
			image.back().reserve(bitmap.bmWidth);
			for (int c = 0; c < bitmap.bmWidth; ++c) {
				image.back().emplace_back();
				image.back().back().reserve(bytes_per_pixel);
				for (int i = 0; i < bytes_per_pixel; ++i) {
					image.back().back().emplace_back(
					    data.at(r * bitmap.bmWidthBytes + c * bytes_per_pixel + i));
				}
			}
		}
		minimum_core_assert(image.size() == bitmap.bmHeight);
		minimum_core_assert(image.back().size() == bitmap.bmWidth);
		minimum_core_assert(image.back().back().size() == bytes_per_pixel);
		return image;
	} else {
		return {};
	}
#else
	return {};
#endif
}

//
// Finds the orange rows at the start of the puzzle. Returns
// a list of pairs (row, column) of the last pixel of each
// orange row.
//
vector<pair<int, int>> get_orange_rows(const Image& image) {
	if (image.empty()) {
		return {};
	}

	map<int, int> rows_map;
	for (auto r : range(image.size())) {
		for (auto c : range(image.back().size())) {
			float delta = 0;
			delta += abs(float(image[r][c][2]) - 255.0f);
			delta += abs(float(image[r][c][1]) - 128.0f);
			delta += abs(float(image[r][c][0]) - 0.0f);
			if (delta <= 3) {
				// cout << "(" << r << ", " << c << ")\n";
				rows_map[r] = max(rows_map[r], int(c));
			}
		}
	}
	vector<pair<int, int>> rows;
	for (auto r : rows_map) {
		if (r.first > 100) {
			rows.emplace_back(r.first, r.second);
		}
	}
	auto e = unique(
	    rows.begin(), rows.end(), [](auto a, auto b) { return abs(a.first - b.first) < 10; });
	rows.erase(e, rows.end());

	if (rows.size() != num_drivers) {
		return {};
	} else {
		return rows;
	}
}

Image compute_integral_image(Image integral) {
	if (integral.empty()) {
		return integral;
	}

	for (auto r : range(integral.size())) {
		for (auto c : range(size_t(1), integral.back().size())) {
			for (auto i : range(integral.back().back().size())) {
				integral[r][c][i] += integral[r][c - 1][i];
			}
		}
	}
	for (auto r : range(size_t(1), integral.size())) {
		for (auto c : range(integral.back().size())) {
			for (auto i : range(integral.back().back().size())) {
				integral[r][c][i] += integral[r - 1][c][i];
			}
		}
	}
	return integral;
}

// Computes a box integral using an integral image.
float box(const Image& integral, int r1, int c1, int r2, int c2, int i) {
	auto sum = float(integral.at(r1).at(c1).at(i)) + float(integral.at(r2).at(c2).at(i))
	           - float(integral.at(r1).at(c2).at(i)) - float(integral.at(r2).at(c1).at(i));
	return sum / ((r2 - r1) * (c2 - c1));
}

// Reads driver qualifications from the captured image.
void parse_driver_qualified(const Image& integral,
                            const vector<pair<int, int>>& rows,
                            int driver_qualified[num_drivers][num_lines],
                            int _,
                            int X) {
	constexpr int square_width = 66;
	const int square_height = rows[1].first - rows[0].first;

	for (auto driver : range(num_drivers)) {
		driver_qualified[driver][0] = X;
		driver_qualified[driver][1] = X;
		driver_qualified[driver][2] = X;

		int r = rows[driver].first - square_height;
		int c = rows[driver].second - num_days * square_width;
		auto red1 = box(integral, r + 13, c - 77, r + 27, c - 61, 2);
		auto red2 = box(integral, r + 13, c - 55, r + 27, c - 38, 2);
		float tol = 10;
		if (abs(red1 - 189) < tol) {
			driver_qualified[driver][0] = _;
		}
		if (abs(red1 - 121) < tol || abs(red2 - 121) < tol) {
			driver_qualified[driver][1] = _;
		}
		if (abs(red1 - 91) < tol || abs(red2 - 91) < tol) {
			driver_qualified[driver][2] = _;
		}
	}
}

// Reads the driver information matrix from the captured image.
void parse_driver_info(const Image& integral,
                       const vector<pair<int, int>>& rows,
                       int driver_info[num_drivers][num_days],
                       int _,
                       int X,
                       int O,
                       int E,
                       int L) {
	constexpr int square_width = 66;
	const int square_height = rows[1].first - rows[0].first;

	for (auto driver : range(num_drivers)) {
		int r = rows[driver].first - square_height;
		int c = rows[driver].second - num_days * square_width;
		for (auto d : range(num_days)) {
			auto red = box(integral,
			               r + 12,
			               c + d * square_width,
			               r + square_height,
			               c + (d + 1) * square_width,
			               2);
			auto blue = box(integral,
			                r + 12,
			                c + d * square_width,
			                r + square_height,
			                c + (d + 1) * square_width,
			                0);
			if (blue < 100) {
				driver_info[driver][d] = X;
			} else if (red < 180) {
				driver_info[driver][d] = O;
			} else {
				auto blue1 = box(integral,
				                 r + 12,
				                 c + d * square_width,
				                 r + square_height,
				                 c + (d + 0.5) * square_width,
				                 0);
				auto blue2 = box(integral,
				                 r + 12,
				                 c + (d + 0.5) * square_width,
				                 r + square_height,
				                 c + (d + 1) * square_width,
				                 0);
				if (abs(blue1 - blue2) < 10) {
					driver_info[driver][d] = _;
				} else if (blue1 > blue2) {
					driver_info[driver][d] = L;
				} else {
					driver_info[driver][d] = E;
				}
			}
		}
	}
}

// Moves the mouse to absolute coordinates.
void mouse_goto(const Image& image, int r, int c) {
#ifdef WIN32
	INPUT input;
	input.type = INPUT_MOUSE;
	MOUSEINPUT& mouse_input = input.mi;
	mouse_input.dx = 65535.0f * float(c) / float(image.back().size());
	mouse_input.dy = 65535.0f * float(r) / float(image.size());
	mouse_input.mouseData = 0;
	mouse_input.dwFlags = MOUSEEVENTF_ABSOLUTE | MOUSEEVENTF_MOVE | MOUSEEVENTF_VIRTUALDESK;
	mouse_input.time = 0;
	mouse_input.dwExtraInfo = 0;
	SendInput(1, &input, sizeof(INPUT));
#endif
}

// Moves the mouse to the unused slots at the top.
void mouse_goto_unused_slot(
    const Image& image, const vector<pair<int, int>>& rows, int d, int line, bool early) {
	minimum_core_assert(0 <= d && d < num_days);
	minimum_core_assert(0 <= line && line < num_lines);
	constexpr int square_width = 66;
	const int square_height = rows[1].first - rows[0].first;

	int r = rows[0].first - square_height - 86 + line * 30;
	int c = rows[0].second - (num_days - d) * square_width + 18;
	if (!early) {
		c += 30;
	}
	mouse_goto(image, r, c);
}

// Moves the mouse to the slots in the driver matrix.
void mouse_goto_driver_slot(
    const Image& image, const vector<pair<int, int>>& rows, int c, int d, bool early) {
	minimum_core_assert(0 <= d && d < num_days);
	minimum_core_assert(0 <= c && c < num_drivers);
	constexpr int square_width = 66;
	const int square_height = rows[1].first - rows[0].first;

	int row = rows[0].first - square_height + c * square_height + square_height / 2;
	int col = rows[0].second - (num_days - d) * square_width + 18;
	if (!early) {
		col += 30;
	}
	mouse_goto(image, row, col);
}

// Generates a mouse click.
void mouse_click() {
#ifdef WIN32
	INPUT input;
	input.type = INPUT_MOUSE;
	MOUSEINPUT& mouse_input = input.mi;
	mouse_input.dx = 0;
	mouse_input.dy = 0;
	mouse_input.mouseData = 0;
	mouse_input.dwFlags = MOUSEEVENTF_LEFTDOWN;
	mouse_input.time = 0;
	mouse_input.dwExtraInfo = 0;
	SendInput(1, &input, sizeof(INPUT));
	Sleep(10);
	mouse_input.dwFlags = MOUSEEVENTF_LEFTUP;
	SendInput(1, &input, sizeof(INPUT));
	Sleep(500);
#endif
}

int main_program() {
	map<int, int> rows_map;

	Timer t("Parsing clipboard");
	auto image = get_clipboard();
	auto rows = get_orange_rows(image);
	auto integral = compute_integral_image(image);
	if (rows.empty()) {
		t.fail();
	} else {
		t.OK();
	}

	int _ = 0;
	int X = 1;  // Not available.
	int driver_qualified[num_drivers][num_lines] = {
	    {_, X, _},  // A
	    {X, _, _},  // B
	    {X, _, _},  // C
	    {_, X, X},  // C
	    {X, X, _},  // E
	    {_, X, X},  // F
	    {_, X, _},  // G
	    {X, _, X},  // H
	    {_, X, _},  // I
	    {_, _, X},  // J
	    {X, _, _}   // K
	};
	if (!rows.empty()) {
		Timer t("Parsing driver qualifications");
		parse_driver_qualified(integral, rows, driver_qualified, _, X);
		t.OK();
	}

	int O = 2;  // Optional day off.
	int E = 3;  // Want an early shift.
	int L = 4;  // Want a late shift.
	int driver_info[num_drivers][num_days] = {
	    {_, _, X, X, _, _, L, O, E, X, X, E, E, O},  // A
	    {L, _, _, X, X, E, O, O, L, _, X, X, _, _},  // B
	    {O, _, _, _, X, X, _, E, _, O, _, X, X, L},  // C
	    {_, _, O, O, _, X, X, _, _, E, O, _, X, X},  // D
	    {X, _, _, _, _, _, X, X, _, O, O, E, E, X},  // E
	    {X, X, L, _, E, O, _, X, X, _, _, O, _, L},  // F
	    {E, X, X, _, O, O, E, _, X, X, _, _, _, O},  // G
	    {O, E, X, X, _, O, O, _, E, X, X, _, _, _},  // H
	    {_, _, _, X, X, L, _, L, _, _, X, X, _, O},  // I
	    {O, L, L, _, X, X, _, O, O, _, E, X, X, _},  // J
	    {_, _, O, _, L, X, X, E, O, _, _, _, X, X}   // K
	};
	if (!rows.empty()) {
		Timer t("Parsing driver matrix");
		parse_driver_info(integral, rows, driver_info, _, X, O, E, L);
		t.OK();
	}

	// Print the input data.
	for (auto c : range(num_drivers)) {
		cout << char('A' + c) << ": ";
		for (auto l : range(num_lines)) {
			if (driver_qualified[c][l] == _) {
				cout << l + 1;
			} else {
				cout << " ";
			}
		}
		cout << " |";
		for (auto d : range(num_days)) {
			auto ch = driver_info[c][d];
			if (ch == _) {
				cout << ".";
			} else if (ch == X) {
				cout << "X";
			} else if (ch == O) {
				cout << "O";
			} else if (ch == L) {
				cout << "L";
			} else if (ch == E) {
				cout << "E";
			} else {
				minimum_core_assert(false, "Unknown driver info.");
			}
		}
		cout << "|\n";
	}

	IP ip;
	auto early = ip.add_boolean_cube(num_drivers, num_days, num_lines);
	auto late = ip.add_boolean_cube(num_drivers, num_days, num_lines);
	Sum objective;

	// Helper variable that tells us if a driver is working at all
	// on a particular day.
	auto working_on = ip.add_boolean_grid(num_drivers, num_days);
	// At most one shift per day.
	for (auto c : range(num_drivers)) {
		for (auto d : range(num_days)) {
			Sum day_sum = 0;
			for (auto l : range(num_lines)) {
				day_sum += early[c][d][l];
				day_sum += late[c][d][l];
			}
			ip.add_constraint(day_sum <= 1);
			ip.add_constraint(day_sum == working_on[c][d]);
		}
	}

	// Driver qualifications.
	for (auto c : range(num_drivers)) {
		for (auto l : range(num_lines)) {
			if (driver_qualified[c][l] == X) {
				for (auto d : range(num_days)) {
					ip.add_constraint(early[c][d][l] == 0);
					ip.add_constraint(late[c][d][l] == 0);
				}
			}
		}
	}

	// Driver force off.
	for (auto c : range(num_drivers)) {
		for (auto d : range(num_days)) {
			if (driver_info[c][d] == X) {
				for (auto l : range(num_lines)) {
					ip.add_constraint(early[c][d][l] == 0);
					ip.add_constraint(late[c][d][l] == 0);
				}
			}
		}
	}

	// Driver want off.
	Sum dayoff_preferences = 0;
	for (auto c : range(num_drivers)) {
		for (auto d : range(num_days)) {
			if (driver_info[c][d] == O) {
				auto slack = ip.add_boolean();
				ip.add_constraint(working_on[c][d] + slack <= 1);
				dayoff_preferences += driver_want_off * slack;
			}
		}
	}

	// Driver shift preferences.
	Sum shift_preferences = 0;
	for (auto c : range(num_drivers)) {
		for (auto d : range(num_days)) {
			if (driver_info[c][d] == E) {
				Sum day_sum = 0;
				for (auto l : range(num_lines)) {
					day_sum += early[c][d][l];
				}
				shift_preferences += shift_preference_granted * day_sum;
			}

			if (driver_info[c][d] == L) {
				Sum day_sum = 0;
				for (auto l : range(num_lines)) {
					day_sum += late[c][d][l];
				}
				shift_preferences += shift_preference_granted * day_sum;
			}
		}
	}

	// Long rests.
	//
	// Count the number of three days off followed by a working day or the
	// end of the period. That should be equal to the number of long rests.
	PseudoBoolean long_rests = Sum(0);
	for (auto c : range(num_drivers)) {
		for (auto d : range(num_days - 2)) {
			PseudoBoolean range_product = Sum(1);
			for (auto d2 : range(d, d + 3)) {
				range_product *= (1 - working_on[c][d2]);
			}
			if (d + 3 < num_days) {
				range_product *= working_on[c][d + 3];
			}
			long_rests += long_rest_bonus * range_product;
		}
	}

	// Cover the shifts.
	for (auto d : range(num_days)) {
		for (auto l : range(num_lines)) {
			Sum early_sum = 0;
			Sum late_sum = 0;
			for (auto c : range(num_drivers)) {
				early_sum += early[c][d][l];
				late_sum += late[c][d][l];
			}
			auto early_slack = ip.add_boolean();
			ip.add_constraint(early_sum + early_slack == 1);
			objective += unassigned_shift_cost * early_slack;

			auto late_slack = ip.add_boolean();
			ip.add_constraint(late_sum + late_slack == 1);
			objective += unassigned_shift_cost * late_slack;
		}
	}

	// Late shift followed by early.
	for (auto c : range(num_drivers)) {
		for (auto d : range(num_days - 1)) {
			Sum late_sum = 0;
			Sum early_sum = 0;
			for (auto l : range(num_lines)) {
				late_sum += late[c][d][l];
				early_sum += early[c][d + 1][l];
			}
			auto slack = ip.add_boolean();
			ip.add_constraint(late_sum + early_sum - slack <= 1);
			objective += late_shift_plus_early * slack;
		}
	}

	// Late shift after three late shifts.
	for (auto c : range(num_drivers)) {
		for (auto d : range(num_days - max_late_shifts)) {
			Sum late_sum = 0;
			for (auto d2 : range(d, d + max_late_shifts + 1)) {
				for (auto l : range(num_lines)) {
					late_sum += late[c][d2][l];
				}
			}
			auto slack = ip.add_boolean();
			ip.add_constraint(late_sum - slack <= max_late_shifts);
			objective += four_consecutive_late_shifts * slack;
		}
	}

	// Late shifts should be four.
	Sum deviation_target_late_shifts = 0;
	for (auto c : range(num_drivers)) {
		Sum late_sum = 0;
		for (auto d : range(num_days)) {
			for (auto l : range(num_lines)) {
				late_sum += late[c][d][l];
			}
		}
		auto slack1 = ip.add_variable(IP::Integer);
		auto slack2 = ip.add_variable(IP::Integer);
		ip.add_constraint(slack1 >= 0);
		ip.add_constraint(slack2 >= 0);
		ip.add_constraint(late_sum + slack1 - slack2 == target_late_shifts);
		deviation_target_late_shifts += shift_assigned_not_four * slack1;
		deviation_target_late_shifts += shift_assigned_not_four * slack2;
	}

	objective += shift_preferences;
	objective += dayoff_preferences;
	objective += deviation_target_late_shifts;
	ip.add_objective(objective);
	ip.add_pseudoboolean_objective(long_rests);

	ip.linearize_pseudoboolean_terms();

	IPSolver solver;
	minimum_core_assert(solver.solutions(&ip)->get());

	for (auto c : range(num_drivers)) {
		cout << char('A' + c) << ": |";
		int num_late = 0;
		for (auto d : range(num_days)) {
			bool found = false;
			for (auto l : range(num_lines)) {
				if (early[c][d][l].bool_value()) {
					cout << l + 1;
					minimum_core_assert(!found);
					found = true;
				}
			}
			if (driver_info[c][d] == X) {
				minimum_core_assert(!found);
				cout << ".";
			} else if (!found) {
				cout << " ";
			}

			found = false;
			for (auto l : range(num_lines)) {
				if (late[c][d][l].bool_value()) {
					cout << l + 1;
					minimum_core_assert(!found);
					found = true;
					num_late += 1;
				}
			}
			if (driver_info[c][d] == X) {
				minimum_core_assert(!found);
				cout << ".";
			} else if (!found) {
				cout << " ";
			}
			cout << "|";
		}
		cout << "  " << num_late << " late.\n";
	}
	cout << "Objective: " << objective.value() + long_rests.value() << "\n";
	cout << "Shift preferences: " << shift_preferences.value() << "\n";
	cout << "Dayoff preferences: " << dayoff_preferences.value() << "\n";
	cout << "Long rests: " << long_rests.value() << "\n";
	cout << "Deviation target late shifts: " << deviation_target_late_shifts.value() << "\n";

	if (!rows.empty()) {
		cout << "\n";
		Timer t("Auto-assigning using mouse");
		for (auto c : range(num_drivers)) {
			for (auto d : range(num_days)) {
				for (auto l : range(num_lines)) {
					if (early[c][d][l].bool_value()) {
						mouse_goto_unused_slot(image, rows, d, l, true);
						mouse_click();
						mouse_goto_driver_slot(image, rows, c, d, true);
						mouse_click();
					}
				}

				for (auto l : range(num_lines)) {
					if (late[c][d][l].bool_value()) {
						mouse_goto_unused_slot(image, rows, d, l, false);
						mouse_click();
						mouse_goto_driver_slot(image, rows, c, d, false);
						mouse_click();
					}
				}
			}
		}
		t.OK();
	}

	return 0;
}

int main() {
	try {
		return main_program();
	} catch (std::exception& err) {
		std::cerr << "ERROR: " << err.what() << std::endl;
		return 1;
	}
}
