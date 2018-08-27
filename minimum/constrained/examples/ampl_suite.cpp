#include <fstream>
#include <iostream>
#include <streambuf>
#include <string>

#include <gflags/gflags.h>

#include <minimum/constrained/ampl.h>
#include <minimum/constrained/successive.h>
#include <minimum/core/check.h>
#include <minimum/core/color.h>
#include <minimum/core/enum.h>
#include <minimum/core/main.h>
#include <minimum/core/numeric.h>
#include <minimum/core/range.h>
#include <minimum/core/string.h>
#include <minimum/nonlinear/constrained_function.h>
#include <minimum/nonlinear/data/util.h>
#include <minimum/nonlinear/solver.h>
using namespace std;
using minimum::core::check;
using minimum::core::range;
using minimum::core::Timer;
using minimum::nonlinear::data::get_directory;

DEFINE_int32(max_variables, 100, "Maximum number of variables to process.");

namespace {
string get_file(const string& filename) {
	ifstream file(filename, ios::binary);
	check(!!file, "Could not open ", filename);
	return {istreambuf_iterator<char>(file), istreambuf_iterator<char>()};
}

#ifdef ERROR
#undef ERROR
#endif
MAKE_ENUM(Status, OK, SKIP, PARSE_ERROR, ERROR, WRONG_RESULT);

Status test_problem(std::string_view data, double expected_objective) {
	minimum::constrained::AmplParser parser(data);
	std::optional<minimum::constrained::AmplProblem> problem;
	try {
		Timer t("Parsing");
		problem.emplace(parser.parse_ampl());
		t.OK();
	} catch (std::exception& err) {
		std::cerr << "\n" << minimum::core::RED << err.what() << minimum::core::NORMAL << "\n";
		return Status::PARSE_ERROR;
	}
	std::cerr << problem->get_variables().size() << " variables.\n";
	if (problem->get_variables().size() > FLAGS_max_variables && problem->has_constraints()) {
		return Status::SKIP;
	}

	if (problem->has_constraints()) {
		Timer t("Creating constrained function.");
		minimum::nonlinear::ConstrainedFunction f;
		problem->create_function(&f);
		t.OK();

		minimum::constrained::SuccessiveLinearProgrammingSolver solver;
		solver.function_improvement_tolerance = 1e-6;
		solver.log_function = [](std::string_view s) { std::cerr << s << "\n"; };

		auto result = solver.solve(&f);
		for (auto i : range(problem->get_variables().size())) {
			std::cerr << problem->get_variables()[i].name << "=" << problem->get_values()[i]
			          << "\n";
		}
		minimum_core_assert(result
		                    == minimum::constrained::SuccessiveLinearProgrammingSolver::Result::OK);
		minimum_core_assert(f.is_feasible());
		if (!std::isnan(expected_objective)) {
			auto objective = f.objective().evaluate();
			std::cerr << "Expected objective: " << expected_objective << "\n";
			std::cerr << "Objective: " << objective << "\n";
			if (minimum::core::relative_error(expected_objective, objective) > 1e-4) {
				return Status::WRONG_RESULT;
			}
		}
	} else {
		Timer t("Creating unconstrained function.");
		minimum::nonlinear::Function f;
		problem->create_function(&f);
		t.OK();

		minimum::nonlinear::LBFGSBSolver solver;
		solver.maximum_iterations = 1000;
		solver.gradient_tolerance = 1e-4;
		solver.log_function = [](std::string_view s) { std::cerr << s << "\n"; };
		minimum::nonlinear::SolverResults results;
		solver.solve(f, &results);
		std::cerr << results << "\n";
		minimum_core_assert(results.exit_success());

		if (!std::isnan(expected_objective)) {
			auto objective = f.evaluate();
			std::cerr << "Expected objective: " << expected_objective << "\n";
			std::cerr << "Objective: " << objective << "\n";
			if (minimum::core::relative_error(expected_objective, objective) > 1e-4) {
				return Status::WRONG_RESULT;
			}
		}
	}
	return Status::OK;
}

Status test_global_problem(std::string problem_name,
                           double expected_objective = std::numeric_limits<double>::quiet_NaN()) {
	std::string data = get_file(get_directory() + "/coconut/global/" + problem_name + ".mod");
	try {
		std::cerr << "\n"
		          << minimum::core::GREEN << problem_name << minimum::core::NORMAL << "\n\n";
		return test_problem(data, expected_objective);
	} catch (std::exception& err) {
		std::cerr << "\n" << minimum::core::RED << err.what() << minimum::core::NORMAL << "\n";
		return Status::ERROR;
	}
}
}  // namespace

int main_program(int num_args, char* args[]) {
	std::map<string, Status> results;
	const double nan = numeric_limits<double>::quiet_NaN();
	std::map<string, double> problems = {
	    {"abel", nan},          {"alkyl", nan},         {"bearing", nan},
	    {"camcge", nan},        {"camcns", nan},        {"camshape100", nan},
	    {"camshape200", nan},   {"camshape400", nan},   {"catmix100", nan},
	    {"catmix200", nan},     {"chain100", nan},      {"chain200", nan},
	    {"chain25", nan},       {"chain400", nan},      {"chain50", nan},
	    {"chakra", nan},        {"chance", nan},        {"chem", nan},
	    {"chenery", nan},       {"circle", nan},        {"demo7", nan},
	    {"dispatch", nan},      {"elec25", nan},        {"elec50", nan},
	    {"etamac", nan},        {"ex14_1_1", nan},      {"ex14_1_2", nan},
	    {"ex14_1_3", nan},      {"ex14_1_4", nan},      {"ex14_1_5", nan},
	    {"ex14_1_6", nan},      {"ex14_1_7", nan},      {"ex14_1_8", nan},
	    {"ex14_1_9", nan},      {"ex14_2_1", nan},      {"ex14_2_2", nan},
	    {"ex14_2_3", nan},      {"ex14_2_4", nan},      {"ex14_2_5", nan},
	    {"ex14_2_6", nan},      {"ex14_2_7", nan},      {"ex14_2_8", nan},
	    {"ex14_2_9", nan},      {"ex2_1_10", nan},      {"ex2_1_1", nan},
	    {"ex2_1_2", nan},       {"ex2_1_3", nan},       {"ex2_1_4", nan},
	    {"ex2_1_5", nan},       {"ex2_1_6", nan},       {"ex2_1_7", nan},
	    {"ex2_1_8", nan},       {"ex2_1_9", nan},       {"ex3_1_1", nan},
	    {"ex3_1_2", nan},       {"ex3_1_3", nan},       {"ex3_1_4", nan},
	    {"ex4_1_1", nan},       {"ex4_1_2", nan},       {"ex4_1_3", nan},
	    {"ex4_1_4", nan},       {"ex4_1_5", nan},       {"ex4_1_6", nan},
	    {"ex4_1_7", nan},       {"ex4_1_8", nan},       {"ex4_1_9", nan},
	    {"ex5_2_2_case1", nan}, {"ex5_2_2_case2", nan}, {"ex5_2_2_case3", nan},
	    {"ex5_2_4", nan},       {"ex5_2_5", nan},       {"ex5_3_2", nan},
	    {"ex5_3_3", nan},       {"ex5_4_2", nan},       {"ex5_4_3", nan},
	    {"ex5_4_4", nan},       {"ex6_1_1", nan},       {"ex6_1_2", nan},
	    {"ex6_1_3", nan},       {"ex6_1_4", nan},       {"ex6_2_10", nan},
	    {"ex6_2_11", nan},      {"ex6_2_12", nan},      {"ex6_2_13", nan},
	    {"ex6_2_14", nan},      {"ex6_2_5", nan},       {"ex6_2_6", nan},
	    {"ex6_2_7", nan},       {"ex6_2_8", nan},       {"ex6_2_9", nan},
	    {"ex7_2_10", nan},      {"ex7_2_1", nan},       {"ex7_2_2", nan},
	    {"ex7_2_3", nan},       {"ex7_2_4", nan},       {"ex7_2_5", nan},
	    {"ex7_2_6", nan},       {"ex7_2_7", nan},       {"ex7_2_8", nan},
	    {"ex7_2_9", nan},       {"ex7_3_1", nan},       {"ex7_3_2", nan},
	    {"ex7_3_3", nan},       {"ex7_3_4", nan},       {"ex7_3_5", nan},
	    {"ex7_3_6", nan},       {"ex8_1_1", nan},       {"ex8_1_2", nan},
	    {"ex8_1_3", nan},       {"ex8_1_4", nan},       {"ex8_1_5", nan},
	    {"ex8_1_6", nan},       {"ex8_1_7", nan},       {"ex8_1_8", nan},
	    {"ex8_2_1", nan},       {"ex8_2_4", nan},       {"ex8_3_10", nan},
	    {"ex8_3_11", nan},      {"ex8_3_12", nan},      {"ex8_3_13", nan},
	    {"ex8_3_14", nan},      {"ex8_3_1", nan},       {"ex8_3_2", nan},
	    {"ex8_3_3", nan},       {"ex8_3_4", nan},       {"ex8_3_5", nan},
	    {"ex8_3_6", nan},       {"ex8_3_7", nan},       {"ex8_3_8", nan},
	    {"ex8_3_9", nan},       {"ex8_4_1", nan},       {"ex8_4_2", nan},
	    {"ex8_4_3", nan},       {"ex8_4_4", nan},       {"ex8_4_5", nan},
	    {"ex8_4_6", nan},       {"ex8_4_7", nan},       {"ex8_4_8", nan},
	    {"ex8_5_1", nan},       {"ex8_5_2", nan},       {"ex8_5_3", nan},
	    {"ex8_5_4", nan},       {"ex8_5_5", nan},       {"ex8_5_6", nan},
	    {"ex8_6_1", nan},       {"ex8_6_2", nan},       {"ex9_1_10", nan},
	    {"ex9_1_1", nan},       {"ex9_1_2", nan},       {"ex9_1_3", nan},
	    {"ex9_1_4", nan},       {"ex9_1_5", nan},       {"ex9_1_6", nan},
	    {"ex9_1_7", nan},       {"ex9_1_8", nan},       {"ex9_1_9", nan},
	    {"ex9_2_1", nan},       {"ex9_2_2", nan},       {"ex9_2_3", nan},
	    {"ex9_2_4", nan},       {"ex9_2_5", nan},       {"ex9_2_6", nan},
	    {"ex9_2_7", nan},       {"ex9_2_8", nan},       {"gancns", nan},
	    {"gancnsx", nan},       {"ganges", nan},        {"gangesx", nan},
	    {"glider100", nan},     {"glider50", nan},      {"gtm", nan},
	    {"harker", nan},        {"haverly", nan},       {"hhfair", nan},
	    {"himmel11", nan},      {"himmel16", nan},      {"house", nan},
	    {"hydro", nan},         {"immun", nan},         {"korcge", nan},
	    {"korcns", nan},        {"launch", nan},        {"least", nan},
	    {"like", nan},          {"linear", nan},        {"lnts100", nan},
	    {"lnts200", nan},       {"lnts50", nan},        {"meanvar", nan},
	    {"mhw4d", nan},         {"minlphi", nan},       {"nemhaus", nan},
	    {"otpop", nan},         {"pindyck", nan},       {"pollut", nan},
	    {"polygon25", nan},     {"polygon50", nan},     {"process", nan},
	    {"prolog", nan},        {"qp1", nan},           {"qp2", nan},
	    {"qp3", nan},           {"qp4", nan},           {"qp5", nan},
	    {"ramsey", nan},        {"rbrock", nan},        {"robot100", nan},
	    {"robot50", nan},       {"rocket100", nan},     {"rocket200", nan},
	    {"rocket50", nan},      {"sambal", nan},        {"sample", nan},
	    {"ship", nan},          {"srcpm", nan},         {"turkey", nan},
	    {"wall", nan},          {"water", nan},         {"weapons", nan},
	    {"worst", nan},
	};
	problems["abel"] = 225.19460000;
	problems["alkyl"] = -1.76500000;
	problems["camshape100"] = -4.28410000;
	problems["camshape200"] = -4.27850000;
	problems["catmix100"] = -0.04810000;
	problems["catmix200"] = -0.04810000;
	problems["chain100"] = 5.06980000;
	problems["chain200"] = 5.06890000;
	problems["chain400"] = 5.06860000;
	problems["chain50"] = 5.07230000;
	problems["chakra"] = -179.13360000;
	problems["chance"] = 29.89440000;
	problems["chem"] = -47.70650000;
	problems["chenery"] = -1058.91990000;
	problems["circle"] = 4.57420000;
	problems["demo7"] = -1589042.38620000;
	problems["dispatch"] = 3155.28790000;
	problems["elec25"] = 243.81270000;
	problems["elec50"] = 1055.18200000;
	problems["etamac"] = -15.29470000;
	problems["ex14_1_1"] = 0.00000000;
	problems["ex14_1_2"] = 0.00000000;
	problems["ex14_1_3"] = 0.00000000;
	problems["ex14_1_4"] = 0.00000000;
	problems["ex14_1_5"] = 0.00000000;
	problems["ex14_1_6"] = 0.00000000;
	problems["ex14_1_7"] = 0.00000000;
	problems["ex14_1_8"] = 0.00000000;
	problems["ex14_1_9"] = 0.00000000;
	problems["ex14_2_1"] = 0.00000000;
	problems["ex14_2_2"] = 0.00000000;
	problems["ex14_2_3"] = 0.00000000;
	problems["ex14_2_4"] = 0.00000000;
	problems["ex14_2_5"] = 0.00000000;
	problems["ex14_2_6"] = 0.00000000;
	problems["ex14_2_7"] = 0.00000000;
	problems["ex14_2_8"] = 0.00000000;
	problems["ex14_2_9"] = 0.00000000;
	problems["ex2_1_1"] = -17.00000000;
	problems["ex2_1_10"] = 49318.01800000;
	problems["ex2_1_2"] = -213.00000000;
	problems["ex2_1_3"] = -15.00000000;
	problems["ex2_1_4"] = -11.00000000;
	problems["ex2_1_5"] = -268.01460000;
	problems["ex2_1_6"] = -39.00000000;
	problems["ex2_1_7"] = -4150.41010000;
	problems["ex2_1_8"] = 15639.00000000;
	problems["ex2_1_9"] = -0.37500000;
	problems["ex3_1_1"] = 7049.20830000;
	problems["ex3_1_2"] = -30665.54000000;
	problems["ex3_1_3"] = -310.00000000;
	problems["ex3_1_4"] = -4.00000000;
	problems["ex4_1_1"] = -7.48730000;
	problems["ex4_1_2"] = -663.50010000;
	problems["ex4_1_3"] = -443.67170000;
	problems["ex4_1_4"] = 0.00000000;
	problems["ex4_1_5"] = 0.00000000;
	problems["ex4_1_6"] = 7.00000000;
	problems["ex4_1_7"] = -7.50000000;
	problems["ex4_1_8"] = -16.73890000;
	problems["ex4_1_9"] = -5.50800000;
	problems["ex5_2_2_case1"] = -400.00000000;
	problems["ex5_2_2_case2"] = -600.00000000;
	problems["ex5_2_2_case3"] = -750.00000000;
	problems["ex5_2_4"] = -450.00000000;
	problems["ex5_2_5"] = -3500.00010000;
	problems["ex5_3_2"] = 1.86420000;
	problems["ex5_3_3"] = 3.23400000;
	problems["ex5_4_2"] = 7512.22590000;
	problems["ex5_4_3"] = 4845.46200000;
	problems["ex5_4_4"] = 10077.77540000;
	problems["ex6_1_1"] = -0.02020000;
	problems["ex6_1_2"] = -0.03250000;
	problems["ex6_1_3"] = -0.35250000;
	problems["ex6_1_4"] = -0.29450000;
	problems["ex6_2_10"] = -3.05200000;
	problems["ex6_2_11"] = 0.00000000;
	problems["ex6_2_12"] = 0.28920000;
	problems["ex6_2_13"] = -0.21620000;
	problems["ex6_2_14"] = -0.69540000;
	problems["ex6_2_5"] = -70.75210000;
	problems["ex6_2_6"] = 0.00000000;
	problems["ex6_2_7"] = -0.16080000;
	problems["ex6_2_8"] = -0.02700000;
	problems["ex6_2_9"] = -0.03410000;
	problems["ex7_2_1"] = 1227.18960000;
	problems["ex7_2_2"] = -0.38880000;
	problems["ex7_2_3"] = 7049.21810000;
	problems["ex7_2_4"] = 3.91800000;
	problems["ex7_2_5"] = 10122.48280000;
	problems["ex7_2_6"] = -83.24990000;
	problems["ex7_2_7"] = -5.73990000;
	problems["ex7_2_8"] = -6.08200000;
	problems["ex7_2_9"] = 1.14360000;
	problems["ex7_2_10"] = 0.10000000;
	problems["ex7_3_1"] = 0.34170000;
	problems["ex7_3_2"] = 1.08990000;
	problems["ex7_3_3"] = 0.81750000;
	problems["ex7_3_4"] = 6.27460000;
	problems["ex7_3_5"] = 1.20360000;
	problems["ex7_3_6"] = 0.00000000;
	problems["ex8_1_1"] = -2.02180000;
	problems["ex8_1_2"] = -1.07090000;
	problems["ex8_1_3"] = 3.00000000;
	problems["ex8_1_4"] = 0.00000000;
	problems["ex8_1_5"] = -1.03160000;
	problems["ex8_1_6"] = -10.08600000;
	problems["ex8_1_7"] = 0.02930000;
	problems["ex8_1_8"] = -0.38880000;
	problems["ex8_2_1"] = -979.17830000;
	problems["ex8_2_4"] = -58307398.39170000;
	problems["ex8_3_1"] = -0.79370000;
	problems["ex8_3_10"] = -0.88450000;
	problems["ex8_3_12"] = -0.58250000;
	problems["ex8_3_13"] = -42.39040000;
	problems["ex8_3_14"] = -2.41240000;
	problems["ex8_3_2"] = -0.41230000;
	problems["ex8_3_3"] = -0.41660000;
	problems["ex8_3_4"] = -3.58000000;
	problems["ex8_3_5"] = -0.06880000;
	problems["ex8_3_6"] = -0.50000000;
	problems["ex8_3_7"] = -1.23260000;
	problems["ex8_3_8"] = -3.24680000;
	problems["ex8_3_9"] = -0.76300000;
	problems["ex8_4_1"] = 0.61860000;
	problems["ex8_4_2"] = 0.48520000;
	problems["ex8_4_3"] = 0.00460000;
	problems["ex8_4_4"] = 0.21250000;
	problems["ex8_4_5"] = 0.00030000;
	problems["ex8_4_6"] = 0.00110000;
	problems["ex8_4_7"] = 28.89800000;
	problems["ex8_4_8"] = 3.32180000;
	problems["ex8_5_1"] = 0.00000000;
	problems["ex8_5_2"] = -0.00020000;
	problems["ex8_5_3"] = -0.00420000;
	problems["ex8_5_4"] = -0.00040000;
	problems["ex8_5_5"] = -0.01080000;
	problems["ex8_5_6"] = -0.00120000;
	problems["ex8_6_1"] = -28.42250000;
	problems["ex8_6_2"] = -31.88860000;
	problems["ex9_1_1"] = -13.00000000;
	problems["ex9_1_10"] = -3.25000000;
	problems["ex9_1_2"] = -16.00000000;
	problems["ex9_1_3"] = -58.00000000;
	problems["ex9_1_4"] = -37.00000000;
	problems["ex9_1_5"] = -1.00000000;
	problems["ex9_1_6"] = -52.00000000;
	problems["ex9_1_7"] = -50.00000000;
	problems["ex9_1_8"] = -3.25000000;
	problems["ex9_1_9"] = 2.00000000;
	problems["ex9_2_1"] = 17.00000000;
	problems["ex9_2_2"] = 99.99950000;
	problems["ex9_2_3"] = 0.00000000;
	problems["ex9_2_4"] = 0.50000000;
	problems["ex9_2_5"] = 5.00000000;
	problems["ex9_2_6"] = -1.00000000;
	problems["ex9_2_7"] = 17.00000000;
	problems["ex9_2_8"] = 1.50000000;
	problems["gtm"] = 543.56510000;
	problems["harker"] = -986.51350000;
	problems["haverly"] = -400.00000000;
	problems["hhfair"] = -87.15900000;
	problems["himmel11"] = -30665.54000000;
	problems["himmel16"] = -0.86600000;
	problems["house"] = -4500.00000000;
	problems["hydro"] = 4366944.00000000;
	problems["immun"] = 0.00000000;
	problems["launch"] = 2257.79760000;
	problems["least"] = 14085.13980000;
	problems["linear"] = 89.00000000;
	problems["lnts100"] = 0.55460000;
	problems["lnts50"] = 0.55470000;
	problems["meanvar"] = 5.24340000;
	problems["mhw4d"] = 0.02930000;
	problems["minlphi"] = 582.23610000;
	problems["nemhaus"] = 31.00000000;
	problems["otpop"] = 0.00000000;
	problems["pindyck"] = -1612.17830000;
	problems["pollut"] = -5353268.62880000;
	problems["polygon25"] = -0.77220000;
	problems["polygon50"] = -0.77570000;
	problems["process"] = -1161.33660000;
	problems["prolog"] = 0.00000000;
	problems["qp1"] = 0.00080000;
	problems["qp2"] = 0.00080000;
	problems["qp3"] = 0.00080000;
	problems["qp4"] = 0.00080000;
	problems["qp5"] = 0.43150000;
	problems["ramsey"] = -2.48750000;
	problems["rbrock"] = 0.00000000;
	problems["robot50"] = 9.14690000;
	problems["rocket100"] = -1.01260000;
	problems["rocket50"] = -1.01280000;
	problems["sambal"] = 3.96820000;
	problems["sample"] = 726.63670000;
	problems["ship"] = 0.00000000;
	problems["srcpm"] = 2109.78180000;
	problems["turkey"] = -29330.15800000;
	problems["wall"] = -1.00000000;
	problems["water"] = 906.35190000;
	for (auto& entry : problems) {
		if (num_args != 2 || args[1] == entry.first) {
			results[entry.first] = test_global_problem(entry.first, entry.second);
		}
	}
	std::cout << "\n\n";
	std::map<Status, minimum::core::Color> colors;
	colors.emplace(Status::OK, minimum::core::GREEN);
	colors.emplace(Status::SKIP, minimum::core::YELLOW);
	colors.emplace(Status::WRONG_RESULT, minimum::core::YELLOW);
	std::map<Status, int> result_counts;
	for (auto& entry : results) {
		std::cout << entry.first << ": ";
		minimum::core::Color color = minimum::core::RED;
		auto itr = colors.find(entry.second);
		if (itr != colors.end()) {
			color = itr->second;
		}
		std::cerr << color << to_string(entry.second) << "\n" << minimum::core::NORMAL;
		result_counts[entry.second] += 1;
	}
	std::cerr << "\n";
	for (auto& entry : result_counts) {
		minimum::core::Color color = minimum::core::RED;
		auto itr = colors.find(entry.first);
		if (itr != colors.end()) {
			color = itr->second;
		}
		std::cerr << color << to_string(entry.first) << minimum::core::NORMAL << ": "
		          << entry.second << "\n";
	}
	return 0;
}

int main(int argc, char* argv[]) { return main_runner(main_program, argc, argv); }
