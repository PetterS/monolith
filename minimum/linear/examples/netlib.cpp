// Solution values to NETLIB problems:
// http://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=AE5378994792586ABA2A0326E95C2415?doi=10.1.1.88.7282&rep=rep1&type=pdf
// and
// http://www.netlib.org/lp/data/readme

#ifndef _MSC_VER
int main() {}  // Need to link with -lstdc++fs with gcc.
#else

#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
using namespace std;

#include <gflags/gflags.h>

#include <minimum/core/check.h>
#include <minimum/core/numeric.h>
#include <minimum/linear/data/util.h>
#include <minimum/linear/first_order_solver.h>
#include <minimum/linear/glpk.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/scs.h>
using namespace minimum::core;
using namespace minimum::linear;

DEFINE_bool(verify_with_simplex,
            false,
            "Verfies the problem solutions with simplex (for testing the MPS reader).");

DEFINE_string(ip_proto_directory,
              "",
              "If set, searches for and opens IP proto files instead of Netlib MPS. Default: not "
              "set. Implies verify_with_simplex.");

DEFINE_int64(max_iterations, 100'000, "Maximum number of iterations to use.");

std::vector<std::string> problem_names = {
    "25fv47",   "80bau3b",  "adlittle", "afiro",    "agg",      "agg2",     "agg3",     "bandm",
    "beaconfd", "blend",    "bnl1",     "bnl2",     "boeing1",  "boeing2",  "bore3d",   "brandy",
    "capri",    "cycle",    "czprob",   "d2q06c",   "d6cube",   "degen2",   "degen3",   "dfl001",
    "e226",     "etamacro", "fffff800", "finnis",   "fit1d",    "fit1p",    "fit2d",    "fit2p",
    "forplan",  "ganges",   "gfrd-pnc", "greenbea", "greenbeb", "grow15",   "grow22",   "grow7",
    "israel",   "kb2",      "lotfi",    "maros-r7", "maros",    "modszk1",  "nesm",     "perold",
    "pilot-ja", "pilot-we", "pilot",    "pilot4",   "pilot87",  "pilotnov", "recipelp", "sc105",
    "sc205",    "sc50a",    "sc50b",    "scagr25",  "scagr7",   "scfxm1",   "scfxm2",   "scfxm3",
    "scorpion", "scrs8",    "scsd1",    "scsd6",    "scsd8",    "sctap1",   "sctap2",   "sctap3",
    "seba",     "share1b",  "share2b",  "shell",    "ship04l",  "ship04s",  "ship08l",  "ship08s",
    "ship12l",  "ship12s",  "sierra",   "stair",    "standata", "standmps", "stocfor1", "stocfor2",
    "stocfor3", "truss",    "tuff",     "vtp-base", "wood1p",   "woodw",    "qap8",     "qap12",
    "qap15"};

std::vector<std::pair<double, int>> problem_solutions = {{0.55018458882867447945812325883916, 4},
                                                         {0.98722419240909025757911711738799, 6},
                                                         {0.22549496316238038228101176621492, 6},
                                                         {-0.46475314285714285714285714285714, 3},
                                                         {-0.35991767286576506712640824319636, 8},
                                                         {-0.20239252355977109024317661926133, 8},
                                                         {0.10312115935089225579061058796215, 8},
                                                         {-0.15862801845012064052174123768736, 3},
                                                         {0.335924858072, 5},
                                                         {-0.30812149845828220173774356124984, 2},
                                                         {0.19776295615228892439564398331821, 4},
                                                         {0.1811236540358545170448413697691, 4},
                                                         {-0.3352135675071266218429697314682, 3},
                                                         {-0.31501872801520287870462195913263, 3},
                                                         {0.13730803942084927215581987251301, 4},
                                                         {0.15185098964881283835426751550618, 4},
                                                         {0.26900129137681610087717280693754, 4},
                                                         {-0.52263930248941017172447233836217, 1},
                                                         {0.21851966988565774858951155947191, 7},
                                                         {0.12278421081418945895739128812392, 6},
                                                         {0.31549166666666666666666666666667, 3},
                                                         {-0.1435178, 4},
                                                         {-0.987294, 3},
                                                         {0.11266396046671392202377652175477, 8},
                                                         {-0.18751929066370549102605687681285, 2},
                                                         {-0.7557152333749133350792583667773, 3},
                                                         {0.555679564817496376532864378969, 6},
                                                         {0.17279106559561159432297900375543, 6},
                                                         {-0.91463780924209269467749025024617, 4},
                                                         {0.91463780924209269467749025024617, 4},
                                                         {-0.68464293293832069575943518435837, 5},
                                                         {0.68464293293832069575943518435837, 5},
                                                         {-0.66421896127220457481235119701692, 3},
                                                         {-0.10958573612927811623444890424357, 6},
                                                         {0.69022359995488088295415596232193, 7},
                                                         {-0.72555248129845987457557870574845, 8},
                                                         {-0.43022602612065867539213672544432, 7},
                                                         {-0.10687094129357533671604040930313, 9},
                                                         {-0.16083433648256296718456039982613, 9},
                                                         {-0.47787811814711502616766956242865, 8},
                                                         {-0.89664482186304572966200464196045, 6},
                                                         {-0.17499001299062057129526866493726, 4},
                                                         {-0.2526470606188, 2},
                                                         {0.14971851664796437907337543903552, 7},
                                                         {-0.58063743701125895401208534974734, 5},
                                                         {0.32061972906431580494333823530763, 3},
                                                         {0.14076036487562728337980641375835, 8},
                                                         {-0.93807552782351606734833270386511, 4},
                                                         {-0.61131364655813432748848620538024, 4},
                                                         {-0.27201075328449639629439185412556, 7},
                                                         {-0.55748972928406818073034256636894, 3},
                                                         {-0.25811392588838886745830997266797, 4},
                                                         {0.3017103473331105277216600201995, 3},
                                                         {-0.44972761882188711430996211783943, 4},
                                                         {-0.266616, 3},
                                                         {-0.52202061211707248062628010857689, 2},
                                                         {-0.52202061211707248062628010857689, 2},
                                                         {-0.64575077058564509026860413914575, 2},
                                                         {-0.7, 2},
                                                         {-0.14753433060768523167790925075974, 8},
                                                         {-0.2331389824330984, 7},
                                                         {0.18416759028348943683579089143655, 5},
                                                         {0.36660261564998812956939504848248, 5},
                                                         {0.54901254549751454623888724925519, 5},
                                                         {0.18781248227381066296479411763586, 4},
                                                         {0.90429695380079143579923107948844, 3},
                                                         {0.86666666743333647292533502995263, 1},
                                                         {0.50500000077144345289985489081569, 2},
                                                         {0.90499999992546440094445846417654, 3},
                                                         {0.141225, 4},
                                                         {0.17248071428571428571428571428571, 4},
                                                         {0.1424, 4},
                                                         {0.157116, 5},
                                                         {-0.7658931857918568112797274346007, 5},
                                                         {-0.41573224074141948654519910873841, 3},
                                                         {0.1208825346, 10},
                                                         {0.17933245379703557625562556255626, 7},
                                                         {0.17987147004453921719954082726049, 7},
                                                         {0.19090552113891315179803819998425, 7},
                                                         {0.1920098210534619710172695474693, 7},
                                                         {0.1470187919329264822702175543886, 7},
                                                         {0.14892361344061291565686421605401, 7},
                                                         {0.1539436218363193, 8},
                                                         {-0.25126695119296330352803637106304, 3},
                                                         {0.12576995, 4},
                                                         {0.14060175, 4},
                                                         {-0.41131976219436406065682760731514, 5},
                                                         {-0.39024408537882029604587908772433, 5},
                                                         {-0.39976783943649587403509204700686, 5},
                                                         {0.45881584718561673163624258476816, 6},
                                                         {0.29214776509361284488226309994302, 0},
                                                         {0.1298314624613613657395984384889, 6},
                                                         {0.14429024115734092400010936668043, 1},
                                                         {0.1304476333084229269005552085566, 1},
                                                         {2.0350000000E+02, 0},
                                                         {5.2289435056E+02, 0},
                                                         {1.0409940410E+03, 0}};

struct Statistics {
	Statistics(string name_) : name(move(name_)) {}

	int problems_OK = 0;
	int problems_nonoptimal = 0;
	int problems_infeasible = 0;

	string name;

	void solve_and_update(double ground_truth, IP& ip, Solver* solver, ostream& log) {
		auto solutions = solver->solutions(&ip);
		solutions->get();
		bool feasible = ip.is_feasible(1e-2);
		if (!feasible) {
			cout << name << " infeasible.\n\n";
			log << name << " INFEASIBLE\n";
			problems_infeasible++;
		} else {
			double calulated_first_order_objective = ip.get_entire_objective();
			if (!is_equal(calulated_first_order_objective, ground_truth, 1e-2)) {
				cout << name << " could not find correct objective.\n";
				cout << "-- Ground truth : " << ground_truth << "\n";
				cout << "-- Calculated   : " << calulated_first_order_objective << "\n\n";

				log << name << " NONOPTIMAL\n";
				problems_nonoptimal++;
			} else {
				log << name << " OK\n";
				problems_OK++;
			}
		}

		log << solutions->log() << endl;
	}

	void print() const {
		cout << name << " problems solved    : " << problems_OK << "\n";
		cout << name << " problems feasible  : " << problems_nonoptimal << "\n";
		cout << name << " problems infeasible: " << problems_infeasible << "\n\n";
	}
};

Statistics primal_dual_stats("Primal-dual");
Statistics primal_dual_org_stats("Primal-dual-org");
Statistics scs_stats("SCS");

unique_ptr<IP> create_ip(string file_name) {
	if (FLAGS_ip_proto_directory.empty()) {
		ifstream fin(file_name);
		return read_MPS(fin);
	} else {
		ifstream fin(file_name, ios::binary);
		fin.seekg(0, std::ios::end);
		auto length = fin.tellg();
		fin.seekg(0, std::ios::beg);
		vector<char> data(length);
		fin.read(data.data(), length);
		return make_unique<IP>(data.data(), length);
	}
}

void process_file(string file_name, string name) {
	cout << "===== " << name << " =====\n";
	ofstream log(file_name + ".log");

	double ground_truth = numeric_limits<double>::quiet_NaN();

	if (FLAGS_ip_proto_directory.empty()) {
		auto index = find(problem_names.begin(), problem_names.end(), name) - problem_names.begin();
		minimum_core_assert(index < problem_names.size(), "Solution not found.");

		auto sol = problem_solutions.at(index);
		ground_truth = sol.first * pow(10.0, sol.second);
	}

	if (FLAGS_verify_with_simplex || !FLAGS_ip_proto_directory.empty()) {
		auto ip = create_ip(file_name);
		// GlopSolver solver(ip.get());
		GlpkSolver solver;
		solver.set_silent(true);
		check(solver.solutions(ip.get())->get(), "LP solver failed.");

		double calulated_objective = ip->get_entire_objective();
		if (FLAGS_ip_proto_directory.empty()) {
			if (!is_equal(calulated_objective, ground_truth, 1e-5)) {
				cout << "LP solver could not solve problem.\n";
				cout << "Ground truth : " << ground_truth << "\n";
				cout << "Calculated   : " << calulated_objective << "\n";
				std::exit(1);
			}
		} else {
			ground_truth = calulated_objective;
		}
	}
	log << "Solution should be " << ground_truth << "\n";

	{
		auto ip = create_ip(file_name);

		stringstream first_order_log;
		FirstOrderOptions options;
		options.maximum_iterations = FLAGS_max_iterations;
		options.tolerance = 1e-7;
		options.rescale_c = true;
		options.log_function = [&](const std::string& str) { first_order_log << str << "\n"; };

		PrimalDualSolver primal_dual_solver;
		primal_dual_solver.options() = options;

		primal_dual_stats.solve_and_update(ground_truth, *ip, &primal_dual_solver, log);

		log << first_order_log.str() << endl;
	}

	{
		auto ip = create_ip(file_name);

		stringstream first_order_log;
		FirstOrderOptions options;
		options.maximum_iterations = FLAGS_max_iterations;
		options.tolerance = 1e-7;
		options.rescale_c = false;
		options.log_function = [&](const std::string& str) { first_order_log << str << "\n"; };

		PrimalDualSolver primal_dual_solver;
		primal_dual_solver.options() = options;

		primal_dual_org_stats.solve_and_update(ground_truth, *ip, &primal_dual_solver, log);

		log << first_order_log.str() << endl;
	}

	{
		auto ip = create_ip(file_name);

		ScsSolver scs_solver;
		scs_solver.set_silent(true);
		scs_solver.set_max_iterations(FLAGS_max_iterations);
		scs_solver.set_convergence_tolerance(1e-5);

		scs_stats.solve_and_update(ground_truth, *ip, &scs_solver, log);
	}
}

int main_program(int argc, char* argv[]) {
	string path;
	if (FLAGS_ip_proto_directory.empty()) {
		path = data::get_directory() + "/netlib";
	} else {
		path = FLAGS_ip_proto_directory;
	}

	for (auto& p : std::experimental::filesystem::directory_iterator(path)) {
		auto name = p.path().stem().string();
		std::transform(name.begin(), name.end(), name.begin(), ::tolower);

		if (FLAGS_ip_proto_directory.empty() && !ends_with(p.path().filename().string(), ".SIF")) {
			continue;
		}
		if (!FLAGS_ip_proto_directory.empty() && !ends_with(p.path().filename().string(), ".ip")) {
			continue;
		}

		try {
			process_file(p.path().string(), name);
		} catch (std::runtime_error& err) {
			cout << err.what() << "\n\n";
			continue;
		}
	}

	primal_dual_stats.print();
	primal_dual_org_stats.print();
	scs_stats.print();
	return 0;
}

int main(int argc, char* argv[]) {
	gflags::ParseCommandLineFlags(&argc, &argv, true);
	try {
		return main_program(argc, argv);
	} catch (std::exception& err) {
		std::cerr << "ERROR: " << err.what() << std::endl;
		return 1;
	}
}
#endif
