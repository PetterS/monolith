// Bin Packing
//
// See e.g.
// http://yetanothermathprogrammingconsultant.blogspot.com/2018/05/multi-start-non-linear-programming-poor.html
// for the problem formulation.
//
// Outputs an SVG to stdout and logs to stderr.
//
#include <iostream>
#include <random>
#include <string_view>
#include <vector>
using namespace std;

#include <absl/strings/substitute.h>
#include <gflags/gflags.h>

#include <minimum/constrained/successive.h>
#include <minimum/core/check.h>
#include <minimum/core/color.h>
#include <minimum/core/main.h>
#include <minimum/core/numeric.h>
#include <minimum/core/random.h>
#include <minimum/core/range.h>
#include <minimum/core/time.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/solver.h>
#include <minimum/nonlinear/auto_diff_term.h>
#include <minimum/nonlinear/constrained_function.h>
#include <minimum/nonlinear/solver.h>
using namespace minimum::core;
using namespace minimum::linear;
using namespace minimum::nonlinear;

DEFINE_int32(n, 50, "Number of circles.");
DEFINE_int32(repetitions,
             1,
             "Number of different initial positions to try (with the same circles).");
DEFINE_int32(seed, 0, "Seed for the random initial position.");
DEFINE_bool(verbose, false, "Print solver logs.");
DEFINE_string(solver,
              "successive",
              "Solver type. Equal to \"successive\" (default), or \"agumented\".");

void write_circles_to_svg(std::ostream* out,
                          double box_width,
                          const vector<vector<double>> p,
                          const vector<double>& r) {
	const double scale = 20;
	*out << R"(<?xml version="1.0" encoding="UTF-8" standalone="no"?>)" << endl;
	*out << absl::Substitute(R"(<svg
		xmlns:dc="http://purl.org/dc/elements/1.1/"
		xmlns:cc="http://creativecommons.org/ns#"
		xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
		xmlns:svg="http://www.w3.org/2000/svg"
		xmlns="http://www.w3.org/2000/svg"
		height="$0"
		width="$0"
		version="1.1">)",
	                         scale * box_width)
	     << endl;
	*out << absl::Substitute(
	    R"(<rect x="0" y="0" width="$0" height="$0" stroke-width="$1" style="fill:gray;stroke:black;" />)",
	    scale * box_width,
	    scale * 0.04);
	for (int i : range(p.size())) {
		*out
		    << absl::Substitute(
		           R"(<circle cx="$0" cy="$1" r="$2" stroke="black" stroke-width="$3" fill="red" />)",
		           scale * p[i].at(0),
		           scale * p[i].at(1),
		           scale * r.at(i),
		           scale * 0.04)
		    << endl;
	}
	*out << "</svg>" << endl;
}

double log_info(string_view iteration, ConstrainedFunction& constrained_function, double time) {
	double objective = constrained_function.objective().evaluate();
	double max_violation = 0;
	for (auto& itr : constrained_function.constraints()) {
		itr.second.evaluate_and_cache();
		max_violation = max(max_violation, itr.second.violation());
	}
	cerr << setw(6) << right << iteration << " " << setw(12) << right << objective << " "
	     << setw(12) << right << max_violation << " " << setw(12) << right << time << "\n";
	return objective;
}

double solve_problem(int initial_position_seed) {
	check(FLAGS_n >= 2, "Need more than two circles.");
	auto rng_initial = repeatably_seeded_engine<mt19937_64>(initial_position_seed, 1, 9, 3, 3, 9);
	auto rng_radius = repeatably_seeded_engine<mt19937_64>(1, 9, 3, 3, 9);

	Timer timer("Creating data");
	vector<vector<double>> p(FLAGS_n, vector<double>{0, 0});
	vector<double> r(FLAGS_n);
	uniform_real_distribution<double> r_dist(0.1, 2);
	uniform_real_distribution<double> pos_dist(0, sqrt((double)FLAGS_n));
	double area = 0;
	for (int i : range(FLAGS_n)) {
		p[i][0] = pos_dist(rng_initial);
		p[i][1] = pos_dist(rng_initial);
		r[i] = r_dist(rng_radius);
		area += 3.141592 * r[i] * r[i];
	}
	double S = sqrt(area);
	timer.OK();
	cerr << "Starting S: " << S << endl;

	timer.next("Setting up function and constraints.");
	ConstrainedFunction f;
	f.feasibility_tolerance = 1e-8;
	f.add_term(make_differentiable<1>([](auto S) { return S[0]; }), &S);

	for (int i : range(FLAGS_n)) {
		double Ri = r[i];
		for (int j : range(i + 1, FLAGS_n)) {
			double Rj = r[j];
			f.add_constraint_term(to_string("radius", i, ",", j),
			                      make_differentiable<2, 2>([Ri, Rj](auto pi, auto pj) {
				                      auto dx = pi[0] - pj[0];
				                      auto dy = pj[1] - pi[1];
				                      auto Rsum = Ri + Rj;
				                      return Rsum * Rsum - dx * dx - dy * dy;
			                      }),
			                      p[i].data(),
			                      p[j].data());
		}

		auto xupper =
		    make_differentiable<2, 1>([Ri](auto pi, auto S) { return pi[0] - S[0] + Ri; });
		auto yupper =
		    make_differentiable<2, 1>([Ri](auto pi, auto S) { return pi[1] - S[0] + Ri; });
		auto xlower = make_differentiable<2>([Ri](auto pi) { return Ri - pi[0]; });
		auto ylower = make_differentiable<2>([Ri](auto pi) { return Ri - pi[1]; });
		f.add_constraint_term(to_string("xupper", i), xupper, p[i].data(), &S);
		f.add_constraint_term(to_string("xlower", i), xlower, p[i].data());
		f.add_constraint_term(to_string("yupper", i), yupper, p[i].data(), &S);
		f.add_constraint_term(to_string("ylower", i), ylower, p[i].data());
	}
	timer.OK();

	if (FLAGS_solver == "augmented") {
		SolverResults results;
		LBFGSSolver solver;
		if (!FLAGS_verbose) {
			solver.log_function = nullptr;
		}
		f.solve(solver, &results);
		cerr << results << endl;
	} else if (FLAGS_solver == "successive") {
		minimum::constrained::SuccessiveLinearProgrammingSolver solver;
		solver.log_function = [](string_view s) { cerr << s << endl; };
		solver.initial_step_size = 1'000;
		solver.solve(&f);
	}

	if (FLAGS_repetitions == 1) {
		write_circles_to_svg(&cout, S, p, r);
	}
	cerr << "S=" << S << endl;
	return S;
}

int main_program(int num_args, char* args[]) {
	vector<double> Ss;
	for (int repetition = 0; repetition < FLAGS_repetitions; ++repetition) {
		Ss.emplace_back(solve_problem(FLAGS_seed + repetition));
	}
	cerr << "S = " << to_string(Ss) << endl;
	return 0;
}

int main(int argc, char* argv[]) { return main_runner(main_program, argc, argv); }
