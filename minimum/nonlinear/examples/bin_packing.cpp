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
#include <vector>
using namespace std;

#include <absl/strings/substitute.h>
#include <gflags/gflags.h>

#include <minimum/core/check.h>
#include <minimum/core/color.h>
#include <minimum/core/main.h>
#include <minimum/core/random.h>
#include <minimum/core/range.h>
#include <minimum/nonlinear/auto_diff_term.h>
#include <minimum/nonlinear/constrained_function.h>
#include <minimum/nonlinear/solver.h>
using namespace minimum::core;
using namespace minimum::nonlinear;

DEFINE_int32(n, 50, "Number of circles.");
DEFINE_int32(seed, 0, "Seed for the random initial position.");
DEFINE_bool(verbose, false, "Print solver logs.");

void write_circles_to_svg(std::ostream* out,
                          double box_width,
                          const vector<double>& x,
                          const vector<double>& y,
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
	for (int i : range(x.size())) {
		*out
		    << absl::Substitute(
		           R"(<circle cx="$0" cy="$1" r="$2" stroke="black" stroke-width="$3" fill="red" />)",
		           scale * x.at(i),
		           scale * y.at(i),
		           scale * r.at(i),
		           scale * 0.04)
		    << endl;
	}
	*out << "</svg>" << endl;
}

int main_program(int num_args, char* args[]) {
	check(FLAGS_n >= 2, "Need more than two circles.");
	auto rng_initial = repeatably_seeded_engine<mt19937_64>(FLAGS_seed, 1, 9, 3, 3, 9);
	auto rng_radius = repeatably_seeded_engine<mt19937_64>(1, 9, 3, 3, 9);

	Timer timer("Creating data");
	vector<double> x(FLAGS_n);
	vector<double> y(FLAGS_n);
	vector<double> r(FLAGS_n);
	uniform_real_distribution<double> r_dist(0.1, 2);
	uniform_real_distribution<double> pos_dist(0, sqrt((double)FLAGS_n));
	double area = 0;
	for (int i : range(FLAGS_n)) {
		x[i] = pos_dist(rng_initial);
		y[i] = pos_dist(rng_initial);
		r[i] = r_dist(rng_radius);
		area += 3.141592 * r[i] * r[i];
	}
	double S = sqrt(area);
	cerr << "Starting S: " << S << endl;

	timer.next("Setting up function and constraints.");
	ConstrainedFunction f;
	f.add_term(make_differentiable<1>([](auto S) { return S[0]; }), &S);

	for (int i : range(FLAGS_n)) {
		double Ri = r[i];
		for (int j : range(i + 1, FLAGS_n)) {
			double Rj = r[j];
			f.add_constraint_term(
			    to_string("radius", i, ",", j),
			    make_differentiable<1, 1, 1, 1>([Ri, Rj](auto xi, auto yi, auto xj, auto yj) {
				    auto dx = xi[0] - xj[0];
				    auto dy = yj[0] - yi[0];
				    auto Rsum = Ri + Rj;
				    return Rsum * Rsum - dx * dx - dy * dy;
			    }),
			    &x[i],
			    &y[i],
			    &x[j],
			    &y[j]);
		}

		auto upper = make_differentiable<1, 1>([Ri](auto xi, auto S) { return xi[0] - S[0] + Ri; });
		auto lower = make_differentiable<1>([Ri](auto xi) { return Ri - xi[0]; });
		f.add_constraint_term(to_string("xupper", i), upper, &x[i], &S);
		f.add_constraint_term(to_string("xlower", i), lower, &x[i]);
		f.add_constraint_term(to_string("yupper", i), upper, &y[i], &S);
		f.add_constraint_term(to_string("ylower", i), lower, &y[i]);
	}
	timer.OK();

	SolverResults results;
	LBFGSSolver solver;
	if (!FLAGS_verbose) {
		solver.log_function = nullptr;
	}
	f.solve(solver, &results);

	cerr << results << endl;
	cerr << "S=" << S << endl;
	write_circles_to_svg(&cout, S, x, y, r);
	return 0;
}

int main(int argc, char* argv[]) { return main_runner(main_program, argc, argv); }
