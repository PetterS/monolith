// Petter Strandmark 2013.
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <random>

#include <badiff.h>
#include <fadiff.h>

#include <catch.hpp>

#include <minimum/curvature/data_term.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace minimum::curvature;

TEST_CASE("Interpolate/Length", "") {
	int M = 5;
	int N = 5;
	int O = 5;
	std::vector<double> un(M * N * O, 1.0);
	std::vector<double> voxeldimensions(3, 1.0);

	PieceWiseConstant data_term(&un[0], M, N, O, voxeldimensions);

	// Makes sure the length regularization is ok
	auto test_length = [&data_term](
	                       double x, double y, double z, double dx, double dy, double dz) -> void {
		double forward_cost = data_term.evaluate_line_integral(x, y, z, x + dx, y + dy, z + dz);
		double backward_cost = data_term.evaluate_line_integral(x + dx, y + dy, z + dz, x, y, z);

		double length = std::sqrt(double(dx * dx + dy * dy + dz * dz));

		CHECK(Approx(forward_cost) == backward_cost);
		CHECK(Approx(forward_cost) == length);
	};

	test_length(0, 0, 0, 4, 4, 4);
	test_length(1, 2, 3, -1, -1, -1);
	test_length(0, 2, 0, 1, -1, 4);
	test_length(3, 2, 1, -1, 1, 1);
	test_length(2, 4, 4, 0, 0, -2);
}

TEST_CASE("Interpolate/Unary", "") {
	int M = 5;
	int N = 5;
	int O = 5;
	int id;
	std::vector<double> un(M * N * O);
	std::vector<double> voxeldimensions(3, 1.0);

	// Manually fill it up
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < O; k++) {
				id = i + j * M + k * M * N;
				un[id] = id;
			}
		}
	}

	PieceWiseConstant data_term(&un[0], M, N, O, voxeldimensions);

	// Makes sure the length regularization is ok
	auto test_unary = [](PieceWiseConstant& data_term,
	                     double x,
	                     double y,
	                     double z,
	                     double dx,
	                     double dy,
	                     double dz,
	                     double value) -> void {
		double forward_cost = data_term.evaluate_line_integral(x, y, z, x + dx, y + dy, z + dz);
		double backward_cost = data_term.evaluate_line_integral(x + dx, y + dy, z + dz, x, y, z);

		CHECK(Approx(forward_cost) == backward_cost);
		CHECK(Approx(forward_cost) == value);
	};

	test_unary(data_term, 2, 2, 0, 0, 2, 0, 34.0);
	test_unary(data_term, 1, 2, 3, 0, 0, 0, 0);
	test_unary(data_term, 0, 0, 0, 4, 4, 4, 4.295486002770816e+02);
	test_unary(data_term, 1, 2, 4, -1, -1, -1, 1.654108574926568e+02);
	test_unary(data_term, 4, 0, 2, -1, 4, 0, 2.618172076554679e+02);
	test_unary(data_term, 4, 4, 2, 0, -4, 0, 256);

	// Infinity tests
	M = 4;
	N = 4;
	O = 1;
	std::vector<double> un_inf(M * N * O);

	id = 0;
	un_inf[id++] = 0.435994902142004;
	un_inf[id++] = 0.025926231827891;
	un_inf[id++] = 0.549662477878709;
	un_inf[id++] = 0.435322392618277;
	un_inf[id++] = 0.420367802087489;
	un_inf[id++] = std::numeric_limits<double>::infinity();
	un_inf[id++] = 0.204648634037843;
	un_inf[id++] = 0.619270966350664;
	un_inf[id++] = 0.299654673674523;
	un_inf[id++] = 0.266827275102867;
	un_inf[id++] = std::numeric_limits<double>::infinity();
	un_inf[id++] = 0.529142094277039;
	un_inf[id++] = 0.134579945344934;
	un_inf[id++] = 0.513578121265746;
	un_inf[id++] = 0.184439865646915;
	un_inf[id] = 0.785335147816673;

	PieceWiseConstant data_inf(&un_inf[0], M, N, O, voxeldimensions);

	test_unary(data_inf, 0, 0, 0, 3, 1, 0, 0.799221196019935);
	test_unary(data_inf, 0, 0, 0, 1, 0, 0, 0.230960566984948);
	test_unary(data_inf, 1, 0, 0, 1, 1, 0, 0.163041051224839);
	test_unary(data_inf, 2, 1, 0, -1, 1, 0, 0.333383812519488);
	test_unary(data_inf, 1, 2, 0, 1, 1, 0, 0.319094055350835);
	test_unary(data_inf, 2, 3, 0, 1, 0, 0, 0.484887506731794);
}

TEST_CASE("Interpolate/Unary inf") {}

template <typename DataTermImplementation>
void perform_stress_test() {
	std::mt19937_64 rng(std::mt19937_64::default_seed);
	auto rand = std::bind(std::uniform_real_distribution<double>(0.0, 1.0), rng);

	for (int iter = 1; iter <= 100; ++iter) {
		int M = 4 + int(6 * rand());
		int N = 4 + int(6 * rand());
		int O = 4 + int(6 * rand());
		std::vector<double> un(M * N * O);
		std::vector<double> voxeldimensions(3, 1.0);

		// Manually fill it up
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {
				for (int k = 0; k < O; k++) {
					auto id = i + j * M + k * M * N;
					un[id] = rand();
				}
			}
		}

		DataTermImplementation data_term(&un[0], M, N, O, voxeldimensions);

		// TODO: Starting in integer points fails.
		double x1 = (1 + rand() * (M - 3));
		double y1 = (1 + rand() * (N - 3));
		double z1 = (1 + rand() * (O - 3));
		double x2 = (1 + rand() * (M - 3));
		double y2 = (1 + rand() * (N - 3));
		double z2 = (1 + rand() * (O - 3));

		int n = 100000;
		double integral = 0;
		integral += data_term.evaluate(x1, y1, z1) / 2.0;
		integral += data_term.evaluate(x2, y2, z2) / 2.0;
		for (int k = 1; k <= n - 1; ++k) {
			double t = double(k) / double(n);
			auto x = (1 - t) * x1 + t * x2;
			auto y = (1 - t) * y1 + t * y2;
			auto z = (1 - t) * z1 + t * z2;
			integral += data_term.evaluate(x, y, z);
		}
		double d = std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
		integral *= d;
		integral /= double(n);

		CAPTURE(M);
		CAPTURE(N);
		CAPTURE(O);
		CAPTURE(x1);
		CAPTURE(y1);
		CAPTURE(z1);
		CAPTURE(x2);
		CAPTURE(y2);
		CAPTURE(z2);
		CAPTURE(integral);
		double real_integral = data_term.evaluate_line_integral(x1, y1, z1, x2, y2, z2);
		CAPTURE(real_integral);
		if (d > 1e-6) {
			CHECK(
			    (std::abs(real_integral - integral) / (0.5 * (integral + real_integral)) <= 1e-4));
		} else {
			CHECK(std::abs(real_integral - integral) <= 1e-4);
		}
	}
}

TEST_CASE("Stress test -- PieceWiseConstant") { perform_stress_test<PieceWiseConstant>(); }

/*
TEST_CASE("Stress test -- TriLinear")
{
    perform_stress_test<TriLinear>();
}
*/

template <typename DataTermImplementation>
void perform_differentiation_test() {
	std::mt19937_64 rng(std::mt19937_64::default_seed);
	auto rand = std::bind(std::uniform_real_distribution<double>(0.0, 1.0), rng);

	for (int iter = 1; iter <= 100; ++iter) {
		int M = 4 + int(6 * rand());
		int N = 4 + int(6 * rand());
		int O = 4 + int(6 * rand());
		std::vector<double> un(M * N * O);
		std::vector<double> voxeldimensions(3, 1.0);

		// Manually fill it up
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {
				for (int k = 0; k < O; k++) {
					auto id = i + j * M + k * M * N;
					un[id] = rand();
				}
			}
		}

		DataTermImplementation data_term(&un[0], M, N, O, voxeldimensions);

		double x1 = (1 + rand() * (M - 3));
		double y1 = (1 + rand() * (N - 3));
		double z1 = (1 + rand() * (O - 3));
		double x2 = (1 + rand() * (M - 3));
		double y2 = (1 + rand() * (N - 3));
		double z2 = (1 + rand() * (O - 3));

		double h = 1e-6;
		double f1, f2;
		f1 = data_term.evaluate_line_integral(x1 - h, y1, z1, x2, y2, z2);
		f2 = data_term.evaluate_line_integral(x1 + h, y1, z1, x2, y2, z2);
		double dIdx1_numeric = (f2 - f1) / (2 * h);
		f1 = data_term.evaluate_line_integral(x1, y1 - h, z1, x2, y2, z2);
		f2 = data_term.evaluate_line_integral(x1, y1 + h, z1, x2, y2, z2);
		double dIdy1_numeric = (f2 - f1) / (2 * h);
		f1 = data_term.evaluate_line_integral(x1, y1, z1 - h, x2, y2, z2);
		f2 = data_term.evaluate_line_integral(x1, y1, z1 + h, x2, y2, z2);
		double dIdz1_numeric = (f2 - f1) / (2 * h);

		f1 = data_term.evaluate_line_integral(x1, y1, z1, x2 - h, y2, z2);
		f2 = data_term.evaluate_line_integral(x1, y1, z1, x2 + h, y2, z2);
		double dIdx2_numeric = (f2 - f1) / (2 * h);
		f1 = data_term.evaluate_line_integral(x1, y1, z1, x2, y2 - h, z2);
		f2 = data_term.evaluate_line_integral(x1, y1, z1, x2, y2 + h, z2);
		double dIdy2_numeric = (f2 - f1) / (2 * h);
		f1 = data_term.evaluate_line_integral(x1, y1, z1, x2, y2, z2 - h);
		f2 = data_term.evaluate_line_integral(x1, y1, z1, x2, y2, z2 + h);
		double dIdz2_numeric = (f2 - f1) / (2 * h);

		typedef fadbad::F<double, 6> F;
		F f_x1 = x1;
		f_x1.diff(0);
		F f_y1 = y1;
		f_y1.diff(1);
		F f_z1 = z1;
		f_z1.diff(2);
		F f_x2 = x2;
		f_x2.diff(3);
		F f_y2 = y2;
		f_y2.diff(4);
		F f_z2 = z2;
		f_z2.diff(5);

		F result = data_term.evaluate_line_integral(f_x1, f_y1, f_z1, f_x2, f_y2, f_z2);
		double dIdx1_auto = result.d(0);
		double dIdy1_auto = result.d(1);
		double dIdz1_auto = result.d(2);
		double dIdx2_auto = result.d(3);
		double dIdy2_auto = result.d(4);
		double dIdz2_auto = result.d(5);

		CAPTURE(M);
		CAPTURE(N);
		CAPTURE(O);
		CAPTURE(x1);
		CAPTURE(y1);
		CAPTURE(z1);
		CAPTURE(x2);
		CAPTURE(y2);
		CAPTURE(z2);
		CAPTURE(dIdx1_auto);
		CAPTURE(dIdx1_numeric);
		CAPTURE(dIdy1_auto);
		CAPTURE(dIdy1_numeric);
		CAPTURE(dIdz1_auto);
		CAPTURE(dIdz1_numeric);

		CAPTURE(dIdx2_auto);
		CAPTURE(dIdx2_numeric);
		CAPTURE(dIdy2_auto);
		CAPTURE(dIdy2_numeric);
		CAPTURE(dIdz2_auto);
		CAPTURE(dIdz2_numeric);

		CHECK(
		    (std::abs(dIdx1_numeric - dIdx1_auto) / (0.5 * (dIdx1_numeric + dIdx1_auto)) <= 1e-4));
		CHECK(
		    (std::abs(dIdy1_numeric - dIdy1_auto) / (0.5 * (dIdy1_numeric + dIdy1_auto)) <= 1e-4));
		CHECK(
		    (std::abs(dIdz1_numeric - dIdz1_auto) / (0.5 * (dIdz1_numeric + dIdz1_auto)) <= 1e-4));
		CHECK(
		    (std::abs(dIdx2_numeric - dIdx2_auto) / (0.5 * (dIdx2_numeric + dIdx2_auto)) <= 1e-4));
		CHECK(
		    (std::abs(dIdy2_numeric - dIdy2_auto) / (0.5 * (dIdy2_numeric + dIdy2_auto)) <= 1e-4));
		CHECK(
		    (std::abs(dIdz2_numeric - dIdz2_auto) / (0.5 * (dIdz2_numeric + dIdz2_auto)) <= 1e-4));
	}
}

TEST_CASE("Differentiation test -- PieceWiseConstant") {
	perform_differentiation_test<PieceWiseConstant>();
}

/*
TEST_CASE("Differentiation test -- TriLinear")
{
  perform_differentiation_test<TriLinear>();
}
*/

#ifdef USE_OPENMP
TEST_CASE("Interpolate/Benchmark", "") {
	int M = 20;
	int N = 20;
	int O = 20;
	std::vector<double> un(M * N * O, 1.0);
	std::vector<double> voxeldimensions(3, 1.0);

	PieceWiseConstant data_term(&un[0], M, N, O, voxeldimensions);

	double start_time = omp_get_wtime();
	for (int i = 0; i < 100000; ++i) {
		data_term.evaluate_line_integral(1.0, 5.0, 3.0, 4.0, 2.0, 2.0);
		data_term.evaluate_line_integral(1.0, 2.0, 2.0, 4.0, 2.0, 3.0);
		data_term.evaluate_line_integral(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		data_term.evaluate_line_integral(1.0, 1.0, 0.0, 1.0, 0.0, 1.0);
	}
	double end_time = omp_get_wtime();
	double elapsed_time = end_time - start_time;
	std::cerr << "Elapsed time: " << elapsed_time << std::endl;
}
#endif
