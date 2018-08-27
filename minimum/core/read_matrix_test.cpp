#include <catch.hpp>

#include <minimum/core/read_matrix.h>
#include <iostream>

using namespace std;
using minimum::core::read_matrix;
using minimum::core::to_string;

TEST_CASE("integer_2_by_2") {
	istringstream sin("  +123,,  -456;, 789\n\r\r\n111; 222; 333");
	CHECK((read_matrix<int>(sin) == vector<vector<int>>{{123, -456, 789}, {111, 222, 333}}));
}

TEST_CASE("integer_2_by_2_positive") {
	istringstream sin(" 123,,  456;, 789\n\r\r\n111; 222; 333\n\n\n");
	CHECK((read_matrix<int>(sin) == vector<vector<int>>{{123, 456, 789}, {111, 222, 333}}));
}

TEST_CASE("empty") {
	istringstream sin("");
	CHECK((read_matrix<int>(sin) == vector<vector<int>>{}));
}

TEST_CASE("three") {
	istringstream sin("\t123,,  456;, 789\n111; 222; 333\n\n11, 22\n\n33 44 55\n66 77 88");
	CHECK((read_matrix<int>(sin, true) == vector<vector<int>>{{123, 456, 789}, {111, 222, 333}}));
	CHECK((read_matrix<int>(sin, true) == vector<vector<int>>{{11, 22}}));
	CHECK((read_matrix<int>(sin, true) == vector<vector<int>>{{33, 44, 55}, {66, 77, 88}}));
}

TEST_CASE("floating") {
	istringstream sin("  +123.4,,  -456.5;, 789.6\n\r\r\n111; 222; 333");
	CHECK((read_matrix<double>(sin)
	       == vector<vector<double>>{{123.4, -456.5, 789.6}, {111, 222, 333}}));
}

TEST_CASE("floating_positive") {
	istringstream sin("  123.4,,  456.5;, 789.6\n\r\r\n111; 222; 333\n\n");
	CHECK((read_matrix<double>(sin)
	       == vector<vector<double>>{{123.4, 456.5, 789.6}, {111, 222, 333}}));
}

TEST_CASE("floating_empty") {
	istringstream sin("");
	CHECK((read_matrix<double>(sin) == vector<vector<double>>{}));
}

TEST_CASE("floating_three") {
	istringstream sin(" 123,,  456;, 789\n111; 222; 333\n\n11, 22\n\n33 44 55\n66 77 88");
	CHECK((read_matrix<double>(sin, true)
	       == vector<vector<double>>{{123, 456, 789}, {111, 222, 333}}));
	CHECK((read_matrix<double>(sin, true) == vector<vector<double>>{{11, 22}}));
	CHECK((read_matrix<double>(sin, true) == vector<vector<double>>{{33, 44, 55}, {66, 77, 88}}));
}

TEST_CASE("nan") {
	istringstream sin("nan nan");
	CHECK(read_matrix<double>(sin).at(0).size() == 2);
}

TEST_CASE("floating_32") {
	istringstream sin("  +123.4,,  -456.5;, 789.6\n\r\r\n111; 222; 333");
	CHECK((read_matrix<float>(sin)
	       == vector<vector<float>>{{123.4f, -456.5f, 789.6f}, {111.f, 222.f, 333.f}}));
}

TEST_CASE("comment") {
	string data =
	    R"(# The matrix.
	   1, 2, 3, // First row
	   4, 5, 6  // Second row.
	   // Pause
	   7, 8, 9 // Last one.
	   
	   )";

	{
		istringstream sin(data);
		CHECK((read_matrix<int>(sin) == vector<vector<int>>{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}));
	}
	{
		istringstream sin(data);
		CHECK(
		    (read_matrix<double>(sin) == vector<vector<double>>{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}));
	}
}
