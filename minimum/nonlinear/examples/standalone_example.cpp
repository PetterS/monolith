// Petter Strandmark.
//
// When everything is built and installed correctly, it is possible
// to build this file with:
//
//     g++ -std=c++14 standalone_example.cpp -lminimum_nonlinear
//
// (on Linux, that is. Windows uses Visual Studio.)
//

#include <iostream>

#include <minimum/nonlinear/function.h>

int main() {
	minimum::nonlinear::Function function;
	std::cout << "OK.\n";
}
