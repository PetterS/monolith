#include <catch.hpp>

#include <minimum/core/openmp.h>
using minimum::core::OpenMpExceptionStore;

TEST_CASE("throw") {
	OpenMpExceptionStore exception_store;
#pragma omp parallel for
	for (int i = 0; i < 100; ++i) {
		try {
			throw std::runtime_error("Fail.");
		} catch (...) {
			exception_store.store();
		}
	}
	CHECK_THROWS(exception_store.throw_if_available());
}

TEST_CASE("nothrow") {
	OpenMpExceptionStore exception_store;
#pragma omp parallel for
	for (int i = 0; i < 100; ++i) {
		try {
		} catch (...) {
			exception_store.store();
		}
	}
	CHECK_NOTHROW(exception_store.throw_if_available());
}

TEST_CASE("invalid") { CHECK_THROWS(OpenMpExceptionStore()); }
