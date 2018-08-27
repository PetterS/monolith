#pragma once

#include <atomic>
#include <exception>
#include <memory>
#include <stdexcept>

namespace minimum {
namespace core {

// Allows capturing and rethrowing exceptions for OpenMP loops.
//
//    OpenMpExceptionStore exception_store;
//    #pragma omp parallel for
//    for (int i = 0; i < 100; ++i) {
//        try {
//            do_something();
//        } catch (...) {
//            exception_store.store();
//        }
//    }
//    exception_store.throw_if_available();
class OpenMpExceptionStore {
   public:
	void store() {
		auto t = std::make_shared<std::exception_ptr>(std::current_exception());
		std::atomic_store(&thrown_exception, t);
	}

	void throw_if_available() {
		has_called_throw_if_available = true;
		if (thrown_exception) {
			std::rethrow_exception(*thrown_exception);
		}
	}

	~OpenMpExceptionStore() noexcept(false) {
		if (!has_called_throw_if_available) {
			throw std::logic_error("OpenMpExceptionStore::throw_if_available must be called.");
		}
	}

   private:
	std::shared_ptr<std::exception_ptr> thrown_exception;
	bool has_called_throw_if_available = false;
};
}  // namespace core
}  // namespace minimum
