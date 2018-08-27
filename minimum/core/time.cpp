
#include <chrono>
#include <ctime>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <minimum/core/time.h>

namespace minimum {
namespace core {

double wall_time() {
#ifdef USE_OPENMP
	return omp_get_wtime();
#else
	using Clock = std::chrono::high_resolution_clock;
	return double(Clock::now().time_since_epoch().count()) * double(Clock::period::num)
	       / double(Clock::period::den);
#endif
}

double cpu_time() { return double(std::clock()) / double(CLOCKS_PER_SEC); }
}  // namespace core
}  // namespace minimum
