#pragma once
#include <minimum/core/check.h>
#include <minimum/linear/ip.h>

class OsiSolverInterface;

namespace minimum {
namespace linear {

// Short-hand mostly for the transition to solver classes.
MINIMUM_LINEAR_API bool solve(IP* ip);
MINIMUM_LINEAR_API bool solve_minisat(IP* ip);

// Iterator-like class for obtaining solutions. After the first call to get,
// the solution is available in the IP object (if successful).
class MINIMUM_LINEAR_API Solutions {
   public:
	virtual ~Solutions() {}

	// Call this function repeatedly to obtain the next solution.
	// Not all solvers support enumerating all optimal solutions.
	virtual bool get() = 0;

	// Replaces the problem with a new, similar one, taking advantage of the already
	// solved problem to warm-start the solution.
	virtual void warm_start(IP* ip) { minimum::core::check(false, "Warm-start not implemented"); }

	// Returns the current solver log, if any. Most solvers do
	// not currently implement this.
	virtual std::string log() { return ""; }
};

// This is just a simple “smart pointer” container for a Solutions object.
// The reason this is used instead of a unique_ptr is only because the get() method
// of the unique_ptr is hidden here to avoid confusion.
class SolutionsPointer {
   public:
	SolutionsPointer(std::unique_ptr<Solutions>&& solutions_) : solutions(std::move(solutions_)) {}
	Solutions* operator->() { return solutions.get(); }

   private:
	std::unique_ptr<Solutions> solutions;
};

// API for solving IP instances.
class MINIMUM_LINEAR_API Solver {
   public:
	virtual ~Solver() {}

	// Start solving the given ip. This typicall does some initialization work,
	// whereas the main time is spent in Solutions::get().
	virtual SolutionsPointer solutions(IP* ip) const = 0;

	void set_silent(bool silent_) { silent = silent_; }
	bool silent = false;
};

// Default solver.
class MINIMUM_LINEAR_API IPSolver : public Solver {
   public:
	virtual SolutionsPointer solutions(IP* ip) const override;

	bool solve_relaxation(IP* ip_to_solve) const;

	typedef std::function<void()> CallBack;
	void set_callback(const CallBack& callback_function);
	IPSolver::CallBack callback_function = nullptr;

	// Creates and returns a pointer to a solver.
	//
	// existing_solver can be empty, in which case it will be created.
	// This allows any interface to be created, e.g. OsiGrbSolverInterface,
	// which this library might not know about.
	void get_problem(const IP& ip, std::unique_ptr<OsiSolverInterface>& existing_solver) const;

	// Saves the integer program to an MPS file.
	void save_MPS(const IP& ip, const std::string& file_name) const;

	double time_limit_in_seconds = 0;
};
}  // namespace linear
}  // namespace minimum
