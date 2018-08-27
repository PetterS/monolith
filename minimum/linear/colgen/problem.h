#pragma once

#include <minimum/linear/colgen/column.h>
#include <minimum/linear/colgen/column_pool.h>
#include <minimum/linear/colgen/export.h>
#include <minimum/linear/colgen/proto.h>

namespace minimum {
namespace linear {
class IP;
namespace colgen {

class MINIMUM_LINEAR_COLGEN_API Problem {
   public:
	Problem(int number_of_rows);
	Problem(const Problem&) = delete;
	~Problem();

	// Called when new columns are needed (every iteration).
	virtual void generate(const std::vector<double>& dual_variables) = 0;

	struct FixInformation {
		int iteration = 0;
		double objective_change = 0;
	};

	// Override and return number of columns fixed. If a negative value is
	// returned, the colgen iteration will stop with a success.
	virtual int fix(const FixInformation& information) = 0;

	// Optionally override to provide a rounded current integer solution
	// value to be displayed in the log.
	// This method could also, e.g., save the best solution to disk.
	//
	// The objective constant set by set_objective_constant() will be
	// added afterwards.
	virtual double integral_solution_value();

	double solve();

	const std::vector<std::size_t>& active_columns() const;
	const Column& column_at(std::size_t i) const { return pool.at(i); }

	// Constant added to the objective function (for display purposes).
	void set_objective_constant(double constant);

   protected:
	ColumnPool pool;

	// Whether logging of additional data is enabled.
	const bool proto_logging_enabled;
	// Optionally override this method to return data to be logged each iteration.
	// Will only call this method if proto_logging_enabled.
	virtual proto::LogEntry create_log_entry() const;
	// Attaches this event to the next proto log entry written.
	void attach_event(proto::Event event);

	// This method is called whenever Ctrl-C is detected during optimization.
	// If this method returns true, optimization continues, otherwise optimization
	// stops. The default behaviour is to throw an std::runtime_error.
	//
	// This function is not called in the actual interrupt handler, so it is
	// safe to use.
	virtual bool interrupt_handler();

	// Returns an IP representing the main linear program for the columns
	// in active_columns().
	std::unique_ptr<IP> create_ip(bool use_integer_variables = false) const;

	void set_row_lower_bound(int row, double lower);
	void set_row_upper_bound(int row, double upper);

   private:
	class Implementation;
	Implementation* impl;
};
}  // namespace colgen
}  // namespace linear
}  // namespace minimum