#ifndef _logging_hpp_INCLUDED
#define _logging_hpp_INCLUDED

/*------------------------------------------------------------------------*/
#ifdef LOGGING
/*------------------------------------------------------------------------*/

#include <vector>

namespace CaDiCaL {

// For debugging purposes and to help understanding what the solver is doing
// there is a logging facility which is compiled in by './configure -l'.  It
// still has to be enabled at run-time though (again using the '-l' option
// in the stand-alone solver).  It produces quite a bit of information.

using namespace std;

struct Clause;
struct Internal;

struct Logger {

static void print_log_prefix (Internal *);

// Simple logging of a C-style format string.
//
static void log (Internal *, const char * fmt, ...);

// Prints the format string (with its argument) and then the clause.  The
// clause can also be a zero pointer and then is interpreted as a decision
// (current decision level > 0) or unit clause (zero decision level) and
// printed accordingly.
//
static void log (Internal *, const Clause *, const char *fmt, ...);

// Same as before, except that this is meant for the global 'clause' stack
// used for new clauses (and not for reasons).
//
static void log (Internal *, const vector<int> &, const char *fmt, ...);

// Another variant, to avoid copying (without logging).
//
static void log (Internal *,
                 const vector<int>::const_iterator & begin,
                 const vector<int>::const_iterator & end,
                 const char *fmt, ...);

static void log_empty_line (Internal *);

};

}

/*------------------------------------------------------------------------*/

// Make sure that 'logging' code is really not included (second case of the
// '#ifdef') if logging code is not included.

#define LOG(...) \
do { \
  if (!internal->opts.log) break; \
  Logger::log (internal, __VA_ARGS__); \
} while (0)

/*------------------------------------------------------------------------*/
#else   // end of 'then' part of 'ifdef LOGGING'
/*------------------------------------------------------------------------*/

#define LOG(...) do { } while (0)

/*------------------------------------------------------------------------*/
#endif // end of 'else' part of 'ifdef LOGGING'
/*------------------------------------------------------------------------*/
#endif
