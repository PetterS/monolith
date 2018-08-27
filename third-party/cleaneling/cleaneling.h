/***************************************************************************
Copyright (C) 2012 - 2014 Armin Biere JKU Linz

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to
deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.
****************************************************************************

CONTENTS OF THIS FILE
=====================

The generic API for the Cleaneling library.  The actual implementation,
defined in the file 'libcleaneling.cc', is hidden in an internal subclass
'Solver' of the 'AbstractSolver' to avoid exposing implementation details.

There is also an API for accessing global time and memory statistics.

--------------------------------------------------------------------------*/

#ifndef libcleaneling_h_INCLUDED
#define libcleaneling_h_INCLUDED

/*------------------------------------------------------------------------*/

#include <cstdlib>			// need 'size_t'

/*------------------------------------------------------------------------*/

namespace Cleaneling {

// Return values for the actual SAT solving routine 'solve'. These have
// have the same encoding as the exit code of SAT solvers as required
// in the SAT competition.
//
enum Status { UNKNOWN = 0, SATISFIABLE = 10, UNSATISFIABLE = 20 };

// Value of a variable as result of the 'value' function for obtaining
// the assignment of a variable after 'solve' returned 'SATISFIABLE'.
//
enum Value { UNASSIGNED = 0, FALSE = -1, TRUE = 1 };

/*------------------------------------------------------------------------*/

class AbstractSolver {

protected:

  AbstractSolver () { }			// hide actual implementation

public:

  virtual ~AbstractSolver () { }	// C++ requires this

  virtual void addlit (int) = 0;	// add literal (0 terminates clause)

  virtual Status solve () = 0;		// solve problem
  virtual Value value (int) = 0;	// value assigned to literal

  virtual void
    option (const char *, int) = 0;	// set option value

  virtual int
    option (const char *) = 0;		// get current value

  virtual void statistics () = 0;	// print statistics
  virtual size_t maximum_bytes () = 0;  // maximum allocated bytes
  virtual size_t current_bytes () = 0;  // currently allocated bytes
};

/*------------------------------------------------------------------------*/

AbstractSolver * new_solver ();		// factory for actual 'Solver'

/*------------------------------------------------------------------------*/

// Helper functions used by 'cleaneling.cc'.

double seconds ();                      // process time

void print_line (char ch = '-');	// for verbose output only

/*------------------------------------------------------------------------*/
}; // end of 'namespace Cleaneling'
/*------------------------------------------------------------------------*/

#endif
