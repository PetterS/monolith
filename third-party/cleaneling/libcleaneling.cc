/**************************************************************************\
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

This file contains the library code for Cleaneling. More precisely it
implements the concrete 'Solver' subclass of 'AbstractSolver'.  The various
parts of this file start with a 'PART' header (search for 'PART').

  PART1   includes, macros, global utilities
  PART2   template code for 'vector' and 'priority_queue'
  PART3   basic internal data structures
  PART4   declaration of the 'Solver' class
  PART5   implemementation of internal member functions of 'Solver'
  PART6   implementation of public / external member functions

In 'PART1' we have the include directives for header files, define macros
and some global utility functions particularly including those for
monitoring memory usage.  The next 'PART2' consists of template code for
a 'vector' and a 'priority_queue' container.

Then the internal data structures of 'Solver' are defined in 'PART3'.  After
the declaration of the 'Solver' class itself in 'PART4' the internal member
functions are implemented in 'PART5'.  Finally in 'PART6' we give the
implementation of the public member functions over-writing their virtual
counter-parts in 'AbstractSolver'.

--------------------------------------------------------------------------*
***** PART1: includes, macros, memory monitoring, utilities ***************
--------------------------------------------------------------------------*/

#include "cleaneling.h"

/*------------------------------------------------------------------------*/

#include <cassert>
#include <climits>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <iomanip>
#include <iostream>

/*------------------------------------------------------------------------*/

// Set operating system dependent configuration options here.

#ifdef __linux__
#define HAVE_GETRUSAGE
#else
// If you happen to have 'getrusage' but your are not on Linux, you can
// still just define 'HAVE_GETRUSAGE' here to support timing satistics.
#endif

/*------------------------------------------------------------------------*/

extern "C" {
#ifdef HAVE_GETRUSAGE
#include <sys/time.h>		// Use 'getrusage' to gather timing stats.
#include <sys/resource.h>	// This is also needed for 'getrusage'.
#endif
#if defined(__linux__)
#include <fpu_control.h>	// Set FPU to double precision on Linux.
#endif
#ifndef NLGL
#include "lglib.h"		// For 'clause_entailed' check to work
  				// we use a Lingeling instance.
#else
struct LGL;			// Declare 'LGL' struct at least.
#endif
};

/*------------------------------------------------------------------------*/
// To enable logging code at compile time use the '-l' option, e.g.
// 'configure.sh -l'.  This is norequired when debugging is included, since
// the '-g' option automatically includes logging code.
//
// If logging code has been compiled in, it still has to be enabled at run
// time using the 'log' option through the 'option' set call, or from the
// command line of the binary using the '-l' command line option
// ('cleaneling -l').  The 'log' option is not legal of course if logging
// code is not compiled in.
//
// This macro can only be used within 'Solver' member functions or where at
// least a references to the 'opts' options and the 'level' of the solver is
// available.
//
#ifndef NLOG
#define LOG(CODE) \
  do { \
    if (!opts.log) break; \
    do { cout << "c @" << level << " " << CODE; } while (0); \
    cout << endl; \
  } while (0)
#define LOGCLAUSE(CLAUSE,CODE) \
  do { \
    if (!opts.log) break; \
    do { \
      cout << "c @" << level << " " << CODE; \
      cout << ' ' << ((CLAUSE)->redundant ? "redundant" : "irredundant"); \
      cout << " clause"; \
      for (const int * p = (CLAUSE)->lits; *p; p++) \
	cout << ' ' << *p; \
      cout << endl; \
    } while (0); \
  } while (0)
#else
#define LOG(...) do { } while (0)
#define LOGCLAUSE(...) do { } while (0)
#endif

/*------------------------------------------------------------------------*/

using namespace std;
using namespace Cleaneling;

/*------------------------------------------------------------------------*/
namespace Cleaneling {
/*------------------------------------------------------------------------*/

static void internal_error (const char * msg) {
  cout << "internal error in 'libcleaneling.a': " << msg << endl << flush;
  abort ();
}

/*------------------------------------------------------------------------*/

// We want to monitor memory usage for each instance of the solver.  So we
// use our own memory manager for allocating and deallocating heap objects.
// It also allows us to use 'realloc' for resizing vectors / stacks and
// checking cheaply for memory leaks.
//
// In principle, the memory manager has the same API as the corresponding C
// library functions 'malloc', 'realloc', and 'free'. These are actually
// used internally anyhow. The difference is to require from the user to
// specify the size of the memory freed in 'free' respectively the current
// size of the allocated block in 'realloc'.
//
// If this argument for 'free' and 'realloc' is incorret this will most
// likely be catched by an assertion, either in 'Mem.dec' or '~Solver'.

struct Mem {

  size_t maximum;		// maximum memory allocated in Bytes
  size_t current;		// currently allocated memory in Bytes

  Mem () : maximum (0), current (0) { }

  // Increment statistics of currently allocated memory.
  //
  void inc (size_t bytes) {
    current += bytes;
    assert (current >= bytes);			// Overflow check.
    if (current > maximum) maximum = current;
  }

  // Decrement statistics of currently allocated memory.
  //
  void dec (size_t bytes) {
    // If this assertion is violated it is most likely due to a wrong
    // 'bytes' argument in 'free' resp. wrong 'old_bytes' argument
    // in 'realloc', but this might have happened way earlier before
    // this assertion is violated.
    //
    assert (current >= bytes);
    current -= bytes;
  }

  // Same semantics as 'C' version.  Needs one additional argument,
  // though, which is described above.
  //
  void * realloc (void * ptr, size_t old_bytes, size_t new_bytes) {
    dec (old_bytes);
    void * res = ::realloc (ptr, new_bytes);
    if (!res) internal_error ("out of memory in 'realloc'");
    inc (new_bytes);
    return res;
  }

  // Same semantics as 'C' version.
  //
  void * malloc (size_t bytes) {
    void * res = ::malloc (bytes);
    if (!res) internal_error ("out of memory in 'malloc'");
    inc (bytes);
    return res;
  }

  // Same semantics as 'C' version.  Needs one additional argument,
  // though, which is described above.
  //
  void free (void * ptr, size_t bytes) { dec (bytes); ::free (ptr); }
};

/*------------------------------------------------------------------------*/

// Process time = user + system time for this process.  This is kind of
// POSIX but not available on Windows.
//
// TODO implement alternatives, e.g. for Windows.
//
double seconds () {
  double res = 0;
#ifdef HAVE_GETRUSAGE
  struct rusage u;
  if (getrusage (RUSAGE_SELF, &u)) return 0;
  res += u.ru_utime.tv_sec + 1e-6 * u.ru_utime.tv_usec;
  res += u.ru_stime.tv_sec + 1e-6 * u.ru_stime.tv_usec;
#endif
  return res;
}

/*------------------------------------------------------------------------*/

// If you do not tell the processor to use double precision all the time,
// you risk diverging behaviour of the solver for different levels of
// optimizations.  We use the same trick as in MiniSAT to force the FPU to
// only use double precision.
//
// TODO This is still not really portable and software floats, as in
// 'PicoSAT' or 'Lingeling' should be preferred instead of 'double'.

static void set_fpu_to_double_precision () {
#ifdef __linux__
  fpu_control_t control;
  _FPU_GETCW (control);
  control &= ~_FPU_EXTENDED;
  control &= ~_FPU_SINGLE;
  control |= _FPU_DOUBLE;
  _FPU_SETCW (control);
#endif
}

static double average (double a, double b) { return b ? a/b : 0; }

static double percent (double a, double b) { return b ? 100*a/b : 0; }

static int sign (int lit) { return lit < 0 ? -1 : 1; }

static int64_t luby (int64_t i) {
  int64_t res = 0, s = 0, k;

  for (k = 1; !res && k < 64; s++, k++) {
    if (i == (1l << k) - 1)
      res = 1l << (k - 1);
  }

  for (k = 1; !res; k++, s++)
    if ((1l << (k - 1)) <= i && i < (1l << k) - 1)
      res = luby (i - (1l << (k-1)) + 1);

  assert (res > 0);

  return res;
}

void print_line (char ch) {
  cout << "c ";
  for (int i = 0; i < 76; i++) cout << ch;
  cout << endl;
}

/*------------------------------------------------------------------------*/

// Provides compact 6 characters formatting for potentially large non
// negative integers using scientific notation.  We avoid exponents as much
// as possible.  This is current only used in 'report'.

class FmtNat {
    int64_t num;
  public:
    FmtNat (int64_t n) : num (n) { assert (n >= 0); }
    friend ostream & operator << (ostream &, const FmtNat &);
};

ostream & operator << (ostream & os, const FmtNat & fn) {
  int64_t n = fn.num;
  assert (n >= 0);						// 123456
  if (n < 10)                    os << "     " << n;		// .....n
  else if (n < 100)              os << "    " << n;		// ....nn
  else if (n < 1000)             os << "   " << n;		// ...nnn
  else if (n < 10000)            os << "  " << n;		// ..nnnn
  else if (n < 100000)           os << " " << n;		// .nnnnn
  else if (n < 1000000)          os << n;			// nnnnnn
  else if (n < 10000000)         os << n/1000          << "e3";	// nnnne3
  else if (n < 100000000)        os << n/10000         << "e4";	// nnnne4
  else if (n < 1000000000)       os << n/100000        << "e5";	// nnnne5
  else if (n < 10000000000l)     os << n/1000000       << "e6";	// nnnne6
  else if (n < 100000000000l)    os << n/10000000      << "e7";	// nnnne7
  else if (n < 1000000000000l)   os << n/100000000     << "e8";	// nnnne8
  else if (n < 10000000000000l)  os << n/1000000000    << "e9";	// nnnne9
  //            dcba987654321             987654321
  //		  ba987654321             987654321
  else if (n < 100000000000000l) os << n/100000000000l << "e11";// nnne11
  //            edcba987654321            ba987654321
  //		   ba987654321            ba987654321
  else {
    int e = 12;
    n /= 1000000000000l;
    while (n >= 1000) n /= 10, e++;
    os << n << 'e' << e;
  }
  return os;
}

/*------------------------------------------------------------------------*/
/**** PART2: template code for 'vector' and 'priority_queue' containers ***/ 
/*------------------------------------------------------------------------*/

// A generic stack implementation very similar to 'std::vector'. We try
// to provide a subset of the interface of the STL standard 'std::vector'
// class.
//
// The main reason for this simplistic replacement of 'std::vector' is to be
// able to monitor memory usage.  We also can enforce calls to 'realloc',
// which sometimes avoids copying memory.  Finally this is way easier to use
// during debugging.
//
// Note that this vector class does not use copy operators for the element
// types during resizing the vector.  So if you rely on copy constructors
// to actually copy elements, this will not work.
//
// The same applies for the destructor.  You have to call the 'release'
// destructor explicity to avoid memory leaks!  The alternative would be to
// add a reference to the memory manager to each instance, which is wasting
// a pointer and in addition needs to be initialized too.

template<class T> class vector {

    T * start;			// Points to first element.
    T * top;			// Points one element past the top element.
    T * allocated;		// Allocated memory for exponential realloc.

    void realloc (Mem & mem, size_t new_size) {
      size_t count = size ();
      size_t old_size = allocated - start;
      size_t old_bytes = old_size * sizeof (T);
      if (new_size) {
	size_t new_bytes = new_size * sizeof (T);
	start = (T*) mem.realloc (start, old_bytes, new_bytes);
	top = start + min (count, new_size);
	allocated = start + new_size;
      } else {
	mem.free (start, old_bytes);
	start = top = allocated = 0;
      }
    }

  public:

    // Most of the functions below are also in the API of 'std::vector'.
    // We also try to use the same semantics.  However, since we handle
    // memory with an explicit memory manager all those function allocating
    // or deallocating heap memory have an additional 'Mem' argument.

    vector () : start (0), top (0), allocated (0) { }

    // We need a reference to 'Mem' to actually delete the memory allocated
    // in this vector.  Since we do not want to waste the space to hold a
    // pointer/reference in each 'vector' seperately, we rely on the user of
    // 'vector' to call 'release' before the destructor is called.

    ~vector () { assert (!start); }  // Fails if above contract is violated.

    bool empty () const { return start == top; }

    size_t size () const { return top - start; }

    T & back () { assert (start < top); return top[-1]; }

    const T & back () const { assert (start < top); return top[-1]; }

    void push_back (Mem & mem, const T & data) { 
      if (top == allocated) realloc (mem, size () ? 2 * size () : 1);
      *top++ = data; 
    }

    T & operator[] (size_t i) { assert (i < size ()); return start[i]; }

    const T & operator[] (const size_t i) const { 
      assert (i < size());
      return start[i];
     }

    void pop_back () { assert (start < top); --top; }

    // This shrinks the allocated space to fit the argument size and is not
    // in the API of the vector class of STL.  It also requires that you
    // kind of call destructors of the elements not used anymore manually
    // before calling this function.
    //
    void shrink (Mem & mem, size_t new_size) { 
      assert (new_size <= size ());
      realloc (mem, new_size);
    }

    // Shrink respectively fit allocated memory to the actual size.
    //
    void shrink (Mem & mem) { shrink (mem, size ()); }

    void release (Mem & mem) {
      size_t bytes = (allocated - start) * sizeof *start;
      mem.free (start, bytes);	// do not call destructors for elements!
      start = top = allocated = 0;
    }

    // In contrast to 'shrink' this function exists in the STL and has the
    // same meaning.  However it also does not call destructors of removed
    // elements.  You need to do this manually before calling this function.
    //
    void resize (size_t new_size) {
      assert (new_size <= size ());
      top = start + new_size;
    }

    // Same as 'resize' but resize the vector to the new top pointer.  This
    // is mainly used for resizing an iterated vector.

    void resize (const T * new_top) {
      assert (start <= new_top && new_top <= allocated);
      size_t new_size = new_top - start;
      top = start + new_size;
    }

    // Same comments as for 'resize' apply here.  Be carefull when calling
    // this function on vectors which contain elements with destructors.
    //
    void clear () { top = start; }

    typedef int (*Cmp)(T*, T*);

    typedef int (*QsortCmp)(const void *, const void *);

    void sort (Cmp cmp) {
      QsortCmp qsort_cmp = (QsortCmp) cmp;
      qsort (start, size (), sizeof (T), qsort_cmp);
    }

    typedef T * iterator;

    iterator begin () { return start; }
    iterator end () { return top; }

    typedef const T * const_iterator;

    const_iterator begin () const { return start; }
    const_iterator end () const { return top; }
};

/*------------------------------------------------------------------------*/

// A specialized heap implementation, which stores priorities and allows
// fast updates of priorities.  Since this does not exist in the STL, we are
// using our own version.  As with the vector class defined above we also
// want to monitor memory usage explicitly.
//
// The template parameter 'P' gives the type for the priorities which should
// be comparable with '<' (only 'less' needed) parameter.  Currently we use
// 'double' and 'int64_t'. This container also serves as permanent store of
// priorities for all elements ever addded.
//
// Elements in the heap are indexed by non negative integers of type 'int'.
// Therefore the heap as well as the index to the position in heap tables
// only have 'int' elements.  This might prevent this template from being
// used in a more general context.
//
// In general, I think it is hard, if not impossible, to provide an
// efficient, yet simple to use non-intrusive priority queue with updates.
// Making it intrusive is an option, by adding 'piority' and 'position'
// fields to the elements.  However, then it will become harder to allow
// the same element to be in multiple priority queues.  Further, since we
// want to put variables into these queues, storing references to them
// in the queues requires to fix these references if the variable stack
// moves the addresses of the variables when growing due to new variables
// being added.  This was the way how 'Heap' worked in PrecoSAT.

template<class P> class priority_queue {

  vector<int> heap;	// This is the actual binary 'heap' data structure.

  vector<P> prior;	// Priorities of type 'P' comparable with '<'.

  vector<int> pos;	// Maps non negative integers to heap position.
  			// Negative position = not pushed on heap.

  // Elements are compared by comparing their priority.

  bool less (int a, int b) { return prior[a] < prior[b]; }

  // Standard 'heap' walking functions to compute the 'left', 'right',
  // resp. 'parent' position of a given position 'p' resp. 'c'.

  static int left (int p) { return 2*p + 1; }

  static int right (int p) { return 2*p + 2; }

  static int parent (int c) { assert (c > 0); return (c - 1)/2; }

  // The policy for updating the heap, by bubbling elements 'up and 'down,
  // is to only swap two elements if '<' (through the 'less' wrapper)
  // holds and the elements are in the wrong order.  In particular, if the
  // priority of two elements is the same, they are not swapped.

  void up (int e) { 				// Precolate element up.
    int epos = pos[e];
    while (epos > 0) {
      int ppos = parent (epos);
      int p = heap [ppos];
      if (!less (p, e)) break;
      heap [epos] = p, heap [ppos] = e;
      pos[p] = epos, epos = ppos;
    }
    pos[e] = epos;
  }

  void down (int e) {				// Precolate element down.
    assert (contains (e));
    int epos = pos[e], size = (int) heap.size ();
    for (;;) {
      int cpos = left (epos);
      if (cpos >= size) break;
      int c = heap [cpos], o;
      int opos = right (epos);
      if (!less (e, c)) {
	if (opos >= size) break;
	o = heap [opos];
	if (!less (e, o)) break;
	cpos = opos, c = o;
      } else if (opos < size) {
	o = heap [opos];
	if (!less (o, c)) { cpos = opos; c = o; }
      }
      heap [cpos] = e, heap [epos] = c;
      pos[c] = epos, epos = cpos;
    }
    pos[e] = epos;
  }

  // Only if an elment has been 'import'ed before (actually 'push'ed
  // or 'update'd) it has a definite priority.  Otherwise the priority is
  // unknown.  This is not precise, since importing 5 for instance also
  // imports implicitly 4 if it was not imported before.
  //
  bool imported (int e) const { 
    assert (0 <= e);
    return e < (int) pos.size (); 
  }

  // Make sure that this element is mapped to a priority.
  //
  void import (Mem & mem, int e) {
    while (!imported (e)) {
      pos.push_back (mem, -1);
      prior.push_back (mem, P());	// 'P' needs default constructor
    }
  }

public:

  void release (Mem & mem) {
    heap.release (mem);
    prior.release (mem);
    pos.release (mem);
  }

  bool empty () const { return heap.empty (); }

  size_t size () const { return heap.size (); }

  bool contains (int e) const {
    if (!imported (e)) return false;
    return pos[abs (e)] >= 0;
  }

  // This container also serves as a map of integer elements to priorities.
  // This map can only be changed by 'update'ing a specific priority of
  // an element.  The domain of the map can not be reduced, e.g. it is not
  // possible to 'unmap' an element.  See also 'update' below.
  //
  const P & priority (int e) const {
    assert (imported (e));
    return prior[abs (e)];
  }

  // Top element on the heap = largest element.
  //
  int top () const { return heap[0]; }

  // Push an element on to the heap either with the default priority 'P()'
  // or if it was on the heap before use the old priority for heap sorting.
  //
  void push (Mem & mem, int e) {
    assert (!contains (e));
    import (mem, e);
    pos[e] = (int) heap.size ();
    heap.push_back (mem, e);
    assert (heap[pos[e]] == e);
    up (e);
  }

  // Update priority of an elmenent.  There are no restrictions on whether
  // this has been 'push'ed or 'pop'ed.  If it is on the heap its position
  // is updated according to the different to the previous priority.  If it
  // is not on the heap or has never been 'push'ed before this function just
  // stores the new priority for this element.
  //
  // This is the only function for changing the priority of an element!
  //
  void update (Mem & mem, int e, const P & p) {
    import (mem, e);
    P q = prior[e];
    if (q == p) return;
    prior[e] = p;
    if (pos[e] < 0) return;
    if (p < q) down (e);
    if (q < p) up (e);
  }

  // Pop the given element from heap.  Its priority is kept as is.
  //
  void pop (int e) {
    assert (contains (e));
    int i = pos[e];
    pos[e] = -1;
    int last = heap.back ();
    heap.pop_back ();
    int j = (int) heap.size ();
    if (i == j) return;
    assert (i < j);
    pos[last] = i;
    heap[i] = last;
    up (last);
    down (last);
  }

  void pop () { pop (top ()); }

  // Rescore all priorities by multiplying them with the same factor.  This
  // is needed by 'Solver::rescore' and can be implemented faster internally
  // to the priority queue instead of going through 'update'.  An
  // alternative would be to pop and then push all elements (very slow) or
  // allow access to the heap, at least read access in the same order as the
  // elments occur on the heap (a constant iterator).
  //
  void rescore (double factor) {
    for (int i = 0; i < (int) prior.size (); i++)
      prior[i] *= factor;
  }
};

/*------------------------------------------------------------------------*/
/**** PART3: internal data structures for 'Solver' ************************/ 
/*------------------------------------------------------------------------*/

struct Option {
  const char * name;		// option name (= opts.#name)
  const char * description;	// usage information for this option
  int val;			// initially default then current value
  int min, max;			// legal option value range

  // Allows to use 'opts.name' instead of 'opts.name.val' for read access to
  // current value of the option.
  //
  operator int () const { return val; }

  Option () { }
  Option (const char * n, int v, int i, int s, const char * d) :
    name (n), description (d), val (v), min (i), max (s) { }
};

// Fore each 'OPTION(NAME,DEFAULT,MIN,MAX,DESCRIPTION)' line in 'options.h'
// add a 'NAME' member to 'Options' with exactly these fields.  Maybe it is
// possible to do this properly in C++.  However, we want to use the
// first argument 'NAME' of 'OPTION' both as string and as C++ symbol and
// this seems to me impossible without further shell scripts to process or
// generate 'options.h'.

struct Options {
  public:
#undef OPTION
#define OPTION(NAME,...) Option NAME;
#include "options.h"
  Options () {
#undef OPTION
#define OPTION(NAME,DEFAULT,MIN,MAX,DESCRIPTION) \
  NAME = Option (#NAME,DEFAULT,MIN,MAX,DESCRIPTION);
#include "options.h"
  }
};

/*------------------------------------------------------------------------*/

struct Clause {
  int glue;			// LBD = glue value as in Glucose

  bool redundant;		// Redundant (so not irredundant clause).
  bool remove;			// Remove this clause during 'reduce'.
  bool important;               // Do not reduce ('glue <= opts.gluekeep').
  bool forcing;			// Acting as reason for a literal.
  bool dumped;			// Dumped but not collected.
  bool satisfied;		// Top-level satisfied.

  int64_t activity;		// Measure of usefulness.

  // 'Inlined' and zero terminated list of literals for this clause.
  // Inlining means that 'lits' is actually of an abitrary size.  A zero
  // literal as sentinel is always added at the end, which avoids to store
  // the size explicitly.

  int lits[2];			// Of size 2 to avoid 'clang' warnings
  				// the actual size is 'size () + 1' to
				// include the zero sentinel.

  // As example, for an empty clause 'lits' is really just of size 1,
  // for a unit clauses 'lits' is of size 2, e.g. the unit and 0,
  // for a binary clause 'lits' is of size 3, etc.

  int size () const;		// Number of literals of this clause
  				// (without the '0' sentinel of course).

  bool large () const;		// Clause is large iff 'size () >= 3'.

  static int cmp (Clause **, Clause **);  // For 'sort' / 'qsort'.

  // The following two functions are used for monitoring memory usage.
  //
  size_t bytes () const { return  bytes (size ()); }
  static size_t bytes (int sz);

  void print () const;		// In DIMACS format for debugging only.
};

/*------------------------------------------------------------------------*/

int Clause::size () const {
  size_t res = 0;
  for (const int * p = lits; *p; p++) res++;
  return res;
}

bool Clause::large () const { return lits[0] && lits[1] && lits[2]; }

size_t Clause::bytes (int size) {
  return sizeof (Clause) + (size-1) * sizeof (int);
}

void Clause::print () const {
  for (const int * p = lits; *p; p++) cout << *p << ' ';
  cout << '0' << endl;
}

// Determines the clause evection strategy.  For more details see
// 'bump_clause' which uses the 'cbump' option, as well as initialization of
// 'glue' and 'activity' in 'new_clause', which also use the 'cbump' 'glured'
//
int Clause::cmp (Clause ** p, Clause ** q) {
  Clause * a = *p, * b = *q;
  if (a->glue < b->glue) return -1;
  if (a->glue > b->glue) return 1;
  if (a->activity < b->activity) return 1;
  if (a->activity > b->activity) return -1;
  return 0;
}

/*------------------------------------------------------------------------*/

struct Watch {			// Watcher of a literal in a clause.

  int blit;			// Blocking literal contained in 'clause'.
				// If this literal is true clause satisfied.

  bool binary;			// Flag is non zero if clause is binary.
				// Then clause 'clause' does not have to be
				// accessed in propagation.

  Clause * clause;		// Points to actual clause (even if binary).

  Watch (int l, int b, Clause * c) : blit (l), binary (b), clause (c) { }
};

typedef vector<Watch> Watches;	// Container of watches for a literal.

/*------------------------------------------------------------------------*/

// Full occurrence lists for dense mode.

struct Occs {
  int64_t count;		// Number of non garbage (dumped) clauses.
  vector<Clause*> clauses;
  Occs () : count (0) { }
  void release (Mem & mem) { clauses.release (mem); count = 0; }
  void dec () { assert (count > 0); count--; }
  void add (Mem & mem, Clause * c) { count++; clauses.push_back (mem, c); }

  // Just syntactic sugar for less code for traversing occurrences:
  //
  vector<Clause*>::const_iterator begin () const { return clauses.begin (); }
  vector<Clause*>::const_iterator end () const { return clauses.end (); }
  vector<Clause*>::iterator begin () { return clauses.begin (); }
  vector<Clause*>::iterator end () { return clauses.end (); }
};

/*------------------------------------------------------------------------*/

struct Var {			// Variable of index 1..maxvar.

  enum State { 
    FREE = 0,			// Not eliminated nor fixed.
    FIXED = 1,			// Fixed to a value on the top-level.
    ELIMINATED = 2,		// Eliminated through variable elimination.
  };

  State state;			// See above.

  int level;			// Decision level if assigned.
  int mark;			// Mark field for various traverals.

  Clause * reason;		// Forcing reason clause if assigned.

  Watches watches[2];		// Watches:     0=pos and 1=neg watches.
  Occs occs[2];			// Occurrences: 0=pos and 1=neg occurrences.

  Var () : 
    state (FREE),
    level (INT_MAX), 		// Needed for sorting watchers.
    mark (0), reason (0)
  { }

  bool free () const { return state == FREE; }

  void release (Mem & mem) {
    watches[0].release (mem);
    watches[1].release (mem);
    occs[0].release (mem);
    occs[1].release (mem);
  }
};

/*------------------------------------------------------------------------*/

struct Frame {			// Each decision opens up a new frame.
  int decision;			// Decision at this decision 'level'.
  int level;			// Decision 'level' of this frame.
  int trail;			// Trail height before this 'decision'.
  
  bool mark;			// Frame pulled into 1st UIP clause?
  				// Also used for calculating LBD = glue.

  Frame (int d = 0, int l = 0, int t = 0) : 
    decision (d), level (l), trail (t), mark (false) { }
};

/*------------------------------------------------------------------------*/

struct ClsPtrLitPair {
  Clause * cls;
  int lit;
  ClsPtrLitPair (Clause * c, int l) : cls (c), lit (l) { }
};

/*------------------------------------------------------------------------*/

struct Stats {

  int64_t conflicts;	// Number of leafs = conflicts during search.
  int64_t decisions;	// Number of decisions during search.
  int64_t levels;	// Sum of decision levels.
  int64_t iterations;	// Number of learned unit clauses during search.
  int64_t propagations;	// Number of propagated assignments.
  int64_t reductions;	// Number of reductions of learned clause cache.
  int64_t reported;	// Number of reported status messages.
  int64_t simplifications;  // Number of simplificaiton during inprocessing.
  int64_t steps;	    // Number of basic steps.

  int64_t sizes;	// Sum of learned clause sizes.

  struct {
    int64_t sum;	// Sum of glues.
    int64_t count;	// Number of glues summed up.
    int64_t updates;	// Number of updated glues.
  } glues;

  struct {
    int64_t count;	// Number of restarts.
    struct {
      int64_t count;	// Number of times part of the trail was reused.
      int64_t sum;	// Sum of the percentages of reused trail levels.
    } reuse;
  } restarts;

  struct {
    int64_t irredundant;// Current number of irredundant clauses.
    int64_t redundant;	// Current number of redundant clauses.
    int64_t collected;	// Accumulated number of collected clauses.
    int64_t reduced;	// Accumulated number of reduced redundant clauses.
    int64_t eliminated;	// Accumulated number of eliminated clauses.
    int64_t blocked;	// Accumulated number of blocked clauses.
  } clauses;

  struct {
    int64_t subsumed;	// Number of backward subsumed clauses.
    int64_t strengthened;	// Number of backward srengthened clauses.
  } backward;

  struct {
    int64_t units;       // Number of units found during distillation.
    int64_t subsumed;    // Number of subsumed clauses during distillation.
    int64_t strengthened;// Number of strengthened clauses in distillation.
  } distill;

  struct {
    int64_t fixed;	// Number of top-level assigned (fixed) variables.
    int64_t equivalent;	// Number of equivalent variables.
    int64_t eliminated;	// Number of eliminated variables.
  } vars;

  struct {		// Monitor effectiveness of clauses minimization.
    int64_t learned;	// Learned literals before minimization.
    int64_t minimized;	// Minimized literals during minimization.
  } lits;

  // Time spent in seconds in various routines for profiling purposes.
  // This is alwways enabled even in release code, but does not seem to
  // produce any measurable overhead (for 'getrusage' on Linux).
  //
  struct {

     struct {
       double analyze;
       double bias;
       double extend;
       double minimize;
       double reduce;
       double restart;
       double total;
     } search;

     struct {
       double block;
       double collect;
       double distill;
       double elim;
       double simplify;
       double total;
     } simp;

     double solve;
  } time;

  Stats () { memset (this, 0, sizeof *this); }	// Clear initially.
};

/*------------------------------------------------------------------------*/

struct Limits {
  int64_t restart;			// Restart after this conflict.
  int64_t max_restart_interval;	// Only 'report' restart if increased.

  struct {
    int64_t redundant;		// Reduce limit on redundant clauses.
    int64_t forcing;		// Number of forcing redundant clauses.
    int64_t important;		// Number of 'glue < = 2' clauses.
  } reduce;

  struct {

    // Gives the absolute number of removed variables of the last
    // simplification phase.  This is used as a metric for the
    // efficiency of simplification.  If simplification is very
    // efficient and removes many variables, then we reduce the
    // interval of conflicts allocated for the next search.
    //
    int removed_vars;

    // Bound simplification phase and stop it after that many steps.
    //
    int64_t steps;		
    int64_t inc;		// Next increment of 'steps' limit
    				// initialized by 'opts.simpint'.
  } simp;

  struct {

    // Bound search phase and stop when this number of conflicts is reached.
    //
    int64_t conflicts;
    int64_t inc;		// Next increment of 'conflicts' limit
    				// initialized by 'opts.searchint'.
  } search;

  Limits () { memset (this, 0, sizeof *this); }	// Clear initially.
};

/*------------------------------------------------------------------------*/

struct Timer {
  double * timeptr;		// Points to profiling statistics to update.
  double started;		// Old time, when this timer was started.

  Timer (double * tp, double s) : timeptr (tp), started (s) { }
};

/*------------------------------------------------------------------------*/

enum Simplifier { NOSIMP = 0, TOPSIMP = 1, BLOCK = 2, ELIM = 3, };

/*------------------------------------------------------------------------*/
/**** PART4: declaration of 'Solver' **************************************/ 
/*------------------------------------------------------------------------*/

class Solver : public AbstractSolver {

    Mem mem;			// Custom memory manager for memory stats.

    int level;			// Decision level (top/root level == 0).

    bool dense;			// Dense mode (full occurrences) or 
				// sparse mode (only watches).

    Simplifier simplifier;	// Points to current simplifier.

    bool schedule;		// Schedule (push) touched literals.

    int maxvar () const;	// Maximum variable index used/imported.

    void import (int lit);	// Make sure all data structures are
    				// ready to use this literal / variable.

    vector<Var> vars;		// Variable vector.
    vector<signed char> vals;	// Values seperated for fast access.
    vector<signed char> phases;	// Saved phases of previous assignment.

    void bias_phases ();	// Compute (initial) phase bias.

    int remaining_vars () const;  // Remaining variables for verbose output
    				  // without fixed, eliminated etc.

    Var & var (int lit);	  // Get variable of a literal.
    const Var & var (int) const;

    Watches & watches (int lit);		// Get literal watches.
    const Watches & watches(int lit) const;

    Occs & occs (int lit);			// Get full occurrences.
    const Occs & occs(int lit) const;

    Value val (int lit) const;	// Literal value ('lit' might be negative).

    priority_queue<double> decisions;

    vector<Frame> control;	// At each decisition a frame is added.
    vector<int> trail;		// To save assigned literals.
    size_t next;		// Position of next literal to propagate.
    size_t touched;		// Touched top-level fixed literal.

    void assign (int, Clause* c = 0);
    void unassign (int);

    void backtrack (int new_level = 0);	// Backtrack to to new level.

    Clause * empty;		// Points to an empty clause if non zero.
    vector<int> original;	// Save original clauses to check assignment.
    vector<int> addedlits;	// Save added literals before adding clause.
    vector<Clause*> clauses;

    LGL * lgl;			// Check entailment ('./configure.sh -c').
    bool clause_entailed ();	// Check clause 'addedlits' entailed.
    bool trivial_clause ();	// Check if clause 'addedlits' is non trivial
    				// and if not remove duplicated literals.

    bool each_variable_occurs_only_once ();	// .. in clause 'addedlits'.

    Clause * new_clause (bool redundant, int glue);	// Allocate.
    void connect_clause (Clause *);			// Connect watches.

    // Combines 'new_clause' and 'connect_clause' but also pushes it.
    //
    void new_push_connect_clause (bool redundant = false, int glue = 0);

    void dump_clause (Clause *);	// Decrease stats but do not delete.
    void delete_clause (Clause *);	// Physically delete clause memory.

    // Add occurrence of 'lit' with blocking literal 'blit'.
    //
    void add_watch (int lit, int blit, bool binary, Clause *);

    void touch (int lit);		// Update and push 'cands.*'.
    void add_occ (int lit, Clause * c);	// Add full occurrence of 'lit'.
    void dec_occ (int lit);		// Decrement full occurrence.

    void touch_fixed ();   // Touch literals in top-level satisfied clauses.

    void connect_occs (Clause * c);	// Connect occurrencees 
    					// all literals in this clause.

    void connect_occs ();	// Go to dense mode and connect
    				// occurrences of all irredundant clauses.

    vector<int> seen;		// Save marked literals as 'seen'.
    void unmark (int level=0);	// Unmark seen literals (after seen level).
    void mark (int);		// Mark literal and save it on 'seen'.
    int marked (int) const;	// Return signed mark (-1,0,1).

    Status simplify ();		// Bounded inprocessing.
    Status search ();		// Bounded search aka CDCL.

    bool propagated ();		// All assignments propagated.
    Clause * bcp ();		// Boolean contstraint propagation.
    void analyze (Clause *);	// Analyze conflict.
    void assume (int decision); // Assume decision / probe.
    bool decide ();		// Returns false if all assigned.

    void minimize_clause ();	// Minimize the 1st UIP clause.
    bool minimize_lit (int);	// Try to remove literal from 1st UIP clause.

    void restart ();		// Backtrack to the top.
    void new_restart_limit ();	// Set new restart limit.
    bool restarting () const;	// Restart limit hit.

    void reduce ();		// Reduce redudant clauses cache.
    bool reducing ();		// Reduce limit hit.

    vector<int> frames;		// Seen/marked frames during analyze.
    int unmark_frames ();	// Unmark all those seen frames.
    bool mark_frame (int lit);	// Mark this frame in analyze.
    bool pull_lit (int lit);	// Analyze this single literal.

    double score_increment;	// Current score increment.
    void bump_lit (int lit);	// Increase score (inverse of decay).
    void bump_clause (Clause *);// Increase activity of clause.
    void rescore (void);	// Rescore score of all literals.

    struct {
      priority_queue<int64_t> elim;	// Candidates to 'elim'inate.
      priority_queue<int64_t> block;	// Candidates to 'block'.
    } cands;

    priority_queue<int64_t> & current_cands ();	// Depends on 'simplifier'.
    void update_cands ();			// Update candidates.

    bool donotelim (int cand);	// Check for too many occurrences etc.
    void elim ();		// Bounded Variable Elimination (BVE).

    void block_lit (int blocking_lit);
    bool block_clause (Clause *, int blocking_lit);

    void block ();		// Blocked Clause Elimination (BCE).

    void backward (int lit);	    // Backward subsume from 'lit' clauses.
    void backward (Clause *, int);  // Backward subsume from clause.

    vector<ClsPtrLitPair> tostrengthen;	// Used in 'backward'.
    void strengthen (Clause *, int);	// Remove literal in clause.

    bool try_elim (int);	// Try BVE on this candidate.
    void do_elim (int);		// Do BVE on this candidate.

    bool try_resolve (Clause * c, int cand, Clause * d);
    void do_resolve (Clause * c, int cand, Clause * d);

    bool contains_eliminated (Clause *) const; // ... variable.
    void dump_eliminated_redundant ();	       // ... redundant clauses.

    vector<int> probes;			// Probes to try.
    void probing (int64_t steps);		// Simple failed literal probing.

    Clause * ignore;		// Ignore this and redundant clauses in BCP.
    int size_penalty () const;	// Reduce simp steps for large formulas.
    size_t distilled;		// Last distilled clause index.
    void distill ();		// Distill units, subsume/strengthen clauses.

    void disconnect_clauses ();		// Disconnect all clauses.
    Clause * reduce_clause (Clause *);	// Remove false literals.
    void collect_clauses ();		// Collect garbage (dumped) clauses.
    void connect_clauses ();		// Connect all clauses.
    void collect ();			// Garbage collection.

    vector<int> extension;	// Saved clauses to extend assignment.
    void push_extension (int);	// Push (and log) to extension.
    void extend ();		// Extend partial to full assignment.

    void push_extension (Clause *, int blit);

    void check_satisfied () const;   // Check assignment satisfies original.

    bool satisfied (Clause * c) const;
    bool satisfied (vector<int>::const_iterator &) const;

    void update_glue (Clause * c);  // Recompute and reduce glue of clause.

    vector<Timer> timers;
    void start (double *);	// Start timer and register time stats.
    void stop ();		// Stop it, update registered time stats.

    void report (char type);	// Progress report if verbosity is enabled.

    Stats stats;
    Limits limits;
    Options opts;

    void print () const;	// Print in DIMACS format to 'stdout'
    				// for debugging purposes.

    void print_time_for (double time, const char *); // Print timing stats.

  public:

    // Exactly the virtual functions in the interface from 'AbstractSolver'
    // make up the public interface of the concrete 'Solver' class.

    Solver ();
    ~Solver ();

    void addlit (int);

    void init_limits ();
    void update_limits ();

    Status solve ();
    Value value (int);

    void option (const char *, int);
    int option (const char *);

    void statistics ();
    size_t maximum_bytes ();
    size_t current_bytes ();
};

/*------------------------------------------------------------------------*/
}; // closing 'namespace Cleaneling'
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/**** PART5: implementation of private member functions of 'Solver' *******/ 
/*------------------------------------------------------------------------*/

int Solver::maxvar () const { 
  int res = (int) vars.size ();
  if (res) { assert (res > 1); res--; }	// Avoid under flow.
  return res;
}

Var & Solver::var (int lit) { 		// Get variable for this literal.
  int idx = abs (lit);
  assert (0 < idx && idx < (int) vars.size ());
  return vars[idx];
}

const Var & Solver::var (int lit) const { 
  int idx = abs (lit);
  assert (0 < idx && idx < (int) vars.size ());
  return vars[idx];
}

// Signed current value, that is, if 'lit' is negative the value of the
// variable is negated to get the value of the literal.  The functions
// 'assign' resp. 'unassign' are used to set resp. reset the value.
//
Value Solver::val (int lit) const {
  int res = vals [abs (lit)];
  if (lit < 0) res = -res;
  return (Value) res;
}

// Signed mark flag, as in 'val' above, i.e. if 'lit' is negative the
// negation of the mark field of the variable is returned.
//
int Solver::marked (int lit) const {
  int res = (int) var (lit).mark;
  return (lit < 0) ? -res : res;
}

// Mark a variable either with '1' resp. '-1' if 'lit' is pos. resp. neg.
// Marked literals are saved on the 'seen' stack, such they can easily be
// unmarked through 'unmark'.
//
void Solver::mark (int lit) {
  Var & v = var (lit);
  assert (!v.mark);
  v.mark = sign (lit);
  seen.push_back (mem, lit);
}

// Unmark all variables previously marked and saved on 'seen'.
//
void Solver::unmark (int level) {
  assert (level <= (int) seen.size ());
  while (level < (int) seen.size ()) {
    int lit = seen.back ();
    seen.pop_back ();
    Var & v = var (lit);
    assert (v.mark == sign (lit));
    v.mark = 0;
  }
}

/*------------------------------------------------------------------------*/

// Compute an initial phase bias based on the Jeroslov-Wang heuristic.

void Solver::bias_phases () {
  start (&stats.time.search.bias);
  size_t bytes = 2 * (maxvar() + 1) * sizeof (double);
  double * score = (double*) mem.malloc (bytes);
  if (!score) internal_error ("out of memory allocating score array");
  memset (score, 0, bytes);
  score += maxvar ();
  vector<Clause*>::const_iterator i;
  for (i = clauses.begin (); i != clauses.end (); i++) {
    Clause * c = *i;
    if (c->redundant || c->satisfied) continue;
    double inc = 1;  // compute 2^-size
    for (int size = c->size (); size; size--) inc /= 2.0;
    for (const int * p = c->lits; *p; p++) score[*p] += inc;
  }
  for (int idx = 1; idx <= maxvar (); idx++)
    phases[idx] = ((score[idx] > score[-idx]) ? 1 : -1);
  score -= maxvar ();
  mem.free (score, bytes);
  stop ();
}

/*------------------------------------------------------------------------*/

// We can use an external SAT solver (currently only Lingeling is
// supported), to check whether learned clauses and resolvents are implied
// by the original formula.  This has to be configured at compile time
// either in debug mode or with './configure.sh -c'.  In any case
// the Lingeling library has to accessible in '../lingeling' and checking
// will only be performed if you set the 'check' option.

bool Solver::clause_entailed () {
#ifndef NLGL
  if (!opts.check) return true;
  vector<int>::const_iterator i;
  for (i = addedlits.begin (); i != addedlits.end (); i++)
    lglassume (lgl, -*i);
  bool res = (lglsat (lgl) == 20);
  bool print = false;
  if (!res) print = true;
#ifndef NLOG
  if (opts.log) print = true;
#endif
  if (print) {
    cout << "c @" << level << " checking entailment of clause";
    for (i = addedlits.begin (); i != addedlits.end (); i++)
      cout << ' ' << *i;
    cout << (res ? " succeeded" : " failed") << endl;
  }
  return res;
#else
  return true;
#endif
}

/*------------------------------------------------------------------------*/

// Used in debugging for checking that the variables in 'addedlits' occur
// at most once.  For external clauses this ensured by 'trivial_clause'.
//
bool Solver::each_variable_occurs_only_once () {
  int lit;
  vector<int>::const_iterator p;
  for (p = addedlits.begin (); p != addedlits.end (); p++)
    assert (!marked (lit = *p)), assert (!marked (-lit));
  bool res = true;
  for (p = addedlits.begin (); res && p != addedlits.end (); p++)
    if (marked (lit = *p)) res = false;
    else if (marked (-lit)) res = false;
    else mark (lit);
  unmark ();
  return res;
}

// Check if clause is non trivial and if not remove duplicated literals.
//
bool Solver::trivial_clause () {
  vector<int>::iterator q = addedlits.begin ();
  vector<int>::const_iterator p;
  int lit, m;
  bool res = false;
  for (p = q; !res && p != addedlits.end (); p++) {
    lit = *p;
    assert (lit);
    m = marked (lit);
    if (m < 0) res = true; 
    else if (!m) *q++ = lit, mark (lit);
  }
  addedlits.resize (q - addedlits.begin ());
  unmark ();
  return res;
}

Watches & Solver::watches (int lit) { return var (lit).watches[lit < 0]; }

const Watches &
Solver::watches (int lit) const { return var (lit).watches[lit < 0]; }

// Add watch for 'lit' in 'clause' to watcher list of 'lit'.  The size,
// clause and blocking literal are stored as is.  See description of 'Watch'
// above or 'new_clause' below for more details.
//
void Solver::add_watch (int lit, int blit, bool binary, Clause * clause) {
  watches (lit).push_back (mem, Watch (blit, binary, clause));
}

Occs & Solver::occs (int lit) { return var (lit).occs[lit < 0]; }

const Occs &
Solver::occs (int lit) const { return var (lit).occs[lit < 0]; }

void Solver::touch (int lit) {
  if (lit < 0) lit = -lit;
  int64_t new_priority = occs (lit).count + occs (-lit).count;
  if (!var (lit).free ()) return;
  priority_queue<int64_t> & queue = current_cands ();
  queue.update (mem, lit, -new_priority);
  if (schedule && !queue.contains (lit)) {
    LOG ("push " << lit << " priority " << -new_priority);
    queue.push (mem, lit);
  }
}

// Further, touch literals in newly top-level satisfied clauses using the
// literals not touched yet in the trail starting at 'touched'.  This global
// trail pointer avoids going through all top-level literals repeatedly.
//
void Solver::touch_fixed (void) {
  assert (dense);
  assert (!level);
  assert (schedule);
  while (touched < trail.size ()) {
    int lit = trail[touched++];
    assert (val (lit) > 0);
    assert (!var (lit).level);
    Occs & os = occs (lit);
    vector<Clause *>::const_iterator i;
    for (i = os.begin (); i != os.end (); i++) {
      Clause * c = *i;
      int other;
      for (const int * p = c->lits; (other = *p); p++) {
	if (val (other)) continue;
	assert (other != lit);
	touch (other);
      }
    }
  }
}

void Solver::add_occ (int lit, Clause * c) {
  assert (dense);
  occs (lit).add (mem, c);
  touch (lit);
}

void Solver::connect_occs (Clause * c) {
  assert (dense);
  for (const int * p = c->lits; *p; p++)
    add_occ (*p, c);
}

void Solver::connect_occs () {
  assert (!dense);
  dense = true;
  Clause * c;
  vector<Clause*>::const_iterator i;
  for (i = clauses.begin (); i != clauses.end (); i++)
    if (!(c = *i)->redundant) connect_occs (c);
}

void Solver::dec_occ (int lit) {
  assert (dense);
  occs (lit).dec ();
  touch (lit);
}

/*------------------------------------------------------------------------*/

// Assign literal 'lit' with 'reason' as reason clause .  This assignment
// marks 'lit' as a decision if 'reason' is zero.  Otherwise 'reason' is
// marked as being used as reason to avoid garbage collection during
// 'collect_clauses' and 'reduce'.  Further save this literal as being
// assigned on the 'trail' for backtracking.
//
void Solver::assign (int lit, Clause * reason) {
  Var & v = var (lit);
  assert (!val (lit));
  int idx = abs (lit);
  Value s = (Value) sign (lit);
  vals[idx] = s;
  if (simplifier == NOSIMP) phases[idx] = s;
  if (!(v.level = level)) {
    LOG ("fixed " << lit);
    stats.vars.fixed++;
    if (v.state == Var::ELIMINATED) {
      assert (stats.vars.eliminated > 0);
      stats.vars.eliminated--;
    } else assert (v.state == Var::FREE);
    v.state = Var::FIXED;
  }
  trail.push_back (mem, lit);
  if ((v.reason = reason)) {
    assert (!reason->forcing);
    reason->forcing = true; // Mark used for 'reduce' and 'collect_clauses'.
    // Then increase the redundant cache size temporarily if necessary.
    if (reason->redundant && !reason->important) limits.reduce.forcing++;
    LOGCLAUSE (reason, "assign " << lit << " forced by");
    update_glue (reason);
  } else LOG ("assign " << lit << " decision");
}

// Unassign a previously assigned literal.  Add it back to 'decisions'.
//
void Solver::unassign (int lit) {
  assert (level > 0);
  Clause * reason;
  Var & v = var (lit);
  assert (val (lit) == TRUE);
  vals [abs (lit)] = UNASSIGNED;
  assert (v.level == level);
  v.level = INT_MAX;
  if ((reason = v.reason)) {
    assert (reason->forcing);
    reason->forcing = false;
    if (reason->redundant && !reason->important) {
      assert (limits.reduce.forcing > 0);
      limits.reduce.forcing--;
    }
  }
  LOG ("unassign " << lit);
  int idx = abs (lit);
  if (!decisions.contains (idx)) decisions.push (mem, idx);
}

// Backtracking involves popping frames from the control stack, unassigning
// all variables in these frames until the last frame still on the control
// stack has the given new decision level.  For decision level 0 we
// initially pushed an empty frame on the control stack, which here serves
// as sentinel.
//
void Solver::backtrack (int new_level) {
  assert (0 <= new_level && new_level <= level);
  if (new_level == level) return;
  LOG ("starting to backtrack to level " << new_level);
  const Frame * f;
  while ((f = &control.back ())->level > new_level) {
    assert (f->level == level);
    assert (f->trail < (int) trail.size ());
    while (f->trail < (int) trail.size ()) {
      int lit = trail.back ();
      assert (var (lit).level == f->level);
      trail.pop_back ();
      unassign (lit);
    }
    assert (level > 0);
    level--;
    LOG ("pop frame decision " << f->decision << " trail " << f->trail);
    trail.resize (f->trail);
    next = f->trail;
    control.pop_back ();
  }
  LOG ("finished backtracking to level " << new_level);
  assert (new_level == level);
}

/*------------------------------------------------------------------------*/

// Allocate, and copy 'addedlits' of previously added literals stored in
// 'addedlits'.  This routine also backtracks to an apropriate decision
// level in which the clause is not empty given the current assignment. If
// it is a unit it will actually also assign the unit literal.  If the
// clause is empty even at the top-level it will then set the global 'empty'
// pointer to this clause if not already an empty clause has been found
// earlier.  Backtracking is actually back jumping to the second largest
// decision level if the clause contains at least two literals.
//
Clause * Solver::new_clause (bool redundant, int glue) {
  size_t size;
  Clause * c;

  assert (!glue || redundant);
  assert (each_variable_occurs_only_once ());

  size = addedlits.size ();

  // Allocate and adjust stats (see comments on 'inlining lits' above).
  {
    size_t bytes = Clause::bytes (size);
    c = (Clause *) mem.malloc (bytes);
    if (!c) internal_error ("out of memory while allocating clause");
#ifdef __clang_analyzer__
    if (!c) exit (1);
#endif
    memset (c, 0, bytes);
    if (opts.gluered) c->glue = glue;	// Glue based reduce, so save glue.
    else assert (!c->glue);		// Otherwise all clauses same glue.

    // Always use 'glue' here and keep all 'glue <= opts.gluekeep' clauses.
    //
    c->important = (glue <= opts.gluekeep);	
					  
    c->redundant = redundant;
    c->activity = stats.conflicts;
  }

  // Copy literals from 'addedlits' to the result.
  {
    int * p = c->lits;
    vector<int>::const_iterator i;
    for (i = addedlits.begin (); i != addedlits.end (); i++) *p++ = *i;
    *p = 0;
    assert (!*p);
  }

  // Update clause statistics and log it.
  //
  if (redundant) {
    stats.clauses.redundant++;
    LOGCLAUSE (c, "glue " << glue << " new");
  } else {
    stats.clauses.irredundant++;
    LOGCLAUSE (c, "new");
  }

#ifndef NLGL
  // Add clause to Lingeling for checking entailment of derived clauses.
  //
  if (opts.check) {
    for (const int * p = c->lits; *p; p++) lgladd (lgl, *p);
    lgladd (lgl, 0);
  }
#endif

  return c;
}

void Solver::new_push_connect_clause (bool redundant, int glue) {
  Clause * c = new_clause (redundant, glue);
  clauses.push_back (mem, c);
  connect_clause (c);
}

/*------------------------------------------------------------------------*/

void Solver::connect_clause (Clause * c) {

  if (c->satisfied) return; // do not connect top level satisfied clauses

  size_t size = c->size ();
  bool binary = (size == 2);

  // Swap two literals with maximum level to the front
  // (unassigned variables have 'infinite' level of 'INT_MAX')
  {
    int lit;
    for (int * p = c->lits; *p && p <= c->lits + 1; p++)
      for (int * q = p + 1; (lit = *q); q++) {
	int cmp = (var (lit).level - var (*p).level);
	if (cmp > 0 || (!cmp && (val (lit) > val (*p))))
	  swap (*p, *q);
      }
  }

  int l0 = c->lits[0], l1 = l0 ? c->lits[1] : 0;

  // Backtrack to the minimum of the levels of these two first literals.
  {
    int new_level = (l0 && l1) ? std::min (var (l0).level, var (l1).level) : 0;
    if (new_level != INT_MAX) backtrack (new_level);
  }

  // Now watch the first two literals if both are non zero, so the clause is
  // at least a binary clause.  See comments to 'Watch' above for details.
  //
  if (size >= 2) {
    add_watch (l0, l1, binary, c);
    add_watch (l1, l0, binary, c);
    if (dense) connect_occs (c);
  }

  // Finally check whether this clause is a unit clause or even empty.
  // If so assign the unit respectively set 'empty' (if still zero).
  { 
    bool ignore = false;
    int lit = 0, other, tmp;
    for (const int * p = c->lits; (other = *p) && p <= c->lits + 1; p++) {
      if ((tmp = val (other)) == TRUE) { ignore = true; break; }
      if (tmp == UNASSIGNED) { if (lit) ignore = true; else lit = other; }
    }
    if (!ignore) {
      if (!lit) {
	assert (!level);
	LOG ("found empty clause");
	if (!empty) empty = c;
      } else assign (lit, c);	// Assign unit forced by this clause.
    }
  }

  // Maintain the number of connected important (currently 'glue <= 2')
  // clauses for delaying 'reduce' in case there are too many of them.
  //
  if (c->redundant && c->important) limits.reduce.important++;
}

/*------------------------------------------------------------------------*/

// Pretend that this clause does not exist anymore. Do not actually delete
// the clause now but delay its deletion until garbage collection.  In other
// solvers like Lingeling, eager deletion is used, which forces costly
// disconnection of  watches and full occurrences.  So this is a space for
// time trade-off here.
//
void Solver::dump_clause (Clause * c) {
  if (c->dumped) return;
  if (c->redundant) {
    assert (stats.clauses.redundant > 0);
    stats.clauses.redundant--;
  } else  {
    assert (stats.clauses.irredundant > 0);
    stats.clauses.irredundant--;
    if (dense)
      for (const int * p = c->lits; *p; p++)
	dec_occ (*p);
  }
  c->dumped = true;
  LOGCLAUSE (c, "dump");
}

/*------------------------------------------------------------------------*/

void Solver::delete_clause (Clause * c) {
  dump_clause (c);
  LOGCLAUSE (c, "delete");
  mem.free (c, c->bytes ());
}

/*------------------------------------------------------------------------*/

// Start timer by saving current time and pointer to the right statistics.
//
void Solver::start (double * time_ptr) {
  timers.push_back (mem, Timer (time_ptr, seconds ()));
}

// Stop the (current = inner most) timer and update its statistics.
//
void Solver::stop () {
  Timer timer = timers.back ();
  timers.pop_back ();
  double delta = seconds () - timer.started;
  if (delta < 0) delta = 0;
  *timer.timeptr += delta;
}

// Imports the given literal and makes sure that all data structures, in
// particularly those involved to represent variables are set up to use this
// literal (and its negation).
//
void Solver::import (int lit) {
  int idx = abs (lit), new_idx;
  assert (lit);
  while (idx >= (new_idx = (int) vars.size ())) {
    vars.push_back (mem, Var ());
    vals.push_back (mem, 0);
    phases.push_back (mem, 0);
    if (!new_idx) continue;
    decisions.push (mem, new_idx);
#ifndef NLGL
    lglfreeze (lgl, new_idx);
#endif
  }
}

// All assigned literals propagated.
//
bool Solver::propagated () { return next == trail.size (); }

// Boolean constraint propagation (BCP) is one of the most expensive
// functions in SAT solvers.  We use blocking literals to reduce memory
// accesses, which also gives us faster propagation for binary clauses.
//
Clause * Solver::bcp () {
  int64_t visits = 0, propagations = 0;
  Clause * conflict = empty;
  while (!conflict && next < trail.size ()) { // BFS over assigned literals.
    propagations++;
    int lit = -trail[next++];                 // Check negative occurrences.
    stats.steps++;
    LOG ("propagating " << -lit);
    Watches & ws = watches (lit);
    vector<Watch>::const_iterator i;
    vector<Watch>::iterator j = ws.begin ();
    for (i = j; !conflict && i != ws.end (); i++) {
      Watch w = *j++ = *i;		// By default copy watch.
      int other = w.blit;		// Invariant: blit occurs in clause.
      Value v = val (other);		// Value of blocking literal.
      if (v == TRUE) continue;		// Blocking literal is true.
      Clause * clause = w.clause;
      if (ignore) {				// If 'ignore is non-zero
	 if (ignore == clause) continue;	// ignore 'ignore' clause
	 if (clause->redundant) { 		// and redundant clauses.
	   if (!w.binary) visits++;		// Count as visited clause.
	   continue;
	 }
      }
      if (w.binary) {				// No need to access clause.
	if (v == FALSE) conflict = clause;	// Other literal false too.
	else assign (other, clause);            // Other literal unit.
      } else {
	visits++;
	if (clause->dumped) { j--; continue; }	// Drop watch if dumped.
	int * lits = clause->lits, * p;
	if (lits[0] == lit) swap (lits[0], lits[1]); // Normalize lit pos:
	assert (lits[1] == lit);                     // simpler code below.
	for (p = lits + 2; (other = *p); p++) // Search first non-false
	  if ((v = val (other)) >= 0) break;  // literal as replacement watch.
	if (other) {
	  *p = lit, lits[1] = other;                 // Swap old with new.
	  add_watch (other, lits[0], false, clause); // Add new watch.
	  j--;					     // Drop old watch.
	} else {                      // No non-false replacement found.
	  v = val (other = lits[0]);                  // Get other watch.
	  if (v == FALSE) conflict = clause;	      // All lits false.
	  else if (v != TRUE) assign (other, clause); // Other watch unit.
          else j[-1].blit = other;                    // Optimistic update.

	}
      }
    }
    // If there is a conflict we might have bailed out of the loop early and
    // thus still have to copy the remaining watches before size reset.
    //
    if (conflict) while (i != ws.end ()) *j++ = *i++;
    ws.resize (j);					// Reset size.
  }
  if (conflict) {
    LOGCLAUSE (conflict, "conflicting");
    if (!simplifier) stats.conflicts++;   // Otherwise counted differently.
  }
  stats.propagations += propagations;
  stats.steps += visits;
  return conflict;
}

// Mark this frame as being used to compute the LBD = glue value for
// better predicting learned clause usefullness as in Glucose.  The mark
// flag is also used for clause minimization (see also the comments in
// 'minimize_lit').
//
bool Solver::mark_frame (int lit) {
  int level = var (lit).level;
  Frame & frame = control[level];
  if (frame.mark) return false;
  frame.mark = true;
  frames.push_back (mem, level);
  return true;
}

// Unmark seen frames and return their number which is the 'glue' (LBD).
//
int Solver::unmark_frames () {
  int res = (int) frames.size ();
  while (!frames.empty ()) {
    Frame & f = control[frames.back ()];
    frames.pop_back ();
    assert (f.mark);
    f.mark = 0;
  }
  return res;
}

// Rescore all scores and the score_increment to be at most 1.
// Make sure to keep the relative order of the priorities.
//
void Solver::rescore () {
  double max_score = score_increment;
  for (int idx = 1; idx < (int) vars.size (); idx++) {
    double p = decisions.priority (idx);
    if (p > max_score) max_score = p;
  }
  double factor = 1/max_score;
  LOG ("rescoring by " << factor);
  decisions.rescore (factor);
  score_increment *= factor;
}

// Bump (increase) the EVSIDS score of a literal using the, the exponential
// variant of VSIDS as first implemented in MiniSAT.
//
void Solver::bump_lit (int lit) {
  const double max_priority = 1e300;
  int idx = abs (lit);
  double old_priority;
  if (score_increment > max_priority ||			// Avoid overflow.
      (old_priority = decisions.priority (idx)) > max_priority) {
    rescore ();
    old_priority = decisions.priority (idx);
  }
  double new_priority = old_priority + score_increment;
  decisions.update (mem, idx, new_priority);
  LOG ("bumping score of " << idx  << " from " 
       << scientific << setprecision (3)
       << old_priority << " to " << new_priority);
}

// If the given literal has already been analyzed just return.  Otherwise
// mark it as having been analyzed and bump its score.  If then the literal
// is assigned on the current level return true to tell the caller that
// another literal on the current level is marked.  Finally if the literal
// has not been marked before and was assigned on a smaller decision level
// add it to the learned clause.
//
bool Solver::pull_lit (int lit) {
  if (val (lit) == TRUE) return false;
  if (marked (lit)) return false;
  mark (lit);
  bump_lit (lit);
  if (var (lit).level == level) return true;
  if (mark_frame (lit)) LOG ("pulling frame " << var (lit).level);
  LOG ("adding literal " << lit);
  addedlits.push_back (mem, lit);
  return false;
}

// Bump (increase) activity of a clause.  This depends on the scheme
// selected by the option 'cbump'.  See the cases below.  If reduction is
// only based on 'activity', by disabling the 'gluered' option, the activity
// updated here and its initialization to the number of conflicts during
// clause construction determine which clauses are evicted from the 'cache'
// of redundant clauses.
//
void Solver::bump_clause (Clause * c) {
  switch (opts.cbump)
    {
      case 0:		// 0=none
      default:

        assert (!opts.cbump);	// Make sure we cover all cases.

	// No Bumping at all.  Clauses are sorted chronologically, since
	// they have their activity set to the number of conflicts at the
	// point they were added.
	break;

      case 1:		// 1=inc

	// Number of times a clauses was involved in conflict analysis.

	c->activity++;
	break;

      case 2:		// 2=lru

	// Least recently used conflict clauses are kicked out.

	c->activity = stats.conflicts;
	break;

      case 3:		// 3=avg

	// Running average of most recently used.

	c->activity = (stats.conflicts + c->activity) / 2;
	break;
    }
}

// Minimize the first UIP clause by trying to remove this literal.  If this
// lieral can be removed true is returned and all literals needed in
// recursive resolutions to remove it are marked.
//
// This test is a simple DFS which starts at the antecedent literals of the
// given starting literal.  It goes recursively backward over reason clauses
// and is successful if DFS eventually reaches marked literals and no
// unmarked decision.
//
bool Solver::minimize_lit (int root) {
  assert (marked (root));
  Clause * reason = var (root).reason;
  if (!reason) return false;
  size_t old_seen_size = seen.size (), next_seen = old_seen_size;
  bool res = true;
  int lit = root;
  for (;;) {
    int other;
    for (const int * p = reason->lits; res && (other = *p); p++) {
      if (other == lit) continue;
      assert (val (other) == FALSE);
      if (marked (other)) continue;
      Var & v = var (other);
      if (!v.reason) res = false;

      // The original MiniSAT code uses Bloom filters, made of 64 bit hashes
      // of the levels in the 1st UIP clause, as was done in Quantor and
      // SatElite for faster subsumption checks by computing signatures over
      // respectively hashing literals of clauses.  However, since we want
      // fast access to the decision literals for the 'reuse trail'
      // technique anyhow, and probably because it is faster too, we use a
      // seperate control stack with frame entries for each decision level
      // These frames have additional mark flags to remember whether a
      // particular decision level occurs in the 1st UIP clause.  The third
      // use of these flags is to compute the Glue/LBD of a learned clause.

      else if (!control[v.level].mark) {
	res = false;
      } else mark (other);
    }
    if (!res || next_seen == seen.size ()) break;
    lit = -seen[next_seen++];
    reason = var (lit).reason;
    assert (reason);
  } 
  if (res) LOG ("minimizing " << lit);
  else unmark (old_seen_size);
  return res;
}

// Minimization  of the learned first UIP clause.
//
void Solver::minimize_clause () {
  LOG ("minimize 1st UIP clause");
  start (&stats.time.search.minimize);
  size_t learned = addedlits.size ();
  stats.lits.learned += learned;

  vector<int>::const_iterator i;
  vector<int>::iterator j = addedlits.begin ();
  for (i = j; i != addedlits.end (); i++)
    if (!minimize_lit (-*i)) *j++ = *i;
  addedlits.resize (j);

  stats.lits.minimized += learned - addedlits.size ();
  stop ();
}

// Main analysis function, which learns a first UIP clause from the
// given conflict.  It traverses in reverse assignment order reasons
// of literals assigned on the current decision level which are backward
// reachable in the implication graph from the conflict.
//
void Solver::analyze (Clause * reason) {
  if (empty) { assert (!level); return; }
  start (&stats.time.search.analyze);
  update_glue (reason);
  assert (addedlits.empty ());
  vector<int>::const_iterator it = trail.end ();
  int lit = 0, open = 0;
  for (;;) {
    LOGCLAUSE (reason, "analyzing");
    bump_clause (reason);
    for (const int * p = reason->lits; (lit = *p); p++)
      if (pull_lit (lit)) open++;
    while (it != trail.begin() && !marked (lit = -*--it))
      assert (var (lit).level == level);
    if (it == trail.begin () || !--open) break;
    reason = var (lit).reason;
    assert (reason);
  }
  assert (lit);
  LOG ("adding UIP " << lit);
  addedlits.push_back (mem, lit);
  minimize_clause ();
  unmark ();			// Flush 'seen' and unmark all literals.
  assert (clause_entailed ());
  int glue = unmark_frames ();
  stats.glues.count++;
  stats.glues.sum += glue;
  stats.sizes += addedlits.size ();
  new_push_connect_clause (true, glue);
  addedlits.clear ();
  score_increment *= opts.scincfact/1000.0;

  // Learning a new unit in the CDCL loop starts a new 'iteration'.
  // Thus the number of iterations is the number of unit clauses learned
  // by the conflict analysis during CDCL.
  //
  if (!simplifier && !level && !empty) {

    // Since a new unit was learned and thus 'search' seems to be effective
    // delay next simplification a little bit.
    //
    limits.search.conflicts += opts.itsimpdel;

    stats.iterations++;				
    report ('i');
  }

  stop ();
}

/*------------------------------------------------------------------------*/

int Solver::remaining_vars () const {		// Actually remaining.
  int res = maxvar ();
  res -= stats.vars.fixed;
  res -= stats.vars.equivalent;
  res -= stats.vars.eliminated;
  return res;
}

/*------------------------------------------------------------------------*/

// Provide progress 'report' in verbose mode.
//
void Solver::report (char type) {
  if (!opts.verbose) return;

  // Print the header every 22nd row such that on a 25 row terminal we
  // keep on seeing the header all the time.
  //
  if (!(stats.reported++ % 20)) {
    print_line ();
    cout << 
"c    seconds     irredundant      limit         height         glue" 
<< endl << 
"c          variables     conflicts    redundant         size           MB" 
<< endl;
    print_line ();
  }

  cout 
  << "c " << type 
  << ' '
  << setprecision(1) << setw (7) << fixed << seconds ()
  << ' '
  << FmtNat (remaining_vars ())
  << ' '
  << FmtNat (stats.clauses.irredundant)
  << ' '
  << FmtNat (stats.conflicts)
  << ' '
  << FmtNat (limits.reduce.redundant)
  << ' '
  << FmtNat (stats.clauses.redundant)
  << ' '
  << setprecision(1) << setw(6) << average (stats.levels, stats.decisions)
  << ' '
  << setprecision(1) << setw(6) << average (stats.sizes, stats.conflicts)
  << ' '
  << setprecision(1) << setw(6) << 
     average (stats.glues.sum, stats.glues.count)
  << ' '
  << setprecision(1) << setw(6) << (current_bytes ()/(double)(1<<20))
  << endl << flush;
}

/*------------------------------------------------------------------------*/

// Compute the next conflict number at which a restart is triggered.  This
// is Luby controlled (rapid restarts) with basic interval 'restartint'.
//
void Solver::new_restart_limit () {
  int64_t new_interval = opts.restartint * luby (stats.restarts.count + 1);
  bool doreport = (stats.restarts.count == 1);

  if (new_interval > limits.max_restart_interval) {
    limits.max_restart_interval = new_interval;
    doreport = true;
  }

  limits.restart = stats.conflicts + new_interval;

  // Do not report every 'restart', but only those where the interval to the
  // next restart reached a new maximum (and do not miss the first one).
  //
  if (doreport) report ('R');
}

bool Solver::restarting () const {
  if (!opts.restart) return false;
  return stats.conflicts >= limits.restart;
}

void Solver::restart () {
  start (&stats.time.search.restart);
  LOG ("restart");
  stats.restarts.count++;

  // Get the next decision by popping from the 'decisions' priority queue as
  // long the top is still assigned.  This is the work 'decide' would need
  // to do anyhoze. We have to keep the first unassigned variable on the
  // decision queue.  This makes sensce for both cases, if 'reusetrail' is
  // enabled or not.
  //
  int next_decision = 0;
  while (!next_decision && !decisions.empty ()) {
    int lit = decisions.top ();
    if (val (lit)) decisions.pop (lit);
    else next_decision = lit;
  }

  if (next_decision) {		// Is there an unassigned variable left?

    int new_level;

    if (opts.reusetrail) {

      // Find the outermost decision with smaller priority.  This decision
      // has to be backtracked over, but the rest of the trail can be
      // reused.  This code is based on the idea of 'Reusing the Trail' by
      // Marijn Heule and his students as published in POS'11.

      double next_decision_priority = decisions.priority (next_decision);
      for (new_level = 0; new_level < level; new_level++) {
	Frame & frame = control[new_level + 1];
	int decision = abs (frame.decision);
	if (decisions.priority (decision) < next_decision_priority) {
	  // At level 'new_level + 1' there is a 'decision' with smaller
	  // priority than the 'next_decision'.  So backtrack to decision
	  // level 'new_level' in order to use 'next_decision' as decision
	  // at that level instead.
	  break;
	}
      }

      if (new_level) {
	stats.restarts.reuse.count++;
	stats.restarts.reuse.sum += (100*new_level) / level;
      }
    } else new_level = 0; // Do not reuse trail.  Back track to the top.

    // If all previous decisions have smaller priority than 'next_decision'
    // then 'new_level' will reach 'level' and no restart happens.
    //
    backtrack (new_level);
  }

  new_restart_limit ();		// Update restart limit.
  stop ();
}

/*------------------------------------------------------------------------*/

// Check if the number of redundant clauses without the forcing redundant
// clauses is larger than the reduce limit.
//
bool Solver::reducing () {
  if (!limits.reduce.redundant) limits.reduce.redundant = opts.redinit;
  int64_t limit = limits.reduce.redundant;
  limit += limits.reduce.forcing;
  limit += limits.reduce.important;
  return limit <= stats.clauses.redundant;
}

// This function removes half of the redundant clauses perdiodically and
// thus helps to speed up BCP, which otherwise would be slowed down by too
// many useless redundant clauses.
//
void Solver::reduce () {

  start (&stats.time.search.reduce);

  LOG ("reduce");
  stats.reductions++;

  {
    Clause * c;

    // First store all redundant and non-forcing clauses in 'candidates'.
    //
    vector<Clause *> candidates;
    vector<Clause*>::const_iterator i;
    for (i = clauses.begin (); i != clauses.end (); i++)
      if ((c = *i)->redundant && !c->important && !c->forcing)
	candidates.push_back (mem, c);

    size_t keep = candidates.size ()/2;		// Target is to remove half.

    LOG ("reduction " << stats.reductions << " keeps " <<
         keep << " clauses " << " out of " << candidates.size ());

    // Now sort with respect to glue and activity.
    //
    candidates.sort (Clause::cmp);

    // Mark those that are going to be removed as garbage.
    //
    for (size_t i = keep; i < candidates.size (); i++)
      candidates[i]->remove = true;

    // Flush references to garbage clauses from watches.
    //
    for (int idx = 1; idx <= maxvar (); idx++)
      for (int sign = -1; sign <= 1 ;sign += 2) {
	Watches & ws = watches (sign * idx);
	vector<Watch>::const_iterator i;
	vector<Watch>::iterator j = ws.begin ();
	for (i = j; i != ws.end (); i++)
	  if (!(*i).clause->remove) *j++ = *i;
	ws.resize (j);
      }

    // Flush references to garbage clauses from 'clauses'.
    {
      size_t i, j = 0;
      for (i = 0; i < clauses.size (); i++) {
	c = clauses[i];
	if (i == distilled) distilled = j;
	if (!c->remove) clauses[j++] = c;
      }
      if (i == distilled) distilled = j;
      clauses.resize (j);
    }

    // Finally delete those removed clauses.
    //
    int64_t reduced = 0;
    for (size_t i = keep; i < candidates.size (); i++)
      delete_clause (candidates[i]), reduced++;
    stats.clauses.reduced += reduced;

    candidates.release (mem);
    LOG ("reduced " << reduced << " clauses");

  }

  // The reduce limit is increased in a simple arithmetic way, which results
  // in the conflict number at which reduction occurs to form a second order
  // arithmetic series.
  //
  limits.reduce.redundant += opts.redinc;

  report ('-');
  stop ();
}

/*------------------------------------------------------------------------*/

void Solver::assume (int decision) {
  assert (propagated ());
  level++;
  int height = (int) trail.size ();
  LOG ("push frame decision " << decision << " trail " << height);
  control.push_back (mem, Frame (decision, level, height));
  assert (level + 1 == (int) control.size ());
  assign (decision);
}

// Check if there are unassigned literals left.  If all are assigned return
// false.  Otherwise return the one with maximium score in the decisions
// priority dequeue.
//
bool Solver::decide () {
  assert (propagated ());
  int decision = 0;
  while (!decision && !decisions.empty ()) {
    int lit = decisions.top ();
    decisions.pop (lit);
    if (!val (lit)) decision = lit;
  }
  if (!decision) return false;
  assert (decision > 0);
  if (phases[decision] < 0) decision = -decision;
  stats.decisions++;
  stats.levels += level;
  assume (decision);
  return true;
}

/*------------------------------------------------------------------------*/

void Solver::strengthen (Clause * c, int remove) {
  if (c->dumped || satisfied (c)) return;
  assert (addedlits.empty ());
  int lit;
  for (const int * p = c->lits; (lit = *p); p++)
    if (lit != remove && !val (lit))
      addedlits.push_back (mem, lit);
  assert (clause_entailed ());
  new_push_connect_clause ();
  addedlits.clear ();
  dump_clause (c);
}

int Solver::size_penalty () const {
  int64_t num_clauses = (stats.clauses.irredundant)/opts.sizepen;
  int logres = 0;
  while (num_clauses && logres < opts.sizemaxpen) num_clauses >>= 1, logres++; 
  return 1 << logres;
}

void Solver::distill () {

  int64_t steps = limits.simp.steps;

  assert (!level);
  assert (!opts.plain);

  if (!opts.distill) return;

  start (&stats.time.simp.distill);
  LOG ("distill");

  int64_t limit = stats.steps + steps / size_penalty ();
  bool changed = false, done = clauses.empty ();
  if (distilled >= clauses.size ()) distilled = 0;

  while (!empty && !done && stats.steps++ < limit) {

    Clause * c = clauses[distilled];

    if (!c->redundant &&		// Distill only irredundant clauses.
        !c->dumped &&
	c->large () &&			// Otherwise 'ignore' does not work.
        !satisfied (c)) {

	Clause * conflict = 0;

	assert (!ignore);
	ignore = c;	   // BCP ignoring this clause.
			   // We also have to ignore all redundant clauses
			   // during BCP (see our IJCAR'12 paper).
	int lit;
	for (const int * p = c->lits; !conflict && (lit = *p); p++) {
	  Value v = val (lit);
	  if (v == FALSE) continue;
	  if (v == TRUE) {
	    // Let 'l_i' be the literals of the clause, then having
	    // 'l_{n+1}' = 'lit' assigned to false assuming all the previous
	    // 'l_i' to be false as well, shows the following implication
	    //
	    //   -l_0 -l_1 ... -l_n -> -lit
	    //
	    // which in turn is equivalent to the following clause
	    //
	    //   l_0  l_1 ... l_n -lit
	    //
	    // Self-subsuming resolution of 'c' with this clause allows
	    // to remove 'lit' from 'c'.  This is the same 'strengthening'
	    // techniques as in distillation and vivification.
	    //
	    LOGCLAUSE (c, "distilliation removes " << lit << " in");
	    stats.distill.strengthened++;
	    strengthen (c, lit);
	    changed = true;
	    assert (!conflict);  // So simply continue with next clause.
	    break;
	  }
	  assume (-lit);
	  conflict = bcp ();
	}

	ignore = 0;	// Do not forget to disable ignoring this clause.

	if (conflict) {	// 'F \ {c} & !c' is inconsistent.

	  // We analyze the conflict in order to detect, whether we actually
	  // produced an empty clause (?) or a unit.  Futher the learned
	  // clause might be shorter than the original clause and thus
	  // actually better suited for future propagation than if we would
	  // have turned 'c' into a learned clause.

	  analyze (conflict);

	  if (empty) LOG ("distilled empty clause");	// TODO reachable?
	  else if (!level) {
	    assert (next < trail.size ());

	    // This code subsumes failed literal probing in particular in
	    // combination with first UIP clause learning in 'analyze'.
	    //
	    LOGCLAUSE (c, "distilled unit from");
	    stats.distill.units++;
	    changed = true;
	    conflict = bcp ();
	    if (conflict) analyze (conflict), assert (empty);
	  } else {

	    // In any case having produced conflict shows that 'c' is an
	    // Asymmetric Tautology (AT) and thus can be removed.   We
	    // already learned a (redundant) clause for this conflict which
	    // might be a shorter version of 'c' actually.  So we do not use
	    // the option discussed in our IJCAR'12 paper to mark 'c' to be
	    // redundant instead of simply removing (dumping) it. 

	    LOGCLAUSE (c, "distillation subsumes");
	    stats.distill.subsumed++;
	    dump_clause (c);
	  }
	}

	backtrack ();
    }

    // Next clause ...
    //
    if (++distilled == clauses.size ()) {
      if (changed) changed = false;		// Reset.
      else done = true;				// One round without change.
      distilled = 0;				// Start over.
    }
  }
  stop ();
  report ('l');
}

/*------------------------------------------------------------------------*/

// Check clause to be satisfied by current assignment.

bool Solver::satisfied (Clause * c) const {
  int lit;
  if (c->satisfied) return true;
  for (const int * p = c->lits; (lit = *p); p++)
    if (val (lit) == TRUE) {
      if (!level) c->satisfied = true;
      return true;			// here: return early
    }
  return false;
}

/*------------------------------------------------------------------------*/

// Update respectively shrink glue value if current value is smaller.

void Solver::update_glue (Clause * c) {
  if (!opts.glueupdate) return;
  if (!opts.gluered) { assert (!c->glue); return; }
  assert (frames.empty ());
  for (const int * p = c->lits; *p; p++) mark_frame (*p);
  int new_glue = unmark_frames ();
  if (new_glue >= c->glue) return;
  LOGCLAUSE (c, "update glue " << c->glue << " to " << new_glue << " of ");
  c->glue = new_glue;
  stats.glues.sum += new_glue;
  stats.glues.count++;
  stats.glues.updates++;
}

/*------------------------------------------------------------------------*/

// Try to resolve two clauses.  Returns true iff resolvent is non trivial.

bool Solver::try_resolve (Clause * c, int pivot, Clause * d) {
  assert (!c->dumped && !c->satisfied);
  assert (!d->dumped && !d->satisfied);
  bool res = true;
  assert (seen.empty ());
  int lit;
  LOGCLAUSE (c, "try resolve 1st antecedent");
  stats.steps++;
  for (const int * p = c->lits; (lit = *p); p++) {
    if (lit == pivot) continue;
    assert (!marked (lit));
    mark (lit);
  }
  LOGCLAUSE (d, "try resolve 2nd antecedent");
  stats.steps++;
  for (const int * p = d->lits; res && (lit = *p); p++) {
    if (lit == -pivot) continue;
    int m = marked (lit);
    if (m > 0) continue;
    if (m < 0) res = false;
    else mark (lit);
  }
  unmark ();
  if (res) LOG ("resolvent would not be tautological");
  else LOG ("resolvent would be tautological");
  return res;
}

/*------------------------------------------------------------------------*/

void Solver::backward (Clause * c, int ignore) {

  LOGCLAUSE (c, "backward subsuming with");

  // First get non-false literal different from 'ignore' with the smallest
  // number of occurences.
  //
  int minlit = 0;
  {
    int minoccs = INT_MAX, lit, litoccs;
    stats.steps++;
    for (const int * p = c->lits; (lit = *p); p++) {

      if (lit == ignore) continue;	// Ignore this literal since
					// we want to remove it's negation
					// through self-subsuming
					// resolution.
      if (val (lit) < 0) continue;
      litoccs = occs (lit).count;
      if (minlit && minoccs >= litoccs) continue;
      minlit = lit;
      minoccs = litoccs;
    }

    // Skip backward check if 'minlit' has too many occurrences.
    //
    if (minoccs >= opts.bwocclim) return;
  }

  assert (minlit);
  LOG ("minimum occurrence literal " << minlit << " ignoring " << ignore);

  // Now mark all literals in 'c'.
  //
  for (const int * p = c->lits; *p; p++) mark (*p);

  // Finally go through all clauses with 'minlit' (different from 'c') and
  // find those that are subsumed by 'c' or have all the literals of 'c'
  // except for one which occurs 'negated' in 'd'.  The former will allow to
  // remove 'd', and the latter to strengthen 'd' by removing the negated
  // literal through self-subsuming resolution.
  {
    Occs & os = occs (minlit);
    vector<Clause *>::const_iterator i;
    for (i = os.begin (); i != os.end (); i++) {
      Clause * d = *i;
      if (d == c) continue;
      int lit, count = seen.size (), negated = 0;
      stats.steps++;
      for (const int  * p = d->lits; count && (lit = *p); p++) {
	int m = marked (lit);
	if (!m) continue;
	assert (count > 0);
	count--;
	if (m > 0) continue;
	assert (m < 0);
	if (negated) {
	  count = INT_MAX;		// Can not (self-)subsume at
	  break;			// least two negated literals.
	}
	negated = lit;			// Remember negated.
      }

      if (count) continue;		// Not enough overlap, or at least
					// two overlapping negated literals.
      if (negated) {
	LOGCLAUSE (d, "backward remove " << negated << " in");
	LOGCLAUSE (c, "self-subsuming");

	// There are essentially two alternatives to remove 'negated' from
	// 'd'.  The first one is to remove 'negated' from 'd->lits'.  This
	// would require to remove the occurrence of 'negated' in 'd' from
	// the occurrence list of 'negated'.  Further, it makes it difficult
	// to properly handle the memory allocated by 'd', since 'd->lits'
	// kind of shrinks (this is doable by spending another bit in the
	// header of clauses and adding say INT_MAX markers at the end of
	// 'd->lits').  Removing the occurrence could also be costly.
	//
	// The second alternative is to dump 'd' and add a new clause with
	// the same literals as 'd'.  However, since in 'backward (int)' we
	// have an iterator over the occurence list/vector of 'ignore' which
	// occurs in 'd', we would need to make sure to fix this iterator if
	// a new clause occurrence is added to this list/vector.  This
	// 'fixing' could be done automatically or explicitly, which in both
	// cases is awkward.
	//
	// Instead we simply remember the clause literal pairs
	// ('d','negated'), where 'negated' can be removed from 'd' and save
	// those pairs on the 'tostrengthen' stack.  The actual removal,
	// e.g. flushing 'tostrengthen', is postponed (currently until all
	// the clauses with 'lit' are handled).
	//
	tostrengthen.push_back (mem, ClsPtrLitPair (d, negated));
	stats.backward.strengthened++;
      } else {
	LOGCLAUSE (d, "subsumed");
	LOGCLAUSE (c, "subsuming");
	stats.backward.subsumed++;
	dump_clause (d);		// Could be moved, i.e. marked as
					// redundant instead.
      }
    }
  }

  unmark ();
}

void Solver::backward (int lit) {
  assert (!level);
  assert (dense);
  assert (tostrengthen.empty ());
  Occs & os = occs (lit);
  vector<Clause *>::const_iterator i;
  for (i = os.begin (); i != os.end (); i++) {
    Clause * c = *i;
    assert (!c->redundant);
    stats.steps++;
    if (c->dumped) continue;
    if (c->size () >= opts.bwclslim) continue;	// too large?
    if (satisfied (c)) continue;
    backward (c, lit);
  }
  while (!tostrengthen.empty ()) {
    ClsPtrLitPair cplp = tostrengthen.back ();
    tostrengthen.pop_back ();
    strengthen (cplp.cls, cplp.lit);
  }
}

/*------------------------------------------------------------------------*/

// Determine whether resolution of all positive with all negative
// occurrences of 'cand' produces at most as many non trivial resolvents as
// the number of positive plus negative occurrences.

bool Solver::try_elim (int cand) {
  assert (var (cand).free ());
  Occs & p = occs (cand), & n = occs (-cand);
  int64_t limit = p.count + n.count;
  vector<Clause *>::const_iterator i;
  for (i = p.begin (); limit >= 0 && i != p.end (); i++) {
    Clause * c = *i;
    assert (!c->redundant);
    stats.steps++;
    if (c->dumped || satisfied (c)) continue;
    vector<Clause *>::const_iterator j;
    for (j = n.begin (); limit >= 0 && j != n.end (); j++) {
      Clause * d = *j;
      assert (!d->redundant);
      stats.steps++;
      if (d->dumped || satisfied (d)) continue;
      if (try_resolve(c, cand, d)) limit--;
    }
  }
  if (limit < 0) LOG ("elim " << cand << " failed");
  else LOG ("elim " << cand << " removes " << limit << " clauses");
  return limit >= 0;
}

/*------------------------------------------------------------------------*/

// Now really do resolve those clauses and add resolvent.

void Solver::do_resolve (Clause * c, int pivot, Clause * d) {
  assert (!c->dumped && !c->satisfied);
  assert (!d->dumped && !d->satisfied);
  assert (addedlits.empty ());
  int lit;
  LOGCLAUSE (c, "do resolve 1st antecedent");
  stats.steps++;
  for (const int * p = c->lits; (lit = *p); p++)
    if (lit != pivot) addedlits.push_back (mem, lit);
  LOGCLAUSE (d, "do resolve 2nd antecedent");
  stats.steps++;
  for (const int * p = d->lits; (lit = *p); p++)
    if (lit != -pivot) addedlits.push_back (mem, lit);
  if (trivial_clause ()) LOG ("resolvent tautological");
  else new_push_connect_clause ();
  addedlits.clear ();
}

/*------------------------------------------------------------------------*/

void Solver::push_extension (int lit) {
  extension.push_back (mem, lit);
  LOG ("push extension " << lit);
}

void Solver::push_extension (Clause * c, int blit) {
  push_extension (0);
  int lit;
  for (const int * p = c->lits; (lit = *p); p++)
    if (lit != blit) push_extension (lit);
  push_extension (blit);
}

// Really do eliminate 'cand' by variable elimination.

void Solver::do_elim (int cand) {

  assert (schedule);
  assert (var (cand).free ());

  // Add all resolvents on 'cand'.
  //
  Occs & p = occs (cand), & n = occs (-cand);
  vector<Clause *>::const_iterator i;
  for (i = p.begin (); i != p.end (); i++) {
    Clause * c = *i;
    stats.steps++;
    if (c->dumped || satisfied (c)) continue;
    vector<Clause *>::const_iterator j;
    for (j = n.begin (); j != n.end (); j++) {
      Clause * d = *j;
      stats.steps++;
      if (d->dumped || satisfied (d)) continue;
      do_resolve (c, cand, d);
    }
  }

  // Save clauses with smaller occurences on extension stack.
  {
    int extend;
    Occs * e;
    if (p.count < n.count) extend = cand, e = &p;
    else extend = -cand, e = &n;
    for (i = e->clauses.begin (); i != e->clauses.end (); i++) {
      Clause * c = *i;
      if (c->dumped || satisfied (c)) continue;
      stats.steps++;
      push_extension (c, extend);
    }
    push_extension (0);
    push_extension (-extend);
  }

  // Dump clauses with positive occurrences.
  //
  while (!p.clauses.empty ()) {
    Clause * c = p.clauses.back ();
    stats.steps++;
    p.clauses.pop_back ();
    if (c->satisfied || c->dumped) continue;
    stats.clauses.eliminated++;
    dump_clause (c);
  }
  p.clauses.release (mem);

  // Dump clauses with negative occurrences.
  //
  while (!n.clauses.empty ()) {
    Clause * c = n.clauses.back ();
    stats.steps++;
    n.clauses.pop_back ();
    if (c->satisfied || c->dumped) continue;
    stats.clauses.eliminated++;
    dump_clause (c);
  }
  n.clauses.release (mem);

  LOG ("eliminated " << cand);
  var (cand).state = Var::ELIMINATED;
  stats.vars.eliminated++;

  // Variable elimination might occasionally produce units, so we flush the
  // BCP queue by default.
  //
  Clause * conflict = bcp ();
  if (conflict) analyze (conflict), assert (empty);

  touch_fixed ();		// Finally touch newly satisfied clauses.
}

/*------------------------------------------------------------------------*/

bool Solver::contains_eliminated (Clause * c) const {
  int lit;
  for (const int * p = c->lits; (lit = *p); p++)
    if (var (lit).state == Var::ELIMINATED) return true;
  return false;
}

void Solver::dump_eliminated_redundant () {
  vector<Clause*>::const_iterator i;
  for (i = clauses.begin (); i != clauses.end (); i++) {
    Clause * c = *i;
    if (!c->redundant || c->satisfied || c->dumped) continue;
    if (contains_eliminated (c)) dump_clause (c);
  }
}

priority_queue<int64_t> & Solver::current_cands (void) {
  switch (simplifier) {
    case BLOCK:
      return cands.block;
    default:
      assert (simplifier == ELIM);
      return cands.elim;
  }
}

/*------------------------------------------------------------------------*/

// Update candidates priority queue for new rounds of simplification.

void Solver::update_cands (void) {

  // Connect all clauses through full occurrence lists.
  //
  if (!dense) connect_occs ();

  // The 'schedule' flag, if set to 'true', makes sure that literals in
  // subsequently added and removed clauses are touched.
  //
  assert (!schedule);
  schedule = true;

  // If there are still variables not tried last time, do not touch all
  // variables, but only use these previously untried as elimination
  // candidates.  Otherwise, e.g. if 'cands.elim.empty ()' returns true,
  // touch all.
  //
  if ((schedule = current_cands ().empty ()))
    for (int idx = 1; idx <= maxvar (); idx++) touch (idx);

  // Further, touch literals in satisfied clauses.
  //
  schedule = true;
  touch_fixed ();
}

/*------------------------------------------------------------------------*/

// (Bounded) Variable Elimination = (B)VE = Clause Distribution = DP 
//
// This is the same bounded variable elimination technique as in the SatElite
// paper except that it performs subsumption and self-subsuming resolution
// on-the-fly while trying to eliminate a candidate variable.  It is also
// limited by the number of resolution steps.
//
// On-the-fly subsumption uses the observation, that, except for scheduling,
// subsumption somehow only helps variable elimination.  Thus it could also
// be done only for those literals we are currently trying to eliminate.
// Subsumption and self-subsumption is performed among the clauses in which
// the elimination candidate variable occurs before trying to do any (other)
// resolution.  Since we try to eliminate all literals eventually, we will
// also perform all the (self-)subsumption checks eventually.

bool Solver::donotelim (int cand) {
  int sign;

  // One sided limit on the number of pivot occurrences, e.g. do not
  // eliminate variables which occurs more than 'elmpocclim1' often
  // either positively or negatively.
  //
  if (occs (cand).count > opts.elmpocclim1) return true;
  if (occs (-cand).count > opts.elmpocclim1) return true;

  // Two-sided limit, e.g. do not eliminate variables which occurs
  // more than 'elmpocclim2' both positively and negatively.
  //
  if (occs (cand).count > opts.elmpocclim2 &&
      occs (-cand).count > opts.elmpocclim2) return true;

  for (sign = -1; sign <= 1; sign += 2) {
    Occs & os = occs (sign * cand);
    vector<Clause *>::const_iterator i;
    for (i = os.begin (); i != os.end (); i++) {
      Clause * c = *i;
      assert (!c->redundant);
      if (c->size () >= opts.elmclslim) return true;
      for (const int * p = c->lits; *p; p++)
	if (occs (*p).count >= opts.elmocclim) return true;
    }
  }
  return false;
}

void Solver::elim () {

  int64_t steps = limits.simp.steps;

  assert (!level);
  assert (!opts.plain);

  if (!opts.elim) return;

  start (&stats.time.simp.elim);

  assert (simplifier == TOPSIMP), simplifier = ELIM;
  update_cands ();

  int64_t limit;

  if (stats.simplifications <= opts.elmrtc) limit = LLONG_MAX;
  else {
    // Scale the allowed number of steps by the size of the formula.
    limit = stats.steps + 10*steps/size_penalty ();
  }

  while (!empty && !cands.elim.empty () && stats.steps++ < limit) {
    int cand = cands.elim.top ();
    int64_t priority = cands.elim.priority (cand);
    cands.elim.pop (cand);
    if (!priority || !var (cand).free () || donotelim (cand)) continue;
    LOG ("next elim cand " << cand << " priority " << priority);
    backward (cand);
    backward (-cand);
    if (try_elim (cand)) do_elim (cand);
  }
  assert (schedule), schedule = false;
  assert (simplifier == ELIM), simplifier = TOPSIMP;

  // Redundant clauses with eliminated variables are dumped.
  //
  dump_eliminated_redundant ();

  report ('e');

  stop ();
}

/*------------------------------------------------------------------------*/

// Check whether all resolvents with 'c' on 'blit' are tautological.

bool Solver::block_clause (Clause * c, int blit) {
  if (c->dumped || satisfied (c)) return false;
  Occs & os = occs (-blit);
  vector<Clause *>::const_iterator i;
  for (i = os.begin (); i != os.end (); i++) {
    Clause * d = *i;
    assert (!d->redundant);
    stats.steps++;
    if (d->dumped || satisfied (d)) continue;
    if (try_resolve (c, blit, d)) return false;
  }
  return true;
}

// Find and remove all blocked clauses blocked on 'blit'.

void Solver::block_lit (int blit) {
  Occs & os = occs (blit);
  vector<Clause *>::const_iterator i;
  for (i = os.begin (); i != os.end (); i++) {
    Clause * c = *i;
    assert (!c->redundant);
    if (!block_clause (c, blit)) continue;
    LOGCLAUSE (c, "blocked on " << blit);
    stats.clauses.blocked++;
    push_extension (c, blit);
    dump_clause (c);
  }
}

// Blocked Clause Elimination (see our TACAS'10 paper and JAR'1X article).
//
void Solver::block () {

  int64_t steps = limits.simp.steps;

  assert (!level);
  assert (!opts.plain);

  if (!opts.block) return;
  if (opts.blkwait >= stats.simplifications) return;

  start (&stats.time.simp.block);

  assert (simplifier == TOPSIMP), simplifier = BLOCK;
  update_cands ();

  int64_t limit;
  if (stats.simplifications <= opts.blkrtc) limit = LLONG_MAX;
  else limit = stats.steps + 10*steps/size_penalty ();

  while (!empty && !cands.block.empty () && stats.steps++ < limit) {
    int cand = cands.block.top ();
    int64_t priority = cands.block.priority (cand);
    cands.block.pop (cand);
    if (!priority || !var (cand).free ()) continue;
    LOG ("next block cand " << cand << " priority " << priority);
    block_lit (cand);
    block_lit (-cand);
  }
  assert (schedule), schedule = false;
  assert (simplifier == BLOCK), simplifier = TOPSIMP;

  report ('b');

  stop ();
}

/*------------------------------------------------------------------------*/

void Solver::disconnect_clauses () {

  // Important clauses can only be deleted during 'collect' and particulary
  // not during 'reduce' (since they are kept).  In 'collect' the clauses
  // are first disconnected with this function and then later connected back
  // with 'connect_clause' through 'connect_clauses' (plural!).  In
  // 'connect_clause' (singular) this counter is increased, thus that at the
  // end of 'collect' we have again the correct number of important clauses.
  //
  limits.reduce.important = 0;

  for (int idx = 1; idx <= maxvar (); idx++)
    for (int sign = -1; sign <= 1; sign += 2) {
      int lit = sign*idx;
      watches (lit).release (mem);
      occs (lit).release (mem);
    }
  dense = false;
}

// Remove false literals from clause which has the main benefit of shorting
// clauses for faster 'bcp', in particular if the clause becomes binary.
// Actually if the clause can not be reduced we return it as is other wise
// we delete the clause and return a new reduced shorter copy.
//
Clause * Solver::reduce_clause (Clause * c) {
  int lit;
  for (const int * i = c->lits; (lit = *i); i++)
    if (val (lit) < 0) break;
  if (!lit) return c;		// Without false literal, nothing to reduce.
  assert (addedlits.empty ());
  for (const int * i = c->lits; (lit = *i); i++)	// Copy non-false
    if (val (lit) >= 0) addedlits.push_back (mem, lit);	// literals only.
  bool redundant = c->redundant;
  int glue = c->glue;
  delete_clause (c);
  Clause * res = new_clause (redundant, glue);
  addedlits.clear ();
  return res;
}

// Collect and delete all top-level satisfied or dumped non-forcing clauses.
//
void Solver::collect_clauses () {
  assert (!level);			    // We are at the top-level.
  size_t i, j = 0;
  int64_t collected = 0;
  for (i = j; i < clauses.size (); i++) {
    if (i == distilled) distilled = j;      // Fix distilled.
    Clause * c = clauses[i];
    if (c->forcing ||			    // Keep forcing and neither
        !(c->dumped || satisfied (c)))	    // dumped nor satisfied clauses.
      clauses[j++] = reduce_clause (c);     // Remove false literals.
    else delete_clause (c), collected++;
  }
  if (i == distilled) distilled = j;	    // Fix distilled.
  clauses.resize (j);
  clauses.shrink (mem);
  stats.clauses.collected += collected;
  LOG ("collected " << collected << " clauses");
}

void Solver::connect_clauses () {
  assert (!dense);
  vector<Clause*>::const_iterator i;
  for (i = clauses.begin (); i != clauses.end (); i++)
    connect_clause (*i);
}

// This is the actual garbage collector.  It simply flushes all the
// occurrences and watches first, and then after deleting and reducing
// clauses, connects them again.  By disconnecting and connecting seperately
// we avoid the quadratic overhead in removing watches and occurrences
// during the actual collecting phase.  In addition, we can reuse the
// 'connect_clause' clause code and thus easier make sure we continue to
// properly watch clauses.
//
void Solver::collect () {
  start (&stats.time.simp.collect);
  LOG ("collect");
  disconnect_clauses ();
  collect_clauses ();
  connect_clauses ();
  report ('c');
  stop ();
}

/*------------------------------------------------------------------------*/

// Inprocessing simplifies the instance interleaved with CDCL search.
// The steps limit given as arguments restricts the amount of work allowed
// in each inprocessing interval.
//
Status Solver::simplify () {

  start (&stats.time.simp.total);

  LOG ("simplify");

  // Go into simplification mode, which disables phase saving, and
  // computes different statistics.
  //
  assert (!simplifier);
  simplifier = TOPSIMP;

  stats.simplifications++;
  backtrack ();			// Simplification restarts implicitly.

  int vars_before = remaining_vars ();

  if (opts.plain) goto COLLECT;

  if (!empty) distill ();	// Distillation.
  if (!empty) block ();		// Blocked clause elimination.
  if (!empty) elim ();		// Bounded variable elimination.

  // TODO add more inprocessors of course ...

COLLECT:

  collect ();			// Collect garbage clauses.

  // Now go back to search mode.
  //
  assert (simplifier == TOPSIMP);
  simplifier = NOSIMP;
  assert (!simplifier);

  int vars_after = remaining_vars ();
  assert (vars_before >= vars_after);
  limits.simp.removed_vars = vars_before - vars_after;

  stop ();

  return empty ? UNSATISFIABLE : UNKNOWN;
}

// This is the CDCL search loop.  It is limited by the number of conflicts,
// since we interleave it with simplification through inprocessing.
//
Status Solver::search () {
  start (&stats.time.search.total);
  LOG ("search");
  int64_t conflicts = 0;
  Clause * conflict;
  Status res = UNKNOWN;
  while (!res)
    if (empty) res = UNSATISFIABLE;
    else if ((conflict = bcp ())) analyze (conflict), conflicts++;
    else if (conflicts >= limits.search.conflicts) break;
    else if (reducing ()) reduce ();
    else if (restarting ()) restart ();
    else if (!decide ()) res = SATISFIABLE;
  stop ();
  return (Status) res;
}

/*------------------------------------------------------------------------*/

void Solver::print () const {
  cout << "p cnf " << maxvar () << ' ' << clauses.size () << endl;
  vector<Clause*>::const_iterator i;
  for (i = clauses.begin (); i != clauses.end (); i++)
    (*i)->print ();
}

/*------------------------------------------------------------------------*/

// Extend current partial assignment to full assignment using the clauses
// saved on the extension stack.  At the same time flip assigned literals
// which would falsify saved clauses.  The extension stack is filled in such
// a way that this will always produce a correct satisfying assignment for
// the original clauses.  This techniques was used first in MiniSAT by
// Niklas Soerensson and is also described in our SAT'10 resp. IJCAR'12
// papers.
//
void Solver::extend () {
  start (&stats.time.search.extend);
  LOG ("extend");
  while (!extension.empty ()) {
    int lit = extension.back (), other;
    bool satisfied = false;
    while ((other = extension.back ())) {
      extension.pop_back ();
      if (val (other) == TRUE) satisfied = true;
    }
    extension.pop_back ();
    if (!satisfied) {
      LOG ("flipping assignment of " << lit);
      vals[abs (lit)] = sign (lit);
    }
  }
  stop ();
}

// Return TRUE if one of the literals traversed is TRUE.  Continue
// traversing and advancing the iterator until a zero sentinel is found.

bool Solver::satisfied (vector<int>::const_iterator & it) const {
  bool res;
  int lit;
  for (res = false; (lit = *it); it++)
    if (val (lit) == TRUE) res = true;	// do not return!
  return res;
}

void Solver::check_satisfied () const {
  vector<int>::const_iterator c;
  for (c = original.begin (); c != original.end (); c++)
    if (!satisfied (c))
      internal_error ("invalid assignment");
}

/*------------------------------------------------------------------------*/
/**** PART6: implemenetation of public member functions *******************/
/*------------------------------------------------------------------------*/

Solver::Solver () : 
  level (0),			// top-level = root-level = 0
  dense (false),		// start in sparse mode only watching clauses
  simplifier (NOSIMP),		// only non zero within 'simplify'
  schedule (false),		// do not schedule touched literals
  next (0),			// nothing to propagate yet
  touched (0),			// no fixed literal touched yet
  empty (0),			// no empty clause found yet
  score_increment (1),		// initial bump score increment
  ignore (0),			// use all clauses during bcp
  distilled (0)			// all before clauses[distilled] distilled
{ 
  if (sizeof (int64_t) != 8) internal_error ("sizeof (int64_t) != 8");

  set_fpu_to_double_precision ();
  mem.inc (sizeof *this);

  // We need an empty frame for the top-level without decision.
  //
  control.push_back (mem, Frame ());

#ifndef NLGL
  lgl = lglinit ();
#endif
}

/*------------------------------------------------------------------------*/

Solver::~Solver () {

  // Beside memory for the vectors in containers declared as members of
  // 'Solver' the only other memory we allocate is for clauses and for the
  // watch and occurrence list containers of variables.

  while (!clauses.empty ()) {
    Clause * c = clauses.back ();
    clauses.pop_back ();
    delete_clause (c);
  }
  while (!vars.empty ()) {	// 'vector' does not call destructors
    Var & v = vars.back ();
    v.release (mem);		// so need to release variables explicitly
    vars.pop_back ();
  }

  vars.release (mem);
  vals.release (mem);
  phases.release (mem);
  control.release (mem);
  trail.release (mem);
  decisions.release (mem);
  original.release (mem);
  addedlits.release (mem);
  clauses.release (mem);
  seen.release (mem);
  cands.elim.release (mem);
  cands.block.release (mem);
  tostrengthen.release (mem);
  frames.release (mem);
  probes.release (mem);
  extension.release (mem);
  timers.release (mem);

#ifndef NLGL
  lglrelease (lgl);
#endif

  mem.dec (sizeof *this);

  // Check that all heap memory has been deallocated to avoid memory leaks.
  // To skip the 'abort' while using tools like 'valgrind' to actual debug
  // the memory leak, if this assertion fails, define the environment
  // variable 'CLEANELINGLEAK', e.g. if this assertion fails on 'test.in',
  // then use the following:
  /*
CLEANELINGLEAK=yes \
valgrind --show-reachable=yes --leak-check=yes --show-reachable \
./cleaneling test.in
  */
  assert (getenv ("CLEANELINGLEAK") || !mem.current);
}

/*------------------------------------------------------------------------*/

void Solver::addlit (int lit) {
  original.push_back (mem, lit);
  if (lit) {
    import (lit);
    addedlits.push_back (mem, lit);
  } else {
    if (!trivial_clause ()) new_push_connect_clause ();
    addedlits.clear ();
  }
}

/*------------------------------------------------------------------------*/
// Both simplification (inprocessing) and search (CDCL) are limited.  The
// first by a virtual metric of 'steps' and the second by the number of
// conflicts.  If the limit is hit the solver switches to the other phase.
// Fist we initialize limits.
//
void Solver::init_limits () {

  new_restart_limit ();

  if (!limits.simp.steps) limits.simp.steps = opts.simpint;

  // It is a good idea, at least for running the SAT solver stand-alone
  // in a kind of 'competition mode' to boost the initial effort
  // for simplification.
  //
  if (!stats.simplifications) {
    assert (opts.boost > 0);
    if (limits.simp.steps >= INT_MAX / opts.boost)
      limits.simp.steps = INT_MAX;
    else limits.simp.steps *= opts.boost;
  }

  if (!limits.search.conflicts && opts.searchfirst)
    limits.search.conflicts = opts.searchint;

  // This has to be initialized in case 'search' is run first.
  //
  limits.simp.removed_vars = 0;
}

// After one round of search and simplification we update both the search
// and simplification limits to allow more time for the next round.
//
void Solver::update_limits () {

  if (opts.simpgeom) {
    // Use a geometric increase of simplification effort as well.
    //
    if (limits.simp.inc >= INT_MAX/2) limits.simp.inc = INT_MAX;
    else if (!limits.simp.inc) limits.simp.inc = opts.simpint;
    else limits.simp.inc *= 2;
  } else {
    // In thise case we use an arithmetic increase for the simplification
    // limit, which means that over time less and less relative effort is
    // spent in simplification, since the search limit is increased
    // geometrically.
    //
    if (limits.simp.inc >= INT_MAX - opts.simpint) limits.simp.inc = INT_MAX;
    else limits.simp.inc += opts.simpint;
  }

  limits.simp.steps = limits.simp.inc;

  // Optionally bound the maximum number of steps in simplification.
  //
  if (opts.stepslim && limits.simp.steps > opts.stepslim)
    limits.simp.steps = opts.stepslim;

  if (!limits.search.inc) limits.search.inc = opts.searchint;
  assert (limits.search.inc);

  // The conflict interval allocated for search increases exponentially.
  //
  if (limits.search.conflicts) {

    int inc = limits.search.inc;// Default is to double the increment.

    // As metric to measure the efficiency of simplification we use
    // the number of removed variables during the last simplification.
    // This in turn is used to reduce the next search interval, which
    // otherwise would grow in geometric fashion.  Thus is simplification
    // is efficient spent more effort in simplification.
    //
    // See also the comments for 'limits.simp.removed_vars'.
    //
    // Note that the only dynamic intensification of 'search' is due
    // to learned units during an 'iteration'.  See the code for reporting
    // and counting iterations in 'analyze'.
    //
    int64_t removed = limits.simp.removed_vars;

    if (removed > 0 && remaining_vars ()) {
	int64_t reduction = (100*removed) / remaining_vars ();
	if (reduction > 1) {
	  // Reduce the next conflict interval for search by the percentage
	  // of reduction in terms of removed variables during the last
	  // simplification.  For instance if 10% variables are removed
	  // the interval is divided by 10, removing 2% variables reduces
	  // the interval by factor of 2 etc.
	  //
	  inc /= reduction;
	}
      }

    if (limits.search.inc >= INT_MAX - inc) limits.search.inc = INT_MAX;
    else limits.search.inc += inc;

  } else assert (!opts.searchfirst);	// Do not update in this case.

  limits.search.conflicts = limits.search.inc;
}

/*------------------------------------------------------------------------*/

Status Solver::solve () {

  if (!addedlits.empty ()) internal_error ("API usage: zero literal missing");

  start (&stats.time.solve);
  LOG ("solve");

  init_limits ();
  bias_phases ();

  Status res;

  //------------------------------------------------------------------------
  // MAIN SOLVING LOOP
  //
  // This is the main solving loop, which consists of limited interleaved
  // CDCL search and simplifications through inprocessing.
  //
  for (;;)
    if ((res = search ())) break;
    else if ((res = simplify ())) break;
    else update_limits ();
  //------------------------------------------------------------------------

  char type;				// for reporting result
  switch (res) {
    default:
      type = '?';
      break;
    case SATISFIABLE:
      type = '1'; 
      extend ();			// extend partial assignemtn to full
      check_satisfied ();		// check assignment
      break;
    case UNSATISFIABLE: 
      type = '0';
#ifndef NLGL
      if (opts.check) assert (lglsat (lgl) == 20);  // check consistency
#endif
      break;
  }

  report (type);
  if (opts.verbose && stats.reported) print_line ();
  stop ();

  return res;
}

/*------------------------------------------------------------------------*/

Value Solver::value (int lit) { import (lit); return val (lit); }

/*------------------------------------------------------------------------*/

int Solver::option (const char * str) {		// get current option value
#undef OPTION
#define OPTION(NAME,DEFAULT,MIN,MAX,DESCRIPTION) \
  if (!strcmp (# NAME, str)) return opts.NAME.val;
#include "options.h"
  return 0;
}

void Solver::option (const char * str, int val) {	// set option value
#undef OPTION
#define OPTION(NAME,DEFAULT,MIN,MAX,DESCRIPTION) \
  if (!strcmp (# NAME, str)) { \
    assert (opts.NAME.min <= opts.NAME.max); \
    if (val > opts.NAME.max) val = opts.NAME.max; \
    if (val < opts.NAME.min) val = opts.NAME.min; \
    opts.NAME.val = val; \
    return; \
  }
#include "options.h"
}

/*------------------------------------------------------------------------*/

// Print the time of one solving phase (or subphase) relative to the total
// time used for solving.  We keep it this part even though it is internal.

void Solver::print_time_for (double time, const char * name) {
  cout 
  << "c "
  << fixed << setprecision (2) << setw (8) << right << time
  << ' ' << setw (10) << left << name
  << ' ' << fixed << setprecision (0) << right << setw (3) 
  << percent (time,stats.time.solve) << '%'
  << endl;
}

#define PRINT_TIME_FOR(FIELD,NAME) \
  print_time_for (stats.time.FIELD, NAME)

void Solver::statistics () {
  print_line ();//--------------------------------------------------------
  cout << "c sizeof (Var) = " << sizeof (Var) << endl;
  cout << "c sizeof (Occs) = " << sizeof (Occs) << endl;
  cout << "c sizeof (Watch) = " << sizeof (Watch) << endl;
  cout << "c sizeof (Watches) = " << sizeof (Watches) << endl;
  cout << "c sizeof (Clause) = " << Clause::bytes (0)
       << " + size * " << sizeof (int) << endl;
  print_line ();//--------------------------------------------------------
  while (!timers.empty ()) stop ();
  double solve = stats.time.solve;
  cout
   << fixed << setprecision (0)
  << "c " 
  << stats.simplifications << " simplifications, "
  << stats.iterations << " iterations, "
  << stats.reductions << " reductions"
  << endl
  << "c " 
  << stats.restarts.count << " restarts, " 
  << stats.restarts.reuse.count << " reused trails "
  << percent (stats.restarts.reuse.count, stats.restarts.count) << "%, "
  << average (stats.restarts.reuse.sum, stats.restarts.reuse.count)
  << "% reused levels"
  << endl
  << "c " 
  << stats.decisions << " decisions, "
  << average (stats.decisions, solve) << " decisions per second"
  << endl
  << "c "
  << stats.conflicts << " conflicts, " 
  << average (stats.conflicts, solve) << " conflicts per second"
  << endl
  << "c " 
  << stats.lits.learned << " learned literals, "
  << stats.lits.minimized << " minimized literals " 
  << percent (stats.lits.minimized, stats.lits.learned) << "%"
  << endl
  << "c "
  << setprecision (2) << average (stats.glues.sum, stats.glues.count)
  << " average glue, " << stats.glues.updates << " updates"
  << endl
  << "c "
  << stats.propagations << " propagations, " << setprecision (2)
  << average (stats.propagations/1e6, solve) 
  << " million propagations per second"
  << endl;
  print_line ();//--------------------------------------------------------
  cout
  << "c "
  << stats.distill.units << " distilled units, "
  << stats.distill.subsumed << " subsumed, "
  << stats.distill.strengthened << " strengthened"
  << endl
  << "c "
  << stats.backward.subsumed << " backward subsumed clauses, "
  << stats.backward.strengthened << " strengthened"
  << endl
  << "c "
  << stats.vars.eliminated << " eliminated variables, "
  << stats.clauses.eliminated << " eliminated clauses"
  << endl
  << "c "
  << stats.clauses.blocked << " blocked clauses"
  << endl
  ;
  print_line ('=');//=====================================================
  PRINT_TIME_FOR (search.analyze, "analyze");
  PRINT_TIME_FOR (search.bias, "bias");
  PRINT_TIME_FOR (search.extend, "extend");
  PRINT_TIME_FOR (search.minimize, "minimize");
  PRINT_TIME_FOR (search.reduce, "reduce");
  PRINT_TIME_FOR (search.restart, "restart");
  print_line ();//--------------------------------------------------------
  PRINT_TIME_FOR (search.total, "search");
  print_line ('=');
  PRINT_TIME_FOR (simp.distill, "block");
  PRINT_TIME_FOR (simp.collect, "collect");
  PRINT_TIME_FOR (simp.distill, "distill");
  PRINT_TIME_FOR (simp.elim, "elim");
  print_line ();//--------------------------------------------------------
  PRINT_TIME_FOR (simp.total, "simp");
  print_line ('#');//#####################################################
  PRINT_TIME_FOR (solve, "solve");
  cout << flush;
}

size_t Solver::maximum_bytes () { return mem.maximum; }

size_t Solver::current_bytes () { return mem.current; }

/*------------------------------------------------------------------------*/

// In order to hide solver details (and potentially allow multiple
// implementations) we use an 'AbstractSolver' interface.  This also allows
// to add all the code in this one file and leave it out of the header file.
//
AbstractSolver * Cleaneling::new_solver () { return new Solver (); }
