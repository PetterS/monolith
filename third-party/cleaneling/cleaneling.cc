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

This file contains the code for the 'cleaneling' application, which is a
simple front-end for the Cleaneling library.  It consists of a parser for
(zipped) DIMACS files, and provides a command line interface for the options
available in the library.  This is achieved by sharing the 'options.h' file.
This option list is pretty printed with the '-h' option.

In addition this application sets signal handlers to catch signals and print
statistics in verbose mode even if the application terminates due to an
external signal.

--------------------------------------------------------------------------*/

#include "cleaneling.h"

/*------------------------------------------------------------------------*/

#include <cassert>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <csignal>

/*------------------------------------------------------------------------*/

using namespace std;
using namespace Cleaneling;

/*------------------------------------------------------------------------*/
namespace Cleaneling { 
/*------------------------------------------------------------------------*/

void banner (); 

/*------------------------------------------------------------------------*/

// Having a 'main' in this file means we can use global/static variables
// without harm.  They will only be used from the global 'main' function.
// This also allows us to set up signal catching to provide statistics in
// verbose mode in case the application is termianted prematurely.

static int lineno, nvars, nclauses, verbose;
static AbstractSolver * solver;
static const char * fname;
static FILE * file;

static void statistics () {
  solver->statistics ();
  print_line ();
  cout << fixed
  << "c " << setprecision (2) 
  << seconds ()
  << " seconds total process time" << endl
  << "c " << setprecision (1)
  << solver->maximum_bytes ()/(double)(1<<20)
  << " MB maximally allocated" << endl
  << flush;
}

/*------------------------------------------------------------------------*/

static int catchedsig;

static void (*sig_int_handler)(int);
static void (*sig_segv_handler)(int);
static void (*sig_abrt_handler)(int);
static void (*sig_term_handler)(int);

static void resetsighandlers (void) {
  (void) signal (SIGINT, sig_int_handler);
  (void) signal (SIGSEGV, sig_segv_handler);
  (void) signal (SIGABRT, sig_abrt_handler);
  (void) signal (SIGTERM, sig_term_handler);
}

static void catchsig (int sig) {
  if (!catchedsig) {
    cout << "s UNKNOWN" << endl << flush;
    catchedsig = 1;
    if (verbose) {
      statistics ();
      cout << "c CAUGHT SIGNAL " << sig << endl << flush;
    }
  }
  resetsighandlers ();
  raise (sig);
}

static void setsighandlers (void) {
  sig_int_handler = signal (SIGINT, catchsig);
  sig_segv_handler = signal (SIGSEGV, catchsig);
  sig_abrt_handler = signal (SIGABRT, catchsig);
  sig_term_handler = signal (SIGTERM, catchsig);
}

/*------------------------------------------------------------------------*/

static void usage () {
  const int W = 26;
  cout 
  << "usage: cleaneling [<option> ...] [<dimacs>]" << endl
  << endl
  << "where <option> is one of the following" << endl
  << endl
  << setw(W) << left << "-h | --help" << "print usage" << endl
  << setw(W) << left << "-e" << "print options" << endl
  << endl
#ifndef NLOG
  << setw(W) << left << "-l" << "shortcut for '--log=1'" << endl
#endif
  << setw(W) << left << "-v" << "shortcut for '--verbose=1'" << endl
  << setw(W) << left << "-n" << "shortcut for '--witness=0'" << endl
#ifdef NLGL
  << setw(W) << left << "-c" << "shortcut for '--check=1" << endl
#endif
  << endl
#undef OPTION
#define OPTION(NAME,DEFAULT,MIN,MAX,DESCRIPTION) \
  << setw(W) << left \
  << "--" # NAME "=" # MIN ".." #MAX \
  << DESCRIPTION << " [" << DEFAULT << "]" << endl
#include "options.h"
  << endl
  << "and <dimacs> a CNF in DIMACS format (optionally compressed)."
  << endl << flush;
}

static void options (bool all = true) {
#undef OPTION
#define OPTION(NAME,DEFAULT,MIN,MAX,DESCRIPTION) \
do  { \
  if (!all && \
      (!strcmp (#NAME, "witness") || \
       !strcmp (#NAME, "verbose"))) break; \
  cout << "c --" #NAME "=" << solver->option (#NAME) << endl; \
} while (0);
#include "options.h"
}

/*------------------------------------------------------------------------*/

// Hast the first argument the second one as suffix?

static int suffix (const char * str, const char * sfx) {
  int l = strlen (str), k = strlen (sfx);
  if (l < k) return 0;
  return !strcmp (str + l - k, sfx);
}

// Open pipe to output of running first argument with name argument.  The
// first argument is a format string for 'sprintf' with exactly one '%s'.

static FILE * cmd (const char * fmt, const char * name) {
  FILE * res;
  char * s = new char[strlen (fmt) + strlen (name) + 1];
  sprintf (s, fmt, name);
  res = popen (s, "r");
  delete [] s;
  return res;
}

/*------------------------------------------------------------------------*/

static void err (const char * msg, const char * strarg) {
  cout << "*** cleaneling: " << msg << ' ' << strarg << endl << flush;
  exit (1);
}

static void perr (const char * msg) {
  assert (fname && lineno >= 1);
  cout << "*** parse error in '" << fname
       << "' at line " << lineno << ": " << msg << endl << flush;
  exit (1);
}

/*------------------------------------------------------------------------*/

static int next (void) {		// next char in input file
  int res = getc (file);
  if (res == '\n') lineno++;
  return res;
}

static int white_space (int ch) {
  return ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r';
}

static bool hasopt (const char *);
static void setopt (const char *);

static void parseopt (const char * line) {
  const char * p = line;
  for (p = line; *p == ' ' || *p == '\t'; p++)
    ;
  if (hasopt (p)) setopt (p);
}

static void parse () {
  int sclauses, ch, sign, lit, last, nlits, szline, nline;
  char * line;
  lineno = 1;
  line = (char*) malloc (szline = 128);
HEADER:
  if ((ch = next ()) == EOF) perr ("unexpected end-of-file");
  if (ch == 'c') {
    nline = 0;
    while ((ch = next ()) != '\n')
      if (ch == EOF) perr ("unexpected end-of-file");
      else {
	if (nline + 1 == szline) line = (char*) realloc (line, szline *= 2);
	assert (nline < szline);
	line[nline++] = ch;
      }
    assert (nline < szline);
    line[nline++] = 0;
    parseopt (line);
    goto HEADER;
  }
  if (ch != 'p') perr ("unexpected character");
  if ((fscanf (file, " cnf %d %d", &nvars, &sclauses)) != 2 ||
      nvars < 0 || sclauses < 0)
    perr ("invalid header");
  if (verbose)
    cout << "c found 'p cnf " << nvars << ' ' << sclauses << "' header"
         << endl << flush;
  nlits = last = 0;
  free (line);
BODY:
  ch = next ();
  if (ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r') goto BODY;
  if (ch == EOF) {
    if (last) perr ("zero missing");
    if (nclauses + 1 < sclauses) perr ("clauses missing");
    else if (nclauses < sclauses) perr ("clause missing");
    if (verbose)
      cout << "c parsed " << nlits << " literals in " 
           << nclauses << " clauses" << endl << flush;
    return;
  }
  if (ch == 'c') {
    while ((ch = next ()) != '\n')
      if (ch == EOF) perr ("unexpected end-of-file");
    goto BODY;
  }
  if (ch == '-') {
    if (!isdigit (ch = next ())) perr ("expected digit after '-'");
    sign = -1;
  } else if (!isdigit (ch)) perr ("expected digit or '-'");
  else sign = 1;
  lit = ch - '0';
  while (isdigit (ch = next ()))
    lit = 10*lit + (ch - '0');
  if (lit > nvars) perr ("variable exceeds maximum");
  assert (nclauses <= sclauses);
  if (nclauses == sclauses) perr ("too many clauses");
  if (lit) lit *= sign, nlits++; else nclauses++;
  if (!white_space (ch)) perr ("expected white space");
  solver->addlit (lit);
  last = lit;
  goto BODY;
}

/*------------------------------------------------------------------------*/

// Print witness to 'stdout'.  Todo more compact output (merge lines).

static void witness () {
  int W = 3;					// care for ' ' and '-'
  for (int64_t E = 10; nvars >= E; E *= 10) W++;
  const int M = 77/W;
  cout << right;
  for (int idx = 1; idx <= nvars; idx++) {
    if (!((idx - 1) % M)) {
      if (idx > 1) cout << endl;
      cout << 'v';
    }
    cout << ' ' << ((solver->value (idx) == TRUE) ? idx : -idx);
  }
  cout << endl << "v 0" << endl << flush;
}

/*------------------------------------------------------------------------*/

static bool hasopt (const char * opt, const char * name) {
  const char * p, * q;
  if (opt[0] != '-' || opt[1] != '-') return false;
  for (p = opt + 2, q = name; *q; p++, q++)
    if (*q != *p) return false;
  return *p == '=';
}

static bool hasopt (const char * str) {
#undef OPTION
#define OPTION(NAME,...) if (hasopt (str, # NAME)) return true;
#include "options.h"
  return false;
}

static void setopt (const char * opt) {
  char * tmp;
  assert (opt[0] == '-' && opt[1] == '-');
  tmp = new char[strlen (opt + 2) + 1];
  strcpy (tmp, opt + 2);
  char * val;
  for (val = tmp; *val != '='; val++) assert (*val);
  *val++ = 0;
  solver->option (tmp, atoi (val));
  delete [] tmp;
}

/*------------------------------------------------------------------------*/
}; // closing 'namespace Cleaneling'
/*------------------------------------------------------------------------*/

int main (int argc, char ** argv) {
  int res;
  solver = new_solver ();
  setsighandlers ();
  for (int i = 1; i < argc; i++) {
    if (!strcmp (argv[i], "-h") || !strcmp (argv[i], "--help")) {
      usage (); exit (0);
    } else if (!strcmp (argv[i], "-e")) { options (false); exit (0); }
    else if (!strcmp (argv[i], "-v")) setopt ("--verbose=1");
    else if (!strcmp (argv[i], "-n")) setopt ("--witness=0");
#ifndef NLGL
    else if (!strcmp (argv[i], "-c")) setopt ("--check=1");
#endif
#ifndef NLOG
    else if (!strcmp (argv[i], "-l")) setopt ("--log=1");
#endif
    else if (hasopt (argv[i])) setopt (argv[i]);
    else if (argv[i][0] == '-') {
      cout << "** cleaneling: invalid option '" << argv[i] << "'" << endl;
      exit (1);
    } else if (fname) {
      cout << "** cleaneling: can not handle multiple files '"
           << fname << "' and '" << argv[i] << "'" << endl;
      exit (1);
    } else fname = argv[i];
  }
  verbose = solver->option ("verbose");
  int clf;
  if (!fname) file = stdin, fname = "<stdin>", clf = 0;
  else if (suffix (fname, ".gz")) file = cmd ("gunzip -c %s", fname), clf = 2;
  else if (suffix (fname, ".bz2")) file = cmd ("bzcat %s", fname), clf = 2;
  else if (suffix (fname, ".7z"))
    file = cmd ("7z x -so %s 2>/dev/null", fname), clf = 2;
  else file = fopen (fname, "r"), clf = 1;
  if (!file) err ("can not read", fname);
  if (verbose) {
    banner ();
    print_line ();
    cout << "c reading " << fname << endl;
  }
  parse ();
  if (clf == 1) fclose (file);
  if (clf == 2) pclose (file);
  if (verbose) print_line (), options ();
  res = solver->solve ();
  cout << "s ";
  switch (res) {
    case UNSATISFIABLE: cout << "UNSATISFIABLE"; break;
    case SATISFIABLE: cout << "SATISFIABLE"; break;
    default: cout << "UNKNOWN"; break;
  }
  cout << endl << flush;
  if (res == SATISFIABLE && solver->option ("witness")) witness ();
  if (verbose) statistics ();
  resetsighandlers ();
  delete solver;
  if (verbose) cout << "c exit " << res << endl << flush;
  return res;
}
