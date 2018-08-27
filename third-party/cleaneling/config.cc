// To break a circular compile time dependency in the library we keep the
// only compile unit which depends on 'config.h' (and thus implicitly) on
// the current 'git' ID and the current release date in a separate file.
// This way we make sure that 'config.h' can be updated whenever one of the
// files (except for this one) changes.  This could be avoided for the
// release version of the library, but still the compilation date would be
// needed to be treated in a special way.

#include "config.h"
#include "cxxflags.h"

#include <iostream>

using namespace std;

namespace Cleaneling {

void print_line (char ch);

void banner () {
  print_line ('=');
  cout << "c Cleaneling SAT Solver" << endl;
  cout << "c Version " << CLEANELING_VERSION << ' ' << CLEANELING_ID << endl;
  cout << "c Copyright (C) 2012 - 2014 Armin Biere JKU Linz" << endl;
  print_line ('=');
  cout << "c Released " << CLEANELING_RELEASED << endl;
  cout << "c Compiled " << CLEANELING_COMPILED << endl;
  cout << "c " << CLEANELING_CXX << endl;
  cout << "c " << CLEANELING_CXXFLAGS << endl;
}

};
