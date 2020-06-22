#include "internal.hpp"

#include <ctime>

extern "C" {
#include <sys/types.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
}

namespace CaDiCaL {

static uint64_t machine_identifier () {
  FILE * file = fopen ("/var/lib/dbus/machine-id", "r");
  if (!file) return 0;
  uint64_t res = 0;
  int ch;
  while ((ch = getc (file)) != '\n' && ch != EOF) {
    uint64_t tmp = res >> 56;
    res <<= 8;
    res += ch;
    if (tmp) res ^= tmp * 123123126951911u;
  }
  fclose (file);
  return res;
}

Random::Random () : state (1) {
  add (machine_identifier ());
  add (std::clock ());
  add (::time (0));
#ifndef _MSC_VER
  add (getpid ());
#endif
}

}
