#include "timer.hh"

// Petter
#ifdef USE_OPENMP
#include <omp.h>
#endif

/*
  Copyright (c) 2003-2015 Tommi Junttila
  Released under the GNU Lesser General Public License version 3.
  
  This file is part of bliss.
  
  bliss is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, version 3 of the License.

  bliss is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with bliss.  If not, see <http://www.gnu.org/licenses/>.
*/

namespace bliss {

Timer::Timer()
{
  reset();
}

void Timer::reset()
{
// Petter
  start_time = 0;
#ifdef USE_OPENMP
start_time = omp_get_wtime();
#endif
}


double Timer::get_duration()
{
	// Petter
	double intermediate = 0;
#ifdef USE_OPENMP
	intermediate = omp_get_wtime();
#endif
  return intermediate - start_time;
}

} // namespace bliss
