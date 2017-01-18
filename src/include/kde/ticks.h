// (C) Copyright 2014-2016, Jaime Ferreira, David Martins de Matos, Ricardo Ribeiro
// Spoken Language Systems Lab, INESC ID, IST/Universidade de Lisboa

// This file is part of XOKDE++.

// XOKDE++ is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.

// XOKDE++ is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

#ifndef __XOKDEPP_TICKS_HPP__
#define __XOKDEPP_TICKS_HPP__

#include <time.h>
#include <string>
#include <iostream>
#include <functional>

namespace xokdepp {

  inline double ticks() {
    struct timespec tv;
    if (clock_gettime(CLOCK_REALTIME, &tv) != 0)
      return 0;
    return tv.tv_sec + tv.tv_nsec / 1.0e9;
  }

  inline void time_code(const std::string &label, std::function<void()> function_to_be_timed) {
    double begin = ticks();
    std::cout << label << ": (begin) " << begin << std::endl;
    function_to_be_timed();
    double end = ticks();
    std::cout << label << ": (elapsed time) " << end - begin << std::endl;
  }

} // namespace xokdepp

#endif /* __XOKDEPP_TICKS_HPP__ */
