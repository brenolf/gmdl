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

#ifndef __XOKDEPP_BASIC_PARAMETERS_H__
#define __XOKDEPP_BASIC_PARAMETERS_H__

#include "kde/basic_types.h"

namespace xokdepp {

  // how close are the same double values?
  const data_type precision = 1e-12;

  const data_type epsilon = 1e-9;
  //const data_type epsilon = 10 * std::numeric_limits<data_type>::epsilon();

  const data_type epsilon12 = 1e-12;
  //const data_type epsilon12 = 10 * std::numeric_limits<data_type>::epsilon();

  const data_type min_coef_regularize = 1e-5;
  const data_type epsilon_regularize = 1e-10;      //epsilon to add to the covariance

  const data_type min_coef_whitening = 1e-7;
  const data_type epsilon_whitening = 1e-7;      //epsilon to add to the covariance

  const data_type MIN_BANDWIDTH = 1e-9; // minimum eigenvalues/singular values (bandwidth)
  //const data_type min_bw = xokdepp::epsilon12; // minimum eigenvalues/singular values (bandwidth)

}

#endif /* __XOKDEPP_PARAMETERS_H__ */
