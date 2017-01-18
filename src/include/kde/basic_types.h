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

#ifndef __XOKDEPP_BASIC_TYPES_H__
#define __XOKDEPP_BASIC_TYPES_H__

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <unsupported/Eigen/MatrixFunctions>

//needed for drawing
#include <unsupported/Eigen/AutoDiff>
#include <Eigen/Eigenvalues>

namespace xokdepp {

  typedef double data_type;

  typedef Eigen::Matrix<data_type, Eigen::Dynamic, 1> vector_type;

  typedef Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> matrix_type;

}

#endif /* __XOKDEPP_BASIC_TYPES_H__ */
