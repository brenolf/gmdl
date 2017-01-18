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

#ifndef __XOKDEPP_BASIC_FUNCTIONS_H__
#define __XOKDEPP_BASIC_FUNCTIONS_H__

#include <random>
#include <iostream>
#include <functional>
#include <algorithm>
#include <cmath>

#include "kde/basic_parameters.h"
#include "kde/basic_types.h"

namespace xokdepp {

  template<typename DerivedA, typename DerivedB>
  bool allclose(const Eigen::DenseBase<DerivedA>& a, const Eigen::DenseBase<DerivedB>& b,
                const typename DerivedA::RealScalar& rtol = xokdepp::epsilon, const typename DerivedA::RealScalar& atol =
                    xokdepp::epsilon12) {
    return ((a.derived() - b.derived()).array().abs() <= (atol + rtol * b.derived().array().abs())).all();
  }

  inline void select_subset_matrix(const matrix_type &src, matrix_type &dst, const std::vector<int> &ok_id) {
    for (size_t i = 0; i < ok_id.size(); i++) {
      for (size_t j = 0; j < ok_id.size(); j++) {
        dst(i, j) = src(ok_id[i], ok_id[j]);
      }
    }
  }

  inline void select_subset_vector(const vector_type &src, vector_type &dst, const std::vector<int> &ok_id) {
    for (size_t i = 0; i < ok_id.size(); i++) {
      dst(i) = src(ok_id[i]);
    }
  }

  //so inneficient...
  inline data_type vector_max_val_and_idx(const std::vector<data_type> &vector, size_t &idx) {
    auto result = std::max_element(vector.begin(), vector.end());
    idx = std::distance(vector.begin(), result);
    return *result;
  }

  inline bool double_equals(double a, double b, double precision = xokdepp::precision) {
    return std::abs(a - b) < precision;
  }

  // compute F'*cov*F (diagonal covs)
  inline data_type FtCF(const vector_type &C, const vector_type &F) {
    /*size_t dims = F.rows();
    matrix_type cov = matrix_type::Zero(dims, dims);
    cov.diagonal() = C;
    return F.transpose() * cov * F;*/
    return (F.array() * C.array() * F.array()).sum();
  }

  // compute F'*cov*F
  inline data_type FtCF(const matrix_type &cov, const vector_type &F) {
    return F.transpose() * cov * F;
  }

  // compute F'*cov*F (diagonal covs)
  inline vector_type FtCF(const vector_type &C, const matrix_type &F) {
    size_t dims = F.rows();
    matrix_type cov = matrix_type::Zero(dims, dims);
    cov.diagonal() = C;
    return (F.transpose() * cov * F).diagonal();
  }

  // compute F'*cov*F
  inline matrix_type FtCF(const matrix_type &cov, const matrix_type &F) {
    return F.transpose() * cov * F;
  }

  // compute F'*cov*F
  inline matrix_type FCFt(const matrix_type &cov, const matrix_type &F) {
    return F * cov * F.transpose();
  }


}

#endif /* __XOKDEPP_BASIC_FUNCTIONS_H__ */
