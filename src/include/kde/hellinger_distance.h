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

#ifndef __XOKDEPP_HELLINGER_DISTANCE_H__
#define __XOKDEPP_HELLINGER_DISTANCE_H__

#include <assert.h>
#include <vector>
#include "kde/basic_parameters.h"
#include "kde/basic_types.h"
#include "kde/basic_functions.h"
#include "kde/mixture.h"

namespace xokdepp {

  inline Eigen::JacobiSVD<matrix_type> get_svd_solver(const matrix_type &covariance) {
    return Eigen::JacobiSVD<matrix_type>(covariance, Eigen::ComputeFullU | Eigen::ComputeFullV);
  }

  inline Eigen::JacobiSVD<matrix_type> get_svd_solver(const vector_type &covariance) {
    size_t dims = covariance.rows();
    matrix_type cov = matrix_type::Zero(dims, dims);
    cov.diagonal() = covariance;
    return Eigen::JacobiSVD<matrix_type>(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
  }

  //oKDE 2011 paper - Unscented Hellinger Distance  (37,38,39)
  //DAVID: mixtures are assumed to have the same dimension (in terms of feature space)
  template<typename _PDF1>
  data_type hellinger_distance(const mixture<_PDF1> &one, const mixture<_PDF1> &other) {
    const size_t dims = one.dims();

    mixture<_PDF1> importance_mix(one);
    importance_mix.addAll(other);
    importance_mix.scale_weights(0.5);

    //as long as it's not negative - since we know nDims will be larger than 3 then k will mostly be 0
    size_t k = std::max(0, 3 - (int)dims);

    data_type nk = dims + k;
    data_type wi0 = k / nk;
    data_type wi = 1 / (2 * nk);

    data_type global_sum = 0.0;

    for (size_t i = 0; i < importance_mix.size(); i++) {
      Eigen::JacobiSVD<matrix_type> svd = get_svd_solver(importance_mix.component(i).covariance());

      vector_type s = nk * svd.singularValues();
      Eigen::DiagonalMatrix<data_type, Eigen::Dynamic, Eigen::Dynamic> sqrt_s;
      sqrt_s.diagonal() = s.array().sqrt();

      matrix_type ucols = svd.matrixU() * sqrt_s;

      //compute j=0
      data_type inner_sum = 0.0;
      if (wi0 > 0) {
        // compute only if wi0 is not 0, otherwise we keep innerSum initialized at 0.0
        data_type p0 = importance_mix.likelihood(importance_mix.mean(i));
        data_type p1 = sqrt(one.likelihood(importance_mix.mean(i)));
        data_type p2 = sqrt(other.likelihood(importance_mix.mean(i)));
        data_type diff = p1 - p2;
        data_type gx = diff * diff / p0;
        inner_sum = gx * wi0;
      }

      //compute remaining j's //we make use of the loop and compute both -1 and 1 sj (39)
      for (int j = 0; j < dims; j++) {       //nDims*2  since we compute 1 and -1
        vector_type means_plus_ucols = importance_mix.mean(i) + ucols.col(j);
        vector_type means_minus_ucols = importance_mix.mean(i) - ucols.col(j);

        //lets compute all probabilities first
        data_type p0_plus = importance_mix.likelihood(means_plus_ucols);
        data_type p1_plus = sqrt(one.likelihood(means_plus_ucols));
        data_type p2_plus = sqrt(other.likelihood(means_plus_ucols));

        //this computes ((sqrt(p1)-sqrt(p2))^2)/p0   (37) //carefull, can divide by zero if p0_plus is suficiently negative
        data_type diff_plus = p1_plus - p2_plus;

        //here for numeric stability
        if(p1_plus == std::numeric_limits<data_type>::infinity() && p2_plus == std::numeric_limits<data_type>::infinity()){
        	//std::cout << "solving inf - inf indetermination" << std::endl;
        	diff_plus = 0;
        }

        data_type gx_plus = diff_plus * diff_plus / p0_plus;
        if (isnan(gx_plus))       // 0*0/0 causes a nan, lets assume result is 0
          gx_plus = 0;
        inner_sum += gx_plus * wi;

        data_type p0_minus = importance_mix.likelihood(means_minus_ucols);
        data_type p1_minus = sqrt(one.likelihood(means_minus_ucols));
        data_type p2_minus = sqrt(other.likelihood(means_minus_ucols));
        data_type diff_minus = p1_minus - p2_minus;

        if(p1_minus == std::numeric_limits<data_type>::infinity() && p2_minus == std::numeric_limits<data_type>::infinity()){
        	//std::cout << "solving inf - inf indetermination" << std::endl;
        	diff_plus = 0;
        }

        data_type gx_minus = diff_minus * diff_minus / p0_minus;
        if (isnan(gx_minus))       // 0*0/0 causes a nan, lets assume result is 0
          gx_minus = 0;
        inner_sum += gx_minus * wi;

        if (isnan(inner_sum)) {
          std::cout << "Tools::hellingerDistance::innerSum::NaN_detected at dim j[" << j << "]" << std::endl;
          std::cout << "p0_plus: " << p0_plus << "  p1_plus: " << p1_plus << "  p2_plus: " << p2_plus << std::endl;
          std::cout << "gx_plus: " << gx_plus << "  wi: " << wi << std::endl;
          std::cout << "p0_minus: " << p0_minus << "  p1_minus: " << p1_minus << "  p2_minus: " << p2_minus << std::endl;
          std::cout << "gx_minus: " << gx_minus << "  wi: " << wi << std::endl;
          std::cout << "inner_sum: " << inner_sum << std::endl;
          exit(0);
        }
      }

      if (isnan(inner_sum)) {
        std::cout << "Tools::hellingerDistance::innerSum::NaN_detected" << std::endl;
        exit(0);
      }

      global_sum += importance_mix.weight(i) * inner_sum;
    }

    data_type distance = sqrt(0.5 * global_sum);

    if (isnan(distance)) {
      std::cout << "Tools::hellingerDistance::finalDistance::NaN_detected" << std::endl;
      exit(0);
    }

    if(distance > 1){
    	std::cout << "Tools::hellingerDistance::finalDistance   hellDistance must be [0,1]" << std::endl;
    	distance = 1;
    }

    return distance;
  }

} // namespace xokdepp

#endif /* __XOKDEPP_HELLINGER_DISTANCE_H__ */
