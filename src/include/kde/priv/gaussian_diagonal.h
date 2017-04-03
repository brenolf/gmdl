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

#ifndef __XOKDEPP_GAUSSIAN_DIAGONAL_H__
#define __XOKDEPP_GAUSSIAN_DIAGONAL_H__

#include "kde/priv/gaussian_base.h"

namespace xokdepp {

  template<>
  class gaussian<vector_type> : public gaussian_base<vector_type> {
  public:

    void reset_covariance() {
      _covariance = vector_type::Zero(_dims);
      _covariance_inverse = vector_type::Zero(_dims);
      _covariance_is_invertible = false;
      _inverse_is_dirty = false;
    }

  public:
    gaussian(size_t dims) :
        gaussian_base<vector_type>(dims) {
      reset_covariance();
    }

    gaussian(const vector_type &mu) :
        gaussian<vector_type>(mu.rows()) {
      mean(mu);
      reset_covariance();
    }

    gaussian(const vector_type &mean, const vector_type &cov) :
        gaussian_base<vector_type>(mean, cov) {
    }

    gaussian(vector_type &&mean, vector_type &&cov) :
        gaussian_base<vector_type>(mean, cov) {
    }

    gaussian(const gaussian<vector_type> &other) :
        gaussian_base<vector_type>(other) {
    }

    // gaussian constructor for mix merge into a single gaussian
    // bandwidth is actually not needed / applicable . just there to have compatible interface
    template <typename SAMPLE_PDF_TYPE>
    gaussian(const mixture<SAMPLE_PDF_TYPE> &mix, const vector_type &/*bandwidth*/) :
        gaussian<vector_type>(mix) {
    }

    // moment matching constructor from mixture
    // gaussian constructor for mix merge into a single gaussian
    template <typename SAMPLE_PDF_TYPE>
    gaussian(const mixture<SAMPLE_PDF_TYPE> &mix) :
        gaussian<vector_type>(mix.dims()) {
      data_type weight;

      if (mix.size() == 1) {      //mix is only composed by a single PDF
        weight = mix.weight(0);
        _mean = mix.mean(0);
        _covariance = mix.component(0).covariance();
        _inverse_is_dirty = true;
        return;
      }

      weight = mix.sum_of_weights();

      for (size_t i = 0; i < mix.size(); i++)
        _mean += mix.weight(i) * mix.mean(i);
      _mean /= weight;

      for (size_t i = 0; i < mix.size(); i++) {
        vector_type uuT = (mix.mean(i) * mix.mean(i).transpose()).diagonal() + mix.component(i).covariance();
        _covariance += mix.weight(i) * uuT;
      }
      _covariance /= weight;
      _covariance -= (_mean * _mean.transpose()).diagonal();
      _inverse_is_dirty = true;
    }

  public:

    //follows the okde paper approach
    void regularize_covariance_okde_approach(vector_type &covariance, Eigen::ArrayXd &ref_eig_vals) const {
      data_type dims = ref_eig_vals.rows(); //this is needed as it will be different from _DIMS when mat is degenerate

      matrix_type C = matrix_type::Zero(_dims, _dims);
      C.diagonal() = covariance + vector_type::Constant(_dims, xokdepp::epsilon_whitening);

      Eigen::EigenSolver<matrix_type> es(C, true);
      Eigen::ArrayXd eig_vals = es.eigenvalues().real().array();
      data_type max_e = -std::numeric_limits<data_type>::infinity();
      data_type min_e = std::numeric_limits<data_type>::infinity();

      for (int i = 0; i < dims; i++) {
        if (eig_vals[i] > max_e)
          max_e = eig_vals[i];
        if (eig_vals[i] < min_e)
          min_e = eig_vals[i];
      }
      Eigen::ArrayXd e((int) dims);
      for (int i = 0; i < dims; i++) {
        e[i] = eig_vals[i] / max_e;
      }

      if (min_e < xokdepp::min_coef_regularize) {        //min(e) < min_val
        std::vector<int> ev_ok, ev_ko;
        double defaultval = xokdepp::min_coef_regularize;
        for (int i = 0; i < dims; i++) {
          if (e[i] > xokdepp::min_coef_regularize)
            ev_ok.push_back(i);
          else
            ev_ko.push_back(i);
        }
        if (ev_ok.size() > 0) {
          double sum_ok = 0;
          for (int idx : ev_ok)
            sum_ok += eig_vals[idx];
          defaultval = (sum_ok / ev_ok.size()) * 1E-2;
        }

        for (int idx : ev_ko)
          eig_vals[idx] = defaultval;

        Eigen::DiagonalMatrix<data_type, Eigen::Dynamic, Eigen::Dynamic> D;
        D.diagonal() = vector_type(eig_vals);
        matrix_type V = es.eigenvectors().real();

        covariance = FCFt(D, V).diagonal(); //FIXME (diagonal should be hidden by FCFt)

        ref_eig_vals = eig_vals;        //this will allow to stop the loop when called from compute_covariance_and_inverse
      }

    }

    void compute_determinant_and_inverse() const {
    	compute_determinant_and_inverse(1);//default behaviour
    }

    void compute_determinant_and_inverse(data_type eigenvalue_relaxation_factor) const {
      if (!_inverse_is_dirty)
        return;

      if (_covariance.isZero()) {
        _covariance_inverse = vector_type::Zero(_dims);
        _covariance_log_determinant = -std::numeric_limits<data_type>::infinity();
        _covariance_is_invertible = false;
        _inverse_is_dirty = false;  //DAVID: FIXME: check
        return;
      }

      Eigen::ArrayXd ref_eig_vals = _covariance.array();
//std::cout << "COMPUTE INVERSE" << std::endl;
//std::cout << "_covariance" << std::endl << _covariance << std::endl;
//std::cout << "ref_eig_vals" << std::endl << ref_eig_vals << std::endl;
      while (!is_positive_definite(ref_eig_vals)) {
        regularize_covariance_okde_approach(_covariance, ref_eig_vals); //preserve original covariance data
        _inverse_is_dirty = true;
      }

      //we now know our covariance is invertible
      _covariance_is_invertible = true;
      _covariance_log_determinant = _covariance.array().log().sum();
      _covariance_inverse = _covariance.array().inverse();

      if (isnan(_covariance_log_determinant)) {
        /*std::cout << "gaussian:" << std::endl << *this << std::endl;
        std::cout << "_covariance" << std::endl << _covariance << std::endl;
        std::cout << "_covariance_inverse" << std::endl << _covariance_inverse << std::endl;*/
        std::cerr << "gaussian::compute_determinant_and_inverse::NaN detected" << std::endl;
        //throw std::exception();
        _covariance_log_determinant = std::numeric_limits<data_type>::infinity(); //FIXME this needs to be tested
        exit(1);
      }

      if (_covariance_log_determinant == std::numeric_limits<data_type>::infinity()
    		  || _covariance_log_determinant == -std::numeric_limits<data_type>::infinity()) {
//    	  std::cout << "gaussian::compute_determinant_and_inverse::Inf detected" << std::endl;
//    	  std::cout << "new relaxation factor = " << eigenvalue_relaxation_factor * 10 << std::endl;
    	  dynamic_covariance_relaxation(eigenvalue_relaxation_factor / 100.0);
    	  compute_determinant_and_inverse(eigenvalue_relaxation_factor * 2);//relax one order of greatness
    	  //throw std::exception();
      }


      _inverse_is_dirty = false;  //DAVID: FIXME: check
    }

    //adds small value to all entries diagonal matrix according to a percentage of the actual corresponding diagonal entry
    void dynamic_covariance_relaxation(data_type factor) const {//factor e the percentage 0.01 for 1%
    	for(int i = 0 ; i < _covariance.rows() ; i++){
    		_covariance(i) += _covariance(i)*factor;
    	}
    }

    void whitening_parameters_eigen(const vector_type &covariance, std::vector<int> &eig_ok, std::vector<int> &eig_ko,
                                    matrix_type &F, matrix_type &iF) const {
//      bool fully_singular_covariance = false;
    	matrix_type C = matrix_type::Zero(_dims, _dims);
      C.diagonal() = covariance + vector_type::Constant(_dims, xokdepp::epsilon_whitening);

      // cout << "DEBUG:COVARIANCE " << covariance << endl << endl;

      // cout << "cov" << _dims << endl;
      // cout << covariance << endl;
      // cout << C << endl << endl;

      Eigen::EigenSolver<matrix_type> es(C, true);
      Eigen::ArrayXd eig_vals = es.eigenvalues().real().array();
      matrix_type V = es.eigenvectors().real();

      //get valid/non-valid eigenvalues indexes
      for (size_t i = 0; i < _dims; i++) {
        if (eig_vals[i] > xokdepp::MIN_BANDWIDTH)
          eig_ok.push_back(i);
        else
          eig_ko.push_back(i);
      }

//      if (eig_ok.size() == 0) {
//    	  std::cerr << "Covariance matrix is completely singular " << std::endl;
//    	  eig_vals = vector_type::Constant(_dims,1).array();
//
//      }

      //get mean(valid_eig_vals)*0.01
      data_type min_eig = 0;
      for (int v_idx : eig_ok)
        min_eig += eig_vals[v_idx];

      data_type sq_min_eig = sqrt((min_eig / eig_ok.size()) * 0.01);

      Eigen::ArrayXd eig_vals_sqrt(eig_ok.size());
      for (size_t i = 0; i < eig_ok.size(); i++)
        eig_vals_sqrt[i] = 1 / sqrt((double)eig_vals[eig_ok[i]]);

      if (eig_ok.size() == 0) {
    	  std::cerr << "Covariance matrix is completely singular " << std::endl;
//    	  eig_vals = vector_type::Constant(_dims,1).array();
//    	  eig_vals_sqrt = vector_type::Constant(_dims,1/(2/xokdepp::MIN_BANDWIDTH)).array();
//    	  eig_ok.resize(0);
//    	  eig_ko.resize(0);
//    	  for(size_t i = 0 ; i < _dims ; i++)
//    		  eig_ok.push_back(i);
//    	  fully_singular_covariance = true;
      }

      Eigen::ArrayXd i_eig_vals_sqrt = eig_vals.sqrt();     //this will include even non-valid eigenvalues
      //so we correct those
      for (int nv_idx : eig_ko)
        i_eig_vals_sqrt[nv_idx] = sq_min_eig;

      //lazy way - could possibly avoid two for loops...
      F = matrix_type::Zero(_dims, eig_ok.size());
      for (size_t i = 0; i < eig_ok.size(); i++) {
        F.col(i) = V.col(eig_ok[i]) * eig_vals_sqrt(i);
      }

      iF = matrix_type(_dims, _dims);
      for (size_t i = 0; i < _dims; i++) {
        iF.col(i) = V.col(i) * i_eig_vals_sqrt(i);
      }
      iF.transposeInPlace();

//      if(fully_singular_covariance)
//    	  eig_ok.resize(0);

    }

    void whitening_parameters_svd(const vector_type &covariance, std::vector<int> &sv_ok, std::vector<int> &sv_ko,
                                  matrix_type &F_trns) const {
      matrix_type C = matrix_type::Zero(_dims, _dims);
      C.diagonal() = covariance + vector_type::Constant(_dims, xokdepp::epsilon_whitening);

      //get valid and non valid eigenvalue indexes
      //Eigen::EigenSolver<matrix_type> es(C, true);
      Eigen::JacobiSVD<matrix_type> svd(C, Eigen::ComputeFullU);
      Eigen::ArrayXd svals = svd.singularValues().real().array();
      Eigen::ArrayXd svals_norm = svals / svals.maxCoeff();

      for (size_t i = 0; i < _dims; i++) {
        if (svals_norm[i] > xokdepp::min_coef_whitening)
          sv_ok.push_back(i);
        else
          sv_ko.push_back(i);
      }

      if (sv_ok.size() == 0) {
        std::cerr << "Covariance matrix is completely singular - handle functionality not implemented" << std::endl;
        //throw std::exception(); //FIXME will not do anything, next part should handle. this should be tested
      }

      matrix_type S_inv = matrix_type::Zero(_dims, _dims);
      for (int ok_idx : sv_ok) {
        S_inv(ok_idx, ok_idx) = 1 / svals(ok_idx);
      }
      for (int ko_idx : sv_ko) {
        //slighly diff from okde code (subspacePrewhitenTransform.m:224)
        S_inv(ko_idx, ko_idx) = 1 / ko_idx;
      }

      matrix_type S_inv_abs = S_inv.array().abs();
      matrix_type S_inv_abs_sqrt = S_inv_abs.array().sqrt();

      F_trns = S_inv_abs_sqrt * svd.matrixU().inverse();
    }

    data_type KL_divergence(const gaussian_base<vector_type> &other) const {
      data_type log_det2_det1 = other.covariance_log_determinant() - covariance_log_determinant();

      matrix_type cov1 = matrix_type::Zero(_dims, _dims);
      cov1.diagonal() = covariance();
      matrix_type invcov2 = matrix_type::Zero(_dims, _dims);
      invcov2.diagonal() = other.covariance_inverse();

      matrix_type invCov2_cov1 = invcov2 * cov1;
      vector_type m1m2 = mean() - other.mean();
      data_type u1u2t_invcov2_u1u2 = FtCF(invcov2, m1m2);
      data_type result = 0.5 * (log_det2_det1 + invCov2_cov1.trace() + u1u2t_invcov2_u1u2 - (data_type)_dims);
      if (isnan(result)) {
        //std::cout << "log_det2_det1 " << log_det2_det1 << std::endl;
        //std::cout << "other.covariance_log_determinant() " << other.covariance_log_determinant() << std::endl;
        //std::cout << "covariance_log_determinant() " << covariance_log_determinant() << std::endl;
        std::cout << "gaussian::KL_divergence::NaN detected" << std::endl;
        //throw std::exception();
        result = -std::numeric_limits<data_type>::infinity();
      }
      return result;
    }

    void decompose_in_two(gaussian_base<vector_type> &g1, gaussian_base<vector_type> &g2) const {
      vector_type M = vector_type::Zero(_dims);
      M(0) = 0.5;
      vector_type C = vector_type::Zero(_dims);
      C(0, 0) = 3.0 / 4.0;

      matrix_type cov = matrix_type::Zero(_dims, _dims);
      cov.diagonal() = covariance();

      Eigen::JacobiSVD<matrix_type, Eigen::NoQRPreconditioner> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
      Eigen::DiagonalMatrix<data_type, Eigen::Dynamic, Eigen::Dynamic> sqrt_s_vals;
      sqrt_s_vals.diagonal() = svd.singularValues().array().sqrt();
      vector_type FM = svd.matrixU() * sqrt_s_vals * M;
      vector_type mu1 = _mean + FM;       //not according to official okde paper
      vector_type mu2 = _mean - FM;       //not according to official okde paper
      //matrix_type sigma = FCFt(C, FM);//not according to official okde paper
      matrix_type sigma = cov + (_mean * _mean.transpose()) - 0.5 * (mu1 * mu1.transpose() + mu2 * mu2.transpose());
      g1 = gaussian<vector_type>(mu1, sigma.diagonal());
      g2 = gaussian<vector_type>(mu2, sigma.diagonal());
    }

  };

} // namespace xokdepp

#endif /* __XOKDEPP_GAUSSIAN_DIAGONAL_H__ */
