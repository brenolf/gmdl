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

#ifndef __XOKDEPP_GAUSSIAN_BASE_H__
#define __XOKDEPP_GAUSSIAN_BASE_H__

#include "kde/basic_parameters.h"
#include "kde/basic_types.h"
#include "kde/basic_functions.h"
#include "kde/mixture.h"

namespace xokdepp {

  template<typename COVARIANCE_TYPE>
  class gaussian_base {
  public:
    size_t dims() const {
      return _dims;
    }

  protected:

    const size_t _dims;
    const double _logPowTwoPI = -log(2.0 * M_PI) * (_dims / 2.0);
    //const double _logPowTwoPI = log(pow(2.0 * M_PI, - (double)_dims / 2.0));

    vector_type _mean;

    mutable COVARIANCE_TYPE _covariance;
    mutable COVARIANCE_TYPE _covariance_inverse;

    mutable data_type _covariance_log_determinant = 0;
    mutable bool _covariance_is_invertible = false;

    mutable bool _inverse_is_dirty = false;

  public:

    void reset_mean() {
      _mean = vector_type::Zero(_dims);
    }

  protected:

    gaussian_base(size_t dims) :
        _dims(dims) {
      reset_mean();
    }

    gaussian_base(const vector_type &mean, const COVARIANCE_TYPE &cov) :
        _dims(mean.rows()), _mean(mean), _covariance(cov), _inverse_is_dirty(true) {
    }

    gaussian_base(vector_type &&mean, COVARIANCE_TYPE &&cov) :
        _dims(mean.rows()), _mean(mean), _covariance(cov), _inverse_is_dirty(true) {
    }

    gaussian_base(const gaussian_base<COVARIANCE_TYPE> &other) :
        _dims(other.dims()), _mean(other.mean()), _covariance(other.covariance()), _inverse_is_dirty(true) {
    }

  public:

    virtual ~gaussian_base() {
    }

    virtual void compute_determinant_and_inverse() const = 0;

    //protected:
    //approach used in R. paper from finance and msc thesis
    virtual void regularize_covariance_R_approach(COVARIANCE_TYPE &covariance, Eigen::EigenSolver<matrix_type> &es,
                                                  Eigen::ArrayXd &eig_vals, data_type epsilon) const {
      for (int i = 0; i < eig_vals.rows(); i++) {
        if (eig_vals[i] <= epsilon)
          eig_vals[i] += 1E-1;      //previously 1E-9
      }
      Eigen::DiagonalMatrix<data_type, Eigen::Dynamic, Eigen::Dynamic> D;
      D.diagonal() = vector_type(eig_vals);
      matrix_type V = es.eigenvectors().real();
      covariance = FCFt(D, V);
    }

    //follows the okde paper approach
    virtual void regularize_covariance_okde_approach(COVARIANCE_TYPE &covariance, Eigen::ArrayXd &ref_eig_vals) const = 0;

    virtual void whitening_parameters_eigen(const COVARIANCE_TYPE &covariance, std::vector<int> &eig_ok, std::vector<int> &eig_ko,
                                            matrix_type &F, matrix_type &iF) const = 0;

    virtual void whitening_parameters_svd(const COVARIANCE_TYPE &covariance, std::vector<int> &sv_ok, std::vector<int> &sv_ko,
                                          matrix_type &F_trns) const = 0;

  public:

    virtual void get_whitening_parameters_eigen(std::vector<int> &eig_ok, std::vector<int> &eig_ko, matrix_type &F,
                                                matrix_type &iF) const {
      whitening_parameters_eigen(covariance(), eig_ok, eig_ko, F, iF);
    }
    virtual void get_whitening_parameters(std::vector<int> &sv_ok, std::vector<int> &sv_ko, matrix_type &F_trns) const {
      whitening_parameters_svd(covariance(), sv_ok, sv_ko, F_trns);
    }

  public:
    inline void update(const vector_type &mean, const COVARIANCE_TYPE &cov) {
      _mean = mean;
      _covariance = cov;
      _inverse_is_dirty = true;
    }

    inline void update(gaussian_base<COVARIANCE_TYPE> &&gaussian) {
      if (_dims != gaussian._dims) {
        std::cerr << "gaussian::operator=  dimension mismatch: " << _dims << " vs " << gaussian.dims() << std::endl;
        //throw std::exception();
      }
      _mean = gaussian._mean;
      _covariance = gaussian._covariance;
      _inverse_is_dirty = true;
    }

    inline void operator=(const gaussian_base<COVARIANCE_TYPE> &gaussian) {
      if (_dims != gaussian._dims) {
        std::cerr << "gaussian::operator=  dimension mismatch: " << _dims << " vs " << gaussian.dims() << std::endl;
        //throw std::exception();
      }
      _mean = gaussian._mean;
      _covariance = gaussian._covariance;
      _inverse_is_dirty = true;
    }

  public:

    const vector_type &mean() const {
      return _mean;
    }
    void mean(vector_type &&mean) {
      _mean = mean;
    }
    void mean(const vector_type &mean) {
      _mean = mean;
    }

    virtual const COVARIANCE_TYPE &covariance() const {
      //compute_determinant_and_inverse();  //we don't need to regularize covariance yet
      return _covariance;
    }
    void covariance(COVARIANCE_TYPE &&cov) {
      _covariance = cov;
      _inverse_is_dirty = true;
    }
    void covariance(const COVARIANCE_TYPE &cov) {
      _covariance = cov;
      _inverse_is_dirty = true;
    }

    bool covariance_is_invertible() const {
      compute_determinant_and_inverse();
      return _covariance_is_invertible;
    }

    data_type covariance_log_determinant() const {
      compute_determinant_and_inverse();
      return _covariance_log_determinant;
    }

    const COVARIANCE_TYPE &covariance_inverse() const {
      compute_determinant_and_inverse();
      return _covariance_inverse;
    }

    inline void convolve(const COVARIANCE_TYPE &bandwidth) {
      _covariance += bandwidth;
      _inverse_is_dirty = true;
    }
    inline void deconvolve(const COVARIANCE_TYPE &bandwidth) {
      _covariance -= bandwidth;
      _inverse_is_dirty = true;
    }

  public:
//NOT USED
//
//    bool is_positive_semidefinite(const COVARIANCE_TYPE &matrix) const {
//      Eigen::EigenSolver<COVARIANCE_TYPE> es(matrix, true);
//      Eigen::ArrayXd eig_vals = es.eigenvalues().real().array();
//      return is_positive_semidefinite(eig_vals);
//    }
//
//    bool is_positive_semidefinite(const Eigen::ArrayXd &eig_vals, data_type epsilon =
//                                      std::numeric_limits<data_type>::epsilon()) const {
//      for (int i = 0; i < eig_vals.rows(); i++) {
//        if (eig_vals(i) < -epsilon)
//          return false;
//      }
//      return true;
//    }

    bool is_positive_definite(const Eigen::ArrayXd &eig_vals, data_type epsilon = xokdepp::epsilon) const {
      for (int i = 0; i < eig_vals.rows(); i++) {
        if (eig_vals(i) < epsilon)
          return false;
      }
      return true;
    }

  protected:
    data_type log_probability(const vector_type &delta) const {
      data_type result = FtCF(covariance_inverse(), delta); //the u'*E*u result
      data_type exponent = -0.5 * result; // -u'*E*u/2
      //data_type expLogDet = -covariance_log_determinant() / 2
      data_type expLogDet = -covariance_log_determinant() / 2.0;
      data_type ans = _logPowTwoPI + expLogDet + exponent;

      /*
      if(ans > 0){

    	  //std::cout << "WARNING :: log prob is positive" << std::endl;
    	  //return -std::numeric_limits<data_type>::infinity();
//    	  std::cout << "covariance_inverse:" << std::endl <<  covariance_inverse() << std::endl;
//    	  std::cout << "covariance_log_determinant: " << covariance_log_determinant() << std::endl;
//    	  std::cout << "delta:" << std::endl << delta << std::endl;
//    	  std::cout << "result: " << result << std::endl;
//    	  std::cout << std::endl;
//    	  std::cout << "_logPowTwoPI: " << _logPowTwoPI << std::endl;
//    	  std::cout << "expLogDet: " << expLogDet << std::endl;
//    	  std::cout << "exponent: " << exponent << std::endl;
//    	  std::cout << "ans: " << ans << std::endl;

//    	  const double DEBUG_logPowTwoPI = log(pow(2.0 * M_PI, - (double)_dims / 2.0));
    	  std::cout << "---" << std::endl;
    	  //data_type result = FtCF(covariance_inverse(), delta); //the u'*E*u result
    	  data_type e_exponent = exp(result * -0.5); // -u'*E*u/2
    	  //data_type expLogDet = -covariance_log_determinant() / 2
    	  data_type e_expLogDet = pow(exp(covariance_log_determinant()) , - 0.5);
    	  data_type e_ans = exp(_logPowTwoPI) * e_expLogDet * e_exponent;
//    	  std::cout << "_dims = " << _dims << std::endl;
//    	  std::cout << " _dims / 2.0 = " << (double) _dims / 2.0 << std::endl;
//    	  std::cout << "2.0 * M_PI  = " << 2.0 * M_PI << std::endl;
//    	  std::cout << "-log(2.0 * M_PI)  = " <<-log(2.0 * M_PI) << std::endl;
//    	  std::cout << "pow(2.0 * M_PI, - _dims / 2.0)  = " << pow(2.0 * M_PI, - _dims / 2.0) << std::endl;
//    	  std::cout << "log(pow(2.0 * M_PI, - _dims / 2.0))  = " << log(pow(2.0 * M_PI, - _dims / 2.0)) << std::endl;
//    	  std::cout << "DEBUG_logPowTwoPI = " << DEBUG_logPowTwoPI << std::endl;

    	  std::cout << "_logPowTwoPI: " << _logPowTwoPI << std::endl;
    	  std::cout << "e_logPowTwoPI: " << exp(_logPowTwoPI) << std::endl;
    	  std::cout << "log(e_logPowTwoPI): " << log(exp(_logPowTwoPI)) << std::endl;//should be redundant

    	  std::cout << "exponent: " << exponent << std::endl;
    	  std::cout << "e_exponent: " << e_exponent << std::endl;
    	  std::cout << "log(e_exponent): " << log(e_exponent) << std::endl;

    	  std::cout << "expLogDet: " << expLogDet << std::endl;
    	  std::cout << "e_expLogDet: " << e_expLogDet << std::endl;
    	  std::cout << "log(e_expLogDet): " << log(e_expLogDet) << std::endl;

    	  std::cout << "ans: " << ans << std::endl;
    	  std::cout << "e_ans: " << e_ans << std::endl;
    	  std::cout << "log(e_ans): " << log(e_ans) << std::endl;
    	  std::cout << "---" << std::endl;

//    	  if(expLogDet > _logPowTwoPI) {
//    		  if(allclose(delta, vector_type::Zero(_dims))){
//    			  ans = 0;
//    		  }
//    		  else {
//    			  ans = -std::numeric_limits<data_type>::infinity();
//    		  }
//    	  }
//    	  else{
//    		  throw std::exception();
//    	  }

    	  throw std::exception();
      } */
      return ans;
    }

  public:
    /**
     * WARNING: template method: "covariance"
     */
    data_type log_likelihood(const vector_type &observation) const {
      return log_probability(_mean - observation);
    }

    //observation nElems MUST be equal to gaussian nElems //should I check that?
    data_type likelihood(const vector_type &observation) const {
      return exp(log_likelihood(observation));
    }

    virtual data_type KL_divergence(const gaussian_base<COVARIANCE_TYPE> &other) const = 0;

    virtual void decompose_in_two(gaussian_base<COVARIANCE_TYPE> &g1, gaussian_base<COVARIANCE_TYPE> &g2) const = 0;

    friend std::ostream &operator<<(std::ostream &o, const gaussian_base<COVARIANCE_TYPE> &g) {
      o << "Mean" << std::endl << g.mean() << std::endl;
      o << "Covariance" << std::endl << g.covariance() << std::endl;
      return o;
    }

    friend bool operator==(const gaussian_base<COVARIANCE_TYPE> &lv, const gaussian_base<COVARIANCE_TYPE> &rv) {
      return allclose(lv.mean(), rv.mean()) && allclose(lv.covariance(), rv.covariance());
    }

  };

  template<typename COVARIANCE_TYPE>
  class gaussian: public gaussian_base<COVARIANCE_TYPE> {
    // EMPTY
  };

} // namespace xokdepp

#endif /* __XOKDEPP_GAUSSIAN_BASE_H__ */

