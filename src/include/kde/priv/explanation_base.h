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

#ifndef __XOKDEPP_EXPLANATION_BASE_H__
#define __XOKDEPP_EXPLANATION_BASE_H__

#include "kde/basic_parameters.h"
#include "kde/basic_types.h"
#include "kde/basic_functions.h"
#include "kde/priv/gaussian_base.h"
#include "kde/mixture.h"
#include "kde/goldberger_k_means.h"

namespace xokdepp {

  template<typename SAMPLE_PDF_TYPE, typename COVARIANCE_TYPE>
  class explanation_base: public gaussian<COVARIANCE_TYPE> {

    typedef gaussian<COVARIANCE_TYPE> gaussian_type;

  protected:

    using gaussian_type::_dims;
    using gaussian_type::_mean;
    using gaussian_type::_covariance;

    // the underlying model associated with the explanation
    mixture<SAMPLE_PDF_TYPE> underlyingModel;

    //the covariance matrix of the component before the convolution with the optimal bandwidth
    COVARIANCE_TYPE _base_covariance;

  protected:

    inline explanation_base(size_t dims) :
        gaussian_type(dims), underlyingModel(dims) {
    }

    /**
     * for now we follow the paper and we will use only two gaussians as the underlying model
     */
    inline explanation_base(const vector_type &mean) :
        gaussian_type(mean), underlyingModel(_dims) {
      underlyingModel.add(SAMPLE_PDF_TYPE(mean));
    }

    /**
     * for now we follow the paper and we will use only two gaussians as the underlying model
     */
    inline explanation_base(const vector_type &mean, const COVARIANCE_TYPE &covariance) :
        gaussian_type(mean, covariance), underlyingModel(_dims), _base_covariance(covariance) {
      underlyingModel.add(SAMPLE_PDF_TYPE(mean, covariance));
    }

    inline explanation_base(const gaussian_type &gaussian) :
        explanation_base(gaussian.mean(), gaussian.covariance()) {
    }

    inline explanation_base(explanation_base &&other) :
        gaussian_type(other.mean(), other.covariance()), underlyingModel(other.detailed_model()), _base_covariance(
            other.base_covariance()) {
    }

    // copy constructor
    inline explanation_base(const explanation_base &other) :
        gaussian_type(other.mean(), other.covariance()), underlyingModel(other.detailed_model()), _base_covariance(
            other.base_covariance()) {
    }

    inline explanation_base(const vector_type &mean, const COVARIANCE_TYPE &base_covariance,
                            const mixture<SAMPLE_PDF_TYPE> &underl_model) :
        gaussian_type(mean, base_covariance), underlyingModel(underl_model), _base_covariance(base_covariance) {
    }

  public:

    virtual ~explanation_base() {
    }

  public:

    const COVARIANCE_TYPE &covariance() const {
      gaussian_type::compute_determinant_and_inverse();
      return _covariance;
    }

    void get_whitening_parameters(std::vector<int> &sv_ok, std::vector<int> &sv_ko, matrix_type &F_trns) const {
//      whitening_parameters_svd(base_covariance(), sv_ok, sv_ko, F_trns);
      gaussian_type::whitening_parameters_svd(covariance(), sv_ok, sv_ko, F_trns);
    }

    inline mixture<SAMPLE_PDF_TYPE> &detailed_model() {
      return underlyingModel;
    }
    inline const mixture<SAMPLE_PDF_TYPE> &detailed_model() const {
      return underlyingModel;
    }

    const COVARIANCE_TYPE &base_covariance() const {
      return _base_covariance;
    }
    void base_covariance(COVARIANCE_TYPE &&cov) {
      _base_covariance = cov;
    }
    void base_covariance(const COVARIANCE_TYPE &cov) {
      _base_covariance = cov;
    }

    inline void update(const vector_type &mean, const COVARIANCE_TYPE &covariance, const COVARIANCE_TYPE &base_covariance,
                       const mixture<SAMPLE_PDF_TYPE> &newUnderlyingModel) {
      gaussian_type::update(mean, covariance);
      _base_covariance = base_covariance;
      underlyingModel = newUnderlyingModel;
    }

    void operator=(const explanation_base &explanation) {
      gaussian_type::operator=(explanation);
      _base_covariance = explanation._base_covariance;
      underlyingModel = explanation.underlyingModel;
    }

    inline void convolve(const COVARIANCE_TYPE &bandwidth) {
      gaussian_type::covariance(_base_covariance + bandwidth);
    }
    inline void deconvolve(const COVARIANCE_TYPE &bandwidth) {
      gaussian_type::covariance(_base_covariance - bandwidth);
    }

    inline void convolve_kde(const COVARIANCE_TYPE &bandwidth) {
      gaussian_type::covariance(covariance() + bandwidth);
    }
    inline void deconvolve_kde(const COVARIANCE_TYPE &bandwidth) {
      gaussian_type::covariance(covariance() - bandwidth);
    }

    data_type detailed_weight(int ix) const {
      return underlyingModel.weight(ix);
    }
    const SAMPLE_PDF_TYPE &detailed_component(int ix) const {
      return underlyingModel.component(ix);
    }
    const vector_type &detailed_mean(int ix) const {
      return underlyingModel.mean(ix);
    }
    const COVARIANCE_TYPE &detailed_covariance(int ix) const {
      return underlyingModel.component(ix).covariance();
    }

    template<typename PDF>
    void revitalize(mixture<PDF> &splitted_mix) const {
      if (underlyingModel.size() == 1) {      //can't split lets represent the component by its own detailed model
        splitted_mix.replace(0, PDF(detailed_mean(0), detailed_covariance(0), underlyingModel), detailed_weight(0));
        return;
      }

      //the number of components we want to split the components that need revitalization.
      //this value should be equal to the number of components of the detailed model
      size_t splitted_mix_size = 2;

      for (size_t i = 0; i < splitted_mix_size; i++) {
        mixture<SAMPLE_PDF_TYPE> detailedModel(_dims);
        //split component using each underlying model component
        //TODO underlying model should be a mixture, so this would be more like a copy

        //test if detailed component is a dirac
        if (detailed_covariance(i).isZero()) {
          detailedModel.add(SAMPLE_PDF_TYPE(detailed_mean(i), detailed_covariance(i)), 1 /*detailed_weight(i)*/);
          splitted_mix.replace(i, PDF(detailed_mean(i), detailed_covariance(i), detailedModel), detailed_weight(i));
        } else {
          SAMPLE_PDF_TYPE split1(_dims), split2(_dims);
          detailed_component(i).decompose_in_two(split1, split2);
          detailedModel.add(std::move(split1), 0.5);
          detailedModel.add(std::move(split2), 0.5);

          splitted_mix.replace(i, PDF(detailed_mean(i), detailed_covariance(i), detailedModel), detailed_weight(i));
        }
      }
    }

    void decompose_in_two(explanation_base &e1, explanation_base &e2) const {
      gaussian_type g1(_dims), g2(_dims);
      gaussian_type::decompose_in_two(g1, g2);
      e1 = explanation_base<SAMPLE_PDF_TYPE, COVARIANCE_TYPE>(g1);
      e2 = explanation_base<SAMPLE_PDF_TYPE, COVARIANCE_TYPE>(g2);
    }

    friend std::ostream &operator<<(std::ostream &o, const explanation_base &e) {
      o << "----- explanation -----" << std::endl;
      o << (gaussian_type&)e << "Base covariance" << std::endl;
      o << e.base_covariance() << std::endl << std::endl;
      o << "Underlying model" << std::endl;
      o << "---" << e.underlyingModel << "---" << std::endl << "-----             -----" << std::endl;
      return o;
    }

  };

  template<typename SAMPLE_PDF_TYPE, typename COVARIANCE_TYPE>
  class explanation: public explanation_base<SAMPLE_PDF_TYPE, COVARIANCE_TYPE> {
    //EMPTY
  };

} // namespace xokdepp

#endif /* __XOKDEPP_EXPLANATION_BASE_H__ */
