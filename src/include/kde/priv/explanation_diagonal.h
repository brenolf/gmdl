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

#ifndef __XOKDEPP_EXPLANATION_DIAGONAL_H__
#define __XOKDEPP_EXPLANATION_DIAGONAL_H__

#include "kde/priv/explanation_base.h"

namespace xokdepp {

  template<typename SAMPLE_PDF_TYPE>
  class explanation<SAMPLE_PDF_TYPE, vector_type> : public explanation_base<SAMPLE_PDF_TYPE, vector_type> {

    typedef gaussian<vector_type> gaussian_type;

  protected:

  protected:

    using explanation_base<SAMPLE_PDF_TYPE, vector_type>::_dims;
    using explanation_base<SAMPLE_PDF_TYPE, vector_type>::_mean;
    using explanation_base<SAMPLE_PDF_TYPE, vector_type>::_covariance;
    using explanation_base<SAMPLE_PDF_TYPE, vector_type>::_base_covariance;
    using explanation_base<SAMPLE_PDF_TYPE, vector_type>::underlyingModel;

  public:

    void reset_base_covariance() {
      _base_covariance = vector_type::Zero(_dims);
    }

  public:

    inline explanation(size_t dims) :
        explanation_base<SAMPLE_PDF_TYPE, vector_type>(dims) {
      gaussian_type::reset_mean();
      gaussian_type::reset_covariance();
      reset_base_covariance();
    }

    /**
     * for now we follow the paper and we will use only two gaussians as the underlying model
     */
    inline explanation(const vector_type &mean) :
        explanation_base<SAMPLE_PDF_TYPE, vector_type>(mean) {
      gaussian_type::reset_covariance();
      reset_base_covariance();
    }

    /**
     * for now we follow the paper and we will use only two gaussians as the underlying model
     */
    inline explanation(const vector_type &mean, const vector_type &covariance) :
        explanation_base<SAMPLE_PDF_TYPE, vector_type>(mean, covariance) {
    }

    inline explanation(const gaussian_type &gaussian) :
        explanation<SAMPLE_PDF_TYPE, vector_type>(gaussian.mean(), gaussian.covariance()) {
    }

    inline explanation(const vector_type &mean, const vector_type &base_covariance, const mixture<SAMPLE_PDF_TYPE> &underl_model) :
        explanation_base<SAMPLE_PDF_TYPE, vector_type>(mean, base_covariance, underl_model) {
    }

    //FIXME: should go away!!!
    // merge components
    explanation(const mixture<explanation<SAMPLE_PDF_TYPE, vector_type>> &mix, const vector_type &bandwidth) :
        explanation<SAMPLE_PDF_TYPE, vector_type>(mix.dims()) {

      mixture<SAMPLE_PDF_TYPE> all_details(_dims);
      for (size_t i = 0; i < mix.size(); i++) {
        for (size_t mmx = 0; mmx < mix.component(i).detailed_model().size(); mmx++) {
          all_details.add(mix.component(i).detailed_component(mmx), mix.component(i).detailed_weight(mmx) * mix.weight(i));
        }
      }

      all_details.normalize_weights_preserve_relative_importance();    //normalize the weights by w_i/sum(w_is)
      all_details.convolve(bandwidth);

      //hardcoded for now - if underlying has just one component we can't do k_means into 2 groups.
      int k_parameter = std::min(all_details.size(), 2UL);
      mixture<SAMPLE_PDF_TYPE> newDetailedModel(all_details.dims(), k_parameter);
      int mapping[all_details.size()];
      goldberger_k_means(all_details, bandwidth, newDetailedModel, mapping);

      newDetailedModel.deconvolve(bandwidth); //needed to restore original covariance instead of its temporary convolution
      all_details.deconvolve(bandwidth); //only needed if we are reusing memory and not doing copys

      explanation<SAMPLE_PDF_TYPE, vector_type> moment_matched = SAMPLE_PDF_TYPE(newDetailedModel);
      this->update(moment_matched.mean(), moment_matched.base_covariance(), moment_matched.base_covariance(), newDetailedModel);
      this->convolve(bandwidth); //FIXME this SHOULD BE OPTIONAL, BUT ITS NOT!!!
    }

    // moment matching constructor
    //okde approach - looks like it gives same results as goldberger one
    //templated to accept all types of distributions - should have implicit cast?
    explanation(const mixture<explanation<SAMPLE_PDF_TYPE, vector_type>> &mix, bool use_base_covariance = true) :
        explanation<SAMPLE_PDF_TYPE, vector_type>(mix.dims()) {

      for (size_t i = 0; i < mix.size(); i++)
        _mean += mix.weight(i) * mix.mean(i);
      _mean /= mix.sum_of_weights();

      for (size_t i = 0; i < mix.size(); i++) {
        const vector_type &cov = use_base_covariance ? mix.component(i).base_covariance() : mix.component(i).covariance();
        vector_type uuT = (mix.mean(i) * mix.mean(i).transpose()).diagonal() + cov;
        _base_covariance += mix.weight(i) * uuT;
      }

      _base_covariance -= (_mean * _mean.transpose()).diagonal();   //wj^-1 E{wi(Ei+uiuiT)}
      gaussian_type::covariance(_base_covariance);
    }

    // whitening constructor
    // FIXME: check need for dirty bit
    explanation(const explanation_base<SAMPLE_PDF_TYPE, vector_type> &to_whiten, const vector_type &smp_mean,
                const std::vector<int> &eig_ok, const std::vector<int> &/* eig_ko */, const matrix_type &F_trns) :
        explanation_base<SAMPLE_PDF_TYPE, vector_type>(eig_ok.size()) {

      for (size_t dtl_i = 0; dtl_i < to_whiten.detailed_model().size(); dtl_i++) {
        vector_type white_mean_full = F_trns * (to_whiten.detailed_mean(dtl_i) - smp_mean);
        vector_type white_mean_valid(eig_ok.size());
        select_subset_vector(white_mean_full, white_mean_valid, eig_ok);

        matrix_type cov = matrix_type::Zero(_dims, _dims);
        cov.diagonal() = to_whiten.detailed_covariance(dtl_i);

        matrix_type white_cov_full = FCFt(cov, F_trns);
        matrix_type white_cov_valid(_dims, _dims);
        select_subset_matrix(white_cov_full, white_cov_valid, eig_ok);

        underlyingModel.add(SAMPLE_PDF_TYPE(white_mean_valid, white_cov_valid.diagonal()), to_whiten.detailed_weight(dtl_i));
      }

      //recompute sample covariance using the the detailed models
      //gaussian type because detailed model is gaussian too
      gaussian_type moment_matched(underlyingModel);
      _base_covariance = moment_matched.covariance();
      gaussian_type::update(std::move(moment_matched));      // white version
    }

  public:

  };

} // namespace xokdepp

#endif /* __XOKDEPP_EXPLANATION_FULL_H__ */
