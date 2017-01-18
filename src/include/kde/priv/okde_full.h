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

#ifndef __XOKDEPP_OKDE_FULL_H__
#define __XOKDEPP_OKDE_FULL_H__

#include "kde/priv/okde_base.h"

namespace xokdepp {

  template<typename SAMPLE_PDF_TYPE>
  class oKDE<SAMPLE_PDF_TYPE, matrix_type> : public oKDE_base<SAMPLE_PDF_TYPE, matrix_type> {

    typedef gaussian<matrix_type> gaussian_type;
    typedef explanation<SAMPLE_PDF_TYPE, matrix_type> explanation_type;

  protected:

  protected:

    using oKDE_base<SAMPLE_PDF_TYPE, matrix_type>::_dims;

    using oKDE_base<SAMPLE_PDF_TYPE, matrix_type>::_components;
    using oKDE_base<SAMPLE_PDF_TYPE, matrix_type>::_weights;
    using oKDE_base<SAMPLE_PDF_TYPE, matrix_type>::_optimal_bandwidth;

    using oKDE_base<SAMPLE_PDF_TYPE, matrix_type>::Dth;
    using oKDE_base<SAMPLE_PDF_TYPE, matrix_type>::_Nt;
    using oKDE_base<SAMPLE_PDF_TYPE, matrix_type>::Nalpha;
    using oKDE_base<SAMPLE_PDF_TYPE, matrix_type>::_Mthc;
    using oKDE_base<SAMPLE_PDF_TYPE, matrix_type>::_forg;

  public:
    oKDE(size_t dims) :
        oKDE_base<SAMPLE_PDF_TYPE, matrix_type>(dims) {
    }

    oKDE(size_t dims, data_type dth) :
        oKDE_base<SAMPLE_PDF_TYPE, matrix_type>(dims, dth) {
    }

    // whitening constructor
    oKDE(const oKDE<SAMPLE_PDF_TYPE, matrix_type> &to_whiten, const vector_type &smp_mean, const std::vector<int> &eig_ok,
         const std::vector<int> &eig_ko, const matrix_type &F_trns) :
        oKDE_base<SAMPLE_PDF_TYPE, matrix_type>(eig_ok.size()) {
      size_t whitened_dims = eig_ok.size();
      Dth = to_whiten.Dth;
      _Nt = to_whiten.Nt();
      Nalpha = to_whiten.Nalpha;
      _Mthc = to_whiten.Mthc();
      _forg = to_whiten.forg();

      matrix_type new_full_h = FCFt(to_whiten.optimal_bandwidth(), F_trns);

      matrix_type new_valid_h(whitened_dims, whitened_dims);
      select_subset_matrix(new_full_h, new_valid_h, eig_ok);
      _optimal_bandwidth = new_valid_h;

      //apply transform to sample model
      for (size_t i = 0; i < to_whiten.size(); i++) {
        this->add(explanation_type(to_whiten.component(i), smp_mean, eig_ok, eig_ko, F_trns), to_whiten.weight(i));
        this->convolve(_optimal_bandwidth);
      }
    }

  public:

    data_type compute_rpfg(const mixture<explanation_type> &mix, const matrix_type &F, data_type N) const {
      const size_t dims = F.rows();
      const data_type f = pow(4.0 / ((dims + 2.0) * N), 2.0 / (dims + 4.0));  //THIS IS DIFFERENT FROM PAPER - FOLLOWS KDE CODE IMPL
      const matrix_type G = F * f;
      const matrix_type Fsq = F * F;

      data_type Rpfg = 0.0; //R(p,F,G) final result
      for (size_t i = 0; i < mix.size(); i++) {
        const data_type wi = mix.weight(i);
        explanation_type e_base(mix.mean(i), mix.component(i).base_covariance());
        e_base.convolve(G); //Egi

        for (size_t j = i; j < mix.size(); j++) {
          const data_type wj = mix.weight(j);
          explanation_type ei(e_base); //copy
          ei.convolve_kde(mix.component(j).base_covariance()); // follows the paper definition //adds Egi + Ej

          const matrix_type &Aij = ei.covariance_inverse();
          const matrix_type mulF_Aij = F * Aij;
          const matrix_type mulFsq_Aijsq = Fsq * (Aij * Aij); //F^2 * Aij^2
          const vector_type dij = mix.mean(i) - mix.mean(j);
          const data_type mij = FtCF(Aij, dij);

          //part 1 = wiwjGAUSS(dij)
          //this would be the elegant way, but _dims may not be known at compile time, so it might not be correct
          const data_type p1 = wi * wj * ei.likelihood(mix.mean(j));
          //part 2 = 2tr(...)(1-2mij)
          const data_type p2 = 2 * mulFsq_Aijsq.trace() * (1 - (2 * mij));
          //part 3 = tr^2(FAij)(1-mij)^2
          const data_type mulF_Aij_tr = mulF_Aij.trace();
          const data_type p3 = (mulF_Aij_tr * mulF_Aij_tr) * ((1 - mij) * (1 - mij));

          //since they are symmetric (ij[0][1] = ij[1][0]) we don't need to compute twice
          //i==j is the special case where ij[0][0] or [1][1]....
          const data_type mult = (i == j) ? 1 : 2;
          Rpfg += p1 * (p2 + p3) * mult;
        }
      }
      return Rpfg;
    }

    //according to the oKDE code implementation
    void estimate_bandwidth() {
      const data_type N = _Nt;

      //H=b^2F   F=Sigma_smp
      //b_opt = d((4PI)^(d/2)) N_alpha R(p,F,G)     (6)

      explanation_type moment_matched(*this, true);
      gaussian_type smp(moment_matched.mean(), moment_matched.base_covariance());

      //compute G 2011-(13)
      matrix_type F, iF;
      std::vector<int> eig_ok, eig_ko;
      smp.get_whitening_parameters_eigen(eig_ok, eig_ko, F, iF);

      matrix_type F_smp_cov = FtCF(smp.covariance(), F);
      gaussian_type F_smp(vector_type::Zero(F_smp_cov.rows()), F_smp_cov);

//      mixture<explanation_type> n_mix(eig_ok.size(), *this, smp, F);
      mixture<explanation_type> n_mix(eig_ok.size(), *this, moment_matched, F);
      data_type Rpfg = compute_rpfg(n_mix, F_smp.covariance(), N); //R(p,F,G) final result

      // inner part of equation (6) //we square now to avoid having negative values when we exponentiate
      data_type dims = F_smp.covariance().rows(); //lets avoid problems with type of _dims (size_t)
      data_type F_smp_cov_det = exp(F_smp.covariance_log_determinant() * -0.5); //FIXME this may give numeric problems!!!
      data_type innerBandwidthOptScale = (F_smp_cov_det / N) / (pow(sqrt(4.0 * M_PI), dims) * dims * Rpfg);
      data_type bandwidthOptScaled = pow(innerBandwidthOptScale, (1 / (dims + 4.0)));
      data_type bandwidthOptScalesquared = bandwidthOptScaled * bandwidthOptScaled;

      //BEGIN VODOO
      matrix_type H = matrix_type::Zero(_dims, _dims);
      H.diagonal() = vector_type::Constant(_dims, 1);

      if (bandwidthOptScalesquared != 0) {
        for (int v_idx : eig_ok)
          H(v_idx, v_idx) = bandwidthOptScalesquared;
      } else {
    	  std::cout << "bandwidthOptScalesquared is zero.";// << std::endl;
    	  //std::cout << " Correcting..." << std::endl;
    	  //Ok, bandwidth is in havoc, lets make something better.
    	  //h=En^(-1/5)
    	  //for (int v_idx : eig_ok)
    	//	  H(v_idx, v_idx) = 1/10;
        //  _optimal_bandwidth = H;
      }

      _optimal_bandwidth = FtCF(H, iF);
      //END V0D00
    }

    //hierarchical_clustering implementation based on okde source code
    void hierarchical_clustering() {
      std::vector<int> ev_ok, ev_ko;
      matrix_type F_trns;

      explanation_type smp(*this, false);
      smp.get_whitening_parameters(ev_ok, ev_ko, F_trns);
      oKDE<SAMPLE_PDF_TYPE, matrix_type> whitened_kde(*this, smp.mean(), ev_ok, ev_ko, F_trns);

      std::vector<size_t> revitalization_list;
      whitened_kde.revitalize(revitalization_list);
      whitened_kde.convolve(whitened_kde.optimal_bandwidth());      //kinda wasteful, should only convolve new ones

      //whitened_kde is now revitalized - lets revitalize the same components in the original kde
      for (int idx : revitalization_list) {
        mixture<explanation_type> splitted_component(_dims, _components[idx].detailed_model().size());
        _components[idx].revitalize(splitted_component);
        splitted_component.convolve(_optimal_bandwidth);
        //replace (reuse memory) and add the second component according to the detailed model
        data_type component_weight = _weights[idx]; //save the original component weight BEFORE replacing with splitted one
        this->replace(idx, splitted_component.component(0), component_weight * splitted_component.weight(0)); //replace old component
        this->add(splitted_component.component(1), component_weight * splitted_component.weight(1)); //add one more component  - this can be done like this because detailed model is of size 2
      }
      this->convolve(_optimal_bandwidth);           //kinda wasteful, should only convolve new ones

      std::vector<std::vector<size_t>> group_indexes = generate_similarity_groups_indexes(whitened_kde,
                                                                                          whitened_kde.optimal_bandwidth(), Dth);

//      std::cout << "group indexes" << std::endl;
//      for (size_t i = 0; i < group_indexes.size(); i++) {
//        std::cout << "[ ";
//        for (size_t j = 0; j < group_indexes[i].size(); j++) {
//          std::cout << " " << group_indexes[i][j] << " ";
//        }
//        std::cout << "] ; ";
//      }
//      std::cout << std::endl;

      //merge according to group_indexes
      mixture<explanation_type> compressed(_dims);
      matrix_type opt_band = (_optimal_bandwidth / 2) * (0.5 * 0.5);  //to be coherent with what was done with the whitened data

      //get sub-mixture
      for (size_t gp = 0; gp < group_indexes.size(); gp++) {
        mixture<explanation_type> sub_mix(_dims);
        for (size_t i = 0; i < group_indexes[gp].size(); i++) {
          sub_mix.add(_components[group_indexes[gp][i]], _weights[group_indexes[gp][i]]);
        }
        compressed.add(explanation_type(sub_mix, opt_band), sub_mix.sum_of_weights());
      }

      _components.swap(compressed.components());
      _weights.swap(compressed.weights());

//      estimate_kernel_density();
    }


  };

} // namespace xokdepp

#endif /* __XOKDEPP_OKDE_FULL_H__ */
