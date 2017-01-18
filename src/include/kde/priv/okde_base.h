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

#ifndef __XOKDEPP_OKDE_BASE_H__
#define __XOKDEPP_OKDE_BASE_H__

#include <vector>
#include <iostream>
#include <iomanip> //debug    needed for std::setprecision()

#include "kde/basic_parameters.h"
#include "kde/basic_types.h"
#include "kde/basic_functions.h"
#include "kde/mixture.h"
#include "kde/gaussian.h"
#include "kde/explanation.h"
#include "kde/similarity_groups.h"

namespace xokdepp {

  template<typename SAMPLE_PDF_TYPE, typename COVARIANCE_TYPE>
  class oKDE_base: public mixture<explanation<SAMPLE_PDF_TYPE, COVARIANCE_TYPE>> {

    typedef gaussian<COVARIANCE_TYPE> gaussian_type;
    typedef explanation<SAMPLE_PDF_TYPE, COVARIANCE_TYPE> explanation_type;

  protected:

    using mixture<explanation_type>::_dims;
    using mixture<explanation_type>::_components;
    using mixture<explanation_type>::_weights;

    //--- error threshold
    //the paper reports 0.1 as the best for online construction of classifiers
    //this value follows the one used in the oKDE paper (2011)
    double Dth = 0.1;

    //the threshold regulating the frequency which the compression routine is called
    //uses the okde implementation "heuristic"
    const int dm;

    // the threshold regulating the frequency which the compression routine is called
    double _Mthc;

    // effective number of samples    -   Nt = N(t-1) * f + 1
    double _Nt = 0;

    //forgetting factor to apply to the effective number of samples
    double _forg = 1;  // DEFAULT: no forgetting

    //the inverse of the (sum of the (squared weights))
    data_type Nalpha = 1;

    COVARIANCE_TYPE _optimal_bandwidth;

    void updateSampleModel(); //TODO

  public:
    //TODO
    data_type compute_sample_log_probability(data_type *sample);

  public:
    void debug_set_Nt(double new_nt) {
      _Nt = new_nt;
    }

    double forg() const {
      return _forg;
    }
    void set_forg(double val) {
      _forg = val;
    }

    data_type Nt() const {
      return _Nt;
    }
    void set_Nt(double val) {
      _Nt = val;
    }

    double Mthc() const {
      return _Mthc;
    }
    void set_Mthc(double val) {
      _Mthc = val;
    }

  private:

    inline void update_Mthc() {
      if (this->size() > _Mthc) {
        _Mthc = 1.5 * _Mthc;
      } else if (this->size() < _Mthc / 2) {
        _Mthc = 0.6 * _Mthc;
      } else {
        //do nothing
      }
    }

    inline void update_Nalpha() {
      data_type w0 = 1 / _Nt;
      Nalpha = 1 / ((1 - w0) * (1 - w0) / Nalpha + w0 * w0); //definition near oKDE 2011 (29)
    }

    inline void compression_routine() {
      estimate_kernel_density();
      hierarchical_clustering();            //compress oKDE
      update_Mthc();                        //update the threshold
    }

  public:

    //TODO should be private, for ease of debug is temporarily public
    //compute optimal bandwidth and convolve the mixture with the optimal bandwidth
    void estimate_kernel_density() {
      estimate_bandwidth();
      this->convolve(_optimal_bandwidth);
    }

  public:

    inline void add_sample(vector_type &&features) {
      _Nt = _Nt * _forg + 1;
      update_Nalpha();
      data_type w0 = 1 / _Nt;
      this->scale_weights(1 - w0);
      add(explanation_type(features), w0);
      if (this->size() > _Mthc)
        compression_routine();
    }

    inline void add_sample(const vector_type &features) {
      _Nt = _Nt * _forg + 1;
      update_Nalpha();
      data_type w0 = 1 / _Nt;
      this->scale_weights(1 - w0);
      this->add(explanation_type(features), w0);
      if (this->size() > _Mthc)
        compression_routine();
    }

    inline void add_samples(std::vector<vector_type> &&samples) {
      for (auto sample : samples)
        add_sample(sample);
    }

    inline void add_samples(const std::vector<vector_type> &samples) {
      for (auto sample : samples)
        add_sample(sample);
    }

  public:
    oKDE_base(size_t dims) :
        mixture<explanation_type>(dims), dm(((_dims * _dims - _dims) / 2) + 2 * _dims + 1), _Mthc(std::min(15, dm)) {
    }

    oKDE_base(size_t dims, data_type dth) :
        mixture<explanation_type>(dims), Dth(dth), dm(((_dims * _dims - _dims) / 2) + 2 * _dims + 1), _Mthc(std::min(15, dm)) {
    }

  public:

    const COVARIANCE_TYPE &optimal_bandwidth() const {
      return _optimal_bandwidth;
    }

    virtual void estimate_bandwidth() = 0;

    virtual void hierarchical_clustering() = 0;

    //DAVID: FIXME: unfortunately, this function creates a new kde, forcing specialization!!!
//    //hierarchical_clustering implementation based on okde source code
//    void hierarchical_clustering() {
//      std::vector<int> ev_ok, ev_ko;
//      matrix_type F_trns;
//
//      explanation_type smp(*this, false);
//      smp.get_whitening_parameters(ev_ok, ev_ko, F_trns);
//      oKDE_base whitened_kde(*this, smp, ev_ok, ev_ko, F_trns);
//
//      std::vector<size_t> revitalization_list;
//      whitened_kde.revitalize(revitalization_list);
//      whitened_kde.convolve(whitened_kde.optimal_bandwidth());      //kinda wasteful, should only convolve new ones
//
//      //whitened_kde is now revitalized - lets revitalize the same components in the original kde
//      for (int idx : revitalization_list) {
//        mixture<explanation_type> splitted_component(_dims, _components[idx].detailed_model().size());
//        _components[idx].revitalize(splitted_component);
//        splitted_component.convolve(_optimal_bandwidth);
//        //replace (reuse memory) and add the second component according to the detailed model
//        data_type component_weight = _weights[idx]; //save the original component weight BEFORE replacing with splitted one
//        this->replace(idx, splitted_component.component(0), component_weight * splitted_component.weight(0)); //replace old component
//        this->add(splitted_component.component(1), component_weight * splitted_component.weight(1)); //add one more component  - this can be done like this because detailed model is of size 2
//      }
//      this->convolve(_optimal_bandwidth);    		//kinda wasteful, should only convolve new ones
//
//      std::vector<std::vector<size_t>> group_indexes = generate_similarity_groups_indexes(whitened_kde,
//                                                                                          whitened_kde.optimal_bandwidth(), Dth);
//
////      std::cout << "group indexes" << std::endl;
////      for (size_t i = 0; i < group_indexes.size(); i++) {
////        std::cout << "[ ";
////        for (size_t j = 0; j < group_indexes[i].size(); j++) {
////          std::cout << " " << group_indexes[i][j] << " ";
////        }
////        std::cout << "] ; ";
////      }
////      std::cout << std::endl;
//
//      //merge according to group_indexes
//      mixture<explanation_type> compressed(_dims);
//      COVARIANCE_TYPE opt_band = (_optimal_bandwidth / 2) * (0.5 * 0.5);  //to be coherent with what was done with the whitened data
//
//      //get sub-mixture
//      for (size_t gp = 0; gp < group_indexes.size(); gp++) {
//        mixture<explanation_type> sub_mix(_dims);
//        for (size_t i = 0; i < group_indexes[gp].size(); i++) {
//          sub_mix.add(_components[group_indexes[gp][i]], _weights[group_indexes[gp][i]]);
//        }
//        compressed.add(explanation_type(sub_mix, opt_band), sub_mix.sum_of_weights());
//      }
//
//      _components.swap(compressed.components());
//      _weights.swap(compressed.weights());
//
////      estimate_kernel_density();
//    }

    void revitalize(std::vector<size_t> &revitalization_list) {
      for (size_t cp = 0; cp < this->size(); cp++) {
        if (_components[cp].detailed_model().size() < 2) {
          // this is a singleton
          continue;
        }

        mixture<explanation_type> p0(_dims);
        p0.add(_components[cp], _weights[cp]);

        mixture<explanation_type> detailed_model(_dims);
        for (size_t i = 0; i < _components[cp].detailed_model().size(); i++) {
          detailed_model.add(explanation_type(_components[cp].detailed_mean(i), _components[cp].detailed_covariance(i)),
                             _weights[cp] * _components[cp].detailed_weight(i));
        }
        detailed_model.convolve(_optimal_bandwidth);

        data_type local_clustering_error = hellinger_distance(detailed_model, p0);

        //if local error too large, revitalize the component
        if (local_clustering_error > Dth) {
          //std::cout << "** local clustering error too large: " << local_clustering_error << std::endl;

          mixture<explanation_type> splitted_component(_dims, _components[cp].detailed_model().size());
          _components[cp].revitalize(splitted_component);
          splitted_component.convolve(_optimal_bandwidth);

          data_type component_weight = _weights[cp];
          this->replace(cp, std::move(splitted_component.component(0)), component_weight * splitted_component.weight(0));
          this->add(std::move(splitted_component.component(1)), component_weight * splitted_component.weight(1));

          //save revitalized component index
          revitalization_list.push_back(cp);
        }

      }
    }

  };

  template<typename SAMPLE_PDF_TYPE, typename COVARIANCE_TYPE>
  class oKDE: public oKDE_base<SAMPLE_PDF_TYPE, COVARIANCE_TYPE> {
    // EMPTY
  };

} // namespace xokdepp

#endif /* __XOKDEPP_OKDE_BASE_H__ */
