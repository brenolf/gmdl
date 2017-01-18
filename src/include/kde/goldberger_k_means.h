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

#ifndef __XOKDEPP_GOLBERGER_K_MEANS_H__
#define __XOKDEPP_GOLBERGER_K_MEANS_H__

#include <assert.h>
#include <vector>
#include "kde/basic_parameters.h"
#include "kde/basic_types.h"
#include "kde/basic_functions.h"
#include "kde/mixture.h"

namespace xokdepp {

  template<typename _PDF>
  bool regroup(const mixture<_PDF> &mix, mixture<_PDF> &compressed_mix, int *mapping, std::vector<int> &visited_mapping,
               bool &forcestop) {
    assert(mix.size() > 0);

    data_type max_dist = -std::numeric_limits<data_type>::infinity();
    data_type max_dist_idx = -1;      //maximum distance index

    //TODO for simplicity I will assume that k-means has k=2
    assert(compressed_mix.size() == 2); //at this point the mix must have this size

    //count the assignments to each class
    int c0 = 0, c1 = 0;

    //regroup
    for (size_t i = 0; i < mix.size(); i++) {

      data_type min_dist_0 = fabs(compressed_mix.component(0).KL_divergence(mix.component(i)));
      data_type min_dist_1 = fabs(compressed_mix.component(1).KL_divergence(mix.component(i)));
      data_type internal_max_dist = std::max(min_dist_0, min_dist_1);

      if (internal_max_dist > max_dist) {
        max_dist = internal_max_dist;
        max_dist_idx = i;
      }

      if (min_dist_1 < min_dist_0) {
        mapping[i] = 1;
        visited_mapping[i] = 1;
        c1++;
      } else {
        mapping[i] = 0;
        visited_mapping[i] = 0;
        c0++;
      }

    }

    if (max_dist_idx == -1) {
      std::cout << "mixture::goldberger_k_means - max_dist_idx has not been correctly found. Aborting..." << std::endl;
      //throw std::exception();
      compressed_mix = mix;
      return true; //error flag
    }

    //check mappings for cases where one class is always selected
    //if detected, correct by assigning the first component to the other class (the "empty" one)
    //this follows the oKDE approach
    if (c0 == 0) {
      //put the less similar (i.e. the maxdist one) to the other class, and leave the rest in this one
      mapping[0] = 0;
      visited_mapping[max_dist_idx] = 0;
      forcestop = true;
    } else if (c1 == 0) {
      mapping[0] = 1;
      visited_mapping[max_dist_idx] = 1;
      forcestop = true;
    }

    return false ;//no errors occured
  }

  template<typename _PDF>
  void refit(const mixture<_PDF> &mix, const matrix_type &optimal_bandwidth, mixture<_PDF> &compressed_mix, int *mapping) {
    for (int i = 0; i < compressed_mix.size(); i++) {
      assert(mix.size() > 0);
      const size_t dims = mix.dims();

      //get the components of gm that are close to the j component of gmC
      data_type compWeight = 0.0;
      int numCloseComponents = 0;
      for (int j = 0; j < mix.size(); j++) {
        if (mapping[j] == i) {
          compWeight += mix.weight(j);
          numCloseComponents++;
        }
      }
      if (compWeight > 0.0) {
        mixture<_PDF> close_mix(dims);
        for (int j = 0, cCSM_idx = 0; j < mix.size(); j++) {
          //cCSM_idx the index of the "head" to write to closeComponentSubMixture
          if (mapping[j] == i) {
            close_mix.add(mix.component(j), mix.weight(j) / compWeight);
            cCSM_idx++;
          }
        }

        //mode 1 KDEcovariance
        //TODO the moment match should use covariance or kde_covariance based on it's own mixture PDF?
        //FIXME this only merges the sample model. We should also be merging the underlying model (if PDF is an explanation) how to make this generic???
        compressed_mix.replace(i, _PDF(close_mix, optimal_bandwidth), compWeight);
        //          compressed_mix.component(i).convolve(optimal_bandwidth);
        if (isnan(compressed_mix.component(i).covariance_log_determinant())) {
          std::cout << "NaN detected" << std::endl;
          //exit(0); //FIXME this should never occurr as there are other methods that prevent this to happen
        }
      } else {
        std::cout << "NO MAPPING!!!" << std::endl;
      }
    }
  }

  // compressed_mix must be initialized with the desired size
  template<typename _PDF>
  void goldberger_k_means(const mixture<_PDF> &mix, const matrix_type &optimal_bandwidth, mixture<_PDF> &compressed_mix,
                          int *mapping, int maxIterations = 20) {
    assert(mix.size() > 0);
    const size_t dims = mix.dims();


    /*// REMOVE
    if (mix.size() < 1) {
      std::cerr << "Assert failed in goldberger_k_means" << std::endl;
    }

    // REMOVE
    if (compressed_mix.size() != 2 && compressed_mix.size() != mix.size()) {
      std::cerr << "Assert failed in goldberger_k_means - compressed_mix size: " << compressed_mix.size() << std::endl;
      //throw std::exception();
    }*/

    //speedup cases
    if (mix.size() == 1) {
      compressed_mix.replace(0, mix.component(0), mix.weight(0));
      mapping[0] = 0;
      return;
    }

    if (mix.size() == 2) {
      compressed_mix.replace(0, mix.component(0), mix.weight(0));
      compressed_mix.replace(1, mix.component(1), mix.weight(1));
      mapping[0] = 0;
      mapping[1] = 1;
      return;
    }

    assert(mix.size() > 2);
    assert(compressed_mix.size() == 2);        //TODO for now we force a maximum goldberger-k-means of k=2

    _PDF moment_matched(mix), split_1(dims), split_2(dims);
    moment_matched.decompose_in_two(split_1, split_2);

    compressed_mix.replace(0, split_1, 0.5);
    compressed_mix.replace(1, split_2, 0.5);

    std::vector<std::vector<int>> visited_mapping;
    visited_mapping.push_back(std::vector<int>(mix.size()));
    for (size_t i = 0; i < mix.size(); i++)
      visited_mapping[0][i] = 0;

    int iteration = 1;
    bool stop = false;
    bool forcestop = false;
    while (!stop) {
      visited_mapping.push_back(std::vector<int>(mix.size()));        //new iteration

      regroup(mix, compressed_mix, mapping, visited_mapping[iteration], forcestop);

      refit(mix, optimal_bandwidth, compressed_mix, mapping);

      //test if the mapping has changed from the previous iteration
      stop = true;
      if (!forcestop) {
        for (int cmi = 0; cmi < mix.size(); cmi++) {
          stop = (visited_mapping[iteration - 1][cmi] == visited_mapping[iteration][cmi]);
          if (!stop)
            break;
        }
        iteration = iteration + 1;
        //      stop = ((distance < epsilon) || (std::abs(prevDistance - distance) / distance < epsilon));
        if (!stop) {
          bool cicle = false;
          for (int iti = 0; iti < iteration - 2; iti++) {
            for (int el = 0; el < mix.size(); el++) {           //if all equal cicle is equal then stop
              cicle = (visited_mapping[iteration - 1][el] == visited_mapping[iti][el]);
              if (!cicle)
                break;
            }
            if (cicle) {
              //std::cerr << "Avoiding a cicle...ufff.." << std::endl;
              break;
            }
          }
          stop = cicle;         //if cicle, stops, else continues
          if (iteration > maxIterations) {
            //std::cerr << "goldberger_k_means :: Converging too slow! MaxIterations exceeded. DEBUG!  - " << std::endl
            //    << "Stopping k-means but still going on with the okde" << std::endl;
          }
        }
      }
    }
  }

} // namespace xokdepp

#endif /* __XOKDEPP_GOLBERGER_K_MEANS_H__ */
