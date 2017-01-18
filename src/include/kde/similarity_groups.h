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

#ifndef __XOKDEPP_SIMILARITY_GROUPS_H__
#define __XOKDEPP_SIMILARITY_GROUPS_H__

#include <assert.h>
#include <vector>
#include "kde/basic_parameters.h"
#include "kde/basic_types.h"
#include "kde/basic_functions.h"
#include "kde/mixture.h"
#include "kde/hellinger_distance.h"

namespace xokdepp {

  //returns the new compressed mixture
  //each entry corresponds to one group to be approximated by one gaussian
  template<typename _PDF1>
  std::vector<std::vector<size_t>> generate_similarity_groups_indexes(const mixture<_PDF1> &mix,
                                                                      const matrix_type &optimal_bandwidth, data_type Dth) {
    _PDF1 one_gauss_approx(mix, true);
    mixture<_PDF1> compressed(mix.dims());
    compressed.add(one_gauss_approx);
    double dist = hellinger_distance(compressed, mix);

    //TODO these two should belong to a class structure which independently maintained the state
    std::vector<double> helldists;
    helldists.push_back(dist);

    //each entry is a group, with a vector of indexes corresponding to the original components
    //FIXME In fact I should be using a map. This would be the most efficient way to do this
    std::vector<std::vector<size_t>> group_indexes;
    group_indexes.push_back(std::vector<size_t>());
    for (size_t i = 0; i < mix.size(); i++)
      group_indexes[0].push_back(i);

    bool stop = false;        //kinda useless but ok...
    while (!stop) {
      size_t max_idx;
      double max_elem = vector_max_val_and_idx(helldists, max_idx);
      if (max_elem <= Dth /*some other condition may be missing*/)
        break;

      mixture<_PDF1> sub_mix(mix.dims());
      for (size_t i = 0; i < group_indexes[max_idx].size(); i++)
        sub_mix.add(mix.component(group_indexes[max_idx][i]), mix.weight(group_indexes[max_idx][i]));

      data_type sub_mix_original_weight_sum = sub_mix.sum_of_weights();        //by JAIME - save original sub_mix weight sum
//      std::cout << "sub_mix weight sum:  " << sub_mix.sum_of_weights() << std::endl;
//      sub_mix.print_weights();
//      std::cout << "---" << std::endl;
      sub_mix.normalize_weights_preserve_relative_importance(); //by JAIME - mixs should always sum to one, normalize the mix
//      sub_mix.print_weights();
//      std::cout << "sub_mix weight sum:  " << sub_mix.sum_of_weights() << std::endl;

      int tmpMapping[sub_mix.size()];
      mixture<_PDF1> compressed_mix(sub_mix.dims(), std::min(2, (int)sub_mix.size()));
      goldberger_k_means(sub_mix, optimal_bandwidth, compressed_mix, tmpMapping);

//      std::cout << "k_means mapping" << std::endl;
//      for(size_t i = 0 ; i < sub_mix.size() ; i++)
//        std::cout << " " << tmpMapping[i] << " " ;
//      std::cout << std::endl;

      //since we are spliting in two - replace current entry, and add a new one
      //if compressed_mix size = 1 - then we only had one component to split so we only have one group (hence size = 1)
      if (compressed_mix.size() == 1) {
        //by JAIME - mix has been normalized, so correct the value back
        compressed.replace(max_idx, compressed_mix.component(0), sub_mix_original_weight_sum * compressed_mix.weight(0));
        std::vector<size_t> new_mappings;
        for (size_t i = 0; i < sub_mix.size(); i++) {         //sub_mix.size() = tmpMapping.size()
          if (tmpMapping[i] == 0)
            new_mappings.push_back(group_indexes[max_idx][i]);
        }
        group_indexes[max_idx].swap(new_mappings); //replace with new mappings
        //could call hellinger distance method, but since we only have one component distance should always be 0
        helldists[max_idx] = 0;
      }
      //in this case we can have had many components to split in two groups
      else if (compressed_mix.size() == 2) {
        compressed.replace(max_idx, compressed_mix.component(0), sub_mix_original_weight_sum * compressed_mix.weight(0)); //by JAIME - mix has been normalized, so correct the value back
        compressed.add(compressed_mix.component(1), sub_mix_original_weight_sum * compressed_mix.weight(1)); //by JAIME - mix has been normalized, so correct the value back
        std::vector<size_t> new_mappings_0, new_mappings_1;
        for (size_t i = 0; i < sub_mix.size(); i++) {             //sub_mix.size() (=) tmpMapping.size()
          if (tmpMapping[i] == 0)
            new_mappings_0.push_back(group_indexes[max_idx][i]);
          else if (tmpMapping[i] == 1)
            new_mappings_1.push_back(group_indexes[max_idx][i]);
          else {
            std::cerr << "mixture::generate_similarity_groups_indexes:: erroneous k-means Mapping" << std::endl;
            exit(0);
          }
        }
        group_indexes[max_idx].swap(new_mappings_0); //replace with new mappings
        group_indexes.push_back(new_mappings_1); //add the other one
        //update hellinger distances
        mixture<_PDF1> sb_mix_0(mix.dims()), sb_mix_1(mix.dims());
        for (size_t i = 0; i < group_indexes[max_idx].size(); i++) {
          sb_mix_0.add(mix.component(group_indexes[max_idx][i]), mix.weight(group_indexes[max_idx][i]));
        }
        int group_indexes_last_elem = group_indexes.size() - 1;
        for (size_t i = 0; i < group_indexes[group_indexes_last_elem].size(); i++) { //we know new one is at last element of group_indexes
          sb_mix_1.add(mix.component(group_indexes[group_indexes_last_elem][i]),
                       mix.weight(group_indexes[group_indexes_last_elem][i]));                //here it copies
        }
        //convert compressed components to mixtures
        mixture<_PDF1> cp_0(compressed.dims()), cp_1(compressed.dims());
        cp_0.add(compressed.component(max_idx), compressed.weight(max_idx));
        cp_1.add(compressed.component(group_indexes_last_elem), compressed.weight(group_indexes_last_elem));
        helldists[max_idx] = hellinger_distance(cp_0, sb_mix_0);
        helldists.push_back(hellinger_distance(cp_1, sb_mix_1));
      } else {
        std::cerr << "mixture::generate_similarity_groups_indexes:: we should be splitting in two components only" << std::endl;
        exit(0);
      }

//      std::cout << "computing group indexes" << std::endl;
//      for (size_t i = 0; i < group_indexes.size(); i++) {
//        std::cout << "[ ";
//        for (size_t j = 0; j < group_indexes[i].size(); j++) {
//          std::cout << " " << group_indexes[i][j] << " ";
//        }
//        std::cout << "] ; ";
//      }
//      std::cout << std::endl;
//      for (size_t i = 0; i < helldists.size(); i++) {
//        std::cout << " " << helldists[i] << " ";
//      }
//      std::cout << std::endl;

    }
    return group_indexes;
  }

} // namespace xokdepp

#endif /* __XOKDEPP_SIMILARITY_GROUPS_H__ */
