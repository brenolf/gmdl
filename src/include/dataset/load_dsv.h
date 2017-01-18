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

#ifndef SRC_INCLUDE_LOAD_DSV_H_
#define SRC_INCLUDE_LOAD_DSV_H_

#include <tuple>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

#include "kde/basic_types.h"
#include "dataset/sample.h"
#include "dataset/dataset.h"

namespace xokdepp {

  template<typename label_type, typename dataset_type>
  dataset_type &&load_dsv(const char *file_path, char delimiter, dataset_type &&dataset, size_t dims) {
    std::cout << "loading database" << std::endl;
    std::ifstream f(file_path);

    std::map<std::string, int> sample_classes;
    int class_number = 0; // class number (serial)

    std::string s;
    while (getline(f, s)) {
      std::istringstream iss(s);
      if (s.empty())
        continue;

      int intra_sample_counter = 0;
      vector_type sample_features(dims);

      for (unsigned long i = 0; i < dims; i++) {
        getline(iss, s, delimiter);
        data_type fieldvalue = 0.0;
        std::istringstream(s) >> fieldvalue;
        sample_features[intra_sample_counter] = fieldvalue;
        intra_sample_counter++;
      }

      getline(iss, s, delimiter);
      auto position = sample_classes.find(s);
      if (position == sample_classes.end()) {
        // class has never been seen
        sample_classes[s] = class_number++;
      }
      dataset.add_sample(label_type(s, sample_classes[s]), sample_features);
    }

    dataset.reset_train_test_partitions();//all samples added, lets set the train/test partitions

    dataset.generate_metadata();

    std::cout << "database loaded;" << std::endl;
    std::cout << "number of classes: " << dataset.number_of_classes() << std::endl;
    std::cout << "number of samples: " << dataset.size() << std::endl;
    std::cout << "Classes label \t numerical_label \t count" << std::endl;

    // for (auto cl : dataset.classes()) {
      // std::cout << "\t";
      // std::cout << cl.first;
      // std::cout << "\t";
      // std::cout << cl.second;
      // std::cout << std::endl;
    // }

    return std::move(dataset);
  }

}

#endif /* SRC_INCLUDE_LOAD_DSV_H_ */
