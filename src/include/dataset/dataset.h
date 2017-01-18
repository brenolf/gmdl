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

#ifndef SRC_INCLUDE_DATASET_H_
#define SRC_INCLUDE_DATASET_H_

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

namespace xokdepp {

  template<typename label_type>
  class dataset {
  public:
    typedef sample<label_type> sample_type;

  protected:
    std::vector<sample_type> _samples;
    // class labels and counts
    std::map<label_type, data_type> _classes;

  public:
    dataset() {
    }

    virtual ~dataset() {
    }

    size_t number_of_classes() const {
      return _classes.size();
    }

    size_t size() const {
      return _samples.size();
    }

    const std::map<label_type, data_type> &classes() const {
      return _classes;
    }

    const std::vector<sample_type> &samples() const {
      return _samples;
    }

    void add_sample(const label_type &label, const vector_type &data) {
      _samples.push_back(sample_type(label, data));
    }

    void shuffle() {
      //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    	unsigned seed = 7;
      std::shuffle(this->_samples.begin(), this->_samples.end(), std::default_random_engine(seed));
      //---- alternative ----
      //std::random_device rd;
      //std::mt19937 g(rd());
      //std::shuffle(this->_samples.begin(), this->_samples.end(), g);
    }

    void generate_metadata() {
      _classes.clear(); // reset
      // FIXME: ensure that map entries are initialized to zero (not sure if operator initializes)
      for (sample_type sample : _samples)
        _classes[sample.label()]++;
    }

    void set(const dataset &dataset, size_t from, size_t to) {
      _samples.clear();
      _classes.clear();
      for (size_t i = from; i < to; i++)
        _samples.push_back(dataset._samples[i]);
      generate_metadata();
    }

    int dims() const {
      if (size())
        return _samples[0].dims();
      return 0;
    }

    friend std::ostream &operator<<(std::ostream &o, const dataset<label_type> &dataset) {
      for (sample_type sample : dataset._samples)
        o << sample << std::endl;
      return o;
    }

  };

}

#endif /* SRC_INCLUDE_DATASET_H_ */
