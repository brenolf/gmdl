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

#ifndef SRC_INCLUDE_EVAL_DATASET_H_
#define SRC_INCLUDE_EVAL_DATASET_H_

#include <vector>
#include <string>
#include <algorithm>    //std::find

#include "kde/basic_types.h"
#include "dataset/dataset.h"

namespace xokdepp {

  template<typename label_type>
  class split_dataset: public dataset<label_type> {
    typedef dataset<label_type> dataset_type;

  private:
    dataset_type _train;
    dataset_type _test;

    double _train_set_fraction;

  public:
    split_dataset(double train_set_fraction) : _train_set_fraction(train_set_fraction) {
      //shuffle();
    }

    void shuffle() {
      dataset_type::shuffle();
//      int train_set_last_idx = this->_samples.size() * (_train_set_fraction / 100.0);
//      _train.set(*this, 0, train_set_last_idx);
//      _test.set(*this, train_set_last_idx, this->_samples.size())
      reset_train_test_partitions();
    }

    void reset_train_test_partitions(){
    	int train_set_last_idx = round( this->_samples.size() * (_train_set_fraction / 100.0));
    	//std::cout << "Train samples: " << train_set_last_idx << std::endl;
    	_train.set(*this, 0, train_set_last_idx);
    	_test.set(*this, train_set_last_idx, this->_samples.size());
    }

    const dataset_type &get_train_dataset() const {
      return _train;
    }

    const dataset_type &get_test_dataset() const {
      return _test;
    }

    double train_size() {
      return _train.size();
    }

    double test_size() {
      return _test.size();
    }

    size_t dims() const {
    	return _train.dims();
    }

    friend std::ostream &operator<<(std::ostream &o, const split_dataset<label_type> &dataset) {
      o << "Train set (size = " << dataset._train.size() << ")" << std::endl;
      for (auto sample : dataset._train.samples()) {
        o << sample << std::endl;
      }
      o << "Test set (size = " << dataset._test.size() << ")" << std::endl;
      for (auto sample : dataset._test.samples()) {
        o << sample << std::endl;
      }
      return o;
    }

  };

}

#endif /* SRC_INCLUDE_EVAL_DATASET_H_ */
