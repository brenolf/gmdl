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

#ifndef SRC_INCLUDE_SAMPLE_H_
#define SRC_INCLUDE_SAMPLE_H_

#include <iostream>
#include "kde/basic_types.h"

namespace xokdepp {

  template<typename label_type>
  class sample {
  private:
    label_type _label;
    vector_type _features;

  public:
    sample(const label_type &label, const vector_type &features) :
        _label(label), _features(features) {
    }

    const label_type &label() const {
      return _label;
    }

    const vector_type &features() const {
      return _features;
    }

    size_t dims() const {
      return _features.rows();
    }

    friend std::ostream &operator<<(std::ostream &o, const sample<label_type> &sample) {
      for (int i = 0; i < sample._features.size(); i++)
        o << sample._features[i] << ",";
      o << sample._label;
      return o;
    }

  };

}

#endif /* SRC_INCLUDE_SAMPLE_H_ */
