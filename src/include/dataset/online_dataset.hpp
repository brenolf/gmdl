#ifndef SRC_INCLUDE_ONLINE_DATASET_H_
#define SRC_INCLUDE_ONLINE_DATASET_H_

#include "dataset.hpp"
#include <iostream>
#include <boost/algorithm/string.hpp>

using namespace std;

static map<string, gmdl::Dataset::SampleType> TYPES = {
  {"<Training>", gmdl::Dataset::SampleType::Training},
  {"<Test>", gmdl::Dataset::SampleType::Test},
  {"<Correction>", gmdl::Dataset::SampleType::Correction}
};

namespace gmdl {
  class OnlineDataset: public Dataset {
  public:
    OnlineDataset(const vector<string> &classes): Dataset(classes) {}

    int get_dimension() {
      return -1;
    }

    bool next(Sample &sample, int *prediction, SampleType *type) {
      string line;
      cin >> line;

      if (cin.eof()) {
        return false;
      }

      (*type) = TYPES[line];
      cin >> line;

      if ((*type) == SampleType::Correction) {
        string aux;
        cin >> aux;
        (*prediction) = find_class_id(aux);
      }

      return read_sample(line, sample);
    }
  };
}

#endif
