#ifndef SRC_INCLUDE_INPUT_DATASET_H_
#define SRC_INCLUDE_INPUT_DATASET_H_

#include "dataset.hpp"
#include <iostream>

using namespace std;

namespace gmdl {
  class InputDataset: public Dataset {
  private:
    string replay = "";
    bool replayed = false;
    unsigned long int lines_read = 0;
    unsigned long int training_size = 0;
    bool trainingSetFinished = false;

    bool iterate(Sample &sample) {
      string line;

      if (replay != "" && !replayed) {
        line = replay;
        replayed = true;
        lines_read++;
      } else {
        cin >> line;

        if (cin.eof()) {
          return false;
        }

        lines_read++;
        replay = line;
      }

      return read_sample(line, sample);
    }

    bool training_samples(Sample &sample) {
      return (lines_read >= training_size) ? false : iterate(sample);
    }

    bool testing_samples(Sample &sample) {
      return (lines_read < training_size) ? false : iterate(sample);
    }

  public:
    InputDataset(const vector<string> &classes): Dataset(classes) {
      cin >> training_size;
    }

    int get_dimension() {
      Sample sample;
      
      replay = "";
      replayed = false;

      iterate(sample);

      lines_read--;

      return sample.first.size();
    }

    bool next(Sample &sample, int *prediction, SampleType *type) {
      (*prediction) = -1;
      (*type) = SampleType::Test;

      if (training_samples(sample)) {
        (*type) = SampleType::Training;
        return true;
      }

      if (!trainingSetFinished) {
        trainingSetFinished = true;
        return false;
      }

      return testing_samples(sample);
    }
  };
}

#endif
