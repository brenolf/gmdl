#ifndef SRC_INCLUDE_FILE_DATASET_H_
#define SRC_INCLUDE_FILE_DATASET_H_

#include "dataset.hpp"
#include <fstream>
#include <boost/algorithm/string.hpp>

using namespace std;

namespace gmdl {
  class FileDataset: public Dataset {
  private:
    const string _root; 
    const string training;
    const string testing;
    bool trainingSetFinished = false;

    ifstream _training;
    ifstream _testing;

    void _close_sets() {
      if (_training.is_open()) {
        _training.close();
      }

      if (_testing.is_open()) {
        _testing.close();
      }
    }

    bool _iterate_over_file(Sample &sample, ifstream &file) {
      string line;

      if (!getline(file, line)) {
        return false;
      }

      vector<string> xs;
      vector<double> attributes;

      split(xs, line, boost::is_any_of(_separator));

      int LABEL_COLUMN = (_label_column == -1) ? xs.size() - 1 : _label_column;

      for (int column = 0; column < xs.size(); column++) {
        if (column == LABEL_COLUMN) {
          continue;
        }

        attributes.push_back(stod(xs[column]));
      }

      sample.first = attributes;
      sample.second = find_class_id(xs[LABEL_COLUMN]);

      return true;
    }

    void open_sets() {
      if (testing == "") {
        throw "no test set supplied";
      }

      _close_sets();

      if (training != "") {
        _training.open(_root + training);
      }

      _testing.open(_root + testing);
    }

    bool training_samples(Sample &sample) {
      return _iterate_over_file(sample, _training);
    }

    bool testing_samples(Sample &sample) {
      return _iterate_over_file(sample, _testing);
    }

  public:
    FileDataset(
      const string root, 
      const string training, 
      const string testing, 
      const vector<string> &classes
    ): Dataset(classes), _root(root), training(training), testing(testing) {
      open_sets();
    }

    ~FileDataset() {
      _close_sets();
    }

    int get_dimension() {
      Sample sample;
      _iterate_over_file(sample, _testing);

      _testing.clear();
      _testing.seekg(0);

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
