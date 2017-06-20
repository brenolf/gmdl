#ifndef SRC_INCLUDE_DATASET_H_
#define SRC_INCLUDE_DATASET_H_

#include <fstream>
#include <map>
#include <vector>
#include <boost/algorithm/string.hpp>

using namespace std;

namespace mdc {
  class Dataset {
  private:
    const string _root;

    string _separator = ",";
    int _label_column = -1;
    int _classes_length = -1;

    map<string, int> _classes;
    map<int, string> _classes_lookup;

    ifstream _training;
    ifstream _testing;

    void _close_sets() {
      if (_training.is_open()) {
        _training.close();
      }

      if (_testing.is_open()) {
        _testing.close();
      }

      _classes.clear();
      _classes_lookup.clear();
    }

    bool _iterate_over_file(pair<vector<double>, int> &sample, ifstream &file) {
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
      sample.second = _classes.at(xs[LABEL_COLUMN]);

      return true;
    }

  public:
    Dataset(const string root) : _root(root) {
    }

    ~Dataset() {
      _close_sets();
    }

    void set_separator(string separator) {
      _separator = separator;
    }

    void set_label_column(int label_column) {
      _label_column = label_column;
    }

    string get_label_name(int i) {
      return _classes_lookup.at(i);
    }

    int get_dimension() {
      pair<vector<double>, int> sample;
      _iterate_over_file(sample, _testing);

      _testing.clear();
      _testing.seekg(0);

      return sample.first.size();
    }

    int get_label_length() {
      return _classes_length;
    }

    void open_sets(const string training, const string testing, const vector<string> &classes) {
      if (testing == "") {
        throw "no test set supplied";
      }

      _close_sets();

      _classes_length = classes.size();

      for (int i = 0; i < _classes_length; i++) {
        _classes[classes[i]] = i;
        _classes_lookup[i] = classes[i];
      }

      if (training != "") {
        _training.open(_root + training);
      }

      _testing.open(_root + testing);
    }

    bool training_samples(pair<vector<double>, int> &sample) {
      return _iterate_over_file(sample, _training);
    }

    bool testing_samples(pair<vector<double>, int> &sample) {
      return _iterate_over_file(sample, _testing);
    }
  };
}

#endif
