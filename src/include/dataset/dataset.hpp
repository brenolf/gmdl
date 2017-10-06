#ifndef SRC_INCLUDE_DATASET_H_
#define SRC_INCLUDE_DATASET_H_

#include <vector>
#include <map>
#include <iterator>

using namespace std;

namespace gmdl {
  class Dataset {
  protected:
    string _separator = ",";
    int _label_column = -1;
    const vector<string> &classes;
    map<string, int> class_ids;

  public:
    Dataset(const vector<string> &classes): classes(classes) {
      for (int i = 0; i < classes.size(); i++) {
        class_ids[classes[i]] = i;
      }
    }

    void set_separator(string separator) {
      _separator = separator;
    }

    void set_label_column(int label_column) {
      _label_column = label_column;
    }

    int get_label_length() {
      return classes.size();
    }

    int find_class_id(const string name) {
      return class_ids[name];
    }

    string get_label_name(int i) {
      return classes.at(i);
    }

    virtual int get_dimension() = 0;
    virtual bool training_samples(pair<vector<double>, int> &sample) = 0;
    virtual bool testing_samples(pair<vector<double>, int> &sample) = 0;
  };
}

#endif
