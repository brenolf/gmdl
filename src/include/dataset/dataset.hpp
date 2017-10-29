#ifndef SRC_INCLUDE_DATASET_H_
#define SRC_INCLUDE_DATASET_H_

#include <vector>
#include <map>
#include <iterator>

using namespace std;

namespace gmdl {
  typedef pair<vector<double>, int> Sample;

  class Dataset {
  protected:
    string _separator = ",";
    int _label_column = -1;
    vector<string> classes;
    map<string, int> class_ids;

  public:
    enum class SampleType {Training, Test};

    Dataset(const vector<string> &local_classes) {
      for (int i = 0; i < local_classes.size(); i++) {
        class_ids[local_classes[i]] = i;
        classes.push_back(local_classes[i]);
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

    string get_header() {
      string header = "";

      for (int i = 0; i < classes.size(); i++) {
        header += classes.at(i) + (i + 1 < classes.size() ? " " : "");
      }

      return header;
    }

    bool next(Sample &sample) {
      SampleType type;
      return next(sample, &type);
    }

    virtual int get_dimension() = 0;
    virtual bool next(Sample &sample, SampleType *type) = 0;
  };
}

#endif
