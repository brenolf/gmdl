#include <iostream>
#include <fstream>
#include "json.hpp"
#include "dataset/dataset.h"
#include "mdc/mdc.h"

using namespace std;

int main() {
  std::ifstream config_file("./config.json");
  nlohmann::json config;
  config_file >> config;

  config_file.close();

  string current_set = config["READ_SET"].get<string>();
  vector<string> classes = config["SETS"][current_set].get<vector<string>>();

  mdc::Dataset d(config["DATASETS"].get<string>());
  d.open_set(current_set, classes);

  pair<vector<double>, int> sample;
  mdc::MDC classifier(d);

  double acc = 0;
  int i = 0;

  while (d.testing_samples(sample)) {
    i++;

    mdc::prediction p = classifier.predict(sample.first);

    if (p.label == sample.second) {
      acc++;
    } else {
      cout << "iteration: " << i << endl;

      cout << "data: ";

      for (auto d : sample.first) {
        cout << d << " ";
      }

      cout << endl << "DLs: ";

      for (auto d : p.description_lengths) {
        cout << d << " ";
      }

      cout << endl << "Probs: ";

      for (auto d : p.description_lengths) {
        cout << (1 - d) << " ";
      }

      cout << endl << "Theta: ";

      for (auto d : classifier.get_Theta()) {
        cout << d << " ";
      }

      double diff = abs(p.description_lengths[p.label] - p.description_lengths[sample.second]);

      cout << endl << "predicted: " << p.label << ", expected: " << sample.second << endl;
      cout << "DL diff: " << diff << endl << endl;
    }

    classifier.train(sample);
  }

  cout << endl << endl << "ACC: " << (acc / i) << endl;

  return 0;
}
