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

  mdc::Dataset d(config["DATASETS"].get<string>());

  d.open_set("iris", {"Iris-setosa", "Iris-versicolor", "Iris-virginica"});

  pair<vector<double>, int> sample;
  mdc::MDC classifier(d);

  double acc = 0;
  int TOTAL = 0;

  while (d.testing_samples(sample)) {
    TOTAL++;

    mdc::prediction p = classifier.predict(sample.first);

    if (p.label == sample.second) {
      acc++;
    } else {
      cout << "iteration: " << TOTAL << endl;

      cout << "data: ";

      for (auto d : sample.first) {
        cout << d << " ";
      }

      cout << endl << "DLs: ";

      for (auto d : p.description_lengths) {
        cout << d << " ";
      }

      cout << "predicted: " << p.label << ", expected: " << sample.second << endl << endl;
    }

    classifier.train(sample);
  }

  cout << endl << endl << "ACC: " << (acc / TOTAL) << endl;

  return 0;
}
