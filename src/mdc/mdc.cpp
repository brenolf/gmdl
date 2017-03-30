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

  double LEARNING_RATE = config["LEARNING_RATE"].get<double>();
  double MOMENTUM = config["MOMENTUM"].get<double>();
  double DELTA = config["DELTA"].get<double>();
  double BETA = config["BETA"].get<double>();
  double OMEGA = config["OMEGA"].get<double>();
  double FORGETTING_FACTOR = config["FORGETTING_FACTOR"].get<double>();
  double SIGMA = config["SIGMA"].get<double>();

  mdc::Dataset d(config["DATASETS"].get<string>());
  d.set_label_column(config["LABEL"].get<int>());
  d.open_set(current_set, classes);

  pair<vector<double>, int> sample;
  mdc::MDC classifier(d);

  classifier.set_learning_rate(LEARNING_RATE);
  classifier.set_momentum(MOMENTUM);
  classifier.set_delta(DELTA);
  classifier.set_beta(BETA);
  classifier.set_omega(OMEGA);
  classifier.set_forgeting_factor(FORGETTING_FACTOR);
  classifier.set_sigma(SIGMA);

  double acc = 0;
  int i = 0;

  xokdepp::matrix_type confusion = xokdepp::matrix_type::Zero(classes.size(), classes.size());

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

      cout << endl << "S: ";

      for (auto d : classifier.get_distances(sample.first)) {
        cout << d << " ";
      }

      double diff = abs(p.description_lengths[p.label] - p.description_lengths[sample.second]);

      cout << endl << "predicted: " << p.label << ", expected: " << sample.second << endl;
      cout << "DL diff: " << diff << endl << endl;
    }

    confusion(sample.second, p.label)++;

    classifier.train(sample, p.label);
  }

  cout << endl << endl << "ACC: " << (acc / i) << endl;
  cout << endl;

  cout << "M. CONFUSAO" << endl << confusion << endl;

  return 0;
}
