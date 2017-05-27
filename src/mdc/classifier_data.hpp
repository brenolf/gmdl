#ifndef INCLUDE_MDC_CLASSIFIER  
#define INCLUDE_MDC_CLASSIFIER  

#include "dataset/dataset.h"
#include "mdc/mdc.h"
#include <iostream>
#include <fstream>
#include <json.hpp>

using namespace std;

typedef struct {
  mdc::MDC *classifier;
  mdc::Dataset *dataset;
} ClassifierData;

ClassifierData get_classifier_data() {
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

  mdc::Dataset *d = new mdc::Dataset(config["DATASETS"].get<string>());
  d->set_label_column(config["LABEL"].get<int>());
  d->open_set(current_set, classes);

  mdc::MDC *classifier = new mdc::MDC(*d);

  classifier->set_learning_rate(LEARNING_RATE);
  classifier->set_momentum(MOMENTUM);
  classifier->set_delta(DELTA);
  classifier->set_beta(BETA);
  classifier->set_omega(OMEGA);
  classifier->set_forgeting_factor(FORGETTING_FACTOR);
  classifier->set_sigma(SIGMA);

  return ClassifierData { classifier, d };
}

#endif
