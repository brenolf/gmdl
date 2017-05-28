#ifndef INCLUDE_MDC_CLASSIFIER  
#define INCLUDE_MDC_CLASSIFIER  

#include "dataset/dataset.h"
#include "mdc/mdc.h"
#include <iostream>
#include <fstream>
#include <json.hpp>
#include <map>

using namespace std;

typedef struct {
  mdc::MDC *classifier;
  mdc::Dataset *dataset;
} ClassifierData;

typedef void (mdc::MDC::*FnPtr)(double);

static const map<string, FnPtr> METAPARAMS = {
  {"learning_rate", &mdc::MDC::set_learning_rate},
  {"momentum", &mdc::MDC::set_momentum},
  {"delta", &mdc::MDC::set_delta},
  {"beta", &mdc::MDC::set_beta},
  {"omega", &mdc::MDC::set_omega},
  {"forgetting_factor", &mdc::MDC::set_forgeting_factor},
  {"sigma", &mdc::MDC::set_sigma}
};

ClassifierData get_classifier_data(cmdline::parser *args) {
  std::ifstream config_file("./config.json");
  nlohmann::json config;
  config_file >> config;

  config_file.close();

  const string set = 
    args->exist("set") ? 
    args->get<string>("set") :
    config["set"].get<string>();

  const int label = 
    args->exist("label") ? 
    args->get<int>("label") :
    config["label"].get<int>();

  vector<string> classes = config["datasets"][set].get<vector<string>>();

  mdc::Dataset *d = new mdc::Dataset(config["datasets_path"].get<string>());
  d->set_label_column(label);
  d->open_set(set, classes);

  mdc::MDC *classifier = new mdc::MDC(*d);

  for (auto const &item : METAPARAMS) {
    if (args->exist(item.first)) {
      (classifier->*item.second)(args->get<double>(item.first));
    } else if (!args->exist("inline")) {
      try {
        (classifier->*item.second)(config[item.first].get<double>());
      } catch(...) {}
    }
  }

  return ClassifierData { classifier, d };
}

#endif
