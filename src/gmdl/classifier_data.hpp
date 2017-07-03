#ifndef INCLUDE_GMDL_CLASSIFIER  
#define INCLUDE_GMDL_CLASSIFIER  

#include "dataset/dataset.hpp"
#include "gmdl/gmdl.hpp"
#include <iostream>
#include <fstream>
#include <json.hpp>
#include <map>
#include <boost/algorithm/string.hpp>

using namespace std;

typedef struct {
  gmdl::GMDL *classifier;
  gmdl::Dataset *dataset;
} ClassifierData;

typedef void (gmdl::GMDL::*FnPtr)(double);

static const map<string, FnPtr> METAPARAMS = {
  {"learning_rate", &gmdl::GMDL::set_learning_rate},
  {"momentum", &gmdl::GMDL::set_momentum},
  {"tau", &gmdl::GMDL::set_tau},
  {"beta", &gmdl::GMDL::set_beta},
  {"omega", &gmdl::GMDL::set_omega},
  {"forgetting_factor", &gmdl::GMDL::set_forgeting_factor},
  {"sigma", &gmdl::GMDL::set_sigma}
};

ClassifierData get_classifier_data(cmdline::parser *args) {
  std::ifstream config_file(args->get<string>("config"));
  nlohmann::json config;
  config_file >> config;

  config_file.close();

  const bool isInline = args->exist("inline");

  const string set = 
    args->exist("set") || isInline ? 
    args->get<string>("set") :
    config["set"].get<string>();

  const int label = 
    args->exist("label") || isInline ? 
    args->get<int>("label") :
    config["label"].get<int>();

  const string datasets_path = 
    args->exist("path") || isInline ? 
    args->get<string>("path") :
    config["datasets_path"].get<string>();

  const string training = 
    args->exist("training") || isInline ? 
    args->get<string>("training") :
    config["datasets"][set]["training"].get<string>();

  const string testing = 
    args->exist("testing") || isInline ? 
    args->get<string>("testing") :
    config["datasets"][set]["testing"].get<string>();

  vector<string> classes;

  if (args->exist("labels") || isInline) {
    split(classes, args->get<string>("labels"), boost::is_any_of(","));
  } else {
    classes = config["datasets"][set]["labels"].get<vector<string>>();
  }

  gmdl::Dataset *d = new gmdl::Dataset(datasets_path);
  d->set_label_column(label);
  d->open_sets(training, testing, classes);

  gmdl::GMDL *classifier = new gmdl::GMDL(*d);

  for (auto const &item : METAPARAMS) {
    if (args->exist(item.first)) {
      (classifier->*item.second)(args->get<double>(item.first));
    } else if (!isInline) {
      (classifier->*item.second)(config[item.first].get<double>());
    } else {
      throw "`inline` option is true but argument `" + item.first + "` is missing";
    }
  }

  return ClassifierData { classifier, d };
}

#endif
