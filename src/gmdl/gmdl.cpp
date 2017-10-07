#include "./args.hpp"
#include "./classifier_data.hpp"
#include "./debugger.hpp"
#include <boost/algorithm/string/join.hpp>

int main(int argc, char *argv[]) {
  cmdline::parser *args = get_parser(argc, argv);
  ClassifierData data = get_classifier_data(args);

  data.classifier->train();

  const int length = data.classifier->get_classes_length();
  double acc = 0;
  int i = 0;

  Eigen::MatrixXf confusion = Eigen::MatrixXf::Zero(length, length);

  pair<vector<double>, int> sample;
  vector<string> predicted_classes;

  while (data.dataset->testing_samples(sample)) {
    i++;

    gmdl::prediction p = data.classifier->predict(sample.first);

    if (args->exist("online")) {
      cout << data.dataset->get_label_name(p.label) << endl;
    } else {
      predicted_classes.push_back(data.dataset->get_label_name(p.label));
    }

    if (p.label == sample.second) {
      acc++;
    } else {
      if (args->exist("online")) {
        data.classifier->train(sample, p.label);
      }
      
      if (!args->exist("quiet")) {
        debugger(i, sample, p, data.classifier);
      }
    }

    confusion(sample.second, p.label)++;
  }

  if (!args->exist("online")) {
    cout << boost::algorithm::join(predicted_classes, ", ") << endl;
  }

  if (args->exist("confusion-matrix")) {
    cout << data.dataset->get_header() << endl;
    cout << confusion << endl;
  }

  return 0;
}
