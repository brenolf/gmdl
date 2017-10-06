#include "./args.hpp"
#include "./classifier_data.hpp"
#include "./debugger.hpp"

int main(int argc, char *argv[]) {
  cmdline::parser *args = get_parser(argc, argv);
  ClassifierData data = get_classifier_data(args);

  data.classifier->train();

  const int length = data.classifier->get_classes_length();
  double acc = 0;
  int i = 0;

  Eigen::MatrixXf confusion = Eigen::MatrixXf::Zero(length, length);

  pair<vector<double>, int> sample;

  while (data.dataset->testing_samples(sample)) {
    i++;

    gmdl::prediction p = data.classifier->predict(sample.first);

    if (p.label == sample.second) {
      acc++;
    } else {
      if (!args->exist("no-incremental-learning")) {
        data.classifier->train(sample, p.label);
      }
      
      if (!args->exist("quiet")) {
        debugger(i, sample, p, data.classifier);
      }
    }

    confusion(sample.second, p.label)++;
  }

  cout << confusion << endl;

  return 0;
}
