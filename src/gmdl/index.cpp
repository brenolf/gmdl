#include "./args.hpp"
#include "./classifier_data.hpp"
#include "./debugger.hpp"
#include <boost/algorithm/string/join.hpp>

int main(int argc, char *argv[]) {
  cmdline::parser *args = get_parser(argc, argv);
  ClassifierData data = get_classifier_data(args);

  if (!args->exist("online")) {
    data.classifier->train(data.dataset);
  }

  const int length = data.classifier->get_classes_length();

  Eigen::MatrixXf confusion = Eigen::MatrixXf::Zero(length, length);

  gmdl::Sample sample;
  vector<string> predicted_classes;
  gmdl::Dataset::SampleType sample_type;

  int index = 0;

  while (data.dataset->next(sample, &sample_type)) {
    index++;

    switch (sample_type) {
      case gmdl::Dataset::SampleType::Training:
        if (args->exist("online")) {
          data.classifier->train(sample);
        }
      break;

      default:
        gmdl::prediction p = data.classifier->predict(sample.first);
      
        if (args->exist("online")) {
          cout << data.dataset->get_label_name(p.label) << endl;
        } else {
          predicted_classes.push_back(data.dataset->get_label_name(p.label));
        }
    
        if (p.label != sample.second && !args->exist("quiet")) {
          debugger(index, sample, p, data.classifier);
        }
    
        confusion(sample.second, p.label)++;
    }
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
