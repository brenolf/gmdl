#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <array>
#include <boost/algorithm/string.hpp>
#include <math.h>
#include "kde/okde.h"
#include "kde/explanation.h"

using namespace xokdepp;
using namespace std;
using namespace boost;

typedef gaussian<vector_type> gaussian_type;
typedef explanation<gaussian_type, vector_type> explanation_type;
typedef oKDE<gaussian_type, vector_type> kde_type;

map<int, array<kde_type, 4>> distributions;

pair<int, array<data_type, 3>> predict(vector<string> v) {
  array<data_type, 3> likelihood = {{ 0, 0, 0 }};

  for (int classe = 0; classe < 3; classe++) {
    for (int attr = 0; attr < 4; attr++) {
      vector_type sample(1);
      sample << stod(v.at(attr));

      data_type l = distributions.at(classe).at(attr).likelihood(sample);

      if (l == 0) {
        l = 0.0001;
      }

      likelihood[classe] += ceil(-log2(l));
    }
  }

  data_type min = INFINITY;
  int j = -1;

  for (int classe = 0; classe < 3; classe++) {
    if (likelihood[classe] < min) {
      min = likelihood[classe];
      j = classe;
    }
  }

  return pair<int, array<data_type, 3>>(j, likelihood);
}

int main() {
  ifstream training("/Users/breno/Dropbox/pesquisa/Mestrado - Breno/code/data/iris.train.data");
  ifstream testing("/Users/breno/Dropbox/pesquisa/Mestrado - Breno/code/data/iris.test.data");

  int label;
  string line;

  // initializes the distributions
  for (int classe = 0; classe < 3; classe++) {
    array<kde_type, 4> attrs = {
      { kde_type(1), kde_type(1), kde_type(1), kde_type(1) }
    };

    distributions.insert(pair<int, array<kde_type, 4>>(classe, attrs));
  }

  // reads training data
  while (getline(training, line)) {
    vector<string> x;
    split(x, line, is_any_of(","));

    label = (int) stod(x.at(4));

    for (int attr = 0; attr < 4; attr++) {
      vector_type sample(1);
      sample << stod(x.at(attr));

      distributions.at(label).at(attr).add_sample(sample);
    }
  }

  training.close();

  // estimates kde
  for (int classe = 0; classe < 3; classe++) {
    for (int attr = 0; attr < 4; attr++) {
      distributions.at(classe).at(attr).estimate_kernel_density();
    }
  }

  int TOTAL = 0;
  double acc = 0;

  // reads testing data
  while (getline(testing, line)) {
    TOTAL++;

    vector<string> x;
    split(x, line, is_any_of(","));

    label = (int) stod(x.at(4));

    pair<int, array<data_type, 3>> prediction = predict(x);

    int predicted = prediction.first;
    array<data_type, 3> DLs = prediction.second;

    if (predicted == label) {
      acc++;
    } else {
      cout << "iteration: " << TOTAL << endl;
      cout << "data: " << x.at(0) << " " << x.at(1) << " " << x.at(2) << " " << x.at(3) << endl;
      cout << "DLs: " << DLs.at(0) << " " << DLs.at(1) << " " << DLs.at(2) << endl;
      cout << "predicted: " << predicted << ", expected: " << label << endl << endl;
    }

    for (int attr = 0; attr < 4; attr++) {
      vector_type sample(1);
      sample << stod(x.at(attr));

      distributions.at(label).at(attr).add_sample(sample);
      distributions.at(label).at(attr).estimate_kernel_density();
    }
  }

  testing.close();

  cout << endl << endl << "ACC: " << (acc / TOTAL) << endl;

  return 0;
}
