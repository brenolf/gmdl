#ifndef SRC_INCLUDE_MDC_H_
#define SRC_INCLUDE_MDC_H_

#include <vector>
#include <map>
#include <math.h>
#include <random>

#include "kde/okde.h"
#include "kde/explanation.h"

using namespace std;

namespace mdc {
  typedef struct _prediction_type {
    int label;
    vector<double> description_lengths;
  } prediction;

  class MDC {
    typedef xokdepp::gaussian<xokdepp::vector_type> gaussian_type;
    typedef xokdepp::explanation<gaussian_type, xokdepp::vector_type> explanation_type;
    typedef xokdepp::oKDE<gaussian_type, xokdepp::vector_type> kde_type;

  private:
    int _classes;
    int _dimension;
    double _omega = pow(2, -32);
    double _forgeting_factor = 1;
    map<int, vector<kde_type>> _distributions;
    vector<double> _Theta;
    double _alpha = 0.1;

  private:
    void __update_Theta(pair<vector<double>, int> &sample) {
      double norm = __L_bar(sample.first, sample.second);

      for (int i = 0; i < _dimension; i++) {
        double partial = __L_hat_attribute(sample.first, sample.second, i);

        double J = partial * log(_Theta.at(i)) * (1 - norm);

        _Theta.at(i) -= _alpha * -J;
      }
    }

    double __L_hat_attribute(vector<double> &attributes, int class_index, int attr) {
      xokdepp::vector_type sample(1);
      sample << attributes.at(attr);

      kde_type &pdf = _distributions.at(class_index).at(attr);

      double density = (double) pdf.likelihood(sample);

      if (density == 0 || isinf(density)) {
        density = _omega;
      }

      density = min(density, 1.0);

      return pow(ceil(-log2(density)), _Theta.at(attr));
    }

    double __L_hat(vector<double> &attributes, int class_index) {
      double DL = 0;

      for (int attr = 0; attr < _dimension; attr++) {
        DL += __L_hat_attribute(attributes, class_index, attr);
      }

      return max(DL, 0.0);
    }

    double __L_total(vector<double> &attributes) {
      double DL = 0;

      for (int c = 0; c < _classes; c++) {
        DL += __L_hat(attributes, c);
      }

      return DL;
    }

    double __L_bar(vector<double> &attributes, int class_index) {
      return __L_hat(attributes, class_index) / __L_total(attributes);
    }

  public:
    MDC(mdc::Dataset &dataset) {
      _classes = dataset.get_label_length();
      _dimension = dataset.get_dimension();
      _Theta = vector<double>(_dimension, 0.9999999999);

      for (int c = 0; c < _classes; c++) {
        vector<kde_type> attrs(_dimension, kde_type(1));

        _distributions.insert(pair<int, vector<kde_type>>(c, attrs));
      }

      pair<vector<double>, int> sample;

      while (dataset.training_samples(sample)) {
        train(sample);
      }
    }

    vector<double> &get_Theta() {
      return _Theta;
    }

    void set_omega(double omega) {
      _omega = omega;
    }

    double get_omega() {
      return _omega;
    }

    void set_forgeting_factor(double forgeting_factor) {
      _forgeting_factor = forgeting_factor;

      for (int c = 0; c < _classes; c++) {
        for (int attr = 0; attr < _dimension; attr++) {
          _distributions.at(c).at(attr).set_forg(forgeting_factor);
        }
      }
    }

    double get_forgeting_factor() {
      return _forgeting_factor;
    }

    void train(pair<vector<double>, int> &sample) {
      for (int attr = 0; attr < _dimension; attr++) {
        xokdepp::vector_type vectorized_sample(1);
        vectorized_sample << sample.first.at(attr);

        kde_type &pdf = _distributions.at(sample.second).at(attr);

        pdf.add_sample(vectorized_sample);

        if (pdf.size() >= 3) {
          pdf.estimate_kernel_density();
        }
      }

      __update_Theta(sample);
    }

    prediction predict(vector<double> &attributes) {
      prediction p;

      for (int c = 0; c < _classes; c++) {
        p.description_lengths.push_back(__L_bar(attributes, c));
      }

      double min = INFINITY;
      p.label = -1;

      for (int c = 0; c < _classes; c++) {
        if (p.description_lengths[c] < min) {
          min = p.description_lengths[c];
          p.label = c;
        }
      }

      return p;
    }
  };
}

#endif
