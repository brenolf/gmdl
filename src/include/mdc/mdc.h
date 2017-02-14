#ifndef SRC_INCLUDE_MDC_H_
#define SRC_INCLUDE_MDC_H_

#include <vector>
#include <map>
#include <math.h>

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
    double _omega = 0.0001;
    double _forgeting_factor = 1;
    map<int, vector<kde_type>> _distributions;

  public:
    MDC(mdc::Dataset &dataset) {
      _classes = dataset.get_label_length();
      _dimension = dataset.get_dimension();

      for (int c = 0; c < _classes; c++) {
        vector<kde_type> attrs = {
          { kde_type(1), kde_type(1), kde_type(1), kde_type(1) }
        };

        _distributions.insert(pair<int, vector<kde_type>>(c, attrs));
      }

      pair<vector<double>, int> sample;

      while (dataset.training_samples(sample)) {
        train(sample);
      }
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
    }

    prediction predict(vector<double> &attributes) {
      prediction p;

      p.description_lengths = {{ 0, 0, 0 }};

       for (int c = 0; c < _classes; c++) {
         for (int attr = 0; attr < _dimension; attr++) {
           xokdepp::vector_type sample(1);
           sample << attributes.at(attr);

           kde_type &pdf = _distributions.at(c).at(attr);

           double l = (double) pdf.likelihood(sample);

           if (l == 0) {
             l = _omega;
           }

           p.description_lengths[c] += ceil(-log2(l));
         }
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
