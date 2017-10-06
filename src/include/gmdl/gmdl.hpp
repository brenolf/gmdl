#ifndef SRC_INCLUDE_GMDL_H_
#define SRC_INCLUDE_GMDL_H_

#define GMDL_DEFAULT_OMEGA -32
#define GMDL_DEFAULT_BETA -32
#define GMDL_DEFAULT_SIGMA 1
#define GMDL_DEFAULT_F 1
#define GMDL_DEFAULT_ETA 0.01
#define GMDL_DEFAULT_ALPHA 0.9
#define GMDL_DEFAULT_TAU 1

#include <vector>
#include <map>
#include <math.h>
#include <random>

#include "dataset/dataset.hpp"
#include "kde/okde.h"
#include "kde/explanation.h"

using namespace std;

namespace gmdl {
  typedef struct _prediction_type {
    int label;
    vector<double> description_lengths;
  } prediction;

  class GMDL {
    typedef xokdepp::gaussian<xokdepp::vector_type> gaussian_type;
    typedef xokdepp::explanation<gaussian_type, xokdepp::vector_type> explanation_type;
    typedef xokdepp::oKDE<gaussian_type, xokdepp::vector_type> kde_type;

  private:
    int _classes;
    int _dimension;
    double _omega = pow(2, GMDL_DEFAULT_OMEGA);
    double _beta = pow(2, GMDL_DEFAULT_BETA);
    double _sigma = 1;
    double _forgeting_factor = 1;
    map<int, vector<kde_type>> _distributions;
    map<int, kde_type> _class_distributions;
    vector<double> _Theta;
    vector<double> _gradients;
    map<int, vector<double>> _means;
    map<int, vector<double>> _variances_acc;
    map<int, vector<long long>> _SAMPLES;
    mt19937 _gen;
    gmdl::Dataset &_initial_dataset;

    double _eta = GMDL_DEFAULT_ETA; // learning rate
    double _alpha = GMDL_DEFAULT_ALPHA; // momentum
    double _tau = GMDL_DEFAULT_TAU; // class prototype distance impact

    double _SEED = 123456789;
    double _MAX_THETA = 0.999999999;
    double _MIN_THETA = pow(2, -32);

  private:
    void __update_Theta(pair<vector<double>, int> &sample, int prediction) {
      double norm = __L(sample.first, prediction);
      int kronecker = (prediction == sample.second) ? 1 : 0;

      for (int i = 0; i < _dimension; i++) {
        double partial = __L_hat_attribute(sample.first, prediction, i);
        double grad = partial * log(_Theta.at(i)) * (1 - norm) * (kronecker - norm);

        _gradients.at(i) = _eta * -grad + _alpha * _gradients.at(i);

        _Theta.at(i) -= _gradients.at(i);

        _Theta.at(i) = min(_Theta.at(i), _MAX_THETA);
        _Theta.at(i) = max(_Theta.at(i), _MIN_THETA);
      }
    }

    double __get_distance_to_prototype(vector<double> &attributes, int class_index) {
      const kde_type &pdf = _class_distributions.at(class_index);
      const double SIZE = pdf.size();

      double variance = 0;

      if (SIZE == 0) {
        return INFINITY;
      }

      xokdepp::vector_type mu_bar(pdf.weight(0) * pdf.component(0).mean());

      for (int i = 1; i < SIZE; i++) {
        mu_bar += (pdf.weight(i) * pdf.component(i).mean());
      }

      for (int i = 0; i < SIZE; i++) {
        variance += pdf.weight(i) *
        (
          (pdf.component(i).mean() - mu_bar) *
          (pdf.component(i).mean() - mu_bar).transpose()
        )(0);
      }

      xokdepp::matrix_type S = xokdepp::matrix_type::Zero(_dimension, _dimension);
      S.diagonal() = pdf.weight(0) * pdf.component(0).covariance();

      for (int i = 1; i < SIZE; i++) {
        xokdepp::matrix_type S_prime = xokdepp::matrix_type::Zero(_dimension, _dimension);
        S_prime.diagonal() = pdf.weight(i) * pdf.component(i).covariance();

        S += S_prime;
      }

      S += (variance * xokdepp::matrix_type::Ones(_dimension, _dimension));

      if (S.determinant() == 0) {
        return INFINITY;
      }

      xokdepp::vector_type x(_dimension);

      for (int attr = 0; attr < _dimension; attr++) {
        x[attr] = attributes.at(attr);
      }

      return sqrt((x - mu_bar).transpose() * S.inverse() * (x - mu_bar));
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

    double __L_hat(vector<double> &attributes, int class_index, vector<double> &S) {
      double DL = 0;

      for (int attr = 0; attr < _dimension; attr++) {
        DL += __L_hat_attribute(attributes, class_index, attr);
      }

      return max(0.0, DL * S[class_index] + __LM(class_index));
    }

    double __L_total(vector<double> &attributes, vector<double> &distances) {
      double DL = 0;

      for (int c = 0; c < _classes; c++) {
        DL += __L_hat(attributes, c, distances);
      }

      return DL;
    }

    double __L(vector<double> &attributes, int class_index) {
      vector<double> S = get_distances(attributes);

      return __L_hat(attributes, class_index, S) / __L_total(attributes, S);
    }

    double __LM(int class_index) {
      double complexity = 0;

      for (int attr = 0; attr < _dimension; attr++) {
        complexity += _distributions.at(class_index).at(attr).size();
      }

      if (_tau != 0) {
        complexity += _class_distributions.at(class_index).size();
      }

      return complexity;
    }

    bool __isCovarianceDegenerate(long long samples, double mean, double variance, double new_sample) {

      if (samples == 1) {
        return false;
      }

      double delta = new_sample - mean;
      mean += delta / samples;

      double delta2 = new_sample - mean;

      double covariance = variance + ((delta * delta2) / samples);

      return covariance <= xokdepp::MIN_BANDWIDTH;
    }

    void __estimate_kernel_density(kde_type &pdf) {
      if (pdf.size() >= 3) {
        pdf.estimate_kernel_density();
      }
    }

    void __estimate_kernel_densities(int class_index) {
      if (_tau != 0) {
        __estimate_kernel_density(_class_distributions.at(class_index));
      }

      for (int attr = 0; attr < _dimension; attr++) {
        __estimate_kernel_density(_distributions.at(class_index).at(attr));
      }
    }

    void train(pair<vector<double>, int> &sample) {
      xokdepp::vector_type vectorized_class_sample(_dimension);
      normal_distribution<> noise_distribution(0, _sigma);

      for (int attr = 0; attr < _dimension; attr++) {
        long long &SAMPLES = _SAMPLES.at(sample.second).at(attr);
        double &MEAN = _means.at(sample.second).at(attr);
        double &VAR_ACC = _variances_acc.at(sample.second).at(attr);

        double value = sample.first.at(attr);

        while (__isCovarianceDegenerate(SAMPLES + 1, MEAN, VAR_ACC, value)) {
          value += noise_distribution(_gen);
        }

        SAMPLES++;

        double delta = value - MEAN;

        MEAN += delta / SAMPLES;

        double delta2 = value - MEAN;

        VAR_ACC += delta * delta2;

        vectorized_class_sample[attr] = value;
      }

      if (_tau != 0) {
        kde_type &class_pdf = _class_distributions.at(sample.second);
        class_pdf.add_sample(vectorized_class_sample);
      }

      for (int attr = 0; attr < _dimension; attr++) {
        xokdepp::vector_type vectorized_sample(1);
        vectorized_sample << vectorized_class_sample[attr];

        kde_type &pdf = _distributions.at(sample.second).at(attr);

        pdf.add_sample(vectorized_sample);
      }
    }

  public:
    GMDL(gmdl::Dataset &dataset) : _initial_dataset(dataset) {
      _classes = dataset.get_label_length();
      _dimension = dataset.get_dimension();
      _Theta = vector<double>(_dimension, _MAX_THETA);
      _gradients = vector<double>(_dimension, 0);
      _gen = mt19937(_SEED);

      for (int c = 0; c < _classes; c++) {
        vector<kde_type> attrs(_dimension, kde_type(1));

        _SAMPLES.insert(pair<int, vector<long long>>(c, vector<long long>(_dimension, 0)));
        _means.insert(pair<int, vector<double>>(c, vector<double>(_dimension, 0)));
        _variances_acc.insert(pair<int, vector<double>>(c, vector<double>(_dimension, 0)));

        _distributions.insert(pair<int, vector<kde_type>>(c, attrs));
        _class_distributions.insert(pair<int, kde_type>(c, kde_type(_dimension)));
      }
    }

    vector<double> get_distances(vector<double> &attributes) {
      vector<double> S(_classes, 1);

      if (_tau == 0) {
        return S;
      }

      S[0] = __get_distance_to_prototype(attributes, 0);

      double S_min = S[0];
      double S_max = S[0];

      for (int c = 1; c < _classes; c++) {
        S[c] = __get_distance_to_prototype(attributes, c);

        S_min = min(S_min, S[c]);
        S_max = max(S_max, S[c]);
      }

      for (int c = 0; c < _classes; c++) {
        double normalized = (S[c] - S_min) / (S_max - S_min);

        S[c] = isinf(S[c]) ? 1 : -log2(0.5 * (1 - normalized + _beta));
        S[c] = pow(S[c], _tau);
      }

      return S;
    }

    int get_classes_length() {
      return _classes;
    }

    vector<double> &get_Theta() {
      return _Theta;
    }

    void set_omega(double omega) {
      _omega = pow(2, -omega);
    }

    void set_beta(double beta) {
      _beta = pow(2, -beta);
    }

    void set_tau(double tau) {
      _tau = tau;
    }

    void set_sigma(double sigma) {
      _sigma = sigma;
    }

    void set_learning_rate(double eta) {
      _eta = eta;
    }

    void set_momentum(double alpha) {
      _alpha = alpha;
    }

    void set_forgeting_factor(double forgeting_factor) {
      _forgeting_factor = forgeting_factor;

      for (int c = 0; c < _classes; c++) {
        _class_distributions.at(c).set_forg(forgeting_factor);

        for (int attr = 0; attr < _dimension; attr++) {
          _distributions.at(c).at(attr).set_forg(forgeting_factor);
        }
      }
    }

    void train() {
      pair<vector<double>, int> sample;

      while (_initial_dataset.training_samples(sample)) {
        train(sample);
      }

      for (int c = 0; c < _classes; c++) {
        __estimate_kernel_densities(c);
      }
    }

    void train(pair<vector<double>, int> &sample, int prediction) {
      train(sample);
      __estimate_kernel_densities(sample.second);
      __update_Theta(sample, prediction);
    }

    prediction predict(vector<double> &attributes) {
      prediction p;

      for (int c = 0; c < _classes; c++) {
        p.description_lengths.push_back(__L(attributes, c));
      }

      double min = p.description_lengths[0];
      p.label = 0;

      for (int c = 1; c < _classes; c++) {
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
