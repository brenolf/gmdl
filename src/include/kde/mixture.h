// (C) Copyright 2014-2016, Jaime Ferreira, David Martins de Matos, Ricardo Ribeiro
// Spoken Language Systems Lab, INESC ID, IST/Universidade de Lisboa

// This file is part of XOKDE++.

// XOKDE++ is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.

// XOKDE++ is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

#ifndef __XOKDEPP_MIXTURE_H__
#define __XOKDEPP_MIXTURE_H__

#include <assert.h>
#include <vector>
#include <time.h>       /* time */
#include "kde/basic_parameters.h"
#include "kde/basic_types.h"
#include "kde/basic_functions.h"

namespace xokdepp {

  template<typename _PDF>
  class mixture {
  public:
    size_t dims() const {
      return _dims;
    }

  protected:
    const size_t _dims;

    std::vector<_PDF> _components;

    std::vector<data_type> _weights;

  public:
    inline mixture(size_t dims) :
        _dims(dims) {
      //EMPTY
    }

    inline mixture(size_t dims, size_t size) :
        _dims(dims), _weights(size) {
      for (int i = 0; i < size; i++) {
        _components.push_back(_PDF(dims));
      }
    }

    inline mixture(mixture &&other) :
        _dims(other.dims()) {
      addAll(other);
    }

    inline mixture(const mixture &other) :
        _dims(other.dims()) {
      addAll(other);
    }

    // whitening constructor
    inline mixture(size_t dims, const mixture &to_whiten, const _PDF &smp, const matrix_type &F) :
        mixture(dims) {
      for (size_t i = 0; i < to_whiten.size(); i++)
        add(_PDF(F.transpose() * (to_whiten.mean(i) - smp.mean()), FtCF(to_whiten.component(i).base_covariance(), F)),
            to_whiten.weight(i));
    }

    virtual ~mixture() {
    }

    void reset() {
      _components.clear();
      _weights.clear();
    }

    void add(const _PDF &pdf, data_type weight = 1) {
      assert(_dims == pdf.dims());
      _components.push_back(pdf);
      _weights.push_back(weight);
    }

    void add(_PDF &&pdf, data_type weight = 1) {
      assert(_dims == pdf.dims());
      _components.push_back(pdf);
      _weights.push_back(weight);
    }

    void add(const std::vector<_PDF> &pdfs) {
      for (auto pdf : pdfs) {
        _components.push_back(pdf);
        _weights.push_back(1); // force weight 1
      }
    }

    void add(std::vector<_PDF> &&pdfs) {
      for (auto pdf : pdfs) {
        _components.push_back(pdf);
        _weights.push_back(1); // force weight 1
      }
    }

    void addAll(mixture &&other) {
      for (size_t cx = 0; cx < other.size(); cx++)
        add(other.component(cx), other.weight(cx));
    }

    void addAll(const mixture &other) {
      for (size_t cx = 0; cx < other.size(); cx++)
        add(other.component(cx), other.weight(cx));
    }

    //TODO could be prettier
    void replace(size_t index, _PDF &&pdf, data_type weight = 1) {
      assert(index < _components.size());
      _components[index] = pdf;
      _weights[index] = weight;
    }

    void replace(size_t index, const _PDF &pdf, data_type weight = 1) {
      assert(index < _components.size());
      _components[index] = pdf;
      _weights[index] = weight;
    }

    void operator=(mixture<_PDF> &&mixture) {
      _components = mixture._components;
      _weights = mixture._weights;
    }

    void operator=(const mixture<_PDF> &mixture) {
      _components = mixture._components;
      _weights = mixture._weights;
    }

    size_t size() const {
      return _components.size();
    }

    const _PDF &component(size_t index) const {
      return _components[index];
    }

    _PDF &component(size_t index) {
      return _components[index];
    }

    //usefull for compression on okde
    const std::vector<_PDF> &components() const {
      return _components;
    }

    std::vector<_PDF> &components() {
      return _components;
    }

    const std::vector<data_type> &weights() const {
      return _weights;
    }

    std::vector<data_type> &weights() {
      return _weights;
    }

    const vector_type &mean(size_t index) const {
      return _components[index].mean();
    }

//    const matrix_type &covariance(size_t index) const {
//      return _components[index].covariance();
//    }

    data_type weight(size_t index) const {
      return _weights[index];
    }

    void set_weight(size_t index, data_type new_weight) {
      _weights[index] = new_weight;
    }

    data_type sum_of_weights() const {
      data_type total_sum = 0;
//      for (data_type w : _weights)
//        total_sum += w;
      for (int i = 0; i < _weights.size(); i++)
        total_sum += _weights[i];
      return total_sum;
    }

    //old approach - preserves relative importance of each weight but normalizes such that it sums to 1
    void normalize_weights_preserve_relative_importance() {
      data_type weights_sum = sum_of_weights();
      for (int i = 0; i < _weights.size(); i++) {
        _weights[i] /= weights_sum;
      }
    }

    //new approach, compatible to oKDE implementation.
    void normalize_weights() {
      for (int i = 0; i < _weights.size(); i++) {
        _weights[i] = 1.0 / _weights.size();
      }
    }

    void scale_weights(data_type factor) {
      for (size_t cx = 0; cx < size(); cx++)
        _weights[cx] *= factor;
    }

    void print_weights() const {
      for (double w : _weights)
        std::cout << w << std::endl;
    }

  private:
    //should be somewhere else!!!!!
    data_type logAdd(data_type logX, data_type logY) const {
      // 1. make X the max
      if (logY > logX) {
        data_type temp = logX;
        logX = logY;
        logY = temp;
      }
      // 2. now X is bigger
      if (logX == -std::numeric_limits<data_type>::infinity()) {
        return logX;
      }
      // 3. how far "down" (think decibels) is logY from logX?
      //    if it's really small (20 orders of magnitude smaller), then ignore
      data_type negDiff = logY - logX;
      if (negDiff < -20) {
        return logX;
      }
      // 4. otherwise use some nice algebra to stay in the log domain
      //    (except for negDiff)
      return logX + log(1.0 + exp(negDiff));
    }

  public:
    data_type likelihood(const vector_type &sample) const {
      data_type likelihood = 0.0;
      for (size_t i = 0; i < size(); i++) {
        likelihood += _weights[i] * exp(_components[i].log_likelihood(sample));
//        std::cout << "------------" << std::endl;
//        std::cout << "likelihood: " << likelihood << std::endl;
//        std::cout << "exp(log_likelihood):" << exp(_components[i].log_likelihood(sample)) << std::endl;
//        std::cout << "log_likelihood:" << _components[i].log_likelihood(sample) << std::endl;
//        std::cout << "------------" << std::endl;
      }
      return likelihood;
    }

    data_type likelihood_debug(const vector_type &sample) const {
          data_type likelihood = 0.0;
          for (size_t i = 0; i < size(); i++) {
            likelihood += _weights[i] * exp(_components[i].log_likelihood(sample));
            std::cout << "------------" << std::endl;
            std::cout << "likelihood: " << likelihood << std::endl;
            std::cout << "exp(log_likelihood):" << exp(_components[i].log_likelihood(sample)) << std::endl;
            std::cout << "log_likelihood:" << _components[i].log_likelihood(sample) << std::endl;
            std::cout << "------------" << std::endl;
          }
          return likelihood;
        }

    //convolve the optimalBandwidth H with all components from the mixture
    template<typename COVARIANCE_TYPE>
    void convolve(const COVARIANCE_TYPE &bandwidth) {
//#ifdef OKDE_PARALLELIZATION
//#pragma omp parallel for
//#endif
      for (size_t cx = 0; cx < _components.size(); cx++)
        _components[cx].convolve(bandwidth);
    }

    template<typename COVARIANCE_TYPE>
    void deconvolve(const COVARIANCE_TYPE &bandwidth) {
      for (size_t cx = 0; cx < _components.size(); cx++)
        _components[cx].deconvolve(bandwidth);
    }

    vector_type compute_sample_probability_vector(const vector_type &sample) {
      vector_type ans(this->size());
      for (int i = 0; i < this->size(); i++) {
        ans[i] = _components[i].likelihood(sample);
      }
      return ans;
    }

  };

  template<typename _PDF>
  std::ostream &operator<<(std::ostream &o, const mixture<_PDF> &mix) {
    o << std::endl << "Dims = " << mix.dims() << std::endl;
    for (int i = 0; i < mix.size(); i++) {
      o << "Mix[" << i << "]" << std::endl;
      o << "weight: " << std::setprecision(9) << mix.weight(i) << std::endl;
      o << mix.component(i);
    }
    return o;
  }

  template<typename _PDF>
  std::ostream &operator<<(std::ostream &o, mixture<_PDF> &&mix) {
    o << std::endl << "ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ" << "Dims = " << mix.dims() << std::endl;
    for (int i = 0; i < mix.size(); i++) {
      o << "Mix[" << i << "]" << std::endl;
      o << "weight: " << mix.weight(i) << std::endl;
      o << mix.component(i);
    }
    return o;
  }

} // namespace xokdepp

#endif /* __XOKDEPP_MIXTURE_H__ */
