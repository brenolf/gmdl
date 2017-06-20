#ifndef INCLUDE_MDC_FSCORE
#define INCLUDE_MDC_FSCORE

#include <Eigen/Dense>
#include <iostream>

using namespace std;

double harmonic_mean(double precision, double recall) {
  return 2 * ((precision * recall) / (precision + recall));
}

double macro_fscore(Eigen::MatrixXf &confusion) {
  const int N = confusion.rows();

  Eigen::RowVectorXd recall(N);
  Eigen::RowVectorXd precision(N);

  for (int i = 0; i < N; i++) {
    const double tp = confusion(i, i);
    const double fn = confusion.row(i).sum() - tp;
    const double fp = confusion.col(i).sum() - tp;

    if ((tp + fn) != 0) {
      recall(i) = tp / (tp + fn);
    } else {
      recall(i) = 0;
    }

    if ((tp + fp) != 0) {
      precision(i) = tp / (tp + fp);
    } else {
      precision(i) = 0;
    }
  }

  return harmonic_mean((precision.sum() / N), (recall.sum() / N));
}

double micro_fscore(Eigen::MatrixXf &confusion) {
  const int N = confusion.rows();

  double tp = 0;
  double fp = 0;
  double fn = 0;

  for (int i = 0; i < N; i++) {
    tp += confusion(i, i);
    fn += confusion.row(i).sum() - confusion(i, i);
    fp += confusion.col(i).sum() - confusion(i, i);
  }

  double precision = (tp + fp) == 0 ? 0 : tp / (tp + fp);
  double recall = (tp + fn) == 0 ? 0 : tp / (tp + fn);

  return harmonic_mean(precision, recall);
}

void fscore(Eigen::MatrixXf &confusion) {
  cout << macro_fscore(confusion) << "/" << micro_fscore(confusion) << endl;
}

#endif
