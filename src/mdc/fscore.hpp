#ifndef INCLUDE_MDC_FSCORE
#define INCLUDE_MDC_FSCORE

#include <Eigen/Dense>
#include <iostream>

using namespace std;

void fscore(Eigen::MatrixXf &confusion) {
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

  const double recall_macro = recall.sum() / N;
  const double precision_macro = precision.sum() / N;

  const double F_macro = 2 * (precision_macro * recall_macro) / (precision_macro + recall_macro);

  cout << F_macro << endl;
}

#endif
