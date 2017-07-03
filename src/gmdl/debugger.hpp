#ifndef INCLUDE_GMDL_DEBUGGER
#define INCLUDE_GMDL_DEBUGGER

void debugger(int i, pair<vector<double>, int> &sample, gmdl::prediction &p, gmdl::GMDL *classifier) {
  cout << "iteration: " << i << endl;
  cout << "data: ";

  for (auto d : sample.first) {
    cout << d << " ";
  }

  cout << endl << "DLs: ";

  for (auto d : p.description_lengths) {
    cout << d << " ";
  }

  cout << endl << "Probs: ";

  for (auto d : p.description_lengths) {
    cout << (1 - d) << " ";
  }

  cout << endl << "Theta: ";

  for (auto d : classifier->get_Theta()) {
    cout << d << " ";
  }

  cout << endl << "S: ";

  for (auto d : classifier->get_distances(sample.first)) {
    cout << d << " ";
  }

  double diff = abs(p.description_lengths[p.label] - p.description_lengths[sample.second]);

  cout << endl << "predicted: " << p.label << ", expected: " << sample.second << endl;
  cout << "DL diff: " << diff << endl << endl;
}

#endif
