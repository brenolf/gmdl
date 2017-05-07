#!/bin/bash

OS=`uname -s | tr '[:upper:]' '[:lower:]'`;

wget https://github.com/nlohmann/json/releases/download/v2.1.1/json.hpp -O /tmp/json.hpp;

if [[ "$OS" == "darwin" ]]; then
  brew install llvm eigen boost;
  ln -s /usr/local/opt/llvm/bin/clang++ /usr/local/bin/clang-omp++
  mv /tmp/json.hpp /usr/local/include/json.hpp;
else
  apt-get install clang++ libomp-dev libboost-dev libeigen3-dev;
  ln -s /usr/bin/clang++ /usr/bin/clang-omp++;
  mv /tmp/json.hpp /usr/include/json.hpp;
fi

echo '\033[92mAll done! Run `make` to compile.\033[0m'
