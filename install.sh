#!/bin/bash

Color_Off="\033[0m"
Green="\033[0;32m"
Cyan="\033[0;36m"
BPurple="\033[1;35m"

OS=`uname -s | tr '[:upper:]' '[:lower:]'`;

echo "$BPurple"'Installing MDC dependencies\n\n'"$Color_Off"

echo "$Cyan"'Downloading C++ headers\n'"$Color_Off";

wget https://github.com/nlohmann/json/releases/download/v2.1.1/json.hpp --show-progress -qO /tmp/json.hpp;

wget https://raw.githubusercontent.com/jarro2783/cxxopts/9db62cb338aeaed1fec5806f6b5d9781f5e19e4c/include/cxxopts.hpp --show-progress -qO /tmp/cxxopts.hpp;

echo "$Cyan"'\nMaking headers available\n'"$Color_Off";

if [[ "$OS" == "darwin" ]]; then
  mv /tmp/json.hpp /usr/local/include/json.hpp;
  mv /tmp/cxxopts.hpp /usr/local/include/cxxopts.hpp;
else
  mv /tmp/json.hpp /usr/include/json.hpp;
  mv /tmp/cxxopts.hpp /usr/include/cxxopts.hpp;
fi

echo "$Cyan"'Installing libraries\n'"$Color_Off"

if [[ "$OS" == "darwin" ]]; then
  brew install llvm eigen boost;
  ln -s /usr/local/opt/llvm/bin/clang++ /usr/local/bin/clang-omp++
else
  apt-get install clang++ libomp-dev libboost-dev libeigen3-dev;
  ln -s /usr/bin/clang++ /usr/bin/clang-omp++;
fi

echo "$Green"'\nAll done! Run `make` to compile.'"$Color_Off"
