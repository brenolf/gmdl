#!/bin/bash

Color_Off="\033[0m"
Green="\033[0;32m"
Cyan="\033[0;36m"
BPurple="\033[1;35m"

OS=`uname -s | tr '[:upper:]' '[:lower:]'`;

echo "$BPurple"'Installing GMDL dependencies\n\n'"$Color_Off"

echo "$Cyan"'Downloading C++ headers\n'"$Color_Off";

wget https://github.com/nlohmann/json/releases/download/v2.1.1/json.hpp --show-progress -qO /tmp/json.hpp;

wget https://raw.githubusercontent.com/tanakh/cmdline/master/cmdline.h --show-progress -qO /tmp/cmdline.hpp;

echo "$Cyan"'\nMaking headers available\n'"$Color_Off";

if test "$OS" = "darwin"; then
  mv /tmp/json.hpp /usr/local/include/json.hpp;
  mv /tmp/cmdline.hpp /usr/local/include/cmdline.hpp;
else
  mv /tmp/json.hpp /usr/include/json.hpp;
  mv /tmp/cmdline.hpp /usr/include/cmdline.hpp;
fi

echo "$Cyan"'Installing libraries\n'"$Color_Off"

if test "$OS" = "darwin"; then
  brew install llvm eigen boost;
  # ln -s /usr/local/opt/llvm/bin/clang++ /usr/local/bin/clang-omp++
else
  sudo apt-get install clang-3.8 libomp-dev libboost-dev libeigen3-dev;
  sudo ln -s /usr/bin/clang++-3.8 /usr/bin/clang++
fi

echo "$Green"'\nAll done! Run `make` to compile.'"$Color_Off"
