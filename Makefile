CXXFILES_CORE = $(wildcard src/include/**/*.cpp)
OFILES_CORE = $(CXXFILES_CORE:%.cpp=%.o)

CXXFILES_GMDL = $(wildcard src/gmdl/*.cpp)
OFILES_GMDL = $(CXXFILES_GMDL:%.cpp=%.o)
PROGRAMS_GMDL = gmdl.app

INCLUDES= -I. -Isrc/include
INCLUDES+= -I/usr/local/include
INCLUDES+= -I/usr/local/include/eigen3
INCLUDES+= -I/usr/include/eigen3
INCLUDES+= -I/usr/lib/llvm-3.8/include
INCLUDES+= -I/usr/lib/clang/3.8.0/include/
INCLUDES+= -I/usr/local/opt/llvm/include -fopenmp

CXXFLAGS = -std=c++11 -m64 -ggdb -O3 $(INCLUDES) -fPIC -fpic -Wall -Wextra -Wno-sign-compare -Wno-overloaded-virtual
# CXXFLAGS+= -lclang-3.8 -lstdc++
# CXXFLAGS+= -ffast-math
CXXFLAGS+= -DNDEBUG -DEIGEN_NO_DEBUG -msse2 -msse3 -mssse3 -msse4 -msse4.1 -msse4.2
CXXFLAGS+= -fopenmp
CXXFLAGS+= -finline -fbuiltin #-fexpensive-optimizations

LDFLAGS= -L/usr/local/opt/llvm/lib

CXX = time clang++

.PHONY: all
all: $(PROGRAMS_GMDL)

%.o:: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(PROGRAMS_GMDL): $(OFILES_CORE) $(OFILES_GMDL)
	$(CXX) -o $(PROGRAMS_GMDL) $^ -fopenmp $(LDFLAGS)

depend:
	$(CXX) -MM $(CXXFLAGS) $(CXXFILES_GMDL) >> .deps

clean:
	rm -f $(OFILES_CORE) $(OFILES_GMDL) $(PROGRAMS_GMDL)

-include .deps
