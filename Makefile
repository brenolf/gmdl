CXXFILES_CORE = $(wildcard src/include/**/*.cpp)
OFILES_CORE = $(CXXFILES_CORE:%.cpp=%.o)

CXXFILES_MDC = $(wildcard src/mdc/*.cpp)
OFILES_MDC = $(CXXFILES_MDC:%.cpp=%.o)
PROGRAMS_MDC = mdc.app

INCLUDES= -I. -Isrc/include
INCLUDES+= -I/usr/local/include
INCLUDES+= -I/usr/local/include/eigen3
INCLUDES+= -I/usr/include/eigen3
INCLUDES+= -I/usr/lib/llvm-3.8/include
INCLUDES+= -I/usr/lib/clang/3.8.0/include/

CXXFLAGS = -std=c++11 -m64 -ggdb -O3 $(INCLUDES) -fPIC -fpic -Wall -Wextra -Wno-sign-compare -Wno-overloaded-virtual
# CXXFLAGS+= -lclang-3.8 -lstdc++
# CXXFLAGS+= -ffast-math
CXXFLAGS+= -DNDEBUG -DEIGEN_NO_DEBUG -msse2 -msse3 -mssse3 -msse4 -msse4.1 -msse4.2
CXXFLAGS+= -fopenmp
CXXFLAGS+= -finline -fbuiltin #-fexpensive-optimizations

CXX = time clang++

.PHONY: all
all: $(PROGRAMS_MDC)

%.o:: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(PROGRAMS_MDC): $(OFILES_CORE) $(OFILES_MDC)
	$(CXX) -o $(PROGRAMS_MDC) $^ -fopenmp

depend:
	$(CXX) -MM $(CXXFLAGS) $(CXXFILES_MDC) >> .deps

clean:
	rm -f $(OFILES_CORE) $(OFILES_MDC) $(PROGRAMS_MDC)

-include .deps
