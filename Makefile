CXXFILES_CORE = $(wildcard src/include/**/*.cpp)
OFILES_CORE = $(CXXFILES_CORE:%.cpp=%.o)

CXXFILES_MDC = $(wildcard src/mdc/*.cpp)
OFILES_MDC = $(CXXFILES_MDC:%.cpp=%.o)
PROGRAMS_MDC = src/mdc/mdc

INCLUDES= -I. -Isrc/include
INCLUDES+= -I/usr/local/include
INCLUDES+= -I/usr/local/include/eigen3
INCLUDES+= -I/usr/include/eigen3

CXXFLAGS = -std=c++1y -m64 -ggdb -O3 -ffast-math $(INCLUDES) -fPIC -fpic -Wall -Wextra -Wno-sign-compare -Wno-overloaded-virtual

CXXFLAGS+= -DNDEBUG -DEIGEN_NO_DEBUG -msse2 -msse3 -mssse3 -msse4 -msse4.1 -msse4.2
CXXFLAGS+= -fopenmp
CXXFLAGS+= -finline -fbuiltin #-fexpensive-optimizations

CXX = clang-omp++

.PHONY: all
all: $(PROGRAMS_MDC)

%.o:: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(PROGRAMS_MDC): $(OFILES_CORE) $(OFILES_MDC)
	$(CXX) -o mdc.app $^ -fopenmp

depend:
	$(CXX) -MM $(CXXFLAGS) $(CXXFILES_MDC) >> .deps

clean:
	rm -f $(OFILES_CORE) $(OFILES_MDC) $(PROGRAMS_MDC)

-include .deps
