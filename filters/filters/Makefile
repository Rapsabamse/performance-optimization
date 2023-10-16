# Author: David Holmqvist <daae19@student.bth.se>

CXX=g++
CXXFLAGS=-std=c++17 -g -Wunused -Wall -Wunused

all: blur threshold threshold_par blur_par

threshold: matrix ppm filters threshold.cpp
	$(CXX) $(CXXFLAGS) threshold.cpp matrix.o ppm.o filters.o -o threshold

blur: matrix ppm filters blur.cpp
	$(CXX) $(CXXFLAGS) blur.cpp matrix.o ppm.o filters.o -o blur

threshold_par: matrix ppm filters threshold_par.cpp
	$(CXX) $(CXXFLAGS) threshold_par.cpp matrix.o ppm.o filters.o -o threshold_par

blur_par: matrix ppm filters blur_par.cpp
	$(CXX) $(CXXFLAGS) blur_par.cpp matrix.o ppm.o filters.o -o blur_par

filters: matrix filters.hpp filters.cpp
	$(CXX) $(CXXFLAGS) -c filters.cpp -o filters.o

matrix: matrix.hpp matrix.cpp
	$(CXX) $(CXXFLAGS) -c matrix.cpp -o matrix.o

ppm: ppm.hpp ppm.cpp
	$(CXX) $(CXXFLAGS) -c ppm.cpp -o ppm.o

clean:
	rm -rf blur threshold *.ppm *.o *.dSYM 2> /dev/null
