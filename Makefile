CXX=g++
STD=-std=c++17
WARN=-Wall -Wextra -Wpedantic -Wconversion
OPT=-O3

# Directory of ampsci.... better way?
AMPSCI=/home/uqbrob14/ampsci/src

x-section: x-section.cpp StandardHaloModel.cpp StandardHaloModel.hpp
	$(CXX) $(STD) -o $@ $@.cpp StandardHaloModel.cpp $(WARN) $(OPT) -I$(AMPSCI)