CXX=g++
STD=-std=c++17
WARN=-Wall -Wextra -Wpedantic -Wconversion
OPT=-O3

EXES = dmex dmex2 eimpact

SRC=./src
# Simple: depends on everything
SOURCES = $(SRC)/*.cpp
DEPS = $(SOURCES) $(SRC)/*.hpp

# just compile + link at same time
COMP=$(CXX) $(STD) -o $@ $@.cpp $(SOURCES) $(WARN) $(OPT) -I$(SRC) -fopenmp

all: $(EXES)

dmex: dmex.cpp $(DEPS)
	$(COMP)

dmex2: dmex2.cpp $(DEPS)
	$(COMP)

eimpact: eimpact.cpp $(DEPS)
	$(COMP)

clean:
	rm -rfv $(EXES)