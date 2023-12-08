PYBIND_INCLUDE = $(shell python3 -m pybind11 --includes)
CFLAGS = -O3 -std=c++11 -fPIC -shared -IMinCompSpin_Greedy/Libraries 
OUTFILE = MinCompSpin$(shell python-config --extension-suffix)

MCM_GREEDY_TARGET = MinCompSpin_Greedy/MCM_Greedy.out
MCM_GREEDY_BUILD_PATH = MinCompSpin_Greedy

CC = g++
CXXFLAGS = ERROR NO DEFAULT RULE

all: $(OUTFILE)

$(OUTFILE): $(MCM_GREEDY_TARGET) library.o
	$(CC) $(CFLAGS) -o $(OUTFILE) $(MCM_GREEDY_TARGET) library.o

library.o: library.cpp
	$(CC) $(PYBIND_INCLUDE) $(CFLAGS) -Wall -Wextra -o library.o -c library.cpp

$(MCM_GREEDY_TARGET):
	cd $(MCM_GREEDY_BUILD_PATH) && $(MAKE)

clean:
	rm $(OUTFILE) || true
	rm library.o  || true
	(cd $(MCM_GREEDY_BUILD_PATH) && $(MAKE) clean) || true
