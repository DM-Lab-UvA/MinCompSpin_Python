PYBIND_INCLUDE = $(shell python3 -m pybind11 --includes)
MCM_INCLUDE = MinCompSpin_Greedy/Libraries/library.hpp
MCM_LINK = MinCompSpin_Greedy/Libraries/MCM/*.o
CFLAGS = -O3 -std=c++11 -fPIC -shared -IMinCompSpin_Greedy/Libraries
OUTFILE = MinCompSpin$(shell python-config --extension-suffix)

MCM_GREEDY_PATH = MinCompSpin_Greedy/Libraries/MCM
# MCM_GREEDY_CPPS = $(shell find $(MCM_GREEDY_PATH) -name "*.cpp" -not \( -name "P_s.cpp" \))
MCM_GREEDY_CPPS = $(shell find $(MCM_GREEDY_PATH) -name "*.cpp") MinCompSpin_Greedy/main.cpp MinCompSpin_Greedy/Libraries/main_routines.cpp
# MCM_GREEDY_CPPS = $(shell find $(MCM_GREEDY_PATH) -name "*.cpp")
MCM_GREEDY_OBJS = $(patsubst %.cpp,%.o,$(MCM_GREEDY_CPPS))

CC = g++

all: $(OUTFILE)

$(OUTFILE): $(MCM_GREEDY_OBJS) library.o
	$(CC) $(CFLAGS) -o $(OUTFILE) $(MCM_GREEDY_OBJS) library.o

library.o: library.cpp
	$(CC) $(PYBIND_INCLUDE) $(CFLAGS) -o library.o -c library.cpp

$(MCM_GREEDY_PATH)/%.o: $(MCM_GREEDY_PATH)/%.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm $(OUTFILE) || true
	rm library.o  || true
	rm $(MCM_GREEDY_OBJS) || true
