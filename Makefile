UNAME := $(shell uname)
PYBIND_INCLUDE = $(shell python3 -m pybind11 --includes)
OUTFILE = MinCompSpin$(shell python-config --extension-suffix)

CFLAGS = -Wall -Wextra -O3 -std=c++11 
# -fvisibility=hidden see:
# https://pybind11.readthedocs.io/en/stable/faq.html#someclass-declared-with-greater-visibility-than-the-type-of-its-field-someclass-member-wattributes
CFLAGS_PYBIND_COMP = -IMinCompSpin_Greedy/Libraries -fvisibility=hidden
CFLAGS_PYBIND_LINK = -shared
ifeq ($(OS), Windows_NT) # Windows
	# I don't know what to do here but from what I saw online
	# it looks like it's the same as on linux.
	CFLAGS_PYBIND_COMP += -fPIC
else ifeq ($(UNAME), Darwin) # Mac
	CFLAGS_PYBIND_LINK += -undefined dynamic_lookup
else ifeq ($(UNAME), Linux) # Linux
	CFLAGS_PYBIND_COMP += -fPIC
else
$(error Couldn't find OS: $(UNAME))
endif

MCM_GREEDY_TARGET = MinCompSpin_Greedy/MCM_Greedy.out
MCM_GREEDY_BUILD_PATH = MinCompSpin_Greedy

CC = g++
CXXFLAGS = ERROR NO DEFAULT RULE

all: $(OUTFILE)

$(OUTFILE): $(MCM_GREEDY_TARGET) library.o
	$(CC) $(CFLAGS) $(CFLAGS_PYBIND_LINK) -o $(OUTFILE) $(MCM_GREEDY_TARGET) library.o

library.o: library.cpp
	$(CC) $(CFLAGS) $(CFLAGS_PYBIND_COMP) $(PYBIND_INCLUDE) -o library.o -c library.cpp

$(MCM_GREEDY_TARGET):
	cd $(MCM_GREEDY_BUILD_PATH) && $(MAKE)

clean:
	rm $(OUTFILE) || true
	rm library.o  || true
	(cd $(MCM_GREEDY_BUILD_PATH) && $(MAKE) clean) || true
