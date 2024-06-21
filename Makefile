
########################################################################################################################
######## ENTER THE FOLLOWING IN YOUR TERMINAL:
#### TO COMPILE C++:  					make
#### TO RUN C++: 						make run
#### TO COMPILE PYTHON MODULE: 			make py_module
#### TO CLEAN:  						make clean    --> to use only when you are completely done
########################################################################################################################

######## C++ compiler:
CC = g++
CXX = g++ 	# Flag for implicit rules: compilation of c++ files
#CXXFLAGS = ERROR NO DEFAULT RULE
CXXFLAGS = -std=c++11 -Wall -O3 #-Wextra

######## For PYTHON Package using PYBIND:
MODULE_NAME = MinCompSpin
OUTFILE = ./packages/$(MODULE_NAME)$(shell py -c "from sysconfig import get_config_var; print(get_config_var('EXT_SUFFIX'))") # Generate the filename python wants for its modules

# OUTFILE = $(MODULE_NAME)$(shell python3-config --extension-suffix) # Generate the module filename

# -fvisibility=hidden see:
# https://pybind11.readthedocs.io/en/stable/faq.html#someclass-declared-with-greater-visibility-than-the-type-of-its-field-someclass-member-wattributes
#CFLAGS_PYBIND_COMP = -IMinCompSpin_Greedy/Libraries -fvisibility=hidden
PYBIND_FLAGS = -shared
PYTHON_INCLUDE = -I/cygdrive/c/Python312/include
PYTHON_LIB = -L/cygdrive/c/Python312/libs -lpython312

UNAME := $(shell uname)
ifeq ($(UNAME), CYGWIN_NT-10.0) # Windows
	CXXFLAGS += -m64 -fPIC -U__STRICT_ANSI__ 
	PYBIND_FLAGS += -m64 -fPIC -U__STRICT_ANSI__
	PYBIND_INCLUDE = $(shell py -m pybind11 --includes)
else ifeq ($(UNAME), Darwin) # Mac
	PYBIND_FLAGS += -undefined dynamic_lookup
	PYBIND_INCLUDE = $(shell python -m pybind11 --includes)
else ifeq ($(UNAME), Linux) # Linux
	PYBIND_FLAGS += -fPIC -Wl,--export-dynamic
	PYBIND_INCLUDE = $(shell python -m pybind11 --includes)
else
$(error Couldn't find OS: $(UNAME))
endif

### MCM Files:
DIR_MCM = src/MinCompSpin
objects_MCM = tools.o LogE_LogL.o Complexity.o MCM_info.o Basis_Choice.o P_s.o Operations_OnData.o info_quant.o BestMCM_GreedySearch.o
OBJS := $(objects_MCM:%=$(DIR_MCM)/%)  # String substitution: add the directory name: # As an example, hello.o turns into ./MCM/hello.o

### Public Libraries:
DIR_LIB = includes

MCM_GREEDY_TARGET = MinCompSpin_Greedy/MCM_Greedy.out
MCM_GREEDY_BUILD_PATH = MinCompSpin_Greedy

########################################################################################################################
### make all:  ##### Compile C++ codes (with main):

MCM_Greedy.out: $(OBJS) $(DIR_LIB)/main.o $(DIR_LIB)/main_routines.o $(DIR_LIB)/library.hpp.gch
	g++ $(CXXFLAGS) $(DIR_LIB)/main.o $(DIR_LIB)/main_routines.o $(OBJS) -o MCM_Greedy.out

$(DIR_LIB)/main.o: main.cpp $(DIR_LIB)/default_datafile.h $(DIR_LIB)/library.hpp.gch
	g++ $(CXXFLAGS) -c main.cpp -include $(DIR_LIB)/library.hpp -o $(DIR_LIB)/main.o   # Compile main.cpp

$(DIR_LIB)/main_routines.o: $(DIR_LIB)/main_routines.cpp $(DIR_LIB)/library.hpp.gch
	g++ $(CXXFLAGS) -c $(DIR_LIB)/main_routines.cpp -include $(DIR_LIB)/library.hpp -o $(DIR_LIB)/main_routines.o

$(DIR_LIB)/library.hpp.gch: $(DIR_LIB)/library.hpp
	g++ $(CXXFLAGS) -c $(DIR_LIB)/library.hpp

#all: $(OUTFILE)

#$(OUTFILE): $(MCM_GREEDY_TARGET) library.o
#	$(CC) $(CFLAGS) $(CFLAGS_PYBIND_LINK) -o $(OUTFILE) $(MCM_GREEDY_TARGET) library.o

#library.o: library.cpp
#	$(CC) $(CFLAGS) $(CFLAGS_PYBIND_COMP) $(PYBIND_INCLUDE) -o library.o -c library.cpp

#$(MCM_GREEDY_TARGET):
#	cd $(MCM_GREEDY_BUILD_PATH) && $(MAKE)

########################################################################################################################
### make py_module:  ##### Compile Python module:

#py_module: $(OBJS)
#	g++ -std=c++11 -O3 -shared -undefined dynamic_lookup -include includes/library.hpp $(PYBIND_INCLUDE) binding/library.cpp $(OBJS) -o $(OUTFILE) includes/main_routines.cpp

py_module: 
	@echo "TEST"
	@echo $(PYBIND_INCLUDE)
	g++ -std=c++11 -O3 $(PYBIND_FLAGS) -include includes/library.hpp includes/main_routines.cpp $(PYBIND_INCLUDE) binding/library.cpp $(OBJS) -o $(OUTFILE) $(PYTHON_LIB)


########################################################################################################################
### make run:  ##### Execute C++ codes (with main):
run: MCM_Greedy.out
	./MCM_Greedy.out $(datafilename) $n

########################################################################################################################
clean:
	rm -f $(OUTFILE)
	rm -f $(OBJS) 
	rm -f $(DIR_LIB)/library.hpp.gch $(DIR_LIB)/main_routines.o $(DIR_LIB)/main.o 
	rm -f MCM_Greedy.out

#	rm $(OUTFILE) || true
#	rm library.o  || true
#	(cd $(MCM_GREEDY_BUILD_PATH) && $(MAKE) clean) || true

