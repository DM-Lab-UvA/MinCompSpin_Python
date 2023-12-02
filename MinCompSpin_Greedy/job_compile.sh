
## COMPILE:

cd Libraries/MCM/
g++ -std=c++11 -O2 -c *.cpp

cd ../../
g++ -std=c++11 -O2 -c Libraries/library.hpp
g++ -std=c++11 -O2 Libraries/main_routines.cpp main.cpp -include Libraries/library.hpp Libraries/MCM/*.o -o MCM_Greedy.out
