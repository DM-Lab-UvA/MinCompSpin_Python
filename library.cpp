#include <pybind11/pybind11.h>
#include <stdio.h>

#include "MinCompSpin_Greedy/Libraries/library.hpp"

int main(int argc, const char *argv[]);

const int argc = 3;
int main_wrapper(char *path, char *n_string_buffer) {
    const char *argv[argc] = {
        "invalid",
        path,
        n_string_buffer,
    };
    main(argc, argv);
    return 0;
}

int add(int i, int j) {
    printf("func: %p res %p\n", Original_Basis);
    return i + j;
}

PYBIND11_MODULE(MinCompSpin, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("add", &add, "A function that adds two numbers");
    m.def("main", &main_wrapper, "the main function");
    
    m.def("Read_BasisOp_BinaryRepresentation", &Read_BasisOp_BinaryRepresentation);
    m.def("Read_BasisOp_IntegerRepresentation", &Read_BasisOp_IntegerRepresentation);
    m.def("PrintTerm_Basis", &PrintTerm_Basis);
    m.def("read_datafile", &read_datafile);
    m.def("build_Kset", &build_Kset);
    // m.def("LogL_SubCM", &LogL_SubCM);
    // m.def("LogE_SubCM", &LogE_SubCM);
    // m.def("GeomComplexity_SubCM", &GeomComplexity_SubCM);
    // m.def("ParamComplexity_SubCM", &ParamComplexity_SubCM);
    m.def("LogL_CM", &LogL_CM);
    m.def("LogL_MCM", &LogL_MCM);
    m.def("LogE_MCM", &LogE_MCM);
    m.def("Complexity_MCM", &Complexity_MCM);
    // m.def("Create_MCM", &Create_MCM);
    // m.def("Read_MCMParts_BinaryRepresentation", &Read_MCMParts_BinaryRepresentation);
    m.def("check_partition", &check_partition);
    // m.def("PrintTerminal_MCM_Info", &PrintTerminal_MCM_Info);
    // m.def("PrintInfo_All_Indep_Models", &PrintInfo_All_Indep_Models);
    // m.def("PrintInfo_All_SubComplete_Models", &PrintInfo_All_SubComplete_Models);
    // m.def("MCM_GivenRank_r", &MCM_GivenRank_r);
    // m.def("MCM_AllRank_SmallerThan_r_Ordered", &MCM_AllRank_SmallerThan_r_Ordered);
    // m.def("MCM_AllRank_SmallerThan_r_nonOrdered", &MCM_AllRank_SmallerThan_r_nonOrdered);
    m.def("Original_Basis", &Original_Basis);
    // m.def("count_set_bits", &count_set_bits);
    // m.def("int_to_bstring", &int_to_bstring);
}

