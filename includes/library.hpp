#include <list>
#include <map>
#include <vector>
#include <string>
//#include <utility>


/******************************************************************************/
/******************************************************************************/
/*********************     SPECIFY the BASIS    *******************************/
/******************************************************************************/
/******************************************************************************/

/*** Original Basis:    *******************************************************/
/******************************************************************************/
std::list<__int128_t> Original_Basis(unsigned int r);   // return the original basis, i.e., {s1, s2, ..., sn}

/*** READ BASIS from a FILE:    ***********************************************/
/******************************************************************************/
std::list<__int128_t> Read_BasisOp_BinaryRepresentation(unsigned int r, std::string Basis_binary_filename); // = basis_BinaryRepresentation_filename);   // filename to specify in data.h
std::list<__int128_t> Read_BasisOp_IntegerRepresentation(std::string Basis_integer_filename); // = basis_IntegerRepresentation_filename); 

/*** Print Basis Info in the Terminal:    *************************************/
/******************************************************************************/
void PrintTerm_Basis(std::list<__int128_t> Basis_li, unsigned int r);


/******************************************************************************/
/******************************************************************************/
/******************     READ and TRANSFORM DATA    ****************************/
/******************************************************************************/
/******************************************************************************/

/*** READ DATA and STORE data in Nset:    *************************************/
/******************************************************************************/
std::vector<std::pair<__int128_t, unsigned int>> read_datafile(unsigned int *N, std::string file, unsigned int r); // filename to specify in data.h


/*** DATA CHANGE of BASIS:    *************************************************/
/******************************************************************************/
// *** Build Kset with the following definitions:
// *** mu_m = states of the systems written in the basis specified in `list<uint32_t> Basis`
// *** Kset[sig_m] = Number of times the state mu_m appears in the transformed dataset
//
// *** Rem: the new basis can have a lower dimension then the original dataset; 
// *** in which case the function will reduce the dataset to the subspace defined by the specified basis.
std::vector<std::pair<__int128_t, unsigned int>> build_Kset(std::vector<std::pair<__int128_t, unsigned int>> Nset, std::list<__int128_t> Basis);


/*** REDUCE Kset:    **********************************************************/
/******************************************************************************/
// *** Remove all the states that occur less than K times in the dataset:
std::vector<std::pair<__int128_t, unsigned int>> Reduce_Kset(std::vector<std::pair<__int128_t, unsigned int>> Kset, unsigned int K, unsigned int *N_new);


/******************************************************************************/
/******************************************************************************/
/***************** Log-LIKELIHOOD (LogL), Log-EVIDENCE (LogE) *****************/
/***************************  and COMPLEXITY   ********************************/
/******************************************************************************/
/******************************************************************************/

/****************   for a sub-Complete Model (SubCM)   ************************/
/**********  restricted to the subspace of the Sub-Complete Model  ************/
/******************************************************************************/
// *** the SubCM is the one specified in Ai;
// *** Ai must be an integer encoded on at least n bits, where each 1 indicates the basis elements included in the part:
// *** For ex. Ai = 01001 is encoded on n=5 basis elements, and element Op1 and Op4 belong to the part;
// *** Rem: Basis elements are ordered from the right to the left.
double LogL_ICC(std::vector<std::pair<__int128_t, unsigned int>> Kset, __int128_t Ai, unsigned int N);
double LogE_ICC(std::vector<std::pair<__int128_t, unsigned int>> Kset, __int128_t Ai, unsigned int N);

// *** Complexity of a SC model based on m basis Operators: m >= 1. Rem: C_geom(m=1) = log(pi):
double GeomComplexity_ICC(unsigned int m);                  // Geometric complexity
double ParamComplexity_ICC(unsigned int m, unsigned int N); // Complexity due to the number of parameters

/******************   for a Complete Model (CM)   *****************************/
/******************************************************************************/
double LogL_CM(std::vector<std::pair<__int128_t, unsigned int>> Kset, unsigned int N);

/****************************    for a MCM     ********************************/
/******************************************************************************/

double Complexity_MCM(std::map<unsigned int, __int128_t> Partition, unsigned int N, double *C_param, double *C_geom);

double LogL_MCM(std::vector<std::pair<__int128_t, unsigned int>> Kset, std::map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r);
double LogE_MCM(std::vector<std::pair<__int128_t, unsigned int>> Kset, std::map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r);

double LogL_MCM_infoICC(std::vector<std::pair<__int128_t, unsigned int>> Kset, std::map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r);
double LogE_MCM_infoICC(std::vector<std::pair<__int128_t, unsigned int>> Kset, std::map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r);


/******************************************************************************/
/******************************************************************************/
/********************   DEFINE MCMs and PRINT INFO   **************************/
/***********************    see  "MCM_info.ccp"   *****************************/
/******************************************************************************/
/******************************************************************************/

// *** Read MCM from a file:
//map<unsigned int, __int128_t> read_MCM_fromfile(string Input_MCM_file = communityfile);

// *** Check that the provided model corresponds to a partition of the basis variables (i.e. properly defines an MCM):
std::pair<bool, unsigned int> check_partition(std::map<unsigned int, __int128_t> Partition);  // the second element is the rank of the partition (dimension of the MCM)

// *** Print information about the MCM specified in `MCM_Partition`:
void Print_MCM_Partition(std::map<unsigned int, __int128_t> partition, unsigned int r);

void PrintFile_MCM_Info(std::list<__int128_t> Basis, std::map<unsigned int, __int128_t> MCM_Partition, unsigned int r, std::string filename = "Result");


/******************************************************************************/
/******************************************************************************/
/********************   Statistical Properties   ******************************/
/******************************************************************************/
/******************************************************************************/

// *** Entropy of a dataset:
double Entropy(std::vector<std::pair<__int128_t, unsigned int>> Kset, unsigned int N);

// *** Distance between two probability distributions:
// double KL_divergence(map<__int128_t, unsigned int> Kset, map<unsigned int, __int128_t> Partition, unsigned int N);
// double JS_divergence(std::map<__int128_t, double> Prob1, std::map<__int128_t, double> Prob2, unsigned int N);

// *** Compare two partitions:
double Var_of_Inf(std::map<unsigned int, __int128_t> Partition1, std::map<unsigned int, __int128_t> Partition2, unsigned int r);
double Norm_Mut_info(std::map<unsigned int, __int128_t> Partition1, std::map<unsigned int, __int128_t> Partition2, unsigned int r);

// *** Check if the MCM in "fp1" is a sub-partition of the MCM in "fp2":
bool is_subset(std::map<unsigned int, __int128_t> fp1, std::map<unsigned int, __int128_t> fp2);


/******************************************************************************/
/******************************************************************************/
/***************************   PRINT TO FILE:  ********************************/
/******************   DATA VS MODEL STATE PROBABILITIES  **********************/
/******************************************************************************/
/******************************************************************************/

void PrintFile_StateProbabilites_NewBasis(std::vector<std::pair<__int128_t, unsigned int>> Kset, std::map<unsigned int, __int128_t> MCM_Partition, unsigned int N, unsigned int r, std::string filename = "Result");
void PrintFile_StateProbabilites_OriginalBasis(std::vector<std::pair<__int128_t, unsigned int>> Nset, std::list<__int128_t> Basis, std::map<unsigned int, __int128_t> MCM_Partition, unsigned int N, unsigned int r, std::string filename = "Result");

// *** Create distributions of empirical data and MCMs:
// std::map<__int128_t, double> emp_dist(std::map<__int128_t, unsigned int> Kset, unsigned int N, unsigned int r);
// std::map<__int128_t, double> MCM_distr(std::map<__int128_t, unsigned int> Kset, std::map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r);

/******************************************************************************/
/******************************************************************************/
/**********************   Metropolis algorithm   ******************************/
/******************************************************************************/
/******************************************************************************/

//list<Interaction> write_interactions(double J, string file);

