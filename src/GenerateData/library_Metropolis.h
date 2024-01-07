#include <set>

/********************************************************************/
/*********************    METROPOLIS.CPP    *************************/
/********************************************************************/
/*******************    SAMPLES A DATASET    ************************/
/********************************************************************/
map<unsigned int, list<Interaction>> write_interactions_metropolis(double J, string file);
double delta_energy(__int128_t state, list<Interaction> edges);
void sample_data_metropolis(double J, unsigned int s, string input_file, string output_filename, unsigned int N = 1000);

/********************************************************************/
/*****************    GENERATE_DATA_EXACT.CPP    ********************/
/********************************************************************/
/*******************    SAMPLES A DATASET    ************************/
/********************************************************************/
// *** list_I = list of operators and parameters of the model
// *** N = nb of datapoints

// *** Samples a dataset with "N" datapoints:  ******************** //
void Sample_dataset(list<Interaction> list_I, string output_filename, unsigned int N = 1000);

// **************************************************************** //
// *** Samples a dataset with "N" datapoints,  ******************** //
// *** and compute the model averages and data averages  ********** //
// *** of all the model operators:                       ********** //
// *** fill this information in "list_I"                 ********** //
// **************************************************************** //
void Sample_dataset_AND_Print_ModelData_Info(list<Interaction>& list_I, string output_filename, unsigned int N = 1000);

/********************************************************************/
/************    PRINT INFORMATION ABOUT THE MODEL    ***************/
/**********    and the GENERATED DATA if there were    **************/
/********************************************************************/
double* Probability_AllStates_Ising(list<Interaction> list_I, double* Z);
int Op_Ising(uint32_t Op, uint32_t state);
void Model_averages_Ising(list<Interaction>& list_I);

/******************************************************************************/
/************************** Best Basis ****************************************/
/******************************************************************************/

set<Operator> PairOp_m1_data_rank(map<uint32_t, unsigned int> Nset, unsigned int N);
set<Operator> AllOp_m1_data_rank(map<uint32_t, unsigned int> Nset, unsigned int N);
void add_new_basisOp(uint32_t new_basisOp, set<uint32_t> &OpSet, vector<unsigned int> &vect_op_used, int layer);
vector<pair<Operator, bool> > select_best_basis(set<Operator> &allOp_data, double Nd, double *LogL, int basis_size = n);
void printfile_BestBasis(vector<pair<Operator, bool>> Op, double Nd, string name);
list<__int128_t> Original_basis();
