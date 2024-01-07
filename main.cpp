// To compile: make
// To run: or: make run
//         or: ./GreedySearch.out data_filename n
//
// #define _USE_MATH_DEFINES 
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <map>
#include <vector>
#include <cmath>       /* tgamma */
#include <random>

#include <ctime> // for chrono
#include <ratio> // for chrono
#include <chrono> // for chrono

using namespace std;
using namespace std::chrono;

/*****************************************************************************************/
/*************************   CONSTANT VARIABLES  *****************************************/
/*****************************************************************************************/
#include "includes/default_datafile.h"

string OutputFile_Add_Location(string filename)
{
    return (OUTPUT_directory + filename);
}

/*****************************************************************************************/
/****************   GREEDY SEARCH:   Useful functions and routines    ********************/
/*****************************************************************************************/
// **** Find the best MCM, Greedy Search:
map<unsigned int, __int128_t> MCM_GreedySearch(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, bool print_it = false);
map<unsigned int, __int128_t> MCM_GreedySearch_AND_printInfo(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, bool print_it = false);

// **** Find the best MCM, Greedy Search starting from the model MCM_0:
map<unsigned int, __int128_t> MCM_GreedySearch_MCM0(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, map<unsigned int, __int128_t> MCM_0, bool print_it = false);

// *** Greedy Search on Reduced dataset:
map<unsigned int, __int128_t> MCM_ReducedGreedySearch_AND_PrintInfo(vector<pair<__int128_t, unsigned int>> Kset, unsigned int K, unsigned int N, unsigned int r, bool print_it = false);

// *** Compare two MCMs:
void compare_two_MCMs_AND_printInfo(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, map<unsigned int, __int128_t> fp1, map<unsigned int, __int128_t> fp2);


/*****************************************************************************************/
/*****************************   IMPORT an MCM from a FILE   *****************************/
/*****************************************************************************************/
// *** Read MCM from a file:
map<unsigned int, __int128_t> read_MCM_fromfile(string Input_MCM_file, unsigned int r);
map<unsigned int, __int128_t> read_MCM_fromfile_AND_printInfo(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, string Input_MCM_file, unsigned int r);


/******************************************************************************/
/***************************   ADD OUTPUT FOLDER    ***************************/
/******************************************************************************/

//// ** location of the output folder:
//string OutputFile_Add_Location(string filename);    // defined in ./src/P_s.cpp
//{
//    return (OUTPUT_directory + filename);
//}

/****************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************/
/**************************************************************************   """ TUTORIAL  """    **************************************************************************/
/****************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************/

void tutorial(vector<pair<__int128_t, unsigned int>> Nset, unsigned int N,  unsigned int n)
{
    cout << endl << "*******************************************************************************************";  
    cout << endl << "******************************  CHOICE OF THE BASIS:  *************************************";
    cout << endl << "*******************************************************************************************" << endl;

// original basis of the data: this is the most natural choice a priori:
//    list<__int128_t> Basis_li = Original_Basis(n);

//// *** The basis can also be read from a file: Ex. the following files contain the best basis for the SCOTUS dataset:
   list<__int128_t> Basis_li = Read_BasisOp_IntegerRepresentation(basis_IntegerRepresentation_filename);
//   list<__int128_t> Basis_li = Read_BasisOp_BinaryRepresentation(n, basis_IntegerRepresentation_filename);

    PrintTerm_Basis(Basis_li, n);

    cout << endl << "*******************************************************************************************";
    cout << endl << "**********************  TRANSFORM the DATA in the CHOSEN BASIS   **************************";
    cout << endl << "**********************************   Build Kset:   ****************************************";
    cout << endl << "*******************************************************************************************" << endl;

//// *** Transform the data in the specified in Basis_SCModel[];

    vector<pair<__int128_t, unsigned int>> Kset = build_Kset(Nset, Basis_li);
    cout << "\t Kset.size() = " << Kset.size() << endl;


    cout << endl << "*******************************************************************************************";
    cout << endl << "***********************  HIERARCHICAL GREEDY MERGING: BY STEPS:  **************************";
    cout << endl << "*******************************************************************************************" << endl;

    bool print_checkpoint = true;  

//// *** Finds the best MCM:
    map<unsigned int, __int128_t> fp1 = MCM_GreedySearch(Kset, N, n, print_checkpoint);

//// *** Print Log-Evidence:  
    double LogE_fp1 = LogE_MCM(Kset, fp1, N, n);
    cout << "Log-Evidence(MCM) = " << LogE_fp1 << "\t = " << LogE_fp1/((double) N)/log(2.) << " bits per datapoint \t" << endl;

//// *** Print max-Log-Likelihood:  
    double LogL_fp1 = LogL_MCM(Kset, fp1, N, n);
    cout << "Max Log-Likelihood(MCM) = " << LogL_fp1 << "\t = " << LogL_fp1/((double) N)/log(2.) << " bits per datapoint \t" << endl;

    cout << endl << "*******************************************************************************************";
    cout << endl << "******************************  READ an MCM from a FILE   *********************************";
    cout << endl << "*******************************************************************************************" << endl;

    cout << "#########  EX. READ a CHOSEN MCM:  #########" << endl;

    // the file communityfile = "INPUT/SCOTUS_Communities_inBestBasis.dat" contains the best MCM in the best basis:
    map<unsigned int, __int128_t> fp2 = read_MCM_fromfile(communityfile, n);
    Print_MCM_Partition(fp2, n);

    cout << endl << "*******************************************************************************************";
    cout << endl << "*******************************  COMPARING TWO MCMs   *************************************";
    cout << endl << "*******************************************************************************************" << endl;

    compare_two_MCMs_AND_printInfo(Kset, N, n, fp1, fp2);

    cout << endl << "*******************************************************************************************";
    cout << endl << "***************************  Decomposition of Log-E   *************************************";
    cout << endl << "*******************************   over each ICC   *****************************************";
    cout << endl << "*******************************************************************************************" << endl;

    double LogE_final = LogE_MCM_infoICC(Kset, fp1, N, n);
    //cout << "Log-Evidence(MCM) = " << LogE_final << "\t = " << LogE_final/((double) N)/log(2.) << " bits per datapoint \t" << endl;

    cout << endl << "*******************************************************************************************";
    cout << endl << "*************************  Decomposition of Max-Log-L   ***********************************";
    cout << endl << "*******************************   over each ICC   *****************************************";
    cout << endl << "*******************************************************************************************" << endl;

    double LogL_final = LogL_MCM_infoICC(Kset, fp1, N, n);
    //cout << "Max-Log-Likelihood(MCM) = " << LogL_final << "\t = " << LogL_final/((double) N)/log(2.) << " bits per datapoint \t" << endl;

    cout << endl << "*******************************************************************************************";
    cout << endl << "***************************  Working with a Reduced Dataset   *****************************";
    cout << endl << "**********   Remove from Kset all the states that occur less than K times:   **************";
    cout << endl << "*******************************************************************************************" << endl;

    // All the states that occur less than K times will be removed from the dataset:
    unsigned int K=2;
    map<unsigned int, __int128_t> fp_reduced = MCM_ReducedGreedySearch_AND_PrintInfo(Kset, K, N, n);

    cout << endl << "*******************************************************************************************";
    cout << endl << "**********************  Print information about the found MCM:  ***************************";
    cout << endl << "*******************************************************************************************" << endl;

    // Prints 1) information about the MCM; 2) the state probabilities P(s) of observed states (in the Data VS MCM); 3) the probability P(k) of observing a state with k values "+1" (in the Data VS MCM) 
    PrintFile_StateProbabilites_OriginalBasis(Nset, Basis_li, fp1, N, n, "Result");

    // Print the state probabilities P(s) of observed states (in the Data VS MCM) using the data transformed in the bew basis:
    PrintFile_StateProbabilites_NewBasis(Kset, fp1, N, n, "Result");
}


/****************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************/
/**************************************************************************     MAIN FUNCTION      **************************************************************************/
/****************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************/

int main(int argc, char *argv[])
{
//// *** Read the arguments:
    string n_string_buffer = "";

    if (argc == 3)
    {
        datafilename = argv[1];
        n_string_buffer = argv[2];
        n = stoul(n_string_buffer);
    }
    else if (argc != 1)
    {
        cout << "The number of arguments must be either 0 or 2" << endl;
        return 0;
    }

    cout << "--->> Create the \"OUTPUT\" Folder: (if needed) ";
    system(("mkdir -p " + OUTPUT_directory).c_str());
    cout << endl;

    cout << endl << "*******************************************************************************************";
    cout << endl << "***********************************  READ THE DATA:  **************************************";
    cout << endl << "*******************************************************************************************" << endl;

    //cout << "Read the dataset: " << datafilename << endl;
    //cout << "Number of variables to read: n = " << n << endl;

    unsigned int N = 0; // will contain the number of datapoints in the dataset
    vector<pair<__int128_t, unsigned int>> Nset = read_datafile(&N, datafilename, n);

    if (N == 0)     // Terminate program if the file can't be found or read, or if it is empty:
        { 
        cout << " --->> Datafile cannot be read, or is empty; Terminate the program." << endl << endl;
        return 0; 
        } 

    cout << endl << "###### File has been read successfully:" << endl;
    cout << "\t Number of datapoints: N = " << N << endl;
    cout << "\t Number of different observed states = " << Nset.size() << endl;

    cout << endl << "*******************************************************************************************";
    cout << endl << "*****************************  HIERARCHICAL GREEDY MERGING:  ******************************";
    cout << endl << "********************************  in the ORIGINAL BASIS  **********************************";
    cout << endl << "*******************************************************************************************" << endl;

    bool print_checkpoint = false;  

//// *** Finds the best MCM and print information about it in the terminal:
    map<unsigned int, __int128_t> fp1 = MCM_GreedySearch_AND_printInfo(Nset, N, n, print_checkpoint);

    cout << endl << "*******************************************************************************************";
    cout << endl << "***************************  Decomposition of Log-E   *************************************";
    cout << endl << "*******************************   over each ICC   *****************************************";
    cout << endl << "*******************************************************************************************" << endl;
    LogE_MCM_infoICC(Nset, fp1, N, n);

/*
    cout << endl << "*******************************************************************************************************************";
    cout << endl << "*******************************************************************************************************************";
    cout << endl << "************************************************  TUTORIAL:  ******************************************************";
    cout << endl << "*******************************************************************************************************************";
    cout << endl << "*******************************************************************************************************************" << endl;

    tutorial(Nset, N, n);
*/

    return 0;
}
