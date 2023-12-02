#include <iostream>
#include <map>
#include <list>
#include <string>
#include <cmath>       /* lgamma */
#include <ctime> // for chrono
#include <ratio> // for chrono
#include <chrono> // for chrono

using namespace std;
using namespace std::chrono;

/******************************************************************************/
/**********************    CONSTANTS AND FUNCTIONS    *************************/
/******************************************************************************/
//#include "library.hpp"

/******************************************************************************/
/**********************   ROUTINES for GREEDY SEARCH   ************************/
/******************************************************************************/
map<unsigned int, __int128_t> MCM_GreedySearch(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, bool print_it = false);

map<unsigned int, __int128_t> MCM_GreedySearch_AND_printInfo(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, bool print_it = false)
{
    cout << "######### START GREEDY SEARCH #########" << endl;
    // *** Calculate the optimal partition
    auto start = chrono::system_clock::now();
    map<unsigned int, __int128_t> fp1 = MCM_GreedySearch(Kset, N, r, print_it);
    auto end = chrono::system_clock::now();

    // *** Time it takes to find partition
    chrono::duration<double> elapsed = end - start;
    cout << "Run time : " << elapsed.count() << "s" << endl << endl;

    cout << "######### ENTROPY OF THE DATA #########" << endl;
    // Entropy of dataset
    double S = Entropy(Kset, N);
    cout << "S = " << S << " bits. Range: [0, " << r << "]" << endl << endl;

    cout << "#########  GREEDY: BEST MCM FOUND:   #########" << endl;
    // Log evidence of MCM
    double LogE_g = LogE_MCM(Kset, fp1, N, r);
    double LogL_g = LogL_MCM(Kset, fp1, N, r);

    Print_MCM_Partition(fp1, r);

    cout << "Log-Evidence(MCM)       : " << LogE_g << "\t = " << LogE_g/( (double) N)/log(2.) << " bits per datapoint \t" << endl;
    cout << "Max-Log-Likelihood(MCM) : " << LogL_g << "\t = " << LogL_g/( (double) N)/log(2.) << " bits per datapoint \t" << endl << endl;
    cout << "Average community size  : " << (double)r / (double)fp1.size() << " variables per community \t" << endl << endl;

    return fp1;
}

/******************************************************************************/
/*****************   GREEDY SEARCH: CALL from ORIGINAL BASIS   ****************/
/******************************************************************************/

map<unsigned int, __int128_t> MCM_GreedySearch_OriginalBasis(vector<pair<__int128_t, unsigned int>> Nset, unsigned int N, unsigned int r)
{
    return MCM_GreedySearch(Nset, N, r);
}

map<unsigned int, __int128_t> MCM_GreedySearch_ChosenBasis(vector<pair<__int128_t, unsigned int>> Nset, unsigned int N, list<__int128_t> Basis_li)
{
    vector<pair<__int128_t, unsigned int>> Kset = build_Kset(Nset, Basis_li);
    return MCM_GreedySearch(Kset, N, Basis_li.size());
}

/******************************************************************************/
/********************   GREEDY SEARCH for REDUCED DATASET  ********************/
/******************************************************************************/

map<unsigned int, __int128_t> MCM_ReducedGreedySearch_AND_PrintInfo(vector<pair<__int128_t, unsigned int>> Kset_Vect, unsigned int K, unsigned int N, unsigned int r, bool print_it = false)
{
    auto start = chrono::system_clock::now();
    cout << "######### REDUCE Kset TO STATES OCCURING AT LEAST K TIMES: K = " << K << "  #########" << endl;
    unsigned int N_reduced = N;
    vector<pair<__int128_t, unsigned int>> Kset_reduced = Reduce_Kset(Kset_Vect, K, &N_reduced); // update N_reduced and return Kset_reduced

    cout << "After reduction: " << endl;
    cout << "\t Number of datapoints = " << N_reduced << endl;
    cout << "\t Number of different observed states = " << Kset_reduced.size() << endl;

    cout << endl << "######### ENTROPY OF THE REDUCED DATA #########" << endl;
    // Entropy of dataset
    double H = Entropy(Kset_reduced, N_reduced);
    cout << "Entropy, S = " << H << ". Range: [0, " << r << "]" << endl;

    cout << endl << "######### START GREEDY SEARCH #########" << endl;
    // *** Calculate the optimal partition
    map<unsigned int, __int128_t> fp_reduced = MCM_GreedySearch(Kset_reduced, N_reduced, r, print_it);
    auto end = chrono::system_clock::now();

    // *** Time it takes to find partition
    chrono::duration<double> elapsed = end - start;
    cout << "Elapsed time      : " << elapsed.count() << "s" << endl;

    cout << endl << "#########  REDUCED GREEDY: BEST MCM FOUND:  #########" << endl;
    // Log evidence of MCM
    double LE_reduced = LogE_MCM(Kset_reduced, fp_reduced, N_reduced, r);
    double LE = LogE_MCM(Kset_Vect, fp_reduced, N, r);
    Print_MCM_Partition(fp_reduced, r);

    cout << "Log-Evidence(MCM) for reduced data    : " << LE_reduced << "\t = " << LE_reduced/((double) N)/log(2.) << " bits per datapoint \t" << endl;
    cout << "Log-Evidence(MCM) for original data   : " << LE << "\t = " << LE/((double) N)/log(2.) << " bits per datapoint \t" << endl;

    cout << endl << "Average comm size : " << (double)r / (double)fp_reduced.size() << endl << endl;

    Kset_reduced.clear();

    return fp_reduced;
}

/******************************************************************************/
/*************************   READING MCM from a FILE   ************************/
/******************************************************************************/
map<unsigned int, __int128_t> read_MCM_fromfile(string Input_MCM_file, unsigned int r);


map<unsigned int, __int128_t> read_MCM_fromfile_AND_printInfo(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, string Input_MCM_file, unsigned int r)
{
    cout << "#########  THEORETICAL   #########" << endl;
    map<unsigned int, __int128_t> fp2 = read_MCM_fromfile(Input_MCM_file, r);

    double LE_t = LogE_MCM(Kset, fp2, N, r);
    Print_MCM_Partition(fp2, r);

    cout << "Log-evidence (MCM)     : " << LE_t << endl;
    cout << "Average community size : " << (double)r / (double)fp2.size() << endl << endl;

    return fp2;
}

/******************************************************************************/
/****************************   COMPARING TWO MCMs   **************************/
/******************************************************************************/
void compare_two_MCMs_AND_printInfo(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, map<unsigned int, __int128_t> fp1, map<unsigned int, __int128_t> fp2)
{
    cout <<  "#########  MCM 1:   #########" << endl;
    Print_MCM_Partition(fp1, r);
    double LogE_fp1 = LogE_MCM(Kset, fp1, N, r);
    double LogL_fp1 = LogL_MCM(Kset, fp1, N, r);

    cout << "      Log-Evidence(MCM) = " << LogE_fp1 << "\t = " << LogE_fp1/((double) N)/log(2.) << " bits per datapoint \t" << endl;
    cout << "Max Log-Likelihood(MCM) = " << LogL_fp1 << "\t = " << LogL_fp1/((double) N)/log(2.) << " bits per datapoint \t" << endl;
    cout << "Average community size  = " << (double)r / (double)fp1.size() << endl << endl;

    cout <<  "#########  MCM 2:   #########" << endl;
    Print_MCM_Partition(fp2, r);
    double LogE_fp2 = LogE_MCM(Kset, fp2, N, r);
    double LogL_fp2 = LogL_MCM(Kset, fp2, N, r);

    cout << "      Log-Evidence(MCM) = " << LogE_fp2 << "\t = " << LogE_fp2/((double) N)/log(2.) << " bits per datapoint \t" << endl;
    cout << "Max Log-Likelihood(MCM) = " << LogL_fp2 << "\t = " << LogL_fp2/((double) N)/log(2.) << " bits per datapoint \t" << endl;
    cout << "Average community size  = " << (double)r / (double)fp2.size() << endl << endl;


    cout << "#########  COMPARATIVE MEASURES   #########" << endl;
    double VOI = Var_of_Inf(fp1, fp2, r);
    double NMI = Norm_Mut_info(fp1, fp2, r);
    string istrue = is_subset(fp1, fp2) ? "Yes" : "No";

    cout << "Is MCM_1 \'subset\' of MCM_2 ?  : " << istrue << endl;
    cout << "Variation of Information      : " << VOI << endl;
    cout << "Normalized Mutual Information : " << NMI << endl;
    cout << "Difference in Log-Evidence    : LogE1-LogE2 = " << LogE_fp1 - LogE_fp2  << " = " << (LogE_fp1 - LogE_fp2)/((double) N)/log(2.) << " bits per datapoint \t" << endl << endl;
}


