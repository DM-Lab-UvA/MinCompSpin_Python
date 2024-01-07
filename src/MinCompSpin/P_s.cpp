#include <iostream>
#include <map>
#include <list>
#include <vector>
#include <fstream>

using namespace std;


/******************************************************************************/
/***************************   ADD OUTPUT FOLDER    ***************************/
/******************************************************************************/
//#include "../../includes/default_datafile.h"

string OutputFile_Add_Location(string filename);  // defined in main.cpp and in py_binding.cpp
/*
string OutputFile_Add_Location(string filename)
{
    return (OUTPUT_directory + filename);
}
*/

/******************************************************************************/
/***************************   Constant variables   ***************************/
/******************************************************************************/
const __int128_t one128 = 1;

/********************************************************************/
/*********************    Proba Structure    ************************/
/********************************************************************/
struct Proba {
    __int128_t s;         // state in the original basis
    __int128_t sig;       // state in the new basis

    double P_D_s = 1.;    // empirical probability of s  
    double P_D_sig = 1.;  // empirical probability of sig --> this should be the same then s if r=n, i.e. if the MCM models all the n spins
    double P_MCM = 1.;  // model probability of s 
};

/******************************************************************************/
/******************   TOOL Functions from "tools.cpp"   ***********************/
/******************************************************************************/
unsigned int Bitset_count(__int128_t bool_nb);
string int_to_bstring(__int128_t bool_nb, unsigned int n);

//check if *Partition* is an actual partition of the basis elements, 
// i.e., that no basis element appears in more than 1 part of the partition.
// i.e., that each basis element only appears in a single part of the partition.
pair<bool, unsigned int> check_partition(map<unsigned int, __int128_t> Partition); 

/******************************************************************************/
/********************   Return Kset over a chosen ICC    **********************/
/******************************************************************************/
map<__int128_t, unsigned int> build_Kset_ICC(vector<pair<__int128_t, unsigned int>> Kset, __int128_t Ai); // Ai = integer indicated the spins included in the ICC

/******************************************************************************/
/*****************   Compute the contribution to P_MCM(s)   *******************/
/*****************  due to the sub-CM defined by Kset_ICC   *******************/
/******************************************************************************/

void update_proba_MCM(map<__int128_t, Proba> &all_P, map<__int128_t, unsigned int> Kset_ICC, __int128_t Ai, unsigned int N)
{
  map<__int128_t, Proba>::iterator it_P;

  __int128_t s, sig;        // states
  unsigned int ks=0;      // number of time state s appear in the dataset
  double Nd = (double) N;

//Contribution of the basis "ba" of spins to the model probability P_MCM(s):
  for (it_P = all_P.begin(); it_P!=all_P.end(); ++it_P)
  {
    s = it_P->first;      // initial state s 
    sig = s & Ai;         // troncated state: take only the bits indicated by Ai
    all_P[s].P_MCM *= Kset_ICC[sig]/Nd;
  }
}

/******************************************************************************/
/******************************************************************************/
/*************      Compute and Print the model probability     ***************/
/*********************  for the MCM constructed from Kset  ********************/
/************************   with the given partition   ************************/
/******************************************************************************/
/******************************************************************************/

/*************      Compute the model probabilities     ***************/

// This function can be used directly on the original basis, by replacing Kset by Nset:
map<__int128_t, Proba> P_sig(vector<pair<__int128_t, unsigned int>> Kset, map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r) // Probabilities in the sigma basis
{
  // Fill in the data probability:
  map<__int128_t, Proba> all_P;
  double Nd = (double) N;

  __int128_t s;        // state
  //unsigned int ks=0; // number of time state s appear in the dataset

  // Check partition:
  pair<bool, unsigned int> Is_partition = check_partition(Partition);
  unsigned int rank = Is_partition.second;

  if (!Is_partition.first) {cout << "Error, the argument is not a partition: the function returned an empty map for P[s]." << endl; }
  else
  { 
    double pre_factor = 1./((double) (one128 << (r-rank))); 

    for (auto const& it : Kset)
    {   
      s = (it).first;      // initial state s 
      //ks = it->second;    // # of times s appears in the data set

      all_P[s].P_D_s = ((double) ((it).second))/Nd;  // (it).second = ks = # of times s appears in the data set
      all_P[s].P_MCM = pre_factor;
    }

    // Compute the Kset over each ICC: Kset_ICC:
    map<__int128_t, unsigned int> Kset_ICC;
    map<unsigned int, __int128_t>::iterator Part;

    for (Part = Partition.begin(); Part != Partition.end(); Part++)
    {
      Kset_ICC = build_Kset_ICC(Kset, (*Part).second);         // (*Part).second) = Ai = integer indicated the spin elements included in b_a
      update_proba_MCM(all_P, Kset_ICC, (*Part).second, N);
      Kset_ICC.clear();
    }  
  }

  return all_P;
}

/*************      Print the model probabilities     ***************/

void PrintFile_StateProbabilites_NewBasis(vector<pair<__int128_t, unsigned int>> Kset, map<unsigned int, __int128_t> MCM_Partition, unsigned int N, unsigned int r, string filename = "Result")
{
  // Probabilities in the sigma basis:
  map<__int128_t, Proba> P_all = P_sig(Kset, MCM_Partition, N, r);
  map<__int128_t, Proba>::iterator it_P;

  string Psig_filename = filename + "_DataVSMCM_Psig.dat";

  cout << "--->> Print the state probabilities P(sig) in the file: \'" << Psig_filename << "\'" << endl << endl;

  //fstream file_P_sig((OUTPUT_directory + Psig_filename), ios::out);
  fstream file_P_sig(OutputFile_Add_Location(Psig_filename), ios::out);

  file_P_sig << "## 1:sig \t 2:P_D(sig) \t 3:P_MCM(sig)" << endl;

  for (it_P = P_all.begin(); it_P!=P_all.end(); ++it_P)
  {   
    file_P_sig << int_to_bstring((it_P->first), r) << "\t " << (it_P->second).P_D_s << "\t " << (it_P->second).P_MCM << endl;
  }

  file_P_sig.close();
}

/******************************************************************************/
/******************************************************************************/
/*******************      Compute the model probability     *******************/
/**************************  for the MCM constructed   ************************/
/***************   in a given basis, with a given partition   *****************/
/******************************************************************************/
/******************************************************************************/
__int128_t transform_mu_basis(__int128_t mu, list<__int128_t> basis);
vector<pair<__int128_t, unsigned int>> build_Kset(vector<pair<__int128_t, unsigned int>> Nset, list<__int128_t> Basis);


//// The partition must be guven in the new basis:
map<__int128_t, Proba> P_s(vector<pair<__int128_t, unsigned int>> Nset, list<__int128_t> Basis, map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r) // Probabilities in the sigma basis
{
  double Nd = (double) N;

  // Fill in the data probability:
  map<__int128_t, Proba> all_P;

  if (!check_partition(Partition).first) {cout << "Error, the argument is not a partition: the function returned an empty map for P[s]." << endl; }
  else
  { 
    // Build Kset:
    //map<__int128_t, unsigned int > Kset_map;
    __int128_t s;        // initial state
    __int128_t sig_m;    // transformed state and to the m first spins
    unsigned int ks=0; // number of time state s appear in the dataset

// ***** Compute P[s] from the original data: *************************************************************************************
    for (auto const& it : Nset)
    {
      s = (it).first;       // state s
      ks = (it).second;     // # of times s appears in the data set
      sig_m = transform_mu_basis(s, Basis);

      // Fill in P[s] empirical and value of transformed state:
      all_P[s].P_D_s = ((double) ks)/Nd;
      all_P[s].sig = sig_m;  

      // Fill in Kset:
      //Kset[sig_m] += ks;
    }

// ***** Build Kset: ***************************************************************************************************************
    vector<pair<__int128_t, unsigned int>> Kset = build_Kset(Nset, Basis);

// ***** Compute P_MCM[s] for the chosen MCM based on the new basis (using "Partition" on Kset): ***********************************
    //cout << "--->> Compute P[s] for the chosen MCM..." << endl << endl;
    map<__int128_t, Proba> all_P_sig = P_sig(Kset, Partition, N, r);   // Compute P[s] for the MCM in the new basis

  // Report the values of P_MCM in the original all_P:
    map<__int128_t, Proba>::iterator it_P;
    for (it_P = all_P.begin(); it_P!=all_P.end(); ++it_P)           // Report P[s] for the MCM in the original basis
    {
      sig_m = (it_P->second).sig;
      (it_P->second).P_MCM = all_P_sig[sig_m].P_MCM;
      (it_P->second).P_D_sig = all_P_sig[sig_m].P_D_s;  // not necessary a probability distribution anymore
    }
  
    all_P_sig.clear();
    Kset.clear();
  }

  return all_P;  
}

/******************************************************************************/
/*****************      PRINT FILE: INFO about an MCM     *********************/
/******************************************************************************/

void PrintFile_MCM_Info(list<__int128_t> Basis, map<unsigned int, __int128_t> MCM_Partition, unsigned int r, string filename = "Result")
{
  //***** PRINT BASIS: 
  //fstream file_MCM_info((OUTPUT_directory + filename + "_MCM_info.dat"), ios::out);
  fstream file_MCM_info(OutputFile_Add_Location(filename + "_MCM_info.dat"), ios::out);

  file_MCM_info << "## sig_vec = states in the chosen new basis (ideally the best basis), defined by the basis operators:" << endl;
  int i = 1;
  for (list<__int128_t>::iterator it = Basis.begin(); it != Basis.end(); it++)
  {
    file_MCM_info << "##\t sig_" << i << " = " << int_to_bstring((*it), r) << endl; i++;
  } file_MCM_info << "##" << endl;

  // Print info about the model -- Print MCM:
  file_MCM_info << "## The MCM Partition is defined on the sig_vec basis by the following Parts:" << endl;

  //***** PRINT MCM: 
  i = 1;
  for (map<unsigned int, __int128_t>::iterator it = MCM_Partition.begin(); it != MCM_Partition.end(); it++)
  {    
    __int128_t Part = (*it).second;

    file_MCM_info << "##\t MCM_Part_" << i << " = " << int_to_bstring(Part, r) << endl; i++;
  }
  file_MCM_info << "##" << endl;

  file_MCM_info.close();
}


/******************************************************************************/
/*************      Print the model probabilities in a file     ***************/
/******************************************************************************/
void PrintFile_StateProbabilites_OriginalBasis(vector<pair<__int128_t, unsigned int>> Nset, list<__int128_t> Basis, map<unsigned int, __int128_t> MCM_Partition, unsigned int N, unsigned int r, string filename = "Result")
{
  // Compute all the state probabilities:
  map<__int128_t, Proba> P_all = P_s(Nset, Basis, MCM_Partition, N, r);

  double *Pk_D = (double *)malloc((r+1)*sizeof(double)); 
  double *Pk_MCM = (double *)malloc((r+1)*sizeof(double)); 

  unsigned int k = 0;
  for(k=0; k<=r; k++)
  {
    Pk_D[k] = 0;
    Pk_MCM[k] = 0;
  }

  string Ps_filename = filename + "_DataVSMCM_Ps.dat";
  string Pk_filename = filename + "_DataVSMCM_Pk.dat";

  cout << "--->> Print information about the MCM in the file: \'" << filename << "_MCM_info.dat\'" << endl;
  cout << "--->> Print the state probabilities P(s) in the file: \'" << Ps_filename << "\'" << endl;
  cout << "--->> Print the probability of a state with k \'+1\' bits: \'" << Pk_filename << "\'" << endl << endl;

  //***** Print info about the model -- Print Basis and MCM:  **************/
  PrintFile_MCM_Info(Basis, MCM_Partition, r, filename);

  //***** Print P(s):  *****************************************************/
  __int128_t s;
  //fstream file_Ps((OUTPUT_directory + Ps_filename), ios::out);
  fstream file_Ps(OutputFile_Add_Location(Ps_filename), ios::out);

  file_Ps << "## s = states in the original basis" << endl;
  file_Ps << "## sig = states in the chosen new basis (ideally the best basis)" << endl;
  file_Ps << "## The chosen \'sig\'-basis and the chosen MCM are printed in the file " << filename << "_MCM_info.dat" << endl;

  // Print P(s): 
  file_Ps << "## " << endl;
  file_Ps << "## 1:s \t 2:P_D(s) \t 3:P_MCM(s) \t 4:sig" << endl;

  for (map<__int128_t, Proba>::iterator it_P = P_all.begin(); it_P!=P_all.end(); ++it_P)
  {   
    s = it_P->first;
    
    file_Ps << int_to_bstring(s, r) << "\t" <<  (it_P->second).P_D_s << "\t" << (it_P->second).P_MCM << "\t" << int_to_bstring((it_P->second).sig, r) << endl;

    k = Bitset_count(s);
    Pk_D[k] += (it_P->second).P_D_s;      // P[k] in the data
    Pk_MCM[k] += (it_P->second).P_MCM;    // P[k] from the MCM
  }
  file_Ps.close();

  //***** Print P(k):   ***************************************************/
  //fstream file_Pk((OUTPUT_directory + Pk_filename), ios::out);
  fstream file_Pk(OutputFile_Add_Location(Pk_filename), ios::out);

  file_Pk << "## 1:k \t 2:P_D(k) \t 3:P_MCM(k)" << endl;

  for(k=0; k<=r; k++)
  {
    file_Pk << k << "\t" << Pk_D[k] << "\t" << Pk_MCM[k] << endl;
  }
  file_Pk.close();
}


