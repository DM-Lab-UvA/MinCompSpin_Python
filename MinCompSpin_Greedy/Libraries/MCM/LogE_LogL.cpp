#include <iostream>
#include <cmath>       /* lgamma */
#include <map>
#include <vector>
#include <iomanip>

using namespace std;

/******************************************************************************/
/******************   TOOL Functions from "tools.cpp"   ***********************/
/******************************************************************************/
unsigned int Bitset_count(__int128_t bool_nb);
string int_to_bstring(__int128_t bool_nb, unsigned int n);

/******************************************************************************/
/*******************************    CONSTANTS     *****************************/
/******************************************************************************/
const __int128_t one128 = 1;

pair<bool, unsigned int> check_partition(map<unsigned int, __int128_t> Partition);

/******************************************************************************/
/************************ Build Kset for a single ICC  ************************/
/******************************************************************************/
map<__int128_t, unsigned int> build_Kset_ICC(vector<pair<__int128_t, unsigned int>> Kset, __int128_t Ai)
{
    map<__int128_t, unsigned int> Kset_ICC;
    __int128_t s;        // state

  //Build Kset_ICC:
    for (auto const& it : Kset)
    {
      s = ((it).first) & Ai;          // troncated state: take only the bits of s (=it.first) indicated by Ai
      Kset_ICC[s] += ((it).second);   // # of times s appears in the data set
    }

    return Kset_ICC;
}

/***********************************************************************************************************************/
/***********************************************************************************************************************/
/**************************************************   LOG-E   **********************************************************/
/***********************************************************************************************************************/
/***********************************************************************************************************************/

/******************************************************************************/
/*********  Log-Evidence (LogE) of an ICC part of a MCM   *********************/
/******************************************************************************/
// Compute the LogE of the ICC defined by Ai;
//    i.e. of the sub-part of an MCM identififed by Ai;
// this function doesn't account of the contribution to LogE due to the non-modeled spins (i.e. N*log(2) per spin)

double LogE_ICC(vector<pair<__int128_t, unsigned int>> Kset, __int128_t Ai, unsigned int N)
{
  map<__int128_t, unsigned int> Kset_ICC = build_Kset_ICC(Kset, Ai);  // Data reduced to the ICC
  unsigned int m = Bitset_count(Ai); // rank of the ICC

  double LogE = 0;
  unsigned int Ncontrol = 0; // for control
  unsigned int Ks = 0;

  map<__int128_t, unsigned int >::iterator it;

  for (it = Kset_ICC.begin(); it!=Kset_ICC.end(); ++it)
  {
    Ks = (it->second);  Ncontrol += Ks;
    if (Ks == 0) {cout << "problem Ks = 0 for some mu_m" << endl; }
    LogE += lgamma(Ks + 0.5);
  }  
  if (Ncontrol != N) { cout << "Error \'LogE_ICC\' function: Ncontrol != N" << endl;  }

  //  return LogE - GeomComplexity_ICC(m) - lgamma( (double)( N + (one128 << (m-1)) ) );
  return LogE + lgamma((double)( one128 << (m-1) )) - (Kset_ICC.size()/2.) * log(M_PI) - lgamma( (double)( N + (one128 << (m-1)) ) ); 
}

/******************************************************************************/
/****************************   LogE of a MCM   *******************************/
/******************************************************************************/
double LogE_MCM(vector<pair<__int128_t, unsigned int>> Kset, map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r)
{
  //if (!check_partition(Partition)) {cout << "Error, the argument is not a partition." << endl; return 0;  }

  //else
  //{
    double LogE = 0; 
    unsigned int rank = 0;
    map<unsigned int, __int128_t>::iterator Part;

    for (Part = Partition.begin(); Part != Partition.end(); Part++)
    {
      LogE += LogE_ICC(Kset, (*Part).second, N);
      rank += Bitset_count((*Part).second);
    }  
    return LogE - ((double) (N * (r-rank))) * log(2.);
  //}
}

double LogE_MCM_infoICC(vector<pair<__int128_t, unsigned int>> Kset, map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r)
{
    double LogE = 0, LogE_SCM = 0, LogE_CM = 0, eta = 0; 
    unsigned int rank = 0;
    double Nd = (double) N;
    cout << setprecision(3) << fixed;

    // Complete model:
    __int128_t Part_CM = (one128 << r)-1;
    LogE_CM = LogE_ICC(Kset, Part_CM, N);
    cout << "Complete model: " << int_to_bstring(Part_CM, r) << "\t  LogE(CM) = " << LogE_CM << "\t = " << LogE_CM/Nd/log(2.) << " bits per datapoint" << endl << endl;

    // Each ICC:
    map<unsigned int, __int128_t>::iterator Part;
    int i=1;
    for (Part = Partition.begin(); Part != Partition.end(); Part++)
    {
      LogE_SCM = LogE_ICC(Kset, (*Part).second, N);
      LogE += LogE_SCM;
      rank += Bitset_count((*Part).second);
      //eta = ((LogE_SCM/Nd/log(2.)) + r) / (r + (LogE_CM/Nd/log(2.)));
      cout << "ICC " << i << ": \t " << int_to_bstring((*Part).second, r) << "\t LogE(ICC) = " << LogE_SCM << "\t = " << LogE_SCM/Nd/log(2.) << " bits per datapoint \t" << endl;
      i++;
    }  
    cout << endl;

    LogE -= ((double) (N * (r-rank))) * log(2.);

    cout << "Log-Evidence(MCM) = " << LogE << "\t = " << LogE/((double) N)/log(2.) << " bits per datapoint \t" << endl;

    return LogE;
}

/***********************************************************************************************************************/
/***********************************************************************************************************************/
/**************************************************   LOG-L   **********************************************************/
/***********************************************************************************************************************/
/***********************************************************************************************************************/

/******************************************************************************/
/**************** Log-likelihood (LogL) of a Complete Model  ******************/
/******************************************************************************/
// Compute the log-likelihood of a Complete Model on Kset:

double LogL_CM(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N)
{
  double LogL = 0;
  double Nd = N;

  for (auto& it : Kset)
  {
    LogL += ((it.second) * log((double) (it.second) / Nd) ); 
  } 

  return LogL;
}

/******************************************************************************/
/*********  Log-Likelihood (LogL) of an ICC part of a MCM   *******************/
/******************************************************************************/
// Compute the LogL of the ICC defined by Ai;
//    i.e. of the sub-part of an MCM identififed by Ai;
// this function doesn't account of the contribution to LogL due to the non-modeled spins (i.e. N*log(2) per spin)

double LogL_ICC(vector<pair<__int128_t, unsigned int>> Kset, __int128_t Ai, unsigned int N)
{
  map<__int128_t, unsigned int> Kset_ICC = build_Kset_ICC(Kset, Ai);

  double LogL = 0;

  map<__int128_t, unsigned int>::iterator it;
  unsigned int Ncontrol = 0; // for control
  unsigned int Ks = 0;
  double Nd = N;

  for (it = Kset_ICC.begin(); it!=Kset_ICC.end(); ++it)
  {
    Ks = (it->second);  Ncontrol += Ks;
    if (Ks == 0) {cout << "problem Ks = 0 for some mu_m" << endl; }
    LogL += (Ks * log((double) Ks / Nd) );
  }
  if (Ncontrol != N) { cout << "Error in function 'LogL_ICC': Ncontrol != N" << endl;  }

  return LogL;
}

/******************************************************************************/
/******************** Log-likelihood (LogL) of a MCM  *************************/
/******************************************************************************/

double LogL_MCM(vector<pair<__int128_t, unsigned int>> Kset, map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r)
{
  //if (!check_partition(Partition)) {cout << "Error, the argument is not a partition." << endl; return 0;  }

  //else
  //{
    double LogL = 0; 
    unsigned int rank = 0;

    map<unsigned int, __int128_t>::iterator Part;
    for (Part = Partition.begin(); Part != Partition.end(); Part++)
    {
      LogL += LogL_ICC(Kset, (*Part).second, N);
      rank += Bitset_count((*Part).second);
    }  
    return LogL - ((double) (N * (r-rank))) * log(2.);
  //}
}

double LogL_MCM_infoICC(vector<pair<__int128_t, unsigned int>> Kset, map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r)
{
  //if (!check_partition(Partition)) {cout << "Error, the argument is not a partition." << endl; return 0;  }

  //else
  //{
    double LogL = 0, LogL_SCM = 0, LogL_C = 0, eta_a = 0, eta = 0, eta_denominator = 0;
    unsigned int rank = 0, rank_a = 0;
    double Nd = (double) N;
    cout << setprecision(3) << fixed;

    // Complete model:
    __int128_t Part_CM = (one128 << r) - 1;
    LogL_C = LogL_CM(Kset, N);
    eta_denominator = LogL_C + ((double) (r * N)) * log(2.);

    cout << "Maximum amount of information that can be extracted = LogL(CM) - LogL(0) = ";
    cout << eta_denominator << " = " << eta_denominator/Nd/log(2.) << " bits per datapoint" << endl << endl;

    cout << "Complete model: " << int_to_bstring(Part_CM, r) << "\t LogL(CM) = " << LogL_C << "\t = " << LogL_C/Nd/log(2.) << " bits per datapoint" << endl;
    cout << "Empty model   : " << int_to_bstring(0, r) << "\t  LogL(0) = " << -((double) (r * N)) * log(2.) << "\t = " << -(int) r << " bits per datapoint" << endl << endl;

    // Each ICC:
    map<unsigned int, __int128_t>::iterator Part;
    int i=1;
    for (Part = Partition.begin(); Part != Partition.end(); Part++)
    {
        LogL_SCM = LogL_ICC(Kset, (*Part).second, N);
        rank_a = Bitset_count((*Part).second);
        eta_a = (LogL_SCM + ((double) (rank_a * N)) * log(2.)) / eta_denominator;

        LogL += LogL_SCM;
        rank += Bitset_count((*Part).second);
        eta += eta_a;

        cout << "ICC " << i << ": \t " << int_to_bstring((*Part).second, r) << "\t LogL(ICC) = " << LogL_SCM << "\t = " << LogL_SCM/Nd/log(2.) << " bits per datapoint; \t";
        cout << "eta(ICC) = " << eta_a << endl;
        i++;
    }  

    // final LogL(MCM)
    LogL -= ((double) (N * (r-rank))) * log(2.);

    cout << endl;
    cout << "Max-LogL(MCM) = " << LogL << "\t = " << LogL/((double) N)/log(2.) << " bits per datapoint \t" << endl;
    cout << "Fraction of deviance explained: eta(MCM) = " << eta << endl;

    return LogL;
  //}
}
