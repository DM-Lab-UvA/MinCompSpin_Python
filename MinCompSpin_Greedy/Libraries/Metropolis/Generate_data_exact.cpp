#include <list>
#include <fstream>
#include <map>
#include <bitset>
#include<cmath>

using namespace std;
/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "data.h"

/******************************************************************************/
/*********************   OPERATOR EVALUATED on a STATE   **********************/
/***********************   depend on the convention    ************************/
/******************************************************************************/

/****************************   !! Convention !!   ****************************/
/****     For Ising:      {1 in data <--> -1 in model}                     ****/
/****                     {0 in data <--> +1 in model}                     ****/
/******************************************************************************/
// s_i in {-1;1}   !! Convention !! {0,1} in data <-> {1,-1} in model
int Op_Ising(uint32_t Op, uint32_t state)
  {  return ( (bitset<n>(Op & state).count())%2 )?(-1):1;   } 

/******************************************************************************/
/************* PARTITION FUNCTION and PROBA of ALL STATES *********************/
/******************************************************************************/
// return the probability (not normalised) of all the states and compute the partition function

double* Probability_AllStates_Ising(list<Interaction> list_I, double *Z)  // Convention {-1;1} 
//!! Convention !!:  {1 in data <--> -1 in model}  and  {0 in data <--> 1 in model} 
{
  double H = 0; // -energy of the state
  int Op_s = 1; // value of the operator for the state s ; \in {-1; 1}

  list<Interaction>::iterator I;
  double* all_P = (double*)malloc((NOp_tot+1)*sizeof(double));
  (*Z) = 0;

  for (uint32_t state = 0; state <= NOp_tot; state++)
  {
    H=0;  // here H is (H = -Hamiltonian)
    for (I = list_I.begin(); I != list_I.end(); I++)
    {
      Op_s = Op_Ising((*I).Op, state);
      H += (*I).g * Op_s;
    }
    all_P[state] = exp(H);
    (*Z) += all_P[state];
  }

  return all_P;
}

/******************************************************************************/
/**************************     MODEL AVERAGES    *****************************/
/******************************************************************************/
void Model_averages_Ising_aux(double *P, double Z, list<Interaction> &list_I) 
{
  int Op_s = 1; // value of the operator for the state s ; \in {-1; 1}

  list<Interaction>::iterator I;

  for (I = list_I.begin(); I != list_I.end(); I++)
  {
    (*I).av_M = 0.;
    for (uint32_t state = 0; state <= NOp_tot; state++)
    {
      Op_s = Op_Ising((*I).Op, state);
      (*I).av_M += Op_s * P[state];
    }
    (*I).av_M = (*I).av_M / Z;
  }
}

void Model_averages_Ising(list<Interaction> &list_I)
{
// Compute un-normalized probability distribution and Partition function Z:
  double Z=0;
  double *P=Probability_AllStates_Ising(list_I, &Z);

  Model_averages_Ising_aux(P, Z, list_I);
}

/******************************************************************************/
/*************************     EMPIRICAL AVERAGES     *************************/
/******************************************************************************/
//Number of time an operator is equal to 1 ( = <phi> in the {0,1} representation )
unsigned int k1_op(map<uint32_t, unsigned int> Nset, uint32_t op)  // Complexity = O(|Nset|)
{
  unsigned int k1=0;
  map<uint32_t, unsigned int>::iterator it;  // iterator on Nset

  for (it = Nset.begin(); it!=Nset.end(); ++it)
    {    k1 += (bitset<n>( ((*it).first) & op ).count() % 2)*((*it).second);   }

  return k1;
}

double op_av_Ising(map<uint32_t, unsigned int> Nset, uint32_t op, unsigned int N)
{
  return ( ((double) N) - 2.*k1_op(Nset, op) ) / ((double) N); // [ [-1 * k1] + [+1 * (N-k1)] ] / N
}
/************************    Empirical averages all op    *************************/
void empirical_averages_Ising(map<uint32_t, unsigned int> Nset, list<Interaction> &list_I, unsigned int N) 
{
  list<Interaction>::iterator I;

  for (I = list_I.begin(); I != list_I.end(); I++)
  {
    (*I).av_D = op_av_Ising(Nset, (*I).Op, N);
  }
}

/******************************************************************************/
/**************************     SAMPLE DATASET    *****************************/
/******************************************************************************/
void Sample_dataset(list<Interaction> list_I, string output_filename, unsigned int N=1000)
{
// Compute un-normalized probability distribution and Partition function Z:
  double Z=0;
  double *p=Probability_AllStates_Ising(list_I, &Z);

// Compute Cumulative distribution:
  double *cumul = (double*)malloc((NOp_tot+1)*sizeof(double));

  cumul[0]=p[0]/Z;
  //cout << 0 << "\t " << cumul[0] << endl;

  for(unsigned int i=1; i <= NOp_tot; i++)
  {
        cumul[i] = cumul[i-1] + p[i]/Z;
        //cout << i << "\t " << cumul[i] << endl;
  }

// Randomly sampled "N" data points using the Cumulative:
  double eps=0;
  __int128_t j=0;

// OUTPUT FILE:
  fstream file(output_filename.c_str(), ios::out);

  for(int i=0; i<N; i++)
  {
    eps=(double)rand()/RAND_MAX;
    j=0;
    while(cumul[j]<eps && j < (1 << n))
    { j++;  }

    bitset<n> hi{ static_cast<unsigned long long>(j >> 64) },
        lo{ static_cast<unsigned long long>(j) },
        bits{ (hi << 64) | lo };
    file << bits << endl;
  }

  file.close();

  delete[] p;
  delete[] cumul;
}

/******************************************************************************/
/**************************     SAMPLE DATASET    *****************************/
/*************  and fill info on model + dataset    ***************************/
/********  i.e. compute Model averages and data averages    *******************/
/***************  and print them with the model information    ****************/
/******************************************************************************/
void Sample_dataset_AND_Print_ModelData_Info(list<Interaction>& list_I, string output_filename, unsigned int N=1000)
{
// Compute un-normalized probability distribution and Partition function Z:
  double Z=0;
  double *P=Probability_AllStates_Ising(list_I, &Z);

// Compute Cumulative distribution:
  double *cumul = (double*)malloc((NOp_tot+1)*sizeof(double));

  cumul[0]=P[0]/Z;

  for(int i=1; i <= NOp_tot; i++)
  {
    cumul[i] = cumul[i-1] + P[i]/Z;
  }

// Randomly sampled "N" data points using the Cumulative:
  double eps=0;
  uint32_t datapt=0;

// OUTPUT FILE:
  fstream file((OUTPUT_directory + output_filename).c_str(), ios::out);

// ***** data is also stored in Nset:  ********************************
  map<uint32_t, unsigned int> Nset; // Nset[mu] = #of time state mu appears in the data set

  for(unsigned int i=0; i<N; i++)
  {
    eps=rand();
    datapt=0;

    while(cumul[datapt]<eps)
    { datapt++;  }

    file << bitset<n>(datapt) << endl;
    Nset[datapt] += 1;
  }

  // Model averages:
  Model_averages_Ising_aux(P, Z, list_I);

  // Empirical averages (from the generated dataset):
  empirical_averages_Ising(Nset, list_I, N);

  file.close();
}

