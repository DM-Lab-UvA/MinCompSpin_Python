#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <map>
#include <vector>
//#include <cstring>

using namespace std;


/******************************************************************************/
/***************************   Constant variables   ***************************/
/******************************************************************************/
const __int128_t one128 = 1;

/******************************************************************************/
/******************   TOOL Functions from "tools.cpp"   ***********************/
/******************************************************************************/
string int_to_bstring(__int128_t bool_nb, unsigned int n);
unsigned int Bitset_count(__int128_t bool_nb);


/******************************************************************************/
/***********************     READ DATA FILE    ********************************/
/******************************************************************************/
/**************    READ DATA and STORE them in Nset    ************************/

vector<pair<__int128_t, unsigned int>> read_datafile(unsigned int *N, std::string filename, unsigned int r) // O(N)  where N = data set size
{
    cout << "Read the dataset: " << filename << endl;
    cout << "Number of variables to read: n = " << r << endl;

    string line, line2;     char c = '1';
    __int128_t state = 0, Op;
    (*N) = 0;            // N = dataset sizes

// ***** Store data in Nset_map:  **********************************************
    map<__int128_t, unsigned int> Nset_map; // Nset[mu] = #of time state mu appears in the data set

    ifstream myfile (filename.c_str());
    if (myfile.is_open())
    {
        while ( getline (myfile,line))
        {
            line2 = line.substr (0,r);          //take the r first characters of line
            Op = one128 << (r - 1);
            state = 0;
            for (auto &elem: line2)     //convert string line2 into a binary integer
            {
                if (elem == c) { state += Op; }
                Op = Op >> 1;
            }
            Nset_map[state] += 1;
            (*N)++;
        }
        myfile.close();
    }
    else
    {
        cout << endl << "                     ########## Unable to open file ##########" << endl << endl;
    }
    //cout << "\t\t data size N = " << (*N) << endl;
        
// ***** Convert map to a vector:  for faster reading later on ********************************
    vector<pair<__int128_t, unsigned int>> Nset(Nset_map.size());    //Nset.resize(Nset_map.size());

    int i=0;
    for (auto& my_pair : Nset_map)
    {
        Nset[i]=my_pair;
        i++;
    }
    return Nset;
}


/******************************************************************************/
/**************************     PRINT Nset   **********************************/
/******************************************************************************/
/*
void Print_File_Nset(map<__int128_t, unsigned int> Nset, unsigned int N, unsigned int r, string OUTPUTfilename)
// map.second = nb of time that the state map.first appears in the data set
{
  map<__int128_t, unsigned int>::iterator it;
  int Ncontrol = 0;
  __int128_t one128 = 1;

  fstream file(OUTPUTfilename.c_str(), ios::out);
  file << "#N = " << N << endl;
  file << "#Total number of accessible states = 2^r-1 = 2^(" << r << ") - 1" << endl;
  file << "#Number of visited states, Nset.size() = " << Nset.size() << endl;
  file << "#" << endl;
  file << "#1: state \t #2: nb of pts in state \t #3: Pba state" << endl;

  for (it = Nset.begin(); it!=Nset.end(); ++it)
  {
    file << int_to_bstring((*it).first, r) << " => " << (*it).second; // << endl;
    file << "  \t  P = " << ((*it).second) / (float) N << endl;
    Ncontrol += (*it).second;
  }

  if (Ncontrol != N) { cout << "Error function \'read_Nset\': Ncontrol != N" << endl;  }

  file.close();
}
*/


/******************************************************************************/
/*********************     CHANGE of BASIS: one datapoint  ********************/
/******************************************************************************/
// Given a choice of a basis (defined by the m-basis list) --> returns the new m-state (i.e. state in the new m-basis)
// Rem: must have m <= n 

// mu = old state
// final_mu = new state

__int128_t transform_mu_basis(__int128_t mu, list<__int128_t> basis)
{
  __int128_t un_i = 1, proj;
  __int128_t final_mu = 0;

  list<__int128_t>::iterator phi_i;

  for(phi_i = basis.begin(); phi_i != basis.end(); ++phi_i)
  {
    proj = (*phi_i) & mu;
    if ( (Bitset_count(proj) % 2) == 1) // odd number of 1, i.e. sig_i = 1
    {
      final_mu += un_i;
    }
    un_i = (un_i << 1);
  }

  return final_mu;
}

/******************************************************************************/
/************************** K_SET *********************************************/
/******************************************************************************/
// Build Kset for the states written in the basis of the m-chosen independent 
// operator on which the SC model is based:

vector<pair<__int128_t, unsigned int>> build_Kset(vector<pair<__int128_t, unsigned int>> Nset, list<__int128_t> Basis)
// sig_m = sig in the new basis and cut on the m first spins 
// Kset[sig_m] = #of time state mu_m appears in the data set
{
    map<__int128_t, unsigned int > Kset_map;
    __int128_t sig_m;    // transformed state and to the m first spins

// ***** Build Kset: *************************************************************************************
    cout << endl << "--->> Build Kset..." << endl;

    for (auto const& it : Nset)
    {
        sig_m = transform_mu_basis((it).first, Basis); // transform the initial state s=(it).first into the new basis
        Kset_map[sig_m] += ((it).second); // ks = (it).second = number of time state s appear in the dataset
    }
    cout << endl;

// ***** Convert map to a vector:  for faster reading later on ********************************************
    vector<pair<__int128_t, unsigned int>> Kset(Kset_map.size());

    int i=0;
    for (auto& my_pair : Kset_map)
    {
        Kset[i]=my_pair;
        i++;
    }

    return Kset;
}


/******************************************************************************/
/****************************   REDUCE K_SET   ********************************/
/******************************************************************************/
// Remove all the states that occur less than a chosen number K of times

vector<pair<__int128_t, unsigned int>> Reduce_Kset(vector<pair<__int128_t, unsigned int>> Kset_Vect, unsigned int K, unsigned int *N_new)
{
    cout << endl << "States removed from Kset:" << endl;

    unsigned int* counter = (unsigned int*)malloc((K+1)*sizeof(unsigned int));
    for(int i=0; i<=K; i++)
        {   counter[i]=0;   }

    vector<pair<__int128_t, unsigned int>> Kset_reduced;
    for (auto const& it : Kset_Vect)
    {
        if (it.second <= K)
        {    
            counter[it.second]++;
        }
        else 
        { 
            Kset_reduced.push_back(it);
        }
    }

    for(int i=1; i<=K; i++)
    {   
        cout << "\t -- there are " << counter[i] << " states that appear ks = " << i << " times;" << endl;   
        (*N_new) -= i*counter[i];   // reduced the total number of states by the number of states removed
    }
    cout << endl;

    return Kset_reduced;
}




