#include <iostream>
#include <sstream>     // for stringstream
#include <fstream>
#include <vector>
#include <map>

using namespace std;

/******************************************************************************/
/******************   TOOL Functions from "tools.cpp"   ***********************/
/******************************************************************************/
string int_to_bstring(__int128_t bool_nb, unsigned int n);

/******************************************************************************/
/***********************    READ an MCM from a FILE  **************************/
/******************************************************************************/
map<unsigned int, __int128_t> read_MCM_fromfile(string Input_MCM_file, unsigned int r)
{
    map<unsigned int, __int128_t> Partition;

    string line, line2;
    __int128_t Op = 1;
    Op <<= r - 1;
    vector<int> comm;

    ifstream myfile(Input_MCM_file.c_str());
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            stringstream ss(line);
            while (getline(ss, line2, '\t'))
            {
                comm.push_back(stoi(line2));
            }
            Partition[comm[1]] += Op;
            Op >>= 1;

            comm.clear();
        }
        myfile.close();
    }
    return Partition;
}

/******************************************************************************/
/************************    PRINT MCM in TERMINAL  ***************************/
/******************************************************************************/
// r = number of variables in the dataset
void Print_MCM_Partition(map<unsigned int, __int128_t> partition, unsigned int r)
{
    map<unsigned int, __int128_t>::iterator it;
    int i = 1;

    for (it = partition.begin(); it != partition.end(); it++)
    {
        cout << "ICC " << i << ": \t " << int_to_bstring((*it).second, r) << endl;
        i++;
    }
    cout << endl;
}
