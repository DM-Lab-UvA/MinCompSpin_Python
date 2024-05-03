#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

/******************************************************************************/
/******************   TOOL Functions from "tools.cpp"   ***********************/
/******************************************************************************/
//string int_to_bstring(__int128_t bool_nb, unsigned int n);

/******************************************************************************/
/************************   Function from "LogE.cpp"   ************************/
/******************************************************************************/
double LogE_ICC(vector<pair<__int128_t, unsigned int>> Kset, __int128_t Ai, unsigned int N);
double LogE_MCM(vector<pair<__int128_t, unsigned int>> Kset, map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r);

/******************************************************************************/
/***************    FIND THE BEST MCM: GREEDY PROCEDURE    ********************/
/******************************************************************************/
unsigned int mat_loc(unsigned int i, unsigned int j)  // for i < j
{
    return j * (j - 1) / 2 + i;
}

//void PrintFile_MCM_FullMerge(map<unsigned int, __int128_t> MCM_Partition, unsigned int r, fstream &file);

/******************************************************************************/
/***************    FIND THE BEST MCM: GREEDY PROCEDURE    ********************/
/***************          Starting from chosen MCM_0       ********************/
/******************************************************************************/
// r = number of communities in MCM_0

map<unsigned int, __int128_t> MCM_GreedySearch_MCM0(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, map<unsigned int, __int128_t> MCM_0, bool print_info = true) //, bool Greedy_Full_merging = false)
{
    bool Greedy_Full_merging = false;
    fstream file_MCMs;  // used only if full merging: print all MCMs along merging

    double Nd = (double) N;
    if(Greedy_Full_merging)
    {
        cout << "**** Performing a full merging procedure until fully merged model: print merging path ****" << endl;
        print_info = true;

        file_MCMs.open("MCMs_FullMerge.dat", ios::out);
        file_MCMs << "## Successive MCM Partitions found along the Greedy hierarchical merging search." << endl;
        file_MCMs << "## \tMCM are given in the specified chosen basis:" << endl;
    }

    // ** Create MCM_Best and initialise it to MCM_0:   
    map<unsigned int, __int128_t> Best_MCM(MCM_0);

    // ** Number of ICCs in MCM_Best:
    unsigned int rank = Best_MCM.size();
    double LogE_tot_best = 0;

    // ** Store values of LogE for each ICC:
    map<unsigned int, double> Best_MCM_LogE_ICC;
    double LogE_tot = 0;

    // ** Final Best MCM recorded along the merging path:
    map<unsigned int, __int128_t> The_Best_MCM; 

/*
    for (unsigned int i = 0; i < r; i++)
    {
        Best_MCM_LogE_ICC[i] = LogE_ICC(Kset, Best_MCM[i], N);
        LogE_tot += Best_MCM_LogE_ICC[i];
    }
*/
    for (auto const& it : Best_MCM)
    {
        Best_MCM_LogE_ICC[it.first] = LogE_ICC(Kset, it.second, N);
        LogE_tot += Best_MCM_LogE_ICC[it.first];
    }
    LogE_tot_best = LogE_tot;

    // ** Create table of LogE values for all pairs of merged ICCs:
    cout << "\t **** Filling in the rxr-triangular Matrix ****" << endl;
    double* logE_merged = (double*)malloc(rank * (rank - 1) / 2 * sizeof(double));

    // logE_merged[mat_loc(i,j)] will contain the values of LogE of merged(ICCs i and j)
    for (unsigned int i = 0; i < rank - 1; i++)
    {
        for (unsigned int j = i + 1; j < rank; j++)
        {
            logE_merged[mat_loc(i, j)] = 0;  // initalized to 0
        }
    }
    cout << endl;

    // buffers:
    map<unsigned int, __int128_t> Buffer_ICC_unmerged;
    double Diff_LogE_best = 0, LogE_merged_best = 0, LogE_unmerged = 0;

    // Best:
    map<unsigned int, __int128_t>::iterator it1, it2;
    bool stop = false;

    unsigned int i_keep, i_erase, i_mat;   // define variables for merging

    // ********* HIERARCHICAL MERGING ****************
    cout << "\t **** Start Merging ****" << endl;

    streamsize ss = cout.precision();
    cout << fixed << setprecision(2); // precision of LogE: print 2 digits after "."
    cout << "   Start: \t\t Nb of communities = " << Best_MCM.size();
    cout << "\t LogE = " << LogE_tot;
    cout << " = " << LogE_tot/log(2.)/Nd << " bits/datapoint"<< endl; // << "\t Check: " << LogE_MCM(Kset, Best_MCM, N, r) << endl;

//    if(Greedy_Full_merging) 
//        {  PrintFile_MCM_FullMerge(Best_MCM, r, file_MCMs); }

    int iteration = 0; 
    while (Best_MCM.size() > 1 && (!stop))
    {
        Diff_LogE_best = 0;
        LogE_unmerged = 0; //LogE_merged = 0, 

        // ****** Find the best pair to merge: *********
        bool start_loop = true;
        for (it1 = Best_MCM.begin(); it1 != (--Best_MCM.end()); ++it1)
        {
            //cout << ">> test for it1 = " << (*it1).first << endl;

            for (it2 = next(it1); it2 != Best_MCM.end(); ++it2)
            {
                //cout << "test for it2 = " << (*it2).first << "\t ";

                // LogE of unmerged ICCs it1 and it2:
                LogE_unmerged = Best_MCM_LogE_ICC[(*it1).first] + Best_MCM_LogE_ICC[(*it2).first];

                // LogE of merged ICC (it1 + it2):
                i_mat = mat_loc((*it1).first, (*it2).first);
                //cout << "i_mat = " << i_mat << "\t ";
                if (logE_merged[i_mat] == 0)
                {
                    //cout << "is 0: ";
                    logE_merged[i_mat] = LogE_ICC(Kset, (*it1).second + (*it2).second, N);
                }
                //cout << "LogE[i_mat] = " << logE_merged[i_mat] << "\t";
                //cout << "LogE_unmerged = " << LogE_unmerged  << endl;

                // if merged ICC is better, then save it:
                //if ( ((logE_merged[i_mat] - LogE_unmerged) > Diff_LogE_best) || ( (it1 == Best_MCM.begin()) && (it2 == next(Best_MCM.begin())) ) )
                if ( ((logE_merged[i_mat] - LogE_unmerged) > Diff_LogE_best) || start_loop )
                {
                    //cout << "Is larger \t";
                    Buffer_ICC_unmerged.clear();
                    Buffer_ICC_unmerged.insert(*it1);
                    Buffer_ICC_unmerged.insert(*it2);

                    Diff_LogE_best = logE_merged[i_mat] - LogE_unmerged;
                    LogE_merged_best = logE_merged[i_mat];
                    //cout << "Diff LogE best = " << Diff_LogE_best << " : " << (*it1).first << "\t " << (*it2).first << endl;
                }
                start_loop = false;
            }
        }

        // ********* STOPPING CRITERIA **************:
        if (Diff_LogE_best <= 0 && !Greedy_Full_merging)  // there is no improvement of Log-E
        { 
            stop = true; 
            if (print_info) {cout << "   Stop: No more improvement of LogE" << endl; } 
            The_Best_MCM = Best_MCM;
            continue;
        }
        else if (Diff_LogE_best <= 0 && LogE_tot > LogE_tot_best)
        {
             The_Best_MCM = Best_MCM; 
             LogE_tot_best = LogE_tot;
        }

        // ********* PERFORM MERGING *****:
        it1 = Buffer_ICC_unmerged.begin();  i_keep = (*it1).first;   //  store the merged ICC in lowest index: i_keep
        it1++;  i_erase = (*it1).first; //  erase the ICC in the largest index: i_erase

        Best_MCM[i_keep] = Best_MCM[i_keep] + Best_MCM[i_erase];
        Best_MCM.erase(i_erase);

        Best_MCM_LogE_ICC[i_keep] = LogE_merged_best;
        Best_MCM_LogE_ICC.erase(i_erase);

        LogE_tot += Diff_LogE_best;

        // ******* For all the indices "i" that are still stored in Best_MCM:
        // ************ reset the values of logE_merged between "i" and "i_keep" to "0"
        for (it1 = Best_MCM.begin(); it1 != Best_MCM.end(); it1++)
        {
            if ((*it1).first > i_keep)
            {
                logE_merged[mat_loc(i_keep, (*it1).first)] = 0;
            }
            else if ((*it1).first < i_keep)
            {
                logE_merged[mat_loc((*it1).first, i_keep)] = 0;
            }
        }

        iteration++;
        if (print_info)
        {   
            cout << "   Iteration " << setw(3) << setfill(' ') << left << iteration << " done: \t Nb of communities = " << setw(3) << setfill(' ') << left << Best_MCM.size();
            cout << "\tLogE = " << LogE_tot; //<< "\t Check: "<< LogE_MCM(Kset, Best_MCM, N, r) << endl;
            cout << " = " << LogE_tot/log(2.)/Nd << " bits/dpt";
            cout << "\t D(LogE) = " << Diff_LogE_best;
            cout << " \t Merged ICCs " << i_keep << " & " << i_erase << " into " << i_keep;
            cout << " (delete " << i_erase << ")" << endl;   
        }
//        if(Greedy_Full_merging) 
//            {  PrintFile_MCM_FullMerge(Best_MCM, r, file_MCMs); }
    } // End While

    cout << "\t **** Greedy Merging Finished ****" << endl;
    cout << endl;
    cout.unsetf(ios_base::fixed);
    cout.precision(ss);

    return The_Best_MCM;
}


/******************************************************************************/
/***************    FIND THE BEST MCM: GREEDY PROCEDURE    ********************/
/*****************      ******    *******    ******      **********************/
/***************    Starting from the INDEPENDENT model    ********************/
/******************************************************************************/

map<unsigned int, __int128_t> MCM_GreedySearch(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, bool print_info = true) //, bool Greedy_Full_merging = false)
{
    cout << "The procedure starts from an independent model with r = " << r << " ICCs, each of size 1" << endl << endl;
    // ** Create initial independent model:
    map<unsigned int, __int128_t> MCM_0;
    __int128_t Op = 1;

    for (unsigned int i = 0; i < r; i++)
    {
        MCM_0[i] = Op;
        Op = Op << 1;
    }

    // ** Call Greedy Algo for MCM_0:
    return MCM_GreedySearch_MCM0(Kset, N, r, MCM_0, print_info); //, Greedy_Full_merging);
}

/******************************************************************************/
/***************    FIND THE BEST MCM: GREEDY PROCEDURE    ********************/
/***************          Starting from chosen MCM_0       ********************/
/******************************************************************************/
// r = number of communities in MCM_0

map<unsigned int, __int128_t> MCM_GreedySearch_MCM0_bis(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, map<unsigned int, __int128_t> MCM_0, bool print_info = true)
{
    // ** Create MCM_Best and initialise it to MCM_0:   
    map<unsigned int, __int128_t> Best_MCM(MCM_0);

    // ** Number of ICCs in MCM_Best:
    unsigned int r = Best_MCM.size();

    // ** Store values of LogE for each ICC:
    map<unsigned int, double> Best_MCM_LogE_ICC;
    double LogE_tot = 0;

/*
    for (unsigned int i = 0; i < r; i++)
    {
        Best_MCM_LogE_ICC[i] = LogE_ICC(Kset, Best_MCM[i], N);
        LogE_tot += Best_MCM_LogE_ICC[i];
    }
*/
    for (auto const& it : Best_MCM)
    {
        Best_MCM_LogE_ICC[it.first] = LogE_ICC(Kset, it.second, N);
        LogE_tot += Best_MCM_LogE_ICC[it.first];
    }

    // ** Create table of LogE values for all pairs of merged ICCs:
    cout << "\t **** Filling in the rxr-triangular Matrix ****" << endl;
    double* logE_merged = (double*)malloc(r * (r - 1) / 2 * sizeof(double));

    // logE_merged[mat_loc(i,j)] will contain the values of LogE of merged(ICCs i and j)
    for (unsigned int i = 0; i < r - 1; i++)
    {
        for (unsigned int j = i + 1; j < r; j++)
        {
            logE_merged[mat_loc(i, j)] = 0;  // initalized to 0
        }
    }
    cout << endl;

    // buffers:
    map<unsigned int, __int128_t> Buffer_ICC_unmerged;
    double Diff_LogE_best = 0, LogE_merged_best = 0, LogE_unmerged = 0;

    // Best:
    map<unsigned int, __int128_t>::iterator it1, it2;
    bool stop = false;

    unsigned int i_keep, i_erase, i_mat;   // define variables for merging

    // ********* HIERARCHICAL MERGING ****************
    cout << "\t **** Start Merging ****" << endl;

    streamsize ss = cout.precision();
    cout << fixed << setprecision(2); // precision of LogE: print 2 digits after "."
    cout << "   Start: \t\t Nb of communities = " << Best_MCM.size();
    cout << "\t LogE = " << LogE_tot << endl; // << "\t Check: " << LogE_MCM(Kset, Best_MCM, N, r) << endl;

    int iteration = 0; 
    while (Best_MCM.size() > 1 && (!stop))
    {
        Diff_LogE_best = 0;
        LogE_unmerged = 0; //LogE_merged = 0, 

        // ****** Find the best pair to merge: *********
        for (it1 = Best_MCM.begin(); it1 != (--Best_MCM.end()); ++it1)
        {
            //cout << ">> test for it1 = " << (*it1).first << endl;

            for (it2 = next(it1); it2 != Best_MCM.end(); ++it2)
            {
                //cout << "test for it2 = " << (*it2).first << "\t ";

                // LogE of unmerged ICCs it1 and it2:
                LogE_unmerged = Best_MCM_LogE_ICC[(*it1).first] + Best_MCM_LogE_ICC[(*it2).first];

                // LogE of merged ICC (it1 + it2):
                i_mat = mat_loc((*it1).first, (*it2).first);
                //cout << "i_mat = " << i_mat << "\t ";
                if (logE_merged[i_mat] == 0)
                {
                    //cout << "is 0: ";
                    logE_merged[i_mat] = LogE_ICC(Kset, (*it1).second + (*it2).second, N);
                }
                //cout << "LogE[i_mat] = " << logE_merged[i_mat] << "\t";
                //cout << "LogE_unmerged = " << LogE_unmerged  << endl;

                // if merged ICC is better, then save it:
                if ((logE_merged[i_mat] - LogE_unmerged) > Diff_LogE_best)
                {
                    //cout << "Is larger \t";
                    Buffer_ICC_unmerged.clear();
                    Buffer_ICC_unmerged.insert(*it1);
                    Buffer_ICC_unmerged.insert(*it2);

                    Diff_LogE_best = logE_merged[i_mat] - LogE_unmerged;
                    LogE_merged_best = logE_merged[i_mat];
                    //cout << "Diff LogE best = " << Diff_LogE_best << " : " << (*it1).first << "\t " << (*it2).first << endl;
                }
            }
        }

        // ********* STOPPING CRITERIA **************:
        if (Diff_LogE_best == 0)  // there is no improvement of Log-E
        { 
            stop = true; 
            if (print_info) {cout << "   Stop: No more improvement of LogE" << endl; } 
            continue; 
        }

        // ********* PERFORM MERGING *****:
        it1 = Buffer_ICC_unmerged.begin();  i_keep = (*it1).first;   //  store the merged ICC in lowest index: i_keep
        it1++;  i_erase = (*it1).first; //  erase the ICC in the largest index: i_erase

        Best_MCM[i_keep] = Best_MCM[i_keep] + Best_MCM[i_erase];
        Best_MCM.erase(i_erase);

        Best_MCM_LogE_ICC[i_keep] = LogE_merged_best;
        Best_MCM_LogE_ICC.erase(i_erase);

        LogE_tot += Diff_LogE_best;

        // ******* For all the indices "i" that are still stored in Best_MCM:
        // ************ reset the values of logE_merged between "i" and "i_keep" to "0"
        for (it1 = Best_MCM.begin(); it1 != Best_MCM.end(); it1++)
        {
            if ((*it1).first > i_keep)
            {
                logE_merged[mat_loc(i_keep, (*it1).first)] = 0;
            }
            else if ((*it1).first < i_keep)
            {
                logE_merged[mat_loc((*it1).first, i_keep)] = 0;
            }
        }

        iteration++;
        if (print_info)
        {   
            cout << "   Iteration " << iteration << " done: \t Nb of communities = " << Best_MCM.size();
            cout << "\t LogE = " << LogE_tot; //<< "\t Check: "<< LogE_MCM(Kset, Best_MCM, N, r) << endl;
            cout << "\t Merged ICCs " << i_keep << " & " << i_erase << " into " << i_keep;
            cout << " (delete " << i_erase << ")" << endl;   
        }
    } // End While

    cout << "\t **** Greedy Merging Finished ****" << endl;
    cout << endl;
    cout.unsetf(ios_base::fixed);
    cout.precision(ss);

    return Best_MCM;
}

/******************************************************************************/
/***************    FIND THE BEST MCM: GREEDY PROCEDURE    ********************/
/*****************      ******    *******    ******      **********************/
/***************    Starting from the INDEPENDENT model    ********************/
/******************************************************************************/
map<unsigned int, __int128_t> MCM_GreedySearch_bis(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, bool print_info = true)
{
    cout << "Define independent model: " << endl;
    // ** Create initial independent model:
    map<unsigned int, __int128_t> MCM_0;
    __int128_t Op = 1;

    for (unsigned int i = 0; i < r; i++)
    {
        MCM_0[i] = Op;
        Op = Op << 1;
    }

    // ** Call Greedy Algo for MCM_0:
    return MCM_GreedySearch_MCM0_bis(Kset, N, MCM_0, print_info);
}

/*
map<unsigned int, __int128_t> MCM_GreedySearch_bis(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, bool print_info = true)
{
    // Create initial independent model:
    map<unsigned int, __int128_t> Best_MCM;  __int128_t Op = 1;
    map<unsigned int, double> Best_MCM_LogE_ICC;

    cout << "\t **** Filling in the rxr-triangular Matrix ****" << endl;
    double* logE_merged = (double*)malloc(r * (r - 1) / 2 * sizeof(double));

    // logE_merged[mat_loc(i,j)] will contain the values of LogE of merged(ICCs i and j)
    for (unsigned int i = 0; i < r - 1; i++)
    {
        for (unsigned int j = i + 1; j < r; j++)
        {
            logE_merged[mat_loc(i, j)] = 0;  // initalized to 0"
        }
    }
    cout << endl;

    // Fill Best_MCM[i] with the choice of the initial MCM_0;
    double LogE_tot = 0;
    for (unsigned int i = 0; i < r; i++)
    {
        Best_MCM[i] = Op;
        Best_MCM_LogE_ICC[i] = LogE_ICC(Kset, Op, N);
        LogE_tot += Best_MCM_LogE_ICC[i];
        Op = Op << 1;
    }

    // buffers:
    map<unsigned int, __int128_t> Buffer_ICC_unmerged;
    double Diff_LogE_best = 0, LogE_merged_best = 0, LogE_unmerged = 0;

    // Best:
    map<unsigned int, __int128_t>::iterator it1, it2;
    bool stop = false;

    unsigned int i_keep, i_erase, i_mat;   // define variables for merging

    int iteration = 0; 
    
    cout << "\t **** Start Merging ****" << endl;

    // ********* HIERARCHICAL MERGING ****************
    streamsize ss = cout.precision();
    cout << fixed << setprecision(2); // precision of LogE: print 2 digits after "."
    cout << "   Start: \t\t Nb of communities = " << Best_MCM.size();
    cout << "\t LogE = " << LogE_tot << endl; // << "\t Check: " << LogE_MCM(Kset, Best_MCM, N, r) << endl;

    while (Best_MCM.size() > 1 && (!stop))
    {
        Diff_LogE_best = 0;
        LogE_unmerged = 0; //LogE_merged = 0, 

        // ****** Find the best pair to merge: *********
        for (it1 = Best_MCM.begin(); it1 != (--Best_MCM.end()); ++it1)
        {
            //cout << ">> test for it1 = " << (*it1).first << endl;

            for (it2 = next(it1); it2 != Best_MCM.end(); ++it2)
            {
                //cout << "test for it2 = " << (*it2).first << "\t ";

                // LogE of unmerged ICCs it1 and it2:
                LogE_unmerged = Best_MCM_LogE_ICC[(*it1).first] + Best_MCM_LogE_ICC[(*it2).first];

                // LogE of merged ICC (it1 + it2):
                i_mat = mat_loc((*it1).first, (*it2).first);
                //cout << "i_mat = " << i_mat << "\t ";
                if (logE_merged[i_mat] == 0)
                {
                    //cout << "is 0: ";
                    logE_merged[i_mat] = LogE_ICC(Kset, (*it1).second + (*it2).second, N);
                }
                //cout << "LogE[i_mat] = " << logE_merged[i_mat] << "\t";
                //cout << "LogE_unmerged = " << LogE_unmerged  << endl;

                // if merged ICC is better, then save it:
                if ((logE_merged[i_mat] - LogE_unmerged) > Diff_LogE_best)
                {
                    //cout << "Is larger \t";
                    Buffer_ICC_unmerged.clear();
                    Buffer_ICC_unmerged.insert(*it1);
                    Buffer_ICC_unmerged.insert(*it2);

                    Diff_LogE_best = logE_merged[i_mat] - LogE_unmerged;
                    LogE_merged_best = logE_merged[i_mat];
                    //cout << "Diff LogE best = " << Diff_LogE_best << " : " << (*it1).first << "\t " << (*it2).first << endl;
                }
            }
        }

        // ********* STOPPING CRITERIA **************:
        if (Diff_LogE_best == 0)  // there is no improvement of Log-E
        { 
            stop = true; 
            if (print_info) {cout << "   Stop: No more improvement of LogE" << endl; } 
            continue; 
        }

        // ********* PERFORM MERGING *****:
        it1 = Buffer_ICC_unmerged.begin();  i_keep = (*it1).first;   //  store the merged ICC in lowest index: i_keep
        it1++;  i_erase = (*it1).first; //  erase the ICC in the largest index: i_erase

        Best_MCM[i_keep] = Best_MCM[i_keep] + Best_MCM[i_erase];
        Best_MCM.erase(i_erase);

        Best_MCM_LogE_ICC[i_keep] = LogE_merged_best;
        Best_MCM_LogE_ICC.erase(i_erase);

        LogE_tot += Diff_LogE_best;

        // ******* For all the indices "i" that are still stored in Best_MCM:
        // ************ reset the values of logE_merged between "i" and "i_keep" to "0"
        for (it1 = Best_MCM.begin(); it1 != Best_MCM.end(); it1++)
        {
            if ((*it1).first > i_keep)
            {
                logE_merged[mat_loc(i_keep, (*it1).first)] = 0;
            }
            else if ((*it1).first < i_keep)
            {
                logE_merged[mat_loc((*it1).first, i_keep)] = 0;
            }
        }

        iteration++;
        if (print_info)
        {   
            cout << "   Iteration " << iteration << " done: \t Nb of communities = " << Best_MCM.size();
            cout << "\t LogE = " << LogE_tot; //<< "\t Check: "<< LogE_MCM(Kset, Best_MCM, N, r) << endl;
            cout << "\t Merged ICCs " << i_keep << " & " << i_erase << " into " << i_keep;
            cout << " (delete " << i_erase << ")" << endl;   
        }
    } // End While

    cout << "\t **** Greedy Merging Finished ****" << endl;
    cout << endl;
    cout.unsetf(ios_base::fixed);
    cout.precision(ss);

    return Best_MCM;
}
*/


