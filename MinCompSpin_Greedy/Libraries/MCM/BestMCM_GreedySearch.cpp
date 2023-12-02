#include <iostream>
#include <map>
#include <vector>
//#include <string>

using namespace std;

/******************************************************************************/
/******************   TOOL Functions from "tools.cpp"   ***********************/
/******************************************************************************/
//string int_to_bstring(__int128_t bool_nb, unsigned int n);

/******************************************************************************/
/************************   Function from "LogE.cpp"   ************************/
/******************************************************************************/
double LogE_ICC(vector<pair<__int128_t, unsigned int>> Kset, __int128_t Ai, unsigned int N);

/******************************************************************************/
/***************    FIND THE BEST MCM: GREEDY PROCEDURE    ********************/
/******************************************************************************/

map<unsigned int, __int128_t> MCM_GreedySearch(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, bool print_it = false)
{
    // Create initial independent model:
    map<unsigned int, __int128_t> Best_MCM_Partition;  __int128_t Op = 1;
    map<unsigned int, double> Best_MCM_LogE_Part;

    cout << "\t **** Filling in the rxr-triangular Matrix ****" << endl;
    double* logE_Mat = (double*)malloc(r * (r - 1) / 2 * sizeof(double));

    for (int i = 0; i < r - 1; i++)
    {
        for (int j = i + 1; j < r; j++)
        {
            logE_Mat[j * (j - 1) / 2 + i] = 0;
        }
    }

    for (unsigned int i = 0; i < r; i++)
    {
        Best_MCM_Partition[i] = Op;
        Best_MCM_LogE_Part[i] = LogE_ICC(Kset, Op, N);
        Op = Op << 1;
    }

    // buffers:
    map<unsigned int, __int128_t> MCM_unmerged, Buffer_Community;
    double Diff_LogE_best = 0, LogE_merged = 0, LogE_unmerged = 0, LogE_merged_best = 0, LogE_tot = 0;

    // Best:
    map<unsigned int, __int128_t>::iterator it1, it2, it1_end, it2_start, it;
    int counter = 0;
    bool stop = false;

    int i_keep, i_erase, i_mat; __int128_t new_part;   // define variables for merging

    int k = 0, iteration = 0;
    
    cout << "\t **** Start Merging ****" << endl;
    // ********* HIERARCHICAL MERGING ****************
    while (Best_MCM_Partition.size() > 1 && (!stop))
    {
        counter = 0;
        Diff_LogE_best = 0;
        LogE_merged = 0, LogE_unmerged = 0;

        // ****** Find the best pair to merge: *********
        //cout << endl << "******** Test combinations: " << endl;
        it1_end = Best_MCM_Partition.end();   it1_end--;
        for (it1 = Best_MCM_Partition.begin(); it1 != it1_end; it1++)
        {
            it2_start = it1;    it2_start++;
            for (it2 = it2_start; it2 != Best_MCM_Partition.end(); it2++)
            {
                i_mat = (*it2).first * ((*it2).first - 1) / 2 + (*it1).first;

                //cout << counter << "\t" << int_to_bstring((*it1).second, r) << "\t " << int_to_bstring((*it2).second, r);

                MCM_unmerged.insert(*it1);
                MCM_unmerged.insert(*it2);
                LogE_unmerged = Best_MCM_LogE_Part[(*it1).first] + Best_MCM_LogE_Part[(*it2).first]; //LogE_MCM(Kset, MCM_unmerged, N);

                if (logE_Mat[i_mat] == 0)
                {
                    logE_Mat[i_mat] = LogE_ICC(Kset, (*it1).second + (*it2).second, N);
                }
                //cout << "\t DlogE = " << (LogE_merged-LogE_unmerged);

                if ((logE_Mat[i_mat] - LogE_unmerged) > Diff_LogE_best)
                {
                    Buffer_Community = MCM_unmerged;
                    Diff_LogE_best = logE_Mat[i_mat] - LogE_unmerged;
                    LogE_merged_best = logE_Mat[i_mat];
                }
                //cout << "\t " << Diff_LogE_best << endl;
                //counter++;
                MCM_unmerged.clear();
            }
        }

        // ********* STOPPING CRITERIA **************:
        if (Diff_LogE_best == 0) { stop = true; continue; }

        // ********* PERFORM MERGING *****:
        it = Buffer_Community.begin();  i_keep = (*it).first;   //  new_part = (*it).second;
        it++;  i_erase = (*it).first; //  new_part += (*it).second;

        for (it1 = Best_MCM_Partition.begin(); it1 != Best_MCM_Partition.end(); it1++)
        {
            if ((*it1).first == i_keep) { continue; }
            if ((*it1).first > i_keep)
            {
                logE_Mat[(*it1).first * ((*it1).first - 1) / 2 + i_keep] = 0;
            }
            else
            {
                logE_Mat[i_keep * (i_keep - 1) / 2 + (*it1).first] = 0;
            }

        }

        Best_MCM_Partition[i_keep] = Best_MCM_Partition[i_keep] + Best_MCM_Partition[i_erase];
        Best_MCM_Partition.erase(i_erase);

        Best_MCM_LogE_Part[i_keep] = LogE_merged_best;
        Best_MCM_LogE_Part.erase(i_erase);

        iteration++;
        if (print_it)
            {   cout << "\t Done with iteration " << iteration << endl;   }
    } // End While

    cout << "\t **** Greedy Merging Finished ****" << endl;
    cout << endl;

    return Best_MCM_Partition;
}


/******************************************************************************/
/***************    FIND THE BEST MCM: GREEDY PROCEDURE    ********************/
/***************          Starting from chosen MCM         ********************/
/******************************************************************************/
// r = number of communities in MCM_0
map<unsigned int, __int128_t> MCM_GreedySearch_MCM0(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, map<unsigned int, __int128_t> MCM_0, bool print_it = false)
{
    // Create initial independent model:
    __int128_t Op = 1;
    map<unsigned int, __int128_t> Best_MCM_Partition;
    map<unsigned int, double> Best_MCM_LogE_Part;

    cout << "**** Filling in the rxr-triangular Matrix ****" << endl;
    double* logE_Mat = (double*)malloc(r * (r - 1) / 2 * sizeof(double));

    for (int i = 0; i < r - 1; i++)
    {
        for (int j = i + 1; j < r; j++)
        {
            logE_Mat[j * (j - 1) / 2 + i] = 0;
        }
    }

    for (unsigned int i = 0; i < r; i++)
    {
        Best_MCM_Partition[i] = MCM_0[i];
        Best_MCM_LogE_Part[i] = LogE_ICC(Kset, MCM_0[i], N);
    }

    // buffers:
    map<unsigned int, __int128_t> MCM_unmerged, Buffer_Community;
    double Diff_LogE_best = 0, LogE_merged = 0, LogE_unmerged = 0, LogE_merged_best = 0, LogE_tot = 0;

    // Best:
    map<unsigned int, __int128_t>::iterator it1, it2, it1_end, it2_start, it;
    int counter = 0;
    bool stop = false;

    int i_keep, i_erase, i_mat; __int128_t new_part;   // define variables for merging

    int k = 0, iteration = 0;
    
    cout << "**** Start Merging ****" << endl;
    // ********* HIERARCHICAL MERGING ****************
    while (Best_MCM_Partition.size() > 1 && (!stop))
    {
        counter = 0;
        Diff_LogE_best = 0;
        LogE_merged = 0, LogE_unmerged = 0;

        // ****** Find the best pair to merge: *********
        //cout << endl << "******** Test combinations: " << endl;
        it1_end = Best_MCM_Partition.end();   it1_end--;
        for (it1 = Best_MCM_Partition.begin(); it1 != it1_end; it1++)
        {
            it2_start = it1;    it2_start++;
            for (it2 = it2_start; it2 != Best_MCM_Partition.end(); it2++)
            {
                i_mat = (*it2).first * ((*it2).first - 1) / 2 + (*it1).first;

                //cout << counter << "\t" << int_to_bstring((*it1).second, r) << "\t " << int_to_bstring((*it2).second, r);

                MCM_unmerged.insert(*it1);
                MCM_unmerged.insert(*it2);
                LogE_unmerged = Best_MCM_LogE_Part[(*it1).first] + Best_MCM_LogE_Part[(*it2).first]; //LogE_MCM(Kset, MCM_unmerged, N);

                if (logE_Mat[i_mat] == 0)
                {
                    logE_Mat[i_mat] = LogE_ICC(Kset, (*it1).second + (*it2).second, N);
                }
                //cout << "\t DlogE = " << (LogE_merged-LogE_unmerged);

                if ((logE_Mat[i_mat] - LogE_unmerged) > Diff_LogE_best)
                {
                    Buffer_Community = MCM_unmerged;
                    Diff_LogE_best = logE_Mat[i_mat] - LogE_unmerged;
                    LogE_merged_best = logE_Mat[i_mat];
                }
                //cout << "\t " << Diff_LogE_best << endl;
                //counter++;
                MCM_unmerged.clear();
            }
        }

        // ********* STOPPING CRITERIA **************:
        if (Diff_LogE_best == 0) { stop = true; continue; }

        // ********* PERFORM MERGING *****:
        it = Buffer_Community.begin();  i_keep = (*it).first;   //  new_part = (*it).second;
        it++;  i_erase = (*it).first; //  new_part += (*it).second;

        for (it1 = Best_MCM_Partition.begin(); it1 != Best_MCM_Partition.end(); it1++)
        {
            if ((*it1).first == i_keep) { continue; }
            if ((*it1).first > i_keep)
            {
                logE_Mat[(*it1).first * ((*it1).first - 1) / 2 + i_keep] = 0;
            }
            else
            {
                logE_Mat[i_keep * (i_keep - 1) / 2 + (*it1).first] = 0;
            }

        }

        Best_MCM_Partition[i_keep] = Best_MCM_Partition[i_keep] + Best_MCM_Partition[i_erase];
        Best_MCM_Partition.erase(i_erase);

        Best_MCM_LogE_Part[i_keep] = LogE_merged_best;
        Best_MCM_LogE_Part.erase(i_erase);

        iteration++;
        if (print_it)
            {   cout << "\t Done with iteration " << iteration << endl;   }
    } // End While
    cout << "**** Greedy Merging Finished ****" << endl;
    cout << endl;

    return Best_MCM_Partition;
}
