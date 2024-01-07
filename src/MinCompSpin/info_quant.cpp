#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <numeric>  // for accumulate()

using namespace std;

/******************************************************************************/
/******************   TOOL Functions from "tools.cpp"   ***********************/
/******************************************************************************/
unsigned int Bitset_count(__int128_t bool_nb);

/******************************************************************************/
/**************************    Check sub-MCM:   *******************************/
/***************   Check if "fp1" is a sub-partition of "fp2":  ***************/
/******************************************************************************/

bool is_subset(map<unsigned int, __int128_t> fp1, map<unsigned int, __int128_t> fp2)
{
    bool flag;
    // Check if fp1 can be merged such that it becomes fp2
    for (auto& c1 : fp1)
    {
        flag = false;
        for (auto& c2 : fp2)
        {
            if ((c1.second & c2.second) == c1.second)
            {
                flag = true;
                break;
            }
        }
        if (!flag) { break; }
    }
    return flag;
}

/******************************************************************************/
/************************    COMPARE TWO PARTITIONS  **************************/
/******************************************************************************/

double Norm_Mut_info(map<unsigned int, __int128_t> Partition1, map<unsigned int, __int128_t> Partition2, unsigned int r)
{
    double I, H, p1, p2, p12;
    I = 0;  H = 0;
    int flag = 0;

    map<unsigned int, __int128_t>::iterator com1, com2;
    for (com1 = Partition1.begin(); com1 != Partition1.end(); com1++)
    {
        /*
        bitset<n> hi1{ static_cast<unsigned long long>((*com1).second >> 64) },
            lo1{ static_cast<unsigned long long>((*com1).second) },
            bits1{ (hi1 << 64) | lo1 };
        p1 = (double)bits1.count() / (double)(n);
        */
        p1 = (double)(Bitset_count((*com1).second)) / (double)(r);
        for (com2 = Partition2.begin(); com2 != Partition2.end(); com2++)
        {
            /*
            bitset<n> hi2{ static_cast<unsigned long long>((*com2).second >> 64) },
                lo2{ static_cast<unsigned long long>((*com2).second) },
                bits2{ (hi2 << 64) | lo2 };

            bitset<n> hi12{ static_cast<unsigned long long>(((*com1).second & (*com2).second) >> 64) },
                lo12{ static_cast<unsigned long long>((*com1).second & (*com2).second) },
                bits12{ (hi12 << 64) | lo12 };

            p2 = (double)bits2.count() / (double)(n);
            p12 = (double)bits12.count() / (double)(n);
            */
            p2 = (double)(Bitset_count((*com2).second)) / (double)(r);
            p12 = (double)(Bitset_count( (*com1).second & (*com2).second )) / (double)(r);            

            if (p12 != 0)
            {
                I += p12 * log(p12 / (p1 * p2));
            }
            if (flag < Partition2.size()) { H += p2 * log(p2); flag++; }
        }
        H += p1 * log(p1);
    }
    if (H == 0) { return 1; }
    else { return -2 * I / H; }
}

double Var_of_Inf(map<unsigned int, __int128_t> Partition1, map<unsigned int, __int128_t> Partition2, unsigned int r)
{
    // Variation of information calculates the distance between two partitions. The regular variation of information
    // is equal to the joint entropy minus the mutual information. However, the normalized version (divide by the joint entropy)
    // is preferred over the regular as this is a true metric, i.e., it satisfies the triangle inequality.
    double I, H, p1, p2, p12;
    I = 0;
    H = 0;
    map<unsigned int, __int128_t>::iterator com1, com2;
    for (com1 = Partition1.begin(); com1 != Partition1.end(); com1++)
    {
        /*bitset<n> hi1{ static_cast<unsigned long long>((*com1).second >> 64) },
            lo1{ static_cast<unsigned long long>((*com1).second) },
            bits1{ (hi1 << 64) | lo1 };
        p1 = (double)(bits1.count()) / (double)(n);
        */
        p1 = (double)(Bitset_count((*com1).second)) / (double)(r);

        for (com2 = Partition2.begin(); com2 != Partition2.end(); com2++)
        {
            /*
            bitset<n> hi2{ static_cast<unsigned long long>((*com2).second >> 64) },
                lo2{ static_cast<unsigned long long>((*com2).second) },
                bits2{ (hi2 << 64) | lo2 };

            bitset<n> hi12{ static_cast<unsigned long long>(((*com1).second & (*com2).second) >> 64) },
                lo12{ static_cast<unsigned long long>((*com1).second & (*com2).second) },
                bits12{ (hi12 << 64) | lo12 };

            p2 = (double)(bits2.count()) / (double)(n);
            p12 = (double)(bits12.count()) / (double)(n);
            */

            p2 = (double)(Bitset_count((*com2).second)) / (double)(r);
            p12 = (double)(Bitset_count( (*com1).second & (*com2).second )) / (double)(r);   

            if (p12 != 0)
            {
                I += p12 * log(p12 / (p1 * p2));
                H += p12 * log(p12);
            }
        }
    }
    if (H == 0) { return 0; }
    else { return 1 + I/H; }
}


/******************************************************************************/
/*******************************    ENTROPY  **********************************/
/******************************************************************************/

double Entropy(vector<pair<__int128_t, unsigned int>> Kset, unsigned int N)
{
    // Entropy of emperical data
    double H = 0;
    double p;

    for (auto& elem : Kset)
    {
        p = (double)elem.second / (double)N;
        H -= p * log2(p);
    }

    return H;
}


/******************************************************************************/
/************************    STATISTICAL DISTANCES  ***************************/
/******************************************************************************/
/*
map<__int128_t, double> cartesianProd(map<__int128_t, double> Map1, map<__int128_t, double> Map2, unsigned int N)
{
    map<__int128_t, double> Prod;

    for (auto& elem1 : Map1)
    {
        for (auto& elem2 : Map2)
        {
            Prod[elem1.first + elem2.first] = elem1.second * elem2.second;
        }
    }

    return Prod;
}

map<__int128_t, double> emp_dist(map<__int128_t, unsigned int> Kset, unsigned int N, unsigned int r)
{
    map<__int128_t, double> EMP;
    if (r > 20)
    {
        return EMP;
    }
    for (auto& elem : Kset)
    {
        EMP[elem.first] = (double)elem.second / (double)N;
    }
    return EMP;
}

map<__int128_t, double> MCM_distr(map<__int128_t, unsigned int> Kset, map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r)
{
    map<__int128_t, map<__int128_t, double>> sub_distr;
    map<__int128_t, double> distr;
    map<__int128_t, double> Map;

    map<unsigned int, __int128_t>::iterator it;

    __int128_t Ai, s;
    unsigned int ks, tot_size = 1;

    if (r > 20) { cout << "Too many nodes! JSD is unreliable!" << endl << endl; return distr; }

    for (auto& comm : Partition)
    {
        Ai = comm.second;
        for (auto& elem : Kset)
        {
            s = elem.first;
            ks = elem.second;
            sub_distr[Ai][Ai & s] += (double)ks;
        }
        tot_size *= sub_distr[Ai].size();
    }

    if (tot_size > 100000000) { cout << "Too many states! JSD is unreliable!" << endl << endl; return distr; }

    it = Partition.begin();
    distr = sub_distr[(*it).second]; it++;
    while (it != Partition.end())
    {
        Ai = (*it).second;
        Map = sub_distr[Ai];
        distr = cartesianProd(distr, Map, N);
        it++;
    }

    for (auto& elem : distr) { elem.second /= pow((double)N, (double)Partition.size()); }

    double sum = accumulate(distr.begin(), distr.end(), 0.0, [](double value, const std::map<__int128_t, double>::value_type& p)
        { return value + p.second; }
    );

    return distr;
}

double KL_divergence(map<__int128_t, unsigned int> Kset, map<unsigned int, __int128_t> Partition, unsigned int N)
{
    map<__int128_t, double> distr;
    map<__int128_t, unsigned int> Kset_new;

    map<__int128_t, unsigned int>::iterator it1, it2;
    map<unsigned int, __int128_t>::iterator it;

    __int128_t s, sig, Ai;
    unsigned int ks, kba;

    for (it1 = Kset.begin(); it1 != Kset.end(); it1++)
    {
        s = it1->first;
        distr[s] = 1.0;
    }

    for (it = Partition.begin(); it != Partition.end(); it++)
    {
        Kset_new.clear();
        Ai = (*it).second;
        for (it1 = Kset.begin(); it1 != Kset.end(); it1++)
        {
            s = it1->first;
            ks = it1->second;

            sig = s & Ai;

            Kset_new[sig] += ks;
        }


        for (it2 = Kset.begin(); it2 != Kset.end(); it2++)
        {
            s = it2->first;
            sig = s & Ai;

            distr[s] *= (double)Kset_new[sig];
        }
    }

    double P, Q, Dpq = 0;
    __int128_t norm = 0;

    for (it1 = Kset.begin(); it1 != Kset.end(); it1++)
    {
        s = it1->first;
        ks = it1->second;
        kba = distr[s];

        norm = pow((double)N, (double)Partition.size());

        P = (double)ks / (double)N;
        Q = (double)kba / norm;

        if (Q == 0) 
        { 
            cout << "No absolute continuity!" << endl; return -1; 
        }

        Dpq += P * log2(P / Q);
    }
    return Dpq;
}

double JS_divergence(map<__int128_t, double> Prob1, map<__int128_t, double> Prob2, unsigned int N)
{
    map<__int128_t, double> mcm_dist, av_dist;
    
    __int128_t s;
    unsigned int ks;

    for (auto& elem : Prob1)
    {
        av_dist[elem.first] += elem.second / 2.0;
    }

    for (auto& elem : Prob2)
    {
        av_dist[elem.first] += elem.second / 2.0;
    }

    double jsd = 0;
    for (auto& elem : Prob1)
    {
        jsd += elem.second * log2(elem.second / av_dist[elem.first]);
    }
    for (auto& elem : Prob2)
    {
        jsd += elem.second * log2(elem.second / av_dist[elem.first]);
    }

    return jsd / 2.0;
}
*/

