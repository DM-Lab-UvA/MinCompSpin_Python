#include <iostream>
#include <list>
#include <fstream>
#include <sstream>
#include <map>
#include <bitset>
#include<cmath>

using namespace std;

#include "data.h"

// Create a bitset for integers up to 2**128
bitset<n> bitset128(__int128_t Op)
{
    bitset<n> hi{ static_cast<unsigned long long>(Op >> 64) },
        lo{ static_cast<unsigned long long>(Op) },
        bits{ (hi << 64) | lo }; // Create binary representation

    return bits;
}

map<unsigned int, list<Interaction>> hypergraph_interactions(double J, unsigned int s, string file)
{
    string line, word;
    __int128_t Op2, Op = 1;
    Op <<= (n - 1);

    char c = '1';

    map<unsigned int, list<Interaction>> I_list;
    Interaction I;

    // Open file and read hypergraph
    ifstream myfile(file.c_str());
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            stringstream ss(line);
            string Operator, Weight;

            if (ss >> Operator >> Weight)
            {
                Op2 = 1;
                Op2 = (Op2 << n - 1);
                I.Op = 0;
                for (auto &elem: Operator)
                {
                    if (elem == c) { I.Op += Op2; }
                    Op2 = Op2 >> 1;
                }


                I.g = (double)stoi(Weight) * J;
            }

            for (unsigned int i = 0; i < n; i++)
            {
                if ((Op >> i) & I.Op)
                {
                    I_list[i].push_back(I);
                }
            }
        }
        myfile.close();
    }

    return I_list;
}

/*
list<Interaction> write_interactions(double J, string file)
{
    string line, line2;
    list<Interaction> I_list;
    Interaction I;

    __int128_t Op = 1;
    Op <<= (n - 1);

    bool flag = false;

    // Open file to read network
    ifstream myfile(file.c_str());
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            stringstream ss(line);
            I.Op = 0;
            I.g = J;
            while (getline(ss, line2, '\t'))
            {
                if ((Op >> stoi(line2) - 1) < I.Op) { flag = true; }
                I.Op += (Op >> (stoi(line2) - 1));
            }
            if (flag) { flag = false; continue; }
            I_list.push_back(I);
        }
        myfile.close();
    }
    return I_list;
}*/

map<unsigned int, list<Interaction>> write_interactions_metropolis(double J, string file)
{
    string line, line2;
    unsigned int node; // First node involved with an interaction
    map<unsigned int, list<Interaction>> I_list; // Iterator for list of interactions
    Interaction I;

    __int128_t Op = 1;
    Op <<= (n - 1);

    int flag;

    // Open file to read network
    ifstream myfile(file.c_str());
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            stringstream ss(line);
            // Create interaction
            I.Op = 0;
            I.g = J;
            // Read nodes involved in interaction
            flag = 2;            
            while (getline(ss, line2, '\t'))
            {
                if (flag > 0) 
                {
                    node = stoi(line2) - 1;
                    I.Op += (Op >> node);
                }
                else
                { 
                    I.g *= stoi(line2); 
                }
                flag--;
            }
            // Add interaction to list corresponding to 2nd node
            I_list[node].push_back(I);

            //bitset<n> bits = bitset128(I.Op);
            //cout << bits << endl;
        }
        myfile.close();
    }
    return I_list;
}

double delta_energy(__int128_t state, list<Interaction> edges, int spin)
{
    //signed int sign = (signed)bitset128((1 << spin) & state).count() * 2 - 1; // Get sign of spin;

    __int128_t Op;
    signed int sisj;
    double diff_E = 0;

    for (auto& edge : edges)
    {
        // Apply operator to state
        Op = state & edge.Op;
        // Calculate s_i * s_j
        bitset<n> bits = bitset128(Op);
        sisj = 2 * ( (1 + (signed)bits.count()) % 2 ) - 1;
        // Add to sum
        diff_E += 2.0 * edge.g * (double)sisj; // Difference in energy
    }
    return diff_E;
}

void sample_data_metropolis(double J, unsigned int s, string input_file, string output_filename, unsigned int N = 1000)
{
    map<unsigned int, list<Interaction>> IA;
    list<Interaction> edges;
    int spin;
    __int128_t Op, state = 0, one = 1;
    double eps, energy, diff_E;

    // INTERACTIONS:
    //IA = hypergraph_interactions(J, s, input_file);
    IA = write_interactions_metropolis(J, input_file);

    /*Interaction I;
    I.Op = 3;
    I.g = 4 * J;

    IA[21].push_back(I);
    IA[20].push_back(I);

    
    IA[21].push_back(I);
    IA[2].push_back(I);
    IA[1].push_back(I);
    IA[0].push_back(I);*/

    // INITIAL ENERGY
    for (auto& elem_map : IA)
    {
        for (auto& edge : elem_map.second)
        {
            // Apply operator to state
            Op = state & edge.Op;
            // Create binary representation
            bitset<n> bits = bitset128(Op);
            // Add to energy
            energy += edge.g * (2 * ((1 + (signed)bits.count()) % 2) - 1);
        }
    }

    // WINDUP-PERIOD
    for (int interval = 0; interval < 100 * n; interval++)
    {
        // Select a node
        spin = rand() % n;
        // Draw uniform number
        eps = (double)rand() / RAND_MAX;
        // Select interactions that involve the selected node
        edges = IA[spin];
        // Calculate the difference in energy
        diff_E = delta_energy(state, edges, spin);
        // Accept/reject new state
        if (diff_E < 0) { state ^= (one << (n - 1 - spin)); state = ~state; }
        else if (eps < exp(-diff_E)) { state ^= (one << (n - 1 - spin)); state = ~state; }
    }

    // OUTPUT FILE:
    fstream file(output_filename.c_str(), ios::out);

    float progress; //cout << endl;

    for (int i = 0; i < N; i++)
    {
        for (int interval = 0; interval < 100 * n; interval++)
        {
            // Select a node
            spin = rand() % n;
            // Draw uniform number
            eps = (double)rand() / RAND_MAX;
            // Select interactions that involve the selected node
            edges = IA[spin];
            // Calculate the difference in energy
            diff_E = delta_energy(state, edges, spin);
            // Accept/reject new state
            if (diff_E < 0) { state ^= (one << (n - 1 - spin)); state = ~state; }
            else if (eps < exp(-diff_E)) { state ^= (one << (n - 1 - spin)); state = ~state; }
        }
        bitset<n> bits = bitset128(state);
        
        /*int barWidth = 70;

        std::cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();

        progress += 1.0/(float)N; // for demonstration only*/


        file << bits << endl;
    }

    //cout << endl;

    file.close();
}
