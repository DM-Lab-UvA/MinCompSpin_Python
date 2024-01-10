#include <iostream>
#include <vector>
#include <map>
#include <chrono>

/* ========================= GREEDY MCM DECLARATIONS ======================= */
std::vector<std::pair<__int128_t, unsigned int>> read_datafile(unsigned int *N, std::string file, unsigned int r); // filename to specify in data.h
std::map<unsigned int, __int128_t> MCM_GreedySearch(std::vector<std::pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, bool print_it = false);
double Entropy(std::vector<std::pair<__int128_t, unsigned int>> Kset, unsigned int N);
double LogL_MCM(std::vector<std::pair<__int128_t, unsigned int>> Kset, std::map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r);
double LogE_MCM(std::vector<std::pair<__int128_t, unsigned int>> Kset, std::map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r);
double LogL_MCM_infoICC(std::vector<std::pair<__int128_t, unsigned int>> Kset, std::map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r);
double LogE_MCM_infoICC(std::vector<std::pair<__int128_t, unsigned int>> Kset, std::map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r);

static struct UnusedVariables {
    double S, LogE_g, LogL_g, LogL_infoICC, LogE_infoICC;
} unused_variables;

void greedy_main(std::string path, unsigned int r) {
    unsigned int N;
    auto Kset = read_datafile(&N, path, r);
    auto fp1 = MCM_GreedySearch(Kset, N, r);
    unused_variables.S = Entropy(Kset, N);
    unused_variables.LogL_g = LogL_MCM(Kset, fp1, N, r);
    unused_variables.LogE_g = LogE_MCM(Kset, fp1, N, r);
    unused_variables.LogL_infoICC =  LogL_MCM_infoICC(Kset, fp1, N, r);
    unused_variables.LogE_infoICC =  LogE_MCM_infoICC(Kset, fp1, N, r);
double LogE_MCM_infoICC(std::vector<std::pair<__int128_t, unsigned int>> Kset, std::map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r);
}

// #define DEBUG

#define REPEAT 5
#ifdef DEBUG
#define NUMBER 200
#else
#define NUMBER 1000
#endif // DEBUG

int main(void) {
    std::vector<long int> results;
    for (size_t i = 0; i < REPEAT; i++) {
        auto t_start = std::chrono::high_resolution_clock::now();
        for (size_t j = 0; j < NUMBER; j++) {
            greedy_main("MinCompSpin_Greedy/INPUT/SCOTUS_n9_N895_Data.dat", 9);
        }
        auto t_end = std::chrono::high_resolution_clock::now();
        auto t_duration = t_end - t_start;
        long micro_duration =
            std::chrono::duration_cast<std::chrono::microseconds>(t_duration)
            .count();
        results.push_back(micro_duration);
    }

    for (size_t i = 0; i < REPEAT; i++) {
        double micro = results[i] / NUMBER;
        std::cout << micro << ", ";
    }
    std::cout << std::endl;
}
