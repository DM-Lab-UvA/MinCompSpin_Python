//#include <cassert>
#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <cstdint>
#include <stdio.h>
#include <vector>
#include <string>
#include <list>
#include <map>

//#include "MinCompSpin_Greedy/Libraries/library.hpp"
#include "../includes/default_datafile.h"

std::string OutputFile_Add_Location(std::string filename)
{
    return (OUTPUT_directory + filename);
}

std::map<unsigned int, __int128_t> MCM_GreedySearch_AND_printInfo(std::vector<std::pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, bool print_it = false);
std::map<unsigned int, __int128_t> MCM_GreedySearch(std::vector<std::pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, bool print_it = false);

typedef std::vector<std::pair<__int128_t, unsigned int>> GreedyKset;
typedef std::map<unsigned int, __int128_t> GreedyPartitions;

/* Frees memory for an array of uint64_t or uint32_t.
 * Used when generating handles for python. */
void array_u64_destructor(void *ptr) {
    uint64_t *buf = reinterpret_cast<uint64_t *>(ptr);
    std::cerr << "Element [0] = " << buf[0] << "\n";
    std::cerr << "freeing memory @ " << ptr << "\n";
    delete[] buf;
}

/*
 * m_array is a numpy array with dtype np.uint64. The shape is (N, 3).
 * Every row is an item in the histogram. The value at index [i, 0]
 * (where i is the index of the item) is the count. The value at [i, 1]
 * are the least significant 64 bits of the original 128 bit value.
 * The value at [i, 2] contains the most significant 64 bits.
 */
class MCM_Kset {
public:
    pybind11::array_t<uint64_t> m_array;
    unsigned int m_N; // Number of points in the data set.

    MCM_Kset(pybind11::array_t<uint64_t> array, unsigned int N)
        : m_array(array), m_N(N) {}
    MCM_Kset(GreedyKset kset, unsigned int N) {
        size_t nr_bins = kset.size();
        uint64_t *buf = new uint64_t[nr_bins * 3];
        size_t i = 0;
        for (auto const& pair : kset) {
            buf[i*3+0] = pair.first & ((__uint128_t)-1 >> 64);
            buf[i*3+1] = pair.first >> 64;
            buf[i*3+2] = pair.second;
            i++;
        }

        // This allows python to free the memory
        pybind11::capsule free_when_done(buf, array_u64_destructor);
        m_array = pybind11::array_t<uint64_t>(
            { nr_bins, (size_t)3 },
            { sizeof(uint64_t) * 3, sizeof(uint64_t) },
            buf,
            free_when_done
        );
        m_N = N;
    }
    GreedyKset to_greedy() {
        size_t nr_bins = m_array.shape(0);
        GreedyKset greedy_Kset = GreedyKset(nr_bins);
        for (size_t i = 0; i < nr_bins; i++) {
            std::pair<__int128_t, unsigned int> item;
            item.first = (__int128_t)m_array.at(i, 0) |
                        ((__int128_t)m_array.at(i, 1) << 64);
            item.second = m_array.at(i, 2);
            greedy_Kset[i] = item;
        }
        return greedy_Kset;
    }
};

static void declare_MCM_Kset(pybind11::handle m) {
    pybind11::class_<MCM_Kset>(m, "Kset")
        .def(pybind11::init<pybind11::array_t<uint64_t>, unsigned int>())
        .def_readwrite("N", &MCM_Kset::m_N)
        .def_readwrite("array", &MCM_Kset::m_array);
}

/*
 * m_array is a numpy array with dtype uint64 and shape (n, 3) where n is the
 * number of partitions.
 */
class MCM_Partitions {
public:
    unsigned int m_r;
    pybind11::array_t<uint64_t> m_array;

    MCM_Partitions(pybind11::array_t<uint64_t> array, unsigned int r)
        : m_r(r), m_array(array) {}
    MCM_Partitions(GreedyPartitions greedy_partitions, unsigned int r) : m_r(r) {
        size_t partition_count = greedy_partitions.size();
        uint64_t *buf = new uint64_t[partition_count * 3];
        size_t i = 0;
        for (auto pair : greedy_partitions) {
            buf[i++] = pair.first;
            buf[i++] = pair.second & ((__int128)-1 >> 64);
            buf[i++] = pair.second >> 64;
        }

        // This allows python to free the memory
        pybind11::capsule free_when_done(buf, array_u64_destructor);
        m_array = pybind11::array_t<uint64_t>(
            { partition_count, (size_t)3 },
            { sizeof(uint64_t) * 3, sizeof(uint64_t) },
            buf,
            free_when_done
        );
    }
    GreedyPartitions to_greedy() {
        GreedyPartitions greedy_partitions;
        ssize_t partition_count = m_array.shape(0);
        for (ssize_t i = 0; i < partition_count; i++) {
            auto pair = std::pair<unsigned int, __int128_t>(
                m_array.at(i, 0),
                 (__int128_t)m_array.at(i, 1) |
                ((__int128_t)m_array.at(i, 2) << 64)
            );
            greedy_partitions.insert(pair);
        }
        return greedy_partitions;
    }
};

static void declare_MCM_Partitions(pybind11::handle m) {
    pybind11::class_<MCM_Partitions>(m, "Partitions")
        .def(pybind11::init<pybind11::array_t<uint64_t>, unsigned int>())
        .def_readwrite("r", &MCM_Partitions::m_r)
        .def_readwrite("array", &MCM_Partitions::m_array);
}

/*
int main(int argc, const char *argv[]);

const int argc = 3;
int main_wrapper(char *path, char *n_string_buffer) {
    const char *argv[argc] = {
        "invalid",
        path,
        n_string_buffer,
    };
    main(argc, argv);
    return 0;
}
*/

MCM_Kset read_datafile_wrapper(
        std::string file,
        unsigned int r) {
    unsigned int N;
    GreedyKset greedy_Kset = read_datafile(&N, file, r);
    MCM_Kset Kset(greedy_Kset, N);
    return Kset;
}

MCM_Partitions MCM_GreedySearch_AND_printInfo_wrapper(
        MCM_Kset Kset,
        unsigned int r,
        bool print_it = false) {
    GreedyKset greedy_Kset = Kset.to_greedy();
    auto greedy_partitions = MCM_GreedySearch_AND_printInfo(
        greedy_Kset, Kset.m_N, r, print_it);
    return MCM_Partitions(greedy_partitions, r);
}

MCM_Partitions MCM_GreedySearch_wrapper(
        MCM_Kset Kset,
        unsigned int r,
        bool print_it = false) {
    GreedyKset greedy_Kset = Kset.to_greedy();
    GreedyPartitions greedy_partitions = MCM_GreedySearch(greedy_Kset, Kset.m_N, r, print_it);
    auto model = MCM_Partitions(greedy_partitions, r);
    return model;
}

double Entropy_wrapper(MCM_Kset Kset) {
    GreedyKset greedy_Kset = Kset.to_greedy();
    return Entropy(greedy_Kset, Kset.m_N);
}

double LogE_MCM_wrapper(
        MCM_Kset Kset,
        MCM_Partitions Partition,
        unsigned int r) {
    GreedyKset greedy_Kset = Kset.to_greedy();
    GreedyPartitions greedy_partitions = Partition.to_greedy();
    return LogE_MCM(greedy_Kset, greedy_partitions, Kset.m_N, r);
}

double LogL_MCM_wrapper(
        MCM_Kset Kset,
        MCM_Partitions Partition,
        unsigned int r) {
    GreedyKset greedy_Kset = Kset.to_greedy();
    GreedyPartitions greedy_partitions = Partition.to_greedy();
    return LogL_MCM(greedy_Kset, greedy_partitions, Kset.m_N, r);
}

void Print_MCM_Partition_wrapper(
        MCM_Partitions Partition,
        unsigned int r) {
    GreedyPartitions greedy_partitions = Partition.to_greedy();
    return Print_MCM_Partition(greedy_partitions, r);
}

double LogE_MCM_infoICC_wrapper(
        MCM_Kset Kset,
        MCM_Partitions Partition,
        unsigned int r) {
    GreedyKset greedy_Kset = Kset.to_greedy();
    GreedyPartitions greedy_partitions = Partition.to_greedy();
    return LogE_MCM_infoICC(greedy_Kset, greedy_partitions, Kset.m_N, r);
}

PYBIND11_MODULE(MinCompSpin, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    declare_MCM_Kset(m);
    declare_MCM_Partitions(m);

//    m.def("main", &main_wrapper, "the main function");
    
    m.def("read_datafile", &read_datafile_wrapper);
    m.def("MCM_GreedySearch_AND_printInfo", &MCM_GreedySearch_AND_printInfo_wrapper,
        pybind11::arg("Kset"),
        pybind11::arg("r"),
        pybind11::arg("print_it") = false
    );
    m.def("MCM_GreedySearch", &MCM_GreedySearch_wrapper,
        pybind11::arg("Kset"),
        pybind11::arg("r"),
        pybind11::arg("print_it") = false
    );
    m.def("Entropy", &Entropy_wrapper);
    m.def("LogE_MCM", &LogE_MCM_wrapper);
    m.def("LogL_MCM", &LogL_MCM_wrapper);
    m.def("Print_MCM_Partition", &Print_MCM_Partition_wrapper);
    m.def("LogE_MCM_infoICC", &LogE_MCM_infoICC_wrapper);
    // m.def("test", &test);
}

