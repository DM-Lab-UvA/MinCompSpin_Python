#include <cassert>
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

#include "MinCompSpin_Greedy/Libraries/library.hpp"

std::map<unsigned int, __int128_t> MCM_GreedySearch_AND_printInfo(std::vector<std::pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, bool print_it = false);
std::map<unsigned int, __int128_t> MCM_GreedySearch(std::vector<std::pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, bool print_it = false);

typedef std::vector<std::pair<__int128_t, unsigned int>> GreedyKset;
typedef std::map<unsigned int, __int128_t> GreedyPartitions;
typedef std::pair<pybind11::array_t<uint64_t>, pybind11::array_t<uint32_t>> NumpyKset;
typedef std::pair<pybind11::array_t<uint32_t>, pybind11::array_t<uint64_t>> NumpyPartitions;

/*
 * m_array is a numpy array with dtype np.uint64. The shape is (N, 3).
 * Every row is an item in the histogram. The value at index [i, 0]
 * (where i is the index of the item) is the count. The value at [i, 1]
 * are the least significant 64 bits of the original 128 bit value.
 * The value at [i, 2] contains the most significant 64 bits.
 */
class MCM_Kset {
public:
    pybind11::array_t<uint64_t> m_bitsets;
    pybind11::array_t<uint32_t> m_counts;

    MCM_Kset(NumpyKset Kset) {
        m_bitsets = Kset.first;
        m_counts = Kset.second;
    }
    MCM_Kset(GreedyKset kset) {
        size_t N = kset.size();
        uint64_t *bitsets_buf = new uint64_t[N * 2];
        uint32_t *counts_buf = new uint32_t[N];
        size_t i = 0;
        for (auto const& pair : kset) {
            bitsets_buf[i*2+0] = pair.first & ((__uint128_t)-1 >> 64);
            bitsets_buf[i*2+1] = pair.first >> 64;
            counts_buf[i] = pair.second;
            i++;
        }

        // This allows python to free the memory
        pybind11::capsule bitsets_free_when_done(bitsets_buf, [](void *ptr) {
            uint64_t *buf = reinterpret_cast<uint64_t *>(ptr);
            std::cerr << "Element [0] = " << buf[0] << "\n";
            std::cerr << "freeing memory @ " << ptr << "\n";
            delete[] buf;
        });
        pybind11::capsule counts_free_when_done(counts_buf, [](void *ptr) {
            uint64_t *buf = reinterpret_cast<uint64_t *>(ptr);
            std::cerr << "Element [0] = " << buf[0] << "\n";
            std::cerr << "freeing memory @ " << ptr << "\n";
            delete[] buf;
        });
        m_bitsets = pybind11::array_t<uint64_t>(
            { N, (size_t)2 },
            { sizeof(uint64_t) * 2, sizeof(uint64_t) },
            bitsets_buf,
            bitsets_free_when_done
        );
        m_counts = pybind11::array_t<uint32_t>(
            { N, },
            { sizeof(uint32_t) },
            counts_buf,
            counts_free_when_done
        );
    }
    NumpyKset to_numpy() {
        return NumpyKset(m_bitsets, m_counts);
    }
    GreedyKset to_greedy() {
        size_t N = m_bitsets.shape(0);
        GreedyKset greedy_Kset = GreedyKset(N);
        for (size_t i = 0; i < N; i++) {
            std::pair<__int128_t, unsigned int> item;
            item.first = ((__int128_t)m_bitsets.at(i, 1) << 64) | (__int128_t)m_bitsets.at(i, 0);
            item.second = m_counts.at(i);
            greedy_Kset[i] = item;
        }
        return greedy_Kset;
    }
};

/*
 * m_array is a numpy array with dtype uint64 and shape (n, 3) where n is the
 * number of partitions.
 */
class MCM_Partitions {
public:
    size_t m_r;
    pybind11::array_t<uint64_t> m_ids;
    pybind11::array_t<uint64_t> m_partitions;

    MCM_Partitions(NumpyPartitions Partitions) {
        m_ids = Partitions.first;
        m_partitions = Partitions.second;
    }
    MCM_Partitions(GreedyPartitions greedy_partitions, unsigned int r) : m_r(r) {
        size_t partition_count = greedy_partitions.size();
        uint64_t *partitions_buf = new uint64_t[partition_count * 2];
        uint32_t *ids_buf = new uint32_t[partition_count];
        size_t i = 0;
        for (auto pair : greedy_partitions) {
            ids_buf[i] = pair.first;
            partitions_buf[i*2+0] = pair.second & ((__int128)-1 >> 64);
            partitions_buf[i*2+1] = pair.second >> 64;
            i++;
        }

        // This allows python to free the memory
        pybind11::capsule partitions_free_when_done(partitions_buf, [](void *ptr) {
            uint64_t *buf = reinterpret_cast<uint64_t *>(ptr);
            std::cerr << "Element [0] = " << buf[0] << "\n";
            std::cerr << "freeing memory @ " << ptr << "\n";
            delete[] buf;
        });
        pybind11::capsule ids_free_when_done(ids_buf, [](void *ptr) {
            uint32_t *buf = reinterpret_cast<uint32_t *>(ptr);
            std::cerr << "Element [0] = " << buf[0] << "\n";
            std::cerr << "freeing memory @ " << ptr << "\n";
            delete[] buf;
        });
        m_partitions = pybind11::array_t<uint64_t>(
            { partition_count, (size_t)2 },
            { sizeof(uint64_t) * 2, sizeof(uint64_t) },
            partitions_buf,
            partitions_free_when_done
        );
        m_ids = pybind11::array_t<uint32_t>(
            { partition_count },
            { sizeof(uint32_t) },
            ids_buf,
            ids_free_when_done
        );
    }
    NumpyPartitions to_numpy() {
        return NumpyPartitions(m_ids, m_partitions);
    }
    GreedyPartitions to_greedy() {
        GreedyPartitions greedy_partitions;
        ssize_t partition_count = m_ids.shape(0);
        for (ssize_t i = 0; i < partition_count; i++) {
            auto pair = std::pair<unsigned int, __int128_t>(
                m_ids.at(i),
                (__int128_t)m_partitions.at(i, 0) + ((__int128_t)m_partitions.at(i, 1) << 64)
            );
            greedy_partitions.insert(pair);
        }
        return greedy_partitions;
    }
};

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

std::pair<NumpyKset, unsigned int> read_datafile_wrapper(std::string file, unsigned int r) {
    unsigned int N;
    GreedyKset greedy_Kset = read_datafile(&N, file, r);
    MCM_Kset Kset(greedy_Kset);
    return std::pair<NumpyKset, unsigned int>(Kset.to_numpy(), N);
}

NumpyPartitions MCM_GreedySearch_AND_printInfo_wrapper(NumpyKset Kset, unsigned int N, unsigned int r, bool print_it = false) {
    GreedyKset greedy_Kset = MCM_Kset(Kset).to_greedy();
    auto greedy_partitions = MCM_GreedySearch_AND_printInfo(greedy_Kset, N, r, print_it);
    return MCM_Partitions(greedy_partitions, r).to_numpy();
}

NumpyPartitions MCM_GreedySearch_wrapper(NumpyKset Kset, unsigned int N, unsigned int r, bool print_it = false) {
    GreedyKset greedy_Kset = MCM_Kset(Kset).to_greedy();
    GreedyPartitions greedy_partitions = MCM_GreedySearch(greedy_Kset, N, r, print_it);
    auto model = MCM_Partitions(greedy_partitions, r);
    return model.to_numpy();
}

double Entropy_wrapper(NumpyKset Kset, unsigned int N) {
    GreedyKset greedy_Kset = MCM_Kset(Kset).to_greedy();
    return Entropy(greedy_Kset, N);
}

double LogE_MCM_wrapper(
        NumpyKset Kset,
        NumpyPartitions Partition,
        unsigned int N, unsigned int r) {
    GreedyKset greedy_Kset = MCM_Kset(Kset).to_greedy();
    GreedyPartitions greedy_partitions = MCM_Partitions(Partition).to_greedy();
    return LogE_MCM(greedy_Kset, greedy_partitions, N, r);
}

double LogL_MCM_wrapper(
        NumpyKset Kset,
        NumpyPartitions Partition,
        unsigned int N, unsigned int r) {
    GreedyKset greedy_Kset = MCM_Kset(Kset).to_greedy();
    GreedyPartitions greedy_partitions = MCM_Partitions(Partition).to_greedy();
    return LogL_MCM(greedy_Kset, greedy_partitions, N, r);
}

void Print_MCM_Partition_wrapper(
        NumpyPartitions Partition,
        unsigned int r) {
    GreedyPartitions greedy_partitions = MCM_Partitions(Partition).to_greedy();
    return Print_MCM_Partition(greedy_partitions, r);
}

double LogE_MCM_infoICC_wrapper(
        NumpyKset Kset,
        NumpyPartitions Partition,
        unsigned int N, unsigned int r) {
    GreedyKset greedy_Kset = MCM_Kset(Kset).to_greedy();
    GreedyPartitions greedy_partitions = MCM_Partitions(Partition).to_greedy();
    return LogE_MCM_infoICC(greedy_Kset, greedy_partitions, N, r);
}

PYBIND11_MODULE(MinCompSpin, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("main", &main_wrapper, "the main function");
    
    m.def("read_datafile", &read_datafile_wrapper);
    m.def("MCM_GreedySearch_AND_printInfo", &MCM_GreedySearch_AND_printInfo_wrapper,
        pybind11::arg("Kset"),
        pybind11::arg("N"),
        pybind11::arg("r"),
        pybind11::arg("print_it") = false
    );
    m.def("MCM_GreedySearch", &MCM_GreedySearch_wrapper,
        pybind11::arg("Kset"),
        pybind11::arg("N"),
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

