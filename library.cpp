#include <cassert>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <stdio.h>
#include <vector>
#include <string>
#include <list>
#include <map>

#include "MinCompSpin_Greedy/Libraries/library.hpp"

std::map<unsigned int, __int128_t> MCM_GreedySearch_AND_printInfo(std::vector<std::pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, bool print_it = false);
std::map<unsigned int, __int128_t> MCM_GreedySearch(std::vector<std::pair<__int128_t, unsigned int>> Kset, unsigned int N, unsigned int r, bool print_it = false);

/* Every item in the CppHistogram maps to 3 uint64s in the PyHistogram
 * [high_bits, low_bits, count]
 */
typedef std::vector<std::pair<__int128_t, unsigned int>> CppKset;
typedef std::map<unsigned int, __int128_t> CppModel;

class MCM_ModelBuffer {
public:
    MCM_ModelBuffer(size_t partitions, unsigned int r) : m_partitions(partitions), m_r(r) {
        m_data = new uint64_t[partitions * (r + 1)];
    }
    uint64_t *data() { return m_data; }
    size_t partitions() const { return m_partitions; }
    size_t r() const { return m_r; }
private:
    size_t m_partitions, m_r;
    uint64_t *m_data;
};

static void declare_MCM_Model(pybind11::handle m) {
    pybind11::class_<MCM_ModelBuffer>(m, "MCM_Model", pybind11::buffer_protocol())
        .def_buffer([](MCM_ModelBuffer &m) -> pybind11::buffer_info {
            return pybind11::buffer_info(
                m.data(),                                        /* pointer to buffer */
                sizeof(float),                                   /* size of one scalar */
                pybind11::format_descriptor<uint64_t>::format(), /* Python format descriptor */
                2,                                               /* Number of dimensions */
                { m.partitions(), m.r() },                       /* Buffer dimensions */
                { sizeof(float) * (m.r() + 1),                   /* Strides in bytes for each dimensions */
                  sizeof(float) },
                true                                             /* Is the data readonly */
            );
        });
}

void test(pybind11::buffer buf) {
    pybind11::buffer_info info = buf.request();
    assert(info.ndim == 2);
    char *ptr = (char *)info.ptr;
    printf("1: %ld, 2: %ld\n", info.shape[0], info.shape[1]);
    for (ssize_t i = 0; i < info.shape[0]; i++) {
        for (ssize_t j = 0; j < info.shape[1]; j++) {
            uint64_t *val = (uint64_t *)(ptr + info.strides[0] * i + info.strides[1] * j);
            printf("%lu, ", *val);
            *val += 1;
        }
        printf("\n");
    }
}

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

static CppKset py_Kset_to_cpp(std::vector<uint64_t> py_Kset) {
    size_t nr_items = py_Kset.size() / 3;
    CppKset cpp_Kset = CppKset(nr_items);
    for (size_t i = 0; i < nr_items; i++) {
        std::pair<__int128_t, unsigned int> item;
        item.first = ((__int128_t)py_Kset[i * 3] << 64) | (__int128_t)py_Kset[i * 3 + 1];
        item.second = py_Kset[i * 3 + 2];
        cpp_Kset[i] = item;
    }
    return cpp_Kset;
}

static std::vector<uint64_t> cpp_Model_to_py(CppModel cpp_Model) {
    std::vector<uint64_t>py_Model = std::vector<uint64_t>(cpp_Model.size() * 3);

    int i = 0;
    for (auto const& item : cpp_Model) {
        auto id = item.first;
        auto partition = item.second;
        py_Model[i * 3 + 0] = id;
        py_Model[i * 3 + 1] = partition >> 64;
        py_Model[i * 3 + 2] = partition & ((__uint128_t)-1 >> 64);
        i++;
    }
    return py_Model;
}

static CppModel py_Model_to_cpp(std::vector<uint64_t> py_Model) {
    CppModel cpp_Model = CppModel();
    for (size_t i = 0; i < py_Model.size(); i += 3) {
        unsigned int id = py_Model[i + 0];        
        uint64_t high_bits = py_Model[i + 1];
        uint64_t low_bits = py_Model[i + 2];
        __uint128_t partition = ((__int128_t)high_bits << 64) | (__int128_t)low_bits;
        auto item = std::pair<unsigned int, __uint128_t>(id, partition);
        cpp_Model.insert(item);
    }
    return cpp_Model;
}

std::pair<std::vector<uint64_t>, unsigned int> read_datafile_wrapper(std::string file, unsigned int r) {
    unsigned int N;
    CppKset cpp_Kset = read_datafile(&N, file, r);
    std::vector<uint64_t> py_Kset = std::vector<uint64_t>(cpp_Kset.size() * 3); // 3 uint64_t per entry

    int i = 0;
    for (auto const& cpp_pair : cpp_Kset) {
        uint64_t high_bits = cpp_pair.first >> 64;
        uint64_t low_bits = cpp_pair.first & ((__uint128_t)-1 >> 64);
        uint64_t count = cpp_pair.second;
        py_Kset[i * 3 + 0] = high_bits;
        py_Kset[i * 3 + 1] = low_bits;
        py_Kset[i * 3 + 2] = count;
        i++;
    }

    return std::pair<std::vector<uint64_t>, unsigned int>(py_Kset, N);
}

std::vector<uint64_t> MCM_GreedySearch_AND_printInfo_wrapper(std::vector<uint64_t> Kset, unsigned int N, unsigned int r, bool print_it = false) {
    CppKset cpp_Kset = py_Kset_to_cpp(Kset);
    auto cpp_Model = MCM_GreedySearch_AND_printInfo(cpp_Kset, N, r, print_it);
    return cpp_Model_to_py(cpp_Model);
}

std::vector<uint64_t> MCM_GreedySearch_wrapper(std::vector<uint64_t> Kset, unsigned int N, unsigned int r, bool print_it = false) {
    CppKset cpp_Kset = py_Kset_to_cpp(Kset);
    CppModel cpp_Model = MCM_GreedySearch(cpp_Kset, N, r, print_it);
    return cpp_Model_to_py(cpp_Model);
}

double Entropy_wrapper(std::vector<uint64_t> Kset, unsigned int N) {
    CppKset cpp_Kset = py_Kset_to_cpp(Kset);
    return Entropy(cpp_Kset, N);
}

double LogE_MCM_wrapper(
        std::vector<uint64_t> Kset,
        std::vector<uint64_t> Partition,
        unsigned int N, unsigned int r) {
    CppKset cpp_Kset = py_Kset_to_cpp(Kset);
    CppModel cpp_Partition = py_Model_to_cpp(Partition);
    return LogE_MCM(cpp_Kset, cpp_Partition, N, r);
}

double LogL_MCM_wrapper(
        std::vector<uint64_t> Kset,
        std::vector<uint64_t> Partition,
        unsigned int N, unsigned int r) {
    CppKset cpp_Kset = py_Kset_to_cpp(Kset);
    CppModel cpp_Partition = py_Model_to_cpp(Partition);
    return LogL_MCM(cpp_Kset, cpp_Partition, N, r);
}

void Print_MCM_Partition_wrapper(
        std::vector<uint64_t> Partition,
        unsigned int r) {
    CppModel cpp_Partition = py_Model_to_cpp(Partition);
    return Print_MCM_Partition(cpp_Partition, r);
}

double LogE_MCM_infoICC_wrapper(
        std::vector<uint64_t> Kset,
        std::vector<uint64_t> Partition,
        unsigned int N, unsigned int r) {
    CppKset cpp_Kset = py_Kset_to_cpp(Kset);
    CppModel cpp_Partition = py_Model_to_cpp(Partition);
    return LogE_MCM_infoICC(cpp_Kset, cpp_Partition, N, r);
}

PYBIND11_MODULE(MinCompSpin, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    declare_MCM_Model(m);

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
    m.def("test", &test);
}

