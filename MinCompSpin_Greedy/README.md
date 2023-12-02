# Greedy search
# for the best Minimally Complex Models (MCM)

This program allows to **detect community structures in binary data**, while taking into account possible **high-order patterns of the data** in the detection of the communities (i.e., possible high-order correlations between the variables). The idea of the algorithm is based on performing statistical inference with a family of spin models (maximum entropy models for binary data) that have minimal information theoretic complexity. These models are called **Minimally Complex Models (MCM)**. We use the fact that each partition of the variables (i.e., a community structure) can be mapped to an MCM.

One huge advantage of Minimally Complex Models is that the model evidence (integration of the likelihood over the model parameter) has a known analytic expression that is easy to compute numerically. As a consequence, the computation of the model evidence doesn’t require fitting the parameters of the model nor computing a numerical integral, which allows to significantly accelerate the comparison between models. The selected model can also be used as a generative model for data. Details can be found in the paper [*Statistical Inference of Minimally Complex Models*](https://arxiv.org/abs/2008.00520) [1].

This repository contains a code initially developed for the paper Ref. [1] and later optimised for the paper Ref.[2]. **The program performs a greedy search** for the best Minimally Complex Spin Model (MCM) on a basis provided by the user. The comparison between models is based on their evidence (i.e., posterior probability that the model generates the data, integrated over its parameter values). The code performs an hierarchical merging procedure to find an optimal MCM in the basis chosen by the user. The code can run for systems with up to `n=128` variables. 

This program is complementary to the program available in [github.com/clelidm/MinCompSpin](https://github.com/clelidm/MinCompSpin), which performs an **exhaustive search** for the best community. We recommand using the greedy search over the exhaustive search when the number of variables exceeds `15`. A **simulated annealing version** of the optimization procedure can also be found in [github.com/ebokai/MinCompSpin_SimulatedAnnealing](https://github.com/ebokai/MinCompSpin_SimulatedAnnealing), which can find solutions closer to the global optimal.


[1]  C. de Mulatier, P. P. Mazza, M. Marsili, *Statistical Inference of Minimally Complex Models*, [arXiv:2008.00520](https://arxiv.org/abs/2008.00520)

[2] E. Peerbooms*, S. Kamphof*, J.-H. Lestang, C. de Mulatier, *Statistical Modeling of Community Structures in Binary Data*

## Different programs for searching for the best MCM:

The comparison between models is based on their evidence (posterior probability that the model generates the data, integrated over its parameter values). All these programs aim at finding the models with the largest log-evidence. 

**Exhaustive search:** The program available [here](https://github.com/clelidm/MinCompSpin) performs an exhaustive search for the best community. The code go through all possible MCMs of a given rank `r`, and finds the one with the largest evidence (exact optimal model). The limiting factor of this code is that the space of possible MCMs increases almost exponentially with the number `r` of variables (the number of MCMs on a chosen basis is given by the Bell number of `r`). We recommand using a different strategy when the number of variables exceeds `~13` to `15`.

**Greedy search VS Simulated annealing:** 
We recommand using the simulated annealing algorithm if you are interested in finding a solution that is as close as possible to the global optimal, or when the local greedy procedure fails to merge the initial spin variables into communities (this could be due to non-pairwise). While we recommand the Greedy algorithm if you are looking for a fast converging algorithm, as it will allow finding solution in reasonable time when the search space becomes too large.

## General information

The code performs an hierarchical merging procedure to find an optimal MCM in the basis chosen by the user.

**Chosen Basis:** To do so, the code start by re-writing the data in the basis provided by the user. The basis is encoded in the variable `Basis_li` and **by default it is the original basis of the data** (i.e., by default, there is no basis transformation). A different basis can be specified in a file by the user and read with the functions `Read_BasisOp_IntegerRepresentation()` or the function `Read_BasisOp_BinaryRepresentation()` (see examples in the `INPUT` folder, the files `SCOTUS_n9_BestBasis_Binary.dat` and `SCOTUS_n9_BestBasis_Integer.dat` which are encoding the best basis for the US Supreme Court data used as an example in Ref.[1]).

**Greedy Search:** The Greedy search can be done using the functions `MCM_GreedySearch()` or `MCM_GreedySearch_AND_printInfo()` (the only difference is the latter function also print some information about the found best MCM). To performs the hierarchical merging procedure, the code starts with `r` initial communities, and then successively merges the two communities that gives the largest increase in the log-Evidence of the model. The code iterates this procedure until no more communities can be merged without decreasing the log-Evidence. **By default, the `n` initial communities are taken to be such that there is only one variable in each community.**  This procedure generates an approximation of the MCM that achieves a maximal value of the evidence along the hierarchical merging process, as the number of communities varies from `n` to `1`.

**Starting from a chosen MCM:** The search can also be started from a different initial MCM, using the function `MCM_GreedySearch_MCM0()`. This can be used for instance to run a short greedy merging at the end of a simulated annealing search.

## Requirements

The code uses the C++11 version of C++.

## Usage with Makefile:

Run the following commands in your terminal, from the main folder (folder containing the `makefile` document):

 - **To compile:** `make`

 - **To Execute:** `make run` . This will use the datafile and variables that are specified in the makefile.

To change datafile: open the makefile and replace the values of the two following variables at the very top of the file (an example is provided):
>  - `datafile`: path to your own datafile;
>  - `n`: number of variables in your file; maximum possible value `n = 128`.

You can also execute the code by running in your terminal the command (from the main folder):
```bash
./GreedySearch.out  datafilename  n
```

where you must replace `datafilename` by the name of your datafile and `n` by your number of variables.

 - **To clean:** `make clean` (to use only once you're done using the code)

## Usage without Makefile:

 - **To compile:**  Type in your Terminal from the main folder:
```bash
cd Libraries/MCM/
g++ -std=c++11 -O2 -c *.cpp
cd ../../
g++ -std=c++11 -O2 -c Libraries/library.hpp
g++ -std=c++11 -O2 Libraries/main_routines.cpp main.cpp -include Libraries/library.hpp Libraries/MCM/*.o -o GreedySearch.out
```
or use the script `job_compile.sh` by running in your Terminal:
```bash
bash job_compile.sh
```

 - **To execute:**
```bash
./GreedySearch.out  datafilename  n
```

where you must replace `datafilename` by the name of your datafile and `n` by your number of variables.

or use the script `job_run.sh` by running in your Terminal:
```bash
bash job_run.sh
```
after replacing `datafilename` by the name of your datafile and `n` by your number of variables in the file `job_run.sh`.


## Examples

All the functions that can be called from `int main()` are declared at the beginning of the `main.cpp` file or in the file `library.hpp`. The most useful functions are described in the section "General information" above. 

In the input folder, we provided as an example the binary dataset `SCOTUS_n9_N895_Data.dat` analysed in the example, which is the dataset of votes of the US Supreme court analysed in Ref.[3] and  used as an example in Ref.[1]. For hands-on and simple tests of the program, please check the examples in the function `int main()` of the `main.cpp` file, and the usage of the program detailed in the function `tutorial()` of the main file (you can call this function by uncommenting the last lines in the function `main()`). 

In the input folder, we also provided as an example the binary dataset `Big5-IPC1_VS3_N5e4.dat`. This is a binarized version of the first `50 000` samples of the Big 5 dataset [4], which has `50` variables. This dataset is the one run as an example by the job file `job_run.sh`. See paper [1] for comments on the results obtained by the analysis using MCMs.

[3] E.D. Lee, C.P. Broedersz, W. Bialek, Statistical Mechanics of the US Supreme Court. [J Stat Phys 160, 275–301 (2015)](https://link.springer.com/article/10.1007/s10955-015-1253-6).

[4] Raw data from [Open-Source Psychometrics Project](https://openpsychometrics.org/_rawdata/) in the line indicated as "Answers to the IPIP Big Five Factor Markers"; [here](https://openpsychometrics.org/_rawdata/IPIP-FFM-data-8Nov2018.zip) is a direct link to the same zipfile.

## License

This code is an open source project under the GNU GPLv3.

----

## Important information:
### Basis change:
To change the basis of the data to a chosen basis and apply the MCM search in this new basis:
 1. Specify the basis elements in a list of integers `list<__int128_t> basis_li = ` using one of the available function.
 2. Transform the dataset `Nset` into the new basis (transformed data is in `Kset`) using the function `build_Kset`. This defines the new transformed data, stored in `Kset`:
```c++
vector<pair<__int128_t, unsigned int>> Kset = build_Kset(Nset, Basis_li);
```

**!! Important!!**
when performing this basis transformation, basis operators are placed from right to left in the new basis, 
i.e. the rightmost bit (lowest bit) corresponds to the first operator in `list<__int128_t> Basis`.

This very important for properly interpreting the output of the MCM algorithm after basis transformation.

