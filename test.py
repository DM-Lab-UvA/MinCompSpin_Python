#! /usr/bin/env python
import MinCompSpin as mod

def greedy_main(path, r):
    Kset, N = mod.read_datafile(path, r)
    print(f"Found {N} samples")
    print("Entropy (S):", mod.Entropy(Kset, N))

    print("Starting search")
    fp1 = mod.MCM_GreedySearch(Kset, N, r)
    print("model:", fp1)
    print("LogE_g:", mod.LogE_MCM(Kset, fp1, N, r))
    print("LogL_g:", mod.LogL_MCM(Kset, fp1, N, r))
    print("MCM_infoICC:", mod.LogE_MCM_infoICC(Kset, fp1, N, r))
    

def test_buf():
    import numpy as np
    dtype = np.dtype('uint64')
    print('dtype', dtype)
    arr = np.array([[1, 2, 3], [5, 6, 7]], dtype=np.uint64)
    print(1)
    print(arr)
    mod.test(arr)
    print(2)
    print(arr)
    mod.test(arr)
    print(3)
    print(arr)
    mod.test(arr)
    return

def test():
    path = "MinCompSpin_Greedy/INPUT/SCOTUS_n9_BestBasis_Binary.dat"
    r = 9
    greedy_main(path, r)
    return

    res = mod.read_datafile(path, r)
    print("read file")
    print(res)
    Kset, N = res

    print("starting search")
    # mod.MCM_GreedySearch_AND_printInfo(Kset, N, r)
    fp1 = mod.MCM_GreedySearch(Kset, N, r)
    print("did search")
    print(fp1)

    S = mod.Entropy(Kset, N)
    print("S", S)

    # print("LogE_g", mod.LogE_MCM(Kset, fp1, N, r))
    # print("LogL_g", mod.LogL_MCM(Kset, fp1, N, r))

    # print("Print_MCM_Partition", mod.Print_MCM_Partition(fp1, r))

    print("LogE_MCM_infoICC", mod.LogE_MCM_infoICC(Kset, fp1, N, r))

    mod.main("MinCompSpin_Greedy/INPUT/SCOTUS_n9_BestBasis_Binary.dat", "9")


if __name__ == '__main__':
    test_buf()
