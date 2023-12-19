#! /usr/bin/env python
import MinCompSpin as mod

def greedy_main(path, r):
    Kset_orig = mod.read_datafile(path, r)
    Kset = mod.Kset(Kset_orig.array, Kset_orig.N)
    print("Kset", type(Kset), Kset)
    print("Kset.array", type(Kset.array))
    print(f"Found {Kset.N} samples")
    print("Entropy (S):", mod.Entropy(Kset))

    print("Starting search")
    fp1 = mod.MCM_GreedySearch(Kset, r)
    print("model:", fp1)
    # Create a new MCM_Partitions instance with an array
    fp2 = mod.Partitions(fp1.array, r)
    print("LogE_g:", mod.LogE_MCM(Kset, fp2, r))
    print("LogL_g:", mod.LogL_MCM(Kset, fp2, r))
    print("MCM_infoICC:", mod.LogE_MCM_infoICC(Kset, fp2, r))

def test():
    path = "MinCompSpin_Greedy/INPUT/SCOTUS_n9_N895_Data.dat"
    r = 9

    greedy_main(path, r)
    # mod.main(path, "9")
    return

if __name__ == '__main__':
    test()
