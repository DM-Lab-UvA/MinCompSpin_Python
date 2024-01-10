#! /usr/bin/env python3
import timeit

debug = True

setup = '''
import MinCompSpin as mod
def greedy_main(path, r):
    Kset = mod.read_datafile(path, r)

    fp1 = mod.MCM_GreedySearch(Kset, r)
    S = mod.Entropy(Kset)
    LogL_g = mod.LogL_MCM(Kset, fp1, r)
    LogE_g = mod.LogE_MCM(Kset, fp1, r)
    LogL_infoICC = mod.LogL_MCM_infoICC(Kset, fp1, r)
    LogE_infoICC = mod.LogE_MCM_infoICC(Kset, fp1, r)
'''

stmt = '''
greedy_main("MinCompSpin_Greedy/INPUT/SCOTUS_n9_N895_Data.dat", 9)
'''
stmt = '''
greedy_main("MinCompSpin_Greedy/INPUT/Big5-reduced.dat", 50)
'''

def main():
    repeat = 5
    number = 1000
    if debug:
         number = 1
    res = timeit.repeat(stmt, setup, repeat=repeat, number=number)
    print("min", min(res) / number * 1000_000, "micro seconds")
    print("mean", sum(res) / len(res) / number * 1000_000, "micro seconds")

if __name__ == '__main__':
    main()
