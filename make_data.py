data_path_prefix = "MinCompSpin_Greedy/INPUT/"
base_files = [
    ["my_data_n100_N2000.dat", 100, 2000],
    ["my_data_n100_N1000.dat", 100, 1000],
    ["my_data_n20_N1000.dat",  20,  1000],
    ["my_data_n40_N1000.dat",  40,  1000],
    ["my_data_n60_N1000.dat",  60,  1000]
]

def make_alt(base, alt_N):
    with open(data_path_prefix + base[0], 'r') as base_file:
        with open(data_path_prefix + f"my_data_n{base[1]}_N{alt_N}.dat", "w") as new_file:
            line_count = 0
            for line in base_file.readlines():
                if line_count == alt_N:
                    break
                new_file.write(line)
                line_count += 1

def main():
    for base in base_files[1:]:
        for alt_N in [750, 500, 250, 100, 10]:
            make_alt(base, alt_N)

if __name__ == '__main__':
    main()
