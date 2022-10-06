import pandas as pd
from settings import *
from src.slice_pairwise import compute_tube_segregation_frequency

# Computes the mean tube segregation frequency of all autosomes (i.e. excluding sexual chromosomes)
# and prints the result on the standard output

print("Computing F_mean genomewide for " + name_root)
print("\n\n")

# Load Segregation data
df_segregation = pd.read_pickle(data_path + name_root + "/segregation_" + name_root + ".pkl")

# Remove chrX from this analysis
chr_list.remove("chrX")

# Array to store all F values
F_arr_genomewide = np.array([])
# For each chromosome, compute F_arr and store it in F_arr_genomewide
for chr in chr_list:
    print("     evaluating " + chr)
    segregation_table = df_segregation[chr]["segregation_table"]
    F_arr = compute_tube_segregation_frequency(segregation_table)
    F_arr_genomewide = np.concatenate((F_arr_genomewide, F_arr))
    del F_arr
print("\n\n")

# Print result
F_mean = np.nanmean(F_arr_genomewide)
print("F_mean = " + str(F_mean))
print("\n\n")

