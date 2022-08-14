from email import generator
from numpy import single
import pandas as pd
import sys
from scipy.spatial.distance import squareform
import multiprocessing as mp
from numpy.random import Generator, MT19937

from settings import name_root, print_info, data_path
from src.slice_pairwise import single_chromosome
from src.slice_pairwise_inter_beta_evaluation import inter_chromosome

# Computes the PI2s for a single chromosome

# Load Segregation data
df_segregation = pd.read_pickle(data_path + name_root + "/segregation_" + name_root + ".pkl")


start_chr = 1;   stop_chr = 19

stop_chr += 1
chromosomes = []
n_chromosomes = stop_chr - start_chr

   
for i in range(start_chr, stop_chr):
    chrA = 'chr' + str(i)
    chromosomes += [chrA]

s = 3545135
rg = Generator(MT19937(s))

print("Computing pairwise pi matrices in intra and in inter for dataset " , name_root, "\n")

for i in range(n_chromosomes):
    for j in range(i, n_chromosomes):
        chrA = chromosomes[i]
        chrB = chromosomes[j]
        if(chrA!=chrB):
            segregation_table_A = df_segregation[chrA]["segregation_table"]
            segregation_table_B = df_segregation[chrB]["segregation_table"]
            inter_chromosome(segregation_table_A, segregation_table_B, chrA, chrB, 95, False,  True, False, generator = rg)
            del segregation_table_A, segregation_table_B
        else:
            segregation_table = df_segregation[chrA]["segregation_table"]
            single_chromosome(chrA, segregation_table,threshold = 95, verbose = False, save = True, ret = False, generator = rg)

del df_segregation