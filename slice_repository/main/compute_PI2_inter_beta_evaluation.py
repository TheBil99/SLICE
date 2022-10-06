from numpy import single
import pandas as pd
from scipy.spatial.distance import squareform
from numpy.random import Generator, MT19937

from settings import name_root, print_info, data_path
from src.slice_pairwise import single_chromosome
from src.slice_pairwise_inter_beta_evaluation import inter_chromosome

# Computes the PI2s in intra and in inter for the given chromosomes
print_info()

# Load Segregation data
df_segregation = pd.read_pickle(data_path + name_root + "/segregation_" + name_root + ".pkl")

# CALL FUNCTIONS HERE

# indicate in the following lines which chromosomes to run interSLICE on

n_chromosomes = 19
chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19"]

# the following lines allow to chose the random generator and a seed for the random generation
# setting rg to None to use the default numpy random generator, with no need to specify the seed 

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