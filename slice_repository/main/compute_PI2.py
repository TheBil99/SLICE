import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from settings import name_root, print_info, data_path
from src.slice_pairwise import single_chromosome
import sys

# Computes the PI2s for a single chromosome
print_info()

# Load Segregation data
df_segregation = pd.read_pickle(data_path + name_root + "/segregation_" + name_root + ".pkl")

# CALL FUNCTIONS HERE

if( len(sys.argv) ==1  ):
    print("No chromosome passed, computing default chromosome: chr19\n")
    segregation_table = df_segregation['chr19']["segregation_table"]
    del df_segregation
    single_chromosome("chr19", segregation_table, verbose = True, save = True, ret = False)
if( len(sys.argv) ==2  ):
    segregation_table = df_segregation["chr" + sys.argv[1]]["segregation_table"]
    del df_segregation
    single_chromosome("chr" + sys.argv[1], segregation_table, verbose = True, save = True, ret = False)
if( len(sys.argv) ==3  ):
    segregation_table = df_segregation["chr" + sys.argv[1]]["segregation_table"]
    del df_segregation
    single_chromosome("chr" + sys.argv[1], segregation_table, threshold = float(sys.argv[2]), verbose = True, save = True, ret = False )


if(len(sys.argv)>3  ):
    print("error: too many arguments passed \n")

