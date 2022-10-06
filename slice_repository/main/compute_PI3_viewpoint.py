import pandas as pd
from src.slice_triplewise_viewpoint import *
from matplotlib import pyplot as plt
from matplotlib import cm
from scipy.spatial.distance import squareform
import warnings

# Computes the PI3s fixing a viewpoint for a viewpoint of interest

print_info()

# Load Segregation data
df_segregation = pd.read_pickle(data_path + name_root + "/segregation_" + name_root + ".pkl")


def chromosome_fixed_i(chr, viewpoint_coordinate):
    """ Computes the PI3s for a specific chromosome fixing a reference window i.
        Attention: the input reference_coordinate must be in bp. """

    print("Computing PIs for " + chr + " fixing the window at the genomic position " + str(viewpoint_coordinate) + " bp...")
    print("\n")

    segregation_table = df_segregation[chr]["segregation_table"]

    # Compute i
    i = int(np.floor(viewpoint_coordinate / resolution))
    print("Reference window i: " + str(i))
    print("\n")

    # Compute PI3s with warnings silenced
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "invalid value encountered in", category=RuntimeWarning)
        pi3_viewpoint_mat, pi3_min_viewpoint_mat, pi3_max_viewpoint_mat = SLICE_TRIPLEWISE_VIEWPOINT(i, segregation_table, chr)

    print("Printing data info before setting all PI3 conditions:")
    print_data_info(pi3_viewpoint_mat)

    # Impose all PI3 conditions
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "invalid value encountered in", category=RuntimeWarning)
        np.putmask(pi3_viewpoint_mat, pi3_viewpoint_mat < pi3_min_viewpoint_mat, pi3_min_viewpoint_mat)
        np.putmask(pi3_viewpoint_mat, pi3_viewpoint_mat > pi3_max_viewpoint_mat, pi3_max_viewpoint_mat)
        np.putmask(pi3_viewpoint_mat, pi3_min_viewpoint_mat == pi3_max_viewpoint_mat, pi3_max_viewpoint_mat)
        pi3_viewpoint_mat[pi3_min_viewpoint_mat > pi3_max_viewpoint_mat] = np.nan
        pi3_viewpoint_mat[pi3_viewpoint_mat == 0] = np.nan
    del pi3_min_viewpoint_mat
    del pi3_max_viewpoint_mat

    print("Printing data info after setting all PI3 conditions:")
    print_data_info(pi3_viewpoint_mat)

    # Save the upper triangular of the matrix as horizontal array
    pi3_viewpoint_arr = squareform(pi3_viewpoint_mat, checks=False)
    pi3_viewpoint_arr = np.reshape(pi3_viewpoint_arr, (1, len(pi3_viewpoint_arr)))
    np.savetxt(data_path + name_root + "/PI3_viewpoint/PI3_viewpoint_" + viewpoint_name + "_" + name_root + ".txt",
               pi3_viewpoint_arr)

    return None


# CALL FUNCTIONS HERE

viewpoint_name = "Lama1"
chr = "chr17"
#position = 68046604        #Lama1 position in genomic assembly mm9
position = 67697259  # in bp    <-lama1 position in genomic assembly mm10

chromosome_fixed_i(chr, position)
