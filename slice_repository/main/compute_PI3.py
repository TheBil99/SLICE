import pandas as pd
from src.slice_triplewise import *
import warnings

# Computes the PI3s for a chromosome

print_info()

# Load Segregation data
df_segregation = pd.read_pickle(data_path + name_root + "/segregation_" + name_root + ".pkl")


def single_chromosome(chr):
    """ Compute the PI3s for a specific chromosome.  """

    print("\nComputing PI3s for " + chr + " ...")

    segregation_table = np.float32(df_segregation[chr]["segregation_table"])

    # Compute PI3s
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "invalid value encountered in", category=RuntimeWarning)

        # Compute the PI3 tensors
        pi3_tensor, pi3_min_tensor, pi3_max_tensor, problematic_tensor = SLICE_TRIPLEWISE(segregation_table, chr)
        del segregation_table

        # Print percentage of triplets from problematics PI2 values
        print_problematic_tensor(problematic_tensor)

        # NO MIN-MAX CONDITIONS IMPOSED ON PI3
        # Print data
        print("Printing PI3 data info before conditions:")
        print_data_info(pi3_tensor)
        print("Percentage of PI3 between -1 and 2 (out of not-NaN): " + str(100 * compute_percentage_of_data_between_values(pi3_tensor, -1, 2)) + " %")
        print("\n\n")

        # IMPOSING THE MIN-MAX CONDITIONS ON PI3
        warnings.filterwarnings("ignore", "invalid value encountered in", category=RuntimeWarning)
        np.putmask(pi3_tensor, pi3_tensor < pi3_min_tensor, pi3_min_tensor)
        np.putmask(pi3_tensor, pi3_tensor > pi3_max_tensor, pi3_max_tensor)
        np.putmask(pi3_tensor, pi3_min_tensor == pi3_max_tensor, pi3_max_tensor)
        pi3_tensor[pi3_min_tensor > pi3_max_tensor] = np.nan
        # Print data
        print("Printing PI3 data info after conditions (including zeros):")
        print_data_info(pi3_tensor)
        print("Printing PI3 min data:")
        print_data_info(pi3_min_tensor)
        print("Printing PI3 max data:")
        print_data_info(pi3_max_tensor)

        # REMOVE THE PI3=0 VALUES
        pi3_tensor[pi3_tensor == 0] = np.nan
        print("Printing PI3 data info after conditions (removing zeros):")
        print_data_info(pi3_tensor)
        # Save the compressed PI3 tensor as compressed numpy array
        pi3_compressed_tensor = compress_symmetric_tensor(pi3_tensor)
        print("Number of good PI3s: " + str(pi3_compressed_tensor.shape[1]) + "\n\n")
        np.savez_compressed(data_path + name_root + "/PI3/PI3_tensor_" + chr + "_" + name_root, pi3_compressed_tensor)
        del pi3_compressed_tensor

    del pi3_tensor
    del pi3_min_tensor
    del pi3_max_tensor
    del problematic_tensor

    return None


# CALL FUNCTIONS HERE

single_chromosome("chr19")
