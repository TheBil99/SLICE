import warnings
import pickle

from scipy.spatial.distance import squareform
from src.utilities import *
#from src.utilities import *

# Functions for pairwise SLICE

def scialdone_equation_for_s_mat(R_mat, K, F):
    """ Equation for computing Scialdone's s_mat using a generic R ( M2/(M1+M2) ) term. """
    
    numerator = 1 - 2 * (1 - F) ** (1 / K) + ((1 - 2 * F + R_mat) / (1 + R_mat)) ** (1 / K)
    denominator = (1 - (1 - F) ** (1 / K)) ** 2
    s_mat = numerator / denominator
    
    return s_mat


def compute_s_mat(F_arr, F_mat, K, F):
    """ Computes the SLICE-normalized Co-Segregation matrix. """

    F_i = np.tile(F_arr, (len(F_arr), 1))
    F_j = np.transpose(F_i)
    F_ij = F_mat
    R_ij = F_ij / (F_i + F_j - F_ij)
    s_mat = scialdone_equation_for_s_mat(R_ij, K, F)
    np.fill_diagonal(s_mat, np.nan)

    return s_mat


def compute_beta_mat(F_arr, F_mat, K, F):
    """ Computes beta matrix using Scialdone's method. """

    F_ij = F_mat
    F_i = np.tile(F_arr, (len(F_arr), 1))
    F_j = np.transpose(F_i)
    F_i_plus_F_j = F_i + F_j

    # Conditions of low statistic to stop the evaluation of beta at high genomic distances
    # If at a given genomic distance (i.e. diagonal) there are less than stop_threshold non-NaN values,
    # the algorithm stops performing the average for the current and next diagonals.
    # Instead, the value of beta for the last non-NaN diagonal is used.
    # Since the number of NaNs of F_ij and F_i + F_j is the same
    # (F_mat[:, np.isnan(F_arr)] = np.nan and F_mat[np.isnan(F_arr), :] = np.nan in the main)
    # we can evaluate the number of NaNs in both matrices.
    stop_threshold = 20  # 100
    stop_average = False
    last_beta = np.nan

    beta_mat = np.nan * np.ones(F_mat.shape)  # initialize with all nans

    # First, we perform a loop filling beta_mat on the upper-diagonals
    for diagonal_index in np.arange(1, len(beta_mat)):

        # First, if stop_average is True, beta doesn't have to be computed. Instead, the last non-NaN beta is used
        if stop_average:
            beta_mat[kth_diag_indices(beta_mat, diagonal_index)] = last_beta
            continue

        # If the diagonal has less than stop_threshold non-NaN values, beta doesn't have to be computed for this
        # diagonal and the next ones, so stop_average is set as True
        if np.sum(~np.isnan(np.diagonal(F_ij, diagonal_index))) < stop_threshold:
            stop_average = True
            beta_mat[kth_diag_indices(beta_mat, diagonal_index)] = last_beta
            continue

        F_ij_diagonal_mean = np.nanmean(np.diagonal(F_ij, diagonal_index))
        F_i_plus_F_j_diagonal_mean = np.nanmean(np.diagonal(F_i_plus_F_j, diagonal_index))
        R_ij_diagonal_mean = F_ij_diagonal_mean / (F_i_plus_F_j_diagonal_mean - F_ij_diagonal_mean)
        beta = scialdone_equation_for_s_mat(R_ij_diagonal_mean, K, F)
        last_beta = beta  # The last non-NaN beta value is updated
        beta_mat[kth_diag_indices(beta_mat, diagonal_index)] = beta

    # Now we fill the lower-diagonals symmetrically using the squareform function
    beta_mat = squareform(beta_mat, checks=False)
    beta_mat = squareform(beta_mat)

    # We set as 0 the entries where beta is negative, and as alpha the entries where beta > alpha
    beta_mat[beta_mat < 0] = 0
    beta_mat[beta_mat > alpha] = alpha

    return beta_mat


def compute_pi(s_mat, beta_mat):
    """ Computes the PI matrix. """

    pi_mat = (s_mat - beta_mat) / (alpha - beta_mat)
    np.fill_diagonal(pi_mat, np.nan)

    return pi_mat


def compute_pi_threshold(beta_mat, n_tubes, K, F, threshold, generator):
    """ Computes the PI threshold using Scialdone's original method. """

    R_threshold_mat = np.nan * np.ones(beta_mat.shape)

    # We loop on the upper-diagonals of R_threshold_mat and fill them using Scialdone's method
    for diagonal_index in np.arange(1, len(R_threshold_mat)):

        # Parameters for the multinomial distribution extraction
        n_multinomial_repetitions = 10000
        #threshold = 95

        # Take the beta value for the diagonal considered
        beta_diagonal = np.diagonal(beta_mat, diagonal_index)
        if np.sum(~np.isnan(beta_diagonal)) == 0:
            continue
        beta = np.nanmean(beta_diagonal)

        # Compute the M-values, i.e. the probabilities for the multinomial extraction
        M2 = -1 + 2 * F + (-1 + 2 * (1 - F) ** (1 / K) + beta * (1 - (1 - F) ** (1 / K)) ** 2) ** K
        M1 = 2 - 2 * F - 2 * (-1 + 2 * (1 - F) ** (1 / K) + beta * (1 - (1 - F) ** (1 / K)) ** 2) ** K
        M0 = 1 - M1 - M2
        if M2 < 0 or M2 > 1 or M1 < 0 or M1 > 1 or M0 < 0 or M0 > 1:
            print("invalid M-value value encountered")
            continue  # if the M values are not good, we skip this diagonal leaving it as NaN

        # Extract the T-values from the multinomial distribution
        # The if below allows the possibility to give a random generator as input to the function: can be useful for parallel applications
        if(generator == None):
            T_values = np.random.multinomial(n=n_tubes, pvals=[M2, M1, M0], size=n_multinomial_repetitions)
        else:
            T_values = generator.multinomial(n=n_tubes, pvals=[M2, M1, M0], size=n_multinomial_repetitions)
            
        
        T2_values = T_values[:, 0]
        T1_values = T_values[:, 1]

        # Compute the threshold R-value and put it in the respective diagonal of the threshold R-matrix
        R_values = T2_values / (T1_values + T2_values)
        R_threshold = np.percentile(R_values, threshold)
        R_threshold_mat[kth_diag_indices(R_threshold_mat, diagonal_index)] = R_threshold

    # Now we fill the lower-diagonals symmetrically using the square-form function
    R_threshold_mat = squareform(R_threshold_mat, checks=False)
    R_threshold_mat = squareform(R_threshold_mat)

    # Finally we compute the PI threshold using the PI functions
    s_threshold_mat = scialdone_equation_for_s_mat(R_threshold_mat, K, F)
    pi_threshold_mat = compute_pi(s_threshold_mat, beta_mat)

    return pi_threshold_mat


# GENERAL PAIRWISE SLICE FUNCTION

def SLICE_PAIRWISE(segregation_table, chr, generator,threshold = 95):
    """ Performs SLICE. """

    n_tubes = np.shape(segregation_table)[1]

    if chr == "chrX" or chr == "chrY":
        K = int(effective_NPs_per_tube / 2)
        F = 1 - (1 - F_mean) ** (1 / 2)
    else:
        K = effective_NPs_per_tube
        F = F_mean

    F_arr = compute_tube_segregation_frequency(segregation_table)
    F_mat = compute_tube_cosegregation_matrix(segregation_table)

    F_mat[:, np.isnan(F_arr)] = np.nan
    F_mat[np.isnan(F_arr), :] = np.nan

    s_mat = compute_s_mat(F_arr, F_mat, K, F)

    beta_mat = compute_beta_mat(F_arr, F_mat, K, F)
    beta_arr = beta_mat[0, :]

    pi_mat = compute_pi(s_mat, beta_mat)
    pi_mat[np.isinf(pi_mat)] = np.nan

    pi_threshold_mat = compute_pi_threshold(beta_mat, n_tubes, K, F, threshold, generator)
    pi_significant_mat = np.copy(pi_mat)
    pi_significant_mat[pi_mat < pi_threshold_mat] = np.nan

    return pi_mat, pi_significant_mat, beta_mat, beta_arr

# FUNCTION THAT EXECUTES SLICE ANALYSIS AND REMOVES VALUES OUTSIDE O-1 RANGE

def single_chromosome(chr, segregation_table,threshold = 95, verbose = True, save = False, ret = True, generator = None):
    """ Saves PI values for a specific chromosome.
        The PIs are saved as an horizontal array containing the upper-triangular elements of the matrix. """
    if(verbose == True):
        print("\nComputing PIs for " + chr + " ...")
    

    # We silence the warning given by the comparison of NaN values
    with warnings.catch_warnings():

        warnings.filterwarnings("ignore", "invalid value encountered in", category=RuntimeWarning)

        # Compute Co-Segregation Matrix, PI matrix and PI significant matrix
        pi_mat, pi_significant_mat, beta_mat, beta_arr = SLICE_PAIRWISE(segregation_table, chr, generator, threshold)
        
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "divide by zero encountered in log", category=RuntimeWarning)      
            npmi_mat = compute_npmi(segregation_table)

        np.savetxt(data_path + name_root + "/beta/beta_arr_" + chr + "_" + name_root + ".txt", beta_arr)

        # Print percentages of PI NaN, < 0 and > 1
        pi_upper_triangular_arr = squareform(pi_mat, checks=False)
        pi_upper_triangular_arr_not_nan = pi_upper_triangular_arr[~np.isnan(pi_upper_triangular_arr)]
        nan_percentage = 100 - 100 * len(pi_upper_triangular_arr_not_nan) / len(pi_upper_triangular_arr)
        pi_less_zero_percentage = 100 * len(pi_upper_triangular_arr_not_nan[pi_upper_triangular_arr_not_nan < 0]) / len(pi_upper_triangular_arr_not_nan)
        pi_bigger_one_percentage = 100 * len(pi_upper_triangular_arr_not_nan[pi_upper_triangular_arr_not_nan > 1]) / len(pi_upper_triangular_arr_not_nan)
        significant_pi_percentage = 100 * np.sum(~np.isnan(pi_significant_mat)) / np.sum(~np.isnan(pi_mat))
        
        if(verbose == True):
            print("% of NaN PI values: " + str(nan_percentage) + " %")
            print("% of PI < 0 (out of non-nans): " + str(pi_less_zero_percentage) + " %")
            print("% of PI > 1 (out of non-nans): " + str(pi_bigger_one_percentage) + " %")
            print("Percentage of significative PIs: " + str(significant_pi_percentage) + " %")

        # Sort out problematic PI values
        pi_mat[pi_mat < 0] = np.nan
        pi_mat[pi_mat > 1] = 1.
        pi_significant_mat[pi_significant_mat > 1] = 1.

        # Print additional info
        mean_pi = np.nanmean(pi_mat)
        std_pi = np.nanstd(pi_mat)
        mean_significant_pi = np.nanmean(pi_significant_mat)
        std_significant_pi = np.nanstd(pi_significant_mat)

        if(verbose == True):
            print("Mean PI value: " + str(mean_pi))
            print("STD PI value: " + str(std_pi))
            print("Mean SIGNIFICATIVE PI value: " + str(mean_significant_pi))
            print("STD SIGNIFICATIVE PI value: " + str(std_significant_pi))

    # Convert matrices to upper-triangular arrays
    pi_arr = squareform(pi_mat, checks=False)
    pi_significant_arr = squareform(pi_significant_mat, checks=False)
    npmi_arr = squareform(npmi_mat, checks=False)
    # Reshape the arrays into a matrix with one row.
    pi_arr = np.reshape(pi_arr, (1, len(pi_arr)))
    pi_significant_arr = np.reshape(pi_significant_arr, (1, len(pi_significant_arr)))
    npmi_arr = np.reshape(npmi_arr, (1, len(npmi_arr)))

    # Save arrays as txt files

    if(save == True):
        if(verbose == True):
            print("Saving data ...")
        np.savetxt(data_path + name_root + "/PI2/PI2_" + chr + "_" + name_root + ".txt", pi_arr)
        np.savetxt(data_path + name_root + "/PI2/PI2_significant_" + str(threshold) + "_" + chr + "_" + name_root + ".txt", pi_significant_arr)
        np.savetxt(data_path + name_root + "/NPMI/NPMI_" + chr + "_" + name_root + ".txt", npmi_arr)

    
    if(verbose == True):
        print("Done\n")

    if(ret == False):
        del pi_mat
        del pi_significant_mat
        del beta_mat
        del beta_arr
        del pi_upper_triangular_arr
        del pi_upper_triangular_arr_not_nan
        del pi_arr
        del pi_significant_arr

        return None
    if(ret == True):
        return pi_mat, pi_significant_mat

