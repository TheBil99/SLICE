import warnings

from src.utilities import *

#from scipy.spatial.distance import squareform


#########################################################################
#########################################################################
# Versione di interSLICE che lavora solo con la parte off diagonal di AB #
#########################################################################
#########################################################################


def equation_for_s_mat(R_mat, K, F):
    """ Equation for computing Scialdone's s_mat using a generic R ( M2/(M1+M2) ) term. """
    
    # numerator =  - (1 - F) ** (1 / K) + ( (   ((1 - 2 * F + R_mat) / (1 + R_mat))  ) )**(1/2*K)
    # denominator = (1 - (1 - F) ** (1 / 2*K)) ** 2
    # s_mat = numerator / denominator

    numerator = 1 - 2 * (1 - F) ** (1 / K) + ((1 - 2 * F + R_mat) / (1 + R_mat)) ** (1 / K)
    denominator = (1 - (1 - F) ** (1 / K)) ** 2
    s_mat = numerator / denominator
    
    
    return s_mat


def compute_s_mat_offdiag(F_arr_A, F_arr_B, F_mat, K, F):
    """ Computes the SLICE-normalized Co-Segregation matrix. """

    F_i = np.tile(F_arr_A, (len(F_arr_B), 1)).T
    F_j = np.tile(F_arr_B, (len(F_arr_A), 1))
    F_ij = F_mat
    R_ij = F_ij / (F_i + F_j - F_ij)
    s_mat = equation_for_s_mat(R_ij, K, F)
    

    return s_mat


def compute_pi_offdiag(s_mat, beta):
    """ Computes the PI matrix. """
    # pi_mat = s_mat / (alpha - 1)
    pi_mat = (1/2) *  (s_mat - beta) / (alpha - beta)

    return pi_mat

def compute_interpi_threshold(n_tubes, K, F, threshold,  beta ,generator):
    """ Computes the PI threshold using Scialdone's original method. """

    R_threshold = np.nan

    # Parameters for the multinomial distribution extraction
    n_multinomial_repetitions = 10000


    # Compute the M-values, i.e. the probabilities for the multinomial extraction

    M2 = -1 + 2 * F + (-1 + 2 * (1 - F) ** (1 / K) + beta*(1 - (1 - F) ** (1 / K)) ** 2) ** K
    M1 = 2 - 2 * F - 2 * (-1 + 2 * (1 - F) ** (1 / K) + beta*(1 - (1 - F) ** (1 / K)) ** 2) ** K
    M0 = 1 - M1 - M2

    if M2 < 0 or M2 > 1 or M1 < 0 or M1 > 1 or M0 < 0 or M0 > 1:
        raise ValueError('Values outside of [0,1] for M0, M1, M2')

    # Extract the T-values from the multinomial distribution
    # The if below allows the possibility to give a random generator as input to the function: can be useful for parallel applications
    if(generator == None):
        T_values = np.random.multinomial(n=n_tubes, pvals=[M2, M1, M0], size=n_multinomial_repetitions)
    else:
        T_values = generator.multinomial(n=n_tubes, pvals=[M2, M1, M0], size=n_multinomial_repetitions)
    
    T2_values = T_values[:, 0]
    T1_values = T_values[:, 1]

    # Compute the threshold R-value 
    R_values = T2_values / (T1_values + T2_values)
    R_threshold = np.percentile(R_values, threshold)



    # Finally we compute the PI threshold using the PI functions
    s_threshold = equation_for_s_mat(R_threshold, K, F)
    pi_threshold = compute_pi_offdiag(s_threshold, beta)

    return pi_threshold

# GENERAL PAIRWISE SLICE FUNCTION BETWEEN TWO CHROMOSOMES

def SLICE_PAIRWISE_INTER(segregation_table_A, segregation_table_B, chrA, chrB,  generator,  threshold = 95):
    """ Performs SLICE. """

    n_tubes = np.shape(segregation_table_A)[1]


    if chrA == "chrX" or chrA == "chrY" or chrB == "chrX" or chrB == "chrY":
        raise ValueError('The interSLICE algorithm has not still extended to sexual chromosomes')
    else:
        K = effective_NPs_per_tube
        F = F_mean

    F_arr_A, F_arr_B = compute_tube_segregation_frequency_offdiag(segregation_table_A, segregation_table_B)
    F_mat = compute_tube_cosegregation_matrix_offdiag(segregation_table_A, segregation_table_B)

    F_mat[:, np.isnan(F_arr_B)] = np.nan
    F_mat[np.isnan(F_arr_A), :] = np.nan

    s_mat = compute_s_mat_offdiag(F_arr_A, F_arr_B, F_mat, K, F)

    ######
    beta = np.nanmean(s_mat)
    ######

    pi_mat = compute_pi_offdiag(s_mat, beta = beta)
    pi_mat[np.isinf(pi_mat)] = np.nan

    pi_threshold = compute_interpi_threshold(n_tubes, K, F, threshold, beta = beta,generator = generator)
    pi_significant_mat = np.copy(pi_mat)
    pi_significant_mat[pi_mat < pi_threshold] = np.nan

    return pi_mat, pi_threshold, pi_significant_mat, F_mat, beta

# FUNCTION THAT EXECUTES SLICE ANALYSIS AND REMOVES VALUES OUTSIDE O-1 RANGE

def inter_chromosome(segregation_table_A, segregation_table_B, chrA, chrB, threshold = 95, verbose = True, save = False, ret = True, generator = None):
    # calculates pis for interchromosomal contacts and removes values outside [0,1] interval
    if(verbose == True):
        print("\nComputing inter PIs for " + chrA + " and " + chrB + " ...")


    # We silence the warning given by the comparison of NaN values
    with warnings.catch_warnings():

        warnings.filterwarnings("ignore", "invalid value encountered in", category=RuntimeWarning)

        # Compute Co-Segregation Matrix, PI matrix and PI significant matrix
        pi_mat, _, pi_significant_mat, cosegregation_mat, beta = SLICE_PAIRWISE_INTER(segregation_table_A, segregation_table_B, chrA, chrB,  generator ,threshold)
        #npmi_mat = compute_npmi(segregation_table)

        # Print percentages of PI NaN, < 0 and > 1
        pi_mat_not_nan = pi_mat[~np.isnan(pi_mat)]

        nan_percentage               = 100 - 100 * pi_mat_not_nan.size / pi_mat.size
        pi_less_zero_percentage      = 100 * pi_mat_not_nan[pi_mat_not_nan < 0].size /  pi_mat_not_nan.size
        pi_bigger_one_percentage     = 100 * pi_mat_not_nan[pi_mat_not_nan > 1].size / pi_mat_not_nan.size
        significant_pi_percentage    = 100 * np.sum(~np.isnan(pi_significant_mat)) / np.sum(~np.isnan(pi_mat))

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
            print("beta: " + str(beta))


    # Save arrays as npy files (binary files)
    if(save == True):
        if(verbose == True):
            print("Saving data ...")

        np.save(data_path + name_root + "/PI2_inter_beta_evaluation/PI2_inter_" + chrA + "_" + chrB + "_" + name_root + ".npy", pi_mat)
        np.save(data_path + name_root + "/PI2_inter_beta_evaluation/PI2_inter_significant_" + str(threshold) + "_" + chrA + "_" + chrB + "_" + name_root + ".npy", pi_significant_mat)
        np.save(data_path + name_root + "/cosegregation_matrix_inter_beta_evaluation/cosegregation_matrix_inter_" + chrA + "_" + chrB + "_" + name_root + ".npy", cosegregation_mat)
        
        #np.save(data_path + name_root + "/NPMI_inter/NPMI_" + chrA + "_" + chrB + "_" + name_root + ".npy", npmi_arr)

    if(verbose == True):
        print("Done\n")

    
    if(ret == False):
        del pi_mat
        del pi_significant_mat
        del pi_mat_not_nan
        del cosegregation_mat
        return None
    if(ret == True):
        return pi_mat, pi_significant_mat
