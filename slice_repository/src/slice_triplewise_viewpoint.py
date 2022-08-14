from src.utilities import *
from src.slice_pairwise import SLICE_PAIRWISE

# Functions for three-way SLICE fixing a viewpoint


def compute_tube_co3segregation_viewpoint_matrix(i, segregation_table):
    """ Computes the tube Co-3-Segregation matrix, fixing a viewpoint window i, from the Segregation table. """

    n_windows = np.shape(segregation_table)[0]
    n_tubes = np.shape(segregation_table)[1]

    s_i = segregation_table[i, :]
    s_i_tile = np.tile(s_i, (n_windows, 1))

    F3_viewpoint_mat = np.matmul(segregation_table * s_i_tile, np.transpose(segregation_table)) / n_tubes

    np.fill_diagonal(F3_viewpoint_mat, np.nan)
    F3_viewpoint_mat[i, :] = np.nan
    F3_viewpoint_mat[:, i] = np.nan

    del s_i
    del s_i_tile

    return F3_viewpoint_mat


def compute_s3_viewpoint_matrix(i, F1_arr, F2_mat, F3_viewpoint_mat, eta_mat, K, F):
    """ Computes the SLICE-Normalized Co-3-Segregation matrix, fixing a viewpoint i. """

    # Quantities related to F1_arr
    p = 1 - (1 - F) ** (1 / K)
    F_i = F1_arr[i]
    F_j = np.tile(F1_arr, (len(F1_arr), 1))
    F_k = np.transpose(F_j)

    # Quantities related to F2_mat
    F_jk = F2_mat
    F_ij = np.tile(F2_mat[i, :], (len(F2_mat[i, :]), 1))
    F_ik = np.transpose(F_ij)

    # Quantity related to F3_viewpoint_mat
    F_ijk = F3_viewpoint_mat

    # Quantities related to eta_mat
    eta_jk = eta_mat
    eta_ij = np.tile(eta_mat[i, :], (len(eta_mat[i, :]), 1))
    eta_ik = np.transpose(eta_ij)

    term_1 = 1 - 3 * p + (p ** 2) * (eta_ij + eta_ik + eta_jk)
    term_2 = - 1 + F_ijk / (F_ij + F_ik + F_jk - F_i - F_j - F_k)
    term_3 = 3 * (1 - p) ** K - (1 - 2 * p + (p ** 2) * eta_ij) ** K - (1 - 2 * p + (p ** 2) * eta_ik) ** K - (1 - 2 * p + (p ** 2) * eta_jk) ** K

    s3_viewpoint_mat = (term_1 - (1 + term_2 * term_3) ** (1 / K)) / (p ** 3)

    del F_i
    del F_j
    del F_k
    del F_ij
    del F_ik
    del F_jk
    del F_ijk
    del eta_ij
    del eta_ik
    del eta_jk
    del term_1
    del term_2
    del term_3

    return s3_viewpoint_mat


def compute_beta3_viewpoint_matrix(i, beta2_mat):
    """ Computes beta3 matrix, fixing a viewpoint i. """

    j = np.tile(np.arange(len(beta2_mat)), (len(beta2_mat), 1))
    k = np.transpose(j)

    beta2_jk = beta2_mat
    beta2_ij = np.tile(beta2_mat[i, :], (len(beta2_mat[i, :]), 1))
    beta2_ik = np.transpose(beta2_ij)

    beta3_viewpoint_mat = np.nan * np.ones(beta2_mat.shape)

    # i is between j and k
    np.putmask(beta3_viewpoint_mat, np.logical_and(j < i, i < k), beta2_ij * beta2_ik)
    np.putmask(beta3_viewpoint_mat, np.logical_and(k < i, i < j), beta2_ij * beta2_ik)
    # j is between i and k
    np.putmask(beta3_viewpoint_mat, np.logical_and(i < j, j < k), beta2_ij * beta2_jk)
    np.putmask(beta3_viewpoint_mat, np.logical_and(k < j, j < i), beta2_ij * beta2_jk)
    # k is between i and j
    np.putmask(beta3_viewpoint_mat, np.logical_and(i < k, k < j), beta2_ik * beta2_jk)
    np.putmask(beta3_viewpoint_mat, np.logical_and(j < k, k < i), beta2_ik * beta2_jk)

    del j
    del k
    del beta2_ij
    del beta2_ik
    del beta2_jk

    return beta3_viewpoint_mat


def compute_gammas_viewpoint_matrices(i, beta2_mat):
    """ Computes gamma_ij, gamma_ik and gamma_jk matrices, fixing a viewpoint i. """

    # The order here is chosen so that the gamma matrices are correct
    k = np.tile(np.arange(len(beta2_mat)), (len(beta2_mat), 1))
    j = np.transpose(k)

    beta2_jk = beta2_mat
    beta2_ik = np.tile(beta2_mat[i, :], (len(beta2_mat[i, :]), 1))
    beta2_ij = np.transpose(beta2_ik)

    # gamma_ij
    gamma_ij_viewpoint_mat = np.nan * np.ones(beta2_mat.shape)
    np.putmask(gamma_ij_viewpoint_mat, np.abs(i - k) <= np.abs(j - k), beta2_ik / v)
    np.putmask(gamma_ij_viewpoint_mat, np.abs(i - k) > np.abs(j - k), beta2_jk / v)

    # gamma ik
    gamma_ik_viewpoint_mat = np.nan * np.ones(beta2_mat.shape)
    np.putmask(gamma_ik_viewpoint_mat, np.abs(i - j) <= np.abs(k - j), beta2_ij / v)
    np.putmask(gamma_ik_viewpoint_mat, np.abs(i - j) > np.abs(k - j), beta2_jk / v)

    # gamma_jk
    gamma_jk_viewpoint_mat = np.nan * np.ones(beta2_mat.shape)
    np.putmask(gamma_jk_viewpoint_mat, np.abs(j - i) <= np.abs(k - i), beta2_ij / v)
    np.putmask(gamma_jk_viewpoint_mat, np.abs(j - i) > np.abs(k - i), beta2_ik / v)

    del j
    del k
    del beta2_ij
    del beta2_ik
    del beta2_jk

    return gamma_ij_viewpoint_mat, gamma_ik_viewpoint_mat, gamma_jk_viewpoint_mat


def compute_viewpoint_matrices(i, pi2_mat, beta2_mat):
    """ Computes beta3_viewpoint_mat, gamma_viewpoint_mat, pi2_sum_viewpoint_mat, omega_viewpoint_mat,
        pi3_min_viewpoint_mat, pi3_max_viewpoint_mat, fixing a viewpoint i, where
                gamma_viewpoint_mat is  as gamma_ij + gamma_ik + gamma_jk,
                pi2_sum_viewpoint_mat is pi_ij + pi_ik + pi_jk,
                omega_viewpoint_mat is pi_ij * gamma_ij + pi_ik * gamma_ik + pi_jk * gamma_jk. """

    # The order here is chosen so that omega_viewpoint_mat is correct
    pi2_jk = pi2_mat
    pi2_ik = np.tile(pi2_mat[i, :], (len(pi2_mat[i, :]), 1))
    pi2_ij = np.transpose(pi2_ik)

    gamma_ij_viewpoint_mat, gamma_ik_viewpoint_mat, gamma_jk_viewpoint_mat = compute_gammas_viewpoint_matrices(i, beta2_mat)

    # Computation of beta3_viewpoint_mat, gamma_viewpoint_mat, pi2_sum_viewpoint_mat, omega_viewpoint_mat

    beta3_viewpoint_mat = compute_beta3_viewpoint_matrix(i, beta2_mat)
    gamma_viewpoint_mat = gamma_ij_viewpoint_mat + gamma_ik_viewpoint_mat + gamma_jk_viewpoint_mat
    pi2_sum_viewpoint_mat = pi2_ij + pi2_ik + pi2_jk
    omega_viewpoint_mat = pi2_ij * gamma_ij_viewpoint_mat + pi2_ik * gamma_ik_viewpoint_mat + pi2_jk * gamma_jk_viewpoint_mat

    # Computation of pi3_min_viewpoint_mat, pi3_max_viewpoint_mat

    pi3_min_viewpoint_mat = (pi2_ij + pi2_ik + pi2_jk - 1) / 2
    pi3_min_viewpoint_mat[pi3_min_viewpoint_mat < 0] = 0

    pi2_ij[np.isnan(pi2_ij)] = np.inf
    pi2_ik[np.isnan(pi2_ik)] = np.inf
    pi2_jk[np.isnan(pi2_jk)] = np.inf
    pi3_max_viewpoint_mat = np.minimum(pi2_ij, pi2_ik)
    pi3_max_viewpoint_mat = np.minimum(pi3_max_viewpoint_mat, pi2_jk)
    pi3_max_viewpoint_mat[np.isinf(pi3_max_viewpoint_mat)] = np.nan

    # These lines of code implement the condition that when one between pi2_ij pi2_ik and pi2_jk is 1,
    # if the difference between the other two pi2s is less than a threshold (set as 0.1),
    # the corresponding pi3_ijk is set as the average value of these two pi2s.
    # To do so, here we impose that when one of the three pi2s is 1 (e.g. pi2_ij = 1)
    # and the difference between the other two is less than the threshold (e.g. |pi2_ik - pi2_jk| < 0.1),
    # then pi3_min_ijk = pi3_max_ijk = (pi2_ik + pi2_jk) / 2.
    # Then, in the main code compute_pi_triplewise, we have the condition that when pi3_min = pi3_max,
    # we set pi3 = pi3_min = pi3_max.
    pi2_ij[np.isinf(pi2_ij)] = np.nan
    pi2_ik[np.isinf(pi2_ik)] = np.nan
    pi2_jk[np.isinf(pi2_jk)] = np.nan
    threshold = 0.1
    np.putmask(pi3_min_viewpoint_mat, np.logical_and(pi2_ij == 1, np.abs(pi2_ik - pi2_jk) < threshold), (pi2_ik + pi2_jk) / 2)
    np.putmask(pi3_max_viewpoint_mat, np.logical_and(pi2_ij == 1, np.abs(pi2_ik - pi2_jk) < threshold), pi3_min_viewpoint_mat)
    np.putmask(pi3_min_viewpoint_mat, np.logical_and(pi2_ik == 1, np.abs(pi2_ij - pi2_jk) < threshold), (pi2_ij + pi2_jk) / 2)
    np.putmask(pi3_max_viewpoint_mat, np.logical_and(pi2_ik == 1, np.abs(pi2_ij - pi2_jk) < threshold), pi3_min_viewpoint_mat)
    np.putmask(pi3_min_viewpoint_mat, np.logical_and(pi2_jk == 1, np.abs(pi2_ij - pi2_ik) < threshold), (pi2_ij + pi2_ik) / 2)
    np.putmask(pi3_max_viewpoint_mat, np.logical_and(pi2_jk == 1, np.abs(pi2_ij - pi2_ik) < threshold), pi3_min_viewpoint_mat)

    del pi2_ij
    del pi2_ik
    del pi2_jk
    del gamma_ij_viewpoint_mat
    del gamma_ik_viewpoint_mat
    del gamma_jk_viewpoint_mat

    return beta3_viewpoint_mat, gamma_viewpoint_mat, pi2_sum_viewpoint_mat, omega_viewpoint_mat, pi3_min_viewpoint_mat, pi3_max_viewpoint_mat


def compute_pi3_viewpoint_matrix(s3_viewpoint_mat, omega_viewpoint_mat, pi2_sum_viewpoint_mat,
                                 beta3_viewpoint_mat, gamma_viewpoint_mat):
    """ Computes the PI3 matrix, fixing a viewpoint i. """

    pi3_viewpoint_mat = (s3_viewpoint_mat - omega_viewpoint_mat + (pi2_sum_viewpoint_mat - 1) * beta3_viewpoint_mat) /\
                        (alpha ** 2 - gamma_viewpoint_mat + 2 * beta3_viewpoint_mat)

    return pi3_viewpoint_mat


def SLICE_TRIPLEWISE_VIEWPOINT(i, segregation_table, chr):
    """ Performs triplewise SLICE, fixing a viewpoint i. """

    if chr == "chrX" or chr == "chrY":
        K = int(effective_NPs_per_tube / 2)
        F = 1 - (1 - F_mean) ** (1 / 2)
    else:
        K = effective_NPs_per_tube
        F = F_mean

    F1_arr = compute_tube_segregation_frequency(segregation_table)
    F2_mat = compute_tube_cosegregation_matrix(segregation_table)
    F3_viewpoint_mat = compute_tube_co3segregation_viewpoint_matrix(i, segregation_table)

    pi2_mat, pi2_significant_mat, beta2_mat, beta2_arr = SLICE_PAIRWISE(segregation_table, chr, generator=None)

    pi2_mat[pi2_mat < 0] = np.nan
    pi2_mat[pi2_mat > 1] = 1.
    pi2_significant_mat[pi2_significant_mat > 1] = 1.

    eta_mat = pi2_mat * (alpha - beta2_mat) + beta2_mat

    s3_viewpoint_mat = compute_s3_viewpoint_matrix(i, F1_arr, F2_mat, F3_viewpoint_mat, eta_mat, K, F)


    beta3_viewpoint_mat, gamma_viewpoint_mat, pi2_sum_viewpoint_mat, omega_viewpoint_mat,\
    pi3_min_viewpoint_mat, pi3_max_viewpoint_mat = compute_viewpoint_matrices(i, pi2_mat, beta2_mat)

    pi3_viewpoint_mat = compute_pi3_viewpoint_matrix(s3_viewpoint_mat, omega_viewpoint_mat, pi2_sum_viewpoint_mat,
                                                     beta3_viewpoint_mat, gamma_viewpoint_mat)
    pi3_viewpoint_mat[np.isinf(pi3_viewpoint_mat)] = np.nan

    del F1_arr
    del F2_mat
    del F3_viewpoint_mat
    del s3_viewpoint_mat
    del pi2_mat
    del pi2_significant_mat
    del beta2_mat
    del beta2_arr
    del eta_mat
    del beta3_viewpoint_mat
    del gamma_viewpoint_mat
    del pi2_sum_viewpoint_mat
    del omega_viewpoint_mat

    return pi3_viewpoint_mat, pi3_min_viewpoint_mat, pi3_max_viewpoint_mat
