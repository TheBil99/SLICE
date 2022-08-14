from src.utilities import *
from src.slice_pairwise import SLICE_PAIRWISE

# Functions for three-way SLICE

# We use the number 1, 2 or 3 to indicates whether the quantity has 1, 2 or 3 indices
# For instance, F2_mat is F_ij, i.e. the tube Co-Segregation frequency matrix.


def compute_tube_co3segregation_tensor(segregation_table):
    """ Computes the 3-way tube Co-3-Segregation tensor from the Segregation Table. """

    n_windows = np.shape(segregation_table)[0]
    n_tubes = np.shape(segregation_table)[1]

    # Computation of F_ijk
    # s_tile and s_tile_t are transposition in 3D-space of the segregation table
    s_tile = np.tile(segregation_table, (len(segregation_table), 1, 1))
    s_tile_t = np.transpose(s_tile, axes=[1, 0, 2])
    F3_tensor = np.matmul(s_tile * s_tile_t, np.transpose(segregation_table)) / n_tubes

    # Set the diagonals of F_ijk to NaN (where i=j, i=k, j=k)
    for i in range(n_windows):
        np.fill_diagonal(F3_tensor[i, :, :], np.nan)
        np.fill_diagonal(F3_tensor[:, i, :], np.nan)
        np.fill_diagonal(F3_tensor[:, :, i], np.nan)

    # Free memory
    del s_tile
    del s_tile_t

    return F3_tensor


def compute_s3_tensor(F1_arr, F2_mat, F3_tensor, eta_mat, K, F):
    """ Computes the SLICE-Normalized Co-3-Segregation tensor. """

    # Quantities related to F1_arr
    p = 1 - (1 - F) ** (1 / K)
    F_i = np.tile(F1_arr, (len(F1_arr), len(F1_arr), 1))
    F_j = np.transpose(F_i, axes=[2, 1, 0])
    F_k = np.transpose(F_i, axes=[1, 2, 0])

    # Quantities related to F2_mat
    # The order here is chosen to match the one in F_i, F_j, F_k for term_3.
    F_ik = np.tile(F2_mat, (len(F2_mat), 1, 1))
    F_jk = np.transpose(F_ik, axes=[1, 2, 0])
    F_ij = np.transpose(F_ik, axes=[2, 0, 1])

    # Quantity related to F3_tensor
    F_ijk = F3_tensor

    # Quantities related to eta_mat
    eta_ik = np.tile(eta_mat, (len(eta_mat), 1, 1))
    eta_jk = np.transpose(eta_ik, axes=[1, 2, 0])
    eta_ij = np.transpose(eta_ik, axes=[2, 0, 1])

    term_1 = 1 - 3 * p + (p ** 2) * (eta_ij + eta_ik + eta_jk)
    term_2 = - 1 + F_ijk / (F_ij + F_ik + F_jk - F_i - F_j - F_k)
    term_3 = 3 * (1 - p) ** K - (1 - 2 * p + (p ** 2) * eta_ij) ** K - (1 - 2 * p + (p ** 2) * eta_ik) ** K - (1 - 2 * p + (p ** 2) * eta_jk) ** K

    s3_tensor = (term_1 - (1 + term_2 * term_3) ** (1 / K)) / (p ** 3)

    # Free memory
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

    return s3_tensor


def compute_beta3_tensor(i, j, k, beta_ik, beta_jk, beta_ij, n):
    """ Computes the beta3 tensor. """

    beta3_tensor = np.nan * np.ones((n, n, n), dtype=np.float32)

    # The order of i, j, k, beta_ij, beta_ik, beta_jk here is chosen so that the appropriate indices are compared.
    # We verify that the order is correct in the testing section.

    # i is between j and k
    np.putmask(beta3_tensor, np.logical_and(j < i, i < k), beta_ij * beta_ik)
    np.putmask(beta3_tensor, np.logical_and(k < i, i < j), beta_ij * beta_ik)
    # j is between i and k
    np.putmask(beta3_tensor, np.logical_and(i < j, j < k), beta_ij * beta_jk)
    np.putmask(beta3_tensor, np.logical_and(k < j, j < i), beta_ij * beta_jk)
    # k is between i and j
    np.putmask(beta3_tensor, np.logical_and(i < k, k < j), beta_ik * beta_jk)
    np.putmask(beta3_tensor, np.logical_and(j < k, k < i), beta_ik * beta_jk)

    return beta3_tensor


def compute_gamma_12_tensor(k, j, i, beta_ik, beta_jk, n):
    """ Computes the gamma tensor where indices 1 and 2 are interacting. """

    # Here for simplicity we write gamma_12 as gamma_ij
    # Again, the order is chosen so that the appropriate indices are compared.
    # Notice that the order is different in the three gamma functions.
    gamma_ij = np.nan * np.ones((n, n, n), dtype=np.float32)
    np.putmask(gamma_ij, np.abs(i - k) <= np.abs(j - k), beta_ik / v)
    np.putmask(gamma_ij, np.abs(i - k) > np.abs(j - k), beta_jk / v)

    return gamma_ij


def compute_gamma_13_tensor(i, k, j, beta_ij, beta_jk, n):
    """ Computes the gamma tensor where indices 1 and 3 are interacting. """

    # Here for simplicity we write gamma_13 as gamma_ik
    gamma_ik = np.nan * np.ones((n, n, n), dtype=np.float32)
    np.putmask(gamma_ik, np.abs(i - j) <= np.abs(k - j), beta_ij / v)
    np.putmask(gamma_ik, np.abs(i - j) > np.abs(k - j), beta_jk / v)

    return gamma_ik


def compute_gamma_23_tensor(j, i, k, beta_ik, beta_ij, n):
    """ Computes the gamma tensor where indices 2 and 3 are interacting. """

    # Here for simplicity we write gamma_23 as gamma_jk
    gamma_jk = np.nan * np.ones((n, n, n), dtype=np.float32)
    np.putmask(gamma_jk, np.abs(j - i) <= np.abs(k - i), beta_ij / v)
    np.putmask(gamma_jk, np.abs(j - i) > np.abs(k - i), beta_ik / v)

    return gamma_jk


def compute_tensors(pi2_mat, beta2_mat):
    """ Computes beta3_tensor, gamma_tensor, pi3_sum_tensor, omega_tensor, pi3_min_tensor, pi3_max_tensor, where
            gamma_tensor is  as gamma_ij + gamma_ik + gamma_jk,
            pi_sum_tensor is pi_ij + pi_ik + pi_jk,
            omega_tensor is pi_ij * gamma_ij + pi_ik * gamma_ik + pi_jk * gamma_jk. """

    n = len(beta2_mat)

    # Computation of generic indices

    # Here the indices are called 1 2 3, without a fixed relation to i, j or k.
    # We use them as inputs to the functions for beta3 and gammas:
    # these functions then convert them into i, j or k so that the appropriate order is always respected.

    index_permutation_1 = np.tile(np.arange(n, dtype=np.float32), (n, n, 1))
    index_permutation_2 = np.transpose(index_permutation_1, axes=[2, 1, 0])
    index_permutation_3 = np.transpose(index_permutation_1, axes=[1, 2, 0])

    beta2_permutation_1 = np.tile(beta2_mat, (len(beta2_mat), 1, 1))
    beta2_permutation_2 = np.transpose(beta2_permutation_1, axes=[1, 2, 0])
    beta2_permutation_3 = np.transpose(beta2_permutation_1, axes=[2, 0, 1])

    # Computation of beta3_tensor, gamma_tensor, pi2_sum_tensor and omega_tensor

    # The order here is chosen so as to match the gamma tensors for the computation of omega_tensor.
    pi2_23_tensor = np.tile(pi2_mat, (n, 1, 1))
    pi2_12_tensor = np.transpose(pi2_23_tensor, axes=[1, 2, 0])
    pi2_13_tensor = np.transpose(pi2_23_tensor, axes=[2, 0, 1])

    beta3_tensor = compute_beta3_tensor(index_permutation_1, index_permutation_2, index_permutation_3,
                                        beta2_permutation_1, beta2_permutation_2, beta2_permutation_3, n)
    gamma_12_tensor = compute_gamma_12_tensor(index_permutation_1, index_permutation_2, index_permutation_3,
                                              beta2_permutation_1, beta2_permutation_3, n)
    gamma_13_tensor = compute_gamma_13_tensor(index_permutation_1, index_permutation_2, index_permutation_3,
                                              beta2_permutation_1, beta2_permutation_2, n)
    gamma_23_tensor = compute_gamma_23_tensor(index_permutation_1, index_permutation_2, index_permutation_3,
                                              beta2_permutation_2, beta2_permutation_3, n)
    gamma_tensor = gamma_12_tensor + gamma_13_tensor + gamma_23_tensor
    pi2_sum_tensor = pi2_12_tensor + pi2_13_tensor + pi2_23_tensor
    omega_tensor = pi2_12_tensor * gamma_12_tensor + pi2_13_tensor * gamma_13_tensor + pi2_23_tensor * gamma_23_tensor

    # Computation of pi3_min_tensor and pi3_max_tensor

    pi3_min_tensor = (pi2_12_tensor + pi2_13_tensor + pi2_23_tensor - 1) / 2
    pi3_min_tensor[pi3_min_tensor < 0] = 0

    pi3_max_tensor = np.copy(pi2_12_tensor)
    np.putmask(pi3_max_tensor, pi3_max_tensor > pi2_13_tensor, pi2_13_tensor)
    np.putmask(pi3_max_tensor, pi3_max_tensor > pi2_23_tensor, pi2_23_tensor)

    # These lines of code implement the condition that when one between pi2_ij pi2_ik and pi2_jk is 1,
    # if the difference between the other two pi2s is less than a threshold (set as 0.1),
    # the corresponding pi3_ijk is set as the average value of these two pi2s.
    # To do so, here we impose that when one of the three pi2s is 1 (e.g. pi2_ij = 1)
    # and the difference between the other two is less than the threshold (e.g. |pi2_ik - pi2_jk| < 0.1),
    # then pi3_min_ijk = pi3_max_ijk = (pi2_ik + pi2_jk) / 2.
    # Then, in the main code compute_pi_triplewise, we have the condition that when pi3_min = pi3_max,
    # we set pi3 = pi3_min = pi3_max.
    threshold = 0.1
    np.putmask(pi3_min_tensor, np.logical_and(pi2_12_tensor == 1, np.abs(pi2_13_tensor - pi2_23_tensor) < threshold), (pi2_13_tensor + pi2_23_tensor) / 2)
    np.putmask(pi3_max_tensor, np.logical_and(pi2_12_tensor == 1, np.abs(pi2_13_tensor - pi2_23_tensor) < threshold), pi3_min_tensor)
    np.putmask(pi3_min_tensor, np.logical_and(pi2_13_tensor == 1, np.abs(pi2_12_tensor - pi2_23_tensor) < threshold), (pi2_12_tensor + pi2_23_tensor) / 2)
    np.putmask(pi3_max_tensor, np.logical_and(pi2_13_tensor == 1, np.abs(pi2_12_tensor - pi2_23_tensor) < threshold), pi3_min_tensor)
    np.putmask(pi3_min_tensor, np.logical_and(pi2_23_tensor == 1, np.abs(pi2_12_tensor - pi2_13_tensor) < threshold), (pi2_12_tensor + pi2_13_tensor) / 2)
    np.putmask(pi3_max_tensor, np.logical_and(pi2_23_tensor == 1, np.abs(pi2_12_tensor - pi2_13_tensor) < threshold), pi3_min_tensor)

    # Free memory
    del index_permutation_1
    del index_permutation_2
    del index_permutation_3
    del beta2_permutation_1
    del beta2_permutation_2
    del beta2_permutation_3
    del gamma_12_tensor
    del gamma_13_tensor
    del gamma_23_tensor
    del pi2_12_tensor
    del pi2_13_tensor
    del pi2_23_tensor

    return beta3_tensor, gamma_tensor, pi2_sum_tensor, omega_tensor, pi3_min_tensor, pi3_max_tensor


def compute_pi3_tensor(s3_tensor, omega_tensor, pi2_sum_tensor, beta3_tensor, gamma_tensor):

    pi3_tensor = (s3_tensor - omega_tensor + (pi2_sum_tensor - 1) * beta3_tensor) / (alpha**2 - gamma_tensor + 2 * beta3_tensor)

    return pi3_tensor


def compute_problematic_tensor(F3_tensor, pi2_mat):
    """ Computes a tensor which gives information on the problematic PI2 values.
        In particular, the entry [i,j,k] is:
            0       if F[i,j,k] is NaN (GAM data are NaN),
            -1      if PI2[i,j] or PI2[i,k] or PI2[j,k] are negative,
            1       otherwise. """

    n = len(pi2_mat)

    problematic_tensor = np.ones((n, n, n), dtype=np.float32)

    pi2_23_tensor = np.tile(pi2_mat, (n, 1, 1))
    pi2_12_tensor = np.transpose(pi2_23_tensor, axes=[1, 2, 0])
    pi2_13_tensor = np.transpose(pi2_23_tensor, axes=[2, 0, 1])

    problematic_tensor[pi2_12_tensor < 0] = -1
    problematic_tensor[pi2_13_tensor < 0] = -1
    problematic_tensor[pi2_23_tensor < 0] = -1

    problematic_tensor[np.isnan(F3_tensor)] = 0

    del pi2_12_tensor
    del pi2_13_tensor
    del pi2_23_tensor

    return problematic_tensor


def SLICE_TRIPLEWISE(segregation_table, chr):

    if chr == "chrX" or chr == "chrY":
        K = int(effective_NPs_per_tube / 2)
        F = 1 - (1 - F_mean) ** (1 / 2)
    else:
        K = effective_NPs_per_tube
        F = F_mean

    F1_arr = np.float32(compute_tube_segregation_frequency(segregation_table))
    F2_mat = np.float32(compute_tube_cosegregation_matrix(segregation_table))
    F3_tensor = compute_tube_co3segregation_tensor(segregation_table)

    F2_mat[:, np.isnan(F1_arr)] = np.nan
    F2_mat[np.isnan(F1_arr), :] = np.nan
    F3_tensor[:, :, np.isnan(F1_arr)] = np.nan
    F3_tensor[:, np.isnan(F1_arr), :] = np.nan
    F3_tensor[np.isnan(F1_arr), :, :] = np.nan

    pi2_mat, pi2_significant_mat, beta2_mat, beta2_arr = SLICE_PAIRWISE(segregation_table, chr, generator=None)
    pi2_mat = np.float32(pi2_mat)
    beta2_mat = np.float32(beta2_mat)
    del pi2_significant_mat
    del beta2_arr

    problematic_tensor = compute_problematic_tensor(F3_tensor, pi2_mat)

    pi2_mat[pi2_mat < 0] = np.nan
    pi2_mat[pi2_mat > 1] = 1.

    eta_mat = pi2_mat * (alpha - beta2_mat) + beta2_mat

    s3_tensor = compute_s3_tensor(F1_arr, F2_mat, F3_tensor, eta_mat, K, F)

    beta3_tensor, gamma_tensor, pi2_sum_tensor, omega_tensor, pi3_min_tensor, pi3_max_tensor = compute_tensors(pi2_mat, beta2_mat)

    pi3_tensor = compute_pi3_tensor(s3_tensor, omega_tensor, pi2_sum_tensor, beta3_tensor, gamma_tensor)
    pi3_tensor[np.isinf(pi3_tensor)] = np.nan

    del F1_arr
    del F2_mat
    del F3_tensor
    del s3_tensor
    del eta_mat
    del beta3_tensor
    del gamma_tensor
    del pi2_sum_tensor
    del omega_tensor

    return pi3_tensor, pi3_min_tensor, pi3_max_tensor, problematic_tensor
