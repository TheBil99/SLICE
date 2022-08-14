from matplotlib import pyplot as plt
from settings import *

# Contains functions needed through the codes

# Functions for segregation/cosegregation/NPMI matrices

def compute_tube_segregation_frequency(segregation_table):
    """ Computes the tube segregation frequency of each DNA window from the Segregation Table. """

    n_tubes = np.shape(segregation_table)[1]

    F_arr = np.sum(segregation_table, axis=1) / n_tubes
    F_arr[F_arr == 0] = np.nan

    return F_arr

def compute_tube_cosegregation_matrix(segregation_table):
    """ Computes the tube Co-Segregation matrix from the Segregation table. """

    n_tubes = np.shape(segregation_table)[1]

    F_mat = np.matmul(segregation_table, np.transpose(segregation_table)) / n_tubes
    np.fill_diagonal(F_mat, np.nan)

    return F_mat

def compute_npmi(segregation_table):
    """ Computes the NPMI normalized tube Co-Segregation matrix. """

    F_arr = compute_tube_segregation_frequency(segregation_table)
    F_mat = compute_tube_cosegregation_matrix(segregation_table)

    F_i = np.tile(F_arr, (len(F_arr), 1))
    F_j = np.transpose(F_i)
    F_ij = F_mat

    npmi_mat = - np.log(F_ij / (F_i * F_j)) / np.log(F_ij)
    np.fill_diagonal(npmi_mat, np.nan)

    return npmi_mat

def compute_npmi_inter(segregation_table_A, segregation_table_B):
    """ Computes the NPMI normalized tube Co-Segregation matrix. """
    F_arr_A = compute_tube_segregation_frequency(segregation_table_A)
    F_arr_B = compute_tube_segregation_frequency(segregation_table_B)

    F_mat = compute_tube_cosegregation_matrix_offdiag(segregation_table_A, segregation_table_B)

    F_i = np.tile(F_arr_A, (len(F_arr_B), 1)).T
    F_j = np.tile(F_arr_B, (len(F_arr_A), 1))
    F_ij = F_mat

    npmi_mat = - np.log(F_ij / (F_i * F_j)) / np.log(F_ij)

    return npmi_mat

def compute_tube_segregation_frequency_offdiag(segregation_table_A, segregation_table_B):
    """ Computes the tube segregation frequency of each DNA window from the Segregation Table. """

    n_tubes = np.shape(segregation_table_A)[1]

    F_arr_A = np.sum(segregation_table_A, axis=1) / n_tubes
    F_arr_A = np.sum(segregation_table_A, axis=1) / n_tubes
    F_arr_A[F_arr_A == 0] = np.nan

    F_arr_B = np.sum(segregation_table_B, axis=1) / n_tubes
    F_arr_B = np.sum(segregation_table_B, axis=1) / n_tubes
    F_arr_B[F_arr_B == 0] = np.nan


    return F_arr_A, F_arr_B

def compute_tube_cosegregation_matrix_offdiag(segregation_table_A, segregation_table_B):
    """ Computes the tube Co-Segregation matrix from the Segregation table. """

    n_tubes = np.shape(segregation_table_A)[1]

    F_mat = np.matmul(segregation_table_A, np.transpose(segregation_table_B)) / n_tubes

    #mettere che nei nan dei F_arr_A e F_arr_B la coseg è nan

    return F_mat

def compute_pi_nan(pi, sign_pi, seg_table, verbose = False):

    """ Function that returns the number of NaN, pi<0 and significative pi for a given pi matrix """
    """ The percentage is carried out only on the upper triangular matrix: this is in order to obtain """
    """ the same values as the one obtained from the function single_chromosome in th script slice_pairwise.py  """

    coseg_zeros = compute_tube_cosegregation_matrix(seg_table)
    F_arr_A= compute_tube_segregation_frequency(seg_table)

    coseg_zeros[:, np.isnan(F_arr_A)] = np.nan
    coseg_zeros[np.isnan(F_arr_A), :] = np.nan
    
    a = ((np.count_nonzero(np.isnan(pi)) - pi.shape[0]   )) /((pi.shape[0]) * ((pi.shape[0]) - 1))
    b = 1 - ((np.count_nonzero(np.isnan(sign_pi)) - sign_pi.shape[0]   )) /((sign_pi.shape[0]) * ((sign_pi.shape[0]) - 1))
    c = (np.count_nonzero(np.isnan(coseg_zeros))  -coseg_zeros.shape[0]  )/((coseg_zeros.shape[0]) * ((coseg_zeros.shape[0]) - 1))
    d = pi[pi == 1].size/((pi.shape[0]) * ((pi.shape[0]) - 1))
    if(verbose == True):
        print("number of GAM nan:\t", c, "\nnumber of pi<0:\t", a-c, "\nnumber of significative pi:\t", b, "\nnumber of pi=1:\t", d)

    return c, a - c, b, d

def compute_pi_nan_inter(pi, sign_pi, seg_tableA, seg_table_B, verbose = False):
    """ Function that returns the number of NaN, pi<0 and significative pi for a given pi matrix in inter """

    coseg_zeros = compute_tube_cosegregation_matrix_offdiag(seg_tableA, seg_table_B)
    F_arr_A, F_arr_B= compute_tube_segregation_frequency_offdiag(seg_tableA, seg_table_B)

    coseg_zeros[:, np.isnan(F_arr_B)] = np.nan
    coseg_zeros[np.isnan(F_arr_A), :] = np.nan
    
    a = np.count_nonzero(np.isnan(pi)) /pi.size
    b = 1 - np.count_nonzero(np.isnan(sign_pi)) /sign_pi.size
    c = np.count_nonzero(np.isnan(coseg_zeros))/coseg_zeros.size
    d = pi[pi == 1].size/pi.shape[0]

    if(verbose == True):
        print("number of GAM nan:\t", c, "\nnumber of pi<0:\t", a-c, "\nnumber of significative pi:\t", b, "\nnumber of pi=1:\t", d)

    return c, a - c, b, d



# Other functions

def kth_diag_indices(mat, k):
    """ Returns the indices of the k-th upper diagonal of a matrix.
        Attention: k=1 is associated to the first super-diagonal. """

    rows, columns = np.diag_indices_from(mat)

    return rows[:-k], columns[k:]


def decide_starting_bin(start):
    """ Decides which bin is the starting one for the locus under consideration. """

    # The bin containing start is found using the floor function, since the bin counting starts from 0. For example:
    #   start = 0.01 Mb, resolution = 0.05 Mb --> bin = floor(0.01 / 0.05) = floor(0.2) = 0.
    bin_containing_start = np.floor(start / resolution)

    # The starting and ending positions of this bin are called head and tail.
    head = bin_containing_start * resolution
    tail = bin_containing_start * resolution + resolution

    # If start is closer to head, the starting bin is the one starting with head.
    # Otherwise, the starting bin is the one starting with tail.
    if np.abs(start - head) <= np.abs(start - tail):
        # The case where start is exactly in the middle of the bin is arbitrarily contained here.
        return int(bin_containing_start)
    if np.abs(start - head) > np.abs(start - tail):
        return int(bin_containing_start) + 1


def locus_bins(start, end):
    """ Computes the bin positions of a particular locus. """

    print("\n")
    print("Computing the locus bins:")
    print("start: " + str(start) + ",   end: " + str(end))

    # To decide where to start, whether from bin i or i+1, we use the function decide_starting_bin.
    start_bin = decide_starting_bin(start)
    # The end bin is found summing to start_bin the length of the region in unit of bins.
    # end_bin is NOT included.
    end_bin = start_bin + int((end - start) / resolution)
    print("start_bin: " + str(start_bin) + ",   end_bin: " + str(end_bin))
    print("\n")

    return start_bin, end_bin


def isolate_locus(x, start, end):
    """ Isolates the locus of interest from the input n-dimensional array x. """

    start_bin, end_bin = locus_bins(start, end)

    # 1D array (e.g. eff_array)
    if len(x.shape) == 1:
        return start_bin, end_bin, x[start_bin:end_bin]

    # 2D array
    if len(x.shape) == 2:
        # Square matrix (e.g. co-segregation matrix)
        if np.shape(x)[0] == np.shape(x)[1]:
            return start_bin, end_bin, x[start_bin:end_bin, start_bin:end_bin]
        # Not-square matrix, i.e. segregation table
        else:
            return start_bin, end_bin, x[start_bin:end_bin, :]


def get_key(dictionary, value):
    """ Gets the key entry of a dictionary for the input value. """

    for k in dictionary.keys():
        if value == dictionary[k]:
            return k


def print_data_info(x):
    """ Given an input n-dimensional array x, gives basic data info,
        i.e. shape, mean, std, percentage of NaN values, <0 values and >1 values. """

    print("shape: " + str(x.shape))
    print("mean: " + str(np.nanmean(x)))
    print("std: " + str(np.nanstd(x)))

    x_flatten = x.flatten()
    x_flatten_not_nan = x_flatten[~np.isnan(x_flatten)]

    if len(x_flatten_not_nan) == 0:
        print("Only NaN values found!")
        print("\n\n")
        del x_flatten
        del x_flatten_not_nan
        return None

    print("% of NaN: " + str(100 - 100 * len(x_flatten_not_nan) / len(x_flatten)) + " %")
    print("% of = 0: " + str(100 * len(x_flatten[x_flatten == 0]) / len(x_flatten)) + " %")
    print("% of = 1: " + str(100 * len(x_flatten[x_flatten == 1]) / len(x_flatten)) + " %")
    print("% of < 0: " + str(100 * len(x_flatten[x_flatten < 0]) / len(x_flatten)) + " %")
    print("% of > 1: " + str(100 * len(x_flatten[x_flatten > 1]) / len(x_flatten)) + " %")

    print("\n\n")

    del x_flatten
    del x_flatten_not_nan

    return None


def remove_symmetric(pi3_arr, i_arr, j_arr, k_arr):
    """ Given as input four n-dimensional arrays pi3_arr, i_arr, j_arr, k_arr,
        it keeps only the entries where i>j and j>k (and thus i>k). """

    pi3_copy_arr = np.copy(pi3_arr)
    i_copy_arr = np.copy(i_arr)
    j_copy_arr = np.copy(j_arr)
    k_copy_arr = np.copy(k_arr)

    # Keep only entries where i > j
    pi3_arr = pi3_copy_arr[i_copy_arr > j_copy_arr]
    i_arr = i_copy_arr[i_copy_arr > j_copy_arr]
    j_arr = j_copy_arr[i_copy_arr > j_copy_arr]
    k_arr = k_copy_arr[i_copy_arr > j_copy_arr]
    pi3_copy_arr = np.copy(pi3_arr)
    i_copy_arr = np.copy(i_arr)
    j_copy_arr = np.copy(j_arr)
    k_copy_arr = np.copy(k_arr)
    # Keep only entries where j > k
    pi3_arr = pi3_copy_arr[j_copy_arr > k_copy_arr]
    i_arr = i_copy_arr[j_copy_arr > k_copy_arr]
    j_arr = j_copy_arr[j_copy_arr > k_copy_arr]
    k_arr = k_copy_arr[j_copy_arr > k_copy_arr]

    del pi3_copy_arr
    del i_copy_arr
    del j_copy_arr
    del k_copy_arr

    return pi3_arr, i_arr, j_arr, k_arr


def compress_symmetric_tensor(tensor):
    """ Given an input n*n*n fully symmetric tensor, produces a compressed version of the tensor.

        The output is a matrix composed by four 1-dimensional arrays stacked horizontally,
        where the first row are the PI3s, the second, third and fourth rows are i, j, k:
                0.23    0.11    0.04    ...
                1       14      44      ...
                4       28      51      ...
                9       33      89      ...
        which means that PI3[1, 4, 9] = 0.23, PI3[14, 28, 33] = 0.11, and so on.
        Notice that NaNs and repetitions are removed in the compressed output. """

    tensor_nonan_arr = tensor[~np.isnan(tensor)]
    i_nonan_arr, j_nonan_arr, k_nonan_arr = np.where(~np.isnan(tensor))

    tensor_nonan_arr, i_nonan_arr, j_nonan_arr, k_nonan_arr = remove_symmetric(tensor_nonan_arr,
                                                                               i_nonan_arr,
                                                                               j_nonan_arr,
                                                                               k_nonan_arr)

    compressed_tensor = np.vstack((tensor_nonan_arr, i_nonan_arr, j_nonan_arr, k_nonan_arr))

    del tensor_nonan_arr
    del i_nonan_arr
    del j_nonan_arr
    del k_nonan_arr

    return compressed_tensor


def uncompress_symmetric_tensor(n, compressed_tensor):
    """ Decompress a compressed fully symmetric tensor (described in the function above) into a n*n*n tensor. """

    tensor = np.nan * np.ones((n, n, n), dtype=np.float32)

    tensor_nonan_arr = compressed_tensor[0, :]
    i_nonan_arr = compressed_tensor[1, :].astype(int)
    j_nonan_arr = compressed_tensor[2, :].astype(int)
    k_nonan_arr = compressed_tensor[3, :].astype(int)

    for tensor_ijk, i, j, k in zip(tensor_nonan_arr, i_nonan_arr, j_nonan_arr, k_nonan_arr):
        tensor[i, j, k] = tensor_ijk
        tensor[i, k, j] = tensor_ijk
        tensor[j, i, k] = tensor_ijk
        tensor[k, i, j] = tensor_ijk
        tensor[j, k, i] = tensor_ijk
        tensor[k, j, i] = tensor_ijk

    del tensor_nonan_arr
    del i_nonan_arr
    del j_nonan_arr
    del k_nonan_arr

    return tensor


def convert_window_arr_to_string(i_arr, chr):
    """ Converts the array containing window bin indices into chr:start-end string format.

        Example: i_arr = [3, 5, 7], resolution = 5 bp, chr = chr1
                 3 --> start = 3 * 5 = 15 bp, end = start + resolution = 20 bp,
                 5 --> start = 5 * 5 = 25 bp, end = start + resolution = 30 bp,
                 7 --> start = 7 * 5 = 35 bp, end = start + resolution = 40 bp,
                 out_str = ['chr1:15-20', 'chr1:25-30', 'chr1:35-40']. """

    start_arr = i_arr * resolution
    end_arr = start_arr + resolution

    start_arr_str = start_arr.astype(str)
    end_arr_str = end_arr.astype(str)

    out_str = np.char.add(chr + ":", start_arr_str)
    out_str = np.char.add(out_str, "-")
    out_str = np.char.add(out_str, end_arr_str)

    return out_str


def compute_top2_to_store(compressed_tensor, chr):
    """ Given an input compressed tensor, computes the top2% tensor values.
        Data are returned as horizontal arrays:
            val1 val2 val3 ...
            i1   i2   i3  ...
            j1   j2   j3  ...
            k1   k2   k3  ...,
        which means that tensor[i1,j1,k1]=val1, tensor3[i2,j2,k3]=pi2, ...
        The coordinates i, j, k in str format and contain
        chromosome and positions in bp (see convert_window_arr_to_string function). """

    tensor_arr = compressed_tensor[0, :]
    i_arr = compressed_tensor[1, :].astype(int)
    j_arr = compressed_tensor[2, :].astype(int)
    k_arr = compressed_tensor[3, :].astype(int)

    threshold = np.nanpercentile(tensor_arr, 98)
    print("Top 2% threshold: " + str(threshold))

    tensor_top2_arr = tensor_arr[tensor_arr >= threshold]
    i_top2_arr = i_arr[tensor_arr >= threshold]
    j_top2_arr = j_arr[tensor_arr >= threshold]
    k_top2_arr = k_arr[tensor_arr >= threshold]

    i_top2_str_arr = convert_window_arr_to_string(i_top2_arr, chr)
    j_top2_str_arr = convert_window_arr_to_string(j_top2_arr, chr)
    k_top2_str_arr = convert_window_arr_to_string(k_top2_arr, chr)

    top2_data_to_save = np.vstack((tensor_top2_arr, i_top2_str_arr, j_top2_str_arr, k_top2_str_arr))

    del tensor_arr
    del i_arr
    del j_arr
    del k_arr
    del tensor_top2_arr
    del i_top2_arr
    del j_top2_arr
    del k_top2_arr

    return top2_data_to_save


def compute_tensor_average_over_genomic_distances(tensor):
    """ Computes the average value of the input tensor for all genomic distances between the three windows.
        In particular, if n is the number of windows, it creates a (n-1)*(n-1) matrix
                genomic_mat = [[g(1,1), g(1,2), ..., g(1,n-1)],
                               [g(2,1), g(2,2), ..., g(2,n-1)],
                                        ...
                               [g(n,1), g(n,2), ..., g(n,n)]]
        where g(a,b) is the average tensor values between all triplets i, j, k with j=i-a and k=i+b:
                g(a,b) = average over all possible i of tensor[i, j=i-a, k=i+b].
        The matrix is returned as a horizontal flatten array. """

    genomic_mat = np.zeros((len(tensor) - 1, len(tensor) - 1))

    # Here a and b range from 1 to n-1,
    # corresponding to the minimum and the maximum genomic distances to be considered
    for a in range(1, len(tensor)):
        for b in range(1, len(tensor)):
            # We store in this list all the tensor values with j=i-a and k=i+b and varying i
            values_at_distance_ab = list()
            for i in range(len(tensor)):
                if i - a in range(len(tensor)) and i + b in range(len(tensor)):
                    values_at_distance_ab.append(tensor[i, i - a, i + b])
            # Convert list into array
            values_at_distance_ab = np.array(values_at_distance_ab)
            # If values_at_distance_ab has only NaN, we set the corresponding entry of genomic_mat as NaN
            if len(values_at_distance_ab[~np.isnan(values_at_distance_ab)]) == 0:
                # Here we have to subtract 1 to a and b,
                # because the entries of a (n-1)*(n-1) numpy matrix go from 0 to n-2
                genomic_mat[a - 1, b - 1] = np.nan
            # If the list is not made up only of NaNs we set the entry as the mean of the values
            else:
                genomic_mat[a - 1, b - 1] = np.nanmean(values_at_distance_ab)
            del values_at_distance_ab

    genomic_mat = genomic_mat.flatten()
    genomic_mat = np.reshape(genomic_mat, (1, len(genomic_mat)))

    return genomic_mat


def compute_zeros_genomic_distance(tensor):
    """ With the same procedure as in 'compute_triplewise_genomic_distance_matrix',
        computes the percentage of 0 values at all fixed genomic distances of the input tensor. """

    zeros_mat = np.zeros((len(tensor) - 1, len(tensor) - 1))

    # Here a and b range from 1 to n-1,
    # corresponding to the minimum and the maximum genomic distances to be considered
    for a in range(1, len(tensor)):
        for b in range(1, len(tensor)):
            count_a_b = 0  # counts the number of entries at fixed genomic distance (a, b)
            for i in range(len(tensor)):
                if i - a not in range(len(tensor)) or i + b not in range(len(tensor)):
                    continue
                if np.isnan(tensor[i, i - a, i + b]):
                    continue
                count_a_b = count_a_b + 1
                if tensor[i, i - a, i + b] == 0:
                    # Here we have to subtract 1 to a and b,
                    # because the entries of a (n-1)*(n-1) numpy matrix go from 0 to n-2
                    zeros_mat[a - 1, b - 1] = zeros_mat[a - 1, b - 1] + 1
            if count_a_b == 0:
                zeros_mat[a - 1, b - 1] = np.nan
            else:
                zeros_mat[a - 1, b - 1] = zeros_mat[a - 1, b - 1] / count_a_b

    zeros_mat = zeros_mat.flatten()
    zeros_mat = np.reshape(zeros_mat, (1, len(zeros_mat)))

    return zeros_mat


def compute_histogram(x, bins, x_min, x_max):
    """ Computes the histogram data for the input n-dimensional array x.
        It produces a value array and a range array, e.g.
                value   range
                1.4     -1.5
                0.7     -1.4
               ...       ...
        which means that the histogram has a value 1.4 from -1.5 to -1.4, 0.7 from -1.4 to -1.3 and so on.
        Note that the len(value array) = bins and len(range array) = bins + 1.
        Data are returned in a unique horizontal array [value array, range array]. """

    x_nonan_arr = x[~np.isnan(x)]

    hist_data = plt.hist(x_nonan_arr, bins=bins, range=(x_min, x_max), histtype='step')

    hist_data = np.hstack((hist_data[0], hist_data[1]))
    hist_data = np.reshape(hist_data, (1, len(hist_data)))

    del x_nonan_arr

    return hist_data


def compute_percentage_of_data_between_values(x, val1, val2):
    """ Computes the percentage (out of not-nan) of entries of x (numpy array or tensor) between two input values. """

    x_flatten = x.flatten()
    x_flatten_not_nan = x_flatten[~np.isnan(x_flatten)]

    if len(x_flatten_not_nan) == 0:
        del x_flatten
        del x_flatten_not_nan
        return np.nan

    x_flatten_not_nan_between_values = np.copy(x_flatten_not_nan)
    x_flatten_not_nan_between_values = x_flatten_not_nan_between_values[x_flatten_not_nan_between_values > val1]
    x_flatten_not_nan_between_values = x_flatten_not_nan_between_values[x_flatten_not_nan_between_values < val2]

    percentage_between_values = len(x_flatten_not_nan_between_values) / len(x_flatten_not_nan)

    del x_flatten
    del x_flatten_not_nan
    del x_flatten_not_nan_between_values

    return percentage_between_values


def print_problematic_tensor(x):
    """ Print the problematic values as they are defined in the slice_triplewise function. """

    x_flatten = x.flatten()

    print("% of triplets from NaN data: " + str(100 * len(x_flatten[x_flatten == 0]) / len(x_flatten)) + " %")
    print("% of triplets from PI2<0: " + str(100 * len(x_flatten[x_flatten == -1]) / len(x_flatten)) + " %")
    print("% of triplets from PI2>0: " + str(100 * len(x_flatten[x_flatten == 1]) / len(x_flatten)) + " %")

    print("\n\n")

    return None
