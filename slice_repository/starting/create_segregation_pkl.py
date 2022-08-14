import pandas as pd

#import sys
#sys.path.insert(1, '/home/federico/Università/Tesi_magistrale/SLICE/slice_repository/')


from src.utilities import *

# Creates a pandas.pkl file containing the segregation table for all chromosomes


def verify_line(input_line):
    """ Returns True if input_line is good, i.e. it contains data of one of the acceptable chromosomes. """

    input_list = input_line.split("\t")
    is_acceptable = input_list[0] in chr_dictionary.keys()

    return is_acceptable


def process_line(input_line):
    """ Process the input line, converting it into separate information.
        Returns: chromosome number, starting genomic position, ending genomic position, segregation data. """

    # Convert line into list of elements (each one is separated by the special character \t)"
    input_list = input_line.split("\t")

    # the first three elements are [n_chromosome, start, stop]
    chr = chr_dictionary[input_list[0]]
    start = int(input_list[1])
    stop = int(input_list[2])

    # segregation_data is [0, 0, 1, ..., 0], where each entry is an int
    segregation_data = list()
    for i in range(3, len(input_list)):
        segregation_data.append(int(input_list[i]))

    return chr, start, stop, segregation_data


print("Creating Segregation pkl file for " + name_root + "\n")

# Load raw data
in_file = open(data_path + name_root + "/rawdata_" + name_root + ".txt", "r")

# These lists contain the data for every line of the input file.
# Which chromosome list
chr_list = list()
# Starting genomic position list
start_list = list()
# Ending genomic position list
stop_list = list()
# Segregation table list
segregation_list = list()

print("     reading the input txt file...\n")
# Loop on the input file, for every row (corresponding to a DNA window)
for line in in_file:
    if verify_line(line):
        chr, start, stop, segregation_data = process_line(line)
        chr_list.append(chr)
        start_list.append(start)
        stop_list.append(stop)
        segregation_list.append(segregation_data)
in_file.close()

# Convert data into numpy array
chr_arr = np.array(chr_list)
start_arr = np.array(start_list)
stop_arr = np.array(stop_list)
segregation_table = np.array(segregation_list)

# The final data is set as a Dictionary.
#   Every key of the dictionary corresponds to a chromosome (as in acceptable chromosomes),
#   and to each chromosome corresponds another dictionary, containing start, stop and segregation data.
#   For instance, data["chr2"] is Dict( { "start_position": np.array([0, 50000, ...]),
#                                        "stop_position": np.array([50000, 100000, ...]),
#                                        "segregation_table": np.array[[0, 1, ...], [0, 0, ...], ...]) } )
print("     creating the pkl file...\n")
data = dict()
for n in chr_dictionary.values():
    chr = get_key(chr_dictionary, n)
    print("         evaluating " + chr + "...")
    data[chr] = dict({"start_position": start_arr[chr_arr == n],
                      "stop_position": stop_arr[chr_arr == n],
                      "segregation_table": segregation_table[chr_arr == n]})

# The data is finally stored as a pandas.DataFrame
df = pd.DataFrame(data)
df.to_pickle(data_path + name_root + "/segregation_" + name_root + ".pkl")

print("\nDone")
