import numpy as np
import os

# Select which dataset to use

#name_root = "iza-mesc_1Mb_420x3"
#name_root = "iza-mesc_150kb_420x3"
name_root = "mesc_46C_1Mb_481x1"
#name_root = "mesc_46C_150kb_481x1"
#name_root = "dopa30_150kb_482x3"


# Insert the path for the pickle containing the segregation table

data_path = os.path.dirname(os.path.abspath(__file__)) + '/data/'


def compute_alpha():

    r_win = r_cell * (resolution / genome_length)**(1/3)
    alpha = (h + 2 * r_win + 2 * r_cell) / (h + 2 * r_win)

    return alpha


human_chr_dictionary = dict({"chr1": 1, "chr2": 2, "chr3": 3, "chr4": 4,
                             "chr5": 5, "chr6": 6, "chr7": 7, "chr8": 8,
                             "chr9": 9, "chr10": 10, "chr11": 11, "chr12": 12,
                             "chr13": 13, "chr14": 14, "chr15": 15, "chr16": 16,
                             "chr17": 17, "chr18": 18, "chr19": 19, "chr20": 20,
                             "chr21": 21, "chr22": 22, "chrX": -1})
mouse_chr_dictionary = dict({"chr1": 1, "chr2": 2, "chr3": 3, "chr4": 4,
                             "chr5": 5, "chr6": 6, "chr7": 7, "chr8": 8,
                             "chr9": 9, "chr10": 10, "chr11": 11, "chr12": 12,
                             "chr13": 13, "chr14": 14, "chr15": 15, "chr16": 16,
                             "chr17": 17, "chr18": 18, "chr19": 19, "chrX": -1})

# Set the parameters for each dataset
if name_root == "iza-mesc_1Mb_420x3":
    effective_NPs_per_tube = 6
    resolution = 1000 * 10**3  # in bp
    r_cell = 4.6  # in um
    h = 0.22  # in um
    genome_length = 5.20 * 10**9  # in bp
    alpha = compute_alpha()
    v = 1 / alpha
    F_mean = 0.4876928549197457
    #F_mean = 0.47544819777316477
    chr_dictionary = mouse_chr_dictionary

if name_root == "iza-mesc_150kb_420x3":
    effective_NPs_per_tube = 6
    resolution = 150 * 10**3  # in bp
    r_cell = 4.6  # in um
    h = 0.22  # in um
    genome_length = 5.20 * 10**9  # in bp
    alpha = compute_alpha()
    v = 1 / alpha
    F_mean = 0.29287174039978625
    
    chr_dictionary = mouse_chr_dictionary

if name_root == "mesc_46C_1Mb_481x1":
    effective_NPs_per_tube = 2
    resolution = 1000 * 10**3  # in bp
    r_cell = 4.6  # in um
    h = 0.22  # in um
    genome_length = 5.20 * 10**9  # in bp
    alpha = compute_alpha()
    v = 1 / alpha
    F_mean = 0.23529319906226276
    #F_mean = 0.22815738408818795
    chr_dictionary = mouse_chr_dictionary

if name_root == "mesc_46C_150kb_481x1":
    effective_NPs_per_tube = 2
    resolution = 150 * 10**3  # in bp
    r_cell = 4.6  # in um
    h = 0.22  # in um
    genome_length = 5.20 * 10**9  # in bp
    alpha = compute_alpha()
    v = 1 / alpha
    F_mean = 0.12871518750761243
    chr_dictionary = mouse_chr_dictionary

if name_root == "dopa30_150kb_482x3":
    effective_NPs_per_tube = 6
    resolution = 150 * 10**3  # in bp
    r_cell = 3.25  # in um
    h = 0.22  # in um
    genome_length = 5.20 * 10**9  # in bp
    alpha = compute_alpha()
    v = 1 / alpha
    F_mean = 0.3160606254627158
    chr_dictionary = mouse_chr_dictionary

if name_root == "xxx":
    effective_NPs_per_tube = np.nan
    resolution = np.nan  # in bp
    r_cell = np.nan  # in um
    h = np.nan  # in um
    genome_length = np.nan  # in bp
    alpha = compute_alpha()
    v = 1 / alpha
    F_mean = np.nan
    chr_dictionary = mouse_chr_dictionary



eff_mean = (1 - (1 - F_mean) ** (1 / effective_NPs_per_tube)) / v

# List of chromosomes of interest, keys from the previous dictionary

chr_list = list(chr_dictionary.keys())

# Print info

def print_info():
    print("Parameters:")
    print("name_root: " + name_root)
    print("Effective number of NPs per tube: " + str(effective_NPs_per_tube))
    print("Resolution: " + str(resolution) + " bp")
    print("Cell nucleus radius: " + str(r_cell) + " um")
    print("Slice thickness h: " + str(h) + " um")
    print("Genome length: " + str(genome_length) + " bp")
    print("alpha = " + str(alpha))
    print("v = " + str(v))
    print("F_mean = " + str(F_mean))
    print("eff_mean = " + str(eff_mean))
    print("\n\n")
