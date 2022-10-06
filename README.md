# SLICE 

This repository coontains the library to implement SLICE and interSLICE analysis on Genome Architecture Mapping data ( Beagrie RA et Al. , Nature 543, 519–524 (2017) ).  
The library is fully developed on Python version 3.9.7 .

## Table of Contents
* [Dependencies](#dependencies)
* [Data preparation for a new dataset](#data-preparation-for-a-new-dataset)
* [Run SLICE](#run-slice)
* [Analysis of SLICE matrices](#analysis-of-slice-matrices)
* [Acknowledgements](#acknowledgements)
* [Contacts](#contacts)


## Dependencies

- numpy version 1.21.2
- matplotlib version 3.4.2
- pandas version 1.3.3

## How to install the library 

To install the library download this repository. It is necessary to add the path of the repository to the enviroment variable PYTHONPATH to allow the import of the library in Python.

### Linux and Mac

Add the following line to _.bashrc_ file:

_export PYTHONPATH="${PYTHONPATH}:/path_to_the_repository/slice_repository"_

### Windows

## Data preparation for a new dataset

The following procedure has to be run only once, when a new segregation table is added to the repository.
- GAM segregation table is provided a .table format. The .table file has to be renamed:
  * set the extension to .txt;
  * the chosen name preceeding the .txt extension, has to start with the word "rawdata". We call name_root this name;
  * the convention for name_root is: rawdata_identification-name_resolution_number-of-tubes. For instance a GAM segregation table at 500kb resolution relative to mESC cells, done with 300 tubes and 3 NP per tube will be called _rawdata_mesc_500kb_300x3.txt_;

- Run the script _prepare_data.sh_.

- Open _settings.py_:
  - add a line where the chosen name for name_root is assigned to the variable name_root. Following the example of point 1): _name_root = rawdata_mesc_500kb_300x3_
  - insert the experimental parameter relative to the new segregation table in the following format, substituting the _xxx_ with the actual values:  
  if name_root == "XXX":  
    effective_NPs_per_tube = xxx  
    resolution = xxx  
    r_cell = xxx  
    h = xxx   
    genome_length = xxx    
    alpha = compute_alpha()  
    v = 1 / alpha  
    F_mean = np.nan  
    chr_dictionary = mouse_chr_dictionary  
  
   - effective_NPs_per_tube is the number of NPs multiplied by 2 if the experiment is unphased, if the experiment is phased it is the number of NPs;
   - resolution is the segregation table resolution in bp (if it is 500kb, insert 500000);
   - r_cell is the cell radius in um;
   - h is the NP thickness 
   - genome_length is the lenght of the whole genome, counting twice the diploid chromosomes and including chromosomes X and Y;   
   - F_mean has to be left np.nan initially;
   - chr_dictionary is mouse_chr_dictionary or human_chr_dictionary if the segregation table is relative to mouse or human, respectively; 

- In _slice_repository/starting_: run _create_segregation_pkl.py_; run _compute_F_mean.py_ and insert the printed output in the variable F_mean in _settings.py_. 

## Run SLICE 

Run SLICE and interSLICE on a given dataset.

In _slice_repository/settings.py_, uncomment only the line which assigns to the variable _name_root_ the name of the selected dataset.

### SLICE for intra-chromosomal pairwise contacts
Run _slice_repository/main/compute_PI2.py_. The chromosome to run SLICE on and the significance threshold can be passed as arguments in the command line. For instance, to run SLICE pairwise on chromosome 10 and to produce significant pi matrices with 95% as significance threshold, run the following command:  
_python3 compute_PI2.py 10 95_

The output of this run will be a pi matrix, a significant pi matrix and the beta array. The matrices are saved in the folder _slice_repository/data/name_root/PI2_, beta array is saved in _slice_repository/data/name_root/beta_. Informations about the pi calculation are printed in the standard output.  
  
PI2 matrices are symmetric matrices with _NaN_ along the diagonal. Only the upper triangular matrix is saved as horizontal array.

### SLICE for intra-chromosomal threewise contacts

Run _slice_repository/main/compute_PI3.py_. The desired chromosome has to be indicated as argument of the function _single_chromosome_ in line l-72. PI3 are saved in _slice_repository/data/name_root/PI3_. Informations about the pi calculation are printed in the standard output.  
  
PI3 are completely symmetrical tensors, saved with a compression procedure for optimization purpose.

### SLICE for intra-chromosomal threewise contacts with a fixed viewpoint

Run _slice_repository/main/compute_PI3_viewpoint.py_. The desired chromosome and the position of the viewpoint in bp have to be indicated as arguments of the _function chromosome_fixed_i_ in line l-68. PI3 respect to the selected viewpoint are saved in _slice_repository/data/name_root/PI3_viewpoint_. Informations about the pi calculation are printed in the standard output.  

PI3 matrices with a fixed viewpoint are symmetric matrices with _NaN_ along the diagonal. Only the upper triangular matrix is saved as horizontal array.

### interSLICE for inter-chromosomal pairwise contacts

Run _slice_repository/main/compute_PI2_inter_beta_evaluation.py_. By default this scripts computes intra- and inter-chromosomal pi matrices on the whole genome, to restrict the computation only on specific chromosomes modify the lines l-31 to l-42. Intra-chromosomal PI2 are saved in _slice_repository/data/name_root/PI2_, inter-chromosomal PI2 are saved in _slice_repository/data/name_root/PI2_inter_beta_evaluation_.  

## Analysis of SLICE matrices

Inside the folder _slice_repository/data_analysis_, the script _contact_matrix_utils.py_ contains functions used to analyze the results produced by SLICE and interSLICE. 
The same folder contains the notebooks _chr12_chr18_1NP_1Mb_beta_comparison.ipynb_,  _chr12_chr18_1NP_150kb.ipynb_, and _averages_genomewide_1NP_1Mb.ipynb_. 

## Acknowledgements

Luca Fiorillo and Francesco Musella developed the scrips for the SLICE analysis of intra-chromosomal contacts. Federico Billeci developed the scripts for the interSLICE analysis of inter-chromosomal contacts.  

## Contacts

Feel free to contact me at _federicobilleci99@gmail.com_.
