# SLICE 

This repository coontains the library developed by Federico Billeci, Luca Fiorillo and Francesco Musella to implement SLICE and interSLICE analysis on Genome Architecture Mapping data. Luca Fiorillo and Francesco Musella developed the scrips for the SLICE analysis of intra-chromosomal contacts. Federico Billeci developed the scripts for the interSLICE analysis of inter-chromosomal contacts.

## How to install the library in linux

To install the library download this repository and add its path to the enviroment variable PYTHONPATH:

export PYTHONPATH=$PYTHONPATH:/place_with_the_repository/slice_repository

## Data preparation for a new dataset

The following procedure has to be run only once, when a new segregation table is added to the repository.
1) GAM segregation table is provided a .table format. The .table file has to be renamed:
  * set the extension to .txt;
  * the chosen name preceeding the .txt extension, has to start with the word "rawdata". We call name_root this name;
  * the convention for name_root is: rawdata_identification-name_resolution_number-of-tubes. For instance a GAM segregation table at 500kb resolution relative to mESC cells, done with 300 tubes and 3 NP per tube will be called _rawdata_mesc_500kb_300x3.txt_;

2) Run the script _prepare_data.sh_.

3) Open _settings.py_:
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

4) in _slice_repository/starting_: run _create_segregation_pkl.py_; run _compute_F_mean.py_ and insert the printed output in the variable F_mean in _settings.py_. 

## Run SLICE 

Run SLICE and interSLICE on a given dataset.

In _slice_repository/settings.py_, uncomment only the line which assigns to the variable _name_root_ the name of the selected dataset.

### SLICE for intra-chromosomal pairwise contacts


### SLICE for intra-chromosomal threewise contacts


### interSLICE for inter-chromosomal pairwise contacts




## Analysis of SLICE matrices











