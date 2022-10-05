import numpy as np
import gzip
import re
from scipy.spatial.distance import squareform
import pickle

#Library containing some functions to extract info from usual 
#cosegregation matrices in text (or gz) format

def count_lines(path, file_name, Zip = True):
    #conto il numero di righe nel file
    if(Zip):
        with gzip.open(path + file_name,'rt') as f:
            for i, line in enumerate(f):
                pass
            return i+1
    else:
        with open(path + file_name,'rt') as f:
            for i, line in enumerate(f):
                pass 
            return i+1

def check_file(path, file_name, Zip = True):                                 #check se il file effettivamente contiene dei numeri 
    count = 0
    if(Zip):
        with gzip.open(path + file_name,'rt') as f:
            for i, line in enumerate(f):
                numbers = re.findall('[0-9]+', line)
                if(len(numbers )> 3):
                    count += 1
    else:
        with open(path + file_name,'rt') as f:
            for i, line in enumerate(f):
                numbers = re.findall('[0-9]+', line)
                if(len(numbers )> 3):
                    count += 1
        return count

def get_line(path, file_name, line_n, Zip = True):                                   #restituisco l'i-esima riga
    
    if(Zip):
        with gzip.open(path + file_name,'rt') as f:
            for i, line in enumerate(f):
                if(i == line_n-1):
                    selected_line = line
    else:
        with open(path + file_name,'rt') as f:
            for i, line in enumerate(f):
                if(i == line_n-1):
                    selected_line = line
        
    
    
    return selected_line

def select_locus(path, file_name, start, stop, resolution, skip_rows,Zip = False):
    #start e stop contengono la prima e l'ultima posizione (in bp's) del locus a cui sono interessato
    #prima inclusa, ultima esclusa
    inf_idx = round(start/resolution) + skip_rows
    if(stop!=-1):
        sup_idx = round(stop/resolution) + skip_rows   #passo a posizioni espresse in windows
    else:
        sup_idx = count_lines(path, file_name, Zip)
    l = []
    if(Zip):
        with gzip.open(path + file_name,'rt') as f:
            for i, line in enumerate(f):
                
                if(i>= sup_idx):						#scritto male ma se lo modifico mi da errore di indentazione, si potrebbero riassorbire questi due if
                    break
                if(i>= inf_idx):
                    l.append(line)
                
        return l

    else:
        with open(path + file_name,'rt') as f:
            for i, line in enumerate(f):
                
                if(i>= sup_idx):
                    break
                if(i>= inf_idx):
                    l.append(line)
                
        return l

def convert(a):
    try: 
        return float(a)
    except:
        return np.nan

def square_selected_locus(path, file_name, start, stop, resolution, skip_rows ,skip_chars ,separator = "\t", Zip = False): 
    #funzione raramente utilizzata, il più delle volte si può rimpiazzare con np.loadtxt, ogni tanto però loadtxt
    #può essere scomoda se ci sono colonne da saltare o in generale il parsing è un poò più antipatico
    #in tal caso usare questa

    #skip_chars è il numero di caratteri non informativi presenti all'inizio di ogni riga                                       
    #rendo il selected_locus una matrice triangolare superiore
    
    sel_locus = select_locus(path, file_name, start, stop, resolution, skip_rows, Zip)
    
    inf_idx = round(start/resolution);    
    
    if(stop!=-1):
        sup_idx = round(stop/resolution)    #passo a posizioni espresse in windows
    else:
        sup_idx = count_lines(path, file_name, Zip)
    
    square_selected_locus = [];    count = 0
    for i in sel_locus:
        temp = i.split(separator)[count+ skip_chars + inf_idx+1:sup_idx + skip_chars]			#il +1 perchè così escludo gli zeri della diagonale
        square_selected_locus += list(map(convert, temp))
        count += 1
    return squareform(np.array(square_selected_locus))

def txt_to_matrix(path, file_name, skip_cols = 0, skip_rows = 0, separator = " ", Zip = False ):
    l = []
    if(Zip):
        with gzip.open(path + file_name,'rt') as f:
            for line in f:
                l.append(line)    
        

    else:
        with open(path + file_name,'rt') as f:
            for line in f:  
                l.append(line)    
     
    del l[:skip_rows]   
       
    matrix = []
    for i in l:
        temp = i.split(separator)[skip_cols:-1]						#-1 perchè l'ultimo termine è il \n	
        matrix.append(list(map(convert, temp)))
        
        
    return np.array(matrix)
   
def txt_to_matrix_single_line(path, file_name, skip_cols = 0, skip_rows = 0, separator = " ", Zip = False ):
    l = ""
    if(Zip):
        with gzip.open(path + file_name,'rt') as f:
            for line in f:
                l+=line    
        

    else:
        with open(path + file_name,'rt') as f:
            for line in f:  
                l+=line
    temp = re.sub("\n", " ", l)
    temp = re.split(separator, temp)
    
    matrix = list(map(convert, temp))
        
        
    return np.array(matrix)
    
def Import(path, file_name, Type = "pkl"):
	if Type == "pkl":
		with open(path + file_name, 'rb') as f:
			mat = pickle.load(f)
		return mat
	#elif( Type == "txt"):
	#	with open(path + file_name, 'rt') as f:
	#		mat = pickle.load(f)
	#	return mat
	#elif( Type == "txt.gz"):
	#	with gzip.open(path + file_name, 'rt') as f:
	#		mat = pickle.load(f)
	#	return mat
	else:
		print("insert a valid format type\n")
		return

def save_pkl(mat, path, file_name):
	with open(path + file_name, 'wb') as f:
		pickle.dump(mat, f)
	














