In questo README scrivo indicazioni aggiuntive da seguire per configurare la libreria che esegue SLICE.
In seguito bisognerà produrre un unico README per fare una guida univoca per configurare la lib, al momento
si da per scontato che tutti i passaggi indicati da francesco vengano eseguiti. 


1. Nel file settings.py modificare il path in cui bisogna andare a cercare il pickle della segregation table
francesco questo path lo indicava come ../data/blabla ma questo path non è comodo nel momento in cui 
si intende runnare le funzioni separatamente o anche da dei notebook esterni, per questo ho messo
la variabile globale data_path, messa nel file settings.py, in modo che venga vista da tutti gli script

2. Le matrici non le salvo in formato txt ma in formato .npy, è un formato binary ottimizzato per salvare 
array numpy ed è più accessibile, per il resto non cambia nulla e se si vuole comunque esportare un .txt basta
cambiare una riga di codice

3. Le matrici di cosegregazione vengono prodotte dalla funzione compute_tube_segregation_matrix in intra e da 
compute_tube_segregation_matrix_offdiag in inter. Queste matrici escono senza alcun valore nan. Questo è perchè
queste funzioni dialogano solamente con le segregation tables che non contengono nan. La funzione 
compute_tube_segregation_frequency mette a nan gli zeri. Poi quando viene calcolata la matrice di coseg in slice_pairwise.py
o slice_pairwise_inter.py vengono messi a nan tutti i valori corrispondenti a nan nelle frequenze di segregazione. 
Questo essenzialmente produce una coseg in cui non sono presenti zeri (perché l'unico modo per avere uno zero nella coseg è
se è già presente nella seg) ma solo nan
