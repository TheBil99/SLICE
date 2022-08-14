La directory ha 1 file python (settings.py) e 5 cartelle (data, main, miscellaneous, src, starting).

- settings.py contiene i parametri relativi a ciascun dataset;
- data contiene sia i file di input che i file di output;
- src contiene le funzioni alla base di SLICE;
- main contiene i programmi da lanciare per ottenere gli output di SLICE;
- starting contiene due file python necessari per poter lanciare i programmi di main;
- miscellaneous contiene diverse funzioni che usano i risultati di SLICE per analisi dati.


COME FARE I RUN DI SLICE

1) Scaricare i file dai Berlinesi

I berlinesi danno un file in formato .table. Va rinominato e va cambiata l'estensione in .txt. Chiamiamo name_root il nome scelto che precede l'estensione .txt. Lo standard per chiamare name_root è di scrivere prima rawdata, poi un nome identificativo, seguito dalla risoluzione, e poi dal numero di tubi. Per esempio, se stiamo usando un dataset fatto da Iza su cellule mESC a 500kb contenente 300 tubi con 3 nuclear profile per tubo, scrivo rawdata_iza-mesc_500kb_300x3.txt. Sottolineo che l'unico vincolo nello scegliere il nome è di cominciare con rawdata e di usare .txt come estensione, il resto può essere scelto a piacere.

2) Creare cartelle in data

Creare una cartella di nome name_root in data e inserire rawdata_name_root.txt qui dentro. Devono poi essere aggiunte, in data/name_root/, 5 cartelle vuote: beta, NPMI, PI2, PI3, PI3_viewpoint. Dentro PI2 e PI3 devono poi essere create due cartelle vuote dal nome visualization.

3) Aggiornare settings.py

Aprire il file settings.py. Nella linea di codice numero 4 (l-4) scrivere il nome scelto per name_root nella corrispondente variabile. Questa operazione permette ai codici di capire che si vuole lavorare su questo nuovo dataset. Dopodiché, aggiungere un if statement per il nuovo name_root a partire da l-28 (si può copia-incollare uno già presente e cambiare i dati). In questa fase, vanno inseriti i seguenti parametri: effective_NPs_per_tube, resolution, r_cell, h, genome_length, chr_dictionary (rimane F_mean che va aggiunto dopo, per ora lasciarlo a NaN):
	effective_NPs_per_tube corrisponde al numero di Nuclear Profile per tubo se il protocollo è phased, e al doppio del numero di Nuclear Profile per tubo se il protocollo è unphased.
	resolution è la risoluzione del dataset in bp
	r_cell è il raggio nucleare della cellula in micrometri
	h è lo spessore di una Nuclear Profile in micrometri
	genome_length è la lunghezza dell'intero genoma in bp (inclusi chrX e chrY)
	chr_dictionary è mouse_chr_dictionary se il dataset si riferisce ai topi, human_chr_dictionary se si riferisce agli umani
Tutti questi parametri sono forniti dai berlinesi (nota però che effective_NPs_per_tube è una nomenclatura che usiamo solo noi). I parametri alpha e v sono calcolati automaticamente, non bisogna fare niente.

4) Lanciare i programmi in starting

Lanciare create_segregation_pkl.py in starting. Questo crea un file segregation_name_root.pkl in data. Successivamente, lanciare compute_F_mean.py e copiare il risultato (stampato sullo standard output) nella variabile F_mean in settings.py (dell'if statement corrispondente al name_root scelto).


A questo punto la fase preparativa è conclusa, e si possono lanciare i codici in main per il calcolo di SLICE.



5) Lanciare SLICE pairwise

Aprire il programma compute_PI2.py in main e, nella l-91, far partire il comando single_chromosome(chr) con il chr desiderato. Il programma salva le PI2, le PI2 significative e il beta array per il chr scelto. I file sono salvati in data/name_root/beta e data/name_root/PI2. Sullo standard output sono poi stampate info sul calcolo delle PI2 (percentuale di NaN, di valori <0, etc).

6) Lanciare SLICE three-way

Aprire il programma compute_PI3.py in main e, nella l-71, far partire il comando single_chromosome(chr) con il chr desiderato. Il programma salva le PI3 in data/name_root/PI3. Di nuovo, sullo standard output sono stampate info sul calcolo delle PI3.

7) Lanciare SLICE three-way fixing a viewpoint

Aprire il programma compute_PI3_viewpoint.py. In l-62, scrivere il nome della viewpoint (viewpoint_name), il cromosoma di appartenenza (chr) e la posizione genomica della viewpoint in bp (position). Il programma salva le PI3 calcolate rispetto a questa viewpoint in data/PI3_viewpoint. Ancora, lo standard output stampa alcune info sul calcolo.



PROGRAMMI AGGIUNTIVI: MISCELLANEOUS

Nella cartella miscellaneous ho inserito dei programmi che usano le PI2 e PI3 calcolate per ricavare dati e plot utili per l'analisi (per esempio la PI2 media al variare della distanza genomica). I programmi sono commentati, quindi con un po' di pazienza si può ricostruire cosa fanno. I dati vengono salvati in PI2/visualization e PI3/visualization.



FORMATO DI OUTPUT PER PI2 E PI3
Le PI2 (e le PI3 con viewpoint fissato) sono delle matrici simmetriche con diagonale a NaN. Per salvare in modo ottimizzato, salvo solo la triangolare superiore come array orizzontale. Le PI3 sono tensori completamente simmetrici . Per salvarli in modo ottimizzato, ho ideato una procedura di compressione (spiegata in src/utilities, compress_symmetric_tensor e uncompress_symmetric_tensor).


NOTA SULLE PI3
Il calcolo delle PI3 diventa molto pesante a risoluzioni alte. Per questo motivo, i rispettivi codici usano il formato dati a 32bit. Ho verificato che su ENEA cresco la risoluzione massima alla quale si può lavorare è 150kb.
