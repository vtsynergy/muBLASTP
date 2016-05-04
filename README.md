# muBLASTP

##muBLASTP: Database-indexed Protein Sequence Search for Multi-core CPUs

###Usage: 
* 1 Format FASTA database
```
	$ ./formatdb database/uniprot_sport.fasta
```
* 2 Sort FASTA database by sequence length
```
    $ ./sortdb database/uniprot_sport.fasta database/uniprot_sport_sort.fasta
```
* 3 Format sorted database 
```
    $ ./formatdb database/uniprot_sport_sort.fasta
```
* 3 Build database index with block of 128k bases (recommand 128K for multithreading, can use larger block for single thread for better performance)
```
    $ ./indexdb database/uniprot_sport_sort.fasta 128
```
* 4(a). Run muBLASTP with a single thread 
```
    $ ./mublastp -i query/unisprot_sprot/query_100/query_0 -d database/uniprot_sport_sort.fasta
```
* 4(b). Run muBLASTP with multiple threads (e.g., 12 threads) 
```
    $ ./mublastp -i query/unisprot_sprot/query_100/query_batch_1000 -d database/uniprot_sport_sort.fasta -t 12
```
###License
Please refer to the included LICENSE file.
