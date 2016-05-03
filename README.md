# muBLASTP

##muBLASTP: Database-indexed Protein Sequence Search for Multi-core CPUs

###Usage: 
* 1 Format FASTA database
```
	$ ./formatdb database/uniprot_sport.fasta
```
* 2 Sort FASTA database by sequence length (would take hours for large databases; 
may use another fasta sort tool, e.g., https://github.com/jimhester/fasta_utilities)
```
    $ ./sortdb database/uniprot_sport.fasta database/uniprot_sport_sort.fasta
```
* 3 Format sorted database 
```
    $ ./formatdb database/uniprot_sport_sort.fasta
```
* 3 Build database index with block of 2k bases
```
    $ ./indexdb database/uniprot_sport_sort.fasta 2000
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
