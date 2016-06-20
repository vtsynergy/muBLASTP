# muBLASTP

##muBLASTP: Database-indexed Protein Sequence Search for Multi-core CPUs

###Usage: 
* 1 Format FASTA database
```
	$ ./formatdb -i <Database> 
```
* 2 Sort FASTA database by sequence length
```
    $ ./sortdb -i <Database> -o <Sorted database>
```
* 3 Build database index with a configurable block size. By default, block size is 128K bases, which is optimal for multithreading. For single thread, enlarging block size (e.g. 1024K bases) can achieve better performance. To enable compressed index, define "COMPRESSED_INDEX" in include/blast.h.
```
    $ ./indexdb -i <Sorted database> [-s Block size in K bases, default 128 (K)]
```
* 4 Run muBLASTP 
```
    $ ./mublastp -i <Query> -d <Sorted database> [-t number of threads]
```
###License
Please refer to the included LICENSE file.

###Acknowledgement
This research was supported in part by the NSF BIGDATA Program via IIS-1247693
and the NSF XPS Program via CCF-1337131. We also acknowledge Advanced Research
Computing at Virginia Tech for the use of their computational resources.
