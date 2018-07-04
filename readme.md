# Protein database search tool (PDBS) 

## INSTALATION

### LINUX

#### BUILD

To build PDBS run the following commands:

    git clone https://github.com/pzuljevic/protein_database_search.git
    cd protein_database_search/
    make

Example:

    make && ./pdbs src/test_data/db/uniprot_sprot.fasta2 50 3 5 1 1 100 0.35 0.90

After running make, an executable named 'pdbs' will appear in the current directory.

#### CLEAN-UP

To clean up the build run:

    make clean
  

## EXAMPLES

### Run:

    ./pdbs <input_file_path> <num threads> <K> <W> <hash type> <MOD> <numLengthClusters> <looseThreshold> <tightThreshold>  <outputFile> <queryInFile>



### Example - generate test data and plot it
>
> make install
> 
> python scripts/gen_data.py > test_data
>
> ./pdbs test_data 1 3 10 1 0 1 0.20 0.90 log.out query.in
>
> python scripts/plot.py
>
> make clean
>


#### Parameter description:
> 
>   K - Kgram length
>
>   W - sliding window size
>
>   hash type - 
>
>               0 (MD5)
>
>               1 (protein alphabet as the number in base 24) 
>
>               2 (same as previous, but with MOD applied)
>
>   MOD - used only for option 2 for hash type
> 
>   numLengthClusters - number of length buckets (when doing phase 1 - clustering by length)
>
>   looseThreshold - threshold to create loose clusters
> 
>   tightThreshold - threshold to create tight clusters
>
>   outputFile - name of the output file (used as an input for plot script)
>
>   queryInFile - file which will be used for online queries, as a result
>                 pdbs will generate 1 file for each sequence query with up to the
>                 10 tight cluster represents from the best loose cluster,
>                 this may be used for the further poa analysis


#### Expected input file description:

Input file must be in the FASTA form and follow the next format:

    line-n:          <FASTA HEADER>
    line-n+1 to n+k: <PROTEIN DATA>
    E.g.:
	>sp|Q6GZX3|002L_FRG3G Uncharacterized protein 002L OS=Frog virus 3 (isolate Goorha) GN=FV3-002L PE=4 SV=1
	MSIIGATRLQNDKSDTYSAGPCYAGGCSAFTPRGTCGKDWDLGEQTCASGFCTSQPLCAR
	IKKTQVCGLRYSSKGKDPLVSAEWDSRGAPYVRCTYDADLIDTQAQVDQFVSMFGESPSL

 
### Example1 (2 threads):

```
>  $ ./pdbs ../src/db/uniprot_sprot.fasta2 2 3 10 0 0 200 0.4 0.9 log.out
>  Reading input file: ../src/db/uniprot_sprot.fasta2, size 0.0247606 GB
>  Using 2 io threads
>  Enqueuing reading file: worker=1/2
>   start_seek_gigabytes=0	end_seek_gigabytes=0.0123803	delta_seek_gigabytes=0.0123803
>   start_seek_bytes=0	end_seek_bytes=13293241	delta_seek_bytes=13293241
>  Enqueuing reading file: worker=2/2
>   start_seek_gigabytes=0.0123803	end_seek_gigabytes=0.0247606	delta_seek_gigabytes=0.0123803
>   start_seek_bytes=13293241	end_seek_bytes=26586482	delta_seek_bytes=13293241
>  Stopping thread pool executor...
>  Processing data for worker: 2 from: 13293241 to: 26586482 took 48.0993 seconds,  read 26299 samples
>  Processing data for worker: 1 from: 0 to: 13293241 took 48.2817 seconds,  read 25324 samples
```

### Example2 (10 threads):

```
>  ./pdbs ../src/db/uniprot_sprot.fasta2 10 3 10 0 0 200 0.4 0.9 log.out
>  Reading input file: ../src/db/uniprot_sprot.fasta2, size 0.0247606 GB
>  Using 10 io threads
>  Enqueuing reading file: worker=1/10
>   start_seek_gigabytes=0	end_seek_gigabytes=0.00247606	delta_seek_gigabytes=0.00247606
>   start_seek_bytes=0	end_seek_bytes=2658648	delta_seek_bytes=2658648
>  Enqueuing reading file: worker=2/10
>   start_seek_gigabytes=0.00247606	end_seek_gigabytes=0.00495212	delta_seek_gigabytes=0.00247606
>   start_seek_bytes=2658648	end_seek_bytes=5317296	delta_seek_bytes=2658648
>  Enqueuing reading file: worker=3/10
>   start_seek_gigabytes=0.00495212	end_seek_gigabytes=0.00742818	delta_seek_gigabytes=0.00247606
>   start_seek_bytes=5317296	end_seek_bytes=7975944	delta_seek_bytes=2658648
>  Enqueuing reading file: worker=4/10
>   start_seek_gigabytes=0.00742818	end_seek_gigabytes=0.00990424	delta_seek_gigabytes=0.00247606
>   start_seek_bytes=7975944	end_seek_bytes=10634592	delta_seek_bytes=2658648
>  Enqueuing reading file: worker=5/10
>   start_seek_gigabytes=0.00990424	end_seek_gigabytes=0.0123803	delta_seek_gigabytes=0.00247606
>   start_seek_bytes=10634592	end_seek_bytes=13293240	delta_seek_bytes=2658648
>  Enqueuing reading file: worker=6/10
>   start_seek_gigabytes=0.0123803	end_seek_gigabytes=0.0148564	delta_seek_gigabytes=0.00247606
>   start_seek_bytes=13293240	end_seek_bytes=15951888	delta_seek_bytes=2658648
>  Enqueuing reading file: worker=7/10
>   start_seek_gigabytes=0.0148564	end_seek_gigabytes=0.0173324	delta_seek_gigabytes=0.00247606
>   start_seek_bytes=15951888	end_seek_bytes=18610536	delta_seek_bytes=2658648
>  Enqueuing reading file: worker=8/10
>   start_seek_gigabytes=0.0173324	end_seek_gigabytes=0.0198085	delta_seek_gigabytes=0.00247606
>   start_seek_bytes=18610536	end_seek_bytes=21269184	delta_seek_bytes=2658648
>  Enqueuing reading file: worker=9/10
>   start_seek_gigabytes=0.0198085	end_seek_gigabytes=0.0222845	delta_seek_gigabytes=0.00247606
>   start_seek_bytes=21269184	end_seek_bytes=23927832	delta_seek_bytes=2658648
>  Enqueuing reading file: worker=10/10
>   start_seek_gigabytes=0.0222845	end_seek_gigabytes=0.0247606	delta_seek_gigabytes=0.00247606
>   start_seek_bytes=23927832	end_seek_bytes=26586480	delta_seek_bytes=2658648
>  Stopping thread pool executor...
>  Processing data for worker: 8 from: 18610536 to: 21269184 took 9.95228 seconds,  read 6580 samples
>  Processing data for worker: 1 from: 0 to: 2658648 took 10.2405 seconds,  read 5416 samples
>  Processing data for worker: 5 from: 10634592 to: 13293240 took 10.335 seconds,  read 4969 samples
>  Processing data for worker: 6 from: 13293240 to: 15951888 took 10.358 seconds,  read 5488 samples
>  Processing data for worker: 7 from: 15951888 to: 18610536 took 10.6365 seconds,  read 4594 samples
>  Processing data for worker: 9 from: 21269184 to: 23927832 took 10.8011 seconds,  read 4978 samples
>  Processing data for worker: 3 from: 5317296 to: 7975944 took 10.9953 seconds,  read 4820 samples
>  Processing data for worker: 10 from: 23927832 to: 26586480 took 11.0448 seconds,  read 4657 samples
>  Processing data for worker: 2 from: 2658648 to: 5317296 took 11.2619 seconds,  read 5030 samples
>  Processing data for worker: 4 from: 7975944 to: 10634592 took 11.9189 seconds,  read 5088 samples
```
