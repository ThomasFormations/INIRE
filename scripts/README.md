
### In sillico mutations 

#### Python script to generate mutations

In order to use the python script, the following python packages must be installed
- pysam
- Biopython Bio Seq, SeqRecord and SeqIo  

Th mutate.py script enables to mutate a genome sequence by provinding a list of intervals, and for each interval a given mutation.
3 different mutations are implemented currently : shuffle, inversion and mask.  
The list of intervals should be provided as a bed file with the following format :
```
chr  start  end  name  score  strand mutation
```
- chr, start, end, name, score and strand as in the standard bed specification.  
- mutation can have one of the three defined mutations: shuffle, inversion or mask.

Script usage
```
python muate.py --bed ctcf.bed --genome bosTau9_chr1_1_32Mb.fa --output bosTau9_chr1_1_32Mb_mutated.fa
```
- bed: the bed file
- genome: the path to the fasta file of the input genome
- output: the path to the fasta file for the mutated genome
  
The data files to test the script can be downloaded from here:  
https://genoweb.toulouse.inra.fr/~faraut/INIRE/data


