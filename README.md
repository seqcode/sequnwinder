# SeqUnwinder

SeqUnwinder is a framework for characterizing class-discriminative motifs in a collection of genomic loci that have several (overlapping) annotation labels.  


Downloading Executables
--------------
The following webpage will maintain executable JAR files for major versions: 
http://mahonylab.org/software/sequnwinder

Building from Source
--------------
If you want to build the code yourself, you will need to first download and build the seqcode-core library (https://github.com/seqcode/seqcode-core) and add its build/classes and lib directories to your CLASSPATH.

Dependencies:
--------------
1. SeqUnwinder requires Java 8+. 
2. SeqUnwinder implements a multi-threaded version of ADMM to train the model. Hence, when using large datasets (tens of thousands of genomic sites), it is advisable to run in a system that allows multiprocessing.
3. SeqUnwinder depends on [MEME](http://meme-suite.org/) (tested with MEME version 4.10.2).

Citation:
--------------
Kakumanu, Akshay, et al. "Deconvolving sequence features that discriminate between overlapping regulatory annotations." bioRxiv (2017): [100511](http://biorxiv.org/content/early/2017/01/15/100511).

Running SeqUnwinder
--------------
On a typical dataset (~20,000 sites and ~8 annotation labels) SeqUnwinder takes a couple of hours to run.

Running from a jar file:

```{r, engine='sh', count_lines}
java -Xmx20G -jar sequnwinder.jar <options - see below>
```

In the above, the “-Xmx20G” argument tells java to use up to 20GB of memory. If you have installed source code from github, and if all classes are in your CLASSPATH, you can run SeqUnwinder as follows:

```{r, engine='sh', count_lines}
java -Xmx20G org.seqcode.projects.sequnwinder.SeqUnwinder <options - see below>
```

Options (Required/important options are in __bold__.)

1. General:

  * --__out__ \<prefix>: Output file prefix. All output will be put into a directory with the prefix name. 
  * --threads \<n\>:  Use n threads to train SeqUnwinder model. Default is 5 threads.
  * --debug: Flag to run in debug mode; prints extra output.
  * --memepath \<path\>: path to the meme bin dir (default: meme is in $PATH).

2. Specifying the Genome:

  * --__geninfo__ \<genome info file\>:  This file should list the lengths of all chromosomes on separate lines using the format chrName\<tab\>chrLength. You can generate a suitable file from UCSC 2bit format genomes using the UCSC utility “twoBitInfo”. The chromosome names should be exactly the same as those used in your input list of genomic regions. 
   
      The genome info files for some UCSC genome versions:  
      | [hg18](http://lugh.bmb.psu.edu/software/multigps/support/hg18.info) | [hg19](http://lugh.bmb.psu.edu/software/multigps/support/hg19.info) | [hg38](http://lugh.bmb.psu.edu/software/multigps/support/hg38.info) | [mm8](http://lugh.bmb.psu.edu/software/multigps/support/mm8.info) | [mm9](http://lugh.bmb.psu.edu/software/multigps/support/mm9.info) | [mm10](http://lugh.bmb.psu.edu/software/multigps/support/mm10.info) | [rn4](http://lugh.bmb.psu.edu/software/multigps/support/rn4.info) | [rn5](http://lugh.bmb.psu.edu/software/multigps/support/rn5.info) | [danRer6](http://lugh.bmb.psu.edu/software/multigps/support/danRer6.info) | [ce10](http://lugh.bmb.psu.edu/software/multigps/support/ce10.info) | [dm3](http://lugh.bmb.psu.edu/software/multigps/support/dm3.info) | [sacCer2](http://lugh.bmb.psu.edu/software/multigps/support/sacCer2.info) | [sacCer3](http://lugh.bmb.psu.edu/software/multigps/support/sacCer3.info) |
  * --__seq__ \<path\> : A directory containing fasta format files corresponding to every named chromosome is required.

3. Input Genomic Regions:

  * --__genregs__ \<file\>: Genomic regions with annotations filename OR --__genseqs__\<file\>: Sequences with annotations filename. A tab delimited file of a list of genomic points/sequences and corresponding annotations/labels. A simple example :
      ```{r, engine='sh', count_lines}
	GenRegs file:
	chr10:100076604	enhancer;shared
	chr6:100316177	promoter;celltypeA

	GenSeqs file:
	ATTGC....TTA	enhancer;shared
	CGTAA....GGT	promoter;celltypeA
      ```
  * --win \<int\>:  Size of the genomic regions in bp. Default = 150.
  * --makerandregs: Flag to make random genomic regions as an extra outgroup class in classification (only applicable when genome is provided).

4. SeqUnwinder Model Options:

  * --mink \<int\>: Minimum length of *K*-mer to consider. Default = 4.
  * --maxk \<int\>: Maximim length of *K*-mer to consider. Default = 5.
   
     For most SeqUnwinder analysis described in the manuscript, *K*-mers of lengths 4 and 5 showed optimal performance. However, with larger datasets (with more data instances for training), maxk can be increased to 6 or 7. 
  * --r \<value\>: Regularization co-efficient in the model. For most SeqUnwinder applications, with ~20k genomic sites and ~6 labels and *K*-mers of 4 and 5, a value of 10.0 has been very effective. However, the optimal value could change with datasets. One might want to use a range of values and choose the one that performs best (in terms of test accuracy).
  * --x \<int\>: Number of folds for cross validation. Default = 3.
  * --mergelow: Flag to merge subclasses with fewer than 200 sites with other relevant classes. By default, all subclasses with less than 200 sites are removed.
  
5. Other SeqUnwinder options (Highly recommend using defaul options):

  * --minscanlen \<value\>: Minimum length of the window to scan *K*-mer models. Default=8.
  * --maxscanlen \<value\>: Maximum length of the window to scan *K*-mer models. Default=14.
  * --hillsthresh \<value\>: Scoring threshold to identify hills. Default=0.1.
  * --mememinw \<value\>: minw arg for MEME. Default=6.
  * --mememaxw \<value\>: maxw arg for MEME. Default=13. This value should always be less than "maxScanLen".
  * --memenmotifs \<int\>: Number of motifs MEME should find in each condition (default=3)
  * --memeargs \<args\> : Additional args for MEME (default:  -dna -mod zoops -revcomp -nostatus)
  * --memesearchwin \<value\>: Window around hills to search for discriminative motifs. Default=16
  * --a \<int\>: Maximum number of allowed ADMM iterations. Default=400.


Example
--------------
This example runs SeqUnwinder v0.1.2 on simulated dataset used in the SeqUnwinder manuscript (Figure 2B, right panel). Simulated sequences to run this example can be found [here](http://lugh.bmb.psu.edu/software/sequnwinder/simulateOverlap.seqs).

Command:
```{r, engine='sh', count_lines}
java -Xmx20G -jar sequnwinder.jar --out simulateOverlap_sequnwinder_04212017 --threads 10 --debug --memepath path-to-meme --geninfo mm10.info --seq path-to-genomes/mm10/ --GenSeqs simulateOverlap.seqs --win 150 --minK 4 --maxK 5 --R 10 --X 3 --maxScanLen 15
```

Results can be found [here](http://lugh.bmb.psu.edu/software/sequnwinder/example/SeqUnwinder_results.html)

Contact
--------------

For queries, please contact Akshay (auk262@psu.edu) or Shaun Mahony (mahony@psu.edu).

Major History:
--------------  

Version 0.1.2 (2017-05-08): Several minor updates over the previous version. SeqUnwinder now automatically estimates "k" for k-means clustering of hills. Additional option provided to deal with sub-classes with very few training instances (see --mergelow). Several option names have been reformatted for consistency.

Version 0.1 (2016-12-09): Initial release to support manuscript submission.
