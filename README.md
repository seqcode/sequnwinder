# SeqUnwinder

SeqUnwinder is a framework for characterizing class-discriminative motifs in a collection of genomic loci that have several (overlapping) annotation labels.  


## Download
--------------
The following webpage will maintain executable JAR files for major versions: 
http://mahonylab.org/software/sequnwinder

### Dependencies:
--------------
1. SeqUnwinder requires Java 8+. To build SeqUnwinder, you will also need to download and build the seqcode-core library (https://github.com/seqcode/seqcode-core), and place the build and lib directories on your paths. 
2. SeqUnwinder implements a multi-threaded version of ADMM to train the model. Hence, when using large datasets (tens of thousands of genomic sites), it is advisable to run in a system that allows multiprocessing.
3. SeqUnwinder depends on [MEME](http://meme-suite.org/) (tested with MEME version 4.10.2).
4. if you want to build the code yourself, you will need to first download build the seqcode-core library (https://github.com/seqcode/seqcode-core) and add its build/classes and lib directories to your CLASSPATH.

## Citation:
--------------
TBD

## Running SeqUnwinder
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

1. Specifying the Genome:

  * --__geninfo__ \<genome info file\>:  This file should list the lengths of all chromosomes on separate lines using the format chrName\<tab\>chrLength. You can generate a suitable file from UCSC 2bit format genomes using the UCSC utility “twoBitInfo”. The chromosome names should be exactly the same as those used in your input list of genomic regions. 
   
      The genome info files for some UCSC genome versions:
   | [hg18](http://lugh.bmb.psu.edu/software/multigps/support/hg18.info) | [hg19](http://lugh.bmb.psu.edu/software/multigps/support/hg19.info) | [hg38](http://lugh.bmb.psu.edu/software/multigps/support/hg38.info) | [mm8](http://lugh.bmb.psu.edu/software/multigps/support/mm8.info) | [mm9](http://lugh.bmb.psu.edu/software/multigps/support/mm9.info) | [mm10](http://lugh.bmb.psu.edu/software/multigps/support/mm10.info) | [rn4](http://lugh.bmb.psu.edu/software/multigps/support/rn4.info) | [rn5](http://lugh.bmb.psu.edu/software/multigps/support/rn5.info) | [danRer6](http://lugh.bmb.psu.edu/software/multigps/support/danRer6.info) | [ce10](http://lugh.bmb.psu.edu/software/multigps/support/ce10.info) | [dm3](http://lugh.bmb.psu.edu/software/multigps/support/dm3.info) | [sacCer2](http://lugh.bmb.psu.edu/software/multigps/support/sacCer2.info) | [sacCer3](http://lugh.bmb.psu.edu/software/multigps/support/sacCer3.info) |
  * --__seq__ <fasta seq directory> : A directory containing fasta format files corresponding to every named chromosome is required.

2. Input Genomic Regions:

  *--__GenRegs__ \<Genomic regions with annotations filename\> OR --__GenSeqs__\<Sequences at Genomic regions with annotations filename\> : A tab delimited file of a list of genomic regions/sequences and corresponding annotations/labels. A simple example :
   ```{r, engine='sh', count_lines}
	chr:4487-7865	enhancer;shared
```
   next point 
Major History:
--------------  

Version 0.1 (2016-12-09): Initial release to support manuscript submission.
