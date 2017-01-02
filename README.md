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
TBD

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

1. Specifying the Genome:

  * --__geninfo__ \<genome info file\>:  This file should list the lengths of all chromosomes on separate lines using the format chrName\<tab\>chrLength. You can generate a suitable file from UCSC 2bit format genomes using the UCSC utility “twoBitInfo”. The chromosome names should be exactly the same as those used in your input list of genomic regions. 
   
      The genome info files for some UCSC genome versions:  
      | [hg18](http://lugh.bmb.psu.edu/software/multigps/support/hg18.info) | [hg19](http://lugh.bmb.psu.edu/software/multigps/support/hg19.info) | [hg38](http://lugh.bmb.psu.edu/software/multigps/support/hg38.info) | [mm8](http://lugh.bmb.psu.edu/software/multigps/support/mm8.info) | [mm9](http://lugh.bmb.psu.edu/software/multigps/support/mm9.info) | [mm10](http://lugh.bmb.psu.edu/software/multigps/support/mm10.info) | [rn4](http://lugh.bmb.psu.edu/software/multigps/support/rn4.info) | [rn5](http://lugh.bmb.psu.edu/software/multigps/support/rn5.info) | [danRer6](http://lugh.bmb.psu.edu/software/multigps/support/danRer6.info) | [ce10](http://lugh.bmb.psu.edu/software/multigps/support/ce10.info) | [dm3](http://lugh.bmb.psu.edu/software/multigps/support/dm3.info) | [sacCer2](http://lugh.bmb.psu.edu/software/multigps/support/sacCer2.info) | [sacCer3](http://lugh.bmb.psu.edu/software/multigps/support/sacCer3.info) |
  * --__seq__ <fasta seq directory> : A directory containing fasta format files corresponding to every named chromosome is required.

2. Input Genomic Regions:

  * --__GenRegs__ \<Genomic regions with annotations filename\> OR --__GenSeqs__\<Sequences at Genomic points with annotations filename\> : A tab delimited file of a list of genomic points/sequences and corresponding annotations/labels. A simple example :
      ```{r, engine='sh', count_lines}
	GenRegs file:
	chr10:100076604	enhancer;shared
	chr6:100316177	promoter;celltypeA

	GenSeqs file:
	ATTGC....TTA	enhancer;shared
	CGTAA....GGT	promoter;celltypeA
      ```
  * --win \<integer value\>:  Size of the genomic regions in bp. Default = 150.
  * --makerandregs: Flag to make random genomic regions as an extra outgroup class in classification.

3. SeqUnwinder Model Options:

  * --minK \<value\>: Minimum length of *K*-mer to consider. Default = 4.
  * --maxK \<value\>: Maximim length of *K*-mer to consider. Default = 5.
   
     For most SeqUnwinder analysis described in (TBD), *K*-mers of lengths 4 and 5 showed optimal performace. However, with larger datsets (with more data instances for training), maxk can be increased to 6 or 7. 
  * --R \<value\>: Regularization co-efficient in the model. For most SeqUnwinder applications, with on an average of ~20k genomic sites and ~6 labels and *K*-mers of 4 and 5, a value of 10.0 has been very effective. However, the optimal value could change with datasets. One might want to use a range of values and choose the one that performs best (in terms of test accuracy).
  * --X <value>: Number of folds for cross validation. Default = 3.
  * --minScanLen \<value\>: Minimum length of the window to scan *K*-mer models. Default=8.
  * --maxScanLen \<value\>: Maximum length of the window to scan *K*-mer models. Default=14.
  * --hillsThresh \<value\>: Scoring threshold to identify hills. Default=0.1.
  * --mememinw \<value\>: minw arg for MEME. Default=6.
  * --mememaxw \<value\>: maxw arg for MEME. Default=13. This value should always be less than "maxScanLen".
  * --memenmotifs \<value\>: Number of motifs MEME should find for each subclass/lable. Default=3.
  * --memeSearchWin \<value\>: Window around hills to search for discriminative motifs. Default=16
  * --memepath \<path\>: Path to the meme bin dir (default: meme in $PATH) 
  * --numClusters \<value\>: Number of clusters to split k-mer hills using k-means. Default=3.
   
      Following are some training parameters, which we higly recommend not to change.
  * --PHO \<value\>: Augmented Lagrangian parameter for training using ADMM framework. Defaule=1.7.
  * --A \<value\>: Maximum number of allowed ADMM iterations. Default=500.
  * --S \<value\>: Maximum number of allowed SeqUnwinder iterations. Default=15.

4. General Options

  * --out \<String\>: Perfix for the output files.
  * --threads \<value\>: Number of threads to use. Default=4.
  * --debug : Flag to run in debug mode.

Example
--------------
This example runs SeqUnwinder v0.1 on simulated dataset used in the SeqUnwinder manuscript (Figure 2B). Simulated sequences to run this example can be found [here](http://lugh.bmb.psu.edu/software/sequnwinder/simulateOverlap.fa).

Command:
```{r, engine='sh', count_lines}
java -Xmx20G -jar sequnwinder.jar ‒‒geninfo ~/genomes/mm9/mm9.info --seq ~/genomes/mm10/ --seqs simulateOverlap.fa --win 150 --minK 4 --maxK 5 --R 10 --PHO 1.7 --A 500 --S 15 --X 3 --out example --numClusters 2 --memepath ~/software/meme_4.10.2/bin --memenmotifs 3 --mememinw 6 --mememaxw 10 --minScanLen 10 --maxScanLen 15 --hillsThresh 0.1 --threads 5
```

Results can be found [here](http://lugh.bmb.psu.edu/software/sequnwinder/example/SeqUnwinder_results.html)

Contact
--------------

For queries, please contact Akshay (auk262@psu.edu) or Shaun Mahony (mahony@psu.edu).

Major History:
--------------  

Version 0.1 (2016-12-09): Initial release to support manuscript submission.
