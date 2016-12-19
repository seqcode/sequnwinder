# SeqUnwinder

SeqUnwinder is a framework for characterizing class-discriminative motifs in a collection of genomic loci that have several (overlapping) annotation labels.  


## Download
--------------
The following webpage will maintain executable JAR files for major versions: 
http://mahonylab.org/software/sequnwinder

Dependencies:
--------------
1. SeqUnwinder requires Java 8+. To build SeqUnwinder, you will also need to download and build the seqcode-core library (https://github.com/seqcode/seqcode-core), and place the build and lib directories on your paths. 
2. SeqUnwinder implements a multi-threaded version of ADMM to learn the model. Hence, when using large datasets (tens of thousands of genomic sites), it advisable to run in a system that allows multiprocessing.
3. SeqUnwinder depends on [MEME](http://meme-suite.org/) (tested with MEME version 4.10.2).

Citation:
--------------
TBD


Major History:
--------------  

Version 0.1 (2016-12-09): Initial release to support manuscript submission.
