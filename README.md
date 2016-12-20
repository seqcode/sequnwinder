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
java -Xmx20G -jar sequnwinder.jar <options - see below>

In the above, the “-Xmx20G” argument tells java to use up to 20GB of memory. If you have installed source code from github, and if all classes are in your CLASSPATH, you can run SeqUnwinder as follows:

java -Xmx20G org.seqcode.projects.sequnwinder.SeqUnwinder <options - see below>

Options
Required/important options are __underscored__.

--__species__




Major History:
--------------  

Version 0.1 (2016-12-09): Initial release to support manuscript submission.
