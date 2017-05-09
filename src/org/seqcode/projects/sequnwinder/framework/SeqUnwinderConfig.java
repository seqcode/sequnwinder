package org.seqcode.projects.sequnwinder.framework;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.RepeatMaskedRegion;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.gsebricks.verbs.location.ChromRegionIterator;
import org.seqcode.gsebricks.verbs.location.ChromosomeGenerator;
import org.seqcode.gsebricks.verbs.location.RepeatMaskedGenerator;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.gseutils.Pair;


/**
 * @author akshaykakumanu
 * @twitter ikaka89
 * @email auk262@psu.edu
 */
public class SeqUnwinderConfig implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public static String version = "0.1.2";
	
	// General options
	protected GenomeConfig gcon;
	protected SequenceGenerator<Region> seqgen = null;
	protected RepeatMaskedGenerator<Region> repMask;
	protected String[] args;
	private double repPropLimit=0.5;
	
	// Feature options (for making the arff file and beyond)
	protected int minK=4;
	protected int maxK=5;
	/** Number K-mers in the model */
	protected int numK;
	protected int win=150;
	protected String outArffFileName="out.arff";
	protected String designFileName="SeqUnwinder.design";
	
	// Peaks features
	/** Peak annotations as given by the user; Enhancer;Shared*/
	protected List<String> annotations = new ArrayList<String>();
	/** List of peaks; this has the same order as the "subGroupsAtPeaks" */
	protected List<Point> peaks = new ArrayList<Point>();
	/** List of regions (win/2 around the peak mid points) */
	protected List<Region> regions = new ArrayList<Region>();
	/** List of sequences around peaks*/
	protected List<String> seqs = new ArrayList<String>();
	/** Flag to include random/background regions in the model */
	protected boolean generateRandRegs=false;
	/** Flag to mask repeats while generating random regions */
	protected boolean screenReps=false;
	/** Random regions generated by SeqUnwinder */ 
	protected List<Region> randRegs = new ArrayList<Region>();
	/** If a subclass has less that this number of instances; notify the user to merge with other subclasses */ 
	public static final int minSubClassSizeNotify = 200;
	/** A fasta file to score the seqs using the k-mer model */
	protected String fastaFname;
	
	
	// Classifier options
	protected double m_Ridge=10;
	/** Augmented Lagrangian parameter rho */
	protected double m_ADMM_pho=1.7;
	protected int m_ADMM_numThreads=5;
	protected int m_SeqUnwinder_MaxIts=15;
	protected int m_ADMM_MaxIts=400;
	protected int sm_NumLayers=2;
	protected int numCrossValidation=3;
	/** Relaxation parameter (to help faster convergence) */
	public static final double ADMM_ALPHA=1.9;
	/** Absolute feasibility tolerance for the primal and dual feasibility conditions */
	public static final double ADMM_ABSTOL = 1E-2; 
	/** Relative  feasibility tolerance for the primal and dual feasibility conditions */
	public static final double ADMM_RELTOL = 1E-2;
	/** Tolerence for internal Nodes convergence */
	public static final double NODES_tol = 1E-2;
	/** The maximum allowed value for pho */
	public static double ADMM_pho_max = 1000;
	/** The options string that is consistent with the WEKA framework */
	protected String[] wekaOptsString;
	protected boolean debugMode = false;
	
	// Discrim options
	/** Minimum length to consider for motif finding */
	protected int minM=6;
	/** Maximum length for motif finding */
	protected int maxM=10;
	/** The base name of the output directory */
	protected String outbase;
	/** The output directory file */
	protected File outdir;
	/** The minimum value of model scan score to consider form motif finding */
	protected double thresold_hills = 0.1;
	/** The number of hills in each cluster to consider for MEME motif finding */
	public static final int NUM_HILLS = 750;
	/** No of iterations of clustering */
	public static final int ITRS_CLUS=10;
	
	// Meme parameters
	protected String MEMEpath;
	protected String MEMEargs = " -dna -mod zoops -revcomp -nostatus ";
	protected int MEMEminw = 6;
	protected int MEMEmaxw = 11;
	protected int MEMEnmotifs = 3;
	protected int MEMEwin = 16;
	public static final double MOTIF_FINDING_ALLOWED_REPETITIVE = 0.2;
	
	// Following are filled during the course of a SeqUnwinder run
	/** Model names **/
	protected List<String> modelNames = new ArrayList<String>();
	/** SubGroup names */
	protected List<String> kmerSubGroupNames = new ArrayList<String>();
	/** Stores the weights of the K-mers
	 *  Keys are the model name
	 * 	Values are the K-mers weights for a given model */ 
	protected HashMap<String,double[]> kmerweights = new HashMap<String,double[]>();
	protected List<WeightMatrix> discrimMotifs = new ArrayList<WeightMatrix>();
	protected HashMap<String,double[]> discrimMotifScores = new HashMap<String,double[]>();
	protected HashMap<String, Double> discrimMotifRocs = new HashMap<String,Double>(); 
	protected String classifier_out ="";
	protected HashMap<String,double[]> testSetStats = new HashMap<String,double[]>();
	protected HashMap<String,double[]> trainSetStats = new HashMap<String,double[]>();
	
	//Settors
	public void resetPeakAnnotations(List<String> modifiedAnnotations){annotations.clear();annotations.addAll(modifiedAnnotations);}
	public void setSubGroupNames(LinkedHashSet<String> sGNs){kmerSubGroupNames.addAll(sGNs);}
	public void setModelNames(List<String> modNames){modelNames.addAll(modNames);}
	public void setNumLayers(int n){
		sm_NumLayers = n;
		wekaOptsString[7] = Integer.toString(sm_NumLayers);
	}
	public void setWeights(HashMap<String,double[]> wts){kmerweights.putAll(wts);}
	public void setDiscrimMotifs(List<WeightMatrix> mots){discrimMotifs.addAll(mots);}
	public void setDiscrimMotifScores(HashMap<String,double[]> scrs){discrimMotifScores.putAll(scrs);}
	public void setClassifierOut(String s){classifier_out = s;}
	public void setDiscrimMotifsRocs(HashMap<String, Double> dscrimRocs){discrimMotifRocs.putAll(dscrimRocs);}
	public void setTrainSetStats(HashMap<String,double[]> trainStats){trainSetStats.putAll(trainStats);}
	public void setTestSetStats(HashMap<String,double[]> tesetStats){testSetStats.putAll(tesetStats);}
	public void setSeqUnwinderMaxIts(int maxItrs){
		m_SeqUnwinder_MaxIts = maxItrs;
		wekaOptsString[13]= Integer.toString(m_SeqUnwinder_MaxIts);
	}
	
	//Gettors
	public List<Point> getPeaks(){return peaks;}
	public List<String> getPeakAnnotations(){return annotations;}
	public List<Region> getRegions(){return regions;}
	public List<String> getSeqs(){return seqs;}
	public Genome getGenome(){return gcon.getGenome();}
	public SequenceGenerator<Region> getSeqGen(){return seqgen;}
	public int getWin(){return win;}
	public String getArffOutName(){return outArffFileName;}
	public int getKmin(){return minK;}
	public int getKmax(){return maxK;}
	public boolean getHasRandRegs(){return generateRandRegs;}
	public List<Region> getGenRandRegs(){return randRegs;}
	public String getDesignFileName(){return designFileName;}
	public String[] getWekaOptions(){return wekaOptsString;}
	public int getNumK(){return numK;}
	public HashMap<String,double[]> getKmerWeights(){return kmerweights;}
	public File getOutDir(){return outdir;}
	public String getOutbase(){return outbase;}
	public List<String> getSubGroupNames(){return kmerSubGroupNames;}
	public List<String> getMNames(){return modelNames;}
	public String getMemeArgs(){String memeargs = MEMEargs+" -nmotifs "+MEMEnmotifs + " -minw "+MEMEminw+" -maxw "+MEMEmaxw; return memeargs;}
	public String getMemePath(){return MEMEpath;}
	public int getMemeSearchWin(){return MEMEwin;}
	public int getMinM(){return minM;}
	public int getMaxM(){return maxM;}
	public int getKmerBaseInd(String kmer){
		int baseInd = 0;
		for(int k=minK; k<kmer.length(); k++){
			baseInd += (int)Math.pow(4, k);
		}
		return baseInd;
	}
	public double getHillsThresh(){return thresold_hills;}
	public List<WeightMatrix> getDiscrimMotifs(){return discrimMotifs;}
	public String getArgs(){
		String a="";
		for(int i=0; i<args.length; i++)
			a = a+" "+args[i];
		return a;
	}
	public String getVersion(){return version;}
	public HashMap<String,Double> getDiscrimMotifRocs(){return discrimMotifRocs;}
	public String getClassifierOutput(){return classifier_out;}
	public HashMap<String,double[]> getTrainSetStats(){return trainSetStats;}
	public HashMap<String,double[]> getTestSetStats(){return testSetStats;}
	public HashMap<String, double[]> getDiscrimMotsScore(){return discrimMotifScores;}
	public String getFastaFname(){return fastaFname;}

	public SeqUnwinderConfig(String[] arguments) throws IOException {
		// Loading general options
		args = arguments;
		ArgParser ap = new ArgParser(args);
		if(ap.hasKey("h") || ap.hasKey("help") || args.length == 0){
			System.err.println(SeqUnwinderConfig.getSeqUnwinderArgsList());
			System.exit(1);
		}
		gcon = new GenomeConfig(args);
		seqgen = gcon.getSequenceGenerator();
		
		// First load all options needed for making the arff file
		
		// Set windown size around peaks
		win = Args.parseInteger(args, "win", 150);

		minK = Args.parseInteger(args, "mink", 4);
		maxK = Args.parseInteger(args, "maxk", 5);

		// Get outdir and outbase and make them; delete dirs that exist with the same
		outbase = Args.parseString(args, "out", "seqUnwinder_out");
		outdir = new File(outbase);
		if(!ap.hasKey("fasta"))
			makeOutPutDirs();

		// use the base to name the arff file
		outArffFileName = outbase+".arff";
		
		numK = 0;
		for(int k=minK; k<=maxK; k++ ){
			numK += (int)Math.pow(4, k);
		}

		// Load peaks and annotations
		if(!ap.hasKey("genregs") && !ap.hasKey("genseqs") && !ap.hasKey("fasta")){
			System.err.println("Please provide genomic locations with labels/annotations and try again !!");
			SeqUnwinderConfig.getSeqUnwinderArgsList();
			System.exit(1);
		}else if(ap.hasKey("genregs")){
			// Reading peaks files and storing annotations
			String peaksFile = ap.getKeyValue("genregs");

			FileReader fr = new FileReader(peaksFile);
			BufferedReader br = new BufferedReader(fr);
			String line;
			while((line = br.readLine()) != null){
				if(!line.startsWith("#") && (line.contains(":") || line.contains("-"))){
					String[] pieces = line.split("\t");
					annotations.add(pieces[1]);
				}
			}
			br.close();
			
			// Loading peaks
			peaks.addAll(RegionFileUtilities.loadPeaksFromPeakFile(gcon.getGenome(), peaksFile, -1));
			
			// Converting peaks to regions
			regions.addAll(RegionFileUtilities.loadRegionsFromPeakFile(gcon.getGenome(), peaksFile, win));
		}else if(ap.hasKey("genseqs")){
			// Reading peaks files and storing annotations
			String seqsFile = ap.getKeyValue("genseqs");

			FileReader fr = new FileReader(seqsFile);
			BufferedReader br = new BufferedReader(fr);
			String line;
			while((line = br.readLine()) != null){
				if(!line.startsWith("#")){
					String[] pieces = line.split("\t");
					annotations.add(pieces[1]);
					seqs.add(pieces[0]);
				}
			}
			br.close();
		}else if(ap.hasKey("fasta")){
			fastaFname = ap.getKeyValue("fasta");
		}
		
		if(ap.hasKey("screenrepeats")){
			screenReps = true;
			repMask = new RepeatMaskedGenerator<Region>(gcon.getGenome()); 

			if(peaks.size() > 0){

				Iterator<Point> peakItr = peaks.iterator();
				Iterator<Region> regionItr = regions.iterator();
				Iterator<String> annotationItr = annotations.iterator();

				while(peakItr.hasNext()){

					Point currPeak = peakItr.next();
					Region currRegion = regionItr.next();
					String currAnnotation = annotationItr.next();

					double repLen = 0;
					Iterator<RepeatMaskedRegion> repItr = repMask.execute(currRegion);
					while (repItr.hasNext()) {
						RepeatMaskedRegion currRep = repItr.next();
						if (currRep.overlaps(currRegion)) {
							repLen += (double) currRep.getWidth();
						}
					}
					if (repLen / (double) currRegion.getWidth() > repPropLimit){
						peakItr.remove();
						regionItr.remove();
						annotationItr.remove();
					}
						

				}
			}


		}		

		//Check if any subclasses have less the minimum number of allowed training instances
		if(ap.hasKey("mergelow")){
			mergeLowCountSubclasses();
		}else{
			removeLowCountSubclasses();
		}

		// Check id random regions are needed to be created 
		if(ap.hasKey("makerandregs")){
			if(ap.hasKey("geninfo") || ap.hasKey("genome") || ap.hasKey("seq")){
				generateRandRegs = true;
			}else{
				System.err.println("Please provide genome files!!!");
				System.exit(1);
			}
		}
		
		
		// Now create randregs if needed
		if(generateRandRegs){

			// First find how many rand regions are needed
			int numRand = Integer.MAX_VALUE; // This should be the size of the subgroup with minimum no of. instances
			Set<String> subgroupNames = new HashSet<String>();
			for(String s : annotations){
				if(subgroupNames.add(s)){}
			}

			for(String s : subgroupNames){
				if(Collections.frequency(annotations, s) < numRand)
					numRand = Collections.frequency(annotations, s);
			}

			RandRegionsGenerator randGenerator = new RandRegionsGenerator(screenReps,numRand);
			randRegs.addAll(randGenerator.execute());
			regions.addAll(randRegs);
			// Now add "Random" to annotations
			for(Region r :  randRegs){
				peaks.add(r.getMidpoint());
				annotations.add("Random");
			}
			
		}
		
		// Now, load all options needed to run the classifier
		// the regularization constant
		debugMode = ap.hasKey("debug");
		m_Ridge = Args.parseDouble(args, "r", 10);
		m_ADMM_pho = Args.parseDouble(args, "pho", 1.7);
		m_ADMM_MaxIts = Args.parseInteger(args, "a", 400);
		m_SeqUnwinder_MaxIts = Args.parseInteger(args, "s", 15);
		m_ADMM_numThreads = Args.parseInteger(args, "threads", 5);
		numCrossValidation = Args.parseInteger(args, "x", 3);
		
		// Now make the weka options string
		if(!debugMode)
			wekaOptsString=new String[20];
		else
			wekaOptsString=new String[21];
		
		wekaOptsString[0] = "-t"; wekaOptsString[1] = outdir.getAbsolutePath()+File.separator+outArffFileName;
		wekaOptsString[2] = "-x"; wekaOptsString[3] = Integer.toString(numCrossValidation);
		wekaOptsString[4] = "-CLS"; wekaOptsString[5] = outdir.getAbsolutePath()+File.separator+designFileName;
		wekaOptsString[6] = "-NL";wekaOptsString[7] = Integer.toString(sm_NumLayers);
		wekaOptsString[8] = "-TY";wekaOptsString[9] = "L1";
		wekaOptsString[10] = "-A";wekaOptsString[11] = Integer.toString(m_ADMM_MaxIts);
		wekaOptsString[12] = "-S";wekaOptsString[13]= Integer.toString(m_SeqUnwinder_MaxIts);
		wekaOptsString[14] = "-PHO";wekaOptsString[15]=Double.toString(m_ADMM_pho);
		wekaOptsString[16] = "-threads";wekaOptsString[17]=Integer.toString(m_ADMM_numThreads);
		wekaOptsString[18] = "-R";wekaOptsString[19]=Double.toString(m_Ridge);
		if(debugMode)
			wekaOptsString[20] = "-DEBUG";

		// Load all MEME arguments
		// Path to MEME binary
		MEMEpath = Args.parseString(args, "memepath", "");
		MEMEargs = Args.parseString(args, "memeargs", MEMEargs);
		
		MEMEminw = Args.parseInteger(args, "mememinw", 6);
		MEMEmaxw = Args.parseInteger(args, "mememaxw", 13);
		//Size of the focussed meme search win
		MEMEwin = Args.parseInteger(args, "memesearchwin", 16);
		MEMEnmotifs = Args.parseInteger(args, "memenmotifs", 3);

		// Load arguments for Discrim analysis
		minM = Args.parseInteger(args, "minscanlen", 6);
		maxM = Args.parseInteger(args, "maxscanlen", 14);
		thresold_hills = Args.parseDouble(args, "hillsthresh", 0.1);
		
		

	}
	
	/**
	 * Removes subclass regions that have a total of less than "minSubClassSizeNotify" train instances
	 */
	public void removeLowCountSubclasses(){
		HashMap<String,Double> subclassCount = new HashMap<String,Double>();

		for(String a : annotations){
			if(subclassCount.containsKey(a)){
				subclassCount.put(a, subclassCount.get(a)+1);
			}else{
				subclassCount.put(a, 1.0);
			}
		}

		List<String> sgroupsToRemove = new ArrayList<String>();

		for(String sname :  subclassCount.keySet()){
			if(subclassCount.get(sname)<minSubClassSizeNotify){
				System.err.println("Wrarning!! -->	"+sname+" has less than "+ Integer.toString(minSubClassSizeNotify)+
						" sites. Removing this subgroup. Use mergeLow option to merge these sites." );
				sgroupsToRemove.add(sname);
			}
		}
		
		//Now remove these annotations peaks,regions and seqs
		for(String toRemove : sgroupsToRemove){
			Iterator<String> a_itr = annotations.iterator();
			Iterator<Point> p_itr = peaks.iterator();
			Iterator<Region> r_itr = regions.iterator();
			Iterator<String> s_itr = seqs.iterator();
			
			while(a_itr.hasNext()){
				String currAnnotation = a_itr.next();
				Point currPeak;
				Region currReg;
				String currSeq;
				if(peaks.size()>0)
					currPeak = p_itr.next();
				if(regions.size()>0)
					currReg = r_itr.next();
				if(seqs.size()>0)
					currSeq = s_itr.next();
				if(currAnnotation.equals(toRemove)){
					a_itr.remove();
					if(peaks.size()>0)
						p_itr.remove();
					if(regions.size()>0)
						r_itr.remove();
					if(seqs.size()>0)
						s_itr.remove();
				}
			}
		}

	}

	/**
	 * Merges all subclass regions that have less than "minSubClassSizeNotify" train instances
	 * This is how the merge happens:-
	 * The regions (of the low number subclass) are randomly assigned to subclasses that atleast share one label in common.
	 */
	public void mergeLowCountSubclasses(){
		HashMap<String,Double> subclassCount = new HashMap<String,Double>();
		
		for(String a : annotations){
			if(subclassCount.containsKey(a)){
				subclassCount.put(a, subclassCount.get(a)+1);
			}else{
				subclassCount.put(a, 1.0);
			}
		}
		
		HashMap<String,List<String>> mergeInds = new HashMap<String,List<String>>();
		
		for(String sname :  subclassCount.keySet()){
			if(subclassCount.get(sname) < minSubClassSizeNotify){
				System.err.println("Wrarning!! -->	"+sname+" has less than "+ Integer.toString(minSubClassSizeNotify)+
						" sites. Merging this subgroup with other relevant subgroups." );
				mergeInds.put(sname,new ArrayList<String>());
				String[] tmpLabs = sname.split(";");
				for(int tl =0; tl<tmpLabs.length; tl++){
					for(String s : subclassCount.keySet()){
						if(s.startsWith(tmpLabs[tl]+";") || s.endsWith(";"+tmpLabs[tl]) || s.contains(";"+tmpLabs[tl]+";") || s.equals(tmpLabs[tl])){
							if(subclassCount.get(s)> minSubClassSizeNotify)
								mergeInds.get(sname).add(s);
						}
					}
				}

			}
		}
		
		// Now update annotations
		if(mergeInds.size()>0){
			Random rand = new Random();
			for(String sname : mergeInds.keySet()){
				if(mergeInds.get(sname).size()>0){
					int maxRand = mergeInds.get(sname).size();
					for(int s = 0; s< annotations.size(); s++){
						if(sname.equals(annotations.get(s))){
							int rind = rand.nextInt(maxRand);
							annotations.set(s, mergeInds.get(sname).get(rind));
						}
					}
				}
			}
		}
		
	}

	public void makeOutPutDirs(){
		//Test if output directory already exists. If it does,  recursively delete contents
		if(outdir.exists())
			deleteDirectory(outdir);
		//(re)make the output directory
		outdir.mkdirs();
		
		// Make the intermediate directory
		File interDir = new File(outdir.getAbsoluteFile()+File.separator+"intermediate");
		interDir.mkdirs();

	}

	public static String getSeqUnwinderArgsList(){
		return(new String("" +
				"Copyright (C) Akshay Kakumanu 2016-2017\n" +
				"\n" +
				"SeqUnwinder comes with ABSOLUTELY NO WARRANTY. This is free software, and you\n"+
				"are welcome to redistribute it under certain conditions.  See the MIT license \n"+
				"for details.\n"+
				"\n OPTIONS:\n" +
				" General:\n"+
				"\t--out <prefix>: Ouput file prefix. All output will be put into a directory with the prefix name\n" +
				"\t--threads <n>: Use n threads to train SeqUnwinder model. Default is 5 threads\n" +
				"\t--debug: Flag to run in debug mode; prints extra output\n" +
				"\t--memepath <path>: path to the meme bin dir (default: meme is in $PATH)\n" +
				" Specify the genome:\n" +
				"\t--geninfo <genome info file> This file should list the lengths of all chromosomes on separate lines using the format chrName<tab>chrLength\n" + 
				"\t\tAND\n" +  
				"\t--seq <path>: A directory containing fasta format files corresponding to every named chromosome is required\n" +
				" Input Genomic Regions:\n" +
				"\t--genregs <List of TF binding sites with annotations; eg: chr1:151736000  Shared;Proximal>\n" +
				"\t\tOR\n" + 
				"\t--genseqs <DNA sequences around at TF binding sites; eg: ATGC...TGC	Shared;Proximal>\n"+
				"\t--win <int>: Size of the genomic regions in bp. Default = 150.\n" +
				"\t--makerandregs: Flag to make random genomic regions as an extra outgroup class in classification (Only applicable when genome is provide.) \n" +
				" SeqUnwinder modelling options \n" +
				"\t--mink <int>: Minimum length of k-mer (default = 4)\n" + 
				"\t--maxk <int>: Maximum length of k-mer (default = 5)\n" + 
				"\t--r <value>: Regularization constant (default = 10)\n" +
				"\t--x <int>: Number of folds for cross validation, default = 3.\n" +
				"\t--mergelow: Flag to merge subclasses with less than 200 sites with other relevant classes. By default, all subclasses with less that 200 sites are removed. \n" +
				" Other SeqUnwinder options (Highly recommend using defaul options): \n"+
				"\t--minscanlen <value>: Minimum length of the window to scan K-mer models. Default=8.\n"+
				"\t--maxscanlen <value>: Maximum length of the window to scan K-mer models. Default=14.\n"+
				"\t--hillsthresh <value>: Scoring threshold to identify hills. Default=0.1.\n"+
				"\t--mememinw <value>: minw arg for MEME. Default=6.\n"+
				"\t--mememaxw <value>: maxw arg for MEME. Default=13. This value should always be less than \"maxscanlen\".\n"+
				"\t‒‒memenmotifs <int>: Number of motifs MEME should find in each condition (default=3)\n" +
				"\t‒‒memeargs <args> : Additional args for MEME (default:  -dna -mod zoops -revcomp -nostatus)\n"+
				"\t--memesearchwin <value>: Window around hills to search for discriminative motifs. Default=16. (Only applicable when run with \"genregs\").\n"+
				"\t--a <int>: Maximum number of allowed ADMM iterations. Default=500.\n"+
				""));
	}


	public boolean deleteDirectory(File path) {
		if( path.exists() ) {
			File[] files = path.listFiles();
			for(int i=0; i<files.length; i++) {
				if(files[i].isDirectory()) {
					deleteDirectory(files[i]);
				}
				else {
					files[i].delete();
				}
			}
		}
		return( path.delete() );
	}

	public class RandRegionsGenerator{

		private int numSamples = 1000;
		private int validSamples=0;

		//private RepeatMaskedGenerator<Region> repMask;
		private double genomeSize=0;
		private long [] chromoSize;
		private String [] chromoNames;
		private int numChroms=0;
		private Random rand = new Random();
		private double repPropLimit=0.5;
		private boolean screenRepeats=false;

		//Settors
		public void setNum(int n){numSamples=n;}
		public void setScreenRepeats(boolean s){screenRepeats=s;}

		public RandRegionsGenerator(boolean screenReps, int num) {
			setScreenRepeats(screenReps);
			setNum(num);
		}


		public List<Region> execute() {

			List<Region>regList = new ArrayList<Region>();
			// First see how big the genome is:
			chromoSize = new long[gcon.getGenome().getChromList().size()];
			chromoNames = new String[gcon.getGenome().getChromList().size()];
			Iterator<NamedRegion> chroms = new ChromRegionIterator(gcon.getGenome());
			while (chroms.hasNext()) {
				NamedRegion currentChrom = chroms.next();
				genomeSize += (double) currentChrom.getWidth();
				chromoSize[numChroms] = currentChrom.getWidth();
				chromoNames[numChroms] = currentChrom.getChrom();
				// System.out.println(chromoNames[numChroms]+"\t"+chromoSize[numChroms]);
				numChroms++;
			}// System.out.println(genomeSize);

			// Now, iteratively generate random positions and check if they are
			// valid
			while (validSamples < numSamples) {
				Region potential;
				long randPos = (long) (1 + (rand.nextDouble() * genomeSize));
				// find the chr
				boolean found = false;
				long total = 0;
				for (int c = 0; c < numChroms && !found; c++) {
					if (randPos < total + chromoSize[c]) {
						found = true;
						if (randPos + win < total + chromoSize[c]) {
							potential = new Region(gcon.getGenome(), chromoNames[c], (int) (randPos - total), (int) (randPos+ win - total - 1));
							boolean regionOK = true;

							// screen repeats
							if (screenRepeats) {
								// is this overlapping a repeat?
								double repLen = 0;
								Iterator<RepeatMaskedRegion> repItr = repMask.execute(potential);
								while (repItr.hasNext()) {
									RepeatMaskedRegion currRep = repItr.next();
									if (currRep.overlaps(potential)) {
										repLen += (double) currRep.getWidth();
									}
								}
								if (repLen / (double) potential.getWidth() > repPropLimit)
									regionOK = false;

								// Is the sequence free from N's?
								String potSeq = seqgen.execute(potential);
								if (potSeq.indexOf('N') >= 0) {
									regionOK = false;
								}
							}
							// Screen dupicates
							for (Region r : regList) {
								if (potential.overlaps(r))
									regionOK = false;
							}

							// Screen for any exclude regions provided
							if (regions.size() != 0) {
								for (Region ex : regions) {
									if (potential.overlaps(ex)) {
										regionOK = false;
									}
								}
							}

							if (regionOK) {
								validSamples++;
								regList.add(potential);
								System.out.println(potential.getChrom() + ":"
										+ potential.getStart() + "-"
										+ potential.getEnd());
							}
						}
					}
					total += chromoSize[c];
				}
			}
			return (regList);
		}


	}

}
