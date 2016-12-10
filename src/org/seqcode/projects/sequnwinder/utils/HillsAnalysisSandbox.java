package org.seqcode.projects.sequnwinder.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.genome.sequence.SequenceUtils;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.gseutils.Pair;
import org.seqcode.motifs.FreqMatrixImport;
import org.seqcode.motifs.MemeER;
import org.seqcode.projects.sequnwinder.clusterkmerprofile.ClusterProfiles;
import org.seqcode.projects.sequnwinder.framework.SeqUnwinderConfig;
import org.seqcode.viz.metaprofile.BinningParameters;


public class HillsAnalysisSandbox {
	
	protected GenomeConfig gCon;
	
	// SeqUnwinder model parameters
	protected int Kmin =4;
	protected int Kmax =5;
	public HashMap<String,double[]> kmerweights = new HashMap<String,double[]>();
	/** Model names */
	protected List<String> kmerModelNames = new ArrayList<String>();
	public int numK;
	
	// Hills clustering parameters
	int numClusters = 3;
	int numClusterItrs = 10;
	// Meme parameters
	protected String MEMEpath;
	protected String MEMEargs = " -dna -mod zoops -revcomp -nostatus ";
	protected int MEMEminw = 6;
	protected int MEMEmaxw = 11;
	protected int MEMEnmotifs = 3;
	protected int MEMEwin = 16;
	public static final double MOTIF_FINDING_ALLOWED_REPETITIVE = 0.2;
	
	
	// Hills
	protected List<Region> hills = new ArrayList<Region>();
	protected List<WeightMatrix> motifs = new ArrayList<WeightMatrix>();
	protected double motifScoreThresh = 0.2;
	
	//CG 
	protected double CGPercCutoff = 0.25;
	
	
	// SeqGenerator
	SequenceGenerator<Region> seqgen = null;
	
	//Settors
	public void setMEMEpath(String s){MEMEpath = s;}
	public void setMEMEminw(int w){MEMEminw = w;}
	public void setMEMEmaxw(int w){MEMEmaxw = w;}
	public void setMEMEnmotifs(int n){MEMEnmotifs = n;}
	public void setMEMEwin(int w){MEMEwin = w;}
	public void setNumClus(int nClus){numClusters = nClus;}
	public void setMaxClusItrs(int maxItrs){numClusterItrs = maxItrs;}
	public void setCGpercCutoff(double cg){CGPercCutoff = cg;}
	public void setMotifThresh(double motThresh){motifScoreThresh =motThresh;}
	public void setKmin(int kmin){Kmin = kmin;}
	public void setKmax(int kmax){Kmax = kmax;}
	public void setNumK(){
		numK = 0;
		for(int k=Kmin; k<=Kmax; k++ ){
			numK += (int)Math.pow(4, k);
		}
	}
	public void setHills(List<Region> modHills){hills = modHills;}
	// Load freq matrices
	public void loadMotifsFromFile(String filename) throws NumberFormatException, IOException {
		motifs.addAll(FreqMatrixImport.readTransfacMatricesAsFreqMatrices(filename));

	}
	public void setKmerWeights(String weightFileName) throws NumberFormatException, IOException{
		BufferedReader reader = new BufferedReader(new FileReader(weightFileName));
		String line;
		boolean header = false;
		while ((line = reader.readLine()) != null) {
			line = line.trim();
			String[] words = line.split("\t");
			if(words[0].charAt(0) == '#' || words[0].contains("Variable") || words[0].contains("Class")){
				header = true;
				for(int i=1; i<words.length; i++){
					kmerModelNames.add(words[i]);
					kmerweights.put(words[i], new double[numK]);
				}
			}else{
				if(!header){
					System.err.println("Please provide a header in the K-mer weight file");
					System.exit(1);
				}
				int ind = getKmerBaseInd(words[0]) + RegionFileUtilities.seq2int(words[0]);
				for(int i = 1; i < words.length; i++ ){
					kmerweights.get(kmerModelNames.get(i-1))[ind] = Double.parseDouble(words[i]);
				}
			}

		}reader.close();
	}

	// Gettors
	public int getKmerBaseInd(String kmer) {
		int baseInd = 0;
		for (int k = Kmin; k < kmer.length(); k++) {
			baseInd += (int) Math.pow(4, k);
		}
		return baseInd;
	}
	public String getMemeArgs(){String memeargs = MEMEargs+" -nmotifs "+MEMEnmotifs + " -minw "+MEMEminw+" -maxw "+MEMEmaxw; return memeargs;}
	
	
	public HillsAnalysisSandbox(GenomeConfig gcon) {
		gCon =gcon;
		seqgen = gCon.getSequenceGenerator();
	}
	
	
	public void printHillStats(String modName){
		System.out.println("#Hill\t"+"HillSeq\t"+modName+"_MotifScore\t"+modName+"_AA-KmerScore\t"+
				modName+"AT-KmerScore\t"+
				modName+"AC-KmerScore\t"+
				modName+"AG-KmerScore\t"+
				modName+"CA-KmerScore\t"+
				modName+"CC-KmerScore\t"+
				modName+"CG-KmerScore\t"+
				modName+"GA-KmerScore\t"+
				modName+"GC-KmerScore\t"+
				modName+"TA-KmerScore\t");

		//First get the k-mer for the given list of motifs
		HashSet<String> motifKmers = new HashSet<String>();
		for(WeightMatrix mot : motifs){
			motifKmers.addAll(WeightMatrix.getConsensusKmerList(mot, Kmin, Kmax));
		}
		
		for(Region hillreg : hills){
			String hillSeq = seqgen.execute(hillreg).toUpperCase();
			double motifScore=0.0, AAScore=0.0, ATScore=0.0, ACScore =0.0, AGScore=0.0, CAScore=0.0, CCScore=0.0, CGScore=0.0, GAScore=0.0, GCScore=0.0, TAScore=0.0;
			for(int k=Kmin; k<=Kmax; k++){
				for(int j=0; j<hillSeq.length()-k+1; j++){
					String currk = hillSeq.substring(j, j+k);
					String revcurrk = SequenceUtils.reverseComplement(currk);
					int  currKInt = RegionFileUtilities.seq2int(currk);
					int  revCurrKInt = RegionFileUtilities.seq2int(revcurrk);
					int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
					int baseInd = this.getKmerBaseInd(currk);
					
					//Motif-Score
					if(motifKmers.contains(currk) || motifKmers.contains(revcurrk) || kmerweights.get(modName)[baseInd+kmer] < 0){
						motifScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//AA-Score
					if(currk.contains("AA") || revcurrk.contains("AA") || kmerweights.get(modName)[baseInd+kmer] < 0){
						AAScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//AT-Score
					if(currk.contains("AT") || revcurrk.contains("AT") || kmerweights.get(modName)[baseInd+kmer] < 0){
						ATScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//AC-Score
					if(currk.contains("AC") || revcurrk.contains("AC") || kmerweights.get(modName)[baseInd+kmer] < 0){
						ACScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//AG-Score
					if(currk.contains("AG") || revcurrk.contains("AG") || kmerweights.get(modName)[baseInd+kmer] < 0){
						AGScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//CA-Score
					if(currk.contains("CA") || revcurrk.contains("CA") || kmerweights.get(modName)[baseInd+kmer] < 0){
						CAScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//CC-Score
					if(currk.contains("CC") || revcurrk.contains("CC") || kmerweights.get(modName)[baseInd+kmer] < 0){
						CCScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//CG-Score
					if(currk.contains("CG") || revcurrk.contains("CG") || kmerweights.get(modName)[baseInd+kmer] < 0){
						CGScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//GA-Score
					if(currk.contains("GA") || revcurrk.contains("GA") || kmerweights.get(modName)[baseInd+kmer] < 0){
						GAScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//GC-Score
					if(currk.contains("GC") || revcurrk.contains("GC") || kmerweights.get(modName)[baseInd+kmer] < 0){
						GCScore += kmerweights.get(modName)[baseInd+kmer];
					}
					//TA-Score
					if(currk.contains("TA") || revcurrk.contains("TA") || kmerweights.get(modName)[baseInd+kmer] < 0){
						TAScore += kmerweights.get(modName)[baseInd+kmer];
					}
				}
			}
			System.out.println(hillreg.toString()+"\t"+hillSeq+"\t"+Double.toString(motifScore)
					+"\t"+Double.toString(Math.round(AAScore*100)/100.0)+
					"\t"+Double.toString(Math.round(ATScore*100)/100.0)+
					"\t"+Double.toString(Math.round(ACScore*100)/100.0)+
					"\t"+Double.toString(Math.round(AGScore*100)/100.0)+
					"\t"+Double.toString(Math.round(CAScore*100)/100.0)+
					"\t"+Double.toString(Math.round(CCScore*100)/100.0)+
					"\t"+Double.toString(Math.round(CGScore*100)/100.0)+
					"\t"+Double.toString(Math.round(GAScore*100)/100.0)+
					"\t"+Double.toString(Math.round(GCScore*100)/100.0)+
					"\t"+Double.toString(Math.round(TAScore*100)/100.0));
		}

	}
	
	
	// 
	public void printClustersAndMotifs(String modName, File outdir) throws IOException{
		// First generate profiles at the hills
		ArrayList<int[]> profiles = new ArrayList<int[]>(); 
		profiles.addAll(getProfilesAtHills(hills));
		List<Double> scores = new ArrayList<Double>();
		scores.addAll(scoreHills(modName));
		HashMap<Integer,String> profileIndices = new HashMap<Integer,String>();
		HashMap<Integer,Double> profileScore = new HashMap<Integer,Double>();
		int ind = 0;
		for(Region r: hills){
			profileIndices.put(ind, r.getLocationString());
			profileScore.put(ind, scores.get(ind));
			ind++;
		}
		
		List<Integer> clusterAssignment = new ArrayList<Integer>();
		ClusterProfiles clusterManager = new ClusterProfiles(numClusterItrs,numClusters,profiles,profileIndices,Kmin,Kmax,profileScore,outdir);
		clusterAssignment = clusterManager.execute();
		String[] randomSequences = new String[MemeER.MOTIF_FINDING_NEGSEQ];
		List<Region> randomRegions = MemeER.randomRegionPick(gCon.getGenome(), null, MemeER.MOTIF_FINDING_NEGSEQ,150);
		for(int i=0; i<randomRegions.size(); i++){
			randomSequences[i] = seqgen.execute(randomRegions.get(i));
		}
		
		List<WeightMatrix> discrimMotifs = new ArrayList<WeightMatrix>();
		
		// Now do meme search on each clusters separately
		String memeargs = getMemeArgs();
		MemeER meme = new MemeER(MEMEpath, memeargs);
		for(int c=0; c<numClusters; c++){ // Over each cluster
			System.err.println("Loading sequences for meme analysis : "+modName+ "Cluster"+c);
			int numHillsLoaded = 0;
			List<String> seqs = new ArrayList<String>();
			for(int p=0; p<hills.size(); p++){
				if(clusterAssignment.get(p) == c && numHillsLoaded < SeqUnwinderConfig.NUM_HILLS){
					numHillsLoaded++;
					String s = seqgen.execute(hills.get(p).getMidpoint().expand(MEMEwin/2));
					if(MemeER.lowercaseFraction(s)<= MOTIF_FINDING_ALLOWED_REPETITIVE){
						seqs.add(s);
					}
				}
				if(numHillsLoaded >= SeqUnwinderConfig.NUM_HILLS)
					break;
			}
			//List<WeightMatrix> selectedMotifs = new ArrayList<WeightMatrix>();
			File meme_outFile = new File(outdir.getAbsolutePath()+File.separator+"Cluster-"+Integer.toString(c)+"_meme");

			System.err.println("Running meme now");
			Pair<List<WeightMatrix>,List<WeightMatrix>> matrices = meme.execute(seqs, meme_outFile, false);
			System.err.println("Finished running meme, Now evaluating motif significance");
			List<WeightMatrix> wm = matrices.car();
			List<WeightMatrix> fm = matrices.cdr();
			// Index for all the selected motifs
			int motInd = 0;

			if(wm.size()>0){
				//Evaluate the significance of the discovered motifs
				double rocScores[] = meme.motifROCScores(wm,seqs,randomSequences);
				System.err.println("MEME results for:" );
				for(int w=0; w<fm.size(); w++){
					if(fm.get(w)!=null){
						System.err.println("\t"+fm.get(w).getName()+"\t"+ WeightMatrix.getConsensus(fm.get(w))+"\tROC:"+String.format("%.2f",rocScores[w]));
					}
					if(rocScores[w] > MemeER.MOTIF_MIN_ROC){
						//selectedMotifs.add(fm.get(w));
						fm.get(w).setName(modName+"_c"+Integer.toString(c)+"_"+Integer.toString(motInd));
						discrimMotifs.add(fm.get(w));
					}
					motInd++;
				}
			}
		}
		
		//First write a transfac file with all the identified discrim motifs
		File motifsTransfac = new File(outdir.getAbsoluteFile()+File.separator+"Discrim_motifs.transfac");
		FileWriter fw = new FileWriter(motifsTransfac);
		BufferedWriter bw = new BufferedWriter(fw);

		for(int m =0; m<discrimMotifs.size(); m++){
			String out = WeightMatrix.printTransfacMatrix(discrimMotifs.get(m),discrimMotifs.get(m).getName().replaceAll("#", ""));
			bw.write(out);
		}
		bw.close();
		
		// Now print the logos
		File motifLogos = new File(outdir.getAbsolutePath()+File.separator+"motif_logos");
		motifLogos.mkdirs();
		// Finally, draw the motif logos
		for(WeightMatrix fm : discrimMotifs){
			File motifFileName = new File(outdir.getAbsolutePath()+File.separator+"motif_logos"+File.separator+fm.getName().replaceAll("#", "")+".png");
			org.seqcode.motifs.DrawMotifs.printMotifLogo(fm, motifFileName, 75, fm.getName(), true);
			File motifFileNameRC = new File(outdir.getAbsolutePath()+File.separator+"motif_logos"+File.separator+fm.getName().replaceAll("#", "")+"_rc.png");
			org.seqcode.motifs.DrawMotifs.printMotifLogo(WeightMatrix.reverseComplement(fm), motifFileNameRC, 75, fm.getName().replaceAll("#", ""), true);
		}
		
	}
	
	// Prints the CG hills that do not score highly for the given list of motifs (usually primary motifs)
	public void printCGhills(String modName){
		System.out.println("#Hill\t"+"HillSeq\t"+modName+"_MotifScore\t"+"fractionCG");
		//First get the k-mer for the given list of motifs
		HashSet<String> motifKmers = new HashSet<String>();
		for(WeightMatrix mot : motifs){
			motifKmers.addAll(WeightMatrix.getConsensusKmerList(mot, Kmin, Kmax));
		}
		
		for(Region hillreg : hills){
			String hillSeq = seqgen.execute(hillreg).toUpperCase();
			double score=0.0;
			for(int k=Kmin; k<=Kmax; k++){
				for(int j=0; j<hillSeq.length()-k+1; j++){
					String currk = hillSeq.substring(j, j+k);
					String revcurrk = SequenceUtils.reverseComplement(currk);
					if(motifKmers.contains(currk) || motifKmers.contains(revcurrk)){
						int  currKInt = RegionFileUtilities.seq2int(currk);
						int  revCurrKInt = RegionFileUtilities.seq2int(revcurrk);
						int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
						int baseInd = this.getKmerBaseInd(currk);
						score = score+kmerweights.get(modName)[baseInd+kmer];
					}
				}
			}
			//Now if th motif score is less check for CGs
			if(score < motifScoreThresh){
				BinningParameters params = new BinningParameters(hillreg.getWidth(), hillreg.getWidth());
				CGScorer scorer = new CGScorer(params,seqgen);
				int MaxCG = (int)hillreg.getWidth()/2;
				CGScoreProfile profile = scorer.execute(hillreg);
				int total = profile.getCGcount();
				double perc = total/(double)MaxCG;
				if(perc > CGPercCutoff)
					System.out.println(hillreg.toString()+"\t"+hillSeq+"\t"+Double.toString(score)+"\t"+Double.toString(perc));
			}
			
		}
	
	}
	
	public void scoreHillsWithCGKmers(String modName){
		System.out.println("#Hill\t"+"HillSeq\t"+modName+"_CGKmerScore");
		
		for(Region hillreg : hills){
			double score=0.0;
			String hillSeq = seqgen.execute(hillreg).toUpperCase();
			for(int k=Kmin; k<=Kmax; k++){
				for(int j=0; j<hillSeq.length()-k+1; j++){
					String currk = hillSeq.substring(j, j+k);
					String revcurrk = SequenceUtils.reverseComplement(currk);
					int  currKInt = RegionFileUtilities.seq2int(currk);
					int  revCurrKInt = RegionFileUtilities.seq2int(revcurrk);
					int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
					int baseInd = this.getKmerBaseInd(currk);
					
					if(currk.contains("CG") || revcurrk.contains("CG") || kmerweights.get(modName)[baseInd+kmer]<0){
						score = score+kmerweights.get(modName)[baseInd+kmer];
					}
				}
			}
			System.out.println(hillreg.toString()+"\t"+hillSeq+"\t"+Double.toString(score));
		}
		
		
	}
	
	private List<Double> scoreHills(String modName){
		List<Double> ret = new ArrayList<Double>();
		for(Region hillreg : hills){
			double score=0.0;
			String hillSeq = seqgen.execute(hillreg).toUpperCase();
			for(int k=Kmin; k<=Kmax; k++){
				for(int j=0; j<hillSeq.length()-k+1; j++){
					String currk = hillSeq.substring(j, j+k);
					String revcurrk = SequenceUtils.reverseComplement(currk);
					int  currKInt = RegionFileUtilities.seq2int(currk);
					int  revCurrKInt = RegionFileUtilities.seq2int(revcurrk);
					int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
					int baseInd = this.getKmerBaseInd(currk);
					score += kmerweights.get(modName)[baseInd+kmer];
				}
			}
			ret.add(score);
		}
		return ret;
		
	}
	
	public void printHillScores(String modName){
		System.out.println("#Hill\t"+"HillSeq\t"+modName+"_MotifScore");
		
		List<Double> scores = new ArrayList<Double>();
		scores.addAll(scoreHills(modName));
		
		for(int r=0; r<hills.size(); r++){
			String hillSeq = seqgen.execute(hills.get(r)).toUpperCase();
			System.out.println(hills.get(r).toString()+"\t"+hillSeq+"\t"+Double.toString(scores.get(r)));
		}
		
	}
		

	
	// Prints all the hills that are high scoring for the given list of motifs
	public void printMotifHills(String modName){
		System.out.println("#Hill\t"+"HillSeq\t"+modName+"_Score\t");
		//First get the k-mer for the given list of motifs
		HashSet<String> motifKmers = new HashSet<String>();
		for(WeightMatrix mot : motifs){
			motifKmers.addAll(WeightMatrix.getConsensusKmerList(mot, Kmin, Kmax));
		}

		for(Region hillreg : hills){
			String hillSeq = seqgen.execute(hillreg).toUpperCase();
			double score=0.0;
			for(int k=Kmin; k<=Kmax; k++){
				for(int j=0; j<hillSeq.length()-k+1; j++){
					String currk = hillSeq.substring(j, j+k);
					String revcurrk = SequenceUtils.reverseComplement(currk);
					int  currKInt = RegionFileUtilities.seq2int(currk);
					int  revCurrKInt = RegionFileUtilities.seq2int(revcurrk);
					int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
					int baseInd = this.getKmerBaseInd(currk);
					if(motifKmers.contains(currk) || motifKmers.contains(revcurrk) || kmerweights.get(modName)[baseInd+kmer] < 0){
						score = score+kmerweights.get(modName)[baseInd+kmer];
					}
				}
			}
			
			if(score > motifScoreThresh){
				System.out.println(hillreg.toString()+"\t"+hillSeq+"\t"+Double.toString(score));
			}
		}
	}
	
	private ArrayList<int[]> getProfilesAtHills(List<Region> rs){
		ArrayList<int[]> ret = new ArrayList<int[]>();
		
		for(Region r: rs){
			int[] pfl = new int[numK];
			String seq = seqgen.execute(r).toUpperCase();
			for(int k=Kmin; k<=Kmax; k++){
				for(int i=0; i<(seq.length()-k+1); i++){
					String currk = seq.substring(i, i+k);
					String revcurrk = SequenceUtils.reverseComplement(currk);
					int currKInt = RegionFileUtilities.seq2int(currk);
					int revcurrKInt = RegionFileUtilities.seq2int(revcurrk);
					int kmer = currKInt<revcurrKInt ? currKInt : revcurrKInt;
					int baseind = getKmerBaseInd(currk);	
					pfl[baseind+kmer]++;
				}
			}
			ret.add(pfl);
		}

		return ret;
	}
	
	/**
	 * For a given list of motifs; this method scans all the hills and checks if 
	 * a motif is present in the hill or not. Finally it prints the number of hills a motif is 
	 * present and the average score of the motif in all the hills it is present.
	 */
	public void printMotifsAtHillsStates(String modName, double cutoff){
		for(WeightMatrix wm : motifs){ // For every motif in the list
			int numHits = 0; // Number of hills at which the motifs scores highly for the given model
			double avgMotScore = 0.0; // Average score of the of the motif at the hills containing it
			HashSet<String> motifKmers = new HashSet<String>();
			motifKmers.addAll(WeightMatrix.getConsensusKmerList(wm, Kmin, Kmax));
			for(Region hillreg : hills ){
				String hillSeq = seqgen.execute(hillreg).toUpperCase();
				double score=0.0;
				for(int k=Kmin; k<=Kmax; k++){
					for(int j=0; j<hillSeq.length()-k+1; j++){
						String currk = hillSeq.substring(j, j+k);
						String revcurrk = SequenceUtils.reverseComplement(currk);
						int  currKInt = RegionFileUtilities.seq2int(currk);
						int  revCurrKInt = RegionFileUtilities.seq2int(revcurrk);
						int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
						int baseInd = this.getKmerBaseInd(currk);
						if(motifKmers.contains(currk) || motifKmers.contains(revcurrk) ){
							score = score+kmerweights.get(modName)[baseInd+kmer];
						}
					}
				}
				if(score > cutoff){
					numHits++;
					avgMotScore += score;
				}
				
			}
			avgMotScore = avgMotScore/numHits;
			
			// Now print the stats
			System.out.println(wm.getName()+"\t"+Integer.toString(numHits)+"\t"+Double.toString(avgMotScore));
		}

	}
	
	
	public static void main(String[] args) throws IOException, ParseException{
		ArgParser ap = new ArgParser(args);
		GenomeConfig gcon = new GenomeConfig(args);
		HillsAnalysisSandbox runner = new HillsAnalysisSandbox(gcon);


		// Length of the smallest K-mer in the K-mer models
		int minK = Args.parseInteger(args, "Kmin", 4);
		runner.setKmin(minK);

		// Length of the largest K-mer in the K-mer models
		int maxK = Args.parseInteger(args, "Kmax", 5);
		runner.setKmax(maxK);

		runner.setNumK();

		// K-mer models file / weights file
		String weights = Args.parseString(args, "weights", null);
		if (weights == null) {
			System.err.println("Provide weights file");
			System.exit(1);
		}
		runner.setKmerWeights(weights);
		
		// Load hills
		String hillsFile = ap.getKeyValue("hillsFile");
		if(hillsFile == null){
			System.err.println("Provide hills file");
			System.exit(1);
		}
		List<Region> regs = RegionFileUtilities.loadRegionsFromPeakFile(gcon.getGenome(), hillsFile, -1);
		runner.setHills(regs);

		// Now load K-mer models
		String motifFile = ap.getKeyValue("motiffile");
		if(ap.hasKey("motiffile"))
			runner.loadMotifsFromFile(motifFile);
		
		double cgThresh = Args.parseDouble(args, "cgThresh", 0.25);
		runner.setCGpercCutoff(cgThresh);
		
		double motThresh = Args.parseDouble(args, "motThresh", 0.2);
		runner.setMotifThresh(motThresh);
		
		String modName = ap.getKeyValue("modname");
		if(modName == null){
			System.err.println("Which modelname to scan? Please provide and re-run");
			System.exit(1);
		}
			
		if(ap.hasKey("scoreHills"))
			runner.printHillScores(modName);
		if(ap.hasKey("printCGhills"))
			runner.printCGhills(modName);
		if(ap.hasKey("printMotHills"))
			runner.printMotifHills(modName);
		if(ap.hasKey("scoreHillsWithCGKmers"))
			runner.scoreHillsWithCGKmers(modName);
		if(ap.hasKey("stats"))
			runner.printHillStats(modName);
		if(ap.hasKey("clusterHills")){
			// Load all MEME arguments
			// Path to MEME binary
			String memepath = Args.parseString(args, "memepath", "");
			runner.setMEMEpath(memepath);
			int mememinw = Args.parseInteger(args, "mememinw", 6);
			runner.setMEMEminw(mememinw);
			int mememaxw = Args.parseInteger(args, "mememaxw", 11);
			runner.setMEMEmaxw(mememaxw);
			//Size of the focussed meme search win
			int memewin = Args.parseInteger(args, "memeSearchWin", 16);
			runner.setMEMEwin(memewin);
			int memenmotifs = Args.parseInteger(args, "memenmotifs", 3);
			runner.setMEMEnmotifs(memenmotifs);
			String outbase = Args.parseString(args, "out", "hills_cluster_out");
			File outdir = new File(outbase);
			outdir.mkdirs();
			
			int numClus_CLUS = Args.parseInteger(args, "numClusters", 3);
			runner.setNumClus(numClus_CLUS);
			
			int numClus_Itrs = Args.parseInteger(args, "numClusItrs", 10);
			runner.setMaxClusItrs(numClus_Itrs);
			
			runner.printClustersAndMotifs(modName, outdir);
		}
		

		if(ap.hasKey("motifAtHills")){
			runner.printMotifsAtHillsStates(modName, motThresh);
		}


	}

}
