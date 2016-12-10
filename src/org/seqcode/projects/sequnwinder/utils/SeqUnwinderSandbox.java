package org.seqcode.projects.sequnwinder.utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.genome.sequence.SequenceUtils;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.math.stats.StatUtil;
import org.seqcode.motifs.FreqMatrixImport;


public class SeqUnwinderSandbox {
	/** Learned k-mer weights */
	protected HashMap<String,double[]> weights = new HashMap<String,double[]>();
	protected List<String> kmerModelNames = new ArrayList<String>();
	protected List<WeightMatrix> motifs = new ArrayList<WeightMatrix>();
	protected boolean incRandom = false;
	
	public int minK =4;
	public int maxK =5;
	public int numK = 0;
	
	public static final String[] dimers = {"AA", "AC", "AG", "AT", "CA", "CC", "CG", "GA", "GC", "TA"}; 
	
	//Settors
	public void includeRandom(){incRandom = true;}
	public void setKmin(int mK){minK = mK;}
	public void setKmax(int mK){maxK = mK;}
	public void setNumK(){
		numK = 0;
		for(int k=minK; k<=maxK; k++ ){
			numK += (int)Math.pow(4, k);
		}
		
	}
	public String getKmerName(int ind){
		int currKmerLen = 0;
		ArrayList<Integer> baseinds = new ArrayList<Integer>();
		for(int k=minK; k <= maxK; k++){
			baseinds.add(getKmerBaseInd(k));
		}
		int search = Collections.binarySearch(baseinds, ind);

		currKmerLen = search >= 0 ? minK + search : -1*(search + 1) - 1 + minK;
		String kmerName = RegionFileUtilities.int2seq(ind- getKmerBaseInd(currKmerLen), currKmerLen);

		return kmerName;
	}
	// Load freq matrices
	public void loadMotifsFromFile(String filename) throws NumberFormatException, IOException {
		FreqMatrixImport motifImport = new FreqMatrixImport();
		motifs.addAll(motifImport.readTransfacMatricesAsFreqMatrices(filename));
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
            		weights.put(words[i], new double[numK]);
            	}
            }else{
            	if(!header){
            		System.err.println("Please provide a header in the K-mer weight file");
            		System.exit(1);
            	}
            	int ind = getKmerBaseInd(words[0].length()) + RegionFileUtilities.seq2int(words[0]);
            	for(int i = 1; i < words.length; i++ ){
            		weights.get(kmerModelNames.get(i-1))[ind] = Double.parseDouble(words[i]);
            	}
            }
        }reader.close();

        // Remove random model
        if(!incRandom){
        	// from kmerModelNames
        	Iterator<String> its = kmerModelNames.iterator();
        	while(its.hasNext()){
        		if(its.next().equals("Random")){
        			its.remove();
        		}
        	}
        	weights.remove("Random");
        }
	}
	
	public int getKmerBaseInd(int s){
		int baseInd = 0;
		for(int k=minK; k<s; k++){
			baseInd += (int)Math.pow(4, k);
		}
		return baseInd;
	}
	
	public void scoreMotifs(){
		Set<String> modNames = weights.keySet();
		StringBuilder sb = new StringBuilder();
		sb.append("Motif");sb.append("\t");
		// Print header
		for(String s: modNames){
			sb.append(s);sb.append("\t");
		}
		sb.deleteCharAt(sb.length()-1);
		sb.append("\n");
		for(WeightMatrix mot : motifs){
			sb.append(mot.getName());sb.append("\t");
			List<String> motKmers = WeightMatrix.getConsensusKmerList(mot, minK, maxK);
			
			// Print the k-mers
			System.out.println(mot.getName());
			for(String s : motKmers){
				System.out.println(s);
				System.out.println(SequenceUtils.reverseComplement(s));
			}
			
			for(String modName : modNames){
				double score=0.0;
				for(String s : motKmers){
					String revS = SequenceUtils.reverseComplement(s);
					int indS = RegionFileUtilities.seq2int(s);
					int indRevS = RegionFileUtilities.seq2int(revS);
					int KmerInd  = indS<indRevS ? indS : indRevS;
					int ind = getKmerBaseInd(s.length()) + KmerInd;
					score = score + weights.get(modName)[ind];
				}
				sb.append(score);sb.append("\t");
			}
			sb.deleteCharAt(sb.length()-1);
			sb.append("\n");
		}

		System.out.println(sb.toString());
	}
	
	public void printDimerZscores(){
		// Let model name index be :- m
		// Let dimer name index be :- d
		// (m*dimers.length)+d, will be the value of dimer-score of dimer[d] and model-name kmerModelNames.get(m)
		double[] scores = new double[kmerModelNames.size()*dimers.length];
		double[] zscore = new double[kmerModelNames.size()*dimers.length];
		
		for(int m=0; m<kmerModelNames.size(); m++){ // over k-mer model names
			for(int d=0; d<dimers.length; d++){ // over dimers
				// Calculate the dimer containing k-mer weights
				double score = 0.0;
				for(int kind=0; kind<weights.get(kmerModelNames.get(m)).length; kind++){
					String kmer = getKmerName(kind);
					if(kmer.contains(dimers[d]) || kmer.contains(SequenceUtils.reverseComplement(dimers[d]))){
						score += weights.get(kmerModelNames.get(m))[kind];
					}
				}
				
				scores[m*dimers.length+d] = score;
			}
		}
		
		// Now calculate the z-scores 
		for(int s=0; s<scores.length; s++){
			zscore[s] = (scores[s] - StatUtil.mean(scores))/StatUtil.std(scores);
		}
		
		// Now make the output string to write
		StringBuilder sb = new StringBuilder();
		
		// First print the z-scores
		sb.append("##################### PRINTING Z-SCORES #######################\n");
		// First, make the header
		sb.append("Dimer");sb.append("\t");
		// Now, append all model names
		for(int m=0; m<kmerModelNames.size(); m++){
			sb.append(kmerModelNames.get(m));sb.append("\t");
		}
		sb.deleteCharAt(sb.length()-1);
		sb.append("\n");
		
		for(int d=0; d<dimers.length; d++){
			sb.append(dimers[d]); sb.append("\t");
			for(int m=0; m<kmerModelNames.size(); m++){
				sb.append(zscore[m*dimers.length+d]); sb.append("\t");
			}
			sb.deleteCharAt(sb.length()-1);
			sb.append("\n");
		}
		
		// Second, print actual scores
		sb.append("#\n");
		sb.append("##################### PRINTING DIMER CONTAINING K-MER SUMS #######################\n");
		sb.append("Dimer");sb.append("\t");
		// Now, append all model names
		for(int m=0; m<kmerModelNames.size(); m++){
			sb.append(kmerModelNames.get(m));sb.append("\t");
		}
		sb.deleteCharAt(sb.length()-1);
		sb.append("\n");
		
		for(int d=0; d<dimers.length; d++){
			sb.append(dimers[d]); sb.append("\t");
			for(int m=0; m<kmerModelNames.size(); m++){
				sb.append(scores[m*dimers.length+d]); sb.append("\t");
			}
			sb.deleteCharAt(sb.length()-1);
			sb.append("\n");
		}
		
		sb.deleteCharAt(sb.length()-1);
		
		System.out.println(sb.toString());
	}
	
	
	public static void main(String[] args) throws NumberFormatException, IOException{
		ArgParser ap = new ArgParser(args);
		SeqUnwinderSandbox sbox = new SeqUnwinderSandbox();
		
		sbox.setKmin(Args.parseInteger(args, "minK", 4));
		sbox.setKmax(Args.parseInteger(args, "maxK", 5));
		sbox.setNumK();
		
		if(ap.hasKey("includeRandom"))
			sbox.includeRandom();
		
		
		String weightsFile = ap.getKeyValue("weights");
		sbox.setKmerWeights(weightsFile);
		
		if(ap.hasKey("scoreMotifs")){
			String motsFilename = ap.getKeyValue("motifs");
			sbox.loadMotifsFromFile(motsFilename);
			sbox.scoreMotifs();
		}
		
		if(ap.hasKey("scoreDimers")){
			sbox.printDimerZscores();
		}
		
	}
	
}
