package org.seqcode.projects.sequnwinder.motifs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.genome.sequence.SequenceUtils;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.motifs.FreqMatrixImport;
import org.seqcode.projects.sequnwinder.framework.SeqUnwinderConfig;



public class ScoreMotif {
	
	protected SeqUnwinderConfig seqConfig;
	protected HashMap<String,double[]> scores = new HashMap<String,double[]>();
	
	public ScoreMotif(SeqUnwinderConfig seqcon) {
		seqConfig = seqcon;
	}
	
	public void execute() throws IOException{
		for(WeightMatrix wm : seqConfig.getDiscrimMotifs()){
			scores.put(wm.getName(), new double[seqConfig.getMNames().size()]);
		}
		// Write the motif associated k-mers to a file
		FileWriter fw = new FileWriter(new File(seqConfig.getOutDir().getAbsolutePath()+File.separator+"kmer_at_motigs.tab"));
		BufferedWriter bw = new BufferedWriter(fw);
		StringBuilder sb = new StringBuilder();
		sb.append("#Kmers used to to score motifs"); sb.append("\n");
		
		for(WeightMatrix mot : seqConfig.getDiscrimMotifs()){
			sb.append(mot.getName());sb.append("\n");
			// Make sure you a freq matrix instead of a log-odds matrix
			mot.toFrequency();
			String consensusMotif = WeightMatrix.getConsensus(mot);
			sb.append(consensusMotif);sb.append("\n");
			sb.append(WeightMatrix.printTransfacMatrix(mot,mot.getName()));
			List<String> motKmers = WeightMatrix.getConsensusKmerList(mot, seqConfig.getKmin(), seqConfig.getKmax());
			
			for(String s: motKmers){
				sb.append(s);sb.append("\n");
				sb.append(SequenceUtils.reverseComplement(s));sb.append("\n");
			}
			
			int j=0;
			for(String modName : seqConfig.getMNames()){
				double weight=0.0;
				for(String s : motKmers){
					String revS = SequenceUtils.reverseComplement(s);
					int indS = RegionFileUtilities.seq2int(s);
					int indRevS = RegionFileUtilities.seq2int(revS);
					int KmerInd  = indS<indRevS ? indS : indRevS;
					int ind = seqConfig.getKmerBaseInd(s) + KmerInd;
					weight = weight + seqConfig.getKmerWeights().get(modName)[ind];
				}
				scores.get(mot.getName())[j] = weight;
				j++;
			}
		}
		
		bw.write(sb.toString());
		bw.close();
		seqConfig.setDiscrimMotifScores(scores);
		
	}
	
	/***
	 * Main method only for testing
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException{
		SeqUnwinderConfig scon = new SeqUnwinderConfig(args);
		ArgParser ap = new ArgParser(args);
		List<WeightMatrix> mots = new ArrayList<WeightMatrix>();
		if(ap.hasKey("motifs")){
			mots.addAll(FreqMatrixImport.readTransfacMatricesAsFreqMatrices(ap.getKeyValue("motifs")));
		}
		
		List<String> modnames = new ArrayList<String>();
		HashMap<String,double[]> kmerweights = new HashMap<String,double[]>();
		
		if(ap.hasKey("kmerWeights")){
			BufferedReader br = new BufferedReader(new FileReader(ap.getKeyValue("kmerWeights")));
			String line = null;
			
			while((line=br.readLine())!=null){
				line.trim();
				String[] pieces = line.split("\t");
				if(pieces[0].contains("Kmer")){ // header
					for(int s=1; s<pieces.length; s++){
						modnames.add(pieces[s]);
						kmerweights.put(pieces[s], new double[scon.getNumK()]);
					}
				}else{
					int ind = scon.getKmerBaseInd(pieces[0]) + RegionFileUtilities.seq2int(pieces[0]);
					for(int i = 1; i < pieces.length; i++ ){
						kmerweights.get(modnames.get(i-1))[ind] = Double.parseDouble(pieces[i]);
					}
				}
			}
			
			br.close();
		}
		
		scon.setDiscrimMotifs(mots);
		scon.setModelNames(modnames);
		scon.setWeights(kmerweights);
		
		ScoreMotif scorer = new ScoreMotif(scon);
		scorer.execute();
		HashMap<String,double[]> scores = scon.getDiscrimMotsScore();
		StringBuilder sb = new StringBuilder();
		sb.append("MotifName");sb.append("\t");
		for(String mname: modnames){
			sb.append(mname);sb.append("\t");
		}
		sb.deleteCharAt(sb.length()-1);sb.append("\n");
		for(WeightMatrix mot : mots){
			double[] scrs = scores.get(mot.getName());
			sb.append(mot.getName());sb.append("\t");
			for(int d=0; d<scrs.length; d++){
				sb.append(scrs[d]);sb.append("\t");
			}
			sb.deleteCharAt(sb.length()-1);sb.append("\n");
		}
			
		System.out.println(sb.toString());
	}
	
}
