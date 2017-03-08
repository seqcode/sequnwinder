package org.seqcode.projects.sequnwinder.framework;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.io.parsing.FASTAStream;
import org.seqcode.genome.sequence.SequenceUtils;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Pair;

public class ScoreSites {

	protected SeqUnwinderConfig seqConfig;
	protected double[][] m_Par; // also has the intercept as the first term; [NumPredictors+Intercept][NumClasses]
	protected List<String> classnames = new ArrayList<String>();
	protected List<Integer[]> data = new ArrayList<Integer[]>();
	protected List<String> seqNames = new ArrayList<String>();

	public ScoreSites(SeqUnwinderConfig sc) {
		seqConfig = sc;
	}

	//Settors
	public void setParams(String file) throws IOException{ // Intercept term should be first and the order should be maintained
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line = null;
		int count=0;
		while((line=br.readLine()) != null){
			String[] pieces = line.split("\t");
			if(line.startsWith("#Kmer") || line.startsWith("#")){
				for(int i=1; i<pieces.length; i++){
					classnames.add(pieces[i]);
				}
				m_Par = new double[seqConfig.getNumK()+1][classnames.size()];
			}else{
				for(int i=1; i<pieces.length; i++){
					m_Par[count][i-1] = Double.parseDouble(pieces[i]);
				}
				count++;
			}
		}
		br.close();
	}


	public void execute(){
		double[] prob = new double[classnames.size()], v = new double[classnames.size()];

		for(int i=0; i<data.size(); i++){
			StringBuilder sb = new StringBuilder();
			sb.append(seqNames.get(i));sb.append("\t");
			// Log-posterior before normalizing
			for (int j = 0; j < classnames.size() ; j++) {
				for (int k = 0; k <= seqConfig.getNumK()+1; k++) {
					v[j] += m_Par[k][j] * data.get(i)[k];
				}
			}

			// Do so to avoid scaling problems
			for (int m = 0; m < classnames.size(); m++) {
				double sum = 0;
				for (int n = 0; n < classnames.size(); n++) {
					sum += Math.exp(v[n] - v[m]);
				}
				prob[m] = 1 / (sum);
				sb.append(prob[m]);sb.append("\t");
			}
			sb.deleteCharAt(sb.length()-1);
			System.out.println(sb.toString());
		}
	}

	public void loadData() throws IOException {

		FASTAStream stream = new FASTAStream(new File(seqConfig.getFastaFname()));
		Integer[] kmerCounts = new Integer[seqConfig.getNumK()];
		
		while(stream.hasNext()){
			for (int i = 0; i < seqConfig.getNumK(); i++)
				kmerCounts[i] = 0;
			Pair<String,String> currSeq = stream.next();


			String seq = currSeq.cdr().toUpperCase();
			// Check if the sequence (seq) contains any N's if present ignore
			// them
			if (seq.contains("N"))
				continue;

			seqNames.add(currSeq.car());

			int ind = 0;
			for (int k = seqConfig.getKmin(); k <= seqConfig.getKmax(); k++) {
				for (int i = 0; i < (seq.length() - k + 1); i++) {
					String currK = seq.substring(i, i + k);
					String revCurrK = SequenceUtils.reverseComplement(currK);
					int currKInt = RegionFileUtilities.seq2int(currK);
					int revCurrKInt = RegionFileUtilities.seq2int(revCurrK);
					int kmer = currKInt < revCurrKInt ? currKInt : revCurrKInt;
					kmerCounts[ind + kmer]++;
				}
				ind = ind + (int) Math.pow(4, k);
			}
			Integer[] kmerCountsWithIntercept = new Integer[kmerCounts.length+1];
			kmerCountsWithIntercept[0] = 1;
			for(int i=0; i<kmerCounts.length; i++){
				kmerCountsWithIntercept[i+1] = kmerCounts[i];
			}

			data.add(kmerCountsWithIntercept);

		}
	}

	public static void main(String[] args) throws IOException{
		ScoreSites runner = new ScoreSites(new SeqUnwinderConfig(args));
		ArgParser ap = new ArgParser(args);
		String model = ap.getKeyValue("model");
		runner.setParams(model);
		runner.loadData();
		runner.execute();
	}


}
