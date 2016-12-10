package org.seqcode.projects.sequnwinder.loaddata;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.seqcode.data.io.BackgroundModelIO;
import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.motifdb.MarkovBackgroundModel;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.sequence.SequenceUtils;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;

import weka.core.Instances;
import weka.core.converters.ArffSaver;
import weka.core.converters.CSVLoader;


public class MakeArffFromFasta {
	
	/** The design file that contains the relationship between subgroups and labels */
	protected String design;
	protected List<String> subGroupNames = new ArrayList<String>();
	protected HashMap<String,Double> subGroupWeights = new HashMap<String,Double>();
	protected List<String> seqs = new ArrayList<String>();
	protected GenomeConfig gcon;
	public MarkovBackgroundModel markov;
	public Random rand = new Random();
	protected String outArffFileName="out.arff";
	protected int len=150;
	protected int Kmin=4;
	protected int Kmax=6;
	
	public MakeArffFromFasta(GenomeConfig gc) {
		gcon = gc;
	}
	
	/**
	 * Use this settor method if you only have subgroups (simple multi-logistic regression)
	 * @param labs
	 * @param randSeqs
	 */
	public void setLabelsSimple(List<String> labs, List<String> randSeqs){
		subGroupNames.addAll(labs);
		// Now make the design file
		StringBuilder designBuilder = new StringBuilder();
		//First get the header
		designBuilder.append("#Index"+"\t"+
				"Name"+"\t"+
				"isSubGroup" +"\t"+
				"AssignedLabs" + "\t"+
				"AssignedSGs" + "\t"+
				"Layer"+"\n");
		// Now subgroups
		HashMap<String,Integer> subgroupNames = new HashMap<String,Integer>();
		for(String s : subGroupNames){
			if(!subgroupNames.containsKey(s))
				subgroupNames.put(s, 1);
		}
		int index = 0;
		for(String s : subgroupNames.keySet()){
			designBuilder.append(index);designBuilder.append("\t"); // subgroup id
			designBuilder.append(s+"\t"); // subgroup 
			designBuilder.append(1);designBuilder.append("\t"); // subgroup indicator
			designBuilder.append("-"+"\t"); // Assigned labels
			designBuilder.append("-"+"\t"); // Assigned subgroups
			designBuilder.append(0);designBuilder.append("\n"); // Layer
			index++;
		}
		
		// Now add the random label subgroup
		designBuilder.append(index);designBuilder.append("\t"); 
		designBuilder.append("RootRandom"+"\t");
		designBuilder.append(1);designBuilder.append("\t");
		designBuilder.append("-"+"\t");
		designBuilder.append("-"+"\t");
		designBuilder.append(0);designBuilder.append("\n");
		
		designBuilder.deleteCharAt(designBuilder.length()-1);
		design = designBuilder.toString();
		
		// Now add the random-seqs to seqs and "RootRandom" label to subGroupNames
		seqs.addAll(randSeqs);
		for(String se : randSeqs){
			subGroupNames.add("RootRandom");
		}
		
		// Finally, find the weights for each subgroup
		for(String s : subGroupNames){
			if(!subGroupWeights.containsKey(s))
				subGroupWeights.put(s, 1.0);
			else
				subGroupWeights.put(s, subGroupWeights.get(s)+1);
		}
		// Now get the sorted keys
		MapKeyComparator<Double> comp_d = new MapKeyComparator<Double>(subGroupWeights);
		List<String> groupWeightsKeyset = comp_d.getKeyList();
		Collections.sort(groupWeightsKeyset, comp_d);

		for(int i=0; i<groupWeightsKeyset.size()-1; i++){
			subGroupWeights.put(groupWeightsKeyset.get(i), 
					subGroupWeights.get(groupWeightsKeyset.get(groupWeightsKeyset.size()-1))/subGroupWeights.get(groupWeightsKeyset.get(i)));
		}
		subGroupWeights.put(groupWeightsKeyset.get(groupWeightsKeyset.size()-1), 1.0);
	}
	/**
	 * Use this settor method, if you have labels in addition to subgroups
	 * @param labs
	 * @param randSeqs
	 */
	public void setLabels(List<String> labs, List<String> randSeqs){
		HashMap<String,Integer> subgroupNames = new HashMap<String,Integer>();
		HashMap<String,Integer> labelNames = new HashMap<String,Integer>();
		int indSG=0;
		int indLab=0;
		for(String s : labs){
			String[] labels = s.split(";");
			String subgroup = labels[0]+"&"+labels[1];
			if(!subgroupNames.containsKey(subgroup)){
				subgroupNames.put(subgroup, indSG);
				indSG++;
			}
			if(!labelNames.containsKey(labels[0])){
				labelNames.put(labels[0], indLab);
				indLab++;
			}
			if(!labelNames.containsKey(labels[1])){
				labelNames.put(labels[1], indLab);
				indLab++;
			}
			subGroupNames.add(subgroup);
		}
		// Now, adjest the indexes of the labels
		for(String s : labelNames.keySet()){
			labelNames.put(s, labelNames.get(s)+indSG+1);
		}

		// Now make the design file
		StringBuilder designBuilder = new StringBuilder();
		//First get the header
		designBuilder.append("#Index"+"\t"+
				"Name"+"\t"+
				"isSubGroup" +"\t"+
				"AssignedLabs" + "\t"+
				"AssignedSGs" + "\t"+
				"Layer"+"\n");
		// Now subgroups
		int index = 0;
		MapKeyComparator<Integer> comp = new MapKeyComparator<Integer>(subgroupNames);
		List<String> subGs = comp.getKeyList();
		Collections.sort(subGs, comp);
		for(String s : subGs){
			designBuilder.append(index);designBuilder.append("\t"); // subgroup id
			designBuilder.append(s+"\t"); // subgroup 
			designBuilder.append(1);designBuilder.append("\t"); // subgroup indicator
			String[] assignedLabs = s.split("&");
			designBuilder.append(labelNames.get(assignedLabs[0])+","+
					labelNames.get(assignedLabs[1])+"\t"); // Assigned labels
			designBuilder.append("-"+"\t"); // Assigned subgroups
			designBuilder.append(0);designBuilder.append("\n"); // Layer
			index++;
		}

		// Now add the random label subgroup
		designBuilder.append(index);designBuilder.append("\t"); 
		designBuilder.append("RootRandom"+"\t");
		designBuilder.append(1);designBuilder.append("\t");
		designBuilder.append("-"+"\t");
		designBuilder.append("-"+"\t");
		designBuilder.append(0);designBuilder.append("\n");

		index++;
		// Now Labels
		comp = new MapKeyComparator<Integer>(labelNames);
		List<String> LABs = comp.getKeyList();
		Collections.sort(LABs, comp);
		for(String s : LABs){
			designBuilder.append(index);designBuilder.append("\t"); // label id
			designBuilder.append("Root"+s+"\t"); // label name
			designBuilder.append(0);designBuilder.append("\t"); // subgroup indicator
			designBuilder.append("-"+"\t"); // Assigned labels

			for(String subS : subgroupNames.keySet()){ 
				if(subS.startsWith(s+"&") || subS.endsWith("&"+s)){
					designBuilder.append(subgroupNames.get(subS)+","); // Assigned subgroups
				}
			}
			designBuilder.deleteCharAt(designBuilder.length()-1);
			designBuilder.append("\t");
			designBuilder.append(1);designBuilder.append("\n");
			index++;
		}
		designBuilder.deleteCharAt(designBuilder.length()-1);
		design = designBuilder.toString();

		// Now add the random-seqs to seqs and "RootRandom" label to subGroupNames
		seqs.addAll(randSeqs);
		for(String se : randSeqs){
			subGroupNames.add("RootRandom");
		}

		// Finally, find the weights for each subgroup
		for(String s : subGroupNames){
			if(!subGroupWeights.containsKey(s))
				subGroupWeights.put(s, 1.0);
			else
				subGroupWeights.put(s, subGroupWeights.get(s)+1);
		}
		// Now get the sorted keys
		MapKeyComparator<Double> comp_d = new MapKeyComparator<Double>(subGroupWeights);
		List<String> groupWeightsKeyset = comp_d.getKeyList();
		Collections.sort(groupWeightsKeyset, comp_d);

		for(int i=0; i<groupWeightsKeyset.size()-1; i++){
			subGroupWeights.put(groupWeightsKeyset.get(i), 
					subGroupWeights.get(groupWeightsKeyset.get(groupWeightsKeyset.size()-1))/subGroupWeights.get(groupWeightsKeyset.get(i)));
		}
		subGroupWeights.put(groupWeightsKeyset.get(groupWeightsKeyset.size()-1), 1.0);


	}
	
	public void setSeqs(List<String> ss){
		seqs.addAll(ss);
	}
	public void setArffOutFileName(String arffOut){outArffFileName = arffOut;}
	public void setKmin(int kmin){Kmin = kmin;}
	public void setKmax(int kmax){Kmax = kmax;}
	public void setLen(int l){len = l;}
	public void setBack(MarkovBackgroundModel b){markov = b;}
	
	//Getters
	public List<String> getRandSeqs(int n){
		List<String> ret = new ArrayList<String>();
		for(int i=0; i<n; i++){
			String s = generateASeq();
			ret.add(s);
		}
		return ret;
	}
	
	public void execute() throws Exception{
		// First count K-mers and generate the .mat file
		KmerCounter counter = new KmerCounter();
		counter.printPeakSeqKmerRange(Kmin, Kmax);
		
		// Now read the mat file and generate the arff file
		generateArff();
		
		// Also, print the desgin file
		BufferedWriter bw = new BufferedWriter(new FileWriter("SeqUnwinder.design"));
		bw.write(design);
		bw.close();
	}
	
	public void generateArff() throws Exception{

		//
		CSVLoader loader = new CSVLoader();
		// Set options
		loader.setNominalAttributes("last");
		loader.setStringAttributes("");
		loader.setMissingValue("?");
		loader.setFieldSeparator("\t");
		StringBuilder sb = new StringBuilder();
		for(String s : subGroupNames){
			sb.append(s);sb.append(",");
		}
		sb.deleteCharAt(sb.length()-1);
		loader.setNominalLabelSpecs(new String[] {sb.toString()});
		loader.setFile(new File("tmpCounts.mat"));
		Instances data = loader.getDataSet();

		//Set subgroup index
		if(data.classIndex() == -1)
			data.setClassIndex(data.numAttributes()-1);

		//First, get weight index
		int wInd = data.numAttributes()-2;
		// Now set weights
		for(int i=0; i<data.numInstances(); i++){
			double weight = data.instance(i).value(wInd);
			data.instance(i).setWeight(weight);
		}
		// Now delete the weight attribute
		data.deleteAttributeAt(wInd);

		//Save the arff file
		ArffSaver saver = new ArffSaver();
		saver.setFile(new File(outArffFileName));
		saver.setInstances(data);
		saver.writeBatch();

	}
	
	public static void main(String[] args) throws Exception{

		// Reading genmome objects
		GenomeConfig gc = new GenomeConfig(args);
		MakeArffFromFasta arffmaker = new MakeArffFromFasta(gc);
		ArgParser ap = new ArgParser(args);

		String backFile =ap.getKeyValue("back");
		if(!ap.hasKey("back")){
			System.err.println("Please provide a background model file!!");
			System.exit(0);
		}

		MarkovBackgroundModel back = BackgroundModelIO.parseMarkovBackgroundModel(backFile, gc.getGenome());
		arffmaker.setBack(back);
		
		// Get peak window size
		int win = Args.parseInteger(args, "len", 150);
		arffmaker.setLen(win);

		// Reading peaks files
		String SeqsFile = ap.getKeyValue("Seqs");

		FileReader fr = new FileReader(SeqsFile);
		BufferedReader br = new BufferedReader(fr);
		String line;
		List<String> labels = new ArrayList<String>();
		List<String> sequences = new ArrayList<String>();
		while((line = br.readLine()) != null){
			if(!line.contains("#")){
				String[] pieces = line.split("\t");
				labels.add(pieces[1]);
				sequences.add(pieces[0]);
			}
		}
		br.close();
		
		

		// Set Seqs
		arffmaker.setSeqs(sequences);

		// generate random seqs
		List<String> randSeqs = new ArrayList<String>();

		
		// First find how many rand regions are needed
		int numRand = Integer.MAX_VALUE; // This should be the size of the subgroup with minimum no of. instances
		Set<String> subgroupNames = new HashSet<String>();
		for(String s : labels){
			if(subgroupNames.add(s)){}
		}

		for(String s : subgroupNames){
			if(Collections.frequency(labels, s) < numRand)
				numRand = Collections.frequency(labels, s);
		}
		randSeqs.addAll(arffmaker.getRandSeqs(numRand));


		// Now set labels, design file string, and ranSeqs to seqs
		// Check if we have labels
		if(ap.hasKey("simple")){
			arffmaker.setLabelsSimple(labels, randSeqs);
		}else{
			arffmaker.setLabels(labels, randSeqs);
		}

		// Now get K-mer params
		int kmin = Args.parseInteger(args, "Kmin", 4);
		int kmax = Args.parseInteger(args, "Kmax", 6);
		arffmaker.setKmin(kmin);
		arffmaker.setKmax(kmax);

		// Name to write the arff file
		String arffOut = Args.parseString(args, "arffOut", "out.arff");
		if(!arffOut.endsWith(".arff"))
			arffOut = arffOut+".arff";

		arffmaker.setArffOutFileName(arffOut);

		// Now execute: Prints the SeqUnwinder design file and the arff file
		arffmaker.execute();


	}
	
	
	
	public class MapKeyComparator<X extends Number> implements Comparator<String>{
		HashMap<String,X> map;

		public MapKeyComparator(HashMap<String,X> m) {
			map = m;
		}

		public List<String> getKeyList(){
			List<String> ret = new ArrayList<String>();
			for(String s : map.keySet()){
				ret.add(s);
			}
			return ret;

		}

		@Override
		public int compare(String o1, String o2) {
			//return (map.get(o1).intValue()).compareTo(map.get(o2).intValue());
			return Integer.compare(map.get(o1).intValue(), map.get(o2).intValue());
		}

	}

	public class KmerCounter{

		protected String outfilename = "tmpCounts.mat";

		public KmerCounter() {}


		// Print the the kmer counts (for a given range of value k) in the sequences
		// for each peak
		public void printPeakSeqKmerRange(int kmin, int kmax) throws IOException {

			int numK = 0;
			for (int k = kmin; k <= kmax; k++) {
				numK = numK + (int) Math.pow(4, k);
			}
			int[] kmerCounts = new int[numK];

			FileWriter of = new FileWriter(outfilename);
			BufferedWriter bw = new BufferedWriter(of);

			StringBuilder sb = new StringBuilder();
			// Printing the header line
			for (int k = kmin; k <= kmax; k++) {
				int N = (int) Math.pow(4, k);
				for (int i = 0; i < N; i++)
					sb.append("\t"+RegionFileUtilities.int2seq(i, k));
			}
			sb.deleteCharAt(0);
			sb.append("\t"+"Weight"+"\t"+"Label"+"\n");
			//bw.write("Weight"+"\t"+"Label"+"\n");
			bw.write(sb.toString());

			//for (Region r : regs) {
			for(int s=0; s<seqs.size(); s++){
				for (int i = 0; i < numK; i++)
					kmerCounts[i] = 0;
				// Check if the sequence (seq) contains any N's if present ignore
				// them
				if (seqs.get(s).contains("N"))
					continue;

				int ind = 0;
				for (int k = kmin; k <= kmax; k++) {
					for (int i = 0; i < (seqs.get(s).length() - k + 1); i++) {
						String currK = seqs.get(s).substring(i, i + k);
						String revCurrK = SequenceUtils.reverseComplement(currK);
						int currKInt = RegionFileUtilities.seq2int(currK);
						int revCurrKInt = RegionFileUtilities.seq2int(revCurrK);
						int kmer = currKInt < revCurrKInt ? currKInt : revCurrKInt;
						kmerCounts[ind + kmer]++;
					}
					ind = ind + (int) Math.pow(4, k);
				}
				sb = new StringBuilder();
				for (int i = 0; i < numK; i++)
					sb.append("\t" + kmerCounts[i]);
				sb.deleteCharAt(0);
				sb.append("\t"+subGroupWeights.get(subGroupNames.get(s))+"\t"+subGroupNames.get(s)+"\n");
				bw.write(sb.toString());
			}
			bw.close();
		}


	}
	
	
	
	
	private String generateASeq(){
		String s = new String();
		//Preliminary bases
		for(int i=1; i<markov.getMaxKmerLen() && i<=len; i++){
			double prob = rand.nextDouble();
			double sum=0; int j=0;
			while(sum<prob){
				String test = s.concat(int2base(j));
				sum += markov.getMarkovProb(test);
				if(sum>=prob){
					s = test;
					break;
				}
				j++;
			}
		}
		//Remaining bases
		for(int i=markov.getMaxKmerLen(); i<=len; i++){
			String lmer = s.substring(s.length()-(markov.getMaxKmerLen() -1));
			double prob = rand.nextDouble();
			double sum=0; int j=0;
			while(sum<prob){
				String test = lmer.concat(int2base(j));
				sum += markov.getMarkovProb(test);
				if(sum>=prob){
					s =s.concat(int2base(j));
					break;
				}
				j++;
			}
		}
		return s;
	}

	private String int2base(int x){
		String c;
		switch(x){
		case 0:
			c="A"; break;
		case 1:
			c="C"; break;
		case 2:
			c="G"; break;
		case 3:
			c="T"; break;
		default:
			c="N";
		}
		return(c);
	}

}
