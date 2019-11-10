package org.seqcode.projects.sequnwinder.loaddata;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.genome.sequence.SequenceUtils;
import org.seqcode.projects.sequnwinder.framework.SeqUnwinderConfig;

import weka.core.Instances;
import weka.core.converters.ArffSaver;
import weka.core.converters.CSVLoader;

public class MakeArff {
	
	protected SeqUnwinderConfig seqConfig;
	/** The design file that contains the relationship between subgroups and labels */
	protected String design;
	/** Weight for each subgroup; estimated based on number of training instances for each subgroup in the training dataset */
	protected HashMap<String,Double> subGroupWeights = new HashMap<String,Double>();
	/** All the subgroup names*/
	protected List<String> subGroupNames = new ArrayList<String>();
	/** subgroup annotations at each each; stored into a list*/
	protected List<String> subGroupsAtPeaks = new ArrayList<String>();

	
	public MakeArff(SeqUnwinderConfig sconf) {
		seqConfig = sconf;
	}
	

	//Setters
	/**
	 * This setter methods does four things :-
	 * 1) Stores subgroup names in "subGroupNames". The format is "label1#label2"
	 * 2) Makes the design file for SeqUnwinder
	 * 3) Appends random regions to peaks and regions and also subGroupNames
	 * 4) Finally, calculates the weights for each subgroup instances.
	 * @param labs
	 * @param randRegs
	 */
	@SuppressWarnings("unchecked")
	public void setLabels(){
		//HashMap<String,Integer> subgroupNames = new HashMap<String,Integer>();
		HashMap<String,Integer> labelNames = new HashMap<String,Integer>();
		for(String s : seqConfig.getPeakAnnotations()){
			String[] labels = s.split(";");
			// Now generate the sub group name
			StringBuilder subgroupSB = new StringBuilder();
			for(int i=0; i<labels.length-1; i++){
				subgroupSB.append(labels[i]);subgroupSB.append("#");
			}
			subgroupSB.append(labels[labels.length-1]);
			String subgroup = subgroupSB.toString();
			subGroupNames.add(subgroup);
			for(String sname : labels){
				if(!labelNames.containsKey(sname)){
					labelNames.put(sname, 1);
				}
			}
			
			subGroupsAtPeaks.add(subgroup);
		}
		
		// ############################################################################
		// Now remove a label that has only one subgroup to it
		HashMap<String,Integer> tmpLabMap = new HashMap<String,Integer>(); 
		for(String sname : subGroupNames){
			for(String piece : sname.split("#")){
				if(tmpLabMap.containsKey(piece)){
					tmpLabMap.put(piece,tmpLabMap.get(piece) + 1);
				}else{
					tmpLabMap.put(piece, 1);
				}
			}
		}
		
		for(String sname : tmpLabMap.keySet()){
			if(tmpLabMap.get(sname) ==1){
				labelNames.remove(sname);
			}
		}
		tmpLabMap.clear();
		
		// Now, adjust the indices of the labels
		int indLab=0;
		for(String s : labelNames.keySet()){
			labelNames.put(s, indLab+subGroupNames.size());
			indLab++;
		}
		// ############################################################################
		
		// Now, if a subgroup has the same name as the label then do this:
		// sname ==> snameOnLy#sname 
		// change in "subGroupsAtPeaks" also
		// And also reset the annotaions in Config. Keep ";" instead of "#"
		for(int i=0; i<subGroupNames.size(); i++){
			String sname = subGroupNames.get(i);
			boolean changeName =false;
			for(String l : labelNames.keySet()){
				if(sname.equals(l)){
					changeName=true;
				}
			}
			if(changeName){
				String modified = sname.concat("Only#").concat(sname);
				subGroupNames.set(i, modified);
				for(int s=0; s< subGroupsAtPeaks.size(); s++){
					if(subGroupsAtPeaks.get(s).equals(sname)){
						subGroupsAtPeaks.set(s, modified);
					}
				}
			}
		}
		
		List<String> modifiedAnnotations = new ArrayList<String>();
		for(String s : subGroupsAtPeaks){
			modifiedAnnotations.add(s.replaceAll("#", ";"));
		}
		
		seqConfig.resetPeakAnnotations(modifiedAnnotations);
		
		// ############################################################################
		// Making design file
		// Set the number of layers to the config file
		if(labelNames.size() == 0){
			seqConfig.setNumLayers(1);
			seqConfig.setSeqUnwinderMaxIts(1);
		}else{
			seqConfig.setNumLayers(2);
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

		for(String s : subGroupNames){
			// Check if this subgroup has parents
			boolean isRoot = true;
			for(String labname: labelNames.keySet()){
				if(s.startsWith(labname+"#") || s.endsWith("#"+labname) || s.contains("#"+labname+"#") || s.equals(labname))
					isRoot = false;
			}
			if(labelNames.size() == 0)
				isRoot = true;

			designBuilder.append(index);designBuilder.append("\t"); // subgroup id
			designBuilder.append(s+"\t"); // subgroup 
			designBuilder.append(1);designBuilder.append("\t"); // subgroup indicator
			if(!isRoot){
				String[] assignedLabs = s.split("#");
				for(String aS : assignedLabs){
					if(labelNames.containsKey(aS))
						designBuilder.append(labelNames.get(aS)+",");// Assigned labels
				}
				designBuilder.deleteCharAt(designBuilder.length()-1);
				designBuilder.append("\t");
			}else{
				designBuilder.append("-"+"\t");
			}
			designBuilder.append("-"+"\t"); // Assigned subgroups
			designBuilder.append(0);designBuilder.append("\n"); // Layer
			index++;
		}
		
		// Now Labels
		MapKeyComparator<Integer> comp = new MapKeyComparator<Integer>(labelNames);
		List<String> LABs = comp.getKeyList();
		Collections.sort(LABs, comp);
		for(String s : LABs){
			designBuilder.append(index);designBuilder.append("\t"); // label id
			designBuilder.append(s+"\t"); // label name
			designBuilder.append(0);designBuilder.append("\t"); // subgroup indicator
			designBuilder.append("-"+"\t"); // Assigned labels

			int subSind = 0;
			for(String subS : subGroupNames){ 
				if(subS.startsWith(s+"#") || subS.endsWith("#"+s) || subS.contains("#"+s+"#") || subS.equals(s)){
					designBuilder.append(subSind+","); // Assigned subgroups
				}
				subSind++;
			}
			designBuilder.deleteCharAt(designBuilder.length()-1);
			designBuilder.append("\t");
			designBuilder.append(1);designBuilder.append("\n");
			index++;
		}
		designBuilder.deleteCharAt(designBuilder.length()-1);
		design = designBuilder.toString();

		// ############################################################################

		// Now, find the weights for each subgroup
		for(String s : subGroupsAtPeaks){
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
		
		// ############################################################################
	}
	//Gettors
	
	
	public void execute() throws Exception{
		
		setLabels();
		// First count K-mers and generate the .mat file
		KmerCounter counter = new KmerCounter();
		counter.printPeakSeqKmerRange(seqConfig.getKmin(), seqConfig.getKmax());
		
		// Now read the mat file and generate the arff file
		generateArff();
		
		// Also, print the desgin file
		BufferedWriter bw = new BufferedWriter(new FileWriter(seqConfig.getOutDir().getAbsolutePath()+File.separator+seqConfig.getDesignFileName()));
		bw.write(design);
		bw.close();
		
		// now up-date cofig with the Subgroup/label names 
		seqConfig.setSubGroupNames(subGroupNames);
	}
	
	public void generateArff() throws Exception{
		
		//
		CSVLoader loader = new CSVLoader();
		// Set options
		loader.setNominalAttributes("last");
		loader.setStringAttributes("");
		loader.setMissingValue("?");
		loader.setFieldSeparator("\t");
		loader.setFile(new File(seqConfig.getOutDir().getAbsolutePath()+File.separator+"tmpCounts.mat"));
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
		saver.setFile(new File(seqConfig.getOutDir().getAbsolutePath()+File.separator+seqConfig.getArffOutName()));
		saver.setInstances(data);
		saver.writeBatch();
		
		
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
			
			FileWriter of = new FileWriter(seqConfig.getOutDir().getAbsolutePath()+File.separator+outfilename);
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
			for(int r=0; r< (seqConfig.getRegions().size() > 0 ? seqConfig.getRegions().size() : seqConfig.getSeqs().size()); r++){
				for (int i = 0; i < numK; i++)
					kmerCounts[i] = 0;

				String seq = (seqConfig.getRegions().size() > 0 ? seqConfig.getSeqGen().execute(seqConfig.getRegions().get(r)).toUpperCase() : seqConfig.getSeqs().get(r));
				// Check if the sequence (seq) contains any N's if present ignore
				// them
				if (seq.contains("N"))
					continue;

				int ind = 0;
				for (int k = kmin; k <= kmax; k++) {
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
				sb = new StringBuilder();
				for (int i = 0; i < numK; i++)
					sb.append("\t" + kmerCounts[i]);
				sb.deleteCharAt(0);
				sb.append("\t"+subGroupWeights.get(subGroupsAtPeaks.get(r))+"\t"+subGroupsAtPeaks.get(r)+"\n");
				bw.write(sb.toString());
			}
			bw.close();
		}
		

	}

}
