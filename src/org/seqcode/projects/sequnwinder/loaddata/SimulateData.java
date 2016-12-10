package org.seqcode.projects.sequnwinder.loaddata;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.seqcode.data.io.BackgroundModelIO;
import org.seqcode.data.motifdb.MarkovBackgroundModel;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.motifs.FreqMatrixImport;
import org.seqcode.projects.sequnwinder.utils.SimulateBindingSite;


public class SimulateData {
	/** This class/code assumes that the first type always has two labels (2 categories; X and Y).
	However, the second type could have either 2 or 3 categories/labels, which is indicated by this variable.  ({A, B, C} or {A, B}) */
	public int num_categories = 3;
	/** Total number of examples in type 2 categories. A, B and C in that order*/
	public int[] num_type2_expls;
	/** Do we need to randomly sample simulated scenarios. If so, generate <code>num_simulations</code> scenarios/datasets  */
	public boolean random_sample_insRate = false;
	public boolean random_sample_overlap = false;
	public Random rand = new Random();
	// Overlap fraction; when using fixed scenarios 
	public double[] overlap_frac;
	public int num_simulations = 100;
	public int[] subclasses;
	public String[] classNames;
	/** Motifs to inset */
	public Map<String,WeightMatrix> motifs = new HashMap<String,WeightMatrix>();
	public GenomeConfig gcon;
	public MarkovBackgroundModel markov;
	public int win=150;
	/**Either this or randomly sample insert-rate freq from 0.3 to 0.9 */
	public double insertRate = 0.7;
	
	public static final double OVERLAP_MIN = 0.5;
	public static final double OVERLAP_MAX = 0.98;
	
	public static final double INSERT_MIN = 0.3;
	public static final double INSERT_MAX = 0.8;
	
	
	public SimulateData(GenomeConfig g) {
		gcon = g;
	}
	
	//Settors
	/**
	 * Sets number of categories for type 2 (<code>num_categories</code>) . Must be either 2 or 3 
	 * @param n
	 */
	public void setNumCategories(int n){
		if(n != 2 && n !=3){
			System.err.println("Invalid number of categories for type II!! Must be either 2 or 3");
			System.exit(1);
		}else{
			num_categories = n;
		}
	}
	/** Sets the number of examples for type 2. This will automatically determine number of exaples in type 1, based on the overlap fraction. */
	public void setNum_type2_expls(int n){
		num_type2_expls = new int[num_categories];
		for(int i=0; i<num_categories; i++){
			num_type2_expls[i] = n;
		}
	}
	/**
	 * f is a double array with fraction overlaps of type 2 with X. 
	 * the length of f should be equal to num_categories.
	 * @param f 
	 */
	public void setSubClasses(double[] f){
		if(f.length != num_categories){
			System.err.println("Make sure the number of categories for Type II are properly set. Chose between 2 or 3");
			System.exit(1);
		}else{
			overlap_frac = new double[f.length];
			for(int i=0; i<f.length; i++){
				overlap_frac[i] = f[i];
			}
			subclasses = new int[2*num_categories];
			classNames = new String[2*num_categories];
			for(int i=0; i<f.length; i++){
				subclasses[i] = Math.round((float)f[i]*num_type2_expls[i]);
				subclasses[i+num_categories] = num_type2_expls[i] - subclasses[i];
			}
			
			if(num_categories == 2){
				classNames = new String[]{"X;A","X;B","Y;A","Y;B"};
			}else if(num_categories == 3){
				classNames = new String[]{"X;A","X;B","X;C","Y;A","Y;B","Y;C"};
			}
			
		}
	}
	public void setRandomSampleInsrtRate(boolean sample){random_sample_insRate = sample;}
	public void setRandomSampleOverlap(boolean sample){random_sample_overlap = sample;}
	public void setMotifs(Map<String,WeightMatrix> ms){motifs = ms;}
	public void setBack(MarkovBackgroundModel back){markov = back;}
	public void setInsrtRate(double rate){insertRate = rate;}
	public void setNumSims(int n){num_simulations = n;}
	public void setWin(int w){win=w;}
	
	public void execute() throws IOException{
		
		// First, check if simulation is needed 
		if(random_sample_overlap && !random_sample_insRate){ // simulate various overlap scenarios
			// generate a random overlap fraction.
			double[] f = new double[num_categories];
			
			for(int s=0; s<num_simulations; s++){
				
				if(num_categories == 2){
					f[0] = OVERLAP_MIN + rand.nextDouble()*(OVERLAP_MAX - OVERLAP_MIN);
					f[1] = 1 - f[0];
				}else if(num_categories == 3){
					f[0] = OVERLAP_MIN + rand.nextDouble()*(OVERLAP_MAX - OVERLAP_MIN);
					f[1] = rand.nextDouble()*(1-f[0]);
					f[2] = 1-f[0]-f[1];
				}
				this.setSubClasses(f);
				StringBuilder sb = new StringBuilder(); sb.append("simulateOverlap_");
				for(int i=0; i<f.length; i++){
					sb.append(classNames[i].replaceAll(";", "")); sb.append("_o");sb.append(Double.toString(f[i]).substring(0, 4));sb.append("_");
				}
				sb.deleteCharAt(sb.length()-1);
				sb.append(".fa");
				String filename = sb.toString();
				
				BufferedWriter bw = new BufferedWriter(new FileWriter(filename));

				bw.write(getData().toString());
				bw.close();

			}	
		}else if(random_sample_insRate && !random_sample_overlap){ // simulate various insert rate scenarios
			
			// Set to default overlap fraction 
			double[] f = new double[num_categories];
			if(num_categories ==  2){
				f[0] = 0.75;
				f[1] = 0.25;
			}else{
				f[0] = 0.75;
				f[1] = 0.12;
				f[2] = 0.13;
			}
			setSubClasses(f);
			for(int s=0; s<num_simulations; s++){
				// generate a random insert rate
				insertRate = INSERT_MIN + rand.nextDouble()*(INSERT_MAX - INSERT_MIN);
			
				StringBuilder sb = new StringBuilder(); sb.append("simulateOverlap_");sb.append("inrtRate_");sb.append(Double.toString(insertRate).substring(0, 4)); sb.append(".fa");
				String filename = sb.toString();
				BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
			
				bw.write(getData().toString());
				bw.close();
			}
		}else if(random_sample_insRate && random_sample_overlap){
			// generate a random overlap fraction.
			double[] f = new double[num_categories];

			for(int s=0; s<num_simulations; s++){

				if(num_categories == 2){
					f[0] = OVERLAP_MIN + rand.nextDouble()*(OVERLAP_MAX - OVERLAP_MIN);
					f[1] = 1 - f[0];
				}else if(num_categories == 3){
					f[0] = OVERLAP_MIN + rand.nextDouble()*(OVERLAP_MAX - OVERLAP_MIN);
					f[1] = rand.nextDouble()*(1-f[0]);
					f[2] = 1-f[0]-f[1];
				}
				this.setSubClasses(f);
				StringBuilder sb = new StringBuilder(); sb.append("simulateOverlap_");
				
				insertRate = INSERT_MIN + rand.nextDouble()*(INSERT_MAX - INSERT_MIN);
				
				for(int i=0; i<f.length; i++){
					sb.append(classNames[i].replaceAll(";", "")); sb.append("_o");sb.append(Double.toString(f[i]).substring(0, 6));sb.append("_");
				}
				sb.deleteCharAt(sb.length()-1);
				
				// Now insert the random insert rate 
				sb.append("_insrtRate_");sb.append(Double.toString(insertRate).substring(0, 4));
				
				sb.append(".fa");
				String filename = sb.toString();

				BufferedWriter bw = new BufferedWriter(new FileWriter(filename));

				bw.write(getData().toString());
				bw.close();

				}
		}else{ // Generate only one data set with a given inser rate and overlap fraction


			// Set to default overlap fraction 
			double[] f = new double[num_categories];
			if(num_categories ==  2){
				f[0] = 0.75;
				f[1] = 0.25;
			}else{
				f[0] = 0.75;
				f[1] = 0.12;
				f[2] = 0.13;
			}
			setSubClasses(f);
			String filename = "simulateOverlap.fa";
			BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
			bw.write(getData().toString());
			bw.close();
		}
	}
	
	
	private StringBuilder getData(){
		// Now generate simulated binding sites for each class
		StringBuilder output = new StringBuilder();
		for(int i=0; i<classNames.length; i++){
			String[] labels = classNames[i].split(";");
			// Add all the motifs for this class into a list
			SimulateBindingSite bs = new SimulateBindingSite(gcon);

			List<WeightMatrix> mots = new ArrayList<WeightMatrix>();

			for(String lab : labels){
				if(motifs.containsKey(lab)){
					mots.add(motifs.get(lab));
				}
			}

			bs.setBack(markov);
			bs.setInsertRate(insertRate);
			bs.setLen(win);
			bs.setMotifs(mots);
			bs.setN(subclasses[i]);

			List<String> seq = bs.execute();
			for(String se : seq){
				output.append(se);output.append("\t");output.append(classNames[i]);output.append("\n");
			}

		}
		output.deleteCharAt(output.length()-1);
		output.append("\n");
		
		return output;
	}
	
	
	
	@SuppressWarnings("static-access")
	public static void main(String[] args) throws IOException, ParseException{
		GenomeConfig gc = new GenomeConfig(args);
		ArgParser ap = new ArgParser(args);
		SimulateData sd = new SimulateData(gc);
		
		// Set motifs
		String motifFilename = ap.getKeyValue("motiffile");
		if(!ap.hasKey("motiffile")){
			System.err.println("Please provide motifs that need to be inserted!!");
			System.exit(0);
		}
		
		String backFile =ap.getKeyValue("back");
		if(!ap.hasKey("back")){
			System.err.println("Please provide a background model file!!");
			System.exit(0);
		}
		
		MarkovBackgroundModel back = BackgroundModelIO.parseMarkovBackgroundModel(backFile, gc.getGenome());
		
		FreqMatrixImport motifImport = new FreqMatrixImport();
		List<WeightMatrix> mots = new ArrayList<WeightMatrix>();
		mots.addAll(motifImport.readTransfacMatricesAsFreqMatrices(motifFilename));
		Map<String,WeightMatrix> motsWithLabs = new HashMap<String,WeightMatrix>();
		
		if(ap.hasKey("sampleOverlap"))
			sd.setRandomSampleOverlap(true);
		if(ap.hasKey("sampleInsrtrate"))
			sd.setRandomSampleInsrtRate(true);
		
		double insertRate = Args.parseDouble(args, "insrtRate", 0.7);
		sd.setInsrtRate(insertRate);
		int numType2 = Args.parseInteger(args, "numType2", 3);
		sd.setNumCategories(numType2);
		int numExamples = Args.parseInteger(args, "N", 3000);
		sd.setNum_type2_expls(numExamples);
		int numSims = Args.parseInteger(args, "numSims", 100);
		sd.setNumSims(numSims);
		int w = Args.parseInteger(args, "win", 150);
		sd.setWin(w);
		
		// Now set motif labels
		if(numType2 == 2){
			motsWithLabs.put("X", mots.get(0)); // Ebox
			motsWithLabs.put("Y", mots.get(1)); //Ebf2
			motsWithLabs.put("A", mots.get(3)); //Oct4
			motsWithLabs.put("B", mots.get(4)); //Onecut2
			
		}else if(numType2 == 3){
			motsWithLabs.put("X", mots.get(0)); // Ebox
			motsWithLabs.put("Y", mots.get(1)); //Ebf2
			motsWithLabs.put("A", mots.get(3)); //Oct4
			motsWithLabs.put("B", mots.get(4)); //Onecut2
			motsWithLabs.put("C", mots.get(5)); // GATA
		}
		
		sd.setMotifs(motsWithLabs);
		sd.setBack(back);
		sd.execute();
		
	}
	
	
}
