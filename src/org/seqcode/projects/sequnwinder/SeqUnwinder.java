package org.seqcode.projects.sequnwinder;

import java.util.HashMap;

import org.seqcode.projects.sequnwinder.framework.Classifier;
import org.seqcode.projects.sequnwinder.framework.SeqUnwinderConfig;
import org.seqcode.projects.sequnwinder.loaddata.MakeArff;
import org.seqcode.projects.sequnwinder.motifs.Discrim;
import org.seqcode.projects.sequnwinder.motifs.ScoreMotif;
import org.seqcode.projects.sequnwinder.utils.Outputwriter;

import weka.classifiers.Evaluation;


public class SeqUnwinder {
	SeqUnwinderConfig seqConfig;
	
	public SeqUnwinder(SeqUnwinderConfig seqCon) {
		seqConfig = seqCon;
	}
	
	// Getters
	public SeqUnwinderConfig getConfig(){return seqConfig;}
	
	public static void main(String[] args) throws Exception{
		
		System.err.println("SeqUnwinder version "+SeqUnwinderConfig.version+"\n\n");
		SeqUnwinderConfig seqCon = new SeqUnwinderConfig(args);
		SeqUnwinder sequnwinder = new SeqUnwinder(seqCon);
		Outputwriter scribe = new Outputwriter(seqCon);
		
		// Now make arff file
		System.err.println("Making Arff file..." );
		MakeArff arffmaker = new MakeArff(sequnwinder.getConfig());
		arffmaker.execute();
		
		//Debug line
		System.exit(1);
		// end
		
		
		// Now run the classifier
		System.err.println("Training SeqUnwinder's core framework..." );
		Classifier classifier = new Classifier(sequnwinder.getConfig().getKmin(),sequnwinder.getConfig().getNumK());
		String classifierout = "";
		try {
			classifierout = Evaluation.evaluateModel(classifier, sequnwinder.getConfig().getWekaOptions());
		} catch (Exception e) {
			if (((e.getMessage() != null) && (e.getMessage().indexOf(
					"General options") == -1))
					|| (e.getMessage() == null)) {
				e.printStackTrace();
			} else {
				System.err.println(e.getMessage());
			}
		}
		
		sequnwinder.getConfig().setModelNames(classifier.getModelNames());
		sequnwinder.getConfig().setClassifierOut(classifierout);
		sequnwinder.getConfig().setWeights(classifier.getWeights());
		scribe.writeKmerWeights();
		
		//Now parse the classifier output and extract the classifier statistics
		HashMap<String,double[]> trainStats = new HashMap<String,double[]>();
		HashMap<String,double[]> testStats = new HashMap<String,double[]>();
		String[] classifierOutputLine = sequnwinder.getConfig().getClassifierOutput().split("\n");
		boolean loadedTrainSetStats=false;
		boolean beginLoadingTrainSetStats = false;
		boolean beginLoadingTestSetStats = false;
		boolean loadedTestSetStats = false;
		for(int s=0; s<classifierOutputLine.length; s++){
			if(!loadedTrainSetStats && classifierOutputLine[s].contains("Detailed Accuracy By Class")){
				beginLoadingTrainSetStats = true;
			}else if(beginLoadingTrainSetStats && !loadedTrainSetStats){
				if(!classifierOutputLine[s].trim().startsWith("Weighted")){
					if(classifierOutputLine[s].equals("") || classifierOutputLine[s].contains("TP Rate"))
						continue;
					String[] pieces = classifierOutputLine[s].trim().split("\\s\\s+");
					trainStats.put(pieces[8], new double[8]);
					for(int p=0; p<pieces.length-1; p++){
						trainStats.get(pieces[8])[p] = Double.parseDouble(pieces[p]);
					}
				}else{
					loadedTrainSetStats=true;
				}
			}else if(loadedTrainSetStats && classifierOutputLine[s].contains("Detailed Accuracy By Class")){
				beginLoadingTestSetStats = true;
			}else if(beginLoadingTestSetStats && !loadedTestSetStats){
				if(!classifierOutputLine[s].trim().startsWith("Weighted")){
					if(classifierOutputLine[s].equals("") || classifierOutputLine[s].contains("TP Rate"))
						continue;
					String[] pieces = classifierOutputLine[s].trim().split("\\s\\s+");
					testStats.put(pieces[8], new double[8]);
					for(int p=0; p<pieces.length-1; p++){
						testStats.get(pieces[8])[p] = Double.parseDouble(pieces[p]);
					}
				}else{
					loadedTestSetStats = true;
				}
			}
			
		}
		sequnwinder.getConfig().setTrainSetStats(trainStats);
		sequnwinder.getConfig().setTestSetStats(testStats);
		
		scribe.writeClssifierOutput();
		
		System.err.println("Converting K-mer models into class-discriminative motifs ..." );

		// Now run Discrim
		Discrim discrim = new Discrim(sequnwinder.getConfig());
		discrim.execute();
		
		scribe.writeDiscrimMotifs();
		
		// Now, score motifs
		ScoreMotif scoremotif = new ScoreMotif(sequnwinder.getConfig());
		scoremotif.execute();
		
		scribe.writeDiscrimScore();
		//scribe.makeDiscrimRHeatmaps();
		scribe.makeDiscrimJavaHeatmaps();
		// Finally, make a html file for output
		scribe.writeHTMLfile();
		
		System.err.println("Finished running SeqUnwinder...");
		
	}

}
