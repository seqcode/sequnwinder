package org.seqcode.projects.sequnwinder.framework;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.projects.sequnwinder.framework.LOne;
import org.seqcode.projects.sequnwinder.framework.Optimizer;
import org.seqcode.projects.sequnwinder.framework.Classifier.ClassRelationStructure.*;

import weka.classifiers.AbstractClassifier;
import weka.classifiers.pmml.producer.LogisticProducerHelper;

import weka.core.*;
import weka.core.Capabilities.Capability;
import weka.core.pmml.PMMLProducer;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.NominalToBinary;
import weka.filters.unsupervised.attribute.RemoveUseless;
import weka.filters.unsupervised.attribute.ReplaceMissingValues;
/**
 * Classifier: Extends Weka's Abstract Classifier class to implements the core of SeqUnwinder's Multi-class structured
 * classification framework. 
 * @author akshaykakumanu
 * @version	%I%, %G%
 */

public class Classifier extends AbstractClassifier implements OptionHandler, WeightedInstancesHandler, PMMLProducer{

	// Model parameters compatible with a normal Multinomial logit model
	
	//protected SeqUnwinderConfig seqConfig;
	
	protected int minK = 4;
	
	protected int numK = 0;

	private static final long serialVersionUID = -6260376163872744102L;

	/** The coefficients of the optimized Structured Multinomial logit model. However, this only stores the leaf node parameters*/
	protected double[][] m_Par; // [NumPredictors+Intercept][NumClasses]

	/** The training data saved as a matrix */
	protected double[][] m_Data; // [NumInstances][NumPredictors+Intercept]

	/** The number of attributes in the model */
	protected int m_NumPredictors;

	/** The number of class lables or the number of leaf nodes */
	protected int m_NumClasses;

	/** The regularization parameter, Is multiplied to the logit part of the loss function */
	protected double m_Ridge = 10;

	/** The index of the class attribute. Usually the last attribute */
	protected int m_ClassIndex;

	/** An attribute filter */
	private RemoveUseless m_AttFilter;

	/** The filter used to make attributes numeric. */
	private NominalToBinary m_NominalToBinary;

	/** The filter used to get rid of missing values. */
	private ReplaceMissingValues m_ReplaceMissingValues;

	/** Log-likelihood of the searched model */
	protected double m_LL;

	/** The maximum number of iterations. */
	private int m_BGFS_MaxIts = -1;

	private Instances m_structure;

	private String m_OptimizationType = "L1";

	// Model parameters of the structured multinomial logit not compatible with the weka logit class

	/** All model parameters including the internal non-leaf nodes */
	protected double[][] sm_Par; //[NumPredictors+Intercept][NumNodes]

	/** Total number of nodes in the structured representation of the classes. Including the root node */
	protected int sm_NumNodes;

	protected int sm_NumLayers;

	/** Class Structure relationships */
	public ClassRelationStructure sm_ClassStructure;

	protected boolean sm_Debug;

	protected int m_ADMM_MaxIts = 1000;

	protected double m_ADMM_pho=0.001;

	protected int m_ADMM_numThreads = 5;

	protected int m_SeqUnwinder_MaxIts = 10;
	

	// Setters
	/**
	 * Sets the Ridge parameter
	 * @param ridge
	 */
	public void setRidge(double ridge){
		m_Ridge = ridge;
	}
	/**
	 * Sets the maximum nuber of training rounds while lerning
	 * @param newMaxIts
	 */
	public void setBGFSMaxIts(int newMaxIts){
		m_BGFS_MaxIts = newMaxIts;
	}
	public void serADMMMaxItrs(int ADMMmax){m_ADMM_MaxIts = ADMMmax;}

	public void setSeqUnwinderMaxItrs(int SeqUnwinderMaxIts){m_SeqUnwinder_MaxIts = SeqUnwinderMaxIts;}

	public void setSmDebug(boolean smD){
		sm_Debug = smD;
	}	

	//Gettors

	/**
	 * Returns the ridge co-eff of the model
	 * @return
	 */
	public double getRidge(){
		return m_Ridge;
	}
	/**
	 * Return the value of m_MaxIts
	 * @return
	 */
	public int getBGFSMaxIts(){
		return m_BGFS_MaxIts;
	}
	public int getADMMMaxItrs(){
		return m_ADMM_MaxIts;
	}
	public boolean getSmDebug(){
		return sm_Debug;
	}

	public int getKmerBaseInd(String kmer){
		int baseInd = 0;
		for(int k=minK; k<kmer.length(); k++){
			baseInd += (int)Math.pow(4, k);
		}
		return baseInd;
	}
	
	public List<String> getModelNames(){
		List<String> ret = new ArrayList<String>();
		for(Node n : sm_ClassStructure.allNodes.values()){
			ret.add(n.nodeName);
		}
		return ret;
	}

	public HashMap<String,double[]> getWeights(){
		HashMap<String,double[]> ret = new HashMap<String,double[]>();
		for(int n=0; n<sm_ClassStructure.allNodes.size(); n++){
			ret.put(sm_ClassStructure.allNodes.get(n).nodeName, new double[numK]);
		}
		
		int j=1;
		for (int i = 0; i < m_structure.numAttributes(); i++) {
			if (i != m_structure.classIndex()) {

				int base = getKmerBaseInd(m_structure.attribute(i).name());
				int kind =  base + RegionFileUtilities.seq2int(m_structure.attribute(i).name());
				for (Node n : sm_ClassStructure.allNodes.values()) {
					ret.get(n.nodeName)[kind] = sm_Par[j][n.nodeIndex];
				}
				j++;
			}
		}

		return ret;
	}

	/**
	 * Returns the coefficients of the leaf nodes of the model. The first dimension
	 * indexes the attributes, and the second the classes.
	 *  
	 * @return the coefficients for this logistic model
	 */
	public double[][] coefficients() {
		return m_Par;
	}
	/**
	 * Returns the coefficients of all the nodes of the trained model
	 * @return
	 */
	public double[][] coefficientsFull(){
		return sm_Par;
	}
	
	public Classifier(int mK, int nK) {
		minK= mK;
		numK= nK;
	}
	
	public Classifier() {
	}


	/**
	 * Returns default capabilities of the classifier.
	 * @return the capabilities of this classifier
	 */
	@Override
	public Capabilities getCapabilities() {
		Capabilities result = super.getCapabilities();
		result.disableAll();

		// attributes
		result.enable(Capability.NOMINAL_ATTRIBUTES);
		result.enable(Capability.NUMERIC_ATTRIBUTES);
		result.enable(Capability.DATE_ATTRIBUTES);
		result.enable(Capability.MISSING_VALUES);

		// class
		result.enable(Capability.NOMINAL_CLASS);
		result.enable(Capability.MISSING_CLASS_VALUES);

		return result;
	}


	@Override
	public double[] distributionForInstance(Instance instance) throws Exception {
		m_ReplaceMissingValues.input(instance);
		instance = m_ReplaceMissingValues.output();
		m_AttFilter.input(instance);
		instance = m_AttFilter.output();
		m_NominalToBinary.input(instance);
		instance = m_NominalToBinary.output();

		// Extract the predictor columns into an array
		double[] instDat = new double[m_NumPredictors + 1];
		int j = 1;
		instDat[0] = 1;
		for (int k = 0; k <= m_NumPredictors; k++) {
			if (k != m_ClassIndex) {
				instDat[j++] = instance.value(k);
			}
		}

		double[] distribution = evaluateProbability(instDat);
		return distribution;
	}

	private double[] evaluateProbability(double[] data) {
		double[] prob = new double[m_NumClasses], v = new double[m_NumClasses];

		// Log-posterior before normalizing
		for (int j = 0; j < m_NumClasses ; j++) {
			for (int k = 0; k <= m_NumPredictors; k++) {
				v[j] += m_Par[k][j] * data[k];
			}
		}

		// Do so to avoid scaling problems
		for (int m = 0; m < m_NumClasses; m++) {
			double sum = 0;
			for (int n = 0; n < m_NumClasses; n++) {
				sum += Math.exp(v[n] - v[m]);
			}
			prob[m] = 1 / (sum);
		}

		return prob;
	}


	@Override
	public void buildClassifier(Instances train) throws Exception {

		// can classifier handle the data?
		getCapabilities().testWithFail(train);

		// remove instances with missing class
		train = new Instances(train);
		train.deleteWithMissingClass();

		// Replace missing values
		m_ReplaceMissingValues = new ReplaceMissingValues();
		m_ReplaceMissingValues.setInputFormat(train);
		train = Filter.useFilter(train, m_ReplaceMissingValues);

		// Remove useless attributes
		m_AttFilter = new RemoveUseless();
		m_AttFilter.setInputFormat(train);
		train = Filter.useFilter(train, m_AttFilter);

		// Transform attributes
		m_NominalToBinary = new NominalToBinary();
		m_NominalToBinary.setInputFormat(train);
		train = Filter.useFilter(train, m_NominalToBinary);

		// Save the structure for printing the model
		m_structure = new Instances(train, 0);

		// Extract data
		m_ClassIndex = train.classIndex();
		m_NumClasses = train.numClasses();

		int nK = m_NumClasses; // All K classes needed unlike Weka's logistic framework
		int nR = m_NumPredictors = train.numAttributes() - 1; // Minus 1 to remove the class attribute
		int nC = train.numInstances();

		m_Data = new double[nC][nR + 1]; // Data values
		int[] Y = new int[nC]; // Class labels
		double[] xMean = new double[nR + 1]; // Attribute means
		double[] xSD = new double[nR + 1]; // Attribute stddev's
		double[] sY = new double[nK]; // Number of classes
		double[] sm_sY = new double[sm_NumNodes]; // Number of instances under each node
		double[] weights = new double[nC]; // Weights of instances
		double totWeights = 0; // Total weights of the instances
		m_Par = new double[nR + 1][nK]; // Optimized parameter values of the leaf nodes
		sm_Par = new double[nR+1][sm_NumNodes]; // Optimized parameter values of all the nodes

		for (int i = 0; i < nC; i++) {
			// initialize X[][]
			Instance current = train.instance(i);
			Y[i] = (int) current.classValue(); // Class value starts from 0
			weights[i] = current.weight(); // Dealing with weights
			totWeights += weights[i];

			m_Data[i][0] = 1;
			int j = 1;
			for (int k = 0; k <= nR; k++) {
				if (k != m_ClassIndex) {
					double x = current.value(k);
					m_Data[i][j] = x;
					xMean[j] += weights[i] * x;
					xSD[j] += weights[i] * x * x;
					j++;
				}
			}
			// Class count
			sY[Y[i]]++;
		}

		if ((totWeights <= 1) && (nC > 1)) {
			throw new Exception(
					"Sum of weights of instances less than 1, please reweight!");
		}
		xMean[0] = 0;
		xSD[0] = 1;
		for (int j = 1; j <= nR; j++) {
			xMean[j] = xMean[j] / totWeights;
			if (totWeights > 1) {
				xSD[j] = Math.sqrt(Math.abs(xSD[j] - totWeights * xMean[j] * xMean[j])/ (totWeights - 1));
			} else {
				xSD[j] = 0;
			}
		}

		// Normalise input data
		for (int i = 0; i < nC; i++) {
			for (int j = 0; j <= nR; j++) {
				if (xSD[j] != 0) {
					m_Data[i][j] = (m_Data[i][j] - xMean[j]) / xSD[j];
				}
			}
		}

		double x[] = new double[(nR + 1) * nK];
		double[][] b = new double[2][x.length]; // Boundary constraints, N/A here

		// Fill sm_sY
		// First fill the leaf node
		for(Node n : sm_ClassStructure.leafs){
			sm_sY[n.nodeIndex] = sY[n.nodeIndex];
		}
		for(int l=1; l<sm_ClassStructure.numLayers; l++){
			for(Node n :sm_ClassStructure.layers.get(l)){
				for(Integer cind : n.children){
					sm_sY[n.nodeIndex] = sm_sY[n.nodeIndex]+sm_sY[cind];
				}
			}
		}

		// Initialize
		for (int p = 0; p < nK; p++) {
			int offset = p * (nR + 1);
			x[offset] = Math.log(sY[p] + 1.0) - Math.log(sY[nK-1] + 1.0); // Null model
			b[0][offset] = Double.NaN;
			b[1][offset] = Double.NaN;
			for (int q = 1; q <= nR; q++) {
				x[offset + q] = 0.0;
				b[0][offset + q] = Double.NaN;
				b[1][offset + q] = Double.NaN;
			}
		}

		double[] sm_x = new double[(nR+1)*sm_NumNodes];
		// Initialize the sm_x parameters
		for(int p=0; p<sm_NumNodes; p++) {
			int offset = p*(nR+1);
			sm_x[offset] = Math.log(sm_sY[p] + 1.0) - Math.log(sm_sY[nK-1] + 1.0); // Null model
			for (int q = 1; q <= nR; q++) {
				sm_x[offset+q] = 0.0;
			}
		}

		Optimizer opt = null;

		if(m_OptimizationType.equals("L1"))
			opt = new LOne(x,sm_x,m_Data);
		
		opt.setRidge(m_Ridge);
		if(opt instanceof LOne){
			((LOne) opt).setADMMmaxItrs(m_ADMM_MaxIts);
		}
		opt.setBGFSmaxItrs(m_BGFS_MaxIts);
		opt.setSeqUnwinderMaxIts(m_SeqUnwinder_MaxIts);
		opt.setClassStructure(sm_ClassStructure);
		opt.setClsMembership(Y);
		opt.setInstanceWeights(weights);
		opt.setNumClasses(m_NumClasses);
		opt.setNumPredictors(m_NumPredictors);
		if(opt instanceof LOne){
			((LOne) opt).setPho(m_ADMM_pho);
			((LOne) opt).set_numThreads(m_ADMM_numThreads);
			((LOne) opt).initUandZ();
		}
		opt.setDebugMode(sm_Debug);
		opt.execute();
		m_Data = null;

		x=opt.getX();
		sm_x= opt.getsmX();
		opt.clearOptimizer();

		//Now update the m_Par and sm_Par with the learned parameters
		for (int i = 0; i < nK; i++) {
			m_Par[0][i] = x[i * (nR + 1)];
			for (int j = 1; j <= nR; j++) {
				m_Par[j][i] = x[i * (nR + 1) + j];
				if (xSD[j] != 0) {
					m_Par[j][i] /= xSD[j];
					m_Par[0][i] -= m_Par[j][i] * xMean[j];
				}
			}
		}

		// Don't need data matrix anymore
		m_Data = null;


		//Now update the sm_Par with the learned parameters
		for(Node n : sm_ClassStructure.allNodes.values()){
			sm_Par[0][n.nodeIndex] = sm_x[n.nodeIndex*(nR+1)];
			for(int j=1; j<=nR; j++){
				sm_Par[j][n.nodeIndex] = sm_x[n.nodeIndex*(nR+1)+j];
				if(xSD[j] != 0 ){
					sm_Par[j][n.nodeIndex] /=xSD[j];
					sm_Par[0][n.nodeIndex] -= sm_Par[j][n.nodeIndex]*xMean[j];
				}
			}
		}
	}


	@Override
	public String toPMML(Instances train) {
		return LogisticProducerHelper.toPMML(train, m_structure, m_Par,
				m_NumClasses);
	}



	public String toString(){
		StringBuffer temp = new StringBuffer();
		String result = "";
		temp.append("Structured MutiLogistic Regression with ridge parameter of " + m_Ridge);
		if (sm_Par == null) {
			return result + ": No model built yet.";
		}

		// find longest attribute name
		int attLength = 0;
		for (int i = 0; i < m_structure.numAttributes(); i++) {
			if (i != m_structure.classIndex() && m_structure.attribute(i).name().length() > attLength) {
				attLength = m_structure.attribute(i).name().length();
			}
		}

		if ("Intercept".length() > attLength) {
			attLength = "Intercept".length();
		}

		if ("Variable".length() > attLength) {
			attLength = "Variable".length();
		}
		attLength += 2;

		int colWidth = 0;
		// check length of class names
		for (Node n : sm_ClassStructure.allNodes.values()) {
			if (n.nodeName.length() > colWidth) {
				colWidth = n.nodeName.length();
			}
		}

		// check against coefficients and odds ratios
		for (int j = 1; j <= m_NumPredictors; j++) {
			for (Node n : sm_ClassStructure.allNodes.values() ) {
				if (Utils.doubleToString(sm_Par[j][n.nodeIndex], 12, 4).trim().length() > colWidth) {
					colWidth = Utils.doubleToString(sm_Par[j][n.nodeIndex], 12, 4).trim().length();
				}
				double ORc = Math.exp(sm_Par[j][n.nodeIndex]);
				String t = " "
						+ ((ORc > 1e10) ? "" + ORc : Utils.doubleToString(ORc, 12, 4));
				if (t.trim().length() > colWidth) {
					colWidth = t.trim().length();
				}
			}
		}

		if ("Class".length() > colWidth) {
			colWidth = "Class".length();
		}
		colWidth += 2;

		temp.append("\nCoefficients...\n");
		temp.append(Utils.padLeft(" ", attLength)
				+ Utils.padLeft("Class", colWidth) + "\n");
		temp.append(Utils.padRight("Variable", attLength));
		
		for(int n=0; n<sm_ClassStructure.allNodes.size(); n++){
			String className = sm_ClassStructure.allNodes.get(n).nodeName;
			temp.append(Utils.padLeft(className, colWidth));
		}

		//for (Node n : sm_ClassStructure.allNodes.values()) {
		//	String className = n.nodeName;
		//	temp.append(Utils.padLeft(className, colWidth));
		//}
		temp.append("\n");
		int separatorL = attLength + ((sm_NumNodes) * colWidth);
		for (int i = 0; i < separatorL; i++) {
			temp.append("=");
		}
		temp.append("\n");

		int j = 1;
		for (int i = 0; i < m_structure.numAttributes(); i++) {
			if (i != m_structure.classIndex()) {
				temp.append(Utils.padRight(m_structure.attribute(i).name(), attLength));
				for(int n=0; n<sm_ClassStructure.allNodes.size(); n++){
					temp.append(Utils.padLeft(Utils.doubleToString(sm_Par[j][sm_ClassStructure.allNodes.get(n).nodeIndex], 12, 4).trim(), colWidth));
				}
				//for (Node n : sm_ClassStructure.allNodes.values()) {
				//	temp.append(Utils.padLeft(Utils.doubleToString(sm_Par[j][n.nodeIndex], 12, 4).trim(), colWidth));
				//}
				temp.append("\n");
				j++;
			}
		}

		temp.append(Utils.padRight("Intercept", attLength));
		for(int n=0; n<sm_ClassStructure.allNodes.size(); n++){
			temp.append(Utils.padLeft(Utils.doubleToString(sm_Par[0][sm_ClassStructure.allNodes.get(n).nodeIndex], 10, 4).trim(), colWidth));
		}
		//for (Node n : sm_ClassStructure.allNodes.values()) {
		//	temp.append(Utils.padLeft(Utils.doubleToString(sm_Par[0][n.nodeIndex], 10, 4).trim(), colWidth));
		//}
		temp.append("\n");

		temp.append("\n\nOdds Ratios...\n");
		temp.append(Utils.padLeft(" ", attLength)
				+ Utils.padLeft("Class", colWidth) + "\n");
		temp.append(Utils.padRight("Variable", attLength));
		
		for(int n=0; n<sm_ClassStructure.allNodes.size(); n++){
			String className = sm_ClassStructure.allNodes.get(n).nodeName;
			temp.append(Utils.padLeft(className, colWidth));
		}
		//for (Node n: sm_ClassStructure.allNodes.values()) {
		//	String className = n.nodeName;
		//	temp.append(Utils.padLeft(className, colWidth));
		//}
		temp.append("\n");
		for (int i = 0; i < separatorL; i++) {
			temp.append("=");
		}
		temp.append("\n");

		j = 1;
		for (int i = 0; i < m_structure.numAttributes(); i++) {
			if (i != m_structure.classIndex()) {
				temp.append(Utils.padRight(m_structure.attribute(i).name(), attLength));
				
				for(int n=0; n<sm_ClassStructure.allNodes.size(); n++){
					double ORc = Math.exp(sm_Par[j][sm_ClassStructure.allNodes.get(n).nodeIndex]);
					String ORs = " "+ ((ORc > 1e10) ? "" + ORc : Utils.doubleToString(ORc, 12, 4));
					temp.append(Utils.padLeft(ORs.trim(), colWidth));
				}
				
				//for (Node n : sm_ClassStructure.allNodes.values()) {
				//	double ORc = Math.exp(sm_Par[j][n.nodeIndex]);
				//	String ORs = " "+ ((ORc > 1e10) ? "" + ORc : Utils.doubleToString(ORc, 12, 4));
				//	temp.append(Utils.padLeft(ORs.trim(), colWidth));
				//}
				temp.append("\n");
				j++;
			}
		}

		return temp.toString();

	}



	protected class ClassRelationStructure implements Serializable{

		/**
		 * 
		 */
		private static final long serialVersionUID = -4392515471314248045L;
		protected List<Node> leafs = new ArrayList<Node>();
		protected List<List<Node>> layers = new ArrayList<List<Node>>();
		protected Map<Integer,Node> allNodes = new HashMap<Integer,Node>();
		protected int numLayers;





		public ClassRelationStructure(File structureFile,int nLayers) throws IOException {
			BufferedReader br = new BufferedReader(new FileReader(structureFile));
			String line;
			numLayers=nLayers;


			for(int i=0; i<numLayers; i++){
				layers.add(new ArrayList<Node>());
			}

			while ((line = br.readLine()) != null) {
				if(!line.startsWith("#")){
					String[] pieces = line.split("\t");
					if(pieces.length <6){System.err.println("Incorrect class structure format!!!"); System.exit(1);}
					Node currNode = new Node(Integer.parseInt(pieces[0]), pieces[1],Integer.parseInt(pieces[2]) == 1? true : false);
					// Add parents indexes
					//if(!currNode.nodeName.contains("Root")){
					if(currNode.isLeaf && !pieces[3].equals("-")){
						String[] pInds = pieces[3].split(",");
						for(String s : pInds){
							currNode.addParent(Integer.parseInt(s));
						}
					}
					// Add children indexes
					if(!currNode.isLeaf && !pieces[4].equals("-")){
						String[] cInds = pieces[4].split(",");
						for(String s : cInds){
							currNode.addChild(Integer.parseInt(s));
						}
					}
					allNodes.put(currNode.nodeIndex, currNode);
					if(currNode.isLeaf)
						leafs.add(currNode);
					layers.get(Integer.parseInt(pieces[5])).add(currNode);	
				}
			}br.close();

			// Print the loaded class structure
			if(getSmDebug()){
				StringBuilder structure = new StringBuilder();
				for(int l=0; l<layers.size(); l++){

					for(Node n : layers.get(l)){
						structure.append(n.nodeName); structure.append("\t");
					}
					structure.deleteCharAt(structure.length()-1);
					structure.append("\n");
				}

				System.err.println(structure.toString());
			}
		}


		protected class Node implements Serializable,Comparable<Node>{

			/**
			 * 
			 */
			private static final long serialVersionUID = 9057711985309713653L;
			protected int nodeIndex;
			protected String nodeName;
			protected boolean isLeaf=false;
			protected int level;
			protected List<Integer> parents = new ArrayList<Integer>();
			protected List<Integer> children = new ArrayList<Integer>();

			protected Node(int nInd, String nName, boolean leaf) {
				nodeIndex = nInd;
				nodeName = nName;
				isLeaf = leaf;
			}

			protected void addParent(Integer pInd){
				parents.add(pInd);
			}
			protected void addChild(Integer cInd){
				children.add(cInd);
			}

			@Override
			public int compareTo(Node o) {
				if(o.nodeIndex == nodeIndex){return 0;}
				else if(nodeIndex > o.nodeIndex){return 1;}
				else{return -1;}
			}

		}
	}

	@Override
	public Enumeration<Option> listOptions() {
		Vector<Option> newVector = new Vector<Option>(7);

		newVector.addElement(new Option(
				"\tUse conjugate gradient descent rather than BFGS updates.", "C", 0,
				"-C"));
		newVector.addElement(new Option("\tSet the ridge in the log-likelihood.",
				"R", 1, "-R <ridge>"));
		newVector.addElement(new Option("\t Pho dual co-efficient for ADMM.",
				"PHO", 1, "-PHO <default is 0.0001>"));
		newVector.addElement(new Option("\tSet the maximum number of iterations for the BGFS x-update"
				+ " (default -1, until convergence).", "M", 1, "-M <number>"));
		newVector.addElement(new Option("\tSet the maximum number of iterations for the ADMM method"
				+ " (default 1000).", "A", 1, "-A <number>"));
		newVector.addElement(new Option("\tSet the maximum number of iterations for SeqUwinder"
				+ " (default 10).", "S", 1, "-S <number>"));
		newVector.addElement(new Option("\tSet the number of threads to run ADMM part in SeqUnwinder"
				+ " (default 5).", "threads", 1, "-threads <number>"));
		newVector.addElement(new Option("\tProvide the class structure file for the multiclasses","CLS",0,"-CLS <File name>"));
		newVector.addElement(new Option("\tProvide the number of layers in the class structur","NL",0,"-NL <Num of Layers>"));
		newVector.addElement(new Option("\tL1 or L2 norm Logit model","TY",0,"-TY <L1 or L2>"));
		newVector.addElement(new Option("\tFlag to run in debug mode","DEBUG",0,"-DEBUG"));
		newVector.addAll(Collections.list(super.listOptions()));

		return newVector.elements();
	}



	@Override
	public void setOptions(String[] options) throws Exception {

		setSmDebug(Utils.getFlag("DEBUG", options));

		String ridgeString = Utils.getOption('R', options);
		if (ridgeString.length() != 0) {
			m_Ridge = Double.parseDouble(ridgeString);
		} else {
			m_Ridge = 10.0;
		}

		String threadString = Utils.getOption("threads", options);
		if(threadString.length() != 0){
			m_ADMM_numThreads = Integer.parseInt(threadString);
		}

		String maxItsString = Utils.getOption('M', options);
		if (maxItsString.length() != 0) {
			m_BGFS_MaxIts = Integer.parseInt(maxItsString);
		} else {
			m_BGFS_MaxIts = -1;
		}

		String AmaxItsString = Utils.getOption('A', options);
		if (AmaxItsString.length() != 0) {
			m_ADMM_MaxIts = Integer.parseInt(AmaxItsString);
		} else {
			m_ADMM_MaxIts = 1000;
		}

		String PhoString = Utils.getOption("PHO", options);
		if(PhoString.length() !=0){
			m_ADMM_pho = Double.parseDouble(PhoString);
		}else{
			m_ADMM_pho = 0.001;
		}

		String SmaxItsString = Utils.getOption('S', options);
		if (SmaxItsString.length() != 0) {
			m_SeqUnwinder_MaxIts = Integer.parseInt(SmaxItsString);
		} else {
			m_SeqUnwinder_MaxIts = 10;
		}

		String numLayerString = Utils.getOption("NL", options);
		if(numLayerString == ""){
			System.err.println("Please provide number of layers!!");
			System.exit(1);
		}
		sm_NumLayers = Integer.parseInt(numLayerString);

		String Type = Utils.getOption("TY", options);
		if (Type.length() != 0) {
			m_OptimizationType = Type;
		}else{
			m_OptimizationType = "L1";
		}

		String classStructureFilename = Utils.getOption("CLS", options);
		if(classStructureFilename == ""){
			System.err.println("Please provide the class structure file!!");
			System.exit(1);
		}else{
			File f = new File(classStructureFilename);
			sm_ClassStructure = new ClassRelationStructure(f,sm_NumLayers);
			sm_NumNodes = sm_ClassStructure.allNodes.size();
		}


		super.setOptions(options);

		Utils.checkForRemainingOptions(options);

	}



	@Override
	public String[] getOptions() {
		Vector<String> options = new Vector<String>();

		options.add("-R");
		options.add("" + m_Ridge);
		options.add("-PHO");
		options.add("" + m_ADMM_pho);
		options.add("-M");
		options.add("" + m_BGFS_MaxIts);
		options.add("-A");
		options.add("" + m_ADMM_MaxIts);
		options.add("-S");
		options.add("" + m_SeqUnwinder_MaxIts);
		options.add("-CLS");
		options.add(""+ " ");
		options.add("-NL");
		options.add("" + 3);
		options.add("-TY");
		options.add("" + m_OptimizationType);
		options.add("-DEBUG");
		options.add("" + sm_Debug);
		options.add("-threads");
		options.add("" + m_ADMM_numThreads);

		Collections.addAll(options, super.getOptions());

		return options.toArray(new String[0]);
	}


	/**
	 * Main method for testing this class.
	 * 
	 * @param argv should contain the command line arguments to the scheme (see
	 *          Evaluation)
	 */
	public static void main(String[] argv) {
		runClassifier(new Classifier(), argv);
	}
}
