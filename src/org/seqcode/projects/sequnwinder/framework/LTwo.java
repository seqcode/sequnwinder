package org.seqcode.projects.sequnwinder.framework;

import org.seqcode.projects.sequnwinder.framework.Classifier.ClassRelationStructure;
import org.seqcode.projects.sequnwinder.framework.Classifier.ClassRelationStructure.Node;

import weka.core.Optimization;
import weka.core.RevisionUtils;
import weka.core.Utils;

public class LTwo extends Optimizer{

	// Fixed SeqUnwinder parameters
	
	/** Tolerence for internal Nodes convergence */
	public final double NODES_tol = 1E-1;
	
	// BGFS parameters 
	
	/** The maximum number of allowed iterations for the BGFS algorithm */
	public int BGFS_maxIts=-1;
	
	// SeqUnwinder parameters
	/** The total number of predictors/features in the model (does not include the intercept term) */
	protected int numPredictors;
	/** Total number of classes to be predicted */
	protected int numClasses;
	/** Relationships between the different nodes in SeqUnwinder */
	protected ClassRelationStructure classStructure;
	/** Total number of nodes (internal nodes and classes) in SeqUnwinder */
	protected int numNodes;
	/** The L1 regularization parameter in the SeqUnwinder model */
	protected double regularization;
	/** Maximum number of iterations to update internal nodes. For small number of node levels (usually 3 to 4), we might need very few iterations*/
	protected int NODES_maxItr=10;
	
	/** Current feature weights for all the nodes in the SeqUnwinder (also contains the intercept term). Dimension :- (numPredictors+1)*numNodes */
	public double[] sm_x;
	/** Current feature weights for all the leaf nodes (classes) in the SeqUnwinder (also contains the intercept term). Dimension :- (numPredictors+1)*numClasses */
	public double[] x;
	/** Boundry conditions needed for the the newton's method. However, null in this case */
	public double[][] b;
	
	// SeqUnwinder training data
	
	/** Training data (Instances)*/
	public double[][] data;
	/** Weights of the instances */
	protected double[] weights;
	/** Instance class membership */
	protected int[] cls;
		
	// Misc
	
	/** Boolean variable indicating debug mode */
	protected boolean sm_Debug;
	
	//Settors
	public void setBGFSmaxItrs(int m){BGFS_maxIts=m;}
	public void setSeqUnwinderMaxIts(int m){NODES_maxItr = m;}
	public void setInstanceWeights(double[] w){weights=w;}
	public void setClsMembership(int[] c){cls=c;}
	public void setNumPredictors(int p){numPredictors=p;}
	public void setNumClasses(int c){numClasses = c;}
	public void setClassStructure(ClassRelationStructure rel){classStructure = rel; setNumNodes(rel.allNodes.size());}
	public void setNumNodes(int n){numNodes = n;}
	public void setRidge(double r){regularization = r;}
	public void setDebugMode(boolean debug){sm_Debug =debug;}
	
	//gettors
	public double[] getX(){return x;}
	public double[] getsmX(){return sm_x;}
	
	public LTwo(double[] xinit, double[] sm_xinit, double[][] d, double[][] bc) {
		x = xinit;
		sm_x = sm_xinit;
		data=d;
		b=bc;
	}
	@Override
	public void execute() throws Exception {
		for(int it=0; it<NODES_maxItr; it++){
			// First, run admm on leaf nodes
			leafUpdate();
			
			// Now update the internal nodes
			//updateInternalNodes();
			
			// Check Convergence
			//boolean converged = true;
			//for(Node n : classStructure.allNodes.values()){
			//	if( getL2NormX(n.nodeIndex) > Math.sqrt(numPredictors)*NODES_tol){
			//		converged=false;
			//		break;
			//	}
			//}
			
			//if(converged)
			//	break;
		}
		
	}
	
	public void leafUpdate() throws Exception{
		// First, create and opt object for the BGFS algorithm
		OptObject oO = new OptObject();
		Optimization opt = new OptEng(oO);
		
		if(sm_Debug){
			System.err.println("Today X values before the Newtown method");
			System.err.println(x[20]);
			System.err.println(x[numPredictors+1+20]);
			System.err.println(x[2*(numPredictors+1)+20]);
		}
		
		 if ( BGFS_maxIts == -1) { // Search until convergence
			 x = opt.findArgmin(x, b);
			 while (x == null) {
				 x = opt.getVarbValues();
				 x = opt.findArgmin(x, b);
			 }
		 } else {
			 opt.setMaxIteration(BGFS_maxIts);
			 x = opt.findArgmin(x, b);
			 if (x == null) {
				 x = opt.getVarbValues();
			 }
		 }
		 
		 if(sm_Debug){
				System.err.println("Today X values after the Newtown method");
				System.err.println(x[20]);
				System.err.println(x[numPredictors+1+20]);
				System.err.println(x[2*(numPredictors+1)+20]);
		 }
		 
		// Now copy the leaf node weights (i.e x) to sm_x
		 for(Node n: classStructure.leafs){
			 int dim = numPredictors+1;
			 int nOffset = n.nodeIndex*dim;
			 for(int w=0; w<dim; w++){
				 sm_x[nOffset+w] = x[nOffset+w];
			 }
		 }
	}
	
	
	/**
	 * This class implements two things:-
	 * It calculates the gradient (The BGFS method will need this)
	 * It calculates the overall objective function. (The BGFS method will need this)
	 * @author akshaykakumanu
	 *
	 */
	public class OptObject {
		
		public double[] evaluateGradient(double[] x){
			
			double[] grad = new double[x.length];
			int dim = numPredictors + 1; // Number of variables per class

			for (int i = 0; i < cls.length; i++) { // ith instance
				double[] num = new double[numClasses]; // numerator of
		                                                     // [-log(1+sum(exp))]'
		        int index;
		        for (int offset = 0; offset < numClasses; offset++) { // Which
		                                                                    // part of 
		        	double exp = 0.0;
		        	index = offset * dim;
		        	for (int j = 0; j < dim; j++) {
		        		exp += data[i][j]*x[index + j];
		        	}
		        	num[offset] = exp;
		        }

		        double max = num[Utils.maxIndex(num)];
		       
		        double denom=0.0;
		        for (int offset = 0; offset < numClasses; offset++) {
		        	num[offset] = Math.exp(num[offset] - max);
		        	denom += num[offset];
		        }
		        Utils.normalize(num, denom);

		        // Update denominator of the gradient of -log(Posterior)
		        double firstTerm;
		        for (int offset = 0; offset < numClasses; offset++) { // Which
		                                                                    // part of x
		        	index = offset * dim;
		        	firstTerm = weights[i] * num[offset];
		        	for (int q = 0; q < dim; q++) {
		        		grad[index + q] += firstTerm * data[i][q];
		        	}
		        }

		        for (int p = 0; p < dim; p++) {
		            grad[cls[i] * dim + p] -= weights[i] * data[i][p];
		        }
		        
			}
		      

			for(Node n : classStructure.leafs){
				int nOffset = n.nodeIndex*dim;
				if(n.parents.size() > 0){
					for(int pid : n.parents){
						int pOffset = pid*dim;
						for(int w=0; w<dim; w++){
							grad[nOffset+w] += regularization*(x[nOffset+w]-sm_x[pOffset+w]);
						}
					}
				}else{
					for(int w=0; w<dim; w++){
						grad[nOffset+w] += regularization*(x[nOffset+w]);
					}
				}
			}	
		      
			return grad;
					
		}
		
		public double objectiveFunction(double[] x){
			double nll=0.0;
			int dim = numPredictors+1;
			
			for (int i = 0; i < cls.length; i++) { // ith instance
				double[] exp = new double[numClasses];
				int index;
				for (int offset = 0; offset < numClasses; offset++) {
					index = offset * dim;
					for (int j = 0; j < dim; j++) {
						exp[offset] += data[i][j] * x[index + j];
					}
				}
				double max = exp[Utils.maxIndex(exp)];
		        double denom = 0;
		        double num = exp[cls[i]] - max;
		        
		        for (int offset = 0; offset < numClasses; offset++) {
		        	denom += Math.exp(exp[offset] - max);
		        }

		        nll -= weights[i] * (num - Math.log(denom)); // Weighted NLL
			}
			
			for(Node n : classStructure.leafs){
				int nOffset = n.nodeIndex*dim;
				if(n.parents.size() >0){
					for(int pid : n.parents){
						int pOffset = pid*dim;
						for(int w=0; w<dim; w++){
							nll += (regularization/2)*(x[nOffset+w]-sm_x[pOffset+w])*(x[nOffset+w]-sm_x[pOffset+w]);
						}
					}
				}else{
					for(int w=0; w<dim; w++){
						nll += (regularization/2)*(x[nOffset+w]*x[nOffset+w]);
					}
				}
			}
		
			if(sm_Debug){
				//System.err.println("Objective:"+nll);
			}
			
			return nll;
		}
		
	}
	
	public class OptEng extends Optimization {
		
		OptObject m_oO = null;
		
		public OptEng(OptObject oO){
			m_oO = oO;
		}
		 
		@Override
		public String getRevision() {
			return RevisionUtils.extract("$Revision: 11247 $");
		}

		@Override
		protected double objectiveFunction(double[] x) throws Exception {
			return m_oO.objectiveFunction(x);
		}

		@Override
		protected double[] evaluateGradient(double[] x) throws Exception {
			return m_oO.evaluateGradient(x);
		}
	}
	
	// To clear memory
	public void clearOptimizer(){
		data=null;
		sm_x=null;
		x=null;
	}
	
}
