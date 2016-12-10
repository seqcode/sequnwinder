package org.seqcode.projects.sequnwinder.framework;

import org.seqcode.projects.sequnwinder.framework.Classifier.ClassRelationStructure;

public abstract class Optimizer {
	
	public abstract void setBGFSmaxItrs(int m);
	
	public abstract void setSeqUnwinderMaxIts(int m);
	public abstract void setInstanceWeights(double[] w);
	public abstract void setClsMembership(int[] c);
	public abstract void setNumPredictors(int p);
	public abstract void setNumClasses(int c);
	public abstract void setClassStructure(ClassRelationStructure rel);
	public abstract void setNumNodes(int n);
	public abstract void setRidge(double r);
	public abstract void setDebugMode(boolean debug);
	
	public abstract void execute() throws Exception;
	
	public abstract double[] getX();
	public abstract double[] getsmX();
	
	public abstract void clearOptimizer();
	

}
