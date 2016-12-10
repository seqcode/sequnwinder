package org.seqcode.projects.sequnwinder.utils;

public class CGScoreProfile {
	private int[] cg;
	private double[] pCG;
		 
	public CGScoreProfile(int[] c, double[] pc) {
		cg=c;
		pCG= pc;
		
	}
	
	public int[] getCGprofile(){
		return cg;
	}
	
	public double[] getCGpercProfile(){
		return pCG;
	}
	
	public int getCGcount(){
		int count=0;
		for(int i=0;i<cg.length;i++){
			count =count+cg[i]; 
		}
		return count;
	}
		 
}


