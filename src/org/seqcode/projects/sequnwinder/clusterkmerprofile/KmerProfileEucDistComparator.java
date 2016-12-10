package org.seqcode.projects.sequnwinder.clusterkmerprofile;

import org.seqcode.ml.clustering.PairwiseElementMetric;

public class KmerProfileEucDistComparator implements PairwiseElementMetric<int[]> {

	@Override
	public double evaluate(int[] e1, int[] e2) {
		double euclDistace=0;
		for(int i=0; i<e1.length; i++){
			euclDistace = euclDistace + Math.pow(e1[i]-e2[i], 2);
		}
		
		euclDistace = Math.sqrt(euclDistace);
		
		return euclDistace;
	}

}
