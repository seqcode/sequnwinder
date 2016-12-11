package org.seqcode.projects.sequnwinder.clusterkmerprofile;

import java.util.Set;

import org.seqcode.ml.clustering.Cluster;
import org.seqcode.ml.clustering.ClusterRepresentative;


/**
 * KmerProfileAvgDistRep: Cluster representative of hills in K-mer space
 * @author akshaykakumanu
 * @version	%I%, %G%
 */
public class KmerProfileAvgDistRep implements ClusterRepresentative<int[]> {
	
	KmerProfileEucDistComparator comp;
	
	public  KmerProfileAvgDistRep(KmerProfileEucDistComparator c) {
		comp = c;
	}
	
	@Override
	public int[] getRepresentative(Cluster<int[]> c) {
		Set<int[]> profiles = c.getElements();
		int[] bestProfile = null;
		double bestDist = Double.MAX_VALUE;
		for(int[] i : profiles){
			double sum=0;
			for(int[] j : profiles){
				sum = sum + comp.evaluate(i, j);
			}
			sum = sum/profiles.size();
			if(sum < bestDist){
				bestDist = sum;
				bestProfile = i;
			}
		}
		
		return bestProfile;
	}

}
