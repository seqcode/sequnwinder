package org.seqcode.projects.sequnwinder.clusterkmerprofile;

import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.Vector;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.ml.clustering.Cluster;
import org.seqcode.ml.clustering.ClusterRepresentative;
import org.seqcode.ml.clustering.DefaultCluster;
import org.seqcode.ml.clustering.PairwiseElementMetric;
import org.seqcode.ml.clustering.kmeans.KMeansClustering;
import org.seqcode.ml.clustering.vectorcluster.DefaultVectorClusterElement;
import org.seqcode.ml.clustering.vectorcluster.EuclideanDistance;
import org.seqcode.ml.clustering.vectorcluster.Mean;
import org.seqcode.ml.clustering.vectorcluster.VectorClusterElement;
import org.tc33.jheatchart.HeatChart;


/**
 * ClusterProfiles: Takes a list of data-points (hills in this case) and removes low penetrance K-mers to generate a
 *  sparse representation of hills in the k-mer space. Next, does a k-means clustering with Euclidean Distance metric.
 * 
 * @author akshaykakumanu
 * @version	%I%, %G%
 */

public class ClusterProfiles {
	private PairwiseElementMetric<VectorClusterElement> metric = new EuclideanDistance<VectorClusterElement>();
	private int K;
	private int numClusItrs;
	private ArrayList<VectorClusterElement> sparse_profiles = new ArrayList<VectorClusterElement>();
	private ArrayList<String> sparse_colnames = new ArrayList<String>();
	private File outdir;
	private HashMap<Integer,String> indToLocation;
	private HashMap<Integer,Double> intToLogitScore;
	private int minK=4;
	private int maxK=5;
	
	private int numBootstraps = 50;
	private double bsSizeFraction = 0.3;
	private int minC = 2;
	private int maxC = 5;

	
	// Minimum penetrance of a K-mer in a cluster to be considered
	//public final double minKmerProp_clus = 0.2;
	public final double minKmerProp_global = 0.04;
	
	// Heatmap options
	public final int W_MARGIN=80;
	public final int H_MARGIN=60;
	public final int W=800;
	public final int H=800;
	
	//Gettors
	public int getKmerBaseInd(int len){
		int baseInd = 0;
		for(int k=minK; k<len; k++){
			baseInd += (int)Math.pow(4, k);
		}
		return baseInd;
	}
	public String getKmerName(int ind){
		int currKmerLen = 0;
		ArrayList<Integer> baseinds = new ArrayList<Integer>();
		for(int k=minK; k <= maxK; k++){
			baseinds.add(getKmerBaseInd(k));
		}
		int search = Collections.binarySearch(baseinds, ind);
		
		currKmerLen = search >= 0 ? minK + search : -1*(search + 1) - 1 + minK;
		String kmerName = RegionFileUtilities.int2seq(ind- getKmerBaseInd(currKmerLen), currKmerLen);
		
		return kmerName;
	}
	
	//Settors
	public void setNumClusters(int nc){K=nc;}
	public void setProfileInds(HashMap<Integer,String> plfsinds){indToLocation=plfsinds;}
	public void setProfileScores(HashMap<Integer,Double> plfscore){intToLogitScore = plfscore;}
	public void setKmerModLenMin(int k){minK = k;}
	public void setKmerModLenMax(int k){maxK = k;}
	public void setOutdir(File f){outdir = f;}
	
	
	public void setSparcedProfiles(ArrayList<int[]> pfs){
		
		// First find the penetrance of each feature
		double[] feature_penetrance = new double[pfs.get(0).length];
		for(int[] p : pfs){
			for(int i=0; i<p.length; i++){
				if(p[i] > 0)
					feature_penetrance[i]++;
			}
		}
		
		// Find the colnames of the once that will be kept and store the colnames in a list
		for(int i=0; i<feature_penetrance.length; i++){
			feature_penetrance[i] = feature_penetrance[i]/pfs.size();
			if(feature_penetrance[i] > minKmerProp_global){
				sparse_colnames.add(getKmerName(i));
			}
		}

		// Now make the sparse profiles
		int sparce_length = sparse_colnames.size();
		//System.err.println(sparce_length);

		for(int[] p : pfs){
			double[] sparce_p = new double[sparce_length];
			int count=0;
			for(int i=0; i<p.length; i++){
				if(feature_penetrance[i] > minKmerProp_global){
					sparce_p[count] = (double)p[i];
					count++;
				}		
			}
			DefaultVectorClusterElement v = new DefaultVectorClusterElement(sparce_p);
			sparse_profiles.add(v);
		}
	}
	
	/**
	 * The method that should be executed after initiating the class object
	 * @throws IOException 
	 */
	public List<Integer> execute() throws IOException{
		for(int c = minC; c<=maxC; c++){
			ArrayList<VectorClusterElement> currBS = generateBootStrapSample();
			double currSH = getSilhouette(c,currBS);
			StringBuilder sb = new StringBuilder();
			sb.append(c);sb.append("\t");sb.append(currSH);
			System.out.println(sb.toString());
		}
		
		
		
		
		//Print the clusters
		List<Integer> clusAssignment = null;
		//List<Integer> clusAssignment = writeClusters(clustermeans);
		
		return clusAssignment;
	}
	
	public boolean added(Vector<VectorClusterElement> added, VectorClusterElement curr){
		boolean ret = false;
		for(VectorClusterElement a : added){
			boolean match = true;
			for(int i=0; i<a.dimension(); i++){
				if(curr.getValue(i) != a.getValue(i))
					match = false;
			}
			if(match){
				ret=true;
				break;
			}
		}
		return ret;
		
	}
	
	
	// Slave methods
	private List<Integer> writeClusters(Vector<VectorClusterElement> clusMeans) throws IOException{
		List<Integer> clusterAssignment = new ArrayList<Integer>();
		File clusout = new File(outdir.getAbsolutePath()+File.separator+"ClusterAssignment.list");
		FileWriter ow = new FileWriter(clusout);
		BufferedWriter bw = new BufferedWriter(ow);
		for(int p=0; p<sparse_profiles.size(); p++){
			int memebership = getClusterAssignment(sparse_profiles.get(p),clusMeans);
			clusterAssignment.add(memebership);
			bw.write(indToLocation.get(p)+"\t"+Integer.toString(memebership)+"\t"+Double.toString(intToLogitScore.get(p))+"\n");
		}
		bw.close();
		return clusterAssignment;
	}
	
	private int getClusterAssignment(VectorClusterElement pfl, Vector<VectorClusterElement> clusMeans){
		int minCluster = -1;
        double minDist = 0.0;
        for(int i = 0; i < clusMeans.size(); i++) { 
            double clustDist = metric.evaluate(pfl, clusMeans.get(i));
            if(minCluster == -1 || clustDist < minDist) { 
                minDist = clustDist;
                minCluster = i;
            }
        }
        return minCluster;
	}

	
	private ArrayList<VectorClusterElement> generateBootStrapSample(){
		ArrayList<VectorClusterElement> ret = new ArrayList<VectorClusterElement>();
		Random rand = new Random();
		int maxInd = sparse_profiles.size();
		int count = 0;
		while(count <(int)(bsSizeFraction*maxInd)){
			ret.add(sparse_profiles.get(rand.nextInt(maxInd)));
			count++;
		}
		return ret;
	}
	
	private double getSilhouette(int nC, ArrayList<VectorClusterElement> data){
		// Initialize
		ClusterRepresentative<VectorClusterElement> crep = new Mean();
		//Random starts
		Random generator = new Random();
		Vector<VectorClusterElement> starts = new Vector<VectorClusterElement>();
		for(int s=0; s<nC; s++){
			boolean found = false;
			while(!found){
				int r = generator.nextInt(data.size());
				if(!added(starts,data.get(r))){
					starts.add(data.get(r));
					found = true;
				}
			}

		}

		//Initialize clustering
		KMeansClustering<VectorClusterElement> kmc = new KMeansClustering<VectorClusterElement>(metric, crep, starts);
		kmc.setIterations(nC);
		kmc.clusterElements(data,0.01);
		return kmc.silhouette();
		
	}

	/**
	 * Constructor that sets up the k-means object
	 * @param itrs
	 * @param k
	 * @param pfls
	 * @param otag
	 */
	public ClusterProfiles(int itrs, int k, ArrayList<int[]> pfls, HashMap<Integer,String> pflsIndsMap, int mink, int maxk, HashMap<Integer,Double> pflscores,File odir) {
		numClusItrs = itrs;
		setNumClusters(k);
		setKmerModLenMin(mink);
		setKmerModLenMax(maxk);
		setSparcedProfiles(pfls);
		setProfileInds(pflsIndsMap);
		setProfileScores(pflscores);
		
		this.setOutdir(odir);
		
	}
	

}
