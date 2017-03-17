package org.seqcode.projects.sequnwinder.clusterkmerprofile;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import java.util.Collections;

import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.Vector;


import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.gseutils.Pair;

import org.seqcode.ml.clustering.ClusterRepresentative;

import org.seqcode.ml.clustering.PairwiseElementMetric;
import org.seqcode.ml.clustering.kmeans.KMeansClustering;
import org.seqcode.ml.clustering.vectorcluster.DefaultVectorClusterElement;
import org.seqcode.ml.clustering.vectorcluster.EuclideanDistance;
import org.seqcode.ml.clustering.vectorcluster.Mean;
import org.seqcode.ml.clustering.vectorcluster.VectorClusterElement;



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
	private double bsSizeFraction = 0.1;
	private int minC = 2;
	private int maxC = 10;

	
	// Minimum penetrance of a K-mer in a cluster to be considered
	public final double minKmerProp_global = 0.04;
	

	
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
	
	public int getNumClusters(){
		return K;
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
		StringBuilder sb = new StringBuilder();
		double[][] sh = new double[maxC-minC+1][numBootstraps];
		double[][] ssd = new double[maxC-minC+1][numBootstraps];
		for(int c = minC; c<=maxC; c++){
			for(int b=0; b<numBootstraps; b++){
				ArrayList<VectorClusterElement> currBS = generateBootStrapSample();
				Pair<Pair<Double,Double>,Vector<VectorClusterElement>> currOut = doClustering(c, currBS);
				double currSH = currOut.car().car();
				double currSSD = currOut.car().cdr();
				sh[c-minC][b] = currSH;
				ssd[c-minC][b] = currSSD;
				sb.append(c);sb.append("\t");sb.append(currSH);sb.append("\t");sb.append("Silhouette");sb.append("\n");
				sb.append(c);sb.append("\t");sb.append(currSSD);sb.append("\t");sb.append("SSD");sb.append("\n");
			}
		}
		
		BufferedWriter bwQual = new BufferedWriter(new FileWriter(outdir.getAbsolutePath()+File.separator+"ClusterQuality.tab"));
		bwQual.write(sb.toString());
		bwQual.close();

		double[] mediansSh = new double[maxC-minC+1];
		for(int c=0; c<sh.length; c++){
			Arrays.sort(sh[c]);
			int middle = ((numBootstraps) / 2);
			if(numBootstraps % 2 == 0){
				double medianA = sh[c][middle];
				double medianB = sh[c][middle-1];
				mediansSh[c] = (medianA + medianB) / 2;
			} else{
				mediansSh[c] = sh[c][middle];
			}
		}
		
		double maxMedian = Double.MIN_VALUE;
		int bestC = 0;
		
		for(int c=0; c<mediansSh.length; c++){
			if(mediansSh[c]>maxMedian){
				maxMedian = mediansSh[c];
				bestC = c+minC;
			}
		}
		
		K= bestC;
		
		Vector<VectorClusterElement>  bestCmeans = null;
		double bestSilhouette = Double.MIN_VALUE;
		
		for(int b=0; b<numBootstraps; b++){
			Pair<Pair<Double,Double>,Vector<VectorClusterElement>> currOut = doClustering(bestC, sparse_profiles);
			double currSH = currOut.car().car();
			if(currSH>bestSilhouette){
				bestSilhouette = currSH;
				bestCmeans = currOut.cdr();
			}
			
		}
		
		List<Integer> clusAssignment = writeClusters(bestCmeans);
		
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
	
	/**
	 * 
	 * @param nC
	 * @param data
	 * @return Pair<Pair<silhouette,ssd>,clustermeans>>
	 */
	private Pair<Pair<Double,Double>,Vector<VectorClusterElement>> doClustering(int nC, ArrayList<VectorClusterElement> data){
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
		kmc.setIterations(numClusItrs);
		kmc.clusterElements(data,0.01);
		Pair<Double,Double> clusterquality = new Pair<Double,Double>(kmc.silhouette(),kmc.sumOfSquaredDistance());
		
		return new Pair<Pair<Double,Double>,Vector<VectorClusterElement>>(clusterquality,kmc.getClusterMeans());
		
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
