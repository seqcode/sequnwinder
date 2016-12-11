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

	
	// Minimum penetrance of a K-mer in a cluster to be considered
	public final double minKmerProp_clus = 0.2;
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
		// Initialize
		ClusterRepresentative<VectorClusterElement> crep = new Mean();
		//Random starts
		Random generator = new Random();
		Vector<VectorClusterElement> starts = new Vector<VectorClusterElement>();
		for(int s=0; s<K; s++){
			boolean found = false;
			while(!found){
				int r = generator.nextInt(sparse_profiles.size());
				if(!added(starts,sparse_profiles.get(r))){
					starts.add(sparse_profiles.get(r));
					found = true;
				}
			}
				
		}
		
		//Initialize clustering
		KMeansClustering<VectorClusterElement> kmc = new KMeansClustering<VectorClusterElement>(metric, crep, starts);
		kmc.setIterations(numClusItrs);
		
		Collection<Cluster<VectorClusterElement>> clusters = kmc.clusterElements(sparse_profiles,0.01);
		Vector<VectorClusterElement> clustermeans = kmc.getClusterMeans();
		
		//Print the clusters
		List<Integer> clusAssignment = writeClusters(clustermeans);
		//Plot the clusters
		Mappable orderedClusters = reorderKmerProfileMaps(clusters);
		//drawClusterHeatmap(orderedClusters);
		printMatrix(orderedClusters);
		
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
	
	
	private Mappable reorderKmerProfileMaps(Collection<Cluster<VectorClusterElement>> clus){
		Mappable ret = null; 
		
		//Mappable features
		double[][] matrix;
		String[] rnames;
		String[] cnames;
		
		// Which colums to retain while drawing the heatmap
		boolean[] keepCol = new boolean[sparse_profiles.get(0).dimension()];
		for(int i=0; i<keepCol.length; i++){
			keepCol[i] = false;
		}
		
		//Which clusters to these columns belong to...
		ArrayList<Integer> colCluster = new ArrayList<Integer>();
		for(int i=0; i<keepCol.length; i++){
			colCluster.add(K+1);
		}
		
		
		int clusID=1;
		for(Cluster<VectorClusterElement> c : clus){
			for(int i=0; i<keepCol.length; i++){
				double colPerc=0;
				for(VectorClusterElement elems : c.getElements()){
					if(elems.getValue(i) > 0){
						colPerc++;
					}
				}
				colPerc = colPerc/c.getElements().size();
				if(colPerc > minKmerProp_clus){
					keepCol[i] = true;
					if(clusID < colCluster.get(i)){
						colCluster.set(i, clusID);
					}
				}
			}
			clusID++;
		}
		
		// Now reorder 
		ArrayIndexComparator comp = new ArrayIndexComparator(colCluster);
		Integer[] indexes = comp.createIndexArray();
		Arrays.sort(indexes, comp);
		
		int sparseLenght = 0;
		for(int i=0;i<indexes.length; i++){
			if(keepCol[indexes[i]])
				sparseLenght++;
			else
				break;
		}
		
		matrix = new double[sparse_profiles.size()][sparseLenght];
		rnames = new String[sparse_profiles.size()];
		cnames = new String[sparseLenght];
		
		//fill the cnames
		for(int j=0; j<sparseLenght; j++){
			cnames[j] = sparse_colnames.get(indexes[j]);
		}
		
		int rowInd = 0;
		for(Cluster<VectorClusterElement> c : clus){
			for(VectorClusterElement elems : c.getElements()){
				for(int j=0; j<sparseLenght; j++){
					matrix[rowInd][j] = elems.getValue(indexes[j]);
				}
				rnames[rowInd] = indToLocation.get(rowInd);
				rowInd++;
			}
		}
		
		ret = new Mappable(matrix, rnames, cnames);
		return ret;
	}
	
	private void printMatrix(Mappable mat) throws IOException{
		StringBuilder sb = new StringBuilder();
		sb.append("Region"+"\t");
		for(int c=0; c<mat.colnames.length; c++){
			sb.append(mat.colnames[c]+"\t");
		}
		sb.deleteCharAt(sb.length()-1);sb.append("\n");
		for(int r=0; r<mat.rownmanes.length; r++){
			sb.append(mat.rownmanes[r]+"\t");
			for(int c=0; c<mat.colnames.length;c++){
				sb.append(mat.matrix[r][c]);sb.append("\t");
			}
			sb.deleteCharAt(sb.length()-1);sb.append("\n");
		}
		
		File matout = new File(outdir.getAbsoluteFile()+File.separator+"K-mer.mat");
		FileWriter ow = new FileWriter(matout);
		BufferedWriter bw = new BufferedWriter(ow);
		bw.write(sb.toString());
		bw.close();
		
		
	}
	
	public void drawClusterHeatmap(Mappable plotMat) throws IOException{
		double[][] matrix = plotMat.matrix;
		HeatChart map = new HeatChart(matrix);
		map.setHighValueColour(new Color(10));
		map.setLowValueColour(new Color(20));
		map.setChartMargin(100);
		map.setAxisLabelsFont(new Font("Ariel",Font.PLAIN,55));
		map.setXValues(plotMat.rownmanes);
		map.setYValues(plotMat.colnames);
		File f = new File(outdir.getAbsoluteFile()+File.separator+"Clusters.png");
		map.saveToFile(f);
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
	
	public class Mappable{
		public double[][] matrix;
		public String[] rownmanes;
		public String[] colnames;
		public Mappable(double[][] m, String[] rnames, String[] cnames) {
			matrix = m;
			rownmanes = rnames;
			colnames = cnames;
		}
	}
	
	public class ArrayIndexComparator implements Comparator<Integer>{
		ArrayList<Integer> list;
		
		public ArrayIndexComparator(ArrayList<Integer> ls) {
			list = ls;
		}
		
		public Integer[] createIndexArray(){
			Integer[] indexes = new Integer[list.size()];
			for(int i=0; i<indexes.length; i++){
				indexes[i] = i;
			}
			return indexes;
		}

		@Override
		public int compare(Integer o1, Integer o2) {
			return list.get(o1).compareTo(list.get(o2));
		}
		
	}

}
