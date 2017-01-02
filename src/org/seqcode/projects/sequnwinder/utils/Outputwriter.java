package org.seqcode.projects.sequnwinder.utils;

import java.awt.Image;
import java.awt.image.RenderedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.projects.sequnwinder.framework.SeqUnwinderConfig;



public class Outputwriter {
	SeqUnwinderConfig seqConfig;
	
	public Outputwriter(SeqUnwinderConfig seqcon) {
		seqConfig = seqcon;
	}
	
	public void writeDiscrimMotifs() throws IOException{
		//First write a transfac file with all the identified discrim motifs
		File motifsTransfac = new File(seqConfig.getOutDir().getAbsoluteFile()+File.separator+"Discrim_motifs.transfac");
		FileWriter fw = new FileWriter(motifsTransfac);
		BufferedWriter bw = new BufferedWriter(fw);

		for(int m =0; m<seqConfig.getDiscrimMotifs().size(); m++){
			String out = WeightMatrix.printTransfacMatrix(seqConfig.getDiscrimMotifs().get(m),seqConfig.getDiscrimMotifs().get(m).getName().replaceAll("#", ""));
			bw.write(out);
		}
		bw.close();

		// Next, print the ROC values
		File motifsROC = new File(seqConfig.getOutDir().getAbsoluteFile()+File.separator+"Discrim_motifs_roc.tab");
		fw = new FileWriter(motifsROC);
		bw = new BufferedWriter(fw);
		for(String s : seqConfig.getDiscrimMotifRocs().keySet()){
			bw.write(s.replaceAll("#", "")+"\t"+Double.toString(seqConfig.getDiscrimMotifRocs().get(s))+"\n");
		}
		bw.close();

		// Now print the logos
		File motifLogos = new File(seqConfig.getOutDir().getAbsolutePath()+File.separator+"motif_logos");
		motifLogos.mkdirs();
		// Finally, draw the motif logos
		for(WeightMatrix fm : seqConfig.getDiscrimMotifs()){
			File motifFileName = new File(seqConfig.getOutDir().getAbsolutePath()+File.separator+"motif_logos"+File.separator+fm.getName().replaceAll("#", "")+".png");
			org.seqcode.motifs.DrawMotifs.printMotifLogo(fm, motifFileName, 75, fm.getName(), true);
			File motifFileNameRC = new File(seqConfig.getOutDir().getAbsolutePath()+File.separator+"motif_logos"+File.separator+fm.getName().replaceAll("#", "")+"_rc.png");
			org.seqcode.motifs.DrawMotifs.printMotifLogo(WeightMatrix.reverseComplement(fm), motifFileNameRC, 75, fm.getName().replaceAll("#", ""), true);
		}
	}

	public void writeClssifierOutput() throws IOException{
		File classifierOutFile = new File(seqConfig.getOutDir().getAbsolutePath()+File.separator+"classifier.out");
		FileWriter fw = new FileWriter(classifierOutFile);
		BufferedWriter bw= new BufferedWriter(fw);
		bw.write(seqConfig.getClassifierOutput());
		bw.close();
	}
	
	public void writeKmerWeights() throws IOException{
		File weightFile = new File(seqConfig.getOutDir().getAbsolutePath()+File.separator+"kmer_weights.mat");
		FileWriter fw = new FileWriter(weightFile);
		StringBuilder weightSB = new StringBuilder();
		weightSB.append("#Kmer");weightSB.append("\t");
		for(String s : seqConfig.getMNames()){
			weightSB.append(s);weightSB.append("\t");
		}
		weightSB.deleteCharAt(weightSB.length()-1);
		weightSB.append("\n");
		for(int k=seqConfig.getKmin(); k<=seqConfig.getKmax(); k++){
			int N = (int) Math.pow(4, k);
			for (int i = 0; i < N; i++){
				weightSB.append(RegionFileUtilities.int2seq(i, k));weightSB.append("\t");
				int base = seqConfig.getKmerBaseInd(RegionFileUtilities.int2seq(i, k));
				int kind = base + RegionFileUtilities.seq2int(RegionFileUtilities.int2seq(i, k)); 
				for(int l=0; l<seqConfig.getMNames().size(); l++){
					weightSB.append(seqConfig.getKmerWeights().get(seqConfig.getMNames().get(l))[kind]);weightSB.append("\t");
				}
				weightSB.deleteCharAt(weightSB.length()-1);
				weightSB.append("\n");
			}
		}
		BufferedWriter bw= new BufferedWriter(fw);
		bw.write(weightSB.toString());
		bw.close();
	}
	
	public void writeDiscrimScore() throws IOException{
		File discrimScoreFile = new File(seqConfig.getOutDir().getAbsolutePath()+File.separator+"Discrim_motifs.scores");
		FileWriter fw = new FileWriter(discrimScoreFile);
		BufferedWriter bw = new BufferedWriter(fw);
		StringBuilder sb = new StringBuilder();
		sb.append("MotifName");sb.append("\t");
		//First header
		for(String s : seqConfig.getMNames()){
			sb.append(s.replaceAll("#", ""));sb.append("\t");
		}
		sb.deleteCharAt(sb.length()-1);
		sb.append("\n");
		
		List<String> orderdMotifNames = new ArrayList<String>(); 
		for(String modname : seqConfig.getMNames()){
			for(String motname : seqConfig.getDiscrimMotsScore().keySet()){
				if(motname.startsWith(modname+"_c"))
					orderdMotifNames.add(motname);
			}
		}
		
		for(String s : orderdMotifNames){
			sb.append(s.replaceAll("#", ""));sb.append("\t");
			for(double scr : seqConfig.getDiscrimMotsScore().get(s)){
				sb.append(scr);sb.append("\t");
			}
			sb.deleteCharAt(sb.length()-1);
			sb.append("\n");
		}
		bw.write(sb.toString());
		bw.close();
		
		// Now simpler version
		File discrimScoreFileSimple = new File(seqConfig.getOutDir().getAbsolutePath()+File.separator+"Discrim_motifs_simple.scores");
		fw = new FileWriter(discrimScoreFileSimple);
		bw = new BufferedWriter(fw);
		sb = new StringBuilder();
		sb.append("MotifName");sb.append("\t");
		//First header
		for(String s:seqConfig.getMNames()){
			String[] pieces = s.split("#");
			if(pieces.length < 2){
				sb.append(s.replaceAll("#", ""));sb.append("\t");
			}
		}
		sb.deleteCharAt(sb.length()-1);
		sb.append("\n");
		
		orderdMotifNames = new ArrayList<String>(); 
		for(String modname : seqConfig.getMNames()){
			String[] pieces = modname.split("#");
			if(pieces.length > 1)
				continue;
			for(String motname : seqConfig.getDiscrimMotsScore().keySet()){
				String[] piecesMotname = modname.split("#");
				if(piecesMotname.length >1)
					continue;
				if(motname.startsWith(modname+"_c"))
					orderdMotifNames.add(motname);
			}
		}
		
		
		for(String s : orderdMotifNames){
			String[] piecesMotName = s.split("#");
			if(piecesMotName.length>1)
				continue;
			sb.append(s.replaceAll("#", ""));sb.append("\t");
			int i=0;
			for(double scr : seqConfig.getDiscrimMotsScore().get(s)){
				String[] piecesModName = seqConfig.getMNames().get(i).split("#");
				i++;
				if(piecesModName.length>1)
					continue;
				sb.append(scr);sb.append("\t");
			}
			sb.deleteCharAt(sb.length()-1);
			sb.append("\n");
		}
		bw.write(sb.toString());
		bw.close();	
	}
	
	public void makeDiscrimJavaHeatmaps() throws NumberFormatException, IOException{
		// First load the discrim motifs
		String motfilename_full = seqConfig.getOutDir().getAbsolutePath()+"/Discrim_motifs.scores";
		
		Map<String,ArrayList<Double>> vals = new HashMap<String,ArrayList<Double>>();
		ArrayList<String> rnames = new ArrayList<String>();
		ArrayList<String> cnames = new ArrayList<String>();
		
		BufferedReader br = new BufferedReader(new FileReader(motfilename_full));
		String line=null;
		while((line=br.readLine()) != null){
			String[] pieces = line.split("\t");
			if(line.contains("MotifName")){
				for(int p=1; p<pieces.length; p++){
					cnames.add(pieces[p]);
				}
				continue;
			}else{
				rnames.add(pieces[0]);
				vals.put(pieces[0], new ArrayList<Double>());
				for(int p=1; p<pieces.length; p++){
					vals.get(pieces[0]).add(Double.parseDouble(pieces[p]));
				}

			}
		}
		br.close();
		
		HeatMapMaker hmapMaker = new HeatMapMaker(vals,cnames,rnames,false);
		// get motifs
		hmapMaker.loadMotifsFromFile(seqConfig.getOutDir().getAbsolutePath()+File.separator+"Discrim_motifs.transfac");
		Image im = hmapMaker.getImage(hmapMaker.getImageWidth(),hmapMaker.getImageHeight());
		ImageIO.write((RenderedImage) im, "png", new File(seqConfig.getOutDir().getAbsolutePath()+"/Discrim_motifs_heatmap.png"));
		
		// Now the simple version
		String motfilename_simple = seqConfig.getOutDir().getAbsolutePath()+"/Discrim_motifs_simple.scores";
		vals = new HashMap<String,ArrayList<Double>>();
		rnames = new ArrayList<String>();
		cnames = new ArrayList<String>();
		br = new BufferedReader(new FileReader(motfilename_simple));
		line=null;
		while((line=br.readLine()) != null){
			String[] pieces = line.split("\t");
			if(line.contains("MotifName")){
				for(int p=1; p<pieces.length; p++){
					cnames.add(pieces[p]);
				}
				continue;
			}else{
				rnames.add(pieces[0]);
				vals.put(pieces[0], new ArrayList<Double>());
				for(int p=1; p<pieces.length; p++){
					vals.get(pieces[0]).add(Double.parseDouble(pieces[p]));
				}

			}
		}
		br.close();

		hmapMaker = new HeatMapMaker(vals,cnames,rnames,false);
		// get motifs
		hmapMaker.loadMotifsFromFile(seqConfig.getOutDir().getAbsolutePath()+File.separator+"Discrim_motifs.transfac");
		im = hmapMaker.getImage(hmapMaker.getImageWidth(),hmapMaker.getImageHeight());
		ImageIO.write((RenderedImage) im, "png", new File(seqConfig.getOutDir().getAbsolutePath()+"/Discrim_motifs_simple_heatmap.png"));

	}
	
	public void writeHTMLfile() throws IOException{
		// Write the HTML file
		File htmlOut = new File(seqConfig.getOutDir().getAbsoluteFile()+File.separator+"SeqUnwinder_results.html");
		FileWriter fout = new FileWriter(htmlOut);
    	fout.write("<html>\n" +
    			"\t<head><title>SeqUnwinder results ("+seqConfig.getOutbase()+")</title></head>\n" +
    			"\t<style type='text/css'>/* <![CDATA[ */ table, th{border-color: #600;border-style: solid;} td{border-color: #600;border-style: solid;} table{border-width: 0 0 1px 1px; border-spacing: 0;border-collapse: collapse;} th{margin: 0;padding: 4px;border-width: 1px 1px 0 0;} td{margin: 0;padding: 4px;border-width: 1px 1px 0 0;} /* ]]> */</style>\n" +
    			"\t<script language='javascript' type='text/javascript'><!--\nfunction motifpopitup(url) {	newwindow=window.open(url,'name','height=75');	if (window.focus) {newwindow.focus()}	return false;}// --></script>\n" +
    			"\t<script language='javascript' type='text/javascript'><!--\nfunction fullpopitup(url) {	newwindow=window.open(url,'name');	if (window.focus) {newwindow.focus()}	return false;}// --></script>\n" +
    			"\t<body>\n" +
    			"\t<h1>SeqUnwinder results ("+seqConfig.getOutbase()+")</h1>\n" +
    			"");
    	DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
    	Date date = new Date();
    	fout.write("\t<p>SeqUnwinder version "+SeqUnwinderConfig.version+" run completed on: "+dateFormat.format(date));
    	fout.write(" with arguments:\n "+seqConfig.getArgs()+"\n</p>\n");

    	// Write the error on training dataset
    	fout.write("\t<h2>Error on training data</h2>\n");
    	fout.write("\t<table>\n");
    	fout.write("\t\t<tr>" +
    			"\t\t<th>Class</th>\n" +
    			"\t\t<th>TP Rare</th>\n" +
    			"\t\t<th>FP Rare</th>\n" +
    			"\t\t<th>Precision</th>\n" +
    			"\t\t<th>Recall</th>\n" +
    			"\t\t<th>F-Measure</th>\n" +
    			"\t\t<th>MCC</th>\n" +
    			"\t\t<th>ROC Area</th>\n" +
    			"\t\t<th>PRC Area</th>\n");
    	fout.write("\t\t</tr>\n");
    	for(String s : seqConfig.getTrainSetStats().keySet()){
    		fout.write("\t\t<tr>" +
    				"\t\t<td>"+s+"</th>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTrainSetStats().get(s)[0])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTrainSetStats().get(s)[1])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTrainSetStats().get(s)[2])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTrainSetStats().get(s)[3])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTrainSetStats().get(s)[4])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTrainSetStats().get(s)[5])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTrainSetStats().get(s)[6])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTrainSetStats().get(s)[7])+"</td>\n");
    		fout.write("\t\t</tr>\n");
    	}
    	fout.write("\t</table>\n");

    	// Cross-validation error
    	fout.write("\t<h2>Stratified cross-validation</h2>\n");
    	fout.write("\t<table>\n");
    	fout.write("\t\t<tr>" +
    			"\t\t<th>Class</th>\n" +
    			"\t\t<th>TP Rare</th>\n" +
    			"\t\t<th>FP Rare</th>\n" +
    			"\t\t<th>Precision</th>\n" +
    			"\t\t<th>Recall</th>\n" +
    			"\t\t<th>F-Measure</th>\n" +
    			"\t\t<th>MCC</th>\n" +
    			"\t\t<th>ROC Area</th>\n" +
    			"\t\t<th>PRC Area</th>\n");
    	fout.write("\t\t</tr>\n");
    	for(String s : seqConfig.getTestSetStats().keySet()){
    		fout.write("\t\t<tr>" +
    				"\t\t<td>"+s+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTestSetStats().get(s)[0])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTestSetStats().get(s)[1])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTestSetStats().get(s)[2])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTestSetStats().get(s)[3])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTestSetStats().get(s)[4])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTestSetStats().get(s)[5])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTestSetStats().get(s)[6])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTestSetStats().get(s)[7])+"</td>\n");
    		fout.write("\t\t</tr>\n");
    	}
    	fout.write("\t</table>\n");
    	
    	fout.write("\t<p><a href='kmer_weights.mat'> Weighted K-mer model.</a>\n");
    	
    	// Print the de novo identified motifs
    	fout.write("\t<h2>De novo motifs</h2>\n");
    	fout.write("\t<table>\n");
    	fout.write("\t\t<tr>" +
    			"\t\t<th>Class</th>\n" +
    			"\t\t<th>Motif</th>\n");
    	fout.write("\t\t</tr>\n");

    	for(String s: seqConfig.getMNames()){
    		fout.write("\t\t<tr>" +
    				"\t\t<td>"+s+"</td>\n");
    		List<String> selectedMotifNames = new ArrayList<String>();
    		for(String motName : seqConfig.getDiscrimMotifRocs().keySet()){
    			if(motName.startsWith(s+"_c"))
    				selectedMotifNames.add(motName);
    		}
    		if(selectedMotifNames.size()>0){
    			fout.write("\t\t<td>");
    			for(String selecMotName : selectedMotifNames){
    				fout.write("\t\t<img src='motif_logos/"+selecMotName.replaceAll("#", "")+".png'><a href='#' onclick='return motifpopitup(\"motif_logos/"+selecMotName.replaceAll("#", "")+"_rc.png\")'>rc</a>");
    				fout.write("<br>\n");
    			}
    			fout.write("</td>\n");
    		}else{
    			fout.write("\t\t<td>No motif found</td>\n");
    		}
    		fout.write("\t\t</tr>\n");
    	}
    	fout.write("\t</table>\n");
    	
    	// Print motif score heatmaps for both complete and simple versions
    	fout.write("\t<h2>Label specific scores of de novo motifs</h2>\n");
    	fout.write("\t<table>\n");
    	fout.write("\t\t<tr>\n");
    	fout.write("\t\t<th>Extended version</th>\n");
    	fout.write("\t\t<th>Compact version</th>\n");
    	fout.write("\t\t</tr>\n");
    	fout.write("\t\t<tr>\n");
		fout.write("\t\t<td><a href='#' onclick='return fullpopitup(\"Discrim_motifs_heatmap.png\")'><img src='Discrim_motifs_heatmap.png' height='450'></a></td>\n");
		fout.write("\t\t<td><a href='#' onclick='return fullpopitup(\"Discrim_motifs_simple_heatmap.png\")'><img src='Discrim_motifs_simple_heatmap.png' height='450'></a></td>\n");
		fout.write("\t\t</tr>\n");
		fout.write("\t</table>\n");
		
		fout.write("\t<p><a href='Discrim_motifs.transfac'> Discriminative motifs.</a> Try inputting these motifs into <a href='http://www.benoslab.pitt.edu/stamp/'>STAMP</a> for validation.</p>\n");
		
		fout.write("\t</body>\n</html>\n");
    	fout.close();
	}

}
