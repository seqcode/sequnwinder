package org.seqcode.projects.sequnwinder.utils;


import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.geom.AffineTransform;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.UnaryOperator;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.data.motifdb.WeightMatrixPainter;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.motifs.FreqMatrixImport;
import org.seqcode.viz.paintable.AbstractPaintable;
import org.seqcode.viz.paintable.PaintableFrame;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

public class HeatMapMaker extends AbstractPaintable {

	protected Map<String,ArrayList<Double>> matrix = new HashMap<String,ArrayList<Double>>(); // matrix to draw a heat-map
	protected List<String> colnames = new ArrayList<String>(); // col-names; labels or sub-groups
	protected List<String> rownames = new ArrayList<String>(); // row-names; motif names
	protected List<WeightMatrix> motifs = new ArrayList<WeightMatrix>(); // motifs; should be in the same order as rownames

	private int ScreenSizeX=1000;
	private int ScreenSizeY=1000;
	private final int topBorder=50, bottomBorder=150;
	private final int leftBorder=25, rightBorder=25,colLabSize=25;
	private int topBound, bottomBound, leftBound, rightBound, baseLine, midLine;
	private static int numVals =4;
	private int BoxHeight=40, BoxWidth=40, colorbarWidth = 120, colorbarHeight=15;

	private final double maxExp=0.4, midExp=0, minExp=-0.4;

	private Color expMaxColor = Color.yellow;
	private Color expMidColor = Color.black;
	private Color expMinColor = Color.blue;
	
	private boolean drawMotiflabs = true;
	
	
	public Image getImage(int w, int h){
		Image im = createImage(w,h);
		return im;
	}
	
	public int getImageWidth(){return ScreenSizeX;}
	public int getImageHeight(){return ScreenSizeY;}
	


	public HeatMapMaker(Map<String,ArrayList<Double>> m, List<String> cn, List<String> rn) {
		matrix.putAll(m);
		colnames.addAll(cn);
		rownames.addAll(rn);
		numVals = colnames.size();
		ScreenSizeY = rownames.size()*BoxHeight + 200;
		ScreenSizeX = colnames.size()*BoxWidth + 300;
		
		topBound = topBorder+50;
		bottomBound = ScreenSizeY-bottomBorder;
		leftBound = leftBorder;
		rightBound = ScreenSizeX - rightBorder-200;
		
	}
	
	public void setDrawMotiflabs(boolean b){drawMotiflabs = b;}

	// Load freq matrices
	public void loadMotifsFromFile(String filename) throws NumberFormatException, IOException {
		FreqMatrixImport motifImport = new FreqMatrixImport();
		List<WeightMatrix> tmp = new ArrayList<WeightMatrix>();
		tmp.addAll(motifImport.readTransfacMatricesAsFreqMatrices(filename));
		for(String motname : rownames){
			for(WeightMatrix wm : tmp){
				if(wm.getName().equals(motname)){
					motifs.add(wm);
					break;
				}
			}
		}
	}
	
	private void drawExpColorBar(Graphics2D g2d, int x, int y){
		//Draw colors 
		GradientPaint colorbar = new GradientPaint(x, y, expMinColor, x+colorbarWidth/2, y, expMidColor, false);
		g2d.setPaint(colorbar);
		g2d.fillRect(x, y, colorbarWidth/2, colorbarHeight);
		colorbar = new GradientPaint(x+colorbarWidth/2, y, expMidColor, x+colorbarWidth, y, expMaxColor, false);
		g2d.setPaint(colorbar);
		g2d.fillRect(x+(colorbarWidth/2), y, colorbarWidth/2, colorbarHeight);
		
		//Draw border
		g2d.setPaint(Color.black);
		g2d.setColor(Color.black);
		g2d.setStroke(new BasicStroke(1.0f));
		g2d.drawRect(x, y, colorbarWidth, colorbarHeight);
		
		//Legend
		g2d.setFont(new Font("Ariel", Font.PLAIN, 12));
		FontMetrics metrics = g2d.getFontMetrics();
		int textY = y+colorbarHeight+ (metrics.getHeight());
		g2d.drawString("0", x+(colorbarWidth/2)-(metrics.stringWidth("0")/2), textY);
		g2d.drawString(String.format("%.1f",minExp), x-(metrics.stringWidth(String.format(".1f",minExp))/2), textY);
		g2d.drawString(String.format("%.1f",maxExp), x+colorbarWidth-(metrics.stringWidth(String.format(".1f",maxExp))/2), textY);
		
		//Title
		g2d.setFont(new Font("Ariel", Font.ITALIC, 12));
		metrics = g2d.getFontMetrics();
		g2d.drawString("Model specific discriminative score", x+(colorbarWidth/2)-(metrics.stringWidth("Model specific discriminative score")/2), y- (metrics.getHeight())/2);
	}

	@Override
	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		//motifs.replaceAll(new WMopr());
		Graphics2D g2d = (Graphics2D)g;
		// get the screen sizes
		int screenSizeX = x2-x1;
		int screenSizeY = y2-y1;
		// make the background
		g2d.setColor(Color.white);
		g2d.fillRect(0, 0, screenSizeX, screenSizeY);
		
		baseLine = (topBound+bottomBound)/2;
		midLine = (leftBound+rightBound)/2;
		int xPos = (leftBound+rightBound)/2;
		
		// paint the column names
		//AffineTransform old = g2d.getTransform();
		//g2d.rotate(-Math.toRadians(90));
		//Font collabelFont = new Font("Arial",Font.PLAIN,colLabSize);
		//g2d.setFont(collabelFont);
		//g2d.setColor(Color.BLACK);
		//for(int c=0; c<colnames.size(); c++){
		//	g2d.drawString(colnames.get(c),topBound-5,xPos-(BoxWidth*numVals/2)+(c*BoxWidth)+5);
		//}
		//g2d.setTransform(old);

		

		//Background rounded boxes
		//int j=0;
		for(int r=0; r<rownames.size(); r++){
			String motifName = rownames.get(r);

			boolean found =false;
			int boxX = xPos-BoxWidth*numVals/2;
			int boxY = topBound+BoxHeight*r;
			g2d.setColor(Color.gray);
			g2d.setStroke(new BasicStroke(1.0f));
			g2d.fillRect(boxX, boxY, BoxWidth*numVals,BoxHeight);
			g2d.setColor(Color.darkGray);
			g2d.fillRect(boxX, boxY, BoxWidth*numVals, BoxHeight);

			//Default expression boxes
			int eY = boxY;
			for(int c=0; c<numVals; c++){
				int eX = xPos - BoxWidth*numVals/2 +(c*BoxWidth);
				g2d.setColor(Color.lightGray);
				g2d.fillRect(eX, eY,BoxWidth, BoxHeight);
				g2d.setColor(Color.white);
				g2d.setStroke(new BasicStroke(1.0f));
				g2d.drawRect(eX, eY,BoxWidth, BoxHeight);
			}

			//Find appropriate expression values

			for(int c=0; c<numVals; c++){
				Double v = matrix.get(motifName).get(c);
				Color currCol = expColor(v);
				int eX = xPos - BoxWidth*numVals/2 +(c*BoxWidth);
				g2d.setColor(currCol);
				g2d.fillRect(eX, eY,BoxWidth, BoxHeight);
				g2d.setColor(Color.white);
				g2d.setStroke(new BasicStroke(1.0f));
				g2d.drawRect(eX, eY,BoxWidth, BoxHeight);					
			}

			// now paint the motif
			WeightMatrixPainter wmp = new WeightMatrixPainter();
			if(drawMotiflabs){
				wmp.paint(motifs.get(r),g2d,boxX+BoxWidth*numVals+30,boxY+2,boxX+BoxWidth*numVals+30+120,boxY+2+BoxHeight-4, rownames.get(r), false);
			}else{
				wmp.paint(motifs.get(r),g2d,boxX+BoxWidth*numVals+30,boxY+2,boxX+BoxWidth*numVals+30+120,boxY+2+BoxHeight-4, " ", false);
			}
		}
		// draw the legend
		drawExpColorBar(g2d,xPos-60,topBound+BoxHeight*(1+matrix.keySet().size()));

	}



	public class WMopr implements UnaryOperator<WeightMatrix>{

		@Override
		public WeightMatrix apply(WeightMatrix t) {
			String s1 = WeightMatrix.getMaxLetters(t);
			s1.replace(' ', 'N');
			String s2 = WeightMatrix.getMaxLetters(WeightMatrix.reverseComplement(t));
			s2.replace(' ', 'N');
			int s1Ind = seq2int(s1);
			int s2Ind = seq2int(s2);
			return s1Ind <= s2Ind ? t : WeightMatrix.reverseComplement(t);
		}

	}

	public  int base2int(char base) {
		int intVal = -1;
		switch (base) {
		case 'A':
			intVal = 0;
			break;
		case 'C':
			intVal = 1;
			break;
		case 'G':
			intVal = 2;
			break;
		case 'T':
			intVal = 3;
			break;
		case 'N':
			intVal = 4;
			break;
		default:
			throw new IllegalArgumentException("Invalid character: " + base);
		}
		return intVal;
	}

	public  int seq2int(String seq) {
		int intVal = 0;
		int len = seq.length();

		for (int i = 0; i < len; i++) {
			long currInt = base2int(seq.charAt(i));
			if (currInt == -1) {
				return -1;
			}
			intVal = intVal * 5;
			intVal += currInt;
		}
		return intVal;
	}


	private Color expColor(double v){
		Color c;
		if(v>midExp){
			Color maxColor = expMaxColor;
			Color minColor = expMidColor;
			double sVal = v>maxExp ? 1 : (v-midExp)/(maxExp-midExp);

			int red = (int)(maxColor.getRed() * sVal + minColor.getRed() * (1 - sVal));
			int green = (int)(maxColor.getGreen() * sVal + minColor.getGreen() * (1 - sVal));
			int blue = (int)(maxColor.getBlue() *sVal + minColor.getBlue() * (1 - sVal));
			c = new Color(red, green, blue);
		}else{
			Color maxColor = expMidColor;
			Color minColor = expMinColor;
			double sVal = v<minExp ? 0 : ((midExp-minExp)-(midExp-v))/(midExp-minExp);
			int red = (int)(maxColor.getRed() * sVal + minColor.getRed() * (1 - sVal));
			int green = (int)(maxColor.getGreen() * sVal + minColor.getGreen() * (1 - sVal));
			int blue = (int)(maxColor.getBlue() *sVal + minColor.getBlue() * (1 - sVal));
			c = new Color(red, green, blue);				
		}
		return(c);
	}
	/**
	 * Main method only for testing 
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException{
		ArgParser ap = new ArgParser(args);
		String valsFname = ap.getKeyValue("vals");
		ArrayList<String> rnames = new ArrayList<String>();
		ArrayList<String> cnames = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new FileReader(valsFname));
		Map<String, ArrayList<Double>> vals = new HashMap<String,ArrayList<Double>>();
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
		
		HeatMapMaker runner = new HeatMapMaker(vals,cnames,rnames);
		runner.loadMotifsFromFile(ap.getKeyValue("motifs"));
		runner.setDrawMotiflabs(!ap.hasKey("nolabs"));
		if(ap.hasKey("raster")){
			runner.saveImage(new File("tmpHeatmap.png"), runner.getImageWidth(), runner.getImageHeight(), true);
		}else{
			runner.saveImage(new File("tmpHeatmap.svg"), runner.getImageWidth(), runner.getImageHeight(), false);
		}
	}

}
