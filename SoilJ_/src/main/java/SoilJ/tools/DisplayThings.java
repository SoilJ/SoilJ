package SoilJ.tools;

/**
 *SoilJ.tools is a collection of classes for SoilJ, 
 *a collection of ImageJ plugins for the semi-automatized processing of 3-D X-ray images of soil columns
 *Copyright 2014 2015 2016 2017 John Koestel
 *
 *This program is free software: you can redistribute it and/or modify
 *it under the terms of the GNU General Public License as published by
 *the Free Software Foundation, either version 3 of the License, or
 *(at your option) any later version.
 *
 *This program is distributed in the hope that it will be useful,
 *but WITHOUT ANY WARRANTY; without even the implied warranty of
 *MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *GNU General Public License for more details.
 *
 *You should have received a copy of the GNU General Public License
 *along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;

import org.apache.commons.math3.stat.StatUtils;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.plugin.ContrastEnhancer;
import ij.plugin.PlugIn;
import ij.plugin.Slicer;
import ij.process.ByteProcessor;
import ij.process.EllipseFitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import sc.fiji.analyzeSkeleton.Point;

/** 
 * DisplayThings is a SoilJ class to display things (no, really!). It is mostly used for debugging purposes when developing SoilJ.
 * 
 * @author John Koestel
 *
 */

public class DisplayThings implements PlugIn {

	public void run(String arg) {
				//ok, this is not needed..
	}

	public ImagePlus createTiffOfShortestPath(ImagePlus distTiff, MorphologyAnalyzer.FloydWarshallReturn mFWR) {
		
		ImagePlus sTiff = new ImagePlus();
		ImageStack sStack = new ImageStack(distTiff.getWidth(), distTiff.getHeight());
		
		for (int i = 0 ; i < distTiff.getNSlices() ; i++) {
			
			ImageProcessor blankIP = new FloatProcessor(distTiff.getWidth(), distTiff.getHeight());
			blankIP.setColor(0);
			blankIP.fill();
			
			sStack.addSlice(blankIP);
		}
		sTiff.setStack(sStack);
		
		ImageStack drawStack = new ImageStack(distTiff.getWidth(), distTiff.getHeight());
		for (int i = 0 ; i < distTiff.getNSlices() ; i++) {
			
			ArrayList<Point> points2Add = new ArrayList<Point>();
			ArrayList<Float> dist2Add = new ArrayList<Float>();
			
			for (int j = 0 ; j < mFWR.shortest.pathPoints.size() ; j++) {
				int z = mFWR.shortest.pathPoints.get(j).z;
				if (z == i + 1) {
					points2Add.add(mFWR.shortest.pathPoints.get(j));
					double nowDist0 = mFWR.shortest.pathDistances.get(i);
					float nowDist = (float)nowDist0;
					if (nowDist > 1000) nowDist = 0;
					dist2Add.add(nowDist);
				}
			}
			
			sTiff.setPosition(i + 1);
			ImageProcessor sIP = sTiff.getProcessor();
			
			for (int j = 0 ; j < points2Add.size() ; j++) {							
				sIP.putPixelValue(points2Add.get(j).x, points2Add.get(j).y, dist2Add.get(j));
			}
			
			drawStack.addSlice(sIP);
		}
		
		ImagePlus outTiff = new ImagePlus();
		outTiff.setStack(drawStack);
		
		return outTiff;
	}
	
	
	public void showMeMyRoi(String name, ImageProcessor myIP, PolygonRoi pRoi, int colorCode) {
	
		ImagePlus nowTiff = new ImagePlus(name, myIP);
		
		if (pRoi != null) {
			Overlay myO = new Overlay(pRoi);		
			myO.setStrokeColor(Color.YELLOW); // if 1
			if (colorCode == 0) myO.setStrokeColor(Color.WHITE);		
			if (colorCode == 2) myO.setStrokeColor(Color.BLUE);
			if (colorCode == 3) myO.setStrokeColor(Color.RED);		
			nowTiff.setOverlay(myO);
		} else nowTiff.setHideOverlay(true);			
		
		nowTiff.setTitle(name);
		nowTiff.updateAndDraw();
		nowTiff.show();
		
		IJ.wait(2500);
		
		//nowTiff.hide();	
	}
	
	public void showMeMyRoi(String name, ImageProcessor myIP) {
		
		ImagePlus nowTiff = new ImagePlus(name, myIP);
		
		nowTiff.setTitle(name);
		nowTiff.updateAndDraw();
		nowTiff.show();
		
		IJ.wait(2500);
		
		//nowTiff.hide();	
	}
	
	public void displayColumnOutlinesByZ(ImagePlus nowTiff, ObjectDetector.EggShapedColCoords3D colCoords, InputOutput.MyFileCollection mFC, ImagePlus surfaceTiff) {
		
		String pathSep = "/";
				
		InputOutput jIO = new InputOutput();
		
		//standard number of angles is 72	
		int[] iallAngels = {0, 72/8 - 1, 72/4 - 1, 3*72/8 - 1, 72/2 - 1, 5*72/8 - 1, 3*72/4 - 1, 7*72/8 - 1};
		double[] myAngels = {0, Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4};
		double[] allAngels = {0, Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4, Math.PI, 5 * Math.PI / 4, 3 * Math.PI / 2, 7 * Math.PI / 4};
		int radius = (nowTiff.getWidth() + nowTiff.getHeight()) / 4;
		double midX = nowTiff.getWidth() / 2;
		double midY = nowTiff.getHeight() / 2;
		
		float[] x = new float[8];
		float[] y = new float[8];
		float[] xo = new float[8];
		float[] yo = new float[8];
		
		//init output images
		ImageProcessor[] outIP = new ImageProcessor[4];
		ImageProcessor standardIP = new ShortProcessor(2 * radius - 1, nowTiff.getNSlices());
		for (int i = 0 ; i < outIP.length ; i++) outIP[i] = standardIP.duplicate();
		
		float[][][] inner = new float[nowTiff.getNSlices()][4][2];
		float[][][] outer = new float[nowTiff.getNSlices()][4][2];
		
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			
			nowTiff.setPosition(z + 1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			IJ.showStatus("Preparing visualization of column outlines for image-slice " + z + " / " + nowTiff.getNSlices());
					
			//if z within the detected column
			if (z >= colCoords.topOfColumn && z < colCoords.bottomOfColumn) {
				
				//inner edge
				for (int j = 0 ; j < allAngels.length ; j++) {
					x[j] = (float)colCoords.xID[z][iallAngels[j]];			
					y[j] = (float)colCoords.yID[z][iallAngels[j]];
				}
				
				//outer edge
				for (int j = 0 ; j < allAngels.length ; j++) {
					xo[j] = (float)colCoords.xOD[z][iallAngels[j]];
					yo[j] = (float)colCoords.yOD[z][iallAngels[j]];
				}				
				
			}
			else {
				for (int j = 0 ; j < allAngels.length ; j++) {
					x[j]=0;				
					y[j]=0;
					xo[j]=0;					
					yo[j]=0;
				}
			}
			
			for (int angle = 0 ; angle < myAngels.length ; angle++) {
				
				double alpha = myAngels[angle];
				float[][] nowPX = new float[2 * radius - 1][4];
				float[][] nowPY = new float[2 * radius - 1][4];
				
				//get profile				
				int cc = 0;				
				double[] minDist = {radius, radius};
				double[] edgeCoord = new double[2];
				double[] monDist = {radius, radius};
				double[] odgeCoord = new double[2];
				
				for (int r = -radius + 1; r < radius ; r++) {
					
					nowPX[cc][angle] = (float)(midX + r * Math.sin(alpha));
					nowPY[cc][angle] = (float)(midY + r * Math.cos(alpha));					
					
					//if z within the detected column
					if (z >= 0 & z < colCoords.zmid.length) {
					
						//check whether one of the inner edges in crossed by the profile
						double distance1 = Math.sqrt(Math.pow(nowPX[cc][angle] - x[angle],2) + Math.pow(nowPY[cc][angle] - y[angle],2));
						double distance2 = Math.sqrt(Math.pow(nowPX[cc][angle] - x[angle + 4],2) + Math.pow(nowPY[cc][angle] - y[angle + 4],2));
						
						if (distance1 < minDist[0]) {
							minDist[0] = distance1;
							edgeCoord[0] = cc;
						}
						if (distance2 < minDist[1]) {
							minDist[1] = distance2;
							edgeCoord[1] = cc;
						}
						
						//check whether one of the outer edges in crossed by the profile
						double oistance1 = Math.sqrt(Math.pow(nowPX[cc][angle] - xo[angle],2) + Math.pow(nowPY[cc][angle] - yo[angle],2));
						double oistance2 = Math.sqrt(Math.pow(nowPX[cc][angle] - xo[angle + 4],2) + Math.pow(nowPY[cc][angle] - yo[angle + 4],2));
						
						if (oistance1 < monDist[0]) {
							monDist[0] = oistance1;
							odgeCoord[0] = cc;
						}
						if (oistance2 < monDist[1]) {
							monDist[1] = oistance2;
							odgeCoord[1] = cc;
						}
					}
					
					cc++;				
					
				}
				
				//transfer detected edges to inner and outer
				if (z >= 0 & z < colCoords.zmid.length) {
					inner[z][angle][0] = (float)edgeCoord[0];				
					inner[z][angle][1] = (float)edgeCoord[1];
					outer[z][angle][0] = (float)odgeCoord[0];
					outer[z][angle][1] = (float)odgeCoord[1];
				}
				
				//remember profile				
				for (int i = 0 ; i < nowPX.length ; i++) {
					int nowPixelValue = nowIP.getPixel((int)nowPX[i][angle], (int)nowPY[i][angle]);
					outIP[angle].putPixel(i, z, nowPixelValue);
				}				
				
			}
		}
				
		ImageStack outStack = new ImageStack(outIP[0].getWidth(), outIP[0].getHeight());
		ImagePlus outTiff = new ImagePlus();
		
		for (int j = 0 ; j < myAngels.length ; j++) {
			
			//construct inner overlay
			float[] Z = new float[2 * (colCoords.bottomOfColumn - colCoords.topOfColumn) - 1];
			for (int i = 0 ; i < colCoords.bottomOfColumn - colCoords.topOfColumn ; i++) Z[i] = i + colCoords.topOfColumn;
			for (int i = 0 ; i < colCoords.bottomOfColumn - colCoords.topOfColumn ; i++) Z[i + colCoords.bottomOfColumn - colCoords.topOfColumn - 1] = colCoords.bottomOfColumn - i;
			float[] inin = new float[Z.length];
			for (int i = 0 ; i < Z.length / 2 ; i++) inin[i] = inner[colCoords.topOfColumn + i][j][0];
			for (int i = Z.length / 2 ; i < Z.length ; i++) {
				int ii = i - Z.length / 2;
				inin[i] = inner[colCoords.bottomOfColumn - ii - 1][j][1];			
			}
			
			//construct outer Overlays
			float[] ZO = new float[colCoords.bottomOfColumn - colCoords.topOfColumn];
			float[] l1 = new float[ZO.length];
			float[] l2 = new float[ZO.length];
			float[] sU = new float[2 * radius - 1];
			float[] sB = new float[2 * radius - 1];
			float[]	sY = new float[2 * radius - 1];
			for (int i = 0 ; i < ZO.length ; i++) {
				ZO[i] = i + colCoords.topOfColumn;
				l1[i] = outer[i + colCoords.topOfColumn][j][0];
				l2[i] = outer[i + colCoords.topOfColumn][j][1];
			}
			for (int r = 0 ; r < 2*radius - 1 ; r++) sY[r] = r;
			
			//construct soil surfaces in case that it has been specified			
			if (surfaceTiff != null) {
				
				surfaceTiff.setPosition(1);
				ImageProcessor suIP = surfaceTiff.getProcessor();

				surfaceTiff.setPosition(2);
				ImageProcessor sdIP = surfaceTiff.getProcessor();
				
				for (int r = -radius + 1; r < radius ; r++) {
					
					float nowPX = (float)(midX + r * Math.sin(myAngels[j]));
					float nowPY = (float)(midY + r * Math.cos(myAngels[j]));					
					
					int nowPixelValue = suIP.getPixel((int)nowPX, (int)nowPY);
					sU[r + radius] = nowPixelValue;
					
					nowPixelValue = sdIP.getPixel((int)nowPX, (int)nowPY);
					sB[r + radius] = nowPixelValue;
				}
			}			
			
			//enhance contrast and save it
			ContrastEnhancer myCE = new ContrastEnhancer();
			myCE.stretchHistogram(outIP[j], 0.5);
			
			ImageProcessor zIP = outIP[j].duplicate();
			ImageProcessor rgbIP = zIP.convertToRGB();
			rgbIP.setColor(Color.YELLOW);
			
			PolygonRoi pRoi = new PolygonRoi(inin, Z, Roi.POLYGON);
			
			rgbIP.setColor(Color.RED);
			PolygonRoi oRoi1 = new PolygonRoi(l1, ZO, Roi.FREELINE);
			PolygonRoi oRoi2 = new PolygonRoi(l2, ZO, Roi.FREELINE);
						
			rgbIP.setLineWidth(2);
			pRoi.drawPixels(rgbIP);
			oRoi1.drawPixels(rgbIP);
			oRoi2.drawPixels(rgbIP);
			
			//ImagePlus check = new ImagePlus("", rgbIP);
			//check.updateAndDraw();
			//check.show();
			
			//draw also surface in case it was specified
			if (surfaceTiff != null) {
				rgbIP.setColor(Color.YELLOW);
				PolygonRoi uRoi = new PolygonRoi(sY, sU, Roi.FREELINE);
				PolygonRoi dRoi = new PolygonRoi(sY, sB, Roi.FREELINE);
				
				uRoi.drawPixels(rgbIP);
				dRoi.drawPixels(rgbIP);				
			}
			
			outStack.addSlice(rgbIP);
		}
		
		outTiff.setStack(outStack);
				
		ContrastEnhancer myCE = new ContrastEnhancer();
		myCE.stretchHistogram(outTiff, 0.5);
		
		//save the results
		String myOutFolderName = mFC.myOutFolder;
		if (surfaceTiff == null) myOutFolderName += pathSep + "ReviewFoundOutlines";
		else myOutFolderName += pathSep + "ReviewFoundSurfaces";
		
		new File(myOutFolderName).mkdir();
		jIO.tiffSaver(myOutFolderName, mFC.fileName, outTiff);		
	}
	
	public void displayColumnOutlinesByZ(ImagePlus nowTiff, ObjectDetector.ColCoords3D colCoords, InputOutput.MyFileCollection mFC, ImagePlus surfaceTiff) {
		
		// String pathSep = "/";
		String pathSep =  System.getProperty("file.separator");
		
		TailoredMaths math = new TailoredMaths();
		InputOutput jIO = new InputOutput();
		
		double[] myAngels = {0, Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4};
		double[] allAngels = {0, Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4, Math.PI, 5 * Math.PI / 4, 3 * Math.PI / 2, 7 * Math.PI / 4};
		int radius = (nowTiff.getWidth() + nowTiff.getHeight()) / 4;
		double midX = nowTiff.getWidth() / 2;
		double midY = nowTiff.getHeight() / 2;
		
		float[] x = new float[8];
		float[] y = new float[8];
		float[] xo = new float[8];
		float[] yo = new float[8];
		
		// init output images
		ImageProcessor[] outIP = new ImageProcessor[4];												// four output 2d images (only initiate)
		ImageProcessor standardIP = new ShortProcessor(2 * radius - 1, nowTiff.getNSlices());		// average of height & width of nowTiff (why?)
		for (int i = 0 ; i < outIP.length ; i++) outIP[i] = standardIP.duplicate();					// set this size for all four output 2D images
		
		float[][][] inner = new float[colCoords.bottomOfColumn - colCoords.topOfColumn][4][2];		// 3D array: height as column, 4x2 matrix angle + x/y?
		float[][][] outer = new float[colCoords.bottomOfColumn - colCoords.topOfColumn][4][2];
		
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {											// loop through all slices
			
			nowTiff.setPosition(z + 1);																// go to respective slice
			ImageProcessor nowIP = nowTiff.getProcessor();											// and get the 2D image (get the slice)
			
			IJ.showStatus("Preparing visualization of column outlines for image-slice " + (z + 1) + " / " + nowTiff.getNSlices());
					
			//if z within the detected column
			if (z >= 0 & z < colCoords.bottomOfColumn - colCoords.topOfColumn - 1) {				// if z in detected column
			
				double[][] xy = new double[8][8];													// 2D array for storing x and y at the 8 given angles
				
				//inner edge			
				double majRad = colCoords.innerMajorRadius[z];										// get the INNER Major & Minor radius at this slice
				double minRad = colCoords.innerMinorRadius[z];
				
				xy = math.getXYOfEllipseFromAngle(allAngels, colCoords.ixmid[z], colCoords.iymid[z], majRad, minRad, colCoords.itheta[z]);	// fill xy
				
				for (int j = 0 ; j < allAngels.length ; j++) {										// loop through all angles 
					x[j] = (float)xy[j][0];															// set x according to xy
					y[j] = (float)xy[j][1];
				} // loop through angles
				
				//outer edge				
				majRad = colCoords.outerMajorRadius[z];												// get the OUTER Major & Minor radius at this slice			
				minRad = colCoords.outerMinorRadius[z];
				
				xy = math.getXYOfEllipseFromAngle(allAngels, colCoords.xmid[z], colCoords.ymid[z], majRad, minRad, colCoords.theta[z]);	// fill OUTER xy with coordinates 
				
				for (int j = 0 ; j < allAngels.length ; j++) {										// loop throgh all angles
					xo[j] = (float)xy[j][0];														// set OUTER x (xo) accodring to xy
					yo[j] = (float)xy[j][1];
				} // loop through all angles				
			} // if z in detected column
			
			for (int angle = 0 ; angle < myAngels.length ; angle++) {								// loop through all angles again
				
				double alpha = myAngels[angle];														// current angle alpha
				float[][] nowPX = new float[2 * radius - 1][4];										// PX = 2D array with # rows = average width/height of nowTiff, 4 columns
				float[][] nowPY = new float[2 * radius - 1][4];
				
				//get profile				
				int cc = 0;
				double[] minDist = {radius, radius};												// array of 2
				double[] edgeCoord = new double[2];													// array of 2
				double[] monDist = {radius, radius};												// array of 2
				double[] odgeCoord = new double[2];													// array of 2
				
				for (int r = -radius + 1; r < radius ; r++) {										// for image -half image size to +half image size (0 in the middle of the picture!)
					
					nowPX[cc][angle] = (float)(midX + r * Math.sin(alpha));							// nowPX at current angle and index cc[0...image width/height] is in the middle + sin(current angle)
					nowPY[cc][angle] = (float)(midY + r * Math.cos(alpha));					
					
					//if z within the detected column
					if (z >= 0 & z < colCoords.zmid.length) {										// if z in column:
					
						//check whether one of the INNER edges is crossed by the profile
						double distance1 = Math.sqrt(Math.pow(nowPX[cc][angle] - x[angle],2) + Math.pow(nowPY[cc][angle] - y[angle],2));
						double distance2 = Math.sqrt(Math.pow(nowPX[cc][angle] - x[angle + 4],2) + Math.pow(nowPY[cc][angle] - y[angle + 4],2));
						
						if (distance1 < minDist[0]) {
							minDist[0] = distance1;
							edgeCoord[0] = cc;
						}
						if (distance2 < minDist[1]) {
							minDist[1] = distance2;
							edgeCoord[1] = cc;
						}
						
						//check whether one of the OUTER edges in crossed by the profile
						double oistance1 = Math.sqrt(Math.pow(nowPX[cc][angle] - xo[angle],2) + Math.pow(nowPY[cc][angle] - yo[angle],2));
						double oistance2 = Math.sqrt(Math.pow(nowPX[cc][angle] - xo[angle + 4],2) + Math.pow(nowPY[cc][angle] - yo[angle + 4],2));
						
						if (oistance1 < monDist[0]) {
							monDist[0] = oistance1;
							odgeCoord[0] = cc;
						}
						if (oistance2 < monDist[1]) {
							monDist[1] = oistance2;
							odgeCoord[1] = cc;
						}
					} // if z in column
					
					cc++;				
					
				} // for -radius - +radius
				
				//transfer detected edges to inner and outer
				if (z >= 0 & z < colCoords.bottomOfColumn - colCoords.topOfColumn) {
					inner[z][angle][0] = (float)edgeCoord[0];				
					inner[z][angle][1] = (float)edgeCoord[1];
					outer[z][angle][0] = (float)odgeCoord[0];
					outer[z][angle][1] = (float)odgeCoord[1];
				}
				
				//remember profile				
				for (int i = 0 ; i < nowPX.length ; i++) {
					int nowPixelValue = nowIP.getPixel((int)nowPX[i][angle], (int)nowPY[i][angle]);
					outIP[angle].putPixel(i, z, nowPixelValue);
				} //for outline pixels				
				
			} // loop through angles
		} // loop through all slices
				
		//prepare output images
		ImageStack outStack = new ImageStack(outIP[0].getWidth(), outIP[0].getHeight());
		ImagePlus outTiff = new ImagePlus();
		
		for (int j = 0 ; j < myAngels.length ; j++) {												// loop through all angles
			
			//construct inner overlay
			float[] Z = new float[2 * (colCoords.bottomOfColumn - colCoords.topOfColumn)];			// array with 2x column height length
			for (int i = 0 ; i < colCoords.bottomOfColumn - colCoords.topOfColumn ; i++) Z[i] = i + colCoords.topOfColumn;	// fill with z values
			for (int i = 0 ; i < colCoords.bottomOfColumn - colCoords.topOfColumn ; i++) Z[i + colCoords.bottomOfColumn - colCoords.topOfColumn] = colCoords.bottomOfColumn - i;
			float[] inin = new float[Z.length];														// inin with same length as Z 
			for (int i = 0 ; i < Z.length / 2 ; i++) inin[i] = inner[i][j][0];						// first half: fill with inner ellipse
			for (int i = Z.length / 2 ; i < Z.length ; i++) {										// second half: also fill with inner ellipse
				int ii = i - Z.length / 2;
				int endNumber = colCoords.bottomOfColumn - colCoords.topOfColumn - 1;
				inin[i] = inner[endNumber - ii][j][1];			
			}
			
			//construct outer Overlays
			float[] ZO = new float[colCoords.bottomOfColumn - colCoords.topOfColumn];				// array with column height length	
			float[] l1 = new float[ZO.length];
			float[] l2 = new float[ZO.length];
			float[] sU = new float[2 * radius - 1];													// array with output image width/height
			float[] sB = new float[2 * radius - 1];
			float[]	sY = new float[2 * radius - 1];
			for (int i = 0 ; i < ZO.length ; i++) {													// loop through column height
				ZO[i] = i + colCoords.topOfColumn;													// get Z values in "original" tiff (not cut off at column top/Bottom)
				l1[i] = outer[i][j][0];																
				l2[i] = outer[i][j][1];
			} // loop through column height	
			for (int r = 0 ; r < 2*radius - 1 ; r++) sY[r] = r;										// sY is [0...image width/height]
			
			//construct soil surfaces in case that it has been specified
			
			if (surfaceTiff != null) {																	
				
				surfaceTiff.setPosition(1);
				ImageProcessor suIP = surfaceTiff.getProcessor();

				surfaceTiff.setPosition(2);
				ImageProcessor sdIP = surfaceTiff.getProcessor();
				
				for (int r = -radius + 1; r < radius ; r++) {
					
					float nowPX = (float)(midX + r * Math.sin(myAngels[j]));
					float nowPY = (float)(midY + r * Math.cos(myAngels[j]));					
					
					int nowPixelValue = suIP.getPixel((int)nowPX, (int)nowPY);
					sU[r + radius] = nowPixelValue;
					
					nowPixelValue = sdIP.getPixel((int)nowPX, (int)nowPY);
					sB[r + radius] = nowPixelValue;
				} // for  0-centered image width/height
			} // if surfaceTiff null			
			
			//enhance contrast and save it
			ContrastEnhancer myCE = new ContrastEnhancer();
			myCE.stretchHistogram(outIP[j], 0.5);													// set contrast
			
			ImageProcessor zIP = outIP[j].duplicate();												// (j >> angle-loop); duplicate the image at teh current angle
			ImageProcessor rgbIP = zIP.convertToRGB();												// and convert it to RGB (better for export)
																	// set drawing color (outline) to red										
	
			PolygonRoi pRoi = new PolygonRoi(inin, Z, Roi.POLYGON);									// outlines "left & right"	
			PolygonRoi oRoi1 = new PolygonRoi(l1, ZO, Roi.FREELINE);								// outline top and bottom (???)
			PolygonRoi oRoi2 = new PolygonRoi(l2, ZO, Roi.FREELINE);
						
			rgbIP.setLineWidth(2);																	// line width of ROI (colCoord3D)
			rgbIP.setColor(Color.YELLOW);	
			pRoi.drawPixels(rgbIP);			// Polygon around soil													// actually draw the outlines (not only store them >> really manipulate image pRoi)
			rgbIP.setColor(Color.RED);	
			oRoi1.drawPixels(rgbIP);		// outer line right
			rgbIP.setColor(Color.CYAN);
			oRoi2.drawPixels(rgbIP);		// outer line left
			
			//ImagePlus check = new ImagePlus("", rgbIP);
			//check.updateAndDraw();
			//check.show();
			
			//draw also surface in case it was specified
			if (surfaceTiff != null) {
				rgbIP.setColor(Color.YELLOW);
				PolygonRoi uRoi = new PolygonRoi(sY, sU, Roi.FREELINE);
				rgbIP.setColor(Color.RED);
				PolygonRoi dRoi = new PolygonRoi(sY, sB, Roi.FREELINE);
				
				uRoi.drawPixels(rgbIP);
				dRoi.drawPixels(rgbIP);				
			}
			
			outStack.addSlice(rgbIP);																			
		} // loop through all angles
		
		outTiff.setStack(outStack);																	// add to imagePlus
				
		ContrastEnhancer myCE = new ContrastEnhancer();
		myCE.stretchHistogram(outTiff, 0.5);														// adapt contrast again
		
		//save the results
		String myOutFolderName = mFC.myOutFolder;													// get the Folder where to save out of the passed mFC object
		if (surfaceTiff == null) myOutFolderName += pathSep + "ReviewFoundOutlines";				// no surfaceTiff >> add folder for the review of the outlines
		else myOutFolderName += pathSep + "ReviewFoundSurfaces";
		
		new File(myOutFolderName).mkdir();															// actually create the folder
		jIO.tiffSaver(myOutFolderName, mFC.fileName, outTiff);										// save the tiffs
	
		//MV to better understand
		// IJ.saveAsTiff(outTiff,mFC.myOutFolder);
	}
	
	// only for Maria to understand the angle! basically same goal as displayColumnOutlinesByZ
	// Strange in displayColumnOutlinesByZ: the 90° and 270° seem to be mirrored. Small columns in only one quarter of the image are still displayed from all angles
	// for me easier to understand if I play around a bit...
	// I don't understand why I cannot plot the lines in yellow or so
	// This method is of course slower than displayColumnOutlinesByZ, as it draws each ellipse and takes each cross-section. but for me it's more tubelisicher
	public void displayInnerCircleInTiff(ImagePlus nowTiff, ObjectDetector.ColCoords3D colCoords, InputOutput.MyFileCollection mFC, ImagePlus surfaceTiff) {
		
		String pathSep =  System.getProperty("file.separator");
		InputOutput jIO = new InputOutput();
		
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());		// initiate stack for adding all sliced with imprinted colunn
		ImagePlus outTiff = new ImagePlus();												// initiate imagePlus for creating tiff out of stack
	
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {									// loop through all slices of original image
			nowTiff.setPosition(z + 1);														// set position ("height") > which slice
			ImageProcessor zIP = nowTiff.getProcessor().duplicate();						// duplicate the slice
			ImageProcessor rgbIP = zIP.convertToRGB();										// convert this slice to RGB
			rgbIP.setLineWidth(4);															// for plotting the outline
			rgbIP.setColor(Color.YELLOW);													// "

			if (z >= 0 & z >= colCoords.topOfColumn & z <= colCoords.bottomOfColumn) {		// if there is a column recorded in this slice		
				int k = z - colCoords.topOfColumn;												// k: colCoords refer to first recorded slice, not to original slice number
				
				EllipseFitter f = new EllipseFitter();										// for creating an allipse
				f.theta = (colCoords.itheta[k]) + Math.PI / 2.;								// theta seems to refer to deviation from vertical line, not horizontal line... >> add 90°
				f.major = colCoords.innerMajorRadius[k] * 2;								// get major, minor axis, center (x,y) of ellipse)
				f.minor = colCoords.innerMinorRadius[k] * 2;								
				f.xCenter = colCoords.ixmid[k];
				f.yCenter = colCoords.iymid[k];
				f.makeRoi(rgbIP);															// create the ROI in the RGB ImageProcessor	
				PolygonRoi pr = new PolygonRoi(f.xCoordinates, f.yCoordinates, f.nCoordinates, Roi.POLYGON);  // create a polygon Roi
				pr.setStrokeColor(Color.WHITE);												// try to set all colors to most visible (problem: after rotation colors (red/yellow, etc) get lost >> must be grey/white/black)
				// pr.setColor(Color.WHITE);
				rgbIP.setColor(Color.WHITE);	
				rgbIP.setRoi(pr);															// add the ellipse to teh slice
				rgbIP.draw(pr);																// and actually draw it (change pixel values)
				
				if (z == colCoords.bottomOfColumn | z == colCoords.topOfColumn) {rgbIP.fill(pr);} // draw bottom and top
				
				// outer circle
				EllipseFitter fo = new EllipseFitter();										// for creating an allipse
				fo.theta = (colCoords.theta[k]) + Math.PI / 2.;								// theta seems to refer to deviation from vertical line, not horizontal line... >> add 90°
				fo.major = colCoords.outerMajorRadius[k] * 2;								// get major, minor axis, center (x,y) of ellipse)
				fo.minor = colCoords.outerMinorRadius[k] * 2;								
				fo.xCenter = colCoords.xmid[k];
				fo.yCenter = colCoords.ymid[k];
				fo.makeRoi(rgbIP);															// create the ROI in the RGB ImageProcessor	
				PolygonRoi prO = new PolygonRoi(fo.xCoordinates, fo.yCoordinates, fo.nCoordinates, Roi.POLYGON);  // create a polygon Roi
				prO.setStrokeColor(Color.WHITE);												// try to set all colors to most visible (problem: after rotation colors (red/yellow, etc) get lost >> must be grey/white/black)
				// pr.setColor(Color.WHITE);
				rgbIP.setColor(Color.WHITE);	
				rgbIP.setRoi(prO);															// add the ellipse to teh slice
				rgbIP.draw(prO);																// and actually draw it (change pixel values)
				
				if (z == colCoords.bottomOfColumn | z == colCoords.topOfColumn) {rgbIP.fill(prO);} // draw bottom and top
					
			} // if z in detected column
			outStack.addSlice(rgbIP);														// add the slice (potentially with "imprinted" ellipse) to the stack
		}
		outTiff.setStack(outStack);															// make tiff out of stack

		ImageStack someAngles = new ImageStack((nowTiff.getWidth()+nowTiff.getHeight())/2, nowTiff.getNSlices()); // new stack for view from different angles
		ImagePlus anglesOut = new ImagePlus();												// new imagePlus for tiff with view from different angles
		
		Slicer slAngles = new Slicer();														// initiate slicer for reslicing (side-view) 
		int [] alphas = {0, 72, 144, 216, 288};  // uneven number, as 0 and pi are the same (or pi/2 and 3pi/2)... and I want not the same in mirrored! For the run... it must be degrees
		// int [] alphas = {0, 120, 240};
		for (int i = 0; i < alphas.length; i ++) {
			int alpha = alphas[i];
			if (alpha > 0) {IJ.run(outTiff, "Rotate... ", "angle=" + alpha + " grid=1 interpolation=Bilinear stack");}	// if necessary rotate image
			ImagePlus reslicedOut = slAngles.reslice(outTiff);								// reslice the image (view not from top, but from side)
			reslicedOut.setSlice(nowTiff.getHeight()/2);  									// where to take the cross-section: here it would make sense to calculate!! as image is a square the size increases
			ImageProcessor angleIP = reslicedOut.getProcessor();							// take the image with teh cross-section
			ImageProcessor angleIPresized = angleIP.resize((nowTiff.getWidth()+nowTiff.getHeight())/2, nowTiff.getNSlices()); // all images in the stack must have the same size (crop may be the better alternative...)
			someAngles.addSlice(angleIPresized);											// add side-view to stack
		}
		anglesOut.setStack(someAngles);														// image from stack
		
		//save the results
		String myOutFolderName = mFC.myOutFolder;												// get the Folder where to save out of the passed mFC object
		if (surfaceTiff == null) myOutFolderName += pathSep + "ReviewFoundOutlines";			// no surfaceTiff >> add folder for the review of the outlines
		else myOutFolderName += pathSep + "ReviewFoundSurfaces";
		
		new File(myOutFolderName).mkdir();														// actually create the folder
		//jIO.tiffSaver(myOutFolderName, mFC.fileName, outTiff);								// save the full tiff
		jIO.tiffSaver(myOutFolderName, "ell_" + mFC.fileName, anglesOut);						// save the side-views
	}
	
	
	
	
	public void displayColumnSurfacesByZ(ImagePlus nowTiff, InputOutput.MyFileCollection mFC, ImagePlus surfaceTiff) {
		
		String pathSep = "/";	
		
		InputOutput jIO = new InputOutput();
		
		double[] myAngels = {0, Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4};
		int radius = (nowTiff.getWidth() + nowTiff.getHeight()) / 4;
		double midX = nowTiff.getWidth() / 2;
		double midY = nowTiff.getHeight() / 2;
	
		//init output images
		ImageProcessor[] outIP = new ImageProcessor[4];
		ImageProcessor verticalIP = new ByteProcessor(2 * radius - 1, nowTiff.getNSlices());		
		for (int i = 0 ; i < outIP.length ; i++) outIP[i] = verticalIP.duplicate();
		
		for (int z = 1 ; z < nowTiff.getNSlices() ; z++) {
			
			nowTiff.setPosition(z + 1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			for (int angle = 0 ; angle < myAngels.length ; angle++) {
				
				double alpha = myAngels[angle];
				int[][] nowPX = new int[2 * radius - 1][4];
				int[][] nowPY = new int[2 * radius - 1][4];
				
				//get profile				
				int cc = 0;
				
				for (int r = -radius + 1; r < radius ; r++) {					
					nowPX[cc][angle] = (int)Math.round(midX + r * Math.sin(alpha));
					nowPY[cc][angle] = (int)Math.round(midY + r * Math.cos(alpha));
					cc++;
				}
				
				//remember profile				
				for (int i = 0 ; i < 2 * radius - 1 ; i++) {
					int nowPixelValue = nowIP.getPixel(nowPX[i][angle], nowPY[i][angle]);
					outIP[angle].putPixel(i, z, nowPixelValue);					
				}
			}
		}
		
		ImageStack outStack = new ImageStack(outIP[0].getWidth(), outIP[0].getHeight());
		ImagePlus outTiff = new ImagePlus();
		
		for (int j = 0 ; j < myAngels.length ; j++) {
			
			//construct outer Overlays	
			float[] sU = new float[2 * radius - 1];
			float[] sB = new float[2 * radius - 1];
			float[]	sY = new float[2 * radius - 1];			
			for (int r = 0 ; r < 2*radius - 1 ; r++) sY[r] = r;
		
			//construct soil surfaces in case that it has been specified			
			surfaceTiff.setPosition(1);
			ImageProcessor suIP = surfaceTiff.getProcessor();
			
			//ImagePlus sTiff = new ImagePlus(" ", suIP);
			//sTiff.updateAndDraw();
			//sTiff.show();

			//surfaceTiff.setPosition(2);
			//ImageProcessor sdIP = surfaceTiff.getProcessor();
			
			for (int r = -radius + 1; r < radius - 1 ; r++) {
				
				int nPX = (int)Math.round(midX + (r+2) * Math.sin(myAngels[j]));
				int nPY = (int)Math.round(midY + (r+2) * Math.cos(myAngels[j]));					
				
				int nowPixelValue = suIP.getPixel(nPX, nPY);
				sU[r + radius] = nowPixelValue;
				
				//nowPixelValue = sdIP.getPixel(nPX, nPY);
				//sB[r + radius] = nowPixelValue;
			}			
			
			//convert to rgb
			ImageProcessor zIP = outIP[j].duplicate();
			zIP.invert();
			ImageProcessor rgbIP = zIP.convertToRGB();
			rgbIP.setColor(Color.RED);			
			rgbIP.setLineWidth(4);
			
			//draw also surface in case it was specified
			PolygonRoi uRoi = new PolygonRoi(sY, sU, Roi.FREELINE);
			//PolygonRoi dRoi = new PolygonRoi(sY, sB, Roi.FREELINE);
				
			uRoi.drawPixels(rgbIP);
			//dRoi.drawPixels(rgbIP);
						
			//ImagePlus shows = new ImagePlus("sadfa",rgbIP);
			//shows.updateAndDraw();
			//shows.show();
			
			outStack.addSlice(rgbIP);
		}
		
		outTiff.setStack(outStack);
		
		//save the results
		String myOutFolderName = mFC.myOutFolder;
		myOutFolderName += pathSep + "ReviewFoundSurfaces";
		
		new File(myOutFolderName).mkdir();
		jIO.tiffSaver(myOutFolderName, mFC.fileName, outTiff);
		
	}
	
//	/**
//	 * Saves resliced "side-view" of the column with the entered column outlines
//	 * @param particleLabels
//	 * @param nowTiff full extent 3D tiff to which Elliptic Outlines exist [Image Plus]
//	 * @param colCoords elliptic outlines of the column wall (inner + outer) [SoilJ.tools.ObjectDetector.colCoords]
//	 * @param file organizational infos [SoilJ.tools.InputOutput.MyFileCollection]
//	 * @return
//	 */
//	public void displayColumnSurfacesEachSlice(ImagePlus nowTiff, ObjectDetector.ColCoords3D colCoords, 
//			InputOutput.MyFileCollection mFC) {
//		// Operating system dependent string to separate folders etc. 
//		String pathSep = System.getProperty("file.separator");
//		// Initiate a new InputOutput with all its functions (eg for writing file)
//		InputOutput jIO = new InputOutput();
//		TailoredMaths math = new TailoredMaths();
//		// The angles from which do display the output (tiff with outlines)
//		double[] myAngles = {0, Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4};
//		double[] allAngles = TailoredMaths.seqfixedLength(0., Math.PI, 360/2);
//		
//		ImagePlus nowTiff4drawing = nowTiff;
//		
//		ImageProcessor[] outIP = new ImageProcessor[4];						// four imageProcessors to store "results"
//		ImageStack imgStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight() );
//		ImagePlus outTiff = new ImagePlus();
//		
//		for (int z = 1; z < nowTiff4drawing.getNSlices(); z++) {					// loop through the image downwards (z direction)
//			nowTiff4drawing.setPosition(z);											// set nowTiff "height" to z 
//			ImageProcessor currentSlice = nowTiff4drawing.getProcessor();			// get the respective ImageProcessor (2D Image)
//			int start = colCoords.topOfColumn;
//			int end = colCoords.bottomOfColumn;
//			
//			if (z >= start & z <= end) { 
//			
//				double h = currentSlice.getHeight();							// variables for height and width of image
//				double w = currentSlice.getWidth();
//				
//				double x[] = new double[360/2];									// array for storing x values at every 2°
//				double y[] = new double[360/2];
//				double xmidO = colCoords.xmid[z-1];
//				double ymidO = colCoords.ymid[z-1];
//				double xmidI = colCoords.ixmid   [z-1];
//				double ymidI = colCoords.iymid[z-1];
//				double majO = colCoords.outerMajorRadius[z-1];
//				double minO = colCoords.outerMinorRadius[z-1];
//				double majI = colCoords.innerMajorRadius[z-1];
//				double minI = colCoords.innerMinorRadius[z-1];
//				double thetaO = colCoords.theta[z-1];
//				double thetaI = colCoords.itheta[z-1];
//	
//				double [][] outerXY = math.getXYOfEllipseFromAngle(allAngles, xmidO, ymidO, majO, minO, thetaO);
//				double [][] innerXY = math.getXYOfEllipseFromAngle(allAngles, xmidI, ymidI, majI, minI, thetaI);
//				
//				int[] outerX = new int[outerXY.length];
//				for (int i = 0; i < outerXY.length; i++){outerX[i] = (int) outerXY[i][0] + (int) w/2;}
//				int[] outerY = new int[outerXY.length];
//				for (int i = 0; i < outerXY.length; i++){outerY[i] = (int) outerXY[i][1] + (int) h/2;}
//				//PolygonRoi roi = PolygonRoi(outerX, outerY, outerXY.length, Roi.POLYGON);
//				PolygonRoi oRoi2 = new PolygonRoi(outerX, outerY, outerX.length, Roi.FREELINE);
//				ImageProcessor currentSliceDup = currentSlice.duplicate();												// (j >> angle-loop); duplicate the image at teh current angle
//				ImageProcessor currentSliceRGB = currentSliceDup.convertToRGB();
//				currentSliceRGB.setColor(Color.YELLOW);
//				currentSliceRGB.setLineWidth(z);
//				oRoi2.drawPixels(currentSliceRGB);
//				imgStack.addSlice(currentSliceRGB);																// stack again of the individual RGB image processors(?)			
//			} else {									// if z in defined column
//				// otherwise just take the image without drawing an outline, convert to RGB and add to stack
//				ImageProcessor currentSliceDup = currentSlice.duplicate();												// (j >> angle-loop); duplicate the image at teh current angle
//				ImageProcessor currentSliceRGB = currentSliceDup.convertToRGB();
//				imgStack.addSlice(currentSliceRGB);			
//			}
//		} // for each slice (z-direction downwards
//		
//		outTiff.setStack(imgStack);
//		// jIO.tiffSaver2(mFC, outTiff);
//		IJ.saveAsTiff(outTiff,mFC.myOutFolder);
//		for (int i = 0; i < myAngles.length; i++) {
//			
//		} // for "output"-angles
//		
//		System.out.println("done =)");
//	} // displayColumnSurfacesEachSlice
	
	
	
	
	public ImagePlus jDisplayParticleLabels(int[][] particleLabels, ImagePlus imp) {  //this object was written by Doube.. but private in BoneJ.. therefore the copy and past.
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getImageStackSize();
		final int wh = w * h;
		ImageStack stack = new ImageStack(w, h);
		double max = 0;
		for (int z = 0; z < d; z++) {
			float[] slicePixels = new float[wh];
			for (int i = 0; i < wh; i++) {
				slicePixels[i] = (float) particleLabels[z][i];
				max = Math.max(max, slicePixels[i]);
			}
			stack.addSlice(imp.getImageStack().getSliceLabel(z + 1),
					slicePixels);
		}
		ImagePlus impParticles = new ImagePlus(imp.getShortTitle() + "_parts",
				stack);
		impParticles.setCalibration(imp.getCalibration());
		impParticles.getProcessor().setMinAndMax(0, max);
		return impParticles;
	}
	
	public void displayIP(ImageProcessor myIP, String label) {		
		
		ImagePlus test = new ImagePlus(label, myIP);		
		test.draw();test.show();
		
	}

	public void plotDepthProfiles(ObjectDetector.RadialModes myRM, int[] donuts2Plot) {
		
		double[] z = new double[myRM.maskingThreshold.length];
		double[] myData1 = new double[myRM.maskingThreshold.length];
		double[] myData2 = new double[myRM.maskingThreshold.length];
				
		for (int i = 0 ; i < myRM.maskingThreshold.length ; i++) {
					
			z[i] = myRM.maskingThreshold.length - i - 1;	
			myData1[i] = myRM.maskedRadialMinima[i][donuts2Plot[0]];	
			myData2[i] = myRM.maskedRadialMinima[i][donuts2Plot[1]];			
		}
		
		//find maxima and minima
		double[] mini = new double[2];mini[0]=StatUtils.min(myData1);mini[1]=StatUtils.min(myData2);
		double[] maxi = new double[2];maxi[0]=StatUtils.max(myData1);maxi[1]=StatUtils.max(myData2);
		
		Plot rmo = new Plot(donuts2Plot[0] + " and " + donuts2Plot[1],"brightness", "depth");
		
		rmo.setLimits(StatUtils.min(mini)-50, StatUtils.max(maxi)+50, 0, myRM.maskingThreshold.length);
		rmo.setColor(Color.BLUE);
		rmo.addPoints(myData1, z, Plot.LINE);
		rmo.setColor(Color.RED);
		rmo.addPoints(myData2, z, Plot.LINE);		
		rmo.draw();rmo.show();
		
	}
	
	public void plotVerticalProfile(double[] thing2Plot, String title, String xLabel) {
		
		double[] z = new double[thing2Plot.length];
		
		for (int i = 0 ; i < thing2Plot.length ; i++) {
					
			z[i] = thing2Plot.length - i - 1;	
			
		}
		
		//find maxima and minima
		double mini = StatUtils.min(thing2Plot);
		double maxi = StatUtils.max(thing2Plot);
		double std = Math.sqrt(StatUtils.variance(thing2Plot));
		
		Plot rmo = new Plot(title, xLabel, "depth");
		
		rmo.setLimits(mini - 0.1 * std, maxi + 0.1 * std, 0, z.length);
		rmo.setColor(Color.BLUE);
		rmo.addPoints(thing2Plot, z, Plot.LINE);
		rmo.draw();rmo.show();
		
	}

	public void plotRadialProfiles(ObjectDetector.RadialModes myRM, int[] depth2Plot) {
		
		double[] myData1 = new double[myRM.radius.length];
		double[] myData2 = new double[myRM.radius.length];
		
		for (int i = 0 ; i < myRM.radius.length ; i++) {	
			myData1[i] = myRM.maskedRadialMinima[depth2Plot[0]][i];	 
			myData2[i] = myRM.maskedRadialMinima[depth2Plot[1]][i];
		}
		
		//find maxima and minima
		double[] mini = new double[2];mini[0]=StatUtils.min(myData1);mini[1]=StatUtils.min(myData2);
		double[] maxi = new double[2];maxi[0]=StatUtils.max(myData1);maxi[1]=StatUtils.max(myData2);
		
		Plot rmo = new Plot(depth2Plot[0] + " and " + depth2Plot[1], "radial distance", "brightness");	
		
		rmo.setLimits(0, StatUtils.max(myRM.radius), StatUtils.min(mini)-50, StatUtils.max(maxi)+50);
		rmo.setColor(Color.BLUE);
		rmo.addPoints(myRM.radius, myData1, Plot.CIRCLE);
		rmo.setColor(Color.RED);
		rmo.addPoints(myRM.radius, myData2, Plot.CIRCLE);		
		rmo.draw();rmo.show();		
	}
	
	public boolean plotSpecialWallFinderRadii(double[] x, double[] y, double[] ygrad, double[] x2) {
		

		RollerCaster rC = new RollerCaster();
		
		///////////////////////////////////////////////////
		// plot grey values
		////////////////////////////////////////////////////
		
		
		//find maxima and minima
		double mini=StatUtils.min(y);				;
		double maxi=StatUtils.max(y);
		double minmax = maxi-mini;
		
		double[] ylim = {mini-0.05*minmax, maxi+0.05*minmax}; 
		double[] outerWallPerimeter = {x2[0], x2[0]};
		double[] innerWallPerimeter = {x2[1], x2[1]};
		
		
		Plot rgv = new Plot("median filtered grey values", "radial coordinate (vx)", "grey value");	
		
		rgv.setLimits(StatUtils.min(x), StatUtils.max(x), ylim[0], ylim[1]);
		rgv.setColor(Color.BLUE);
		rgv.addPoints(rC.castDouble2Float(x), rC.castDouble2Float(y), Plot.LINE);
		rgv.setColor(Color.RED);
		rgv.addPoints(rC.castDouble2Float(outerWallPerimeter), rC.castDouble2Float(ylim), Plot.LINE);
		rgv.addPoints(rC.castDouble2Float(innerWallPerimeter), rC.castDouble2Float(ylim), Plot.LINE);
		rgv.draw();
		PlotWindow win1 = rgv.show();	

		//find position for second plot
		int newPosX = win1.getX();
		int newPosY = win1.getY() + win1.getHeight();
		
		///////////////////////////////////////////////////
		//plot Gradient
		////////////////////////////////////////////////////
		
		//find maxima and minima
		mini=StatUtils.min(ygrad);				;
		maxi=StatUtils.max(ygrad);
		minmax = maxi-mini;
		
		ylim[0] = mini-0.05*minmax;
		ylim[1] = maxi+0.05*minmax; 
		
		Plot rmo = new Plot("radial gradients of median filtered grey values", "radial coordinate (vx)", "gradient");	
		
		rmo.setLimits(StatUtils.min(x), StatUtils.max(x), mini-3, maxi+3);
		rmo.setColor(Color.GREEN);
		rmo.addPoints(rC.castDouble2Float(x), rC.castDouble2Float(ygrad), Plot.LINE);
		rmo.setColor(Color.RED);
		rmo.addPoints(rC.castDouble2Float(outerWallPerimeter), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.addPoints(rC.castDouble2Float(innerWallPerimeter), rC.castDouble2Float(ylim), Plot.LINE);		
		rmo.draw();
		
		PlotWindow win2 = rmo.show();	
		win2.setLocation(newPosX, newPosY);
				
		GenericDialog gd = new GenericDialog("");
		gd.showDialog();
	    if (gd.wasCanceled()) return false;
	    else new WaitForUserDialog("press any key");
	    gd.removeAll();
	    	    
	    win1.close();
	    win2.close();
	    rgv.dispose();
	    rmo.dispose();
	    
	    IJ.freeMemory();IJ.freeMemory();
	    
	    return true;
		
	}
	
	public boolean plotHistogram(double[] x, double[] y, String title, String xLabel, String yLabel) {
		
		RollerCaster rC = new RollerCaster();		
		
		//find maxima and minima
		double mini=StatUtils.min(y);
		double maxi=StatUtils.max(y);
		
		Plot rmo = new Plot(title, xLabel, yLabel);	
		
		rmo.setLimits(StatUtils.min(x), StatUtils.max(x), mini-3, maxi+3);
		rmo.setColor(Color.BLACK);
		rmo.addPoints(rC.castDouble2Float(x), rC.castDouble2Float(y), Plot.LINE);
				
		//add thresholds
		double[] ylim = {mini, maxi};
	
		//show histogram and thresholds..
		rmo.draw();
		
		PlotWindow win1 = rmo.show();
		GenericDialog gd = new GenericDialog("");
		
/*		gd.showDialog();
	    if (gd.wasCanceled()) return false;
	    else new WaitForUserDialog("press any key");
	    gd.removeAll();
	    	    
	    win1.close();	    
	    rmo.dispose();*/
	    
	    IJ.freeMemory();IJ.freeMemory();
	    
	    IJ.freeMemory();IJ.freeMemory();
	    
	    return true;
		
	}
	
	public boolean plotHistogramWithThresholds(double[] x, double[] y, double[] threshes, String title, String xLabel, String yLabel) {
		
		RollerCaster rC = new RollerCaster();		
		
		//find maxima and minima
		double mini=StatUtils.min(y);
		double maxi=StatUtils.max(y);
		
		Plot rmo = new Plot(title, xLabel, yLabel);	
		
		rmo.setLimits(StatUtils.min(x), StatUtils.max(x), mini-3, maxi+3);
		rmo.setColor(Color.BLACK);
		rmo.addPoints(rC.castDouble2Float(x), rC.castDouble2Float(y), Plot.LINE);
				
		//add thresholds
		double[] ylim = {mini, maxi};
		double[] meth1 = {threshes[0], threshes[0]};
		double[] meth2 = {threshes[1], threshes[1]};
		double[] meth3 = {threshes[2], threshes[2]};
		double[] meth4 = {threshes[3], threshes[3]};
		double[] meth5 = {threshes[4], threshes[4]};
		double[] meth6 = {threshes[5], threshes[5]};
		double[] meth7 = {threshes[6], threshes[6]};
		double[] meth8 = {threshes[7], threshes[7]};
		double[] meth9 = {threshes[8], threshes[8]};
		
		rmo.setColor(Color.RED); //default (also some isodata)
		rmo.addPoints(rC.castDouble2Float(meth1), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.setColor(Color.BLUE); // Otsu			
		rmo.addPoints(rC.castDouble2Float(meth2), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.setColor(Color.GREEN); //Huang
		rmo.addPoints(rC.castDouble2Float(meth3), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.setColor(Color.CYAN);  //max entropy			
		rmo.addPoints(rC.castDouble2Float(meth4), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.setColor(Color.MAGENTA);  //minimum
		rmo.addPoints(rC.castDouble2Float(meth5), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.setColor(Color.ORANGE);   //minError
		rmo.addPoints(rC.castDouble2Float(meth6), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.setColor(Color.YELLOW);   //Renyi Entropy
		rmo.addPoints(rC.castDouble2Float(meth7), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.setColor(Color.PINK);   ///triangle
		rmo.addPoints(rC.castDouble2Float(meth8), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.setColor(Color.GRAY);   //IsoData
		rmo.addPoints(rC.castDouble2Float(meth9), rC.castDouble2Float(ylim), Plot.LINE);		
	
		//show histogram and thresholds..
		rmo.draw();
		
		PlotWindow win1 = rmo.show();
		GenericDialog gd = new GenericDialog("");
		
/*		gd.showDialog();
	    if (gd.wasCanceled()) return false;
	    else new WaitForUserDialog("press any key");
	    gd.removeAll();
	    	    
	    win1.close();	    
	    rmo.dispose();*/
	    
	    IJ.freeMemory();IJ.freeMemory();
	    
	    IJ.freeMemory();IJ.freeMemory();
	    
	    return true;
		
	}
	
	public boolean plotHistogramAdvanced(float[] hist, int myThresh) {
		
		RollerCaster rC = new RollerCaster();		
		HistogramStuff mHS = new HistogramStuff();
		
		//find maxima and minima
		double mini=StatUtils.min(rC.castFloat2Double(hist));
		double maxi=StatUtils.max(rC.castFloat2Double(hist));
		
		//get cumulative histogram
		double xLeft = mHS.findQuantileFromHistogram(rC.castFloat2Int(hist), 0.001);
		double xRight = mHS.findQuantileFromHistogram(rC.castFloat2Int(hist), 0.999);
		
		//create x-axis
		int[] x = new int[hist.length];
		for (int i = 0 ; i < hist.length ; i++) x[i] = i;
		
		Plot rmo = new Plot("Joint Histogram", "Gray Value", "Frequency");	
		
		rmo.setLimits(xLeft, xRight, mini-5, maxi+5);
		rmo.setColor(Color.BLACK);
		rmo.addPoints(rC.castInt2Float(x), hist, Plot.LINE);
				
		//add threshold
		double[] ylim = {mini, maxi};
		float[] thresh = {myThresh, myThresh};
		
		rmo.setColor(Color.RED);
		rmo.addPoints(thresh, rC.castDouble2Float(ylim), Plot.LINE);
	
		//show histogram and thresholds..
		rmo.draw();		
		PlotWindow win1 = rmo.show();
		GenericDialog gd = new GenericDialog("");
		
		//save figure
		
		
	    return true;
		
	}
	
	public boolean plotXYXY(double[] x, double[] y, double[] x2, double[] y2,String title, String xLabel, String yLabel) {
		
		RollerCaster rC = new RollerCaster();	
		
		//find maxima and minima
		//find maxima and minima
		double[] mini = new double[2];mini[0]=StatUtils.min(y);mini[1]=StatUtils.min(y2);
		double[] maxi = new double[2];maxi[0]=StatUtils.max(y);maxi[1]=StatUtils.max(y2);
		
		Plot rmo = new Plot(title, xLabel, yLabel);		
		
		rmo.setLimits(StatUtils.min(x), StatUtils.max(x), 0.99 * StatUtils.min(mini), 1.01 * StatUtils.max(maxi));
		if (StatUtils.min(mini) < 0) rmo.setLimits(StatUtils.min(x), StatUtils.max(x), 1.01 * StatUtils.min(mini), 1.01 * StatUtils.max(maxi));
		if (StatUtils.min(maxi) < 0) rmo.setLimits(StatUtils.min(x), StatUtils.max(x), 1.01 * StatUtils.min(mini), 0.99 * StatUtils.max(maxi));
		//rmo.setLimits(0, 400, 7000, 17000);
		rmo.setColor(Color.BLUE);
		rmo.addPoints(rC.castDouble2Float(x), rC.castDouble2Float(y), Plot.LINE);
		rmo.setColor(Color.RED);
		rmo.addPoints(rC.castDouble2Float(x2), rC.castDouble2Float(y2), Plot.LINE);
		rmo.draw();
		PlotWindow pW = rmo.show();		

		
//		GenericDialog gd = new GenericDialog("");
//		gd.showDialog();
//	    if (gd.wasCanceled()) return false;
//	    else new WaitForUserDialog("press any key");
//	    gd.removeAll();
	    
		IJ.wait(500);
		
	    pW.removeAll();
	    pW.dispose();
	    rmo.dispose();
	    
	    IJ.freeMemory();IJ.freeMemory();
	    
	    return true;
		
	}

	public void plotRadialProfilesWithSpecialFunctionFit(ObjectDetector.RadialModes myRM, FitStuff.FittingResults myResults, int depth2Plot) {
		
		FitStuff fittings = new FitStuff();
		
		FitStuff.HyperbolicFunction2Params mLF = fittings.new HyperbolicFunction2Params(); 
		double[] myData = new double[myRM.radius.length];
		double[] myParams = new double[myResults.numberOfParams];
		double[] realRadius = new double[(int)StatUtils.max(myRM.radius)];
		double[] myFit = new double[realRadius.length];
		
		mLF.setup(StatUtils.max(myRM.radius), myRM.maskedRadialMinima[depth2Plot][myRM.radius.length - 1]);
		
		for (int i = 0 ; i < myRM.radius.length ; i++) myData[i] = myRM.maskedRadialMinima[depth2Plot][i];
		for (int i = 0 ; i < myResults.numberOfParams ; i++) myParams[i] = myResults.params[depth2Plot][i];		
		for (int i = 0 ; i < myFit.length ; i++) {		
			realRadius[i] = i;
			myFit[i] = mLF.value(realRadius[i], myParams);
		}
		
		//find maxima and minima
		double mini = StatUtils.min(myData);
		double maxi = StatUtils.max(myData);
		
		Plot rmo = new Plot("" + depth2Plot, "radial distance", "brightness");	
		
		rmo.setLimits(0, StatUtils.max(myRM.radius), mini - 50, maxi + 50);
		rmo.setColor(Color.BLUE);
		rmo.addPoints(myRM.radius, myData, Plot.CIRCLE);
		rmo.setColor(Color.RED);
		rmo.addPoints(realRadius, myFit, Plot.LINE);		
		rmo.draw();rmo.show();		
	}
}

