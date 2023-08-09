package SoilJ.tools;

import java.awt.Color;
import java.awt.Polygon;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Vector;

import org.apache.commons.math3.analysis.function.Logistic;
import org.apache.commons.math3.analysis.function.Logistic.Parametric;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import Utilities.Counter3D;
import Utilities.Object3D;
import features.TubenessProcessor;

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

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.OvalRoi;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.Filters3D;
import ij.plugin.GaussianBlur3D;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.plugin.Resizer;
import ij.plugin.Slicer;
import ij.plugin.filter.GaussianBlur;
import ij.process.AutoThresholder;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import inra.ijpb.binary.BinaryImages;
import inra.ijpb.binary.conncomp.FloodFillComponentsLabeling3D;
import inra.ijpb.data.image.Images3D;
import inra.ijpb.label.LabelImages;
import inra.ijpb.plugins.MarkerControlledWatershed3DPlugin;
import inra.ijpb.watershed.ExtendedMinimaWatershed;
import mcib3d.image3d.ImageByte;
import mcib3d.image3d.processing.FillHoles3D;
import process3d.Dilate_;
import process3d.Erode_;

/** 
 * ImageManipulator is one of the two main SoilJ classes (the other being ObjectDetector). 
 * It contains subroutines which changes image contents. 
 * 
 * @author John Koestel
 *
 */

public class ImageManipulator implements PlugIn {
	
	// carries out all sorts of image manipulations (illumination correction, changes tilt of object, applying filters..)
	
	public void run(String arg) {
		//nothing will be run here..		
	}
	
	///////////////////////////////////////////
	// classes
	///////////////////////////////////////////
		
	public class GrayValueCalibrationParameters {
		
		public String[] fileNames;
		public double[] lowerBounds;
		public double[] upperBounds;	
		public int[] index;
		
	}
	
	public class PhaseOfInterestInfo {
		
		public boolean applyThreshold = false;
		public String smallerLarger = null;
		public int threshold = 0;
		public int myPhaseValue = 0;
		public String outName = null;
		
	}
	
	boolean applyThreshold = false;
	String smallerLarger = null;
	int threshold = 0;
	int myPhaseValue = 0;
	String outName = null;
	
	public ImagePlus cutImageInXYPlane(ImagePlus nowTiff, PolygonRoi[] pRoi, boolean cutCanvas) {
		
		RollerCaster rC = new RollerCaster();
		ImagePlus outTiff = new ImagePlus();
		MenuWaiter menu = new MenuWaiter();
		MenuWaiter.PoreSpaceAnalyzerOptions mPSA = menu.new PoreSpaceAnalyzerOptions();
		RoiHandler roi = new RoiHandler();
		
		//cut image in XY-plane
		int minX = nowTiff.getWidth();
		int minY = nowTiff.getHeight();
		int maxX = 0;
		int maxY = 0;
		for (int i = 0 ; i < pRoi.length ; i++) {
			int[] nowX = pRoi[i].getXCoordinates();
			int[] nowY = pRoi[i].getXCoordinates();
			
			double nmiX = StatUtils.min(rC.castInt2Double(nowX)) + pRoi[i].getXBase();
			double nmiY = StatUtils.min(rC.castInt2Double(nowY)) + pRoi[i].getYBase();
			double nmaX = StatUtils.max(rC.castInt2Double(nowX)) + pRoi[i].getXBase();
			double nmaY = StatUtils.max(rC.castInt2Double(nowY)) + pRoi[i].getYBase();
			
			if (nmiX < minX) minX = (int)Math.round(nmiX);
			if (nmiY < minY) minY = (int)Math.round(nmiY);
			if (nmaX > maxX) maxX = (int)Math.round(nmaX);
			if (nmaY > maxY) maxY = (int)Math.round(nmaY);			
		}
	
		//do the cut
		int[] imageDimensions = {maxX-minX, maxY-minY};
		
		mPSA.mRSO.choiceOfRoi = "Cuboid";
		mPSA.mRSO.cubeX1 = minX;
		mPSA.mRSO.cubeY1 = minY;
		mPSA.mRSO.cubeX2 = maxX;
		mPSA.mRSO.cubeY2 = maxY;
		PolygonRoi cutRoi = roi.makeMeAnIndependentRoi(imageDimensions, mPSA.mRSO);			
		ImageStack cutStack = new ImageStack(maxX-minX, maxY-minY);		
	
		for (int i = 0 ; i < nowTiff.getNSlices() ; i++) {
			
			nowTiff.setPosition(i+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			nowIP.setRoi(pRoi[i]);
			nowIP.setColor(255);
			nowIP.fillOutside(pRoi[i]);
			nowIP.resetRoi();
			
			//nowTiff.updateAndDraw();
			//nowTiff.show();
			
			if (cutCanvas) {
				nowIP.setRoi(cutRoi);			
				ImageProcessor cutIP = nowIP.crop();
			
				cutStack.addSlice(cutIP);	
			}
			else cutStack.addSlice(nowIP);
					
		}
		
		outTiff.setStack(cutStack);
		
		return outTiff;		
	}
	
	public ImagePlus calibrateGrayValuesFromList(ImagePlus nowTiff, InputOutput.CaliRoiLocations myCRL) {
		
		HistogramStuff myHS = new HistogramStuff();
		
		ImageStack newStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImagePlus outTiff = new ImagePlus();
		
		double lTarget = 5000;
		double uTarget = 20000;
		
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			
			nowTiff.setSlice(z + 1);
			ImageProcessor myIP = nowTiff.getProcessor();
			
			IJ.showStatus("Correction gray value in slice " + (z + 1) + "/" + nowTiff.getNSlices() + " ...");
			
			//sample rod gray values
			OvalRoi myRodRoi = new OvalRoi(myCRL.rod.x, myCRL.rod.y, myCRL.rod.width, myCRL.rod.height);
			myIP.setRoi(myRodRoi);
			int[] myHist = myIP.getHistogram();
			double upper = myHS.findMedianFromHistogram(myHist);
			
			//sample insu1 gray values
			OvalRoi myInsu1Roi = new OvalRoi(myCRL.insu1.x, myCRL.insu1.y, myCRL.insu1.width, myCRL.insu1.height);
			myIP.setRoi(myInsu1Roi);
			myHist = myIP.getHistogram();
			double l1 = myHS.findMedianFromHistogram(myHist);
			
			//sample insu2 gray values
			OvalRoi myInsu2Roi = new OvalRoi(myCRL.insu2.x, myCRL.insu2.y, myCRL.insu2.width, myCRL.insu2.height);
			myIP.setRoi(myInsu2Roi);
			myHist = myIP.getHistogram();
			double l2 = myHS.findMedianFromHistogram(myHist);			

			double lower = (l1 + l2) / 2;
						
			for (int x = 0 ; x < nowTiff.getWidth() ; x++) {
				for (int y = 0 ; y < nowTiff.getHeight() ; y++) {
					
					double nowPix = myIP.getPixel(x, y);					
					double newPix = (nowPix - lower) / (upper - lower) * (uTarget - lTarget) + lTarget;
					
					myIP.putPixel(x, y, (int)Math.round(newPix));
					
				}
			}
			
			newStack.addSlice(myIP);
		}
		
		outTiff.setStack(newStack);
		return outTiff;
		
	}
	
	public ImagePlus makeSureItIsTheCorrectBinary(ImagePlus nowTiff, int myImagePhase) {
				
		ImageStack nowStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			
			nowTiff.setSlice(z + 1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			ImageProcessor outIP = nowIP.duplicate();
			
			for (int x = 0 ; x < nowStack.getWidth() ; x++) {
				for (int y = 0 ; y < nowStack.getWidth() ; y++) {
					
					int nowPix = nowIP.getPixel(x, y);
					
					if (nowPix == myImagePhase)	outIP.putPixel(x, y, 255);
					else outIP.putPixel(x, y, 0);						
					
				}
			}
			
			nowStack.addSlice(outIP);
			
		}
		
		nowTiff.setStack(nowStack);
		
		return nowTiff;
	}
	
	public void createROIFromInnerCircle(InputOutput.MyFileCollection mFC, PolygonRoi[] pRoi) {
		
		InputOutput jIO = new InputOutput();
		
		ImageStack nowStack = new ImageStack(mFC.nowWidth, mFC.nowHeight);
		
		for (int z = 0 ; z < mFC.nOfSlices ; z++) {
						
			ImageProcessor nowIP = new ByteProcessor(mFC.nowWidth, mFC.nowHeight);
			
			nowIP.setRoi(pRoi[z]);
			nowIP.setBackgroundValue(0);
			
			nowIP.setRoi(pRoi[z]);
			
			nowIP.setColor(255);
			nowIP.fill(pRoi[z]);
			
			nowIP.setColor(0);
			nowIP.fillOutside(pRoi[z]);
			
			nowStack.addSlice(nowIP);
			
		}
			
		ImagePlus outTiff = new ImagePlus("",nowStack);
		
		jIO.tiffSaver(mFC.myOutFolder, mFC.fileName, outTiff);
	}
	
	public ImagePlus cutByMask(InputOutput.MyFileCollection mFC, ImagePlus nowTiff, ImagePlus maskTiff) {
		
		RoiHandler roi = new RoiHandler();
		MenuWaiter menu = new MenuWaiter();
		
		int beforeRightX = 0;
		int beforeBotY = 0;
		
		int topZ = -1;
		int botZ = -1;
		int leftX = -1;
		int rightX = -1;
		int topY = -1;
		int botY = -1;
		
		ImageStack cutStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImagePlus cutTiff = new ImagePlus();
		
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			
			IJ.showStatus("Preparing cutout in layer " + (z + 1) + "/" + nowTiff.getNSlices());
			
			nowTiff.setSlice(z + 1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			ImageProcessor cutIP = nowIP.duplicate();
			
			maskTiff.setSlice(z + 1);
			ImageProcessor maskIP = maskTiff.getProcessor();
							
			for (int x = 0 ; x < nowTiff.getWidth() ; x++) {				
				
				for (int y = 0 ; y < nowTiff.getHeight() ; y++) {
					
					int nowPix = nowIP.getPixel(x, y);				
					int maskPix = maskIP.getPixel(x, y);
					
					if (x > 0) beforeRightX = maskIP.getPixel(x - 1, y);
					else beforeRightX = 0;
					if (y > 0) beforeBotY = maskIP.getPixel(x, y - 1);
					else beforeBotY = 0;
					
					if (maskPix > 0) {						
						cutIP.putPixel(x, y, nowPix);						
						if (topZ < 0) topZ = z + 1;
						if (leftX < 0) leftX = x;
						else if (x < leftX) leftX = x;
						if (topY < 0) topY = y; 
						else if (y < topY) topY = y;  
					}	
					else {
						cutIP.putPixel(x,y,0);
						if (topZ > 0 & botZ < 0) botZ = z + 1;
						else if (z > botZ) botZ = z + 1;
						if (beforeRightX > 0) {
							if(rightX < 0) rightX = x;
							else if (x > rightX) rightX = x;
						}
						if (beforeBotY > 0) {
							if(botY < 0) botY = y;
							else if (y > botY) botY = y;
						}
					}
				}						
			}
			
			cutStack.addSlice(cutIP);
			
		}	
		
		cutTiff.setStack(cutStack);
		
		//cutTiff.updateAndDraw();
		//cutTiff.show();
		
		//clip away unused canvas
		MenuWaiter.ROISelectionOptions mRSO = menu.new ROISelectionOptions();
		
		//catch negative entries
		if (leftX < 0) leftX = 0;
		if (rightX > nowTiff.getWidth() - 1) rightX = nowTiff.getWidth() - 1;
		if (topY < 0) topY = 0;
		if (botY > nowTiff.getHeight() - 1) botY = nowTiff.getHeight() - 1;
		if (topZ < 1) topZ = 1;
		if (botZ > nowTiff.getNSlices()) botZ = nowTiff.getNSlices();
		
		//assign to output
		mRSO.cubeX1 = leftX;
		mRSO.cubeX2 = rightX;
		mRSO.cubeY1 = topY;
		mRSO.cubeY2 = botY;
		mRSO.cubeZ1 = topZ;
		mRSO.cubeZ2 = botZ;
		mRSO.choiceOfRoi = "Cuboid";
		mRSO.cutCanvas = true;
	
		RoiHandler.ColumnRoi colRoi = roi.prepareDesiredRoi(mFC, cutTiff, mRSO, "255", "");
		
		return colRoi.nowTiff;		
		
	}
	
	public ImagePlus cutByMask(ImagePlus nowTiff, PolygonRoi[] myRoi) {
			
		int maxX = (int) myRoi[0].getBounds().getMaxX();
		int minX = (int) myRoi[0].getBounds().getMinX();
		int maxY = (int) myRoi[0].getBounds().getMaxY();
		int minY = (int) myRoi[0].getBounds().getMinY();
		
		ImageStack cutStack = new ImageStack(maxX-minX,maxY-minY);
		ImagePlus cutTiff = new ImagePlus();
		
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			
			IJ.showStatus("Preparing cutout in layer " + (z + 1) + "/" + nowTiff.getNSlices());
			
			nowTiff.setSlice(z + 1);
			ImageProcessor nowIP = nowTiff.getProcessor().duplicate();

			nowIP.setRoi(myRoi[z]);
			ImageProcessor cutIP = nowIP.crop();
			
//			ImagePlus showme = new ImagePlus("test", nowIP);
//			showme.updateAndDraw();showme.show();
			
			cutStack.addSlice(cutIP);
			
		}	
		
		cutTiff.setStack(cutStack);

		return cutTiff;		
		
	}
	
	public ImagePlus removeHoles(ImagePlus nowTiff) {
				
		ImageByte imgByte = new ImageByte(nowTiff);
		FillHoles3D.process(imgByte, 255, Runtime.getRuntime().availableProcessors(), false);
		nowTiff.updateImage();	
		
		return nowTiff;
		
	}
	
	public ImagePlus deleteWatershedsAtEdges(ImagePlus nowTiff) {
			
		ArrayList<Integer> myList = new ArrayList<Integer>();
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			
			nowTiff.setPosition(z + 1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			if (z == 0 | z == nowTiff.getNSlices() - 1) {				
				for (int x = 0 ; x < nowTiff.getWidth() ; x++) {
					for (int y = 0 ; y < nowTiff.getHeight() ; y++) {						
						int nowPixel = (int)Math.round(nowIP.getPixelValue(x, y));
						if (nowPixel > 0 & !myList.contains(nowPixel)) myList.add(nowPixel);
					}
				}				
			}
			else {
				for (int x = 0 ; x < nowTiff.getWidth() ; x++) {					
					//top boundary
					int nowPixel = (int)Math.round(nowIP.getPixelValue(x, 0));
					if (nowPixel > 0 & !myList.contains(nowPixel)) myList.add(nowPixel);
					
					//bottom boundary
					nowPixel = (int)Math.round(nowIP.getPixelValue(x, nowTiff.getHeight() - 1));
					if (nowPixel > 0 & !myList.contains(nowPixel)) myList.add(nowPixel);
				}
					
				for (int y = 0 ; y < nowTiff.getHeight() ; y++) {			
					
					//left boundary
					int nowPixel = (int)Math.round(nowIP.getPixelValue(0, y));
					if (nowPixel > 0 & !myList.contains(nowPixel)) myList.add(nowPixel);
					
					//right boundary
					nowPixel = (int)Math.round(nowIP.getPixelValue(nowTiff.getWidth() - 1, y));
					if (nowPixel > 0 & !myList.contains(nowPixel)) myList.add(nowPixel);		
				}					
			}
		}		
		
		//sweep through the whole image and set all canvas edge touching watersheds to 0 
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			
			nowTiff.setPosition(z + 1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			ByteProcessor outIP = new ByteProcessor(nowIP.getWidth(), nowIP.getHeight());
			
			for (int x = 0 ; x < nowTiff.getWidth() ; x++) {
				for (int y = 0 ; y < nowTiff.getHeight() ; y++) {						
					int nowPixel = (int)Math.round(nowIP.getPixelValue(x, y));
					if (myList.contains(nowPixel)) outIP.putPixel(x, y, 255);
					else outIP.putPixel(x, y, 0);
				}
			}
			
			outStack.addSlice(outIP);
		}
		
		ImagePlus outTiff = new ImagePlus("", outStack);
		
		return outTiff;
	}
	
	public ImagePlus fuseMasks(ImagePlus virginTiff, ImagePlus eroTiff) {
		
		ImageStack fuseStack = new ImageStack(eroTiff.getWidth(), eroTiff.getHeight());
		
		for (int z = 0 ; z < eroTiff.getNSlices() ; z++) {
			
			virginTiff.setPosition(z + 1);			
			eroTiff.setPosition(z + 1);
						
			ImageProcessor vigIP = virginTiff.getProcessor();			
			ImageProcessor eroIP = eroTiff.getProcessor();
			
			ByteProcessor outIP = new ByteProcessor(vigIP.getWidth(), vigIP.getHeight());
			
			for (int x = 0 ; x < outIP.getWidth() ; x++) {
				for (int y = 0 ; y < outIP.getHeight() ; y++) {
		
					int vPix = vigIP.getPixel(x, y);
					int ePix = eroIP.getPixel(x, y);			
					
					int outPix = 0;
					if (vPix > 0 | ePix > 0) outPix = 255;
					
					outIP.putPixel(x, y, outPix);					
				}
			}			
			fuseStack.addSlice(outIP);
		}
		
		ImagePlus outTiff = new ImagePlus("", fuseStack);
		
		return outTiff;
		
	}
	
	public ImagePlus distanceTransformWatershed(ImagePlus nowTiff, MenuWaiter.TobiasWatershedOptionMenu aMO) {
		
		ImagePlus distTiff = new ImagePlus();
		
		int connectivity = 26;
					
		// Compute distance on specified image
		//copied from "DistanceTransformWatershed3D"
		//MorphoLibJ --- Ignacio Arganda-Carreras (ignacio.arganda@ehu.eus)			
		
		ImageStack dist = BinaryImages.distanceMap(nowTiff.getImageStack(), aMO.weights.getShortWeights(), aMO.normalize);
		
		// invert distance map
		Images3D.invert(dist);

		ImageStack result = ExtendedMinimaWatershed.extendedMinimaWatershed(dist, nowTiff.getImageStack(), aMO.dynamic, connectivity, false);
				
		distTiff = new ImagePlus(nowTiff.getShortTitle() + "_dw", result);
		distTiff.setCalibration(nowTiff.getCalibration());
		
		return distTiff;
		
	}
	
	public ImagePlus calc3D(ImagePlus A, String nowGauge, MenuWaiter.Calc3DMenuReturn m3D, ImagePlus B) {
		
		StackCalculator sC= new StackCalculator();
		RoiHandler roi = new RoiHandler();
		InputOutput jIO = new InputOutput();
		
		ImagePlus outTiff0, outTiff = new ImagePlus();
		
		outTiff0 = null;		
		if (m3D.operation.equalsIgnoreCase("+")) outTiff0 = sC.add(A, B);
		if (m3D.operation.equalsIgnoreCase("-")) outTiff0 = sC.subtract(A, B);
		
		//clear outside 
		if (m3D.useInnerCircle) {	
			
			//read InnerCircle file
			ObjectDetector jOD = new ObjectDetector();
			ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();
			int versio = jIO.checkInnerCircleFileVersion(nowGauge);			
			if (versio == 0) jCO = jIO.readInnerCircleVer0(nowGauge);	
			else jCO = jIO.readInnerCircleVer1(nowGauge);	
			
			ImageStack outStack = new ImageStack(outTiff0.getWidth(), outTiff0.getHeight()); 
			PolygonRoi[] nowRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, -1);
			
			for (int i = 0 ; i < outTiff0.getNSlices() ; i++) {
				
				outTiff0.setPosition(i+1);
				ImageProcessor nowIP = outTiff0.getProcessor();
				ImageProcessor outIP = nowIP.duplicate();
				
				outIP.setBackgroundValue(0d);
				outIP.fillOutside(nowRoi[i]);
				
				outStack.addSlice(outIP);
			}
			
			outTiff.setStack(outStack);
			
		}
		else outTiff = outTiff0;
		
		return outTiff;		
	}
	
	public ImagePlus calc3D(ImagePlus A, String operation, ImagePlus B) {
		
		StackCalculator sC= new StackCalculator();
		
		ImagePlus outTiff = new ImagePlus();
				
		if (operation.equalsIgnoreCase("+")) outTiff = sC.add(A, B);
		if (operation.equalsIgnoreCase("-")) outTiff = sC.subtract(A, B);
		
		return outTiff;		
	}
	
	public class FloatIPCalculator {
		
		public ImageProcessor subtract(ImageProcessor aIP, ImageProcessor bIP) {
			
			ImageProcessor outIP = aIP.duplicate();
			
			for (int x = 0 ; x < aIP.getWidth() ; x++) {
				for (int y = 0 ; y < aIP.getHeight() ; y++) {
					float avox = aIP.getPixelValue(x, y);
					float bvox = bIP.getPixelValue(x, y);
					outIP.putPixelValue(x, y, avox - bvox);
				}
			}
			
			return outIP;
			
		}
		
		public ImageProcessor add(ImageProcessor aIP, ImageProcessor bIP) {
			
			ImageProcessor outIP = aIP.duplicate();
			
			for (int x = 0 ; x < aIP.getWidth() ; x++) {
				for (int y = 0 ; y < aIP.getHeight() ; y++) {
					float avox = aIP.getPixelValue(x, y);
					float bvox = bIP.getPixelValue(x, y);
					outIP.putPixelValue(x, y, avox + bvox);
				}
			}
			
			return outIP;
			
		}
		
	}
	
	
	public class StackCalculator {
		
		public ImagePlus add(ImagePlus a, ImagePlus b) {
			
			ImageStack outStack = new ImageStack(a.getWidth(), a.getHeight()); 
			ImagePlus outTiff = new ImagePlus();
			
			for (int z = 0 ; z < a.getNSlices() ; z++) {
				
				a.setPosition(z + 1);
				b.setPosition(z + 1);
				
				ImageProcessor aIP = a.getProcessor();
				ImageProcessor bIP = b.getProcessor();
				
				ImageProcessor outIP = aIP.duplicate();
				
				for (int x = 0 ; x < a.getWidth() ; x++) {
					for (int y = 0 ; y < a.getHeight() ; y++) {
						int avox = aIP.getPixel(x, y);
						int bvox = bIP.getPixel(x, y);
						outIP.putPixel(x, y, avox + bvox);
					}
				}
				
				outStack.addSlice(outIP);
			}
			
			outTiff.setStack(outStack);
			
			return outTiff;
		}
		
		public ImagePlus subtract(ImagePlus a, ImagePlus b) {
			
			ImageStack outStack = new ImageStack(a.getWidth(), a.getHeight()); 
			ImagePlus outTiff = new ImagePlus();
			
			for (int z = 0 ; z < a.getNSlices() ; z++) {
				
				a.setPosition(z + 1);
				b.setPosition(z + 1);
				
				ImageProcessor aIP = a.getProcessor();
				ImageProcessor bIP = b.getProcessor();
				
				ImageProcessor outIP = aIP.duplicate();
				
				for (int x = 0 ; x < a.getWidth() ; x++) {
					for (int y = 0 ; y < a.getHeight() ; y++) {
						int avox = aIP.getPixel(x, y);
						int bvox = bIP.getPixel(x, y);
						outIP.putPixel(x, y, avox - bvox);
					}
				}
				
				outStack.addSlice(outIP);
			}
			
			outTiff.setStack(outStack);
			
			return outTiff;
		}
		
	}
	
	public class SkeletonizerOptions {
		
		public ImagePlus nowTiff;
		
		public int numberOfSlicesAdded;
		
	}
		
	///////////////////////////////////////////////
	// Melitta
	//////////////////////////////////////////////
	
	public ImagePlus applyMedianFilterAndUnsharpMask(ImagePlus nowTiff, MenuWaiter.MedianFilterAndUnsharpMaskReturn mMUS) {				
		
		StackCalculator mSC = new StackCalculator();
		
		//construct some objects
		ImageStack zStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImageStack readyStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());		
		
		//define filter codes			
		int median3DFilter = 11; // 11 is the code for the 3D median filter..	
	
		//apply 3-D median filter		
		readyStack = nowTiff.getStack();
		if (mMUS.medianFilterSizeXDir > 1 & mMUS.medianFilterSizeYDir > 1 & mMUS.medianFilterSizeZDir > 1 ) {
			zStack = Filters3D.filter(readyStack, median3DFilter, mMUS.medianFilterSizeXDir, mMUS.medianFilterSizeYDir, mMUS.medianFilterSizeZDir);
		}
		else {
			zStack = readyStack;
		}
		ImagePlus zTiff = new ImagePlus();
		zTiff.setStack(zStack);		
	
		//apply 3-D unsharp mask
		ImagePlus blurTiff = zTiff.duplicate();
		GaussianBlur3D.blur(blurTiff, mMUS.uMaskStandardDeviationXDir, mMUS.uMaskStandardDeviationYDir, mMUS.uMaskStandardDeviationZDir);
		
		//calculate difference between blurTiff and original		
		ImagePlus diffTiff = mSC.subtract(zTiff,blurTiff);
		
		//and weighted mask
		ImageStack weightedStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		for (int i = 0 ; i < diffTiff.getNSlices() ; i++) {
			diffTiff.setPosition(i+1);
			ImageProcessor nowIP = diffTiff.getProcessor();
			nowIP.multiply(mMUS.uMaskSharpeningWeight);
			weightedStack.addSlice(nowIP);
		}
		ImagePlus unsharpMask = new ImagePlus();
		unsharpMask.setStack(weightedStack);
		
		//sharpened Image
		ImagePlus filtTiff = mSC.add(zTiff, unsharpMask);
		
		//return filtered 3D image
		return filtTiff;
	}
		
	///////////////////////////////////////////////
	// geo and illu correct
	//////////////////////////////////////////////
	
	public ImagePlus cutOffTopAndBottomSteel(ImagePlus nowTiff, int top, int bot) {
	
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImagePlus outTiff = new ImagePlus();
	
		for (int z = top ; z < bot + 1 ; z++) { 
			
			IJ.showStatus("Correcting illumination of slice " + (z - top) + "/" + (bot - top));
			
			//set stack position to the correct depth
			nowTiff.setPosition(z + 1);
			
			//get and correct current slice
			ImageProcessor localIP = nowTiff.getProcessor();			
			outStack.addSlice(localIP);
		}
		
		outTiff.setStack(outStack);		
		
		return outTiff;
	
	}
	
	public ImagePlus beamDeHardenThis(ImagePlus nowTiff, InputOutput.MyFileCollection mFC, MenuWaiter.BeamDeHardeningReturn mBDH) {
			
		InputOutput jIO = new InputOutput();
		ObjectDetector jOD = new ObjectDetector();
		
		int standardRadius = 0;
		ObjectDetector.ColCoordsEssentials3D jEC = jOD.new ColCoordsEssentials3D();
		ImagePlus correctionMaps = new ImagePlus();
		
		//load gauge file	
		float radialMappingFactor = 1.2f;
		if (mBDH.isSteelColumn) {
			
			//read column coordinates
			
			ObjectDetector.EggShapedColCoords3D jCO = jIO.readInnerCircleSteel(mFC);
			
			//sample radial brightness profiles		
			double[][][] radialGrayValues = jOD.getRadialSteelGrayValues(nowTiff, jCO, radialMappingFactor); //layer, angle, radial distance
			standardRadius = (int)Math.round(radialMappingFactor * (StatUtils.percentile(jCO.innerRadius, 50)));
			
			boolean isSteel = true;
			correctionMaps = getBeamHardeningCorrectionMaps(jCO.heightOfColumn, jCO.anglesChecked, standardRadius, radialGrayValues, isSteel);
						
			jEC = jOD.extractColCoordsEssentials3D(jCO);
			
		}
		else {
			
			ObjectDetector.ColCoords3D jCO = jIO.readInnerCircleAlu(mFC);
			
			//sample radial brightness profiles			
			double[][][] radialGrayValues = jOD.getRadialAluGrayValues(nowTiff, jCO, radialMappingFactor); //layer, angle, radial distance
			standardRadius = (int)Math.round(radialMappingFactor * (StatUtils.percentile(jCO.innerMinorRadius, 50)));
				
			boolean isSteel = false;
			correctionMaps = getBeamHardeningCorrectionMaps(jCO.heightOfColumn, jCO.anglesChecked, standardRadius, radialGrayValues,isSteel);
			
			jEC = jOD.extractColCoordsEssentials3D(jCO);
						
		}
		
		ImagePlus outTiff = correct4BeamHardening(nowTiff, jEC, correctionMaps, standardRadius);
		
		return outTiff;
		
	}
	
	public ImagePlus getBeamHardeningCorrectionMaps(int heightOfColumn, int anglesChecked, int standardRadius, double[][][] radialGrayValues, boolean isSteel) {

		ObjectDetector jOD = new ObjectDetector();
		FitStuff jFS = new FitStuff();
		
		//get two different percentile radial illumination profiles..
		double lowerPerc = 60;
		double upperPerc = 80;
		ObjectDetector.RadialFacts lowerRF = jOD.getRadialFacts(heightOfColumn, anglesChecked, standardRadius, radialGrayValues, lowerPerc);			
		ObjectDetector.RadialFacts upperRF = jOD.getRadialFacts(heightOfColumn, anglesChecked, standardRadius, radialGrayValues, upperPerc);				
		
		//smoothen central profile parts
		IJ.showStatus("Smoothing correction functions ...");
		lowerRF = jFS.fitRadialProfile(lowerRF, isSteel);
		upperRF = jFS.fitRadialProfile(upperRF, isSteel);
		
		//create correction maps
		ImageProcessor lowerCorrectionMap = createBeamHardeningCorrectionMap4Steel(lowerRF);
		ImageProcessor upperCorrectionMap = createBeamHardeningCorrectionMap4Steel(upperRF);
				
		ImageStack outStack = new ImageStack(upperCorrectionMap.getWidth(), upperCorrectionMap.getHeight());
		outStack.addSlice(lowerCorrectionMap);
		outStack.addSlice(upperCorrectionMap);
		
		ImagePlus outTiff = new ImagePlus();
		outTiff.setStack(outStack);
		
		//outTiff.updateAndDraw();outTiff.show();
		
		return outTiff;
		
	}
	
	public ImagePlus scaleWithoutMenu(ImagePlus nowTiff, int newWidth, int newHeight, int newDepth, int interpolationMethod) {
		
		ImageProcessor ip = nowTiff.getProcessor();
		boolean averageWhenDownsizing = true;		
		
		int nSlices = nowTiff.getStackSize();
		int w=nowTiff.getWidth(), h=nowTiff.getHeight(), d=nowTiff.getNSlices();
		
		ImagePlus outTiff = nowTiff.createImagePlus();
		
		ImageStack stack1 = nowTiff.getStack();
		ImageStack stack2 = new ImageStack(newWidth, newHeight);
		boolean virtualStack = stack1.isVirtual();
		double min = nowTiff.getDisplayRangeMin();
		double max = nowTiff.getDisplayRangeMax();
		ImageProcessor ip1, ip2;
		int method = interpolationMethod;
		if (w==1 || h==1)
			method = ImageProcessor.NONE;
		for (int i=1; i<=nSlices; i++) {
			IJ.showStatus("Scale: " + i + "/" + nSlices);
			ip1 = stack1.getProcessor(i);
			String label = stack1.getSliceLabel(i);
			ip1.setInterpolationMethod(method);
			ip2 = ip1.resize(newWidth, newHeight, averageWhenDownsizing);
			if (ip2!=null)
				stack2.addSlice(label, ip2);
			IJ.showProgress(i, nSlices);
		}
		outTiff.setStack("", stack2);
		
		IJ.showProgress(1.0);
		int[] dim = nowTiff.getDimensions();
		outTiff.setDimensions(dim[2], dim[3], dim[4]);
		
		if (newDepth != d) { 
			Resizer resizer = new Resizer();
			resizer.setAverageWhenDownsizing(averageWhenDownsizing);
			outTiff = resizer.zScale(outTiff, newDepth, interpolationMethod);
		}
		
		return outTiff;
		
	}
	
	
	/**
	 * 
	 * @param mFC			[InputOutplut.MyFileCollection]
	 * @param nowTiff		[ImagePlus]
	 * @param mBEO			[MenuWaiter.BioPoreExtractionOptions]
	 * @return	iC.run("OR create stack", roundOneTiff, miniTiff0);
	 */
	public ImagePlus extractBioPores(InputOutput.MyFileCollection mFC, ImagePlus nowTiff, MenuWaiter.BioPoreExtractionOptions mBEO) {
		
		ImageCalculator iC = new ImageCalculator();
		Dilate_ dil = new Dilate_();
		MorphologyAnalyzer jMA = new MorphologyAnalyzer();
		
		//init
		int lengthThreshold=2;
		double threshold=60;
		
		//check how large the image is and scale down in case it is too large.. because tubeness 
		// crashes if they are too large (comment out if you have 1000^3 voxels or less) 
		int w = nowTiff.getWidth(), h = nowTiff.getHeight(), d = nowTiff.getNSlices(); 
		double bulkVol = (double) w * (double) h * (double) d;
		int wnew = w;
		int hnew = h;
		int dnew = d;
		if (bulkVol > 1200000000) {
			IJ.showStatus("Reducing image size..");
			nowTiff = binaryScale2HalfSize(nowTiff);
			wnew = nowTiff.getWidth();					// why not overwrite w, h, d?
			hnew = nowTiff.getHeight();
			dnew = nowTiff.getNSlices();
		}
		
		//*************************CalcTubenessOnLargeImage************************************
		double sigma = 1;
		int n_dilations = 1;
		
		TubenessProcessor tp = new TubenessProcessor(sigma, false);
		ImagePlus original = tp.generateImage(nowTiff.duplicate());	// gives scores for how tube-like a pore is >> non-binary!
		
		original = binarize3DFloatImage(original, threshold);		// binarize (0 / 255) the generated image with the above-defined threshold (60)
				
		n_dilations=1;
		for (int i = 1; i < n_dilations + 1 ; i++) original = dil.dilate(original, 255, true);	// add 1 voxel (dilation) to the 255-phase
		
		//double minimum = 3.14 * (sigma + n_dilations) * (sigma + n_dilations) * (sigma + n_dilations) * lengthThreshold;
		//BoneJParticles myParty = jMA.parallelParticleAnalyzer(original, minimum , Double.POSITIVE_INFINITY);
		//original = myParty.myPartyPic;
		
		//original = binarize3DFloatImage(original,1);
				
		//*************************initTubenessOndownscaledImage************************************		
		IJ.showStatus("Reducing image size..");						// >> do the same again, but on the image that has a reduced size!
		nowTiff = binaryScale2HalfSize(nowTiff);					// reduce image size (even if not reduced before)
		
		sigma=1;
		n_dilations = 1;		
		
		tp = new TubenessProcessor(sigma, false);
		ImagePlus miniTiff0 = tp.generateImage(nowTiff);
		
		miniTiff0 = binarize3DFloatImage(miniTiff0, threshold);
				
		n_dilations=1;
		for (int i = 1; i < n_dilations + 1 ; i++) miniTiff0 = dil.dilate(miniTiff0, 255, false); // what does the false refer to??'
				
		//minimum = 3.14 * (sigma + n_dilations) * (sigma + n_dilations) * (sigma + n_dilations) * lengthThreshold;
		//myParty = jMA.parallelParticleAnalyzer(miniTiff0, minimum , Double.POSITIVE_INFINITY);
		
		//miniTiff0 = myParty.myPartyPic;
		//miniTiff0 = binarize3DFloatImage(miniTiff0,1);
		
		//*************************loop1************************************
		double[] sigmas =   {2,4,6,8,
							6,8,11,15};  		// ???
		
		for (int j = 0 ; j < 4 ; j++) {
			
			sigma=sigmas[j];	
			n_dilations = (int)(sigma);
			
			tp = new TubenessProcessor(sigma, false);
			ImagePlus miniTiff = tp.generateImage(nowTiff);
			
			miniTiff = binarize3DFloatImage(miniTiff,threshold);

			for (int i = 1; i < n_dilations + 1 ; i++) miniTiff = dil.dilate(miniTiff, 255, false);
			
			//minimum = Math.pow(n_dilations, 3) * lengthThreshold;
			//myParty = jMA.parallelParticleAnalyzer(miniTiff, minimum , Double.POSITIVE_INFINITY);
			
			//miniTiff = myParty.myPartyPic;
			//miniTiff = binarize3DFloatImage(miniTiff,1);
			
			//miniTiff.updateAndDraw();miniTiff.show();
			
			miniTiff0 = iC.run("OR create stack", miniTiff0, miniTiff);		
			
			//miniTiff0.updateAndDraw();miniTiff0.show();
			
		}
		
		//*************************Rescale, merge with initial round************************************
		
		ImagePlus roundOneTiff = scaleWithoutMenu(miniTiff0, wnew, hnew, dnew, ImageProcessor.BILINEAR);
		
		roundOneTiff = binarize3DFloatImage(roundOneTiff,128);
		
		roundOneTiff = iC.run("OR create stack", original, roundOneTiff);
		
		//*************************TubenessOn2TimesdownscaledImage************************************		
		IJ.showStatus("Reducing image size ince more ..");
		nowTiff = binaryScale2HalfSize(nowTiff);
		
		for (int j = 4 ; j < mBEO.numberOfSigmas - 2 ; j++) {
			
			sigma=sigmas[j];				
			n_dilations = (int)(sigma / 2);
		
			tp = new TubenessProcessor(sigma, false);
			ImagePlus miniTiff = tp.generateImage(nowTiff);
			
			miniTiff = binarize3DFloatImage(miniTiff,threshold);

			for (int i = 1; i < n_dilations + 1 ; i++) miniTiff = dil.dilate(miniTiff, 255, false);
			
			//minimum = Math.pow(n_dilations, 3) * lengthThreshold;
			//myParty = jMA.parallelParticleAnalyzer(miniTiff, minimum , Double.POSITIVE_INFINITY);
			
			//miniTiff = myParty.myPartyPic;
			//miniTiff = binarize3DFloatImage(miniTiff,1);
			
			//miniTiff.updateAndDraw();miniTiff.show();
			
			if (j == 4) miniTiff0 = miniTiff;
			else miniTiff0 = iC.run("OR create stack", miniTiff0, miniTiff);		
			
			//miniTiff0.updateAndDraw();miniTiff0.show();
		}
		
		//************************* Re-scale final round, merge with previous round************************************
		
		miniTiff0 = scaleWithoutMenu(miniTiff0, wnew, hnew, dnew, ImageProcessor.BILINEAR);
		
		miniTiff0 = binarize3DFloatImage(miniTiff0,128);
		
		return iC.run("OR create stack", roundOneTiff, miniTiff0);    // take roundOneTiff and miniTiff0 and combine them with the OR-Operator, create a new image
		
	}
	
	
	
	public ImagePlus extractGravel(InputOutput.MyFileCollection mFC, ImagePlus nowTiff, MenuWaiter.GravelExtractionOptions mGEO) {
		
		Dilate_ dil = new Dilate_();
		Erode_ ero = new Erode_();
		MorphologyAnalyzer jMA = new MorphologyAnalyzer();
		
		//calculate image resolution
		double resolution = 2 * mGEO.voxelsize;
		
		//check how large the image is and scale down in case it is too large.. because tubeness crashes if they are too large (comment out if you have 1000^3 voxels or less) 
		int w = nowTiff.getWidth(), h = nowTiff.getHeight(), d = nowTiff.getNSlices(); 
		double bulkVol = (double)w * (double)h * (double)d;
		int wnew = w;
		int hnew = h;
		int dnew = d;
				
		if (bulkVol > 1200000000) {
			IJ.showStatus("Reducing image size..");
			
			wnew = nowTiff.getWidth() / 2;
			hnew = nowTiff.getHeight() / 2;
			dnew = nowTiff.getNSlices() / 2;
			nowTiff = scaleWithoutMenu(nowTiff, wnew, hnew, dnew, ImageProcessor.BILINEAR);
			
			resolution = 2 * resolution;
		}
						
		//threshold image
		ImagePlus zwischiTiff = binarize3DImage(nowTiff, mGEO.threshold);

		//erode one time..
		zwischiTiff = ero.erode(zwischiTiff, 255, false);
		
		//apply median filter
		if (mGEO.medianFilter > 0) {
			int median3DFilter = 11; // 11 is the code for the 3D median filter..	
			float filterSize = (float)mGEO.medianFilter;
			ImageStack readyStack = zwischiTiff.getStack();
			readyStack = Filters3D.filter(readyStack, median3DFilter, filterSize, filterSize, filterSize);	
			zwischiTiff.setStack(readyStack);					
		}
		
		//dilate one time again..
		zwischiTiff = dil.dilate(zwischiTiff, 255, false);	
		
		//filter out all fractions with equivalent diameters < 2 mm. 
		ImagePlus outTiff = zwischiTiff;
		
		if (mGEO.filterOutSandFraction) {
			double radius = 2d / resolution;
			double minVol = 4 * 3 * Math.PI * Math.pow(radius, 3);
			
					
			FloodFillComponentsLabeling3D myFCCL = new FloodFillComponentsLabeling3D(26, 32);
			ImageStack myLabelStack = myFCCL.computeLabels(outTiff.getStack());					
				
			myLabelStack = LabelImages.volumeOpening(myLabelStack, (int)Math.round(minVol));
			
			outTiff.setStack(myLabelStack);
				
			//rebinarize
			outTiff = binarize3DFloatImage(outTiff,1);			
		}
				
		return outTiff;
		
	}
	
	public ImagePlus extractNonBioPores(InputOutput.MyFileCollection mFC, ImagePlus nowTiff, ImagePlus bioTiff, MenuWaiter.BioPoreExtractionOptions mBEO) {
		
		//check if image sizes of nowTiff and outTiff are identical and if not scale and resegment
		if (nowTiff.getWidth() > bioTiff.getWidth()) {
			nowTiff = scaleWithoutMenu(nowTiff, bioTiff.getWidth(), bioTiff.getHeight(), bioTiff.getNSlices(), ImageProcessor.BILINEAR);
			nowTiff = binarize3DImage(nowTiff, 128);
		}
		
		return calc3D(nowTiff, "-", bioTiff);
		
	}
		
//	public ImagePlus beamDeHardenThisPVCDEPRECATED(ImagePlus nowTiff, ObjectDetector.ColCoords3D jPCO, MenuWaiter.BeamDeHardeningReturn mBDH) {
//	
//		//init objects	
//		ObjectDetector jOD = new ObjectDetector(); 
//		FitStuff fittings = new FitStuff();
//		ImageManipulator jIM = new ImageManipulator();
//		
//		//init variables
//		ImagePlus outTiff = new ImagePlus();			
//		
//		//get modes of radial histograms
//		ObjectDetector.RadialModes myModes = jOD.getRadialPVCIllumination(nowTiff, jPCO, mBDH);
//		
//		//fit a curve to the matrix brightness data..
//		//FitStuff.FittingResults myLogisticFit = fittings.fitGeneralizedLogisticFunction(myModes, mBDH.maxEval);
//					
//		//fill gaps and bad fits
//		//double[][] filteredGLFParams = jIM.fillGapsInFittedBeamHardeningCorrectionMap(myLogisticFit, myModes, mBDH);
//		
//		//assemble beam hardening correction map
//		//ImageProcessor blurIP = assembleCorrectionMapForBeamHardening(myModes, filteredGLFParams);
//		
////		DisplayThings dT  = new DisplayThings();
////		dT.displayIP(blurIP, "blurIP");			
//		
//		//apply brightness correction to collected minima..
//		//myModes = applyBeamHardeningCorrection2AirPhaseGammaData(myModes, blurIP);		
//		
//		//fit a curve to the air-phase brightness data..
//		//FitStuff.FittingResults myHyperbolicFit = fittings.fitHyperbolicFunction2Params(myModes, mBDH.maxEval);
//		
//		//fill gaps and bad fits
//		//ImageProcessor gammaIP = assembleCorrectionMapForAirPhaseGamma(myLogisticFit, myHyperbolicFit, myModes, mBDH);
//		
////		dT.displayIP(gammaIP, "gammaIP");	
//
//		//correct image for beam hardening and air-phase gamma
//		//outTiff = correct4PVCBeamHardening(nowTiff, blurIP, gammaIP, jPCO, mBDH);
//		
//		//return corrected tiffs..
//		return outTiff;
//		
//	}	
	
	public ObjectDetector.RadialModes applyBeamHardeningCorrection2AirPhaseGammaData(ObjectDetector.RadialModes myModes, ImageProcessor blurIP) { 
		
		for (int i = 0 ; i < myModes.maskingThreshold.length ; i++) {
			for (int j = 0 ; j < myModes.radius.length - 1 ; j++) {
				double nowCorrFac = blurIP.getPixelValue((int)Math.floor(myModes.radius[j]) + 1, i);
				double nowRef = blurIP.getPixelValue(blurIP.getWidth() - 1, i);
				myModes.maskedRadialMinima[i][j] = myModes.maskedRadialMinima[i][j] * nowRef / nowCorrFac;
			}
		}
		
		return myModes;
	}	
	
	public ImageProcessor assembleCorrectionMapForBeamHardening(ObjectDetector.RadialModes myModes, double[][] filteredGLFParams) {
		
		TailoredMaths maths = new TailoredMaths();
		
		//create a correction function map
		double[] radius = new double[(int)maths.max(myModes.radius)];
		for (int i = 0 ; i < radius.length ; i++) radius[i] = i;
		Parametric myLogL = new Logistic.Parametric();
		FloatProcessor corrIP = new FloatProcessor(radius.length, myModes.maskingThreshold.length);
		for (int y = 0 ; y < myModes.maskingThreshold.length ; y++) {
			double[] nowGLFP = new double[6];
			for (int i = 0 ; i < 6 ; i++) nowGLFP[i] = filteredGLFParams[y][i];
			for (int x = 0 ; x < radius.length ; x++) {	
				double myValue2Put = myModes.maskedRadialModes[y][myModes.radius.length - 1];
				try {myValue2Put = myLogL.value(x, nowGLFP);} 
				catch (NotStrictlyPositiveException e) {} 
				corrIP.putPixelValue(x, y, myValue2Put);
			}
		}
		
		//ImagePlus letsSee = new ImagePlus("nonFiltered", corrIP);
		//letsSee.updateAndDraw();letsSee.show();
		
		GaussianBlur myGB = new GaussianBlur();
		double accuracy = 0.01d;
		FloatProcessor blurIP = (FloatProcessor)corrIP.duplicate();
		myGB.blurFloat(blurIP, 20, 20, accuracy);			
		//ImagePlus letsSeeM = new ImagePlus("Filtered",blurIP);
		//letsSeeM.updateAndDraw();letsSeeM.show();
		
		return blurIP;
	}
	
	public ImageProcessor assembleCorrectionMapForAirPhaseGamma(ObjectDetector.RadialModes myModes, double[][] myGapFilledHFParameters) {
		
		FitStuff fittings = new FitStuff();
		TailoredMaths maths = new TailoredMaths();
		
		//create a correction function map
		double[] radius = new double[(int)maths.max(myModes.radius)];
		for (int i = 0 ; i < radius.length ; i++) radius[i] = i;
		FitStuff.HyperbolicFunction2Params myHF = fittings.new HyperbolicFunction2Params();
		FloatProcessor corrIP = new FloatProcessor(radius.length, myModes.maskingThreshold.length);
		for (int y = 0 ; y < myModes.maskingThreshold.length ; y++) {
			double[] nowHFparams = new double[2];
			myHF.setup(StatUtils.max(myModes.radius), myModes.maskedRadialMinima[y][myModes.radius.length - 1]);
			for (int i = 0 ; i < 2 ; i++) nowHFparams[i] = myGapFilledHFParameters[y][i];
			for (int x = 0 ; x < radius.length ; x++) {	
				double myValue2Put = myModes.maskedRadialModes[y][myModes.radius.length - 1];
				try {myValue2Put = myHF.value(x, nowHFparams);} 
				catch (NotStrictlyPositiveException e) {} 
				corrIP.putPixelValue(x, y, myValue2Put);
			}
		}
		
		//ImagePlus letsSee = new ImagePlus("nonFiltered", corrIP);
		//letsSee.updateAndDraw();letsSee.show();
		
		GaussianBlur myGB = new GaussianBlur();
		double accuracy = 0.01d;
		FloatProcessor blurIP = (FloatProcessor)corrIP.duplicate();
		myGB.blurFloat(blurIP, 20, 20, accuracy);			
		//ImagePlus letsSeeM = new ImagePlus("Filtered",blurIP);
		//letsSeeM.updateAndDraw();letsSeeM.show();
		
		return blurIP;
	}
	
	public SkeletonizerOptions addSlices4Skeletonization(ImagePlus nowTiff, MenuWaiter.PoreSpaceAnalyzerOptions mPSA, RoiHandler.ColumnRoi colRoi) {
		
		SkeletonizerOptions mSO = new SkeletonizerOptions();
		
		ImagePlus outTiff = new ImagePlus();
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());

		mSO.numberOfSlicesAdded = mPSA.numberOfCopiedSlices2Add + 1; 
		
		ImageProcessor newIP = new ByteProcessor(nowTiff.getWidth(), nowTiff.getHeight());
		
		//see that the phase to be analyzed is 255
		nowTiff = extractPhaseOfInterest(nowTiff, mPSA.imagePhase2BeAnalyzed, mPSA.nameOfAnalyzedPhase);
		
		newIP.setColor(0);
		newIP.setRoi(0, 0, nowTiff.getWidth(), nowTiff.getHeight());
		newIP.fill();
		ImageProcessor blackIP = newIP.duplicate(); 
				
		newIP.setColor(255);
		newIP.setRoi(1, 1, nowTiff.getWidth() - 2, nowTiff.getHeight() - 2);
		newIP.fill();
		
		nowTiff.setPosition(1);
		ImageProcessor firstIP = nowTiff.getProcessor().duplicate();
		
		nowTiff.setPosition(nowTiff.getNSlices());
		ImageProcessor lastIP = nowTiff.getProcessor().duplicate();
				
		//add a black slice on top
		outStack.addSlice(blackIP);
	
		//add first slices
		for (int i = 0 ; i < mPSA.numberOfCopiedSlices2Add ; i++) {
			outStack.addSlice(firstIP);
		}
		
		//add image
		for (int i = 0 ; i < nowTiff.getNSlices() ; i++) {
			nowTiff.setPosition(i+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			outStack.addSlice(nowIP);
		}
			
		//add last slices
		for (int i = 0 ; i < mPSA.numberOfCopiedSlices2Add ; i++) {
			outStack.addSlice(lastIP);
		}
		
		//add a black slices on bottom	
		outStack.addSlice(blackIP);
		
		outTiff.setStack(outStack);
	
		mSO.nowTiff = outTiff;
		
		return mSO;		
		
	}
	
	public double[][] fillGapsInFittedBeamHardeningCorrectionMap(FitStuff.FittingResults myFit, ObjectDetector.RadialModes myModes, MenuWaiter.BeamDeHardeningReturn mBDH) {
		
		RollerCaster cast = new RollerCaster();
		
		double[] R2 = myFit.R2;
		int standardRadius = (int)Math.round(StatUtils.max(myModes.radius));
		double[] realRadius = new double[standardRadius];
		ArrayList<Integer> testOnes = new ArrayList<Integer>();
		ArrayList<Integer> goodOnes = new ArrayList<Integer>();		
		
		double[][] gapFilledParameters = new double[R2.length][myFit.numberOfParams];
		
		//sample the maximal beamhardening artifact in the soil
		for (int i = R2.length/3 ; i < 2 * R2.length / 3 ; i++) if (R2[i] > mBDH.severeGoodnessCriterion) testOnes.add(i);		
		double[] A = new double[testOnes.size()];			//minimum of radial brightness profile
		double[] K = new double[testOnes.size()];			//maximum of radial brightness profile
		for (int i = 0 ; i < testOnes.size() ; i++) {
			A[i] = myFit.params[testOnes.get(i)][4];
			K[i] = myFit.params[testOnes.get(i)][0];
		}
		double expectedBrightnessDifference = StatUtils.percentile(K, 50) - StatUtils.percentile(A, 50);
		
		//set all R2 of fits with more than the 2 expected brightness difference to 0;
		for (int i = 0 ; i < R2.length ; i++) if (myFit.params[i][0] - myFit.params[i][4] > expectedBrightnessDifference | myFit.params[i][0] - myFit.params[i][4] < 0) {
			R2[i] = 0;
		}
		
		//find first and last good fit
		int firstGood = 0;
		int lastGood = 0;
		for (int i = 0 ; i < R2.length ; i++) {
			if (firstGood == 0 & R2[i] > mBDH.severeGoodnessCriterion) firstGood = i;
			if (R2[i] > mBDH.severeGoodnessCriterion) lastGood = i;
		}
		
		//plug in close-to-wall values before firstGood and after lastGood
		for (int i = 0 ; i < firstGood ; i++) {
			myFit.params[i][0] = myModes.maskedRadialModes[i][myModes.radius.length - 1];
			myFit.params[i][1] = myFit.params[firstGood][1];
			myFit.params[i][2] = myFit.params[firstGood][2];
			myFit.params[i][3] = myFit.params[firstGood][3];
			myFit.params[i][4] = myModes.maskedRadialModes[i][myModes.radius.length - 1] - 1;
			myFit.params[i][5] = myFit.params[firstGood][5];
			R2[i] = mBDH.goodnessCriterion + 0.01;
		}
		for (int i = lastGood + 1 ; i < R2.length ; i++) {
			myFit.params[i][0] = myModes.maskedRadialModes[i][myModes.radius.length - 1];
			myFit.params[i][1] = myFit.params[lastGood][1];
			myFit.params[i][2] = myFit.params[lastGood][2];
			myFit.params[i][3] = myFit.params[lastGood][3];
			myFit.params[i][4] = myModes.maskedRadialModes[i][myModes.radius.length - 1] - 1;
			myFit.params[i][5] = myFit.params[lastGood][5];
			R2[i] = mBDH.goodnessCriterion + 0.01;
		}
		
		//pick out the good ones and imputed ones and make sure that a value for the first and last cross-section is set
		boolean imputedFirst = false;
		boolean imputedLast = false;
		for (int i = 0 ; i < R2.length ; i++) if (R2[i] > mBDH.goodnessCriterion) {			
			goodOnes.add(i);			
		}
		else if (i == 0 | i == R2.length - 1) {
			goodOnes.add(i);
			if (i == 0) imputedFirst = true;
			if (i == R2.length - 1) imputedLast = true;
		}
		
		//re-mould them in vectors and matrices
		int[] depth = new int[goodOnes.size()];
		if (imputedFirst) {
			for (int j = 0 ; j < myFit.numberOfParams ; j++) myFit.params[0][j] = myFit.params[goodOnes.get(1)][j];
		}
		if (imputedLast) {
			for (int j = 0 ; j < myFit.numberOfParams ; j++) myFit.params[goodOnes.size() - 1][j] = myFit.params[goodOnes.get(goodOnes.size() - 2)][j];
		}		
		
		double[][] data = new double[goodOnes.size()][myFit.numberOfParams];		
		for (int i = 0 ; i < depth.length ; i++) {
			depth[i] = goodOnes.get(i);
			for (int j = 0 ; j < myFit.numberOfParams ; j++) data[i][j] = myFit.params[depth[i]][j];
		}
		
		//init realRadius
		for (int i = 0 ; i < realRadius.length ; i++) realRadius[i] = i;
		
		//interpolate the gaps for all parameters
		for (int i = 0 ; i < myFit.numberOfParams ; i++) {			
			double[] nowData = new double[depth.length];
			for (int j = 0 ; j < nowData.length ; j++) nowData[j] = data[j][i];
			LinearInterpolator myLI = new LinearInterpolator();		
			PolynomialSplineFunction mySF = myLI.interpolate(cast.castInt2Double(depth), nowData);			 
			for (int j = 0 ; j < R2.length ; j++) gapFilledParameters[j][i] = mySF.value(j);						
		}
		
		return gapFilledParameters;
		
	}
	
	public ImageProcessor assembleCorrectionMapForAirPhaseGamma(FitStuff.FittingResults myMatrixBrightnessFit, FitStuff.FittingResults myHFFit, ObjectDetector.RadialModes myModes, MenuWaiter.BeamDeHardeningReturn mBDH) {
		
		FitStuff fittings = new FitStuff();
		RollerCaster cast = new RollerCaster();
		TailoredMaths maths = new TailoredMaths();
				
		double[] R2 = myHFFit.R2;
		double[] matrixR2 = myMatrixBrightnessFit.R2;
		int standardRadius = (int)Math.round(StatUtils.max(myModes.radius));
		double[] realRadius = new double[standardRadius];
		ArrayList<Integer> testOnes = new ArrayList<Integer>();
		ArrayList<Integer> goodOnes = new ArrayList<Integer>();		
		
		double[][] gapFilledParameters = new double[R2.length][myHFFit.numberOfParams];
		
		//sample the maximal beam-hardening artifact in the soil ... using the soil matrix fits here..
		for (int i = matrixR2.length/3 ; i < 2 * matrixR2.length / 3 ; i++) if (matrixR2[i] > mBDH.severeGoodnessCriterion) testOnes.add(i);		
		double[] A = new double[testOnes.size()];			//minimum of radial brightness profile
		double[] K = new double[testOnes.size()];			//maximum of radial brightness profile
		for (int i = 0 ; i < testOnes.size() ; i++) {
			A[i] = myMatrixBrightnessFit.params[testOnes.get(i)][4];
			K[i] = myMatrixBrightnessFit.params[testOnes.get(i)][0];
		}
		double expectedBrightnessDifference = StatUtils.percentile(K, 50) - StatUtils.percentile(A, 50);
		
		//set all R2 of fits with more than the 2 expected brightness difference to 0;
		for (int i = 0 ; i < R2.length ; i++) {
			if (myMatrixBrightnessFit.params[i][0] - myMatrixBrightnessFit.params[i][4] > expectedBrightnessDifference 
					| myMatrixBrightnessFit.params[i][0] - myMatrixBrightnessFit.params[i][4] < 0) {		
				R2[i] = 0;
			}
		}
				
		//perform a similar check for the air-phase gamma correction function
		FitStuff.HyperbolicFunction2Params myHF = fittings.new HyperbolicFunction2Params();		
		for (int i = 0 ; i < R2.length ; i++) {
			myHF.setup(StatUtils.max(myModes.radius), myModes.maskedRadialMinima[i][myModes.radius.length - 1]);
			double refMatrixBrightness = myModes.maskedRadialModes[i][myModes.radius.length - 1]; 
			double refAirPhaseBrightness = myModes.maskedRadialMinima[i][myModes.radius.length - 1];
			double[] nowParams = new double[2];			
			for (int j = 0 ; j < 2 ; j++) nowParams[j] = myHFFit.params[i][j];
			double fittedAirPhaseBrightness = myHF.value(0, nowParams);
			if ( (refMatrixBrightness - refAirPhaseBrightness) / (refMatrixBrightness - fittedAirPhaseBrightness) > 3) R2[i] = 0;			
		}		
		
		//find first and last good fit .. using the matrix fits here..
		int firstGood = 0;
		int lastGood = 0;
		for (int i = 0 ; i < matrixR2.length ; i++) {
			if (firstGood == 0 & matrixR2[i] > mBDH.severeGoodnessCriterion) firstGood = i;
			if (matrixR2[i] > mBDH.severeGoodnessCriterion) lastGood = i;
		}
		
		//plug in close-to-wall values before firstGood and after lastGood
		for (int i = 0 ; i < firstGood ; i++) {			
			myHFFit.params[i][1] = 0;
			R2[i] = mBDH.gammaGoodnessCriterion + 0.01;
		}
		for (int i = lastGood + 1 ; i < matrixR2.length ; i++) {			
			myHFFit.params[i][1] = 0;
			R2[i] = mBDH.gammaGoodnessCriterion + 0.01;
		}
	
		//pick out the good ones and imputed ones and make sure that a value for the first and last cross-section is set
		boolean imputedFirst = false;
		boolean imputedLast = false;
		for (int i = 0 ; i < R2.length ; i++) if (R2[i] > mBDH.gammaGoodnessCriterion) {			
			goodOnes.add(i);			
		}
		else if (i == 0 | i == R2.length - 1) {
			goodOnes.add(i);
			if (i == 0) imputedFirst = true;
			if (i == R2.length - 1) imputedLast = true;
		}
		
		//re-mould them in vectors and matrices
		int[] depth = new int[goodOnes.size()];
		if (imputedFirst) {
			for (int j = 0 ; j < myHFFit.numberOfParams ; j++) myHFFit.params[0][j] = myHFFit.params[goodOnes.get(1)][j];
		}
		if (imputedLast) {
			for (int j = 0 ; j < myHFFit.numberOfParams ; j++) myHFFit.params[goodOnes.size() - 1][j] = myHFFit.params[goodOnes.get(goodOnes.size() - 2)][j];
		}		
		
		//re-mould them in vectors and matrices		
		double[][] data = new double[goodOnes.size()][myHFFit.numberOfParams];		
		for (int i = 0 ; i < depth.length ; i++) {
			depth[i] = goodOnes.get(i);
			for (int j = 0 ; j < myHFFit.numberOfParams ; j++) data[i][j] = myHFFit.params[depth[i]][j];
		}
				
				
		//create a correction function map with gaps
		double[] radius = new double[(int)maths.max(myModes.radius)];
		for (int i = 0 ; i < radius.length ; i++) radius[i] = i;	
		FloatProcessor corrIP = new FloatProcessor(radius.length, myModes.maskingThreshold.length);
		for (int y = 0 ; y < myModes.maskingThreshold.length ; y++) {
			if (goodOnes.contains(y)) {
				double[] nowHFparams = new double[2];				
				myHF.setup(StatUtils.max(myModes.radius), myModes.maskedRadialMinima[y][myModes.radius.length - 1]);
				for (int i = 0 ; i < 2 ; i++) nowHFparams[i] = myHFFit.params[y][i];
				for (int x = 0 ; x < radius.length ; x++) {	
					double myValue2Put = myModes.maskedRadialModes[y][myModes.radius.length - 1];
					try {myValue2Put = myHF.value(x, nowHFparams);} 
					catch (NotStrictlyPositiveException e) {} 
					corrIP.putPixelValue(x, y, myValue2Put);
				}
			}
			else {
				for (int x = 0 ; x < radius.length ; x++) {					
					corrIP.putPixelValue(x, y, 0);
				}
			}
		}
		
//		ImagePlus letsSee = new ImagePlus("nonFiltered", corrIP);
//		letsSee.updateAndDraw();letsSee.show();
		
		//fill in the gaps in the map
		for (int y = 0 ; y < goodOnes.size() - 1 ; y++) {
			
			int lay0 = goodOnes.get(y);
			int lay1 = goodOnes.get(y+1);
			double dist = lay1 - lay0;
			
			if (dist > 1) {
				
				//save neighboring pixel rows
				double[] pa0 = new double[corrIP.getWidth()];
				double[] pa1 = new double[corrIP.getWidth()];
				double[] valueDist = new double[corrIP.getWidth()];
				for (int r = 0 ; r < corrIP.getWidth() ; r++) {
					pa0[r] = corrIP.getPixelValue(r, lay0);
					pa1[r] = corrIP.getPixelValue(r, lay1);
					valueDist[r] = pa1[r] - pa0[r];
				}
				
				//interpolate missing values from neighboring pixel rows by linear interpolation... 
				for (int lay = lay0 + 1 ; lay < lay1 ; lay++) {
					double lilDist = lay - lay0;
					for (int r = 0 ; r < corrIP.getWidth() ; r++) {
						double interpolatedValue = pa0[r] + lilDist / dist * valueDist[r];
						corrIP.putPixelValue(r,lay,interpolatedValue);
					}
				}
			}
		}
				
		//blur the final correction map..
		GaussianBlur myGB = new GaussianBlur();
		double accuracy = 0.01d;
		FloatProcessor blurIP = (FloatProcessor)corrIP.duplicate();
		myGB.blurFloat(blurIP, 20, 20, accuracy);			
		//ImagePlus letsSeeM = new ImagePlus("Filtered",blurIP);
		//letsSeeM.updateAndDraw();letsSeeM.show();
		
		return blurIP;
		
	}
	
	public ImagePlus clipImage(int i, ImagePlus nowTiff, InputOutput.MyFileCollection mFC, MenuWaiter.ClipperMenuReturn mSCM) {
		
		InputOutput jIO = new InputOutput(); 
		ObjectDetector jOD = new ObjectDetector(); 
		RoiHandler roi = new RoiHandler();
		
		ImagePlus outTiff = new ImagePlus();		
		
		String nowGaugePath = null;
		String nowSurfPath = null;
		String[] myGandS = new String[2];
		ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();		
		PolygonRoi[] pRoi = new PolygonRoi[nowTiff.getNSlices()];		
		int xmin = nowTiff.getWidth();
		int xmax = 0;
		int ymin = nowTiff.getHeight();
		int ymax = 0;
		
		int referenceSlice = 1;		
		
		//load gauge and surface files
		if (mSCM.isSoilColumn == true) {			
			if (mSCM.referenceIsTopMostSlice) myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC);			
			if (mSCM.referenceIsSoilSurface) {
				myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC);
				nowSurfPath = mFC.mySurfaceFolder + mFC.pathSep + myGandS[1];
			}
			nowGaugePath = mFC.myInnerCircleFolder + mFC.pathSep + myGandS[0];
			
			//read InnerCircle file
			int versio = jIO.checkInnerCircleFileVersion(nowGaugePath);			
			if (versio == 0) jCO = jIO.readInnerCircleVer0(nowGaugePath);	
			else jCO = jIO.readInnerCircleVer1(nowGaugePath);
			
			pRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, mSCM.clipFromInnerPerimeter);
					
			if (mSCM.referenceIsSoilSurface == true) {
				ImagePlus mySurfs = jIO.openTiff3D(nowSurfPath);	
				referenceSlice = jOD.findMedianSoilSurfacePosition(mySurfs);	
			}
			
			//find min and max values for clip
			for (int j = referenceSlice + mSCM.startAtSlice ; j < referenceSlice + mSCM.stopAtSlice ; j++) {
				
				Polygon p = pRoi[j-2].getNonSplineCoordinates();
				
				// for x			
				int[] nowX = new int[pRoi[j].getNCoordinates()]; 				
				nowX = p.xpoints;
				Arrays.sort(nowX, 0, nowX.length); 
				int xmi = nowX[0] + (int)Math.round(pRoi[j].getXBase());
				int xma = nowX[nowX.length - 1] + (int)Math.round(pRoi[j].getXBase());
				
				// and for y
				int[] nowY = new int[pRoi[j].getNCoordinates()]; 				
				nowY = p.ypoints;
				Arrays.sort(nowY, 0, nowY.length); 
				int ymi = nowY[0] + (int)Math.round(pRoi[j].getYBase());
				int yma = nowY[nowY.length - 1] + (int)Math.round(pRoi[j].getYBase());
				
				// remember the mins and maxes
				if (xmi < xmin) xmin = xmi - mSCM.canvasExceedsBy;
				if (xma > xmax) xmax = xma + mSCM.canvasExceedsBy;
				if (ymi < ymin) ymin = ymi - mSCM.canvasExceedsBy;
				if (yma > ymax) ymax = yma + mSCM.canvasExceedsBy;				
			}			
		} 
		else {			// if the sample is not a soil column..			
			xmin = 0 + mSCM.clipFromCanvasEdge - mSCM.canvasExceedsBy;
			xmax = nowTiff.getWidth() - mSCM.clipFromCanvasEdge + mSCM.canvasExceedsBy;
			ymin = 0 + mSCM.clipFromCanvasEdge - mSCM.canvasExceedsBy;
			ymax = nowTiff.getHeight() - mSCM.clipFromCanvasEdge + mSCM.canvasExceedsBy;			
		}			
		
		//create clip-out-roi
		float[] xpf = {xmin, xmax, xmax, xmin, xmin};
		float[] ypf = {ymin, ymin, ymax, ymax, ymin};	
		PolygonRoi coRoi = new PolygonRoi(xpf, ypf, Roi.POLYGON);
		
		//create simple Roi in case that the image does not contain a soil column
		float[] xps = {xmin + mSCM.canvasExceedsBy, xmax - mSCM.canvasExceedsBy, xmax - mSCM.canvasExceedsBy, xmin + mSCM.canvasExceedsBy, xmin + mSCM.canvasExceedsBy};
		float[] yps = {ymin + mSCM.canvasExceedsBy, ymin + mSCM.canvasExceedsBy, ymax - mSCM.canvasExceedsBy, ymax - mSCM.canvasExceedsBy, ymin + mSCM.canvasExceedsBy};
		PolygonRoi simpleRoi = new PolygonRoi(xps, yps, Roi.POLYGON);
		
		//clip and cut
		ImageStack outStack = null;
		if (!mSCM.preserveOriginialCanvasSize) {
		
			outStack = new ImageStack(xmax - xmin, ymax - ymin);
		
			for (int j = referenceSlice + mSCM.startAtSlice ; j < referenceSlice + mSCM.stopAtSlice ; j++) {
			
				nowTiff.setPosition(j);
				ImageProcessor nowIP = nowTiff.getProcessor();
			
				if (mSCM.isSoilColumn == true) simpleRoi = pRoi[j];
			
				nowIP.setRoi(simpleRoi);
				nowIP.setBackgroundValue(0);
				nowIP.fillOutside(simpleRoi);
			
				nowIP.setRoi(coRoi);
				ImageProcessor cropIP = nowIP.crop();
			
				outStack.addSlice(cropIP);
			}			
		} 
		else {
			
			outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
			
			for (int j = referenceSlice + mSCM.startAtSlice ; j < referenceSlice + mSCM.stopAtSlice ; j++) {
			
				nowTiff.setPosition(j);
				ImageProcessor nowIP = nowTiff.getProcessor();
			
				if (mSCM.isSoilColumn == true) simpleRoi = pRoi[j];
			
				nowIP.setRoi(simpleRoi);
				nowIP.setBackgroundValue(0);
				nowIP.fillOutside(simpleRoi);
			
				outStack.addSlice(nowIP);
			}
		}			 
		
		//add blank canvas to bottom if it is wished..
		if (mSCM.addCanvasExccedance2Bottom == true) {			
			ImageProcessor blankIP = null;
			if (nowTiff.getBitDepth() == 8) blankIP = new ByteProcessor(outStack.getWidth(), outStack.getHeight());
			if (nowTiff.getBitDepth() == 16) blankIP = new ShortProcessor(outStack.getWidth(), outStack.getHeight());
			for (int j = 0 ; j < mSCM.canvasExceedsBy ; j++) {
				outStack.addSlice(blankIP);	
			}
		}
		
		outTiff.setStack(outStack);
		
		//outTiff.updateAndDraw();
		//outTiff.show();
		
		return outTiff;
		
	}
	
//	public ImagePlus correct4SteelBeamHardeningDEPRECACTED(ImagePlus nowTiff, ObjectDetector.EggShapedColCoords3D jCO, ImageProcessor lCM, ImageProcessor uCM, int standardRadius) {
//		
//		//init output variables
//		ImagePlus outTiff = new ImagePlus();
//		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
//		
//		//calculate number of angles needed to paint the whole canvas..
//		double spacing = 2 * Math.PI / jCO.anglesChecked;		
//		float[] sR = new float[standardRadius]; 
//		float[] sA = new float[jCO.anglesChecked + 1];
//		
//		//create array for standard radius
//		for (int i = 0 ; i < standardRadius ; i++) sR[i] = i + 1;
//		for (int i = 0 ; i < jCO.anglesChecked + 1 ; i++) sA[i] = i;		 
//
//		//sweep over slices
//		for (int i = 0 ; i < nowTiff.getNSlices() ; i++) {
//			
//			IJ.showStatus("Correcting for beam hardening in slice #" + (i + 1) + "/" + nowTiff.getNSlices());
//			
//			//set image to next slice
//			nowTiff.setPosition(i+1);
//			ImageProcessor nowIP = nowTiff.getProcessor();
//
//			//calculate normalized correction functions			
//			double loRef = lCM.getPixelValue(0, i); 
//			double hiRef = uCM.getPixelValue(0, i); 
//			 
//			//float[] ftotLower = correctionFunctionCalculator(i, myLowerFits, sR);	
//			//float[] ftotUpper = correctionFunctionCalculator(i, myUpperFits, sR);
//			double[] loFunc = new double[standardRadius];				
//			double[] hiFunc = new double[standardRadius];
//			double[] deltaFunc = new double[standardRadius];
//			for (int r = 0 ; r < standardRadius ; r++) {
//				loFunc[r] = lCM.getPixelValue(r, i);
//				hiFunc[r] = uCM.getPixelValue(r, i);
//				deltaFunc[r] = hiFunc[r] - loFunc[r];				
//			}
//			
//			//create an angle map
//			double[] nowAngle = new double[jCO.anglesChecked];
//			for (int j = 0 ; j < jCO.anglesChecked ; j++) {
//				
//				double x = jCO.xID[i][j];
//				double y = jCO.yID[i][j];
//				
//				double dx = x - jCO.xmid[i];
//				double dy = y - jCO.ymid[i];
//				double nowRadius = Math.sqrt(dx*dx + dy*dy);
//				
//				if (j < 36) nowAngle[j] = Math.acos(dy / nowRadius);				
//				if (j >= 36) nowAngle[j] = 2 * Math.PI - Math.acos(dy / nowRadius);
//				
//			}
//			
//			//do the correction
//			for (int x = 0 ; x < nowIP.getWidth() ; x++) {
//				for (int y = 0 ; y < nowIP.getHeight() ; y++) {			
//					
//					//get distance of pixel to center
//					float dx = x - (float)jCO.xmid[i];
//					float dy = y - (float)jCO.ymid[i];
//					float nowRadius = (float)Math.sqrt((double)(dx*dx) + (double)(dy*dy));
//					
//					//get angle
//					double alpha = 0; 
//					if (dx >= 0) alpha = Math.acos(dy / nowRadius);				
//					if (dx < 0) alpha = 2 * Math.PI - Math.acos(dy / nowRadius);			
//															
//					//get radius at this angle
//					double dx0 = jCO.xID[i][0] - jCO.xmid[i];
//					double dy0 = jCO.yID[i][0] - jCO.ymid[i];
//					double radiusAtThisAngle = Math.sqrt((dx0*dx0) + (dy0*dy0));;
//					for (int j = 0 ; j < jCO.anglesChecked ; j++) {						
//						if (nowAngle[j] > alpha) {
//							
//							double weight = (nowAngle[j] - alpha) / spacing;							
//							
//							int previousAngleId = j - 1;
//							if (j == 0) previousAngleId = 71;
//							
//							double dx1 = jCO.xID[i][previousAngleId] - jCO.xmid[i];
//							double dx2 = jCO.xID[i][j] - jCO.xmid[i];
//							double dy1 = jCO.yID[i][previousAngleId] - jCO.ymid[i];
//							double dy2 = jCO.yID[i][j] - jCO.ymid[i];
//							double r1 = Math.sqrt((dx1*dx1) + (dy1*dy1));
//							double r2 = Math.sqrt((dx2*dx2) + (dy2*dy2));
//							
//							radiusAtThisAngle = weight * r2 + (1 - weight) * r1;
//							
//							break;
//						}
//					}
//					
//					//check if pixel is within column and if yes apply correction
//					if (nowRadius < radiusAtThisAngle) {
//						
//						int renormalizedRadius = (int)Math.floor(nowRadius / radiusAtThisAngle * standardRadius);
//						double loFactor = loFunc[renormalizedRadius];
//						double hiFactor = hiFunc[renormalizedRadius];
//						double deltaFactor = hiFunc[renormalizedRadius];
//						
//						float mygray = nowIP.getPixelValue(x, y);											
//						double newgray = mygray * loRef /loFactor;
//						
//						//double base = mygray - loRef;
//						//double extra = base / deltaFactor * 20000;
//						//double newgray = 10000 + extra;
//						
//						nowIP.putPixel(x, y, (int)newgray);
//						
//					}			
//				}	
//		
//			}
//			
//			//nowTiff.updateAndDraw();
//			//nowTiff.show();
//
//			outStack.addSlice(nowIP);
//		}
//		
//		outTiff.setStack(outStack);
//		
//		//outTiff.updateAndDraw();
//		//outTiff.show();
//		
//		return outTiff;		
//		
//	}
	
	public ImagePlus correct4BeamHardening(ImagePlus nowTiff, ObjectDetector.ColCoordsEssentials3D jCO, ImagePlus correctionMaps, int standardRadius) {
		
		//init output variables
		ImagePlus outTiff = new ImagePlus();
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		//extract individual correctionMaps
		correctionMaps.setSlice(1);
		ImageProcessor lCM = correctionMaps.getProcessor();
		correctionMaps.setSlice(2);
		ImageProcessor uCM = correctionMaps.getProcessor();		
		
		//calculate number of angles needed to paint the whole canvas..
		double spacing = 2 * Math.PI / jCO.anglesChecked;		
		float[] sR = new float[standardRadius]; 
		float[] sA = new float[jCO.anglesChecked + 1];
		
		//create array for standard radius
		for (int i = 0 ; i < standardRadius ; i++) sR[i] = i + 1;
		for (int i = 0 ; i < jCO.anglesChecked + 1 ; i++) sA[i] = i;		 

		//sweep over slices
		for (int i = 0 ; i < nowTiff.getNSlices() ; i++) {
			
			IJ.showStatus("Correcting for beam hardening in slice #" + (i + 1) + "/" + nowTiff.getNSlices());
			
			//set image to next slice
			nowTiff.setPosition(i+1);
			ImageProcessor nowIP = nowTiff.getProcessor();

			//calculate normalized correction functions			
			double loRef = lCM.getPixelValue(0, i); 
			double hiRef = uCM.getPixelValue(0, i); 
			 
			//float[] ftotLower = correctionFunctionCalculator(i, myLowerFits, sR);	
			//float[] ftotUpper = correctionFunctionCalculator(i, myUpperFits, sR);
			double[] loFunc = new double[standardRadius];				
			double[] hiFunc = new double[standardRadius];
			double[] deltaFunc = new double[standardRadius];
			for (int r = 0 ; r < standardRadius ; r++) {
				loFunc[r] = lCM.getPixelValue(r, i);
				hiFunc[r] = uCM.getPixelValue(r, i);
				deltaFunc[r] = hiFunc[r] - loFunc[r];				
			}
			
			//create an angle map
			double[] nowAngle = new double[jCO.anglesChecked];
			for (int j = 0 ; j < jCO.anglesChecked ; j++) {
				
				double x = jCO.xID[i][j];
				double y = jCO.yID[i][j];
				
				double dx = x - jCO.xmid[i];
				double dy = y - jCO.ymid[i];
				double nowRadius = Math.sqrt(dx*dx + dy*dy);
				
				if (j < 36) nowAngle[j] = Math.acos(dy / nowRadius);				
				if (j >= 36) nowAngle[j] = 2 * Math.PI - Math.acos(dy / nowRadius);
				
			}
			
			//do the correction
			for (int x = 0 ; x < nowIP.getWidth() ; x++) {
				for (int y = 0 ; y < nowIP.getHeight() ; y++) {			
					
					//get distance of pixel to center
					float dx = x - (float)jCO.xmid[i];
					float dy = y - (float)jCO.ymid[i];
					float nowRadius = (float)Math.sqrt((double)(dx*dx) + (double)(dy*dy));
					
					//get angle
					double alpha = 0; 
					if (dx >= 0) alpha = Math.acos(dy / nowRadius);				
					if (dx < 0) alpha = 2 * Math.PI - Math.acos(dy / nowRadius);			
															
					//get radius at this angle
					double dx0 = jCO.xID[i][0] - jCO.xmid[i];
					double dy0 = jCO.yID[i][0] - jCO.ymid[i];
					double radiusAtThisAngle = Math.sqrt((dx0*dx0) + (dy0*dy0));;
					for (int j = 0 ; j < jCO.anglesChecked ; j++) {						
						if (nowAngle[j] > alpha) {
							
							double weight = (nowAngle[j] - alpha) / spacing;							
							
							int previousAngleId = j - 1;
							if (j == 0) previousAngleId = 71;
							
							double dx1 = jCO.xID[i][previousAngleId] - jCO.xmid[i];
							double dx2 = jCO.xID[i][j] - jCO.xmid[i];
							double dy1 = jCO.yID[i][previousAngleId] - jCO.ymid[i];
							double dy2 = jCO.yID[i][j] - jCO.ymid[i];
							double r1 = Math.sqrt((dx1*dx1) + (dy1*dy1));
							double r2 = Math.sqrt((dx2*dx2) + (dy2*dy2));
							
							radiusAtThisAngle = weight * r2 + (1 - weight) * r1;
							
							break;
						}
					}
					
					//check if pixel is within column and if yes apply correction
					if (nowRadius < radiusAtThisAngle) {
						
						int renormalizedRadius = (int)Math.floor(nowRadius / radiusAtThisAngle * standardRadius);
						double loFactor = loFunc[renormalizedRadius];
						double hiFactor = hiFunc[renormalizedRadius];
						double deltaFactor = hiFunc[renormalizedRadius];
						
						float mygray = nowIP.getPixelValue(x, y);											
						double newgray = mygray * loRef /loFactor;
						
						//double base = mygray - loRef;
						//double extra = base / deltaFactor * 20000;
						//double newgray = 10000 + extra;
						
						nowIP.putPixel(x, y, (int)newgray);
						
					}			
				}	
		
			}
			
			//nowTiff.updateAndDraw();
			//nowTiff.show();

			outStack.addSlice(nowIP);
		}
		
		outTiff.setStack(outStack);
		
		//outTiff.updateAndDraw();
		//outTiff.show();
		
		return outTiff;		
		
	}
	
//	public ImagePlus correct4PVCBeamHardeningDEPRECATED(ImagePlus nowTiff, ImageProcessor blurIP, ImageProcessor gammaIP, ObjectDetector.ColCoords3D jCO, MenuWaiter.BeamDeHardeningReturn mBDH) {
//		
//		RoiHandler roi = new RoiHandler();
//		
//		//init output variables
//		ImagePlus outTiff = new ImagePlus();
//		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
//		
//		//calculate number of angles needed to paint the whole canvas..
//		double spacing = 2 * Math.PI / mBDH.anglesChecked;
//		int standardRadius = blurIP.getWidth();
//		float[] sR = new float[standardRadius]; 
//		float[] sA = new float[mBDH.anglesChecked + 1];
//		
//		//create array for standardradius
//		for (int i = 0 ; i < standardRadius ; i++) sR[i] = i + 1;
//		for (int i = 0 ; i < mBDH.anglesChecked + 1 ; i++) sA[i] = i;	
//		
//		//load polygon roi of inner perimeter
//		PolygonRoi[] nowRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 2);
//
//		//sweep over slices
//		for (int i = 0 ; i < nowTiff.getNSlices() ; i++) {
//			
//			IJ.showStatus("Correcting for beam hardening in slice #" + (i + 1) + "/" + nowTiff.getNSlices());
//			
//			//extract X and Y of pRoi		
//			nowRoi[i].fitSpline(mBDH.anglesChecked);
//			Polygon myPoly = nowRoi[i].getPolygon();			
//			int xID[] = myPoly.xpoints;
//			int yID[] = myPoly.ypoints;
//			
//			//set image to next slice
//			nowTiff.setPosition(i+1);
//			ImageProcessor nowIP = nowTiff.getProcessor();
//			
//			//get matrix brightness reference values for this depth			
//			double[] corrFunc = new double[standardRadius];
//			double nowReference = blurIP.getPixelValue(standardRadius - 1, i);
//			for (int r = 0 ; r < standardRadius ; r++) corrFunc[r] = nowReference / blurIP.getPixelValue(r, i);
//			
//			//get air-phase gamma reference values for this depth			
//			double[] gammaFunc = new double[standardRadius];
//			double nowRefGamma = gammaIP.getPixelValue(standardRadius - 1, i);
//			double gammaDelta = nowReference - nowRefGamma;
//			for (int r = 0 ; r < standardRadius ; r++) gammaFunc[r] = gammaDelta / (nowReference - gammaIP.getPixelValue(r, i));
//					
//			//do the correction
//			for (int x = 0 ; x < nowIP.getWidth() ; x++) {
//				for (int y = 0 ; y < nowIP.getHeight() ; y++) {			
//					
//					//get distance of pixel to center
//					float dx = x - (float)jCO.xmid[i];
//					float dy = y - (float)jCO.ymid[i];
//					float nowRadius = (float)Math.sqrt((double)(dx*dx) + (double)(dy*dy));
//					
//					//get angle
//					double alpha = Math.atan(dy/dx);
//					if (dx < 0 & dy >= 0) alpha = 2 * Math.PI + alpha;
//					if (dy < 0) alpha = Math.PI + alpha;
//															
//					//get radius at this angle
//					double dx0 = xID[0] - jCO.xmid[i];
//					double dy0 = yID[0] - jCO.ymid[i];
//					double radiusAtThisAngle = Math.sqrt((dx0*dx0) + (dy0*dy0));;
//					for (double j = 0 ; j < mBDH.anglesChecked ; j++) {
//						double nowAngle  = 2 * j / mBDH.anglesChecked * Math.PI;
//						if (nowAngle > alpha) {
//							
//							double weight = (nowAngle - alpha) / spacing;							
//							
//							double dx1 = xID[(int)j-1] - jCO.xmid[i];
//							double dx2 = xID[(int)j] - jCO.xmid[i];
//							double dy1 = yID[(int)j-1] - jCO.ymid[i];
//							double dy2 = yID[(int)j] - jCO.ymid[i];
//							double r1 = Math.sqrt((dx1*dx1) + (dy1*dy1));
//							double r2 = Math.sqrt((dx2*dx2) + (dy2*dy2));
//							
//							radiusAtThisAngle = weight * r2 + (1 - weight) * r1;
//							
//							break;
//						}
//					}
//					
//					//check if pixel is within column and if yes apply correction
//					double corrFactor = 1;
//					double gammaFactor = 1;
//					if (nowRadius < radiusAtThisAngle) {						
//						int renormalizedRadius = (int)Math.floor(nowRadius / radiusAtThisAngle * standardRadius);
//						corrFactor = corrFunc[renormalizedRadius];
//						gammaFactor = gammaFunc[renormalizedRadius];
//					}		
//					
//					//apply beam hardening correction
//					float mygray = nowIP.getPixelValue(x, y);					
//					double newgray = (int)Math.round(corrFactor * mygray);
//					
//					//also apply air-phase gamma correction
//					int newestgray = (int)Math.round(nowReference - (nowReference - newgray) * gammaFactor); 
//					
//					nowIP.putPixel(x, y, newestgray);
//				}	
//		
//			}
//	
//			outStack.addSlice(nowIP);
//		}
//		
//		outTiff.setStack(outStack);
//		
//		return outTiff;		
//		
//	}
	
	public float[] correctionFunctionCalculator(int i, FitStuff.FittingResults myFR, float[] sR) {
		
		float[] ftot = new float[sR.length];
		int standardRadius = sR.length;
		
		double a = myFR.params[i][0];
		double b = myFR.params[i][1];
		double c = myFR.params[i][2];
		double r = myFR.params[i][3];
		double dy = myFR.params[i][4];		
		double f1, f2, f3;
		
		for (int j = 0 ; j < standardRadius ; j++) {
			f1 = dy / Math.exp(a * sR[j] - a) + r;
			f2 = dy / Math.exp(b * sR[j] - b) + r;
			f3 = dy / Math.exp(c * sR[j] - c) + r;
			ftot[standardRadius - j - 1] = (float)((f1 + f2 + f3) / 3);  //switch vector around!!!! important!!
		}
		
		return ftot;
	}
	
	public ImagePlus putColumnUprightInCenter(ImagePlus nowTiff, ObjectDetector.ColCoords3D prelimCC, MenuWaiter.ColumnFinderMenuReturn jCFS) {
		
		ObjectDetector jOD = new ObjectDetector();		
		Slicer sly = new Slicer();
		Median jMed = new Median();
		
		ObjectDetector.ColCoords3D newCC = jOD.new ColCoords3D();
		
		ImagePlus straightTiff = new ImagePlus();
		ImagePlus outTiff = new ImagePlus();
		
		//by-pass the cutting step when re-measuring column wall..
		//outTiff = nowTiff;
		
		if (prelimCC.tiltTotal > 0.01) {   //if tilting angle is too big then put column straight
			
			IJ.showStatus("Putting column straight and moving it to center of canvas, step 1 of 4...");
			
			ImagePlus boTiff0 = new ImagePlus();
			ImageStack boStack0 = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
			//find tilting angle relative to XY-plane
			double alpha = prelimCC.tiltInXZ;
			double beta = prelimCC.tiltInYZ;
			double dz = 1000;
			
			double dx = Math.tan(alpha) * dz;
			double dy = Math.tan(beta) * dz;
			double gamma = Math.atan(dx/dy);
			
			//correct for ambiguity in atan response			
			if (dy < 0) gamma = gamma + Math.PI; 
			
			//rotate so that the tilting is only in y-direction
			for (int i = 1 ; i < nowTiff.getNSlices() + 1 ; i++) {
				
				nowTiff.setPosition(i);
				ImageProcessor nowIP = nowTiff.getProcessor();
				nowIP.setInterpolate(true);
				double gammaInDegree = gamma * 360 / 2 / Math.PI ;
				nowIP.rotate(gammaInDegree + 90);
				
				boStack0.addSlice(nowIP);
			}
			
			boTiff0.setStack(boStack0);			
						
			//reslice
			ImagePlus vertiTiff = sly.reslice(boTiff0);
						
			//free memory ..
			boTiff0.flush();
			IJ.freeMemory();IJ.freeMemory();
			
			IJ.showStatus("Putting column straight and moving it to center of canvas, step 2 of 4...");
			
			//put column upright
			ImagePlus boTiff = new ImagePlus();
			ImageStack boStack = new ImageStack(vertiTiff.getWidth(), vertiTiff.getHeight());			
			for (int i = 1 ; i < vertiTiff.getNSlices() + 1 ; i++) {
				
				vertiTiff.setPosition(i);
				ImageProcessor nowIP = vertiTiff.getProcessor();
				nowIP.setInterpolate(true);
				double deltaInDegree = prelimCC.tiltTotal * 360 / 2 / Math.PI ;
				nowIP.rotate(-deltaInDegree);
				
				boStack.addSlice(nowIP);
			}

			//free memory ..
			vertiTiff.flush();
			IJ.freeMemory();IJ.freeMemory();
			
			//transfer images to boTiff
			boTiff.setStack(boStack);	
			
			//reslice back to XY-plane view
			ImagePlus straightTiff0 = sly.reslice(boTiff);
			
			//free memory from boTiff again..
			boTiff.flush();
			IJ.freeMemory();IJ.freeMemory();
			
			IJ.showStatus("Putting column straight and moving it to center of canvas, step 3 of 4...");
			
			//rotate back to original orientation
			ImageStack straightStack = new ImageStack(straightTiff0.getWidth(), straightTiff0.getHeight());
			for (int i = 1 ; i < straightTiff0.getNSlices() + 1 ; i++) {
				
				straightTiff0.setPosition(i);
				ImageProcessor nowIP = straightTiff0.getProcessor();
				nowIP.setInterpolate(true);
				double gammaInDegree = gamma * 360 / 2 / Math.PI ;
				nowIP.rotate(-gammaInDegree - 90);
				
				straightStack.addSlice(nowIP);
			}
			
			//free memory from straightTiff0 again..
			straightTiff0.flush();
			IJ.freeMemory();IJ.freeMemory();
			
			straightTiff.setStack(straightStack);

		}
		else {
			straightTiff = nowTiff;
		}
				
/*		straightTiff.draw();
		straightTiff.show();*/
		
		//re-find column outlines
		boolean look4PreciseCoords = false;
		newCC = jOD.findOrientationOfPVCOrAluColumn(straightTiff, jCFS, look4PreciseCoords);
		
		//move column into the center of the canvas and cut out unnecessary parts of the canvas
		IJ.showStatus("Putting column straight and moving it to center of canvas, step 4 of 4...");
		double rim = 25;
		double mRadius = jMed.evaluate(newCC.outerMajorRadius);
		double toBeLeft = mRadius + rim;
		
		//check if cut-out Roi is not larger than the image..
		if (2 * toBeLeft > straightTiff.getWidth()) toBeLeft = straightTiff.getWidth() / 2;
		if (2 * toBeLeft > straightTiff.getHeight()) toBeLeft = straightTiff.getHeight() / 2;
			
		//do the cutting
		double diameter = 2 * toBeLeft;
		ImageStack outStack = new ImageStack((int)Math.round(diameter), (int)Math.round(diameter));
		
		//Boolean cutSuccessful = true;
		double XC = straightTiff.getWidth() / 2;
		double YC = straightTiff.getHeight() / 2;
		
		double dx = XC - jMed.evaluate(newCC.xmid);
		double dy = YC - jMed.evaluate(newCC.ymid);
		for (int i = 1 ; i < straightTiff.getNSlices() + 1 ; i++) {
			
			IJ.showStatus("Moving column to center of canvas " + i + "/" + straightTiff.getNSlices());
			
			straightTiff.setPosition(i);
			ImageProcessor nowIP = straightTiff.getProcessor();
			nowIP.setInterpolate(true);
			
			nowIP.translate(dx, dy);
			
			Roi cutRoi = new Roi(XC - toBeLeft, YC - toBeLeft, (int)Math.round(diameter), (int)Math.round(diameter));			
			nowIP.setRoi(cutRoi);
			ImageProcessor cutIP = nowIP.crop();
			
			outStack.addSlice(cutIP);
		}
		
		//free memory from straightTiff again..
		straightTiff.flush();
		IJ.freeMemory();IJ.freeMemory();
		
		outTiff.setStack(outStack);		
		
		//outTiff.updateAndDraw();
		//outTiff.show();		
		
		return outTiff;
		
	}
	
	/*public ImagePlus removeRingArtifacts(ImagePlus nowTiff, String innerCirclePath, MenuWaiter.RingArtifactRemoverOptions rARO) {
				
		ObjectDetector jOD = new ObjectDetector();
		InputOutput jIO = new InputOutput();
		RoiHandler roi = new RoiHandler();
		
		ImagePlus outTiff = new ImagePlus();
		PolygonRoi[] pRoi = new PolygonRoi[nowTiff.getNSlices()];
		
		//read gauge file
		if (rARO.useInnerCircle) {
			ObjectDetector.ColumnCoordinates jCO = jIO.readGaugeFile(innerCirclePath);	
			pRoi = roi.makeMeAPolygonRoiStack("inner", "tight", jCO, 0);
		}
		
		double[] axisOfEvil = jOD.detectRingArtifacts(nowTiff, rARO); //returns tow firstOrder polynomials describing the location of the centers of the artifacts: x = az + b; y = cz + d where a,b,c,d == axisOfEvil[0..3]
		
	
		double[] XisazPLUSbYisczPLUSd = {0, 0, 0, 0,};
		double xCenter = (float)nowTiff.getWidth() / 2;
		double yCenter = (float)nowTiff.getHeight() / 2;
		int angles = 12;
		double aincr = 2 * Math.PI / angles;		
				
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			
			//get informed about column location
			double diam = nowTiff.getWidth();
			if (pRoi[z] != null) pRoi[z].getFloatWidth();
			else if (nowTiff.getWidth() > nowTiff.getHeight()) diam = nowTiff.getHeight();
			int checkPosis = (int)Math.floor(diam / 2);
			
			//scan layer z
			nowTiff.setPosition(z+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			//smooth layer z						
			myRF.rank(nowIP, 2, RankFilters.MEDIAN);
			nowIP.sharpen();		
						
			//sample radial profiles
			double[][] radProf = new double[angles][checkPosis-1];
			int angCount = 0;
			for (double alpha = 0 ; alpha < 2 * Math.PI - aincr; alpha += aincr) {
				for (double r = 1 ; r < checkPosis ; r++) {
					int x = (int) Math.round(xCenter + Math.sin(alpha) * r);
					int y = (int) Math.round(yCenter - Math.cos(alpha) * r);
					radProf[angCount][(int)Math.round(r) - 1] = nowIP.getPixel(x, y);
				}
				angCount++;
			}
			

			
		}
		
		return outTiff;
		
	}*/
	
	
	///////////////////////////////////////////////
	// binarize
	///////////////////////////////////////////////
	
	public int[] multiOtsu(InputOutput.MyFileCollection mFC, int[] nowHist, MenuWaiter.MultiRegionSegmentation mRS) {
		
		//following Liao et al. (2001), Journal of Information Science and Engineering, 17, 713-727
		
		DisplayThings disp = new DisplayThings();
		RollerCaster rC = new RollerCaster();
		HistogramStuff hist = new HistogramStuff();
		
		int[] thresholds = new int[mRS.numberOfThresholds];
		
		//if histogram is in 16 bit, scale to 8 bit
		int[] eightHist = new int[256];
		if (mRS.convert2EightBit) hist.convert16to8BitHistogram(nowHist, mRS.lowerBoundary, mRS.upperBoundary);
		else eightHist = nowHist;
		
		ImageProcessor P = calcLookUpTableP(eightHist);
		ImageProcessor S = calcLookUpTableS(P);
		
/*		ImagePlus test0 = new ImagePlus("P", P);		
		test0.draw();test0.show();
		
		ImagePlus test1 = new ImagePlus("S", S);		
		test1.draw();test1.show();*/
			
		//find maximum...
		double ojMax = 0;
		switch(mRS.numberOfThresholds) {
			
			case 1: 
				
				for (int i = 0 ; i < 256 - mRS.numberOfThresholds ; i++) {
					
					double objectiveFunction = S.getPixelValue(i, 0) * S.getPixelValue(i, 0) / P.getPixelValue(i, 0);					
					objectiveFunction += S.getPixelValue(255, i) * S.getPixelValue(255, i) / P.getPixelValue(255, i);
			
						
					if (objectiveFunction > ojMax) {
						ojMax = objectiveFunction;
						thresholds[0] = i;
					}
						
				}
				break;
			
			case 2: 
				
				for (int i = 1 ; i < 256 - mRS.numberOfThresholds ; i++) {
					for (int j = i + 1 ; j < 256 - 1 ; j++) {										
												
							
						double objectiveFunction = S.getPixelValue(i, 0) * S.getPixelValue(i, 0) / P.getPixelValue(i, 0);
						objectiveFunction += S.getPixelValue(j, i) * S.getPixelValue(j, i) / P.getPixelValue(j, i);
						objectiveFunction += S.getPixelValue(255, j) * S.getPixelValue(255, j) / P.getPixelValue(255, j);
						
						if (objectiveFunction > ojMax) {
							ojMax = objectiveFunction;
							thresholds[0] = i;
							thresholds[1] = j;
						}
							
					}
				}
				break;
				
			case 3: 
			
				for (int i = 0 ; i < 256 - mRS.numberOfThresholds ; i++) {
					for (int j = i + 1 ; j < 256 - 2 ; j++) {
						for (int k = j + 1; k < 256 - 1 ; k++) {						
						
							double objectiveFunction = S.getPixelValue(i, 0) * S.getPixelValue(i, 0) / P.getPixelValue(i, 0);
							objectiveFunction += S.getPixelValue(j, i) * S.getPixelValue(j, i) / P.getPixelValue(j, i);
							objectiveFunction += S.getPixelValue(k, j) * S.getPixelValue(k, j) / P.getPixelValue(k, j);
							objectiveFunction += S.getPixelValue(255, k) * S.getPixelValue(255, k) / P.getPixelValue(255, k);	
							
							if (objectiveFunction > ojMax) {
								ojMax = objectiveFunction;
								thresholds[0] = i;
								thresholds[1] = j;
								thresholds[2] = k;
							}
						}						
					}
				}
				break;
				
			case 4: 
				
				for (int i = 0 ; i < 256 - mRS.numberOfThresholds ; i++) {
					for (int j = i + 1; j < 256 - 3 ; j++) {
						for (int k = j + 1; k < 256 - 2 ; k++) {
							for (int m = k + 1; m < 256 - 1 ; m++) {

								double objectiveFunction = S.getPixelValue(i, 0) * S.getPixelValue(i, 0) / P.getPixelValue(i, 0);
								objectiveFunction += S.getPixelValue(j, i) * S.getPixelValue(j, i) / P.getPixelValue(j, i);
								objectiveFunction += S.getPixelValue(k, j) * S.getPixelValue(k, j) / P.getPixelValue(k, j);
								objectiveFunction += S.getPixelValue(m, k) * S.getPixelValue(m, k) / P.getPixelValue(m, k);	
								objectiveFunction += S.getPixelValue(255, m) * S.getPixelValue(255, m) / P.getPixelValue(255, m);
								
								if (objectiveFunction > ojMax) {
									ojMax = objectiveFunction;
									thresholds[0] = i;
									thresholds[1] = j;
									thresholds[2] = k;
									thresholds[3] = m;
								}
							}
						}						
					}
				}
				break;
				
			case 5: 
			
				for (int i = 0 ; i < 256 - mRS.numberOfThresholds ; i++) {
					for (int j = i + 1; j < 256 - 4 ; j++) {
						for (int k = j + 1; k < 256 - 3 ; k++) {
							for (int m = k + 1; m < 256 - 2 ; m++) {
								for (int n = m + 1 ; n < 256 -1 ; n++) {
							
									double objectiveFunction = S.getPixelValue(i, 0) * S.getPixelValue(i, 0) / P.getPixelValue(i, 0);
									objectiveFunction += S.getPixelValue(j, i) * S.getPixelValue(j, i) / P.getPixelValue(j, i);
									objectiveFunction += S.getPixelValue(k, j) * S.getPixelValue(k, j) / P.getPixelValue(k, j);
									objectiveFunction += S.getPixelValue(m, k) * S.getPixelValue(m, k) / P.getPixelValue(m, k);	
									objectiveFunction += S.getPixelValue(n, m) * S.getPixelValue(n, m) / P.getPixelValue(n, m);
									objectiveFunction += S.getPixelValue(255, n) * S.getPixelValue(255, n) / P.getPixelValue(255, n);
									
									if (objectiveFunction > ojMax) {
										ojMax = objectiveFunction;
										thresholds[0] = i;
										thresholds[1] = j;
										thresholds[2] = k;
										thresholds[3] = m;
										thresholds[4] = n;
									}
								}
							}
						}						
					}
				}				
				break;
			
		}

		//double[] x = new double[256]; for (int i = 0 ; i < x.length ; i++) x[i] = i;
		//disp.plotHistogram(x, rC.castInt2Double(eightHist), "test", "histogram classes", "frequency");
		
		return thresholds;
		
	}
	
	public int[] multiMaxEntropy(InputOutput.MyFileCollection mFC, int[] nowHist, MenuWaiter.MultiRegionSegmentation mRS) {
		
		//following S. Schlueter's code from QuantIM
		
		DisplayThings disp = new DisplayThings();
		RollerCaster rC = new RollerCaster();
		HistogramStuff hist = new HistogramStuff();
		
		int[] thresholds = new int[mRS.numberOfThresholds];
		
		//if histogram is in 16 bit, scale to 8 bit
		int[] eightHist = new int[256];
		if (mRS.convert2EightBit) hist.convert16to8BitHistogram(nowHist, mRS.lowerBoundary, mRS.upperBoundary);
		else eightHist = nowHist;
		
		ImageProcessor P = calcLookUpTableP(eightHist);
		ImageProcessor S = calcLookUpTableS4MaxEntro(P);
		
/*		ImagePlus test0 = new ImagePlus("P", P);		
		test0.draw();test0.show();
		
		ImagePlus test1 = new ImagePlus("S", S);		
		test1.draw();test1.show();*/
		
		//find maximum...
		double ojMax = 0;
		switch(mRS.numberOfThresholds) {
			
			case 1: 
				
				for (int i = 0 ; i < 256 - mRS.numberOfThresholds ; i++) {

					double objectiveFunction = Math.log(P.getPixelValue(i, 0)) - S.getPixelValue(i, 0) / P.getPixelValue(i, 0);					
					objectiveFunction += Math.log(P.getPixelValue(255, i)) - S.getPixelValue(255, i) / P.getPixelValue(255, i);
			
						
					if (objectiveFunction > ojMax) {
						ojMax = objectiveFunction;
						thresholds[0] = i;
					}
						
				}
				break;
			
			case 2: 
				
				for (int i = 1 ; i < 256 - mRS.numberOfThresholds ; i++) {
					for (int j = i + 1 ; j < 256 - 1 ; j++) {										
												
							
						double objectiveFunction = Math.log(P.getPixelValue(i, 0)) - S.getPixelValue(i, 0) / P.getPixelValue(i, 0);
						objectiveFunction += Math.log(P.getPixelValue(j, i)) - S.getPixelValue(j, i) / P.getPixelValue(j, i);
						objectiveFunction += Math.log(P.getPixelValue(255, j)) - S.getPixelValue(255, j) / P.getPixelValue(255, j);
						
						if (objectiveFunction > ojMax) {
							ojMax = objectiveFunction;
							thresholds[0] = i;
							thresholds[1] = j;
						}
							
					}
				}
				break;
				
			case 3: 
			
				for (int i = 0 ; i < 256 - mRS.numberOfThresholds ; i++) {
					for (int j = i + 1 ; j < 256 - 2 ; j++) {
						for (int k = j + 1; k < 256 - 1 ; k++) {						
						
							double objectiveFunction = Math.log(P.getPixelValue(i, 0)) - S.getPixelValue(i, 0) / P.getPixelValue(i, 0);
							objectiveFunction += Math.log(P.getPixelValue(j, i)) - S.getPixelValue(j, i) / P.getPixelValue(j, i);
							objectiveFunction += Math.log(P.getPixelValue(k, j)) - S.getPixelValue(k, j) / P.getPixelValue(k, j);
							objectiveFunction += Math.log(P.getPixelValue(255, k)) - S.getPixelValue(255, k) / P.getPixelValue(255, k);	
							
							if (objectiveFunction > ojMax) {
								ojMax = objectiveFunction;
								thresholds[0] = i;
								thresholds[1] = j;
								thresholds[2] = k;
							}
						}						
					}
				}
				break;
				
			case 4: 
				
				for (int i = 0 ; i < 256 - mRS.numberOfThresholds ; i++) {
					for (int j = i + 1; j < 256 - 3 ; j++) {
						for (int k = j + 1; k < 256 - 2 ; k++) {
							for (int m = k + 1; m < 256 - 1 ; m++) {

								double objectiveFunction = Math.log(P.getPixelValue(i, 0)) - S.getPixelValue(i, 0) / P.getPixelValue(i, 0);
								objectiveFunction += Math.log(P.getPixelValue(j, i)) - S.getPixelValue(j, i) / P.getPixelValue(j, i);
								objectiveFunction += Math.log(P.getPixelValue(k, j)) - S.getPixelValue(k, j) / P.getPixelValue(k, j);
								objectiveFunction += Math.log(P.getPixelValue(m, k)) - S.getPixelValue(m, k) / P.getPixelValue(m, k);	
								objectiveFunction += Math.log(P.getPixelValue(255, m)) - S.getPixelValue(255, m) / P.getPixelValue(255, m);
								
								if (objectiveFunction > ojMax) {
									ojMax = objectiveFunction;
									thresholds[0] = i;
									thresholds[1] = j;
									thresholds[2] = k;
									thresholds[3] = m;
								}
							}
						}						
					}
				}
				break;
				
			case 5: 
			
				for (int i = 0 ; i < 256 - mRS.numberOfThresholds ; i++) {
					for (int j = i + 1; j < 256 - 4 ; j++) {
						for (int k = j + 1; k < 256 - 3 ; k++) {
							for (int m = k + 1; m < 256 - 2 ; m++) {
								for (int n = m + 1 ; n < 256 -1 ; n++) {
							
									double objectiveFunction = Math.log(P.getPixelValue(i, 0)) - S.getPixelValue(i, 0) / P.getPixelValue(i, 0);
									objectiveFunction += Math.log(P.getPixelValue(j, i)) - S.getPixelValue(j, i) / P.getPixelValue(j, i);
									objectiveFunction += Math.log(P.getPixelValue(k, j)) - S.getPixelValue(k, j) / P.getPixelValue(k, j);
									objectiveFunction += Math.log(P.getPixelValue(m, k)) - S.getPixelValue(m, k) / P.getPixelValue(m, k);	
									objectiveFunction += Math.log(P.getPixelValue(n, m)) - S.getPixelValue(n, m) / P.getPixelValue(n, m);
									objectiveFunction += Math.log(P.getPixelValue(255, n)) - S.getPixelValue(255, n) / P.getPixelValue(255, n);
									
									if (objectiveFunction > ojMax) {
										ojMax = objectiveFunction;
										thresholds[0] = i;
										thresholds[1] = j;
										thresholds[2] = k;
										thresholds[3] = m;
										thresholds[4] = n;
									}
								}
							}
						}						
					}
				}				
				break;
			
		}

		//double[] x = new double[256]; for (int i = 0 ; i < x.length ; i++) x[i] = i;
		//disp.plotHistogram(x, rC.castInt2Double(eightHist), "test", "histogram classes", "frequency");
		
		return thresholds;
		
	}
	
	public int[] multiMinError(InputOutput.MyFileCollection mFC, int[] nowHist, MenuWaiter.MultiRegionSegmentation mRS) {
		
		//following S. Schlueter's code from QuantIM which in turn is based on
		//Kittler & Illingworth (1986): Pattern Recognition, 19(1), 41-47.
		
		DisplayThings disp = new DisplayThings();
		RollerCaster rC = new RollerCaster();
		HistogramStuff hist = new HistogramStuff();
		
		int[] thresholds = new int[mRS.numberOfThresholds];
		
		//if histogram is in 16 bit, scale to 8 bit
		int[] eightHist = new int[256];
		if (mRS.convert2EightBit) hist.convert16to8BitHistogram(nowHist, mRS.lowerBoundary, mRS.upperBoundary);
		else eightHist = nowHist;
		
		ImageProcessor P = calcLookUpTableP(eightHist);  //lut0 in QuantIM;
		ImageProcessor logP = calcLookUpTableLogP(eightHist);
		ImageProcessor S = calcLookUpTableS(P);   //lut1 in QuantIM;
		ImageProcessor[] S0 = calcLookUpTableS4MinErr(P, S);
		ImageProcessor logS = S0[1];   // lut2 in QuantIM;
		
/*		ImagePlus test0 = new ImagePlus("logP", logP);		
		test0.draw();test0.show();
		
		ImagePlus test1 = new ImagePlus("logS", logS);		
		test1.draw();test1.show();
*/
		
		//find maximum...
		double ojMax = 1e20;
		switch(mRS.numberOfThresholds) {
			
			case 1: 
				
				for (int i = 0 ; i < 256 - mRS.numberOfThresholds ; i++) {

					double objectiveFunction = P.getPixelValue(i, 0) * (logS.getPixelValue(i, 0) - logP.getPixelValue(i, 0));					
					objectiveFunction += P.getPixelValue(255, i) * (logS.getPixelValue(255, i) - logP.getPixelValue(255, i));
			
						
					if (objectiveFunction < ojMax) {
						ojMax = objectiveFunction;
						thresholds[0] = i;
					}
						
				}
				break;
			
			case 2: 
				
				for (int i = 1 ; i < 256 - mRS.numberOfThresholds ; i++) {
					for (int j = i + 1 ; j < 256 - 1 ; j++) {										
												
							
						double objectiveFunction = P.getPixelValue(i, 0) * (logS.getPixelValue(i, 0) - logP.getPixelValue(i, 0));
						objectiveFunction += P.getPixelValue(j, i) * (logS.getPixelValue(j, i) - logP.getPixelValue(j, i));
						objectiveFunction += P.getPixelValue(255, j) * (logS.getPixelValue(255, j) - logP.getPixelValue(255, j));
						
						if (objectiveFunction < ojMax) {
							ojMax = objectiveFunction;
							thresholds[0] = i;
							thresholds[1] = j;
						}
							
					}
				}
				break;
				
			case 3: 
			
				for (int i = 0 ; i < 256 - mRS.numberOfThresholds ; i++) {
					for (int j = i + 1 ; j < 256 - 2 ; j++) {
						for (int k = j + 1; k < 256 - 1 ; k++) {						
						
							double objectiveFunction = P.getPixelValue(i, 0) * (logS.getPixelValue(i, 0) - logP.getPixelValue(i, 0));
							objectiveFunction += P.getPixelValue(j, i) * (logS.getPixelValue(j, i) - logP.getPixelValue(j, i));
							objectiveFunction += P.getPixelValue(k, j) * (logS.getPixelValue(k, j) - logP.getPixelValue(k, j));
							objectiveFunction += P.getPixelValue(255, k) * (logS.getPixelValue(255, k) - logP.getPixelValue(255, k));	
							
							if (objectiveFunction < ojMax) {
								ojMax = objectiveFunction;
								thresholds[0] = i;
								thresholds[1] = j;
								thresholds[2] = k;
							}
						}						
					}
				}
				break;
				
			case 4: 
				
				for (int i = 0 ; i < 256 - mRS.numberOfThresholds ; i++) {
					for (int j = i + 1; j < 256 - 3 ; j++) {
						for (int k = j + 1; k < 256 - 2 ; k++) {
							for (int m = k + 1; m < 256 - 1 ; m++) {

								double objectiveFunction = P.getPixelValue(i, 0) * (logS.getPixelValue(i, 0) - logP.getPixelValue(i, 0));
								objectiveFunction += P.getPixelValue(j, i) * (logS.getPixelValue(j, i) - logP.getPixelValue(j, i));
								objectiveFunction += P.getPixelValue(k, j) * (logS.getPixelValue(k, j) - logP.getPixelValue(k, j));
								objectiveFunction += P.getPixelValue(m, k) * (logS.getPixelValue(m, k) - logP.getPixelValue(m, k));	
								objectiveFunction += P.getPixelValue(255, m) * (logS.getPixelValue(255, m) - logP.getPixelValue(255, m));
								
								if (objectiveFunction < ojMax) {
									ojMax = objectiveFunction;
									thresholds[0] = i;
									thresholds[1] = j;
									thresholds[2] = k;
									thresholds[3] = m;
								}
							}
						}						
					}
				}
				break;
				
			case 5: 
			
				for (int i = 0 ; i < 256 - mRS.numberOfThresholds ; i++) {
					for (int j = i + 1; j < 256 - 4 ; j++) {
						for (int k = j + 1; k < 256 - 3 ; k++) {
							for (int m = k + 1; m < 256 - 2 ; m++) {
								for (int n = m + 1 ; n < 256 -1 ; n++) {
							
									double objectiveFunction = P.getPixelValue(i, 0) * (logS.getPixelValue(i, 0) - logP.getPixelValue(i, 0));
									objectiveFunction += P.getPixelValue(j, i) * (logS.getPixelValue(j, i) - logP.getPixelValue(j, i));
									objectiveFunction += P.getPixelValue(k, j) * (logS.getPixelValue(k, j) - logP.getPixelValue(k, j));
									objectiveFunction += P.getPixelValue(m, k) * (logS.getPixelValue(m, k) - logP.getPixelValue(m, k));	
									objectiveFunction += P.getPixelValue(n, m) * (logS.getPixelValue(n, m) - logP.getPixelValue(n, m));
									objectiveFunction += P.getPixelValue(255, n) * (logS.getPixelValue(255, n) - logP.getPixelValue(255, n));
									
									if (objectiveFunction < ojMax) {
										ojMax = objectiveFunction;
										thresholds[0] = i;
										thresholds[1] = j;
										thresholds[2] = k;
										thresholds[3] = m;
										thresholds[4] = n;
									}
								}
							}
						}						
					}
				}				
				break;
			
		}

		//double[] x = new double[256]; for (int i = 0 ; i < x.length ; i++) x[i] = i;
		//disp.plotHistogram(x, rC.castInt2Double(eightHist), "test", "histogram classes", "frequency");
		
		return thresholds;
		
	}
	
	public ImageProcessor calcLookUpTableP(int[] nowHist) {
		
		ImageProcessor P = new FloatProcessor(256, 256);
		
		double tiny = 1e-20;
		double sumHist = tiny;
		P.putPixelValue(0, 0, sumHist);
		for (int x = 1 ; x < 256 ; x++) {
			sumHist += nowHist[x];
			if (nowHist[x] == 0) sumHist += tiny;
			P.putPixelValue(x, 0, sumHist);
		}
			
		for (int y = 1  ; y < 256 ; y++) {
			for (int x = y + 1 ; x < 256 ; x++) {
				
				P.putPixelValue(x, y, P.getPixelValue(x, 0) - P.getPixelValue(y, 0));
				
			}	
		}		
		
		//get maximal entry and normalize between 0 and 1
/*		double maxEntry = P.getMax();
		
		for (int y = 0  ; y < 256 ; y++) {
			for (int x = y ; x < 256 ; x++) {
				
				double pNow = P.getPixelValue(x, y) / maxEntry; 
				P.putPixelValue(x, y, pNow);
				
			}	
		}		*/
		
			
		return P;
		
	}
	
	public ImageProcessor calcLookUpTableLogP(int[] nowHist) {
		
		ImageProcessor logP = new FloatProcessor(256, 256);
		
		double tiny = 1e-20;
		double sumHist = 0;
		logP.putPixelValue(0, 0, sumHist);
		for (int x = 1 ; x < 256 ; x++) {
			if (nowHist[x] > 0) sumHist += nowHist[x];
			else sumHist += tiny;
			logP.putPixelValue(x, 0, Math.log(sumHist));
		}
			
		for (int y = 1  ; y < 256 ; y++) {
			for (int x = y + 1 ; x < 256 ; x++) {
				
				double P = Math.exp(logP.getPixelValue(x, 0)) - Math.exp(logP.getPixelValue(y, 0));
				logP.putPixelValue(x, y, Math.log(P));
				
			}	
		}		
		
		return logP;
		
	}
	
	public ImageProcessor calcLookUpTableS(ImageProcessor P) {
		
		ImageProcessor S = new FloatProcessor(256, 256);
		
		double tiny = 1e-20;
		double sumS = tiny;
		S.putPixelValue(0, 0, sumS);
		for (int x = 1 ; x < 256 ; x++) {
			sumS += x * (P.getPixelValue(x, 0) - P.getPixelValue(x - 1, 0));
			S.putPixelValue(x, 0, sumS);
		}
			
		for (int y = 1  ; y < 256 ; y++) {
			for (int x = y + 1 ; x < 256 ; x++) {
				
				S.putPixelValue(x, y, S.getPixelValue(x, 0) - S.getPixelValue(y, 0));
				
			}	
		}		
		
		return S;
		
	}
	
	public ImageProcessor calcLookUpTableS4MaxEntro(ImageProcessor P) {
		
		ImageProcessor S = new FloatProcessor(256, 256);
		double tiny = 1e-20;
		
		double sumS = 0;
		S.putPixelValue(0, 0, sumS);
		for (int x = 1 ; x < 256 ; x++) {
			double numAtX = P.getPixelValue(x, 0) - P.getPixelValue(x - 1, 0);
			if (numAtX > 0) sumS += numAtX * Math.log(numAtX);
			else sumS += tiny;
			S.putPixelValue(x, 0, sumS);
		}
			
		for (int y = 1  ; y < 256 ; y++) {
			for (int x = y + 1 ; x < 256 ; x++) {
				
				S.putPixelValue(x, y, S.getPixelValue(x, 0) - S.getPixelValue(y, 0));
				
			}	
		}		
		
		return S;
		
	}
	
	public ImageProcessor[] calcLookUpTableS4MinErr(ImageProcessor P, ImageProcessor S) {
		
		ImageProcessor[] S0 = new FloatProcessor[2];
		ImageProcessor S1 = new FloatProcessor(256, 256);
		ImageProcessor S2 = new FloatProcessor(256, 256);
		double tiny = 1d/256d;
		
		double sumS = tiny;			
		S1.putPixelValue(0, 0, Math.log(Math.sqrt(sumS)));	
		double[] checkS = new double[256];
		checkS[0] = Math.log(Math.sqrt(tiny));
		for (int x = 1 ; x < 256 ; x++) {	
			S1.putPixelValue(x, 0, S.getPixelValue(x, 0) / P.getPixelValue(x, 0));
			
			sumS = 0; 			//cm2 in QuantIM	
			double Px = P.getPixelValue(x, 0);
			for (int j = 1 ; j <= x ; j++) {
				double deltaPj = (P.getPixelValue(j, 0) - P.getPixelValue(j - 1, 0));				
				if (P.getPixelValue(j - 1, 0) < 1e-10) deltaPj = 0;
				double nowse = S1.getPixelValue(x, 0);
				sumS += (j - nowse) * (j - nowse) * deltaPj;
			}
			sumS /= Px;  //second central moment
			if (sumS == 0) sumS = tiny;
			
			checkS[x] = Math.log(Math.sqrt(sumS));
			
			S2.putPixelValue(x, 0, Math.log(Math.sqrt(sumS)));
			
		}
		
		for (int y = 1  ; y < 256 ; y++) {
			for (int x = y + 1 ; x < 256 ; x++) {
								
				S1.putPixelValue(x, y, (S.getPixelValue(x, 0) - S.getPixelValue(y, 0)) / P.getPixelValue(x, y));
				
				sumS = 0; 
				for (int k = y ; k < x ; k++) {
					double deltaPk = (P.getPixelValue(k, 0) - P.getPixelValue(k - 1, 0));				
					if (P.getPixelValue(k - 1, 0) < 1e-10) deltaPk = 0;
					double nowS1 = S1.getPixelValue(x, y);
					sumS += (k - nowS1) * (k - nowS1) * deltaPk;
				}
				sumS /= P.getPixelValue(x, y);
				
				if (sumS == 0) S2.putPixelValue(x, y, Math.log(Math.sqrt(tiny)));
				else S2.putPixelValue(x, y, Math.log(Math.sqrt(sumS)));
				
			}
		}
			
/*		ImagePlus test0 = new ImagePlus("S1", S1);		
		test0.draw();test0.show();
		
		ImagePlus test1 = new ImagePlus("S2", S2);		
		test1.draw();test1.show();*/
		
		S0[0] = S1;
		S0[1] = S2;
		
		return S0;
		
	}
	
		
	public ImagePlus[] applyAChosenThresholdingMethod(ImagePlus nowTiff, InputOutput.MyFileCollection mFC, MenuWaiter.ThresholderMenuReturn mTMR, int[] myZ) {
		
		String pathSep = "/";
		
		InputOutput jIO = new InputOutput();	
		HistogramStuff hist = new HistogramStuff();
		RoiHandler rH = new RoiHandler();
		ObjectDetector jOD = new ObjectDetector();
		
		ImagePlus[] outImg = {null, null, null};
		ImageStack myStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
					
		int[] myHist = new int[256 * 256];
		int[] newHist;
		int myThresh = 0;

		AutoThresholder myAuto = new AutoThresholder();		
	
		PolygonRoi[] pRoi = null;		
		PolygonRoi[] oRoi = null;
		int[] onePercQuantile = new int[nowTiff.getNSlices()];
		
		//read gauge file AND find illumination of this column	
		if (mTMR.useInnerCircle) {
			
			if (mFC.nowInnerCirclePath.contains("Gauge")) {
				
				ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();				
				int versio = jIO.checkInnerCircleFileVersion(mFC.nowInnerCirclePath);			
				if (versio == 0) jCO = jIO.readInnerCircleVer0(mFC.nowInnerCirclePath);	
				else jCO = jIO.readInnerCircleVer1(mFC.nowInnerCirclePath);	
				
				//create ROI stacks..
				pRoi = rH.makeMeAPolygonRoiStack("inner", "exact", jCO, 0);
				oRoi = rH.makeMeAPolygonRoiStack("inner", "manual", jCO, 10);	
				
			}			
			else {

				//read column coordinates
				ObjectDetector.EggShapedColCoords3D jCO = jOD.new EggShapedColCoords3D();						
				jCO = jIO.readInnerCircleSteel(mFC);
		
				//create ROI stacks..
				pRoi = rH.makeMeAPolygonRoiStack("inner", "exact", jCO, 0);
				oRoi = rH.makeMeAPolygonRoiStack("inner", "manual", jCO, 10);	
			
			}
		}
		
		//find threshold	
		if (!mTMR.useConstantThreshold) {
	
			//get the stacked histogram		
			for (int i = 1 ; i < nowTiff.getNSlices() + 1 ; i++) {
				
				IJ.showStatus("Getting histogram of slice #" + i + "/" + nowTiff.getNSlices());
				
				nowTiff.setPosition(i);		
				
				ImageProcessor myIP = nowTiff.getProcessor();
				ImageProcessor modIP = myIP.duplicate();
				
				//cut out everything outside column
				if (mTMR.useInnerCircle) {
					modIP.setRoi(pRoi[i - 1]);
					modIP.setColor(0);
					modIP.fillOutside(pRoi[i - 1]);
				}
				
				newHist=modIP.getHistogram();
				newHist[0] = 0; //set zero GV to zero
				int gotoc = 256;
				if (nowTiff.getBitDepth() == 16) gotoc = 256 * 256;
				for (int j = 0 ; j < gotoc ; j++) {
					myHist[j] = myHist[j] + newHist[j];
				}			
				
				//save 1% perc and wallOfgray
				double[] cumHist = hist.calcCumulativeHistogram(newHist);			
				onePercQuantile[i-1] = hist.findPercentileFromCumHist(cumHist, 0.001);	
			}		
			
			//prepare possible scaling nowTiff to 8-bit
			double[] cumHist = hist.calcCumulativeHistogram(myHist);
			double lBound = 0.0000001;double uBound = 0.99999;			
			int lowestVal = (int)Math.round(hist.findPercentileFromCumHist(cumHist, lBound));			//set lower gray value
			int largestVal = (int)Math.round(hist.findPercentileFromCumHist(cumHist, uBound));		
			float valSpan = (float)largestVal - lowestVal;
			if (!mTMR.setMaxgray2Wallgray | mTMR.useConstantThreshold == true) valSpan = largestVal - lowestVal;
		
			ImagePlus eightBitTiff = new ImagePlus();
			
			if (nowTiff.getBitDepth() == 8) {
				eightBitTiff = nowTiff.duplicate();
			}
			else {
			
				ImageStack eightBitStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
				
				for (int i = 1 ; i < nowTiff.getNSlices() + 1 ; i++) {
					
					IJ.showStatus("Scaling 16-bit image to 8-bit at slice #" + i + "/" + nowTiff.getNSlices());
					
					nowTiff.setPosition(i);
					ImageProcessor myIP = nowTiff.getProcessor();
					ImageProcessor modIP = myIP.duplicate();
					modIP.max(largestVal);
					modIP.subtract(lowestVal);
					modIP.multiply(256 / valSpan);
					modIP.min(0);
					modIP.max(255);
					ImageProcessor eightIP = modIP.convertToByte(false);
					
					eightBitStack.addSlice(eightIP);
				}
				
				eightBitTiff.setStack(eightBitStack);
			}
	
			//get the 8-bit stacked histogram
			int[] my8Hist = new int[256];
			int[] new8Hist;
			for (int i = 1 ; i < eightBitTiff.getNSlices() + 1 ; i++) {
				
				IJ.showStatus("Getting 8-bit histogram of slice #" + i + "/" + (eightBitTiff.getNSlices()));
				
				eightBitTiff.setPosition(i);			
				
				ImageProcessor myIP = eightBitTiff.getProcessor();
				ImageProcessor modIP = myIP.duplicate();
				
				//cut out everything outside column
				if (mTMR.useInnerCircle) {
					modIP.setRoi(pRoi[i - 1]);				
					modIP.setColor(0);
					modIP.fillOutside(pRoi[i - 1]);	
				}
				
				new8Hist=modIP.getHistogram();	
				for (int j = 0 ; j < my8Hist.length ; j++) {
					my8Hist[j] = my8Hist[j] + new8Hist[j];
				}			
				my8Hist[0] = 0; //set zero GV to zero
				my8Hist[255] = 0; //set last GV to zero	
			}				
		
			//find primary threshold
			myThresh = myAuto.getThreshold(mTMR.myPrimaryMethod, my8Hist);
		
			//find secondary threshold
			if(mTMR.mySecondaryMethod != null) {
				for(int i = myThresh ; i < my8Hist.length ; i++) my8Hist[i] = 0;
				myThresh = myAuto.getThreshold(mTMR.mySecondaryMethod, my8Hist);
			}			
				
			//do binarization
			for (int i = 1 ; i < eightBitTiff.getNSlices() + 1 ; i++) {
				
				IJ.showStatus("Binarizing slice #" + i + "/" + (eightBitTiff.getNSlices()));			
							
				eightBitTiff.setPosition(i);
				ImageProcessor myIP = eightBitTiff.getProcessor();
				
				//apply the threshold
				myIP.threshold(myThresh);
				
				//create a binary where the less dense phase is 255 and everything else 0
				myIP.invert();
				//for (int x = 0 ; x < myIP.getWidth() ; x++) {
				//	for (int y = 0 ; y < myIP.getHeight() ; y++) {						
				//		int thisImgPixel = myIP.getPixel(x, y);					
				//		myIP.putPixel(x, y, thisImgPixel);
				//	}
				//}
				
				//cut out everything outside column
				if (mTMR.useInnerCircle) {
					myIP.setRoi(pRoi[i - 1]);				
					myIP.setColor(0);
					myIP.fillOutside(pRoi[i - 1]);
				}
			
				myStack.addSlice(myIP);
			}
			
			//create binary out image 
			ImagePlus sillyTiff = new ImagePlus();
			sillyTiff.setStack(myStack);
			outImg[0] = sillyTiff;	
						
			//outImg[0].updateAndDraw();
			//outImg[0].show();
			
			//nowTiff.updateAndDraw();
			//nowTiff.show();
			
			//also save threshold		
			String thresholdSaverPath = mFC.myOutFolder + pathSep + "EightBitThresholds.txt";
			jIO.writeThreshold(thresholdSaverPath, nowTiff.getShortTitle(), myThresh);	
			
			//also save threshold comparison images	
			if ((mTMR.save3DImage | mTMR.save4GeoDict) & mTMR.save4Evaluation) jIO.writeSnapshots4ThresholdComparison(mFC, nowTiff, outImg[0], oRoi, myZ);
			if (mTMR.save4Evaluation & !(mTMR.save3DImage | mTMR.save4GeoDict)) {
				for (int i = 0 ; i < myZ.length ; i++) myZ[i] = i + 1;
				jIO.writeSnapshots4ThresholdComparison(mFC,  nowTiff, outImg[0], oRoi, myZ);
			}
		}
		
		else { 	//in case that a constant threshold is used for all images..
	
			if (mTMR.minThreshold > 0 | mTMR.maxThreshold > 0) {
				
				ImageStack nowStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());		
				
				for (int i = 0 ; i < nowTiff.getNSlices() ; i++) {
				
					IJ.showStatus("Binarizing slice " + i + "/" + nowTiff.getNSlices());
				
					nowTiff.setPosition(i+1);
					ImageProcessor myIP = nowTiff.getProcessor();
					ImageProcessor lowIP = myIP.duplicate();
					ImageProcessor highIP = myIP.duplicate();
		
					if (mTMR.useInnerCircle) {
						lowIP.setRoi(oRoi[i]);					
						lowIP.setColor(0);
						lowIP.fillOutside(oRoi[i]);
						
						highIP.setRoi(oRoi[i]);					
						highIP.setColor(0);
						highIP.fillOutside(oRoi[i]);
					}
					
		/*			ImagePlus test = new ImagePlus("",lowIP);
					test.updateAndDraw();
					test.show();*/
					
					lowIP.threshold(mTMR.minThreshold);
					highIP.threshold(mTMR.maxThreshold);
										
					StackCalculator sC = new StackCalculator();	
					ImagePlus l = new ImagePlus("", lowIP);
					ImagePlus h = new ImagePlus("", highIP);
					ImagePlus outI = sC.subtract(l, h);
					ImageProcessor outIP = outI.getProcessor();					
					
					ImageProcessor eightIP = outIP.convertToByte(false);
				
					nowStack.addSlice(eightIP);				
				}
				
				//create binary out image 
				ImagePlus binTiff = new ImagePlus();
				binTiff.setStack(nowStack);
				outImg[0] = binTiff;
				
			/*	binTiff.updateAndDraw();
				binTiff.show();*/
				
			}
			
			//also save threshold comparison images			
			if (mTMR.save4Evaluation) jIO.writeSnapshots4ThresholdComparison(mFC, nowTiff, outImg[0], oRoi, myZ);			

		}
					
		return outImg;
	}
	
	public ImagePlus binarize3DImage(ImagePlus nowTiff, int threshold) {
		
		ImagePlus outTiff = new ImagePlus();
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			
			nowTiff.setSlice(z + 1);			
			ImageProcessor nowIP = nowTiff.getProcessor().duplicate();			
			
			ImageProcessor outIP = new ByteProcessor(nowIP.getWidth(), nowIP.getHeight());
			if (nowTiff.getBytesPerPixel() == 2) {
				
				for (int x = 0 ; x < outIP.getWidth() ; x++) {
					for (int y = 0 ; y < outIP.getHeight() ; y++) {
						
						int nowPix = nowIP.getPixel(x, y);
						
						if (nowPix > threshold) outIP.putPixel(x, y, 255);
						else outIP.putPixel(x, y, 0);					
						
					}					
				}
				
			}
			else {
				nowIP.threshold(threshold);
				outIP = nowIP;
			}			
					
			outStack.addSlice(outIP);
			
		}
	
		outTiff.setStack(outStack);	
			
		return outTiff;
		
	}

	
	/** 
	 * Puts each pixel in each slice to 0 if below threshold, otherwise to 255
	 * @param nowTiff
	 * @param threshold
	 * @return ImagePlus outTiff
	 */
	public ImagePlus binarize3DFloatImage(ImagePlus nowTiff, double threshold) {
		
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImagePlus outTiff = new ImagePlus();
		
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			
			nowTiff.setSlice(z + 1);
			
			ImageProcessor nowIP = nowTiff.getProcessor();
			ImageProcessor newIP = new ByteProcessor(nowTiff.getWidth(), nowTiff.getHeight());
			
			for (int x = 0 ; x < nowTiff.getWidth() ; x++) {
				
				for (int y = 0 ; y < nowTiff.getHeight() ; y++) {
				
					float nowPixel = nowIP.getPixelValue(x, y);
					if (nowPixel < threshold) newIP.putPixelValue(x, y, 0);
					else newIP.putPixelValue(x, y, 255);
					
				}
				
			}
			
			outStack.addSlice(newIP);
			
		}		
		
		outTiff.setStack(outStack);
		
		return outTiff;
		
	}
	
	public void segmentTernaryAndSave(InputOutput.MyFileCollection mFC, ImagePlus nowTiff, ImagePlus gradTiff, ImagePlus pomRegion) {
	
 		InputOutput jIO = new InputOutput();
 		ImageManipulator jIM = new ImageManipulator();
		
		int w = nowTiff.getWidth();
		int h = nowTiff.getHeight();
			
		ImagePlus markerTiff = new ImagePlus();		
		ImagePlus outTiff = new ImagePlus();
		ImageStack markerStack = new ImageStack(w, h);
			
		//scale again to original height
		pomRegion = jIM.scaleWithoutMenu(pomRegion, pomRegion.getWidth(), pomRegion.getHeight() * 2, 1, ImageProcessor.BILINEAR);
		
		//assign POM mask process
		ImageProcessor pomRegIP = pomRegion.getProcessor();
		
		//flip back upright again		
		pomRegIP.flipVertical();
			
		//determine lower and upper threshold
		int myMin = 70000; 
		int myMax = 0;
		for (int y = 0 ; y < 6554 ; y++) {
			for (int x = 0 ; x < 6554 ; x++) {
					
				int nowPix = pomRegIP.getPixel(x, y);
			
				if (nowPix > 0 & x * 10 < myMin) myMin = x * 10;
				if (nowPix > 0 & x * 10 > myMax) myMax = x * 10;
				
			}
			
		}
		
		//get reduced starting range
		int deltaPOM = myMax - myMin;
				
		//do the segmentation
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			
			IJ.showStatus("Segmenting slice " + (z+1) + "/" + nowTiff.getNSlices());
			
			nowTiff.setSlice(z + 1);
			ImageProcessor mP = nowTiff.getProcessor();
			
			gradTiff.setSlice(z + 1);
			ImageProcessor gP = gradTiff.getProcessor().convertToShort(false);
			
			ImageProcessor markerIP = new ByteProcessor(w, h);
			ImageProcessor pomSeedsIP = new ByteProcessor(w, h);
			
			for (int x = 0 ; x < mP.getWidth(); x++) {
				
				for (int y = 0 ; y < mP.getHeight(); y++) {
					
					int mPix = mP.getPixel(x, y);
					int gPix = gP.getPixel(x, y);	
					
					//do segmentation
					markerIP.putPixel(x, y, 0);
					if (mPix < myMin + 0.33 * deltaPOM) markerIP.putPixel(x, y, 1);  
					if (mPix > myMax) markerIP.putPixel(x, y, 3);
					
					//do extra IP for pomRegion Seeds
					if (pomRegIP.getPixel(mPix/10, gPix) > 0) pomSeedsIP.putPixel(x, y, 255);
					
				}
			}			
			
			//ImagePlus test = new ImagePlus("",pomSeedsIP);
			//test.show();
			
			//filter out thin POM seed regions
			pomSeedsIP.dilate();
			
			for (int x = 0 ; x < mP.getWidth(); x++) {
				
				for (int y = 0 ; y < mP.getHeight(); y++) {
					
					float mPix = pomSeedsIP.getPixel(x, y);
					
					if (mPix == 255) markerIP.putPixel(x, y, 2);	
					
				}
			}
			
			//add new markerSlice
			markerStack.addSlice(markerIP);
		
		}	
		
		markerTiff.setStack(markerStack);
		
		IJ.freeMemory();IJ.freeMemory();
		
		//markerTiff.updateAndDraw();
		//markerTiff.show();
		
		//use watershed algorithm
		int connectivity = 26;		
		MarkerControlledWatershed3DPlugin mcWS = new MarkerControlledWatershed3DPlugin();
		MarkerControlledWatershed3DPlugin.binaryMarkers = false;
		MarkerControlledWatershed3DPlugin.getDams = false;		
		outTiff = mcWS.process(gradTiff, markerTiff, null, connectivity);		
		
		//save images		
		jIO.tiffSaver(mFC.myOutFolder, mFC.fileName, outTiff);
		
	}
	
	public PhaseOfInterestInfo parsePhaseOfInterestInfo(String myPhase) {
		
		PhaseOfInterestInfo mPII = new PhaseOfInterestInfo();
		
		//parse myPhase
		try {
			
			if (myPhase.contains(">") | myPhase.contains("<")) {
				
				mPII.applyThreshold = true;
					
				for (int i = 0 ; i < myPhase.length() - 1 ; i++) {
					if (myPhase.substring(i, i + 1).equalsIgnoreCase(">")) {
						mPII.threshold = Integer.parseInt(myPhase.substring(i + 1));
						mPII.smallerLarger = ">";
						mPII.outName = "GT" + mPII.threshold;
						break;
					}
					if (myPhase.substring(i, i + 1).equalsIgnoreCase("<")) {
						mPII.threshold = Integer.parseInt(myPhase.substring(i + 1));
						mPII.smallerLarger = "<";
						mPII.outName = "LT" + mPII.threshold;
						break;
					}
				}
				
			}
			
			else mPII.myPhaseValue = Integer.parseInt(myPhase);
			
		}
		
		catch(Exception e) {    			
		
			IJ.error("Invalid input string!\nTry entering an Integer number between 0 and 255 or something like >0 or <255 or similar.\nNote that >= or <= are not allowed!");
			return null;
			
	    }
	
		return mPII;
		
		
	}
	
	public ImagePlus extractPhaseOfInterest(ImagePlus nowTiff, String myPhase, String myPhaseName) {
		
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		if (myPhase.compareToIgnoreCase("null") != 0) {
		
			PhaseOfInterestInfo mPII = parsePhaseOfInterestInfo(myPhase);		
		
			for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
				
				IJ.showStatus("Creating binary for '" + myPhase + "_" + myPhaseName + "' slice " + (z + 1) + "/" + nowTiff.getNSlices());
				
				nowTiff.setSlice(z + 1);
				ImageProcessor nowIP = nowTiff.getProcessor();			
						
				for (int x = 0 ; x < nowTiff.getWidth() ; x++) {
					for (int y = 0 ; y < nowTiff.getHeight() ; y++) {
						
						int nowPix = nowIP.getPixel(x, y);
						
						if (!mPII.applyThreshold) {
							
							if (nowPix == mPII.myPhaseValue) nowIP.putPixel(x, y, 255);
							else nowIP.putPixel(x, y, 0);
							
						}
						else {
						 
							if (mPII.smallerLarger.equalsIgnoreCase("<")) {
							
								if (nowPix < threshold) nowIP.putPixel(x, y, 255);
								else nowIP.putPixel(x, y, 0);
								
							}
							
							if (mPII.smallerLarger.equalsIgnoreCase(">")) {
								
								if (nowPix > mPII.threshold) nowIP.putPixel(x, y, 255);
								else nowIP.putPixel(x, y, 0);
								
							}
						}
					}
				}
				
				outStack.addSlice(nowIP);				
			}
			
			nowTiff.setStack(outStack);
			
		}
		
		
		
		//nowTiff.updateAndDraw();nowTiff.show();
		
		return nowTiff;
		
	}
	
	public ImagePlus[] removePoresAboveSurface(ImagePlus nowTiff, int maxTopDepression, ImageProcessor topIP, PolygonRoi[] pRoi, PolygonRoi[] iRoi) {
		

		ImagePlus[] outTiff = new ImagePlus[2];
		ImagePlus img1 = new ImagePlus();
		ImagePlus img2 = new ImagePlus();
		
		ImageStack corrected4SurfaceStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight()); 
		ImageStack corr4ThicknessAnalyses = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight()); 
		
		//cut away soil above the soil surface
		for (int i = 0 ; i < maxTopDepression ; i++) {
			
			nowTiff.setPosition(i + 1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			if (i < maxTopDepression) {		
				
				IJ.showStatus("Correcting for soil top surface in slice #" + (i + 1) + "/" + maxTopDepression);				
				
				ImageProcessor modIP = nowIP.duplicate();
				ImageProcessor forThIP = nowIP.duplicate();
				
				for (int x = 0 ; x < nowIP.getWidth() ; x++) {
					for (int y = 0 ; y < nowIP.getHeight() ; y++) {					
						int nowPix = nowIP.getPixel(x, y);
						int surPix = topIP.getPixel(x, y);
						if (i <= surPix) modIP.putPixelValue(x, y, 0);
						else modIP.putPixelValue(x, y, nowPix);
						if (i <= surPix - 100) forThIP.putPixelValue(x, y, 0);		//to prevent underestimation of pore diameters reaching the top..
						else forThIP.putPixelValue(x, y, nowPix);	
					}			
				}			
				
				modIP.setColor(Color.BLACK);
				modIP.setRoi(pRoi[i]);
				modIP.fillOutside(pRoi[i]);
				modIP.resetRoi();
				
				forThIP.setColor(Color.BLACK);
				forThIP.setRoi(pRoi[i]);
				forThIP.fillOutside(pRoi[i]);
				forThIP.resetRoi();
					
				if (iRoi != null) {
					modIP.setColor(Color.BLACK);
					modIP.setRoi(iRoi[i]);
					modIP.fill(iRoi[i]);
					modIP.resetRoi();
					
					forThIP.setColor(Color.BLACK);
					forThIP.setRoi(iRoi[i]);
					forThIP.fill(iRoi[i]);
					forThIP.resetRoi();
				}
										
				corrected4SurfaceStack.addSlice(modIP);
				corr4ThicknessAnalyses.addSlice(forThIP);
				
			}
		}
		
		img1.setStack(corrected4SurfaceStack);
		img2.setStack(corr4ThicknessAnalyses);
		
		outTiff[0] = img1;
		outTiff[1] = img2;
		
		return outTiff;			
	}
	
	public ImagePlus[] removePoresBelowSurface(ImagePlus nowTiff, int maxBotDepression, ImageProcessor botIP, PolygonRoi[] pRoi, PolygonRoi[] iRoi) {
		
		ImagePlus[] outTiff = new ImagePlus[2];
		ImagePlus img1 = new ImagePlus();
		ImagePlus img2 = new ImagePlus();
		
		ImageStack corrected4SurfaceStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImageStack corr4ThicknessAnalyses = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		//cut away soil above the soil surface
		for (int i = maxBotDepression ; i < nowTiff.getNSlices() ; i++) {
			
			IJ.showStatus("Correcting for soil bottom surface in slice #" + (maxBotDepression - (nowTiff.getNSlices() - i)) + "/" + maxBotDepression);				
		
			nowTiff.setPosition(i + 1);
			ImageProcessor nowIP = nowTiff.getProcessor();
		
			ImageProcessor modIP = nowIP.duplicate();
			ImageProcessor forThIP = nowIP.duplicate();
			
			for (int x = 0 ; x < nowIP.getWidth() ; x++) {
				for (int y = 0 ; y < nowIP.getHeight() ; y++) {					
					int nowPix = nowIP.getPixel(x, y);
					int surPix = nowTiff.getNSlices() - botIP.getPixel(x, y);
					if (i >= surPix) modIP.putPixelValue(x, y, 0);
					else modIP.putPixelValue(x, y, nowPix);
					if (i >= surPix + 100) forThIP.putPixelValue(x, y, 0);		//to prevent underestimation of pore diameters reaching the top..
					else forThIP.putPixelValue(x, y, nowPix);	
				}			
			}		
			
			modIP.setColor(Color.BLACK);
			modIP.setRoi(pRoi[i]);
			modIP.fillOutside(pRoi[i]);
			modIP.resetRoi();
			
			forThIP.setColor(Color.BLACK);
			forThIP.setRoi(pRoi[i]);
			forThIP.fillOutside(pRoi[i]);
			forThIP.resetRoi();
			
			corrected4SurfaceStack.addSlice(modIP);
			corr4ThicknessAnalyses.addSlice(forThIP);
		}
		

		img1.setStack(corrected4SurfaceStack); // is empty?
		img2.setStack(corr4ThicknessAnalyses);
		
		outTiff[0] = img1;
		outTiff[1] = img2;
		
		return outTiff;	
		
	}
	
	public ImagePlus[] addColumnRegion2Stack(ImagePlus nowTiff, int maxTopDepression, int maxBotDepression, PolygonRoi[] pRoi, PolygonRoi[] iRoi, boolean hasSurface) {
		
		ImagePlus[] outTiff = new ImagePlus[2];
		ImagePlus img1 = new ImagePlus();
		ImagePlus img2 = new ImagePlus();
		
		ImageStack corrected4SurfaceStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImageStack corr4ThicknessAnalyses = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		//decide whether using max bottom depression makes sense..
		if (maxBotDepression > nowTiff.getNSlices()) maxBotDepression = nowTiff.getNSlices();
		
		if (hasSurface) {
			//cut away soil above the soil surface
		
			for (int i = maxTopDepression ; i < maxBotDepression ; i++) {
			
				nowTiff.setPosition(i + 1);
				ImageProcessor nowIP = nowTiff.getProcessor();
				
				IJ.showStatus("Adding in-between slice " + (i - maxTopDepression + 1) + "/" + (nowTiff.getNSlices() - maxBotDepression -  maxTopDepression));
				
				ImageProcessor modIP = nowIP.duplicate();			
				ImageProcessor forThIP = nowIP.duplicate();
				
				modIP.setColor(Color.BLACK);
				modIP.setRoi(pRoi[i]);
				modIP.fillOutside(pRoi[i]);
				modIP.resetRoi();
				
				forThIP.setColor(Color.BLACK);
				forThIP.setRoi(pRoi[i]);
				forThIP.fillOutside(pRoi[i]);
				forThIP.resetRoi();
				
				/*ImagePlus tester = new ImagePlus("", modIP);
				tester.updateAndDraw();
				tester.show();*/
				
				corrected4SurfaceStack.addSlice(modIP);		
				corr4ThicknessAnalyses.addSlice(forThIP);
			}
			
			img1.setStack(corrected4SurfaceStack);
			img2.setStack(corr4ThicknessAnalyses);
			
			outTiff[0] = img1;
			outTiff[1] = img2;
		}		
		else {
			
			for (int i = maxTopDepression ; i < maxBotDepression ; i++) {
				
				nowTiff.setPosition(i + 1);
				ImageProcessor nowIP = nowTiff.getProcessor();
				
				IJ.showStatus("Adding in-between slice " + (i - maxTopDepression + 1) + "/" + (maxBotDepression -  maxTopDepression));
				
				ImageProcessor modIP = nowIP.duplicate();			
				
				modIP.setColor(Color.BLACK);
				modIP.setRoi(pRoi[i]);
				modIP.fillOutside(pRoi[i]);
				modIP.resetRoi();
			
				corrected4SurfaceStack.addSlice(modIP);		
			
			}
			
			img1.setStack(corrected4SurfaceStack);			
			
			outTiff[0] = img1;
			outTiff[1] = null;
			
		}
		
		return outTiff;	
		
	}
	
	public ImagePlus[] stackImagePlusArrays(ImagePlus[] top, ImagePlus[] bot) {
	
		ImagePlus[] outTiff = new ImagePlus[top.length];
		
		for (int i = 0 ; i < top.length ; i++) {
									
			ImagePlus topI = top[i];
			ImagePlus botI = bot[i];
			ImageStack outStack = new ImageStack(topI.getWidth(), topI.getHeight());
			
			for (int j = 0 ; j < topI.getNSlices() ; j++) {
				topI.setPosition(j+1);
				ImageProcessor nowIP = topI.getProcessor();
				outStack.addSlice(nowIP);				
			}
			
			for (int j = 0 ; j < botI.getNSlices() ; j++) {
				botI.setPosition(j+1);
				ImageProcessor nowIP = botI.getProcessor();
				outStack.addSlice(nowIP);				
			}
			
			ImagePlus newImg = new ImagePlus();
			newImg.setStack(outStack);
			outTiff[i] = newImg;		
			
			//newImg.updateAndDraw();
			//newImg.show();
			
		}
		
		return outTiff;
				
	}
			
	public ImagePlus binarizeInTwoSteps(ImagePlus myImg) {
		
		ImagePlus outImg = new ImagePlus();
						
		ImageStack myStack = new ImageStack(myImg.getWidth(), myImg.getHeight());
		
		ImageProcessor myIP = myImg.getProcessor().convertToByte(true);
				
		int[] myHist = new int[256];
		int[] newHist;
		int myThresh;
		int i, j;
		
		AutoThresholder myAuto = new AutoThresholder();
		String[] myMethods = AutoThresholder.getMethods();
		
		for (i = 1 ; i < myImg.getStackSize() + 1 ; i++) {
			myImg.setPosition(i);			
			myIP = myImg.getProcessor().convertToByte(true);
			newHist=myIP.getHistogram();	
			for (j = 0 ; j < myHist.length ; j++) {
				myHist[j] = myHist[j] + newHist[j];
			}			
		}
		
		myThresh = myAuto.getThreshold(myMethods[7], myHist);  //13 Renyi, 7 mean 
		
		for (i = myThresh; i < myHist.length ; i++) {
			myHist[i]=0;
		}
		
		int[] cutHist = new int[myThresh];
		for (i = 0 ; i < myThresh ; i++) {
			cutHist[i] = myHist[i];
		}
		
		myThresh = myAuto.getThreshold(myMethods[13], myHist);
				
		//do binarization
		for (i = 1 ; i < myImg.getStackSize() + 1 ; i++) {
			IJ.showStatus("Binarizing slice " + i + "/" + myImg.getNSlices());			
			myImg.setPosition(i);	
			myIP = myImg.getProcessor().convertToByte(true);
			myIP.threshold(myThresh);
			myIP.invert();
			myStack.addSlice(myIP);
		}
		
		//create binary outimage 
		outImg.setStack(myStack);
			
		return outImg;
	}
	
	public ImagePlus fuseSkeletonAndDistance(ImagePlus nowTiff, ImagePlus distTiff, int numberOfAddedLayers) {
		
		ImagePlus outTiff = new ImagePlus();
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			
			nowTiff.setPosition(z + 1);
			
			if (z < numberOfAddedLayers + 1) {
				distTiff.setPosition(2);
			}
			else {
				if (z > numberOfAddedLayers + distTiff.getNSlices()) distTiff.setPosition(distTiff.getNSlices() - 1);
				else distTiff.setPosition(z - numberOfAddedLayers);			
			}
		
			if (z >= numberOfAddedLayers + distTiff.getNSlices()) distTiff.setPosition(distTiff.getNSlices());
			
			ImageProcessor nowIP = nowTiff.getProcessor();
			ImageProcessor distIP = distTiff.getProcessor();
			ImageProcessor outIP = new FloatProcessor(distIP.getWidth(), distIP.getHeight());
			
			for (int x = 0 ; x < nowIP.getWidth() ; x++) {
				for (int y = 0 ; y < nowIP.getHeight() ; y++) {
					
					float nowSkel = nowIP.getPixelValue(x, y);
					float nowDist = distIP.getPixelValue(x, y);
					
					if (nowSkel == 255) {
						
						if (nowDist > 0) outIP.putPixelValue(x, y, nowDist);
						else outIP.putPixelValue(x, y, 1); 
						
						int cc = 0;
						
					}										
				}				
			}	
			
			outStack.addSlice(outIP);
			
		}
		
		outTiff.setStack(outStack);
		
		//outTiff.updateAndDraw();
		//outTiff.show();
	
		return outTiff;
		
	}
	
	public ImagePlus selectSlicesFromTiff(ImagePlus nowTiff, int[] loadSlices) {
		
		ImagePlus outTiff = new ImagePlus();
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		for (int z = 0 ; z < loadSlices.length ; z++) {
			
			nowTiff.setPosition(loadSlices[z]);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			outStack.addSlice(nowIP);			
		}
		
		outTiff.setStack(outStack);
		
		return outTiff;
		
	}
	
	public ImagePlus calibrateGrayValues(InputOutput.MyFileCollection mFC, ImagePlus nowTiff, MenuWaiter.CalibrationReferences myNR) {
		
		//init units
		RoiHandler roi = new RoiHandler();
		HistogramStuff hist = new HistogramStuff();
		TailoredMaths maths = new TailoredMaths();		
		InputOutput jIO = new InputOutput();			

		//init vars
		ImagePlus outTiff = new ImagePlus();								
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());						
		int[] myHist = new int[256 * 256];
		int[] wallHist = new int[256 * 256];
		int[] newHist;
		int i, j;
		
		PolygonRoi[] pRoi = null;
		PolygonRoi[] oRoi  = null;
		PolygonRoi[] ooRoi  = null;
		
		//read gauge file
		if (myNR.useInnerCircle) {
			
			//read InnerCircle file
			ObjectDetector jOD = new ObjectDetector();
			
			if (myNR.material.equalsIgnoreCase("aluminium")) {
				
				ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();			
				int versio = jIO.checkInnerCircleFileVersion(mFC.nowInnerCirclePath);			
				if (versio == 0) jCO = jIO.readInnerCircleVer0(mFC.nowInnerCirclePath);	
				else jCO = jIO.readInnerCircleVer1(mFC.nowInnerCirclePath);	
				
				int wallThickness = (int) (0.9 * Math.round(StatUtils.max(jCO.wallThickness)));				
			
				pRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 4);		
				oRoi = roi.makeMeAPolygonRoiStack("outer", "manual", jCO, 4);
				ooRoi = roi.makeMeAPolygonRoiStack("outerFromInner", "manual", jCO, -wallThickness-6);
			}
			
			else {

				//read column coordinates
				ObjectDetector.EggShapedColCoords3D jCO = jOD.new EggShapedColCoords3D();						
				jCO = jIO.readInnerCircleSteel(mFC);
				
				pRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 4);		
				oRoi = roi.makeMeAPolygonRoiStack("outer", "manual", jCO, 4);
				ooRoi = roi.makeMeAPolygonRoiStack("outerFromInner", "manual", jCO, -(int)StatUtils.mean(jCO.wallThickness)-6);
				
			}
		}
				
		//find illumination of this column	
		double[] lowerQuantile = new double[nowTiff.getNSlices()];
		double[] upperQuantile = new double[nowTiff.getNSlices()];
		double[] mode = new double[nowTiff.getNSlices()];
		double[] wall = new double[nowTiff.getNSlices()];	
	
		//get the stacked histogram		
		for (i = 1 ; i < nowTiff.getNSlices() + 1 ; i++) {
			
			IJ.showStatus("Getting 16-bit histogram of slice #" + i + "/" + nowTiff.getNSlices());
			
			nowTiff.setPosition(i);		
			
			ImageProcessor myIP = nowTiff.getProcessor();
			ImageProcessor modIP = myIP.duplicate();
			ImageProcessor outIP = myIP.duplicate();
			ImageProcessor faroutIP = myIP.duplicate();
			ImageProcessor wallIP = myIP.duplicate();
			
			//cut away everything outside column
			if (myNR.useInnerCircle) {
				modIP.setRoi(pRoi[i - 1]);			
				modIP.setColor(0);
				modIP.fillOutside(pRoi[i - 1]);
			}
			
			//get histogram of soil
			newHist=modIP.getHistogram();
			newHist[0] = 0; //set zero GV to zero
			for (j = 0 ; j < myHist.length ; j++) {
				myHist[j] = myHist[j] + newHist[j];
			}
			double[] inCumHist = hist.calcCumulativeHistogram(newHist);
			
			//find mode  (2017/02/15: presently disabled)
			if (myNR.hiRef.equalsIgnoreCase("mode")) {
				mode[i-1] = (double)hist.findModeFromHistogram(newHist);
			}
			
			//cut away soil column including walls and another wall distance of the surroundings
			if (myNR.useInnerCircle) {
				faroutIP.setRoi(ooRoi[i - 1]);			
				faroutIP.setColor(0);
				faroutIP.fill(ooRoi[i - 1]);	
				faroutIP.resetRoi();
			}
			
			//get histogram of soil
			newHist=faroutIP.getHistogram();
			newHist[0] = 0; //set zero GV to zero
			for (j = 0 ; j < myHist.length ; j++) {
				myHist[j] = myHist[j] + newHist[j];
			}
			double[] faroutCumHist = hist.calcCumulativeHistogram(newHist);
			
			//cut away soil column including walls and everything outside another wall distance from the outer wall perimeter
			if (myNR.useInnerCircle) {
				outIP.setRoi(oRoi[i - 1]);				
				outIP.setColor(0);
				outIP.fill(oRoi[i - 1]);	
				outIP.resetRoi();
				
				outIP.setRoi(ooRoi[i - 1]);
				outIP.fillOutside(ooRoi[i - 1]);
				outIP.resetRoi();
			}
			
			//get histogram of soil
			newHist=outIP.getHistogram();
			newHist[0] = 0; //set zero GV to zero
			for (j = 0 ; j < myHist.length ; j++) {
				myHist[j] = myHist[j] + newHist[j];
			}
			double[] outCumHist = hist.calcCumulativeHistogram(newHist);
			
			//find reference quantiles
			if (myNR.lowRef.equalsIgnoreCase("quantile")) {
				lowerQuantile[i-1] = hist.findPercentileFromCumHist(inCumHist, myNR.lowerReference); 
				if (myNR.sampleLowerWithinSoil) lowerQuantile[i-1] = hist.findPercentileFromCumHist(inCumHist, myNR.lowerReference); 
				else {
					if (myNR.lowerTag.equalsIgnoreCase("Outside")) lowerQuantile[i-1] = hist.findPercentileFromCumHist(outCumHist, myNR.lowerReference);
					else lowerQuantile[i-1] = hist.findPercentileFromCumHist(faroutCumHist, myNR.lowerReference);
				}
			}		
			
			if (myNR.hiRef.equalsIgnoreCase("quantile")) {				
				upperQuantile[i-1] = hist.findPercentileFromCumHist(inCumHist, myNR.upperReference);
				if (myNR.sampleUpperWithinSoil) upperQuantile[i-1] = hist.findPercentileFromCumHist(inCumHist, myNR.upperReference);
				else {
					if (myNR.upperTag.equalsIgnoreCase("Outside")) upperQuantile[i-1] = hist.findPercentileFromCumHist(outCumHist, myNR.upperReference);
					else upperQuantile[i-1] = hist.findPercentileFromCumHist(faroutCumHist, myNR.upperReference);
				}
			}
			
			//cut out everything but the wall
			if (myNR.useInnerCircle) {
				wallIP.setRoi(pRoi[i - 1]);			
				wallIP.setColor(0);
				wallIP.fill(pRoi[i - 1]);
			
				wallIP.setRoi(oRoi[i - 1]);
				wallIP.setColor(0);
				wallIP.fillOutside(oRoi[i - 1]);
			}
			
			//get histogram of wall
			newHist=wallIP.getHistogram();
			newHist[0] = 0; //set zero GV to zero
			for (j = 0 ; j < wallHist.length ; j++) {
				wallHist[j] = newHist[j];
				myHist[j] = myHist[j] + newHist[j];
			}	

			//DisplayThings dT = new DisplayThings();
			//dT.showMeMyRoi("pRoi, slice" + i, wallIP, pRoi[i - 1], 1);
					
			//find median wall gray			
			double[] wallCumHist = hist.calcCumulativeHistogram(wallHist);			
			wall[i-1] = hist.findPercentileFromCumHist(wallCumHist, 0.5);		
			
		/*	//normalize lower reference
			if (myNR.lowRef.equalsIgnoreCase("quantile")) {
				if (myNR.upperReference > 0) normLower[i - 1] = lowerQuantile[i-1] / upperQuantile[i-1];
				else normLower[i - 1] = lowerQuantile[i-1] / wall[i-1];				
			}
			else {
				normLower[i - 1] = wall[i-1] / upperQuantile[i-1];					
			}*/
			  
		}		
		
		//write results into file
		int windowHalfSize = 50;  //window half-size for LOESS filter
		if (myNR.lowRef.equalsIgnoreCase("quantile")) {					
			double[] smoothLower = maths.filterAnArray(lowerQuantile, myNR.windowSize, myNR.lowerSmoothingFilter);
			smoothLower = maths.fillStartAndEndZeros(smoothLower);
			myNR.originalLower = maths.LinearLOESSFilter(smoothLower, windowHalfSize);
		}
		if (myNR.lowRef.equalsIgnoreCase("wall")) {
			double[] smoothLower = maths.filterAnArray(wall, myNR.windowSize, myNR.lowerSmoothingFilter);
			smoothLower = maths.fillStartAndEndZeros(smoothLower);
			myNR.originalLower = maths.LinearLOESSFilter(smoothLower, windowHalfSize);	
		}
		
		if (myNR.hiRef.equalsIgnoreCase("quantile")) {
			double[] smoothUpper = maths.filterAnArray(upperQuantile, myNR.windowSize, myNR.upperSmoothingFilter);
			smoothUpper = maths.fillStartAndEndZeros(smoothUpper);
			myNR.originalUpper = maths.LinearLOESSFilter(smoothUpper, windowHalfSize);	
		}
		if (myNR.hiRef.equalsIgnoreCase("wall")) {
			double[] smoothUpper = maths.filterAnArray(wall, myNR.windowSize, myNR.upperSmoothingFilter);
			smoothUpper = maths.fillStartAndEndZeros(smoothUpper);
			myNR.originalUpper = maths.LinearLOESSFilter(smoothUpper, windowHalfSize);	
		}
			
		//String myFileName = nowTiff.getTitle().substring(0, nowTiff.getTitle().length() - 4);
		//jIO.writeIlluminationCorrectionIntoAsciiFile(myNR, mFC.myOutFolder, myFileName);
		
		//apply correction image
		for (i = 0 ; i < nowTiff.getNSlices() ; i++) {
			
			IJ.showStatus("Normalizing slice #" + i + "/" + nowTiff.getNSlices());
			
			nowTiff.setPosition(i + 1);
			
			ImageProcessor myIP = nowTiff.getProcessor();
			
			double nowIntercept = myNR.lowerTarget;
			double nowSlope = (myNR.upperTarget - myNR.lowerTarget) / (myNR.originalUpper[i] - myNR.originalLower[i]);
						
			for (int x = 0 ; x < myIP.getWidth() ; x++) {
				for (int y = 0 ; y < myIP.getHeight() ; y++) {
					
					double nowX = myIP.getPixelValue(x, y) - myNR.originalLower[i];
					double newPixelValue = nowIntercept + nowSlope * nowX;					
									
					myIP.putPixelValue(x, y, newPixelValue);
					
				}
			}
			
			outStack.addSlice(myIP);

		}
		
		outTiff.setStack(outStack);
		
		return outTiff;
		
	}
	
	
	
	
	///////////////////////////////////////////////
	// gimmicks
	//////////////////////////////////////////////
	
	public ImageProcessor invertSoilBinary(ImageProcessor myIP, Roi pRoi) {
	
		ImageProcessor modIP = myIP.duplicate();
		
		modIP.invert();
		
		if (pRoi != null) {
			modIP.setRoi(pRoi);
			modIP.setColor(0);
			modIP.fillOutside(pRoi);
			modIP.resetRoi();
		}
		
		return modIP;
				
	}
	
	public ImagePlus invertImage(ImagePlus nowTiff, PolygonRoi[] pRoi) {
		
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {			
			
			nowTiff.setPosition(z + 1);
			
			PolygonRoi nowRoi = null;
			if (pRoi != null) nowRoi = pRoi[z];				
			ImageProcessor nowIP = invertSoilBinary(nowTiff.getProcessor().duplicate(), nowRoi);
			
			outStack.addSlice(nowIP);
		}		
		
		ImagePlus outTiff = new ImagePlus("", outStack);
		
		return outTiff;
		
	}
	
	public ImagePlus flipTiff(ImagePlus nowTiff) {
		
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImagePlus outTiff = new ImagePlus();
		
		int i;
		
		for (i = nowTiff.getNSlices() ; i > 0 ; i--) {
			
			nowTiff.setPosition(i);
			
			ImageProcessor myIP = nowTiff.getProcessor();
			
			outStack.addSlice(myIP);			
		}
		outTiff.setStack(outStack);
			
		return outTiff;
		
	}
	
	public ImagePlus cutOutRealThickness(ImagePlus nowTiff, ImagePlus rawThickImp) {
		
		ImagePlus outTiff = new ImagePlus();
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		//rawThickImp.updateAndDraw();
		//rawThickImp.show();
		
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			nowTiff.setPosition(z+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			rawThickImp.setPosition(z+1);
			ImageProcessor rawIP = rawThickImp.getProcessor();
			ImageProcessor outIP = rawIP.duplicate();
			for (int x = 0 ; x < nowTiff.getWidth() ; x++) {		
				for (int y = 0 ; y < nowTiff.getHeight() ; y++) {
					int nowPix = nowIP.getPixel(x, y);
					if (nowPix == 0) outIP.putPixel(x, y, 0);
				}
			}
			outStack.addSlice(outIP);
		}
		
		outTiff.setStack(outStack);
		
		//outTiff.updateAndDraw();
		//outTiff.show();
		
		return outTiff;
	}
	
	public ImagePlus binarizeLargerThanZero(ImagePlus nowTiff) {
		
		ImagePlus binTiff = new ImagePlus();
		ImageStack binStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			
			nowTiff.setPosition(z + 1);
			
			ImageProcessor binIP = nowTiff.getProcessor().duplicate().convertToByte(false);
			
			for (int x = 0 ; x < binIP.getWidth() ; x++) {
				for (int y = 0 ; y < binIP.getHeight() ; y++) {
					
					int nowPix = binIP.get(x, y);
					if (nowPix > 0) binIP.putPixel(x, y, 255);
					
				}
			}
			
			binStack.addSlice(binIP);
		}
		
		binTiff.setStack(binStack);
		
		return binTiff;		
	}
	
	public ImagePlus binarizeGradientMask(ImagePlus nowTiff, int myThresh) {
		
		ImagePlus outTiff = new ImagePlus();
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());	
		
		//create mask for cutting out the rest..
		for (int i = 0 ; i < nowTiff.getNSlices() ; i++) {
			nowTiff.setPosition(i+1);		
		
			ImageProcessor binIP = nowTiff.getProcessor();
			binIP.threshold(myThresh);
			//binIP.dilate();			
			binIP.multiply(1/255.0);
			
			outStack.addSlice(binIP);
		}
		outTiff.setStack(outStack);
		
		return outTiff;		
	}
	
	public ImagePlus makeASubstack(ImagePlus nowTiff, int topSlice, int botSlice) {
		
		ImagePlus outTiff = new ImagePlus();
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		for (int i = topSlice ; i < botSlice ; i++) {
			
			nowTiff.setPosition(i + 1);
			ImageProcessor nowIP = nowTiff.getProcessor();
						
			outStack.addSlice(nowIP);			
		}
		
		outTiff.setStack(outStack);
		return outTiff;
		
	}

	public ImageProcessor fillHolesWithNeighboringValues(ImageProcessor surfaceIP, PolygonRoi[] pRoi, String topOrBottom) {
	
		DisplayThings disp = new DisplayThings();
		
		//segment out deep roots			
		ImageProcessor justHoles = surfaceIP.duplicate();
		justHoles.threshold(1);
		justHoles.invert();
		
		justHoles.setBackgroundValue(0);
		if (topOrBottom.equalsIgnoreCase("top")) justHoles.fillOutside(pRoi[0]);
		else justHoles.fillOutside(pRoi[pRoi.length - 1]);
		ImageProcessor justHoles8Bit = justHoles.convertToByte(false);
		
		//identify the locations 
		ImagePlus holeOne = new ImagePlus("", justHoles8Bit);
		int imgVol = holeOne.getWidth() * holeOne.getHeight();
		Counter3D myOC = new Counter3D(holeOne, 128, 0, imgVol, false, false);
		myOC.getObjects();		
		ImagePlus myHoles = myOC.getObjMap();
	
		Vector allObjects = myOC.getObjectsList();
		int numOfObjects = allObjects.size();
		
		//fill Holes
		ImageProcessor outIP = surfaceIP.duplicate();
		for (int i = 0; i < numOfObjects; i++){		         
			Object3D currObj=(Object3D) allObjects.get(i);
			
			int[] tmpArrayInt=currObj.bound_cube_TL;
			
			//fill holes that were spotted
			int xmin = tmpArrayInt[0] - 1;
			int xmax = xmin + currObj.bound_cube_width + 1;
			int ymin = tmpArrayInt[1] - 1;
			int ymax = ymin + currObj.bound_cube_height + 1;
			
			float[] X = new float[]{xmin, xmin, xmax, xmax, xmin};
			float[] Y = new float[]{ymin, ymax, ymax, ymin, ymin};
			
			PolygonRoi mRoi = new PolygonRoi(X, Y, Roi.POLYGON);
	
			//init list for hole-neighboring pixels
			ArrayList<Integer> neighbors = new ArrayList<Integer>();
			
			//cut out canvas around hole #i, remove all holes with other IDs and binarize the image.
			ImageProcessor hoIP = myHoles.getProcessor();
			
			ImageProcessor thIP = new ByteProcessor(xmax - xmin + 1, ymax - ymin + 1); 
			for (int x = xmin ; x < xmax + 1 ; x++) {
				for (int y = ymin ; y < ymax + 1 ; y++) {	
					int nowVal = (int)Math.round(hoIP.getPixelValue(x, y));
					if (nowVal == i) thIP.putPixel(x - xmin, y - ymin, 255); 		
					else thIP.putPixel(x, y, 0); 	
				}
			}
			
			// do the same for the surface elevation map
			ImageProcessor vlIP = surfaceIP.duplicate();
			vlIP.setRoi(mRoi);
			vlIP.crop();
			
			//make a copy of the binarized hole, dilate it and sample gray value of neighbors
			ImageProcessor cpIP = thIP.duplicate();
			cpIP.erode();
	
			for (int x = 0 ; x < cpIP.getWidth() ; x++) {
				for (int y = 0 ; y < cpIP.getHeight() ; y++) {						
					int nowTH = thIP.getPixel(x, y);
					int nowCP = cpIP.getPixel(x, y);
					int nowVL = (int)Math.round(vlIP.getPixelValue(x + xmin, y + ymin));
					if (nowTH == 0 & nowCP > 0 & nowVL > 0) {
						neighbors.add(nowVL);
					}
				}
			}
	
			//calculate median neighbor value
			double[] neighborsAsArray = new double[neighbors.size()]; 
			for (int j = 0 ; j < neighbors.size() ; j++) neighborsAsArray[j] = neighbors.get(j);
			int medianNeighbor = (int)StatUtils.percentile(neighborsAsArray, 10);
			
			//fill hole with neighborvalue
			for (int x = 0 ; x < cpIP.getWidth() ; x++) {
				for (int y = 0 ; y < cpIP.getHeight() ; y++) {
					int nowTH = thIP.getPixel(x, y);
					if (nowTH > 0) {
						outIP.putPixel(x + xmin, y + ymin, medianNeighbor);
					}
				}
			}				
		}
		
		return outIP;
	}
	
	public ImagePlus clearOutside(ImagePlus nowTiff, PolygonRoi[] pRoi) {
		
		ImagePlus outTiff = new ImagePlus();
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
				
		//get the stacked histogram		
		for (int i = 1 ; i < nowTiff.getNSlices() + 1 ; i++) {
			
			nowTiff.setPosition(i);		
			
			ImageProcessor myIP = nowTiff.getProcessor();
			ImageProcessor modIP = myIP.duplicate();
			
			//cut out everything outside column			
			modIP.setRoi(pRoi[i - 1]);
			modIP.setColor(0);
			modIP.fillOutside(pRoi[i - 1]);
		
			outStack.addSlice(modIP);			
		}
		
		outTiff.setStack("", outStack);		
		
		return outTiff;
	}
	
	public ImagePlus binaryScale2HalfSize(ImagePlus binTiff) {
		
		ImagePlus smallTiff = new ImagePlus();
		ImageStack smallStack = new ImageStack(binTiff.getWidth() / 2, binTiff.getHeight() / 2);
		
		for (int z = 0 ; z < binTiff.getNSlices() - 1; z += 2) {
			
			binTiff.setSlice(z + 1);
			ImageProcessor zIP1 = binTiff.getProcessor(); 
			binTiff.setSlice(z + 2);
			ImageProcessor zIP2 = binTiff.getProcessor();
			
			ImageProcessor outIP = new ByteProcessor(smallStack.getWidth(), smallStack.getHeight());
			
			for (int x = 0 ; x < binTiff.getWidth() - 1 ; x += 2) {
				for (int y = 0 ; y < binTiff.getHeight() - 1 ; y += 2) {
					
					int black = 0;
					int white = 0;
					
					int nowPixel = 128;
					
					for (int i = 0 ; i < 2 ; i++) {
						for (int j = 0 ; j < 2 ; j++) {
							nowPixel = zIP1.get(x + i, y + j);
							if (nowPixel == 0) black++;
							if (nowPixel == 255) white++;
						}
					}
					
					for (int i = 0 ; i < 2 ; i++) {
						for (int j = 0 ; j < 2 ; j++) {
							nowPixel = zIP2.get(x + i, y + j);
							if (nowPixel == 0) black++;
							if (nowPixel == 255) white++;
						}
					}
					
					if (black > white) outIP.putPixel(x / 2, y / 2, 0);
					else outIP.putPixel(x / 2, y / 2, 255);
					
				}
			}
			
			smallStack.addSlice(outIP);
			
		}
		
		smallTiff.setStack(smallStack);
		
		return smallTiff;
		
	}

	//create correction maps
	public ImageProcessor createBeamHardeningCorrectionMap4Steel(ObjectDetector.RadialFacts rF) {
		
		TailoredMaths math = new TailoredMaths();
		
		ImageProcessor corrMap = new FloatProcessor(rF.standardRadius, rF.imageHeight);
		double[] smoothInZ = new double[rF.imageHeight];
		
		for (int x = 0 ; x < rF.standardRadius ; x++) {
			
			//apply 1-D median filter in vertical direction
			for (int z = 0 ; z < rF.imageHeight ; z++) smoothInZ[z] = rF.smoothedRadialProfile[z][x];
			smoothInZ = math.oneDMeanFilter(smoothInZ, 15);
			
			//assign values to correction map..
			for (int z = 0 ; z < rF.imageHeight ; z++) corrMap.putPixelValue(x, z, smoothInZ[z]);
		}
		
		//ImagePlus test = new ImagePlus("",corrMap);
		//test.updateAndDraw();test.show();
		
		return corrMap;
		
	}
	
	/*public ImagePlus scaleIsotropicly3D(ImagePlus nowTiff, double scalingFactor) {
		
		//probe size of scaled Image
		nowTiff.setSlice(1);
		ImageProcessor nowIP = nowTiff.getProcessor().duplicate(); 
		nowIP.setInterpolationMethod(ImageProcessor.BILINEAR);
		nowIP.scale(scalingFactor, scalingFactor);
		
		//init intermediate and output images.
		ImagePlus smallTiff = new ImagePlus();
		ImagePlus zwischiTiff = new ImagePlus();
		ImageStack smallStack = new ImageStack(nowIP.getWidth(), nowIP.getHeight());
		ImageStack zwischiStack = new ImageStack(nowIP.getWidth(), nowIP.getHeight());
		
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			
			nowTiff.setSlice(z + 1);
			nowIP = nowTiff.getProcessor().duplicate(); 
			
			nowIP.setInterpolationMethod(ImageProcessor.BILINEAR);
			nowIP.sc(scalingFactor, scalingFactor);
			
			zwischiStack.addSlice(nowIP);
			
		}
		
		zwischiTiff.setStack(smallStack);
		
		for (int z = 0 ; z < Math.floorDiv(zwischiTiff.getNSlices(),2) - 1 ; z += 2) {
			
			zwischiTiff.setSlice(z + 1);
			ImageProcessor nowIP1 = zwischiTiff.getProcessor().duplicate();
			
			zwischiTiff.setSlice(z + 2);
			ImageProcessor nowIP2 = zwischiTiff.getProcessor().duplicate();
			
			nowIP = nowIP1.duplicate();
			
			for (int x = 0 ; x < nowIP1.getWidth() ; x++) {
				for (int y = 0 ; y < nowIP1.getHeight() ; y++) {
					
					int p1 = nowIP1.getPixel(x, y);
					int p2 = nowIP2.getPixel(x, y);
					
					nowIP.putPixel(x, y, Math.floorDiv(p1 + p2, 2));
					
				}
			}
			
			smallStack.addSlice(nowIP);
			
		}
			
		smallTiff.setStack(smallStack);
		
		return smallTiff;
		
	}*/
	
}
	


	


