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
import java.util.ArrayList;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import SoilJ.tools.HistogramStuff.IlluminationInfo;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.measure.CurveFitter;
import ij.measure.SplineFitter;
import ij.plugin.ContrastEnhancer;
import ij.plugin.PlugIn;
import ij.plugin.Selection;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.filter.RankFilters;
import ij.plugin.filter.ThresholdToSelection;
import ij.process.ByteProcessor;
import ij.process.EllipseFitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import inra.ijpb.binary.BinaryImages;
import inra.ijpb.binary.conncomp.FloodFillComponentsLabeling3D;
import inra.ijpb.label.LabelImages;
import inra.ijpb.morphology.Morphology;
import inra.ijpb.morphology.Strel3D;
import inra.ijpb.morphology.strel.CubeStrel;
import process3d.Dilate_;
import process3d.Erode_;

/** 
 * ObjectDetector is one of the two main SoilJ classes (the other being ImageManipulator). 
 * Object detector contains subroutines for object recognition.
 * 
 * @author John Koestel
 *
 */

public class ObjectDetector implements PlugIn {

	
	public void run(String arg) {
				//ok, this is not needed..
	}
	
	public static final long MilliSecondsOfOneDay=86400000L;
	
	public class RadialFacts {
		
		public double[][] radialProfile;		
		public double[][] smoothedRadialProfile;	
		public double[] min;		
		public double[] max;
		public int standardRadius;
		public int anglesChecked;
		public int imageHeight;
		
	}
	
	public class ColCoords2D {
		
		//outer coords
		public double xCenter;
		public double yCenter;
		public double zCenter;
		public double outerMinorRadius;
		public double outerMajorRadius;				
		public double theta;		
		public double outerR2;
		
		//inner coords
		public double ixCenter;
		public double iyCenter;
		public double innerMinorRadius;
		public double innerMajorRadius;				
		public double itheta;		
		public double innerR2;
		
		//wall properties
		public double wallThickness;
		public double CVOfWallThickness;
		public double CVOfWallBrightness;
		public double ratioBetweenOuterAndInnerRadius;
		
		//global props
		public boolean columnIsAtThisDepth;		
	
		//fitting parameters
		public double airWallContrast;	
		public double absoluteWallgrayValue;
		public double wallSoilStdContrastThreshold;		
				
	}
	
	public class ColCoordsEssentials3D {
		
		public int anglesChecked;
		
		public double[][] xID;
		public double[][] yID;
		
		public double[] xmid;
		public double[] ymid;
		
	}
	
	public class BestWallFinds {
		
		public float[] xD;
		public float[] yD;
		public double[] angle;		
		
	}
	
	public class SampleBrightness {
		
		int[][] histogram;
		int[][] centralHistogram;
		int[][] referenceHistogram;
		
		double[] wall;
		
		int[] q01;
		int[] q50;
		int[] q60;
		int[] q80;
		int[] q95;
		int[] central_q01;
		int[] central_q50;
		int[] central_q60;
		int[] central_q80;
		int[] central_q95;
		int[] ref_q01;
		int[] ref_q50;
		int[] ref_q60;
		int[] ref_q80;
		int[] ref_q95;
				
	}
	
	public class ColumnContainer {
		
		public ImagePlus nowTiff;
		public ImagePlus cutTiff;
		public ImagePlus normalizedTiff;
		public ImagePlus binaryTiff;
		
		public int startslice;
		public int stopslice;
		public int increment;		
		public int dAlpha;
		
		public boolean autoTuned;
		
		public ColCoords3D prelimCC;
		public ColCoords3D preciseCC;
		public ColCoords3D mendedCC;
			
		public MenuWaiter.ColumnFinderMenuReturn jCFS;
		
	}
	
	public class ApproximateColumnIllumination {
		
		int outside;
		int wall;
		int inside;
		int globalContrast;
		
	}
	
	public class RadialModes {
		
		double[][] maskedRadialModes;
		double[][] maskedRadialMinima;
		double[] radius;
		int[] maskingThreshold;
		
	}
	
	public class EggShapedColCoords3D {
		
		public double[] xmid;			//x midpoint
		public double[] ymid;			//y midpoint
		public double[] zmid;			//y midpoint
		
		public double[] ixmid;			//x midpoint (inner circle)
		public double[] iymid;			//y midpoint (inner circle)
		
		public double[][] xOD;
		public double[][] yOD;
		public double[][] xID;
		public double[][] yID;
		
		public double[] medianOuterDiameter;
		public double[] innerRadius;
		
		public double[] wallThickness;
		public double[] theta;  //angle of major ellipse axis
		public double[] itheta;  //angle of major ellipse axis (inner circle)
		
		public double[] outerR2;
		public double[] innerR2;
		public boolean[] columnIsAtThisDepth;
		
		public int bevelStartsAt;
		public FitStuff.LinearFunction bevelFunction;
			
		public double tiltInXZ; 
		public double tiltInYZ;
		public double tiltTotal;
			
		public int topOfColumn;
		public int heightOfColumn;
		public int bottomOfColumn;
		public int numberOfImputedLayers;
		
		//fitting parameters		
		public int anglesChecked;
		
		public double[] airWallContrast;
		public double[] absoluteWallgrayValue;
		public double[] wallSoilStdContrastThreshold;		
				
	}
	
	public class ColCoords3D {
	
		public double[] xmid;			//x midpoint
		public double[] ymid;			//y midpoint
		public double[] zmid;			//y midpoint
		
		public double[] ixmid;			//x midpoint (inner circle)
		public double[] iymid;			//y midpoint (inner circle)
		
		public double[] outerMajorRadius;
		public double[] innerMajorRadius;
		public double[] outerMinorRadius;
		public double[] innerMinorRadius;
		public double[] deltaMajorMinorRadius;
		public double[] wallThickness;
		public double[] theta;  //angle of major ellipse axis
		public double[] itheta;  //angle of major ellipse axis (inner circle)
		
		public double[] outerR2;
		public double[] innerR2;
		public boolean[] columnIsAtThisDepth;
		
		public int bevelStartsAt;
		public FitStuff.LinearFunction bevelFunction;
			
		public double tiltInXZ; 
		public double tiltInYZ;
		public double tiltTotal;
			
		public int topOfColumn;
		public int heightOfColumn;
		public int bottomOfColumn;
		public int numberOfImputedLayers;
		public int topCutOff;
		public int botCutOff;
		
		//fitting parameters		
		public int anglesChecked;
		
		public double[] airWallContrast;
		public double[] absoluteWallgrayValue;
		public double[] wallSoilStdContrastThreshold;		
				
	}
	
	public RadialFacts getRadialFacts(int imageHeight, int anglesChecked, int standardRadius, double[][][] radialGrayValues, double percentile) {
		
		RadialFacts myRF = new RadialFacts();
		double[][] radialProfile = new double[imageHeight][standardRadius];
		double[] meanProfileAtThisHeight = new double[standardRadius];
		double[] nowRadii = new double [anglesChecked];
		double[] maxi = new double[imageHeight];
		double[] mini = new double[imageHeight];
		
		myRF.anglesChecked = anglesChecked;
		myRF.standardRadius = standardRadius;
		myRF.imageHeight = imageHeight;
		
		for (int z = 0 ; z < imageHeight ; z++) {
			
			for (int r = 0 ; r < standardRadius ; r++) {
				
				for (int alpha = 0 ; alpha < anglesChecked ; alpha++) {
					
					nowRadii[alpha] = radialGrayValues[z][alpha][r];
					
				}
				
				radialProfile[z][r] = StatUtils.percentile(nowRadii, percentile);
				meanProfileAtThisHeight[r] = radialProfile[z][r]; 				
			}		
			
			maxi[z] = StatUtils.max(meanProfileAtThisHeight);
			mini[z] = StatUtils.min(meanProfileAtThisHeight);
			
		}
		
		//assign results to output structure
		myRF.max = maxi;
		myRF.min = mini;
		myRF.radialProfile = radialProfile;
			
		return myRF;
		
	}	
	
	public double[][][] getRadialSteelGrayValues(ImagePlus nowTiff, EggShapedColCoords3D jCO, float radialMappingFactor) {
		
		RollerCaster rC = new RollerCaster();
		DisplayThings disp = new DisplayThings();
		
		//init variables	
		int standardRadius = (int)Math.round(radialMappingFactor * (StatUtils.percentile(jCO.innerRadius, 50)));
		//double[] sRadius = new double[standardRadius];
		double[][][] radialgrayValues = new double[jCO.heightOfColumn][jCO.anglesChecked][standardRadius];
		int maxAlpha = 360;
		int dAlpha = maxAlpha / jCO.anglesChecked;
		
		
		for (int i = 0 ; i < nowTiff.getNSlices() ; i++) {
			
			IJ.showStatus("Sampling radial illumination of slice #" + (i + 1) + "/" + nowTiff.getNSlices());
			
			//set image to next slice
			nowTiff.setPosition(i+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			//sweep through all checked angles and get radial gray values
			int angleCounter = 0;			
			for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {				
								
				double dx = jCO.xID[i][angleCounter] - jCO.xmid[i];
				double dy = jCO.yID[i][angleCounter] - jCO.ymid[i];
				double nowRadius = Math.sqrt((double)(dx*dx) + (double)(dy*dy));
				int[] x = new int[(int)Math.round(nowRadius)];
				int[] y = new int[(int)Math.round(nowRadius)];
				double[] grayAtThisAngle = new double[(int)Math.round(nowRadius)];
				//double[] splinedGrayAtThisAngle = new double[standardRadius];
				int[] radius = new int[(int)Math.round(nowRadius)];
				
				int distanceCounter = 0;
				for (double checkRadius = 0 ; checkRadius < (int)Math.round(nowRadius); checkRadius++) {					
					x[distanceCounter] = (int)Math.round(Math.sin(angle) * checkRadius + jCO.xmid[i]);
					y[distanceCounter] = (int)Math.round(Math.cos(angle) * checkRadius + jCO.ymid[i]);				
					grayAtThisAngle[distanceCounter] = (int)Math.round(nowIP.getPixelValue(x[distanceCounter], y[distanceCounter]));	
					radius[distanceCounter]	= (int)Math.round(checkRadius);
					distanceCounter++;					
				}	
				
				//stretch gray values to standardized radius
				float[] newRadius = new float[distanceCounter];
				for (int j = 0 ; j < distanceCounter ; j++) newRadius[j] = j * standardRadius / distanceCounter;
				SplineFitter sF = new SplineFitter(newRadius, rC.castDouble2Float(grayAtThisAngle), radius.length);
				for (int j = 0 ; j < standardRadius ; j++) {
					double newgray = sF.evalSpline(j);
					//sRadius[j] = j;
					//splinedGrayAtThisAngle[j] = newgray;
					radialgrayValues[i][angleCounter][j] = newgray;
				}
				
				//display radial gray value profile
				//disp.plotXYXY(rC.castInt2Double(radius), grayAtThisAngle, sRadius, splinedGrayAtThisAngle, "radius", "GV", "splined GV");
				
				angleCounter++;
				 				
			}
			
		}
		
		return radialgrayValues;
	}
	
	public double[][][] getRadialAluGrayValues(ImagePlus nowTiff, ColCoords3D jCO, float radialMappingFactor) {
		
		RollerCaster rC = new RollerCaster();
		DisplayThings disp = new DisplayThings();
		TailoredMaths math = new TailoredMaths();
	
		//init variables
		if (jCO.anglesChecked == 0) jCO.anglesChecked = 72;  //standard value..
		jCO.heightOfColumn = nowTiff.getNSlices();
		
		int standardRadius = (int)Math.round(radialMappingFactor * (StatUtils.percentile(jCO.innerMinorRadius, 50)));
		//double[] sRadius = new double[standardRadius];
		double[][][] radialgrayValues = new double[jCO.heightOfColumn][jCO.anglesChecked][standardRadius];
		int maxAlpha = 360;
		
		int dAlpha = maxAlpha / jCO.anglesChecked;		
		
		for (int i = 0 ; i < nowTiff.getNSlices() ; i++) {
			
			IJ.showStatus("Sampling radial illumination of slice #" + (i + 1) + "/" + nowTiff.getNSlices());
			
			//set image to next slice
			nowTiff.setPosition(i+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
						
			//create array of considered angles 
			int angleCounter = 0;			
			double[] myAngles = new double[jCO.anglesChecked];
			for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {				
				
				myAngles[angleCounter] = angle;
				
				angleCounter++;
			}
			
			//calculate x and y coordinate
			double[][] xy = math.getXYOfEllipseFromAngle(myAngles, jCO.ixmid[i], jCO.iymid[i], jCO.innerMajorRadius[i], jCO.innerMinorRadius[i], jCO.itheta[i]);
					
			//cast it into two arrays..
			double[] x = new double[jCO.anglesChecked];
			double[] y = new double[jCO.anglesChecked];
			for (int j = 0 ; j < myAngles.length ; j++) x[j] = xy[j][0];
			for (int j = 0 ; j < myAngles.length ; j++) y[j] = xy[j][1]; 
				
			//sweep through all checked angles and get radial gray values
			for (int alpha = 0 ; alpha < myAngles.length ; alpha++) {
				
				double dx = x[alpha] - jCO.xmid[i];
				double dy = y[alpha] - jCO.ymid[i];
				double nowRadius = Math.sqrt((double)(dx*dx) + (double)(dy*dy));
				int[] X = new int[(int)Math.round(nowRadius)];
				int[] Y = new int[(int)Math.round(nowRadius)];
				double[] grayAtThisAngle = new double[(int)Math.round(nowRadius)];
				//double[] splinedGrayAtThisAngle = new double[standardRadius];
				int[] radius = new int[(int)Math.round(nowRadius)];
				
				int distanceCounter = 0;
				for (double checkRadius = 0 ; checkRadius < (int)Math.round(nowRadius); checkRadius++) {					
					X[distanceCounter] = (int)Math.round(Math.sin(myAngles[alpha]) * checkRadius + jCO.xmid[i]);
					Y[distanceCounter] = (int)Math.round(Math.cos(myAngles[alpha]) * checkRadius + jCO.ymid[i]);				
					grayAtThisAngle[distanceCounter] = (int)Math.round(nowIP.getPixelValue(X[distanceCounter], Y[distanceCounter]));	
					radius[distanceCounter]	= (int)Math.round(checkRadius);
					distanceCounter++;					
				}	
				
				//stretch gray values to standardized radius
				float[] newRadius = new float[distanceCounter];
				for (int j = 0 ; j < distanceCounter ; j++) newRadius[j] = j * standardRadius / distanceCounter;
				SplineFitter sF = new SplineFitter(newRadius, rC.castDouble2Float(grayAtThisAngle), radius.length);
				for (int j = 0 ; j < standardRadius ; j++) {
					double newgray = sF.evalSpline(j);
					//sRadius[j] = j;
					//splinedGrayAtThisAngle[j] = newgray;
					radialgrayValues[i][alpha][j] = newgray;
				}
				
				//display radial gray value profile
				//disp.plotXYXY(rC.castInt2Double(radius), grayAtThisAngle, sRadius, splinedGrayAtThisAngle, "radius", "GV", "splined GV");
				 				
			}
			
		}
		
		return radialgrayValues;
	}
	
	public ImagePlus createGradientImage(InputOutput.MyFileCollection mFC, ImagePlus myTiff, boolean saveImg) {
		
		InputOutput jIO = new InputOutput();
		
		ImageStack nowStack = myTiff.getStack().duplicate();
		Strel3D se = CubeStrel.fromDiameter(2);
		ImageStack grad = Morphology.gradient(nowStack, se);
		
		ImagePlus gradTiff = new ImagePlus("Gradient3D", grad);
	
		if (saveImg) jIO.tiffSaver(mFC.myGradientFolder, mFC.fileName, gradTiff);
		
		return gradTiff;		
	}
	
	public ImagePlus applyClosing2Image(InputOutput.MyFileCollection mFC, ImagePlus myTiff, boolean saveImg) {
		
		InputOutput jIO = new InputOutput();
		
		ImageStack nowStack = myTiff.getStack().duplicate();
		Strel3D se = CubeStrel.fromDiameter(4);
		ImageStack closed = Morphology.closing(nowStack, se);
		
		ImagePlus closedTiff = new ImagePlus("", closed);
	
		if (saveImg) jIO.tiffSaver(mFC.myGradientFolder, mFC.fileName, closedTiff);
		
		return closedTiff;		
	}
	
	public SampleBrightness getSampleBrightness(ImagePlus nowTiff, ObjectDetector.ColCoords3D jCO) {
		
		//init units
		SampleBrightness myBrightness = new SampleBrightness();
		RoiHandler roi = new RoiHandler();
		HistogramStuff hist = new HistogramStuff();
		
		//init variables
		int[][] outHist = new int[jCO.heightOfColumn][65536];
		int[][] centralHist = new int[jCO.heightOfColumn][65536];
		int[][] referenceHist = new int[jCO.heightOfColumn][65536];
		int[] q01 = new int[jCO.heightOfColumn];
		int[] q50 = new int[jCO.heightOfColumn];
		int[] q60 = new int[jCO.heightOfColumn];
		int[] q80 = new int[jCO.heightOfColumn];
		int[] q95 = new int[jCO.heightOfColumn];
		int[] c01 = new int[jCO.heightOfColumn];
		int[] c50 = new int[jCO.heightOfColumn];
		int[] c60 = new int[jCO.heightOfColumn];
		int[] c80 = new int[jCO.heightOfColumn];
		int[] c95 = new int[jCO.heightOfColumn];
		int[] r01 = new int[jCO.heightOfColumn];
		int[] r50 = new int[jCO.heightOfColumn];
		int[] r60 = new int[jCO.heightOfColumn];
		int[] r80 = new int[jCO.heightOfColumn];
		int[] r95 = new int[jCO.heightOfColumn];
		
		
		//load polygon roi of inner perimeter		
		PolygonRoi[] innerRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 2);
		PolygonRoi[] centralRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, (int)Math.round(StatUtils.percentile(jCO.innerMinorRadius,50) - 100));
		PolygonRoi[] referenceRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, (int)Math.round(StatUtils.percentile(jCO.innerMinorRadius,50) - 200));
		myBrightness.wall = findMedianWallgrayValues(nowTiff, jCO);
				
		for (int i = 0; i < jCO.heightOfColumn ; i++) {
			
			IJ.showStatus("Sampling 1-D brightness profiles ... reached slice " + (i + 1) + "/" + nowTiff.getNSlices());
			
			//set image to next slice
			nowTiff.setPosition(i+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			//get median gray value of slice
			nowIP.setRoi(innerRoi[i]);
			int[] myHist = nowIP.getHistogram();
			for (int j = 0 ; j < myHist.length ; j++) outHist[i][j] = myHist[j];			
			q01[i] = hist.findQuantileFromHistogram(myHist, 0.01);	
			q50[i] = hist.findQuantileFromHistogram(myHist, 0.50);			
			q80[i] = hist.findQuantileFromHistogram(myHist, 0.80);
			nowIP.resetRoi();
			
			//get median from small inner circle
			nowIP.setRoi(centralRoi[i]);
			myHist = nowIP.getHistogram();
			for (int j = 0 ; j < myHist.length ; j++) centralHist[i][j] = myHist[j];			
			c01[i] = hist.findQuantileFromHistogram(myHist, 0.01);	
			c50[i] = hist.findQuantileFromHistogram(myHist, 0.50);			
			c80[i] = hist.findQuantileFromHistogram(myHist, 0.80);
			nowIP.resetRoi();
			
			//get reference illumination
			ImageProcessor copyIP = nowIP.duplicate();
			copyIP.setRoi(centralRoi[i]);
			copyIP.setValue(0);			
			copyIP.fill(centralRoi[i]);
			copyIP.setRoi(referenceRoi[i]);
			copyIP.setValue(0);	
			copyIP.fillOutside(referenceRoi[i]);
			copyIP.resetRoi();
			
			myHist = copyIP.getHistogram();
			myHist[0] = 0;
			for (int j = 0 ; j < myHist.length ; j++) referenceHist[i][j] = myHist[j];			
			r01[i] = hist.findQuantileFromHistogram(myHist, 0.01);	
			r50[i] = hist.findQuantileFromHistogram(myHist, 0.50);			
			r80[i] = hist.findQuantileFromHistogram(myHist, 0.80);

		}
		
		myBrightness.histogram = outHist;
		myBrightness.centralHistogram = centralHist;
		myBrightness.referenceHistogram = referenceHist;
		myBrightness.q01 = q01;
		myBrightness.q50 = q50;		
		myBrightness.q80 = q80;		
		myBrightness.central_q01 = c01;
		myBrightness.central_q50 = c50;		
		myBrightness.central_q80 = c80;
		myBrightness.ref_q01 = c01;
		myBrightness.ref_q50 = c50;		
		myBrightness.ref_q80 = c80;		
		
		return myBrightness;
		
	}
	

	
	public ApproximateColumnIllumination probeApproximateColumnIllumination(ColumnContainer colCon) {
		
		Median jMed = new Median();
		RoiHandler roi = new RoiHandler();
		FitStuff fit = new FitStuff(); 
		HistogramStuff hist = new HistogramStuff();
		ApproximateColumnIllumination medACI = new ApproximateColumnIllumination();
		DisplayThings disp = new DisplayThings();
		
		//define in which part of the image to search for the wall
		int imageHeight = colCon.nowTiff.getNSlices();
		int startCheck = (int)Math.round((double)imageHeight/5);
		int stopCheck = (int)Math.round((double)imageHeight/5*4);
		ApproximateColumnIllumination[] myACI = new ApproximateColumnIllumination[stopCheck - startCheck];
		
		//calculate median column positions
		double oXmid = jMed.evaluate(colCon.prelimCC.xmid);
		double oYmid = jMed.evaluate(colCon.prelimCC.ymid);	
		double oMajR = jMed.evaluate(colCon.prelimCC.outerMajorRadius);
		double oMinR = jMed.evaluate(colCon.prelimCC.outerMinorRadius);
		double otheta = jMed.evaluate(colCon.prelimCC.theta);		
		FitStuff.FittedEllipse fE = fit.new FittedEllipse();
				
		//create rois
		int safetyNet = 10;
		float imageDims = (colCon.nowTiff.getWidth() + colCon.nowTiff.getHeight()) / 2;
		float thickestWallOnEarth = imageDims / 12;
		float thinnestWallOnEarth = imageDims / 30;
		
		//init gray levels
		double[] inside = new double[stopCheck - startCheck];
		double[] wall = new double[stopCheck - startCheck];
		double[] outside = new double[stopCheck - startCheck];
		double[] globalContrast = new double[stopCheck - startCheck];
		
		//mercury
		fE.xCenter = oXmid; fE.yCenter = oYmid; fE.majorRadius = oMajR - thickestWallOnEarth; fE.minorRadius = oMinR - thickestWallOnEarth; fE.theta = otheta;
		PolygonRoi mercury = roi.makeRoiFromFittedEllipse(fE);
		
		//venus
		fE.xCenter = oXmid; fE.yCenter = oYmid; fE.majorRadius = oMajR - thinnestWallOnEarth; fE.minorRadius = oMinR - thinnestWallOnEarth; fE.theta = otheta;
		PolygonRoi venus = roi.makeRoiFromFittedEllipse(fE);
		
		//earth
		fE.xCenter = oXmid; fE.yCenter = oYmid; fE.majorRadius = oMajR - safetyNet; fE.minorRadius = oMinR - safetyNet; fE.theta = otheta;
		PolygonRoi earth = roi.makeRoiFromFittedEllipse(fE);

		//mars
		fE.xCenter = oXmid; fE.yCenter = oYmid; fE.majorRadius = oMajR + safetyNet; fE.minorRadius = oMinR + safetyNet; fE.theta = otheta;
		PolygonRoi mars = roi.makeRoiFromFittedEllipse(fE);
		
		//universe
		Roi all = new Roi(1, 1, colCon.nowTiff.getWidth(), colCon.nowTiff.getHeight()); 
		
		//probe the illuminations..
		for (int i = startCheck ; i < stopCheck ; i++) {
			
			colCon.nowTiff.setPosition(i+1);
			ImageProcessor nowIP = colCon.nowTiff.getProcessor();
			nowIP.setValue(0);
			nowIP.setBackgroundValue(0);
			
			//get soil illumination
			ImageProcessor soilIP = nowIP.duplicate();			
			soilIP.fillOutside(mercury);			
			inside[i-startCheck] = hist.findMedianFromHistogram(soilIP.getHistogram());
			//disp.displayIP(soilIP, "soil");			
			
			//get wall illumination
			ImageProcessor wallIP = nowIP.duplicate();			
			wallIP.fillOutside(earth);
			wallIP.fill(venus);
			wall[i-startCheck] = hist.findMedianFromHistogram(wallIP.getHistogram());
			//disp.displayIP(wallIP, "wall");			
			
			//get outside illumination
			ImageProcessor outsiIP = nowIP.duplicate();			
			outsiIP.fill(mars);	
			outsiIP.setRoi(all);
			outside[i-startCheck] = hist.findMedianFromHistogram(outsiIP.getHistogram());			
			//disp.displayIP(outsiIP, "outside");
			
			//getGlobalContrast			
			IlluminationInfo jII = probeIlluminationLevel(nowIP, median3DCoords(colCon.prelimCC), false);
			globalContrast[i-startCheck] = jII.quantile90 - jII.quantile10;
			
		}
		
		//calculate medians		
		medACI.inside = (int)Math.round(jMed.evaluate(inside));
		medACI.wall = (int)Math.round(jMed.evaluate(wall));
		medACI.outside = (int)Math.round(jMed.evaluate(outside));
		medACI.globalContrast = (int)Math.round(jMed.evaluate(globalContrast));
		
		return medACI;
	}
	
	public HistogramStuff.IlluminationInfo probeIlluminationLevel(ImageProcessor nowIP, ColCoords2D prelimCC, boolean cutAwaySoil) {
		
		//init objects
		HistogramStuff hist = new HistogramStuff();
		RoiHandler roi = new RoiHandler();
		FitStuff fit = new FitStuff();
		
		//init varis
		HistogramStuff.IlluminationInfo jII = hist.new IlluminationInfo();
		ImageProcessor iIP = nowIP.duplicate();
		
		//prepare roi on the wall .. make a roi that is approximately on the center of the column wall..
		if (cutAwaySoil == true) {	//if prelimCC is not empty
			double ratioBetweenInnerAndOuterRadius = 0.97; 		//may need adaptation for very thin columns..
			FitStuff.FittedEllipse fE = fit.new FittedEllipse();		
			
			//if there has been an inner diameter found, use the inner diameter (because it is not affected by a bevel).. otherwise, use an outer diameter
			if (prelimCC.ixCenter > 0) {
				fE.xCenter = prelimCC.ixCenter;
				fE.yCenter = prelimCC.iyCenter;
				fE.majorRadius = 1.005 * prelimCC.innerMajorRadius;
				fE.minorRadius = 1.005 * prelimCC.innerMinorRadius;
				fE.theta = prelimCC.itheta;
			}
			else {
				fE.xCenter = prelimCC.xCenter;
				fE.yCenter = prelimCC.yCenter;
				fE.majorRadius = ratioBetweenInnerAndOuterRadius * prelimCC.outerMajorRadius;
				fE.minorRadius = ratioBetweenInnerAndOuterRadius * prelimCC.outerMinorRadius;
				fE.theta = prelimCC.theta;
			}
					
			PolygonRoi wallRoi = roi.makeRoiFromFittedEllipse(fE);
				
			//cut away soil to get gray values of wall and surrounding, only...			
			iIP.setRoi(wallRoi);
			iIP.setColor(Color.black);
			iIP.fill(wallRoi);
			iIP.resetRoi();
			
			//DisplayThings dT = new DisplayThings();
			//dT.displayIP(iIP, "");			
		}	
		
		//get histogram of image
		int[] myHist = iIP.getHistogram();
		myHist[0] = 0;		//set 0 entries to zero
				
		//assign values	
		jII.min = hist.findQuantileFromHistogram(myHist, 0.00001);
		jII.quantile01 = hist.findQuantileFromHistogram(myHist, 0.01);
		jII.quantile05 = hist.findQuantileFromHistogram(myHist, 0.05);
		jII.quantile10 = hist.findQuantileFromHistogram(myHist, 0.1);
		jII.lowerQuartile = hist.findQuantileFromHistogram(myHist, 0.25);
		jII.median = hist.findQuantileFromHistogram(myHist, 0.5);
		jII.upperQuartile = hist.findQuantileFromHistogram(myHist, 0.75);
		jII.quantile90 = hist.findQuantileFromHistogram(myHist, 0.9);
		jII.quantile95 = hist.findQuantileFromHistogram(myHist, 0.95);
		jII.quantile99 = hist.findQuantileFromHistogram(myHist, 0.99);
		jII.max = hist.findQuantileFromHistogram(myHist, 0.999999);
		
		jII.mean = (int)Math.round(hist.findMeanFromHistogram(myHist));
		
		return jII;		
	}
	
	public ColCoords3D findOrientationOfPVCOrAluColumn(ImagePlus nowTiff, MenuWaiter.ColumnFinderMenuReturn jCFS, boolean isPrecise) {

		ColCoords3D pCC = new ColCoords3D();
		FitStuff fit = new FitStuff();
		TailoredMaths math = new TailoredMaths();
		
		//init outline related variables
		int i, j, k;
		int maxAlpha = 360;
		int dAlpha = 5;		
		int checkRange = (nowTiff.getWidth()/2 + nowTiff.getHeight()/2) ;		
		double xmid = nowTiff.getWidth() / 2;
		double ymid = nowTiff.getHeight() / 2;
		double radius = nowTiff.getHeight() / 2;
		int footprintOfMedianFilter = jCFS.medianFilter;
		if (footprintOfMedianFilter == 0) footprintOfMedianFilter = 5;
		int iCounter = 0;
		
		//find slices corresponding to the 30 and 70 percentiles and 
		double imageHeight = nowTiff.getNSlices();
		int startSlice = (int)(imageHeight * 0.25) + 1;
		int stopSlice = (int)(imageHeight * 0.75) + 1;
		int samplingNumber = 30;
		double dSlice = (stopSlice - startSlice) / samplingNumber;
		
		//change this selection if the precise column coordinates are sought
		if (isPrecise) {			
			
			startSlice = 1;
			samplingNumber = nowTiff.getNSlices();
			dSlice = 1;					
		}
		
		//init outline related vectors
		double[] x = new double[checkRange];
		double[] y = new double[checkRange];	
		double[] grayAtThisAngle = new double[checkRange];
		double[] angleAtThisAngle = new double[maxAlpha/dAlpha];
		float[] xOD = new float[maxAlpha/dAlpha]; // x of outer diameter
		float[] yOD = new float[maxAlpha/dAlpha]; // y of outer diameter
				
		//init out vectors
		double confidenceLevel = 0.99;
		double[] xCenter = new double[samplingNumber];
		double[] yCenter = new double[samplingNumber];
		double[] zCenter = new double[samplingNumber];
		double[] majorRadius = new double[samplingNumber];
		double[] minorRadius = new double[samplingNumber];
		double[] xyAngle = new double[samplingNumber];
		double[] R2 = new double[samplingNumber];
		double[] tilt = new double[3];
									
		for (i = startSlice ; i < startSlice + samplingNumber * dSlice ; i += dSlice) {  
			
			IJ.showStatus("Searching for location and orientation of column..");
							
			int cc=0;
			nowTiff.setPosition(i);				
			ImageProcessor myIP = nowTiff.getProcessor();
			
			//probe illumination level of background
			ColumnContainer colCon = new ColumnContainer();
			colCon.jCFS = jCFS;
			HistogramStuff.IlluminationInfo jII = probeIlluminationLevel(myIP, null, false);
			int myContrast = jII.quantile90 - jII.quantile10;
			int minAirPVCGradient = (int)Math.round(jCFS.airWallContrast * (double)myContrast);
				
			//sample gray values around the column
			for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {				
				
				int cc0=1;
				angleAtThisAngle[cc] = angle;
				for (double checkRadius = radius - 1; checkRadius > radius - checkRange ; checkRadius--) {					
					x[cc0] = Math.sin(angle) * checkRadius + xmid;
					y[cc0] = Math.cos(angle) * checkRadius + ymid;				
					grayAtThisAngle[cc0] = myIP.getPixelValue((int)x[cc0], (int)y[cc0]);					
					cc0++;
				}
				grayAtThisAngle[0] = 65000; 	//set first entry to a very large value so that there is no gradient (it would be 0 otherwise) 
				
				//apply medianFilter
				double[] mgrayAtThisAngle = math.oneDMedianFilter(grayAtThisAngle, footprintOfMedianFilter);
				
				//calculate the first derivative, dMG
				double[] dMG = math.oneDFirstDerivative(mgrayAtThisAngle, footprintOfMedianFilter);
				
				if (jCFS.isAlreadyNormalized) {
					for (j = 0 ; j < checkRange - 3; j++) {
						if (mgrayAtThisAngle[j] > jCFS.fixedWallGrayValue - jCFS.stdFixedWallGrayValue) {
							xOD[cc] = (float)x[j];
							yOD[cc] = (float)y[j];
							break;
						}
					}
				}
				else {
					//look for the outer edge
					for (j = 0 ; j < checkRange - 3; j++) {
						if (dMG[j] > minAirPVCGradient) {
							int myJ = j;
							for (k = j ; k < j + 3 ; k++) {
								if (dMG[k] > dMG[myJ]) myJ = k;
							}		
							xOD[cc] = (float)x[myJ];
							yOD[cc] = (float)y[myJ];
							break;
						}
					}
				}
				
				//increase vector index
				cc++;
			}
			
			//kick out values where the wall was not found
			int zeroCounter = 0;
			cc = 0;
			for (j = 0 ; j < xOD.length ; j++) if (xOD[j] == 0 & yOD[j] == 0) zeroCounter++;
			float[] xOR = new float[xOD.length - zeroCounter]; // x of outer diameter and zeros removed
			float[] yOR = new float[yOD.length - zeroCounter]; // y of outer diameter and zeros removed
			double[] angleR = new double[yOD.length - zeroCounter];  //also make an angle that has equally many entries..
			for (j = 0 ; j < xOD.length ; j++) if (xOD[j] == 0 & yOD[j] == 0) {} // do nothing
			else {
				xOR[cc] = xOD[j];
				yOR[cc] = yOD[j];
				angleR[cc] = angleAtThisAngle[j];
				cc++;
			}		
			
			//fit an elliptic ROI to the x and y
			EllipseFitter jEF = new EllipseFitter();			
			PolygonRoi pRoi = new PolygonRoi(xOR, yOR, Roi.POLYLINE); // create pRoi with outer Wall coordinates
			myIP.setRoi(pRoi); 
			jEF.fit(myIP, myIP.getStatistics());	
			
			//also kick out outliers 
			FitStuff.GoodnessOfFit gOF = fit.calculateR2OfEllipsoidFit(angleR, xOR, yOR, jEF);
			double maxDeviationAllowed = 0.05;
			ArrayList<Integer> badOnes = new ArrayList<Integer>();
			for (j = 0 ; j < xOR.length ; j++) if (Math.abs((gOF.d2m[j] - gOF.median2d2m) / gOF.median2d2m) > maxDeviationAllowed) badOnes.add(j);
			float[] xOS = new float[xOR.length - badOnes.size()]; // x of outer diameter and zeros and outliers removed
			float[] yOS = new float[yOR.length - badOnes.size()]; // y of outer diameter and zeros and outliers removed
			double[] angleS = new double[yOR.length - badOnes.size()];  //also make an angle that has equally many entries..
			cc = 0;
			for (j = 0 ; j < xOR.length ; j++) if (Math.abs((gOF.d2m[j] - gOF.median2d2m) / gOF.median2d2m) < maxDeviationAllowed) {
				xOS[cc] = xOR[j];
				yOS[cc] = yOR[j];
				angleS[cc] = angleR[j];
				cc++;
			}
			
			//do the fit again but without outliers
			PolygonRoi sRoi = new PolygonRoi(xOS, yOS, Roi.POLYLINE); // create pRoi with outer Wall coordinates
			myIP.setRoi(sRoi); 
			EllipseFitter jEFS = new EllipseFitter();
			jEFS.fit(myIP, myIP.getStatistics());			
			gOF = fit.calculateR2OfEllipsoidFit(angleS, xOS, yOS, jEFS);
			
			//transfer fitting results to vectors for export
			xCenter[iCounter] = jEFS.xCenter;
			yCenter[iCounter] = jEFS.yCenter;
			zCenter[iCounter] = i;
			minorRadius[iCounter] = jEFS.minor / 2;
			majorRadius[iCounter] = jEFS.major / 2;				
			xyAngle[iCounter] = jEFS.theta;			
			R2[iCounter] = gOF.R2;
			
			iCounter++;
		}
		
		//kick out fits that did not work
		int cc = 0;
		for (i = 0 ; i < R2.length ; i++) if (R2[i] > confidenceLevel) cc++;
		double[] xFiltered = new double[cc];
		double[] yFiltered = new double[cc];
		double[] zFiltered = new double[cc];
		double[] oMajRF = new double[cc];
		double[] oMinRF = new double[cc];
		double[] thetaF = new double[cc];
		double[] R2F = new double[cc];
		
		cc = 0;
		for (i = 0 ; i < R2.length ; i++) {
			if (R2[i] > confidenceLevel) {
				xFiltered[cc]=xCenter[i];
				yFiltered[cc]=yCenter[i];
				zFiltered[cc]=zCenter[i];
				oMajRF[cc]=majorRadius[i];
				oMinRF[cc]=minorRadius[i];
				thetaF[cc]=xyAngle[i];
				R2F[cc]=R2[i];
				cc++;
			}
		}
		
		//find column Tilt!
		tilt = findColumnTilt(xFiltered, yFiltered, zFiltered);
		
		//save results in struct
		pCC.xmid = xFiltered;
		pCC.ymid = yFiltered;
		pCC.zmid = zFiltered;
		pCC.outerMinorRadius = oMinRF;
		pCC.outerMajorRadius = oMajRF;
		pCC.theta = thetaF;
		pCC.outerR2 = R2F;
		pCC.tiltInXZ = tilt[0];
		pCC.tiltInYZ = tilt[1];
		pCC.tiltTotal = tilt[2];
		
		return pCC;

	}
		
	public ColCoords3D findOrientationSteelColumn(ImagePlus nowTiff) {

		//init objects
		ColCoords3D pCC = new ColCoords3D();
		HistogramStuff hist = new HistogramStuff();
		TailoredMaths math = new TailoredMaths();
		FitStuff fit = new FitStuff();
		
		//init outline related variables		
		int i, j;
		int maxAlpha = 360;
		int dAlpha = 5;		
		int checkRange = nowTiff.getWidth() / 4;
		double xmid = nowTiff.getWidth() / 2;
		double ymid = nowTiff.getHeight() / 2;
		double radius = nowTiff.getHeight() / 2;				
		int footprintOfMedianFilter = 5;
		int iCounter = 0;
		int gradientSearchWindow = 5;		
		double steelgrayThreshold;
		
		//find slices corresponding to the 30 and 70 percentiles and 
		double imageHeight = nowTiff.getNSlices();
		int startSlice = (int)(imageHeight * 0.3) + 1;
		int stopSlice = (int)(imageHeight * 0.7) + 1;
		int samplingNumber = 30;
		int dSlice = (stopSlice - startSlice) / samplingNumber;		
		
		//init outline related vectors		
		double[] x = new double[checkRange];
		double[] y = new double[checkRange];	
		double[] grayAtThisAngle = new double[checkRange];
		double[] angleAtThisAngle = new double[maxAlpha/dAlpha];
		float[] xOD = new float[maxAlpha/dAlpha]; // x of outer diameter
		float[] yOD = new float[maxAlpha/dAlpha]; // y of outer diameter
		
		//init out vectors
		double[] xCenter = new double[samplingNumber];
		double[] yCenter = new double[samplingNumber];
		double[] zCenter = new double[samplingNumber];
		double[] majorRadius = new double[samplingNumber];
		double[] minorRadius = new double[samplingNumber];
		double[] medianOuterRadius = new double[samplingNumber];				
		double[] xyAngle = new double[samplingNumber];
		double[] R2 = new double[samplingNumber];
		double[] tilt = new double[3];	
	
								
		//loop over the some of the slices in the image..
		for (i = startSlice ; i < startSlice + samplingNumber * dSlice ; i = i + dSlice) {  
			
			IJ.showStatus("Searching for location and orientation of column..");

			int cc=0;
			nowTiff.setPosition(i);				
			ImageProcessor myIP = nowTiff.getProcessor();
			
			//determine steelgray Threshold
			int[] myHist = myIP.getHistogram();
			double[] cumHist = hist.calcCumulativeHistogram(myHist);
			steelgrayThreshold = hist.findPercentileFromCumHist(cumHist, 0.95);
			
			//nowTiff.draw();
			//nowTiff.show();
			
			//sample gray values around the column
			for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {				
				
				int cc0=0;
				angleAtThisAngle[cc] = angle;
				for (double checkRadius = radius ; checkRadius > radius - checkRange ; checkRadius--) {					
					x[cc0] = Math.sin(angle) * checkRadius + xmid;
					y[cc0] = Math.cos(angle) * checkRadius + ymid;				
					grayAtThisAngle[cc0] = myIP.getPixelValue((int)x[cc0], (int)y[cc0]);					
					cc0++;
				}
				
				//apply medianFilter
				double[] mgrayAtThisAngle = math.oneDMedianFilter(grayAtThisAngle, footprintOfMedianFilter);
							
				//calculate the first derivative, dMG
				double[] dMG = math.oneDFirstDerivative(mgrayAtThisAngle, footprintOfMedianFilter);
								
				//find first steel
				int myWallPosition = 0;
				for (j = 0 ; j < mgrayAtThisAngle.length ; j++) if (mgrayAtThisAngle[j] > steelgrayThreshold) {
					myWallPosition = j;
					break;
				}
				
				//if first steel was not found..
				if (myWallPosition == 0) for (j = 0 ; j < mgrayAtThisAngle.length ; j++) if (mgrayAtThisAngle[j] > 0.95 * steelgrayThreshold) {
					myWallPosition = j;
					break;
				}
								
				//look for the outer edge
				double maxGradient = -50000;
				int maxGradPosition = 0;
				if (myWallPosition > gradientSearchWindow & myWallPosition < dMG.length - gradientSearchWindow) {
					for (j = myWallPosition - gradientSearchWindow ; j < myWallPosition + gradientSearchWindow; j++) {
						if (dMG[j] > maxGradient) {	
							maxGradient = dMG[j];
							maxGradPosition = j;
						}
					}
				} else {
					maxGradPosition = 0;
				}
				
				xOD[cc] = (float)x[maxGradPosition];
				yOD[cc] = (float)y[maxGradPosition];
				
				//increase vector index
				cc++;
			}
			
			//fit an elliptic ROI to the x and y
			EllipseFitter jEF = new EllipseFitter();			
			PolygonRoi pRoi = new PolygonRoi(xOD, yOD, Roi.POLYLINE); // create pRoi with outer Wall coordinates
			myIP.setRoi(pRoi); 
			jEF.fit(myIP, myIP.getStatistics());	
			
			//transfer fitting results to vectors for export
			xCenter[iCounter] = jEF.xCenter;
			yCenter[iCounter] = jEF.yCenter;
			zCenter[iCounter] = i;
			minorRadius[iCounter] = jEF.minor / 2;
			majorRadius[iCounter] = jEF.major / 2;
			medianOuterRadius[iCounter] = (minorRadius[iCounter] + majorRadius[iCounter]) / 2;
			xyAngle[iCounter] = jEF.theta;
			FitStuff.GoodnessOfFit gOF = fit.calculateR2OfEllipsoidFit(angleAtThisAngle, xOD, yOD, jEF);
			R2[iCounter] = gOF.R2;
			
			iCounter++;		
 
		}
		
		//find column Tilt!
		tilt = findColumnTilt(xCenter, yCenter, zCenter);
		
		//save results in struct
		pCC.xmid = xCenter;
		pCC.ymid = yCenter;
		pCC.zmid = zCenter;
		pCC.tiltInXZ = tilt[0];
		pCC.tiltInYZ = tilt[1];
		pCC.tiltTotal = tilt[2];
		
		return pCC;		
		
	}
	
	public EggShapedColCoords3D findEggShapedWalls(ImagePlus nowTiff, ColCoords3D prelimCC, double wallThickness) {

		//init objects
		Median jMed = new Median();
		EggShapedColCoords3D preciseCC = new EggShapedColCoords3D();
		TailoredMaths math = new TailoredMaths();
		
		//init outline related variables
		int i, j;
		int maxAlpha = 360;
		int dAlpha = 5;
		double xmid = nowTiff.getWidth() / 2;
		double ymid = nowTiff.getHeight() / 2;
		double[] imageEdges = {nowTiff.getHeight(), nowTiff.getWidth()};
		double radius = StatUtils.max(imageEdges) / 2;
		int checkRange = (int)Math.round(radius / 3) ;		
		int footprintOfMedianFilter = 3;
		int imageHeight = nowTiff.getNSlices();
		double steelgrayThreshold = 35000;
		double appliedSteelgrayThreshold = steelgrayThreshold;
		int maximumSearchWindow = 5;
		
		//create angle vector
		double[] myAngle = new double[maxAlpha/dAlpha];
		int cc = 0;	
		for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) 
			{myAngle[cc] = angle; cc++;}
		
		//set starting point for outerColumn search
		double[] x = new double[checkRange];
		double[] y = new double[checkRange];	
		double[] grayAtThisAngle = new double[checkRange];
		double[] angleAtThisAngle = new double[maxAlpha/dAlpha];
		float[] xOD = new float[maxAlpha/dAlpha + 1]; // x of outer diameter
		float[] yOD = new float[maxAlpha/dAlpha + 1]; // y of outer diameter
		float[] xID = new float[maxAlpha/dAlpha + 1]; // x of inner diameter
		float[] yID = new float[maxAlpha/dAlpha + 1]; // y of inner diameter
		double[][] aXOD = new double[nowTiff.getNSlices()][maxAlpha/dAlpha + 1]; // x of outer diameter
		double[][] aYOD = new double[nowTiff.getNSlices()][maxAlpha/dAlpha + 1]; // y of outer diameter
		double[][] aXID = new double[nowTiff.getNSlices()][maxAlpha/dAlpha + 1]; // x of inner diameter
		double[][] aYID = new double[nowTiff.getNSlices()][maxAlpha/dAlpha + 1]; // y of inner diameter
		
		//init pre-out vectors
		double[] xCenter = new double[imageHeight];
		double[] yCenter = new double[imageHeight];
		double[] zCenter = new double[imageHeight];
		double[] outerDiameter = new double[imageHeight];
		double[] accurateWallThickness = new double[imageHeight];
		double[] tilt = new double[3];	
		boolean[] columnDetected = new boolean[imageHeight];
		
		//init out vectors
		PolygonRoi[] oRoi = new PolygonRoi[imageHeight]; 	
		PolygonRoi[] iRoi = new PolygonRoi[imageHeight]; 	
		
		//nowTiff.draw();
		//nowTiff.show();
								
		for (i = 0 ; i < nowTiff.getNSlices() ; i++) {  
		
			IJ.showStatus("Searching for column's wall of slice #" + (i + 1) + "/" + nowTiff.getNSlices());
							
			int coco=0;
			nowTiff.setPosition(i+1);				
			ImageProcessor myIP = nowTiff.getProcessor();
			
			//nowTiff.draw();
			//nowTiff.show();
			
			//sample gray values around the column
			for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {				
				
				int cc0=0;
				angleAtThisAngle[coco] = angle;
				for (double checkRadius = radius ; checkRadius > radius - checkRange ; checkRadius--) {					
					x[cc0] = Math.sin(angle) * checkRadius + xmid;
					y[cc0] = Math.cos(angle) * checkRadius + ymid;				
					grayAtThisAngle[cc0] = myIP.getPixelValue((int)x[cc0], (int)y[cc0]);					
					cc0++;
				}
				
				//apply medianFilter
				double[] mGrayAtThisAngle = math.oneDMedianFilter(grayAtThisAngle, footprintOfMedianFilter);
				
				//calculate the first derivative,dMG 
				//double[] dMG = math.oneDFirstDerivative(mgrayAtThisAngle, footprintOfMedianFilter);
			  
				//find first steel 
				int myWallPosition = 0; 
				appliedSteelgrayThreshold = steelgrayThreshold; 
				for (j = 3 ; j < mGrayAtThisAngle.length ; j++) 
				if (mGrayAtThisAngle[j] > appliedSteelgrayThreshold) { 
					myWallPosition = j;
					break; 
				}		  
				
				//if wall is not found set coordinates to 0;
				if (myWallPosition == 0) {
					xOD[coco] = 0f;
					yOD[coco] = 0f;									
					xID[coco] = 0f;
					yID[coco] = 0f;
				}					
				else {
								
					//find maximum gray value
					int maxValue = 0;
					double maxGray = 0;
					for (j = myWallPosition - maximumSearchWindow ; j < myWallPosition + maximumSearchWindow ; j++) {
						if (j >= mGrayAtThisAngle.length || j < 0);
						else if (mGrayAtThisAngle[j] > maxGray) {
							maxValue = j;
							maxGray = mGrayAtThisAngle[j];		
						}
					}
					
					//look for the outer edge 
					//double maxGradient = -500000; 
					//int maxGradPosition = 0; 
					//if (myWallPosition > gradientSearchWindow & myWallPosition < dMG.length - gradientSearchWindow) { 
					//	for (j = myWallPosition - gradientSearchWindow ; j < myWallPosition + gradientSearchWindow; j++) { 
					//		if (dMG[j] > maxGradient) {
					//			maxGradient = dMG[j]; maxGradPosition = j; 
					//			} 
					//		} 
					//	} 
					//else	maxGradPosition = 0;					
					//}
					
					xOD[coco] = (float)x[maxValue - 1];			//minus 1 to move the detected edge one voxel further to the outside. 
					yOD[coco] = (float)y[maxValue - 1];	
					if (j + (int)Math.round(wallThickness) < x.length) {
						xID[coco] = (float)x[maxValue + (int)Math.round(wallThickness) - 1];
						yID[coco] = (float)y[maxValue + (int)Math.round(wallThickness) - 1];
					}					
				}
				//test approach
				//myIP.setValue(60000);
				//myIP.drawDot((int)xID[coco], (int)yID[coco]);
				//nowTiff.updateAndDraw();
				//nowTiff.show();
				
				//increase vector index
				coco++;
			}
			
			xOD[coco]=xOD[0];
			yOD[coco]=yOD[0];
			xID[coco]=xID[0];
			yID[coco]=yID[0];
			
			//test if ROI was found properly..			 
			/*
			 * if (zeros.size() < 5) { Overlay mO = new Overlay(); mO.add(oRoi[i]);
			 * mO.setStrokeColor(Color.red);
			 * 
			 * nowTiff.setPosition(i); ImageProcessor nowIP =
			 * nowTiff.getProcessor().convertToRGB();
			 * 
			 * nowIP.setRoi(oRoi[i]); nowIP.drawOverlay(mO);
			 * 
			 * ImagePlus testImg = new ImagePlus("", nowIP); testImg.updateAndDraw();
			 * testImg.show();
			 * 
			 * IJ.wait(100); testImg.close(); }
			 */
			
			//calculate diameters
			double[] diameters = new double[maxAlpha/dAlpha/2];
			double[] xcent = new double[maxAlpha/dAlpha/2];
			double[] ycent = new double[maxAlpha/dAlpha/2];
			for (int alpha = 0 ; alpha < maxAlpha/dAlpha/2 ; alpha++) {
				
				double dx = xOD[alpha] - xOD[alpha + maxAlpha/dAlpha/2];
				double dy = yOD[alpha] - yOD[alpha + maxAlpha/dAlpha/2];
				
				diameters[alpha] = Math.sqrt(dx*dx + dy*dy);
				
				double[] xs = {xOD[alpha], xOD[alpha+maxAlpha/dAlpha/2]};
				double[] ys = {yOD[alpha], yOD[alpha+maxAlpha/dAlpha/2]};
				xcent[alpha] = StatUtils.mean(xs);
				ycent[alpha] = StatUtils.mean(ys);
				
			}
			
			outerDiameter[i] = StatUtils.percentile(diameters,50);
			xCenter[i] = StatUtils.percentile(xcent,50);
			yCenter[i] = StatUtils.percentile(ycent,50);			
			zCenter[i] = i;
			
			//test goodness of find
			double[] radius0 = new double[maxAlpha/dAlpha];
			for (int alpha = 0 ; alpha < maxAlpha/dAlpha ; alpha++) {
				
				double dx = xOD[alpha] - xCenter[i];
				double dy = yOD[alpha] - yCenter[i];
				
				radius0[alpha] = Math.sqrt(dx*dx + dy*dy);
			}
			
			double medianRadius = StatUtils.percentile(radius0,50);
			
			
			for (int alpha = 0 ; alpha < maxAlpha/dAlpha ; alpha++) {
				
				double rdr = Math.abs(radius0[alpha] - medianRadius) / medianRadius;
				
				if (rdr > 0.05) xOD[alpha] = 0;  //flag as fishy if the deviation is too large.
				
			}
			
			//write ROI coordinates into array.. 
			for (j = 0 ; j < xOD.length ; j++) aXOD[i][j] = xOD[j]; 
			for (j = 0 ; j < yOD.length ; j++) aYOD[i][j] = yOD[j];
			for (j = 0 ; j < xID.length ; j++) aXID[i][j] = xID[j]; 
			for (j = 0 ; j < yID.length ; j++) aYID[i][j] = yID[j];
		}
		
		//find depth dependent wall thickness
		double unbevelledDiameter = StatUtils.percentile(outerDiameter,50);
		for (i = 0 ; i < outerDiameter.length ; i++) {			
			accurateWallThickness[i] = wallThickness + (outerDiameter[i] - unbevelledDiameter) / 2 ;			
		}
		
		//find inner diameter more accurate
		for (i = 0 ; i < outerDiameter.length ; i++) {
			
			double nowDiamRatio = accurateWallThickness[i] / wallThickness;
			
			if (nowDiamRatio < 0.999 | nowDiamRatio > 1.001) {
				for( j = 0 ; j < xOD.length ; j++) {
					xOD[j] = (float)aXOD[i][j];
					yOD[j] = (float)aYOD[i][j];
					xID[j] = (float)aXID[i][j];
					yID[j] = (float)aYID[i][j];
					
					double dx = xOD[j] - xID[j];
					double dy = yOD[j] - yID[j];
					
					double beta = Math.atan2(dx, dy);
					
					double newIX = Math.sin(beta) * nowDiamRatio * wallThickness;
					double newIY = Math.cos(beta) * nowDiamRatio * wallThickness;
					
					xID[j] = xOD[j] - (float)newIX;
					yID[j] = yOD[j] - (float)newIY;
					aXID[i][j] = xID[j];
					aYID[i][j] = yID[j];
				}				
			}
			
			//kick out zero values
			ArrayList<Integer> zeros = new ArrayList<Integer>();
			for (j = 0 ; j < xOD.length ; j++) if (xOD[j] == 0) zeros.add(j);
		
			//make purified versions of xOD and co.
			float[] pxOD = new float[xOD.length - zeros.size()];
			float[] pyOD = new float[xOD.length - zeros.size()];
			float[] pxID = new float[xOD.length - zeros.size()];
			float[] pyID = new float[xOD.length - zeros.size()];
				
			cc = 0;
			for (j = 0 ; j < xOD.length ; j++) {
				if (!zeros.contains(j)) {
					pxOD[cc] = xOD[j];
					pyOD[cc] = yOD[j];	
					pxID[cc] = xID[j];
					pyID[cc] = yID[j];	
					cc++;
				}
			}
			
			//fit an elliptic ROI to the x and y			
			oRoi[i] = new PolygonRoi(pxOD, pyOD, Roi.POLYLINE); // create pRoi with outer Wall coordinates
			iRoi[i] = new PolygonRoi(pxID, pyID, Roi.POLYLINE); // create pRoi with outer Wall coordinates
			if (zeros.size() < 5) {
				oRoi[i].fitSpline();
				iRoi[i].fitSpline();			
				columnDetected[i] = true;
			}
			else columnDetected[i] = false;
			
		}
		
		//find column Tilt!
		tilt = findColumnTilt(xCenter, yCenter, zCenter);
		
	    //find top and bottom of column
		int colTop = -1; 
		int colBot = -1;
		ArrayList<Integer> missingLayers = new ArrayList<Integer>(); 
		for (i = 0 ; i < columnDetected.length ; i++) {
			
			if (colTop < 0) if (columnDetected[i]) colTop = i;
			
			if (colTop >= 0) {				
				if (columnDetected[i]) colBot = i;
				else missingLayers.add(i);
			}			
		}		
		preciseCC.topOfColumn = colTop;
		preciseCC.bottomOfColumn = colBot;
		preciseCC.heightOfColumn = colBot - colTop + 1;
				
		//check for layers in which the column outlines were not found...
		for (i = 0 ; i < missingLayers.size() ; i++) {
			
			if (i > colBot) missingLayers.remove(i);
			
		}
		
		// kick out fishy XID
		ArrayList<Double> detectedDiams = new ArrayList<Double>();
		ArrayList<Double> detectedWallThicknesses = new ArrayList<Double>();
		ArrayList<Double> detectedXmid = new ArrayList<Double>();
		ArrayList<Double> detectedYmid = new ArrayList<Double>();
		for (i = 0 ; i < imageHeight ; i++) {						
			if (columnDetected[i]) {
				detectedDiams.add(outerDiameter[i]);
				detectedWallThicknesses.add(accurateWallThickness[i]);
				detectedXmid.add(xCenter[i]);
				detectedYmid.add(yCenter[i]);
			}			
		}
		double[] detectedDiameters = new double[detectedDiams.size()];
		double[] detectedXm = new double[detectedDiams.size()];
		double[] detectedYm = new double[detectedDiams.size()];
		for (i = 0 ; i < detectedDiams.size(); i++) {			
			detectedDiameters[i] = detectedDiams.get(i) - 2*detectedWallThicknesses.get(i);
			detectedXm[i] = detectedXmid.get(i);
			detectedYm[i] = detectedYmid.get(i);
		}
		double medianInnerDiameterZ = StatUtils.percentile(detectedDiameters, 50);
		double medianXm = StatUtils.percentile(detectedXm, 50);
		double medianYm = StatUtils.percentile(detectedYm, 50);
		
		for (i = 0 ; i < imageHeight ; i++) {						
			if (columnDetected[i]) {
				double minVal = 1000;
				for (j = 0 ; j < maxAlpha/dAlpha ; j++) {
					if (aXOD[i][j] < minVal) minVal = aXOD[i][j];
					if (aYOD[i][j] < minVal) minVal = aYOD[i][j];					
				}
				if (Math.abs((outerDiameter[i] - 2 * accurateWallThickness[i]) - medianInnerDiameterZ) > 0.01 * medianInnerDiameterZ ||
					Math.abs(xCenter[i] - medianXm) > 0.04 * medianXm ||	
					Math.abs(yCenter[i] - medianYm) > 0.04 * medianYm ||
					minVal < 1) {
						columnDetected[i] = false;
						missingLayers.add(i);					
				}	
			}			
		}
		
		//impute missing layers
		if (missingLayers.size() > 0) {			
			for (i = 0 ; i < missingLayers.size() ; i++) {
				
				int toImpute = missingLayers.get(i);
				int onTop = -1;
				int onBot = -1;
				
				//find next outline above
				for (j = toImpute - 1 ; j >= colTop ; j--) if (columnDetected[j]) onTop = j;
				
				// find next outline below
				for (j = toImpute + 1 ; j <= colBot ; j++) if (columnDetected[j]) onBot = j;
				
				xCenter[toImpute] = xCenter[colTop] + (toImpute - colTop) * (xCenter[colBot] - xCenter[colTop]) / (colBot - colTop);
				yCenter[toImpute] = yCenter[colTop] + (toImpute - colTop) * (yCenter[colBot] - yCenter[colTop]) / (colBot - colTop);
				accurateWallThickness[toImpute] = accurateWallThickness[colTop] + (toImpute - colTop) * (accurateWallThickness[colBot] - accurateWallThickness[colTop]) / (colBot - colTop);
				
				for (int alpha = 0 ; alpha < maxAlpha/dAlpha ; alpha++) {
					
					aXOD[toImpute][alpha] = aXOD[colTop][alpha] + (toImpute - colTop) * (aXOD[colBot][alpha] - aXOD[colTop][alpha]) / (colBot - colTop);
					aYOD[toImpute][alpha] = aYOD[colTop][alpha] + (toImpute - colTop) * (aYOD[colBot][alpha] - aYOD[colTop][alpha]) / (colBot - colTop);
					aXID[toImpute][alpha] = aXID[colTop][alpha] + (toImpute - colTop) * (aXID[colBot][alpha] - aXID[colTop][alpha]) / (colBot - colTop);
					aYID[toImpute][alpha] = aYID[colTop][alpha] + (toImpute - colTop) * (aYID[colBot][alpha] - aYID[colTop][alpha]) / (colBot - colTop);
					
				}				
			}			
		}
		
		//transfer column outlines to output struct
		preciseCC.xOD = aXOD;
		preciseCC.yOD = aYOD;
		preciseCC.xID = aXID;
		preciseCC.yID = aYID;
		
		preciseCC.medianOuterDiameter = outerDiameter;
		
		preciseCC.xmid = xCenter; 
		preciseCC.ymid = yCenter;
		preciseCC.zmid = zCenter;
		
		preciseCC.wallThickness = accurateWallThickness;
		
		preciseCC.anglesChecked = xOD.length - 1;

		preciseCC.tiltInXZ = tilt[0];
		preciseCC.tiltInYZ = tilt[1];
		preciseCC.tiltTotal = tilt[2];
					
		//and innerRadius
		double[] innerRadius = new double[preciseCC.medianOuterDiameter.length];
		for (i = 0 ; i < preciseCC.medianOuterDiameter.length ; i++) innerRadius[i] = preciseCC.medianOuterDiameter[i] / 2 - preciseCC.wallThickness[i];
		preciseCC.innerRadius = innerRadius;
		
		
		/*
		 * ImagePlus showImage = nowTiff.duplicate(); ImageConverter myIC = new
		 * ImageConverter(showImage); myIC.convertToRGB(); for(i = 1 ; i <
		 * nowTiff.getStackSize() ; i++) {
		 * 
		 * if (columnDetected[i]) { //Overlay mO = new Overlay(); //mO.add(oRoi[i]);
		 * 
		 * showImage.setPosition(i); ImageProcessor myIP = showImage.getProcessor();
		 * 
		 * myIP.setRoi(oRoi[i]); myIP.setColor(Color.red); myIP.draw(oRoi[i]);
		 * 
		 * myIP.setColor(Color.yellow); myIP.draw(iRoi[i]);
		 * 
		 * showImage.updateAndDraw(); showImage.show(); } }
		 */
		 
		
		return preciseCC;		
		
	}
	
	public ColCoords2D extract2DCoordsFrom3DCoords(ColCoords3D nowCoords, int i) {
		
		ColCoords2D outCoords = new ColCoords2D();
		
		outCoords.xCenter = nowCoords.xmid[i];
		outCoords.yCenter = nowCoords.ymid[i];
		outCoords.zCenter = nowCoords.zmid[i];
		outCoords.outerMinorRadius = nowCoords.outerMinorRadius[i];
		outCoords.outerMajorRadius = nowCoords.outerMajorRadius	[i];			
		outCoords.theta = nowCoords.theta[i];	
		outCoords.outerR2 = nowCoords.outerR2[i];
		
		outCoords.ixCenter = nowCoords.ixmid[i];
		outCoords.iyCenter = nowCoords.iymid[i];
		outCoords.innerMinorRadius = nowCoords.innerMinorRadius[i];
		outCoords.innerMajorRadius = nowCoords.innerMajorRadius[i];			
		outCoords.itheta = nowCoords.itheta[i];		
		outCoords.innerR2 = nowCoords.innerR2[i];
		
		outCoords.wallThickness = nowCoords.wallThickness[i];
				
		outCoords.columnIsAtThisDepth = nowCoords.columnIsAtThisDepth[i];
	
		outCoords.airWallContrast = nowCoords.airWallContrast[i];
		outCoords.absoluteWallgrayValue = nowCoords.absoluteWallgrayValue[i];
		
		outCoords.wallSoilStdContrastThreshold = nowCoords.wallSoilStdContrastThreshold[i];		
			
		
		return outCoords;
		
	}
	
	public ColCoords3D findColumnsTop(InputOutput.MyFileCollection mFC, ColCoords3D samCoords, MenuWaiter.ColumnFinderMenuReturn jCFS) {
		
		IJ.showStatus("Trying to find top of column ...");
				
		InputOutput jIO = new InputOutput();
		RankFilters rF = new RankFilters();
		
		//find topmost found column part in samCoords
		int soFarOnTop = 99999;
		int findCounter = 0;
		int topNumber = 0;
		for (int i = 0 ; i < samCoords.columnIsAtThisDepth.length ; i++) { // throw away the first find because it is often suboptimal
			if (samCoords.columnIsAtThisDepth[i]) {
				findCounter++;
				if (findCounter > 1) {
					soFarOnTop = (int)Math.round(samCoords.zmid[i]);
					topNumber = i;
					break;
				}
			}
		}
		
		//get 2D coords of slice that is so far on top
		ColCoords2D nowTop = extract2DCoordsFrom3DCoords(samCoords, topNumber);
		ColCoords2D formerTop = nowTop;
		
		ArrayList<Integer> toAdd = new ArrayList<Integer>();
		int validationRange = (int)Math.round(mFC.nOfSlices / 100) + 1; // check the topmost x slices
		int[] valiSeries = new int[validationRange];
		ColCoords2D[] topCoords = new ColCoords2D[validationRange];
		boolean hasBeenTried = false;
		int originallyOnTop = soFarOnTop;
		
		//define starting set of Tiffs
		int nowOnTop = soFarOnTop;
		int[] referenceFrame = {0, soFarOnTop};
		ArrayList<Integer> badTops = new ArrayList<Integer>();
		
		//check it!
		//jCFS.debug = true;
		
		int minAddSize = (int)Math.floor(mFC.nOfSlices / 800) + 1;
		if (minAddSize == 1) minAddSize = 2;
		
		while (toAdd.size() < minAddSize) {	
			
			InputOutput.SampleTiffWrapper sTW = jIO.assembleRepresentativeSample(mFC, referenceFrame);		
			
			int cc = 0;
			while (!sTW.hasConverged){  //loop until convergence
			
				cc++;
				IJ.showStatus("Searching for top of column, iteration " + cc + " ...");			
				
				nowOnTop = referenceFrame[1];
				
				//try find column walls
				for (int i = 0 ; i < sTW.samTiff.getNSlices() ; i++) {
					
					sTW.samTiff.setPosition(i+1);
					ImageProcessor myIP = sTW.samTiff.getProcessor().duplicate();
					
					//fill zero values with background
					ImageProcessor zeroIP = myIP.duplicate(); 
					zeroIP.threshold(1);
					zeroIP.dilate();
					
					ArrayList<Integer> fringeVoxels = new ArrayList<Integer>();
					for (int x = 0 ; x < zeroIP.getWidth() ; x++) {
						for (int y = 0 ; y < zeroIP.getHeight() ; y++) {
							int oldP = (int)Math.round(myIP.getPixelValue(x, y));
							int newP = (int)Math.round(zeroIP.getPixelValue(x, y));
							if (oldP > 0 & newP > 0) fringeVoxels.add((int)Math.round(myIP.getPixelValue(x, y)));					
						}
					}
					double[] fillValue = new double[fringeVoxels.size()];
					for (int j = 0 ; j < fringeVoxels.size() ; j++) fillValue[j] = fringeVoxels.get(j);
					if (!jCFS.isAlreadyNormalized) {
						int meanFillValue = (int)Math.round(StatUtils.mean(fillValue));					
						myIP.min(meanFillValue);
					}
					
					//apply 2-D median filter and go for it...
					rF.rank(myIP, jCFS.medianFilter2D, RankFilters.MEDIAN);
					ColCoords2D nowSlice = findColumnWalls2D(sTW.samSlices[i], myIP, nowTop, jCFS);	
					
					if (nowSlice.columnIsAtThisDepth) {
						nowOnTop = sTW.samSlices[i];
						nowTop = nowSlice;
						break;
					}					
				}
				
				if (nowOnTop < soFarOnTop) {
					referenceFrame[1] = nowOnTop;
					soFarOnTop = nowOnTop;				
				}
				else referenceFrame[0] += (referenceFrame[1] - referenceFrame[0]) / 2;
				
				//assemble new sample image;
				sTW = jIO.assembleRepresentativeSample(mFC, referenceFrame);						
			}
			
			//sample x images from top and check when the outline parameters stabilize..	
			int startOfValiSeries = (int)nowTop.zCenter;
			if (badTops.size() > 0) {
				while (badTops.contains(startOfValiSeries)) {
					startOfValiSeries += mFC.nOfSlices / 150;
				}
			}
			
			for (int i = startOfValiSeries ; i < startOfValiSeries + validationRange ; i++) valiSeries[i - startOfValiSeries] = i;
			
			ImagePlus valiTiff = new ImagePlus();
			if (mFC.imageHasBeenLoaded) {
				ImageStack newStack = new ImageStack(mFC.nowTiff.getWidth(), mFC.nowTiff.getHeight());
		    	for (int j = 0 ; j < valiSeries.length ; j++) {
		    		
		    		mFC.nowTiff.setPosition(valiSeries[j]);	    		
		    		ImageProcessor nowIP = mFC.nowTiff.getProcessor().duplicate();
		    		newStack.addSlice(nowIP);    		
		    		
		    	}
		    	ImagePlus outTiff = new ImagePlus("", newStack);
		    	valiTiff = outTiff;
			}			
			else valiTiff = jIO.openTiff3DSomeSlices(mFC.nowTiffPath, mFC.nowWidth, mFC.nowHeight, valiSeries);			
			
			ArrayList<Integer> z = new ArrayList<Integer>();
			
			for (int i = 0 ; i < valiTiff.getNSlices() ; i++) {
				
				if (!hasBeenTried) IJ.showStatus("Top was found.. verifying find for slice " + (valiSeries[i] + 1) + " ...");
				else IJ.showStatus("Top has been lost again.. searching it in slice " + (valiSeries[i] + 1) + " ...");
				
				valiTiff.setPosition(i + 1);
				ImageProcessor myIP = valiTiff.getProcessor().duplicate();
				
				//fill zero values with background
				ImageProcessor zeroIP = myIP.duplicate(); 
				zeroIP.threshold(1);
				zeroIP.dilate();
				
				ArrayList<Integer> fringeVoxels = new ArrayList<Integer>();
				for (int x = 0 ; x < zeroIP.getWidth() ; x++) {
					for (int y = 0 ; y < zeroIP.getHeight() ; y++) {
						int oldP = (int)Math.round(myIP.getPixelValue(x, y));
						int newP = (int)Math.round(zeroIP.getPixelValue(x, y));
						if (oldP > 0 & newP >0) fringeVoxels.add((int)Math.round(myIP.getPixelValue(x, y)));					
					}
				}
				double[] fillValue = new double[fringeVoxels.size()];
				for (int j = 0 ; j < fringeVoxels.size() ; j++) fillValue[j] = fringeVoxels.get(j);
				int meanFillValue = (int)Math.round(StatUtils.mean(fillValue));
				myIP.min(meanFillValue);
				
				//apply 2-D median filter and go for it...				
				rF.rank(myIP, jCFS.medianFilter2D, RankFilters.MEDIAN);
				topCoords[i] = findColumnWalls2D(valiSeries[i], myIP, formerTop, jCFS); 
				
				if (topCoords[i].columnIsAtThisDepth) {
					z.add(i);			
				}
			}
					
			//check col properties		
			int startCheck = 2 * validationRange / 3;
			ArrayList<Integer> zRef = new ArrayList<Integer>();
			for (int i = startCheck ; i < validationRange; i++) {
				if (z.contains(i)) {
					for (int j = 0 ; j < z.size() ; j++) {
						if (z.get(j) >= startCheck & !zRef.contains(z.get(j))) {
							zRef.add(z.get(j));					
							break;
						}
					}
				}
			}
			
			double[] wT = new double[zRef.size()];
			double[] x = new double[zRef.size()];
			double[] y = new double[zRef.size()];
			double[] wb = new double[zRef.size()];
			
			for (int i = 0 ; i < zRef.size() ; i++) { //TAKE MEDIAN OF LAST values as reference
				wT[i] = topCoords[zRef.get(i)].wallThickness;
				x[i] = topCoords[zRef.get(i)].xCenter;
				y[i] = topCoords[zRef.get(i)].yCenter;
				wb[i] = topCoords[zRef.get(i)].CVOfWallBrightness;	
			}
			
			//check start of column			
			boolean isCool = true;
			for (int i = 0 ; i < validationRange; i++) {
				
				isCool = false;
				if (z.contains(i)) isCool = true;			
				
				double wc = Math.abs(topCoords[i].wallThickness - StatUtils.percentile(wT, 50)) / StatUtils.percentile(wT, 50);
				if (wc > jCFS.percentToleratedDifference / 100) isCool = false;
				
				double xc = Math.abs(topCoords[i].xCenter - StatUtils.percentile(x, 50)) / StatUtils.percentile(x, 50);
				if (xc > jCFS.percentToleratedDifference / 100) isCool = false;
				
				double yc = Math.abs(topCoords[i].yCenter - StatUtils.percentile(y, 50)) / StatUtils.percentile(y, 50);
				if (yc > jCFS.percentToleratedDifference / 100) isCool = false;
				
				double medianWB = StatUtils.percentile(wb, 50);
				double wbc = Math.abs(topCoords[i].CVOfWallBrightness - medianWB) / StatUtils.percentile(wb, 50);
				if (wbc > 10 * jCFS.percentToleratedDifference / 100) isCool = false;  //for PVC columns this criterion needs to be relaxed
				
				if (isCool) toAdd.add(i);
			}
			
			//try again some slices below...
			if (toAdd.size() < minAddSize) {
				referenceFrame[0] = nowOnTop - validationRange + 1;
				referenceFrame[1] = originallyOnTop;
				soFarOnTop = originallyOnTop;
				badTops.add(valiSeries[0]);
				toAdd.clear();
				//jCFS.debug = true;
			}			
		}		
						
		//check or duplicates between topCoords and samCoords		
		int notAdd = 0;
		ArrayList<Integer> noAdd = new ArrayList<Integer>();
		for (int i = 0 ; i < validationRange ; i++) {
			if ((int)Math.round(samCoords.zmid[i]) <= valiSeries[valiSeries.length - 1]) {
				notAdd++;
				noAdd.add(i);	//are entries in the samCoords file that are outside the detected column top
			}
		}
		
		//create output coordinates
		ColCoords3D outCoords = initColCords3D(toAdd.size() + samCoords.zmid.length - notAdd);
		for (int i = 0 ; i < toAdd.size() ; i++) {			
			outCoords = plug2DInto3DCoords(i, topCoords[toAdd.get(i)], outCoords);			
		}		
		int[] goodSlices = new int[samCoords.zmid.length - notAdd];
		int[] placings = new int[samCoords.zmid.length - notAdd];
		int cc = 0;
		for (int i = notAdd ; i < samCoords.zmid.length ; i++) {
			goodSlices[cc] = i;
			placings[cc] = toAdd.size() + i - notAdd;			
			cc++;
		}
			
		outCoords = plug3DInto3DCoords(placings, goodSlices, samCoords, outCoords);
		outCoords.topOfColumn = (int)Math.round(topCoords[toAdd.get(0)].zCenter);
		outCoords.anglesChecked = samCoords.anglesChecked;		
		
		//create ColCords3D for output
		return outCoords;
		
	}
	
	public ColCoords3D findClosestXYSlice2Top(InputOutput.MyFileCollection mFC, ColCoords3D samCoords, MenuWaiter.ColumnFinderMenuReturn jCFS) {
		
		IJ.showStatus("Looking for column coordinates close to the top ... ...");
		
		//jCFS.debug = true;
				
		InputOutput jIO = new InputOutput();
		
		//find topmost found column part in samCoords
		int soFarOnTop = (int)samCoords.zmid[0];
		
		//find topmost found column part in samCoords
		int findCounter = 0;
		int topNumber = 0;
		for (int i = 0 ; i < samCoords.columnIsAtThisDepth.length ; i++) { // throw away the first find because it is often suboptimal
			if (samCoords.columnIsAtThisDepth[i]) {
				findCounter++;
				if (findCounter > 0) {
					soFarOnTop = (int)Math.round(samCoords.zmid[i]);
					topNumber = i;
					break;
				}
			}
		}
		
		//get 2D coords of slice that is so far on top
		ColCoords2D nowTop = extract2DCoordsFrom3DCoords(samCoords, topNumber);
		ColCoords2D formerTop = nowTop;
		
		ArrayList<Integer> toAdd = new ArrayList<Integer>();
		int validationRange = mFC.nOfSlices / 100; // check the topmost x slices
		int[] valiSeries = new int[validationRange];
		ColCoords2D[] topCoords = new ColCoords2D[validationRange];
		
		//define starting set of Tiffs
		int nowOnTop = jCFS.topOfColumn;
		ArrayList<Integer> badTops = new ArrayList<Integer>();
		
		int minAddSize = (int)Math.floor((jCFS.bottomOfColumn - jCFS.topOfColumn) / 800) + 1;
		if (minAddSize == 1) minAddSize = 2;
		
		//check it!
		while (toAdd.size() < minAddSize) {	
						
			//sample x images from top and check when the outline parameters stabilize..	
			int startOfValiSeries = nowOnTop;
			if (badTops.size() > 0) {
				while (badTops.contains(startOfValiSeries)) {
					startOfValiSeries += 1;
				}
			}
			
			for (int i = startOfValiSeries ; i < startOfValiSeries + validationRange ; i++) valiSeries[i - startOfValiSeries] = i;
			
			ImagePlus valiTiff = jIO.openTiff3DSomeSlices(mFC.nowTiffPath, mFC.nowWidth, mFC.nowHeight, valiSeries);
			
			ArrayList<Integer> z = new ArrayList<Integer>();
			
			for (int i = 0 ; i < valiTiff.getNSlices() ; i++) {
				
				valiTiff.setPosition(i + 1);
				ImageProcessor myIP = valiTiff.getProcessor().duplicate();
				topCoords[i] = findColumnWalls2D(valiSeries[i], myIP, formerTop, jCFS); 
				
				if (topCoords[i].columnIsAtThisDepth) {
					z.add(i);			
				}
			}
					
			//check col properties		
			int startCheck = 2 * validationRange / 3;
			ArrayList<Integer> zRef = new ArrayList<Integer>();
			for (int i = startCheck ; i < validationRange; i++) {
				if (z.contains(i)) {
					for (int j = 0 ; j < z.size() ; j++) {
						if (z.get(j) >= startCheck & !zRef.contains(z.get(j))) {
							zRef.add(z.get(j));					
							break;
						}
					}
				}
			}
			
			double[] wT = new double[zRef.size()];
			double[] x = new double[zRef.size()];
			double[] y = new double[zRef.size()];
			double[] r = new double[zRef.size()];
			
			for (int i = 0 ; i < zRef.size() ; i++) { //TAKE MEDAIN OF LAST values as reference
				wT[i] = topCoords[zRef.get(i)].wallThickness;
				x[i] = topCoords[zRef.get(i)].xCenter;
				y[i] = topCoords[zRef.get(i)].yCenter;
				r[i] = (topCoords[zRef.get(i)].outerMajorRadius + topCoords[i].outerMinorRadius) / 2;	
			}
			
			//check start of column			
			boolean isCool = true;
			for (int i = 0 ; i < validationRange; i++) {
				
				isCool = false;
				if (z.contains(i)) isCool = true;			
				
				double wc = Math.abs(topCoords[i].wallThickness - StatUtils.percentile(wT, 50)) / StatUtils.percentile(wT, 50);
				if (wc > jCFS.percentToleratedDifference / 100) isCool = false;
				
				double xc = Math.abs(topCoords[i].xCenter - StatUtils.percentile(x, 50)) / StatUtils.percentile(x, 50);
				if (xc > jCFS.percentToleratedDifference / 100) isCool = false;
				
				double yc = Math.abs(topCoords[i].yCenter - StatUtils.percentile(y, 50)) / StatUtils.percentile(y, 50);
				if (yc > jCFS.percentToleratedDifference / 100) isCool = false;
				
				if (isCool) toAdd.add(i);
			}
			
			//if necessary try again some slices below...
			if (toAdd.size() < minAddSize) {
				nowOnTop += 2 * minAddSize;
				badTops.add(valiSeries[0]);
				toAdd.clear();
				formerTop = topCoords[7];
				//jCFS.debug = true;
			}			
		}
						
		//check or duplicates between topCoords and samCoords		
		int notAdd = 0;
		ArrayList<Integer> noAdd = new ArrayList<Integer>();
		for (int i = 0 ; i < validationRange ; i++) {
			if ((int)Math.round(samCoords.zmid[i]) <= valiSeries[valiSeries.length - 1]) {
				notAdd++;
				noAdd.add(i);
			}
		}
		
		//squeeze in top layer
		if (topCoords[toAdd.get(0)].zCenter != jCFS.topOfColumn) topCoords[toAdd.get(0)].zCenter = jCFS.topOfColumn;
		
		//create output coordinates
		ColCoords3D outCoords = initColCords3D(toAdd.size() + samCoords.zmid.length - notAdd);
		for (int i = 0 ; i < toAdd.size() ; i++) {			
			outCoords = plug2DInto3DCoords(i, topCoords[toAdd.get(i)], outCoords);			
		}
				
		int[] goodSlices = new int[samCoords.zmid.length - notAdd];
		int[] placings = new int[samCoords.zmid.length - notAdd];
		int cc = 0;
		for (int i = notAdd ; i < samCoords.zmid.length ; i++) if (!noAdd.contains(i)){
			goodSlices[cc] = i;
			placings[cc] = toAdd.size() + i - notAdd;			
			cc++;
		}
			
		outCoords = plug3DInto3DCoords(placings, goodSlices, samCoords, outCoords);
		outCoords.topOfColumn = jCFS.topOfColumn;
		outCoords.anglesChecked = samCoords.anglesChecked;		
		
		//create ColCords3D for output
		return outCoords;
		
	}
	
	public ColCoords3D findColumnsBottom(InputOutput.MyFileCollection mFC, ColCoords3D samCoords, MenuWaiter.ColumnFinderMenuReturn jCFS) {
		
		IJ.showStatus("Trying to find bottom of column ...");
				
		//jCFS.debug = true;
		
		InputOutput jIO = new InputOutput();
		RankFilters rF=new RankFilters();
		
		//find topmost found column part in samCoords
		int soFarOnBot = 99999;
		int botNumber = 0;
		for (int i = samCoords.columnIsAtThisDepth.length - 1; i > 0 ; i--) { // throw away the first find because it is often suboptimal
			if (samCoords.columnIsAtThisDepth[i]) {
				soFarOnBot = (int)Math.round(samCoords.zmid[i]);
				botNumber = i;
				break;			
			}
		}
		
		//get 2D coords of slice that is so far on top
		ColCoords2D nowBot = extract2DCoordsFrom3DCoords(samCoords, botNumber);
		ColCoords2D outBot = nowBot;
		int botAtStart = (int)Math.round(nowBot.zCenter);
		
		ArrayList<Integer> toAdd = new ArrayList<Integer>();
		int validationRange = (int)Math.round(mFC.nOfSlices / 100) + 1; // check the bottommost x slices
		int[] valiSeries = new int[validationRange];
		ColCoords2D[] botCoords = new ColCoords2D[validationRange];
		boolean hasBeenTried = false;
		int originallyOnBot = soFarOnBot;
		int[] referenceFrame = {mFC.nOfSlices - 1, soFarOnBot};		
		int nowOnBot = soFarOnBot;
		
		ArrayList<Integer> badBots = new ArrayList<Integer>();	
		
		//jCFS.debug = true;
		//jCFS.showFit = true;
		
		int minAddSize = (int)Math.floor(mFC.nOfSlices / 800) + 1;
		if (minAddSize == 1) minAddSize = 2;
		
		while (toAdd.size() < minAddSize) {
		
			//define starting set of Tiffs			
			InputOutput.SampleTiffWrapper sTW = jIO.assembleRepresentativeSample(mFC, referenceFrame);		
				
			int cc = 0;
			while (!sTW.hasConverged){  //loop until convergence
			
				cc++;
				IJ.showStatus("Searching for bottom of column, iteration " + cc + " ...");			
				
				nowOnBot = referenceFrame[1];
				
				//try find column walls
				for (int i = 0 ; i < sTW.samTiff.getNSlices() ; i++) {
					
					sTW.samTiff.setPosition(i+1);
					ImageProcessor myIP = sTW.samTiff.getProcessor().duplicate();
					//fill zero values with background
					ImageProcessor zeroIP = myIP.duplicate(); 
					zeroIP.threshold(1);
					zeroIP.dilate();
					
					ArrayList<Integer> fringeVoxels = new ArrayList<Integer>();
					for (int x = 0 ; x < zeroIP.getWidth() ; x++) {
						for (int y = 0 ; y < zeroIP.getHeight() ; y++) {
							int oldP = (int)Math.round(myIP.getPixelValue(x, y));
							int newP = (int)Math.round(zeroIP.getPixelValue(x, y));
							if (oldP > 0 & newP >0) fringeVoxels.add((int)Math.round(myIP.getPixelValue(x, y)));					
						}
					}
					double[] fillValue = new double[fringeVoxels.size()];
					for (int j = 0 ; j < fringeVoxels.size() ; j++) fillValue[j] = fringeVoxels.get(j);
					if (!jCFS.isAlreadyNormalized) {
						int meanFillValue = (int)Math.round(StatUtils.mean(fillValue));
						myIP.min(meanFillValue);
					}
					
					//apply 2-D median filter and go for it...
					rF.rank(myIP, jCFS.medianFilter2D, RankFilters.MEDIAN);
					ColCoords2D nowSlice = findColumnWalls2D(sTW.samSlices[i], myIP, nowBot, jCFS);	
					
					if (nowSlice.columnIsAtThisDepth) {
						nowOnBot = sTW.samSlices[i];
						outBot = nowSlice;
						break;
					}					
				}
				
				if (nowOnBot > soFarOnBot) {
					referenceFrame[1] = nowOnBot;
					soFarOnBot = nowOnBot;				
				}
				else referenceFrame[0] -= (referenceFrame[0] - referenceFrame[1]) / 2;
				
				//assemble new sample image;
				sTW = jIO.assembleRepresentativeSample(mFC, referenceFrame);						
			}
			
			//sample x images from bot and check when the gray values stabilize..				
			int stopSlice = (int)outBot.zCenter - validationRange;
			if (badBots.size() > 0) {
				while (badBots.get(badBots.size()-1) - stopSlice < 5) {
					stopSlice -= mFC.nOfSlices / 150;
				}
			}
			for (int i = stopSlice + validationRange ; i > stopSlice ; i--) valiSeries[validationRange - (i - stopSlice)] = i;
			
			ImagePlus valiTiff = new ImagePlus();
			if (mFC.imageHasBeenLoaded) {
				ImageStack newStack = new ImageStack(mFC.nowTiff.getWidth(), mFC.nowTiff.getHeight());
		    	for (int j = 0 ; j < valiSeries.length ; j++) {
		    		
		    		mFC.nowTiff.setPosition(valiSeries[j]);	    		
		    		ImageProcessor nowIP = mFC.nowTiff.getProcessor().duplicate();
		    		newStack.addSlice(nowIP);    		
		    		
		    	}
		    	ImagePlus outTiff = new ImagePlus("", newStack);
		    	valiTiff = outTiff;
			}
			else valiTiff = jIO.openTiff3DSomeSlices(mFC.nowTiffPath, mFC.nowWidth, mFC.nowHeight, valiSeries);
			
			ArrayList<Integer> z = new ArrayList<Integer>();
			
			for (int i = 0 ; i < valiTiff.getNSlices() ; i++) {
				
				
				if (!hasBeenTried) IJ.showStatus("Bottom was found.. verifying find for slice " + (valiSeries[i] + 1) + " ...");
				else IJ.showStatus("Bottom has been lost again.. searching it now in slice " + (valiSeries[i] + 1) + " ...");
								
				valiTiff.setPosition(i + 1);
				ImageProcessor myIP = valiTiff.getProcessor().duplicate();
				
				//fill zero values with background
				ImageProcessor zeroIP = myIP.duplicate(); 
				zeroIP.threshold(1);
				zeroIP.dilate();
				
				ArrayList<Integer> fringeVoxels = new ArrayList<Integer>();
				for (int x = 0 ; x < zeroIP.getWidth() ; x++) {
					for (int y = 0 ; y < zeroIP.getHeight() ; y++) {
						int oldP = (int)Math.round(myIP.getPixelValue(x, y));
						int newP = (int)Math.round(zeroIP.getPixelValue(x, y));
						if (oldP > 0 & newP >0) fringeVoxels.add((int)Math.round(myIP.getPixelValue(x, y)));					
					}
				}
				double[] fillValue = new double[fringeVoxels.size()];
				for (int j = 0 ; j < fringeVoxels.size() ; j++) fillValue[j] = fringeVoxels.get(j);
				if (!jCFS.isAlreadyNormalized) {
					int meanFillValue = (int)Math.round(StatUtils.mean(fillValue));				
					myIP.min(meanFillValue);
				}
				
				//ImagePlus test = new ImagePlus("",myIP);
				//test.show();
				
				//apply 2-D median filter and go for it...
				if (jCFS.medianFilter2D > 0) rF.rank(myIP, jCFS.medianFilter2D, RankFilters.MEDIAN);
				botCoords[i] = findColumnWalls2D(valiSeries[i], myIP, outBot, jCFS); 
				
				//test.show();
				
				if (botCoords[i].columnIsAtThisDepth) {
					z.add(i);
				}
			}
					
			//check col properties
			cc = 0;
			int startCheck = 2 * validationRange / 3;
			for (int i = startCheck ; i < validationRange; i++) if (z.contains(i)) cc++;
			
			//if not enough good values were found there, try a larger range..
			if (cc < 2) {
				cc = 0;
				startCheck = 1 * validationRange / 3;
				for (int i = startCheck ; i < validationRange; i++) if (z.contains(i)) cc++;		
			}
			
			//if still not enough good values were found there, try a larger range..
			if (cc < 2) {
				cc = 0;
				startCheck = 0;
				for (int i = startCheck ; i < validationRange; i++) if (z.contains(i)) cc++;		
			}
					
			//get values in reference series
			ArrayList<Integer> zRef = new ArrayList<Integer>();
			for (int i = startCheck ; i < validationRange; i++) {
				if (z.contains(i)) {
					for (int j = 0 ; j < z.size() ; j++) {
						if (z.get(j) >= startCheck & !zRef.contains(z.get(j))) {
							zRef.add(z.get(j));					
							break;
						}
					}
				}
			}
			
			double[] x = new double[zRef.size()];
			double[] y = new double[zRef.size()];
			double[] wb = new double[zRef.size()];		
			
			for (int i = 0 ; i < zRef.size() ; i++) { //TAKE MEDIAN OF LAST values as reference
				x[i] = botCoords[zRef.get(i)].xCenter;
				y[i] = botCoords[zRef.get(i)].yCenter;
				wb[i] = botCoords[zRef.get(i)].CVOfWallBrightness;	
			}
			
			//valiTiff.updateAndDraw();
			//valiTiff.show();
			
			//check start of column
			boolean isCool = true;
			for (int i = 0 ; i < validationRange; i++) {
				
				isCool = false;		
				if (z.contains(i)) isCool = true;
					
				double xc = Math.abs(botCoords[i].xCenter - StatUtils.percentile(x, 50)) / StatUtils.percentile(x, 50);
				if (xc > jCFS.percentToleratedDifference / 100) isCool = false;
				
				double yc = Math.abs(botCoords[i].yCenter - StatUtils.percentile(y, 50)) / StatUtils.percentile(y, 50);
				if (yc > jCFS.percentToleratedDifference / 100) isCool = false;
								
				double wbc = Math.abs(botCoords[i].CVOfWallBrightness - StatUtils.percentile(wb, 50)) / StatUtils.percentile(wb, 50);
				if (wbc > 10 * jCFS.percentToleratedDifference / 100) isCool = false;  //criterion needs to be relaxed for the wall brightness coefficient of variation
				
				if (isCool) toAdd.add(i);
			}
			
			//if toAdd is empty, try some slice further up
			if (toAdd.size() < minAddSize) {
				referenceFrame[0] = nowOnBot + validationRange - 1;
				referenceFrame[1] = originallyOnBot;
				toAdd.clear();
				soFarOnBot = originallyOnBot;
				badBots.add(valiSeries[valiSeries.length - 1]);
				//jCFS.debug = true;
			}		
		}
		
		//create ColCords3D for output
		int botOfColumn = (int)Math.round(botCoords[toAdd.get(0)].zCenter);
		int colHeight = botOfColumn - samCoords.topOfColumn;
		ColCoords3D outCoords = initColCords3D(colHeight);
				
		//check or duplicates between topCoords and samCoords
		ArrayList<Integer> botDepths = new ArrayList<Integer>();
		for (int i = 0 ; i < botDepths.size() ; i++) botDepths.add((int)Math.round(botCoords[toAdd.get(i)].zCenter));
		int notAdd = 0;
		ArrayList<Integer> noAdd = new ArrayList<Integer>();
		for (int i = 0 ; i < samCoords.zmid.length ; i++) {
			if ((int)Math.round(samCoords.zmid[i]) >= botOfColumn | (int)Math.round(samCoords.zmid[i]) == 0 |
					botDepths.contains(samCoords.zmid[i])) {
				notAdd++;
				noAdd.add(i);
			}
		}
		
		//plug bottom coordinates into output column coordinates
		for (int i = 0 ; i < toAdd.size() ; i++) {			
			int placing = (int)Math.round(botCoords[toAdd.get(i)].zCenter) - samCoords.topOfColumn - 1;
			outCoords = plug2DInto3DCoords(placing, botCoords[toAdd.get(i)], outCoords);			
		}
		
		//also fill in already found column coordinates 
		int[] goodSlices = new int[samCoords.zmid.length - notAdd];
		int[] placings = new int[samCoords.zmid.length - notAdd];
		int cc = 0;
		for (int i = 0 ; i < samCoords.zmid.length ; i++) if (!noAdd.contains(i)) {
			goodSlices[cc] = i;
			placings[cc] = (int)Math.round(samCoords.zmid[i]) - samCoords.topOfColumn; 
			cc++;
		}
			
		outCoords = plug3DInto3DCoords(placings, goodSlices, samCoords, outCoords);
		outCoords.topOfColumn = samCoords.topOfColumn;
		outCoords.bottomOfColumn = botOfColumn;
		outCoords.heightOfColumn = botOfColumn - samCoords.topOfColumn;	
		outCoords.anglesChecked = samCoords.anglesChecked;	
					
		return outCoords;
		
	}
	
	public ColCoords3D findClosestXYSlice2Bottom(InputOutput.MyFileCollection mFC, ColCoords3D samCoords, MenuWaiter.ColumnFinderMenuReturn jCFS) {
		
		IJ.showStatus("Looking for column coordinates close to the bottom ...");
				
		InputOutput jIO = new InputOutput();
		
		//find topmost found column part in samCoords
		int soFarOnBot = jCFS.bottomOfColumn;
		int botNumber = 0;
		for (int i = samCoords.columnIsAtThisDepth.length - 2; i > 0 ; i--) { // throw away the first find because it is often suboptimal
			if (samCoords.columnIsAtThisDepth[i]) {
				botNumber = i;
				break;			
			}
		}
		
		//get 2D coords of slice that is so far on top
		ColCoords2D nowBot = extract2DCoordsFrom3DCoords(samCoords, botNumber);
		ColCoords2D outBot = nowBot;
		
		ArrayList<Integer> toAdd = new ArrayList<Integer>();
		int validationRange = mFC.nOfSlices / 100 + 1; // check the botmost x slices
		int[] valiSeries = new int[validationRange];
		ColCoords2D[] botCoords = new ColCoords2D[validationRange];
		int originallyOnBot = soFarOnBot;	
		int nowOnBot = jCFS.bottomOfColumn;
		if (nowOnBot < 2) nowOnBot = mFC.nOfSlices;
		
		ArrayList<Integer> badBots = new ArrayList<Integer>();	
		
		int minAddSize = (int)Math.floor(jCFS.bottomOfColumn / 800) + 1;
		if (minAddSize == 1) minAddSize = 2;
		
		while (toAdd.size() < minAddSize) {
			
			//sample x images from bot and check when the gray values stabilize..				
			int stopSlice = nowOnBot - validationRange;			
			if (badBots.size() > 0) {
				while (badBots.get(badBots.size()-1) - stopSlice < 5) {
					stopSlice -= 7;
				}
			}
			for (int i = stopSlice + validationRange ; i > stopSlice ; i--) valiSeries[validationRange - (i - stopSlice)] = i;
			
			ImagePlus valiTiff = jIO.openTiff3DSomeSlices(mFC.nowTiffPath, mFC.nowWidth, mFC.nowHeight, valiSeries);
			
			ArrayList<Integer> z = new ArrayList<Integer>();
			
			for (int i = 0 ; i < valiTiff.getNSlices() ; i++) {
							
				valiTiff.setPosition(i + 1);
				ImageProcessor myIP = valiTiff.getProcessor();
				botCoords[i] = findColumnWalls2D(valiSeries[i], myIP, outBot, jCFS); 
				
				if (botCoords[i].columnIsAtThisDepth) {
					z.add(i);
				}
			}
					
			//check col properties
			int cc = 0;
			int startCheck = 2 * validationRange / 3;
			for (int i = startCheck ; i < validationRange; i++) if (z.contains(i)) cc++;
			
			//if not enough good values were found there, try a larger range..
			if (cc < 3) {
				cc = 0;
				startCheck = 1 * validationRange / 3;
				for (int i = startCheck ; i < validationRange; i++) if (z.contains(i)) cc++;		
			}
			
			//if still not enough good values were found there, try a larger range..
			if (cc < 3) {
				cc = 0;
				startCheck = 0;
				for (int i = startCheck ; i < validationRange; i++) if (z.contains(i)) cc++;		
			}
					
			//get values in reference series
			ArrayList<Integer> zRef = new ArrayList<Integer>();
			for (int i = startCheck ; i < validationRange; i++) {
				if (z.contains(i)) {
					for (int j = 0 ; j < z.size() ; j++) {
						if (z.get(j) >= startCheck & !zRef.contains(z.get(j))) {
							zRef.add(z.get(j));					
							break;
						}
					}
				}
			}
			
			double[] x = new double[zRef.size()];
			double[] y = new double[zRef.size()];
			double[] r = new double[zRef.size()];		
			
			for (int i = 0 ; i < zRef.size() ; i++) { //TAKE MEDAIN OF LAST values as reference
				x[i] = botCoords[zRef.get(i)].xCenter;
				y[i] = botCoords[zRef.get(i)].yCenter;
				r[i] = (botCoords[zRef.get(i)].innerMajorRadius + botCoords[i].innerMinorRadius) / 2;	
			}
			
			//check start of column
			boolean isCool = true;
			for (int i = 0 ; i < validationRange; i++) {
				
				isCool = false;		
				if (z.contains(i)) isCool = true;
					
				double xc = Math.abs(botCoords[i].xCenter - StatUtils.percentile(x, 50)) / StatUtils.percentile(x, 50);
				if (xc > jCFS.percentToleratedDifference / 100) isCool = false;
				
				double yc = Math.abs(botCoords[i].yCenter - StatUtils.percentile(y, 50)) / StatUtils.percentile(y, 50);
				if (yc > jCFS.percentToleratedDifference / 100) isCool = false;
				
				if (isCool) toAdd.add(i);
			}
			
			//if toAdd is empty, try some slice further up
			if (toAdd.size() < minAddSize) {
				toAdd.clear();
				soFarOnBot = originallyOnBot;
				badBots.add(valiSeries[valiSeries.length - 1]);
				//jCFS.debug = true;
				nowOnBot -= 7;
				outBot = botCoords[7];
			}		
		}
		
		//create ColCords3D for output
		int botOfColumn = jCFS.bottomOfColumn;
		if (botOfColumn == 1) botOfColumn = mFC.nOfSlices;
		int colHeight = botOfColumn - samCoords.topOfColumn + 1;
		ColCoords3D outCoords = initColCords3D(colHeight);
				
		//check or duplicates between topCoords and samCoords
		ArrayList<Integer> botDepths = new ArrayList<Integer>();
		for (int i = 0 ; i < botDepths.size() ; i++) botDepths.add((int)Math.round(botCoords[toAdd.get(i)].zCenter));
		int notAdd = 0;
		ArrayList<Integer> noAdd = new ArrayList<Integer>();
		for (int i = 0 ; i < samCoords.zmid.length ; i++) {
			if ((int)Math.round(samCoords.zmid[i]) >= botOfColumn | (int)Math.round(samCoords.zmid[i]) == 0 |
					botDepths.contains(samCoords.zmid[i])) {
				notAdd++;
				noAdd.add(i);
			}
		}
		
		if (botCoords[toAdd.get(0)].zCenter < jCFS.bottomOfColumn) {
			botCoords[toAdd.get(0)].zCenter = jCFS.bottomOfColumn;
			botCoords[toAdd.get(0)].outerR2 = 0;
			botCoords[toAdd.get(0)].innerR2 = 0;
		}
		
		for (int i = 0 ; i < toAdd.size() ; i++) {			
			int placing = (int)Math.round(botCoords[toAdd.get(i)].zCenter) - samCoords.topOfColumn;
			outCoords = plug2DInto3DCoords(placing, botCoords[i], outCoords);			
		}
		
		int[] goodSlices = new int[samCoords.zmid.length - notAdd];
		int[] placings = new int[samCoords.zmid.length - notAdd];
		int cc = 0;
		for (int i = 0 ; i < samCoords.zmid.length ; i++) if (!noAdd.contains(i)) {
			goodSlices[cc] = i;
			placings[cc] = (int)Math.round(samCoords.zmid[i]) - samCoords.topOfColumn; 
			cc++;
		}
			
		outCoords = plug3DInto3DCoords(placings, goodSlices, samCoords, outCoords);
		outCoords.bottomOfColumn = botOfColumn;
		outCoords.topOfColumn = samCoords.topOfColumn;
		outCoords.heightOfColumn = botOfColumn - samCoords.topOfColumn;	
		outCoords.anglesChecked = samCoords.anglesChecked;	
					
		return outCoords;
		
	}
	
	public int[] findTopAndBottomOfPVCColumn(ImagePlus nowTiff, ColCoords3D pCC) {
		
		//init objects		
		Median jMed = new Median();
		HistogramStuff hist = new HistogramStuff();

		//init outline related variables		
		int i;
		int[] topAndBot = new int[2];
		double imageHeight = nowTiff.getNSlices();
		double[] myP = new double[(int)imageHeight];
		double[] dP = new double[(int)imageHeight - 1];
		double threshold = 0.99;

		for (i = 0 ; i < imageHeight ; i++) {
		
			IJ.showStatus("Searching slice for column's top and bottom #" + (i + 1) + "/" + imageHeight);
		
			nowTiff.setPosition(i+1);				
			ImageProcessor myIP = nowTiff.getProcessor();
			
			int[] myHist = myIP.getHistogram();						
			double[] cumHist = hist.calcCumulativeHistogram(myHist);
		
			//calculate 99 percentile of gray values at this depth
			myP[i] = hist.findPercentileFromCumHist(cumHist, threshold);
			
			//calculate 1st derivative of 99 percentile profile
			if (i > 0) dP[i-1] = myP[i] - myP[i-1];
			
		}
		
		//search maximum and minimum derivatives
		double minDP = 0;
		double maxDP = 0;
		
		for (i = 0 ; i < imageHeight - 1 ; i++) {
			if (dP[i] < minDP & dP[i] < -0.1 * jMed.evaluate(myP)) {
				minDP = dP[i];
				topAndBot[0] = i;
			}
			if (dP[i] > maxDP & dP[i] > 0.1 * jMed.evaluate(myP)) {
				maxDP = dP[i];
				topAndBot[1] = i;
			}
			
		}
		
		return topAndBot;
		
	}
	
	public int findMedianSoilSurfacePosition(ImagePlus surfTiff) {
		
		HistogramStuff hist = new HistogramStuff();
		
		surfTiff.setPosition(1);
		ImageProcessor topSoilIP = surfTiff.getProcessor();
		int[] topSurfHist = topSoilIP.getHistogram();
		topSurfHist[0] = 0;
		return hist.findMedianFromHistogram(topSurfHist);		
	
	}
	
	/*public ColumnCoordinates findTopAndBottomOfSteelColumn(ImagePlus nowTiff, ColumnCoordinates pCC) {
		
		//init objects
		Median jMed = new Median();
		HistogramStuff hist = new HistogramStuff();
				
		//init outline related variables
		int i;
		double imageHeight = nowTiff.getNSlices();
		double[] myP = new double[(int)imageHeight];
		double[] dP = new double[(int)imageHeight - 1];
		double threshold = 0.99;

		for (i = 0 ; i < imageHeight ; i++) {
		
			IJ.showStatus("Searching slice for column's top and bottom #" + (i + 1) + "/" + imageHeight);
		
			nowTiff.setPosition(i+1);				
			ImageProcessor myIP = nowTiff.getProcessor();
			
			int[] myHist = myIP.getHistogram();						
			double[] cumHist = hist.calcCumulativeHistogram(myHist);
					
			//calculate 99 percentile of gray values at this depth
			myP[i] = hist.findPercentileFromCumHist(cumHist, threshold);
			
			//calculate 1st derivative of 99 percentile profile
			if (i > 0) dP[i-1] = myP[i] - myP[i-1];
			
		}
		
		//search maximum and minimum derivatives
		double minDP = 0;
		double maxDP = 0;
		pCC.bottomOfColumn = myP.length; 
		pCC.topOfColumn = 0;
		
		for (i = 0 ; i < imageHeight - 1 ; i++) {
			if (dP[i] < minDP & myP[i] < -0.1 * jMed.evaluate(myP)) {
				minDP = dP[i];
				pCC.bottomOfColumn = i;
				
				if (pCC.bottomOfColumn - pCC.topOfColumn < 0.1 * pCC.approximateHeight) {
					maxDP = 0;
					minDP = 0;
				}
				
			}
			if (dP[i] > maxDP & myP[i] > 0.1 * jMed.evaluate(myP)) {
				maxDP = dP[i];
				pCC.topOfColumn = i;
			}
			
		}
		
		return pCC;
		
	}*/
	
	public ImagePlus makeAggregateMask(InputOutput.MyFileCollection mFC, RoiHandler.ColumnRoi colRoi, MenuWaiter.AggregateMaskOptionMenu aMO) {
		
		ImageManipulator jIM = new ImageManipulator();
		Dilate_ mD = new Dilate_();
		Erode_ mE = new Erode_();
		FloodFillComponentsLabeling3D myFCCL = new FloodFillComponentsLabeling3D(26, 32);
				
		ImagePlus maskTiff = new ImagePlus();	
		
		double roiVolume = colRoi.area * colRoi.nowTiff.getNSlices();
		
		//invert
		IJ.showStatus("Inverting image ... ");
		ImagePlus aggTiff = jIM.invertImage(colRoi.nowTiff, null);
		IJ.freeMemory();IJ.freeMemory();
		
		//filter out small particles
		int minimalParticleSize = (int)Math.round(aMO.filterSize / 100 * roiVolume);		
		ImageStack myLabelStack = myFCCL.computeLabels(aggTiff.getStack());			
		aggTiff.setStack(BinaryImages.binarize(LabelImages.volumeOpening(myLabelStack, minimalParticleSize)));		
		
		ImagePlus virginTiff = aggTiff.duplicate();		
		IJ.freeMemory();IJ.freeMemory();
		
		//dilation
		IJ.showStatus("Dilating ... (step 1 / " + aMO.closingVoxels + ")");
		aggTiff = mD.dilate(aggTiff, 255, false);
		for (int i = 1 ; i < aMO.closingVoxels ; i++) {
			IJ.showStatus("Dilating ... (step " + (i + 1) + " / " + aMO.closingVoxels + ")");
			aggTiff = mD.dilate(aggTiff, 255, false);
		}
		IJ.freeMemory();IJ.freeMemory();
		
		//fill holes
		IJ.showStatus("Removing 'holes' from the aggregate ... ");
		aggTiff = jIM.removeHoles(aggTiff);
		virginTiff = jIM.removeHoles(virginTiff);
		IJ.freeMemory();IJ.freeMemory();
		
		//invert once more
		IJ.showStatus("Inverting image ... ");
		aggTiff = jIM.invertImage(aggTiff, null);
		
		//watershed
		IJ.showStatus("Distance watershed transform ... ");
		//aggTiff = jIM.distanceTransformWatershed(aggTiff, aMO);
		IJ.freeMemory();IJ.freeMemory();

		aggTiff.updateAndDraw();
		aggTiff.show();
		
		//delete all watersheds that touch the 3D image canvas edges
		IJ.showStatus("Removing spurious watersheds ... ");
		aggTiff = jIM.deleteWatershedsAtEdges(aggTiff);
		
		//invert even once more
		IJ.showStatus("Inverting image ... ");
		aggTiff = jIM.invertImage(aggTiff, null);
		
		//erosion
		int erosionSteps = (int)Math.round(aMO.erosionOvershoot * (double)aMO.closingVoxels);
		IJ.showStatus("Eroding ... (step 1 / " + erosionSteps + ")");
		ImagePlus eroTiff =  mE.erode(aggTiff, 255, false);
		for (int i = 1 ; i < erosionSteps ; i++) {
			IJ.showStatus("Eroding ... (step " + (i + 1) + " / " + erosionSteps + ")");
			eroTiff =  mE.erode(eroTiff, 255, false);
		}
		
		eroTiff.updateAndDraw();
		eroTiff.show();
		
		//fuse with the virgin aggTiff
		ImagePlus fuseTiff = jIM.fuseMasks(virginTiff, eroTiff);
		
		fuseTiff.updateAndDraw();
		fuseTiff.show();
		
		//release memory
		aggTiff.killStack();eroTiff.killStack();virginTiff.killStack();
				
		//cut away small remnants from the dilation-watershed-erosion tribulation..		
		myLabelStack = myFCCL.computeLabels(fuseTiff.getStack());			
		maskTiff.setStack(BinaryImages.binarize(LabelImages.volumeOpening(myLabelStack, minimalParticleSize)));		
			
		return maskTiff;
		
	}	
	
	public ImagePlus goWatershed(InputOutput.MyFileCollection mFC, RoiHandler.ColumnRoi colRoi, MenuWaiter.TobiasWatershedOptionMenu aMO) {

		ImageManipulator jIM = new ImageManipulator();
		Dilate_ mD = new Dilate_();
		Erode_ mE = new Erode_();
		FloodFillComponentsLabeling3D myFCCL = new FloodFillComponentsLabeling3D(26, 32);
		
		ImagePlus aggTiff = colRoi.nowTiff;	
		
		double roiVolume = colRoi.area * colRoi.nowTiff.getNSlices();
		
		if (aMO.invertImage) {
			//invert		
			IJ.showStatus("Inverting image ... ");
			aggTiff = jIM.invertImage(aggTiff, null);
			IJ.freeMemory();IJ.freeMemory();
		}

		//cut away small remnants from the dilation-watershed-erosion tribulation..
		int minParticleSize = (int)Math.round(aMO.filterSize*roiVolume/100);
		ImageStack myLabelStack = myFCCL.computeLabels(aggTiff.getStack());			
		aggTiff.setStack(BinaryImages.binarize(LabelImages.volumeOpening(myLabelStack, minParticleSize)));
		
		//invert		
		IJ.showStatus("Inverting image ... ");
		aggTiff = jIM.invertImage(aggTiff, null);
		IJ.freeMemory();IJ.freeMemory();
		
		//watershed
		IJ.showStatus("Distance watershed transform ... ");
		aggTiff = jIM.distanceTransformWatershed(aggTiff, aMO);
		IJ.freeMemory();IJ.freeMemory();
		
		//delete all watersheds that touch the 3D image canvas edges
		IJ.showStatus("Removing spurious watersheds ... ");
		aggTiff = jIM.deleteWatershedsAtEdges(aggTiff);		
		
		//invert even once more
		IJ.showStatus("Inverting image ... ");
		aggTiff = jIM.invertImage(aggTiff, null);
		
		//erosion
		int erosionSteps = 2;
		IJ.showStatus("Eroding ... (step 1 / " + erosionSteps + ")");
		ImagePlus eroTiff =  mE.erode(aggTiff, 255, false);
		for (int i = 1 ; i < erosionSteps ; i++) {
			IJ.showStatus("Eroding ... (step " + (i + 1) + " / " + erosionSteps + ")");
			eroTiff =  mE.erode(eroTiff, 255, false);
		}
		
		//dilation
		IJ.showStatus("Dilating ... (step 1 / " + 2 + ")");
		aggTiff = mD.dilate(aggTiff, 255, false);
		for (int i = 1 ; i < 2 ; i++) {
			IJ.showStatus("Dilating ... (step " + (i + 1) + " / " + 2 + ")");
			aggTiff = mD.dilate(aggTiff, 255, false);
		}
		IJ.freeMemory();IJ.freeMemory();
		
		//cut away small remnants from the dilation-watershed-erosion tribulation..
		myLabelStack = myFCCL.computeLabels(aggTiff.getStack());			
		aggTiff.setStack(BinaryImages.binarize(LabelImages.volumeOpening(myLabelStack, minParticleSize)));
			
		return aggTiff;
		
	}
	
	public ImagePlus makeAggregateMask(InputOutput.MyFileCollection mFC, RoiHandler.ColumnRoi colRoi, MenuWaiter.ClosingMaskOptionMenu aMO) {
		
		MorphologyAnalyzer morph = new MorphologyAnalyzer();
		ImageManipulator jIM = new ImageManipulator();
		Dilate_ mD = new Dilate_();
		Erode_ mE = new Erode_();
		FloodFillComponentsLabeling3D myFCCL = new FloodFillComponentsLabeling3D(26, 32);
	
		double roiVolume = colRoi.area * colRoi.nowTiff.getNSlices();
		
		//invert
		IJ.showStatus("Inverting image ... ");
		ImagePlus aggTiff = jIM.invertImage(colRoi.nowTiff, null);
		IJ.freeMemory();IJ.freeMemory();
		
		//filter out small particles
		int minimalParticleSize = (int)Math.round(aMO.filterSize / 100 * roiVolume);
		ImageStack myLabelStack = myFCCL.computeLabels(aggTiff.getStack());			
		aggTiff.setStack(BinaryImages.binarize(LabelImages.volumeOpening(myLabelStack, minimalParticleSize)));	
		
		ImagePlus virginTiff = aggTiff.duplicate();	
		IJ.freeMemory();IJ.freeMemory();
		
		//dilation
		IJ.showStatus("Dilating ... (step 1 / " + aMO.closingVoxels + ")");
		aggTiff = mE.erode(aggTiff, 255, false);
		for (int i = 1 ; i < aMO.closingVoxels ; i++) {
			IJ.showStatus("Dilating ... (step " + (i + 1) + " / " + aMO.closingVoxels + ")");
			aggTiff = mD.dilate(aggTiff, 255, false);
		}
		IJ.freeMemory();IJ.freeMemory();
		
		//fill holes
		IJ.showStatus("Removing 'holes' from the aggregate ... ");
		aggTiff = jIM.removeHoles(aggTiff);
		virginTiff = jIM.removeHoles(virginTiff);	
		IJ.freeMemory();IJ.freeMemory();
		
		ImagePlus eroTiff = aggTiff;
		//erosion
		if (!aMO.noErosion){
			int erosionSteps = (int)Math.round(aMO.erosionOvershoot * (double)aMO.closingVoxels);		
			IJ.showStatus("Eroding ... (step 1 / " + erosionSteps + ")");
			eroTiff =  mE.erode(aggTiff, 255, false);
			for (int i = 1 ; i < erosionSteps ; i++) {
				IJ.showStatus("Eroding ... (step " + (i + 1) + " / " + erosionSteps + ")");
				eroTiff =  mE.erode(eroTiff, 255, false);
			}
		}
		
		//fuse with the virgin aggTiff
		eroTiff = jIM.fuseMasks(virginTiff, eroTiff);
		
		return eroTiff;
		
	}	
	

	
	
	
	public BestWallFinds findBestWallDetectionPerSector(int numOfSectors, int i, float[] xOD, float[] yOD, double[] myAngle, ColCoords3D prelimCC) {
				
		BestWallFinds bWF = new BestWallFinds();
		
		int numOfAngles = myAngle.length;
		int increment = numOfAngles / numOfSectors;
		
		double xMid = StatUtils.percentile(prelimCC.xmid, 50);
		double yMid = StatUtils.percentile(prelimCC.ymid, 50);
		
		ArrayList<Integer> bestOnes = new ArrayList<Integer>();
		
		for (int j = 0 ; j < numOfAngles ; j += increment) {
			
			int bestOne = -1;
			double farthestDist = -1;
			
			for (int k = 0 ; k < increment ; k++) {
			
				double xDist = xOD[j + k] - xMid;
				double yDist = yOD[j + k] - yMid;				
				double nowDist = Math.sqrt(xDist * xDist + yDist * yDist);
				
				if (nowDist > farthestDist) {
					
					farthestDist = nowDist;
					bestOne = j + k;
					
				}		
			}
			
			if (bestOne > 1) {
				bestOnes.add(bestOne);
			}
			
		}
		
		float[] xD = new float[bestOnes.size()];
		float[] yD = new float[bestOnes.size()];
		double[] angle = new double[bestOnes.size()];
		
		for (int j = 0 ; j < bestOnes.size() ; j++) {
			
			xD[j] = xOD[bestOnes.get(j)];
			yD[j] = yOD[bestOnes.get(j)];
			angle[j] = myAngle[bestOnes.get(j)];
			
		}
		
		bWF.xD = xD;
		bWF.yD = yD;
		bWF.angle = angle;
		
		return bWF;
	}
	
	public ColCoordsEssentials3D extractColCoordsEssentials3D(ColCoords3D jCO) {
		
		ColCoordsEssentials3D outCoords = new ColCoordsEssentials3D();
		TailoredMaths math = new TailoredMaths();
		
		outCoords.anglesChecked = jCO.anglesChecked;
		
		//create array of considered angles 
		int angleCounter = 0;			
		double[] myAngles = new double[jCO.anglesChecked];
		int maxAlpha = 360;
		int dAlpha = maxAlpha / jCO.anglesChecked;		
		for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {				
			
			myAngles[angleCounter] = angle;
			
			angleCounter++;
		}
		
		//calculate x and y coordinate
		double[][] X = new double[jCO.ixmid.length][jCO.anglesChecked];
		double[][] Y = new double[jCO.ixmid.length][jCO.anglesChecked];
		for (int i = 0 ; i < jCO.ixmid.length ; i++) {
			
			double[][] xy = math.getXYOfEllipseFromAngle(myAngles, jCO.ixmid[i], jCO.iymid[i], jCO.innerMajorRadius[i], jCO.innerMinorRadius[i], jCO.itheta[i]);
						
			//cast it into two arrays..
			for (int j = 0 ; j < myAngles.length ; j++) X[i][j] = xy[j][0];
			for (int j = 0 ; j < myAngles.length ; j++) Y[i][j] = xy[j][1];
			
		}
		
		outCoords.xID = X;
		outCoords.yID = Y;		
		
		outCoords.xmid = jCO.ixmid;
		outCoords.ymid = jCO.iymid;
		
		
		return outCoords;
	}
	
	public ColCoordsEssentials3D extractColCoordsEssentials3D(EggShapedColCoords3D jCO) {
		
		ColCoordsEssentials3D outCoords = new ColCoordsEssentials3D();
		
		outCoords.anglesChecked = jCO.anglesChecked;
		
		outCoords.xID = jCO.xID;
		outCoords.yID = jCO.yID;
		
		outCoords.xmid = jCO.xmid;
		outCoords.ymid = jCO.ymid;
		
		return outCoords;
	}
	

	public ColCoords2D findColumnWalls2D(int sliceNum, ImageProcessor myIP, ColCoords2D prelimCC, MenuWaiter.ColumnFinderMenuReturn jCFS) {
	// takes image, slice number, preliminary Column guesses and menu inputs to return the 2D coordinates	
		
//		ImagePlus test = new ImagePlus("",myIP);
//		test.updateAndDraw();
//		test.show();
//		
//		IJ.wait(500);
//		
//		test.close();
		
		//import objects
		TailoredMaths maths = new TailoredMaths();
		FitStuff fit = new FitStuff();		
		RoiHandler roi = new RoiHandler();
		ColCoords2D i2D = new ColCoords2D();
		AndAllTheRest aa = new AndAllTheRest();
		
		//init variables
		int maxAlpha = 360;
		int dAlpha = 5;
		double[] myAngle = new double[maxAlpha / dAlpha];		
		int cc = 0;	
		for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha / dAlpha)) 
			{myAngle[cc] = angle; cc++;}
	
		int footprintOfMedianFilter = 5; // 1-D median filter		// settings for median filter: smaller "window" for smaller images
		if (myIP.getWidth() < 1200) footprintOfMedianFilter = 3;
		
		double xmid = prelimCC.xCenter;
		double ymid = prelimCC.yCenter;
		double radius = (prelimCC.outerMajorRadius + prelimCC.outerMinorRadius) / 2;
		
		//catch if prelimCC did not work
		if (xmid != xmid) xmid = myIP.getWidth() / 2;				// if there are no preliminary Oulines, take the center of the image as ROI center. and assume radius of 
		if (ymid != ymid) ymid = myIP.getHeight() / 2;
		if (radius != radius) radius = 0.85 * (myIP.getWidth() + myIP.getHeight()) / 4;
		
		double startingPoint = Math.round(radius) + 40;
		int myWindowSize = (int)Math.ceil(0.003 * radius) + 1;

		//init coordinates of found wall / not-wall interfaces
		float[] xOD = new float[maxAlpha/dAlpha]; 			// x of outer diameter (what is alpha?? some angles of 360 and 5 but I dont get it)
		float[] yOD = new float[maxAlpha/dAlpha]; 			// y of outer diameter
		int[] locationOD = new int[maxAlpha/dAlpha]; 		// y of outer diameter ?????
		float[] xID = new float[maxAlpha/dAlpha]; 			// x of inner diameter
		float[] yID = new float[maxAlpha/dAlpha]; 			// y of inner diameter
		double[] locationID = new double[maxAlpha/dAlpha]; 	// y of inner diameter ?????
	
		//init startin values for finding the wall
		double oRelAirWallContrast = prelimCC.airWallContrast;
		double oRDStdThresh = prelimCC.wallSoilStdContrastThreshold;	
			
		//init variables to tamper with..
		double relativeAirWallContrast = oRelAirWallContrast; //standard contrast for outer Wall find
		double rDStdThresh = oRDStdThresh;	//standard contrast for inner Wall find
		double adaptedCVWallThicknessThresh = jCFS.maxCVOfWallThickness; 	
		double stdThreshold = jCFS.stdThreshold;
		
		//init parameter bounds
		double rAWC_Top = 3;  //relativeAirWallContrast
		double rAWC_Bot = 0.1;
		double dContrast = 0.25;  //updater for relative air wall contrast
				
		//if the wall is PVC..
/*		if (jCFS.isPVC) {
			relativeAirWallContrast = 0.1;	
			rAWC_Top = 0.4;
			rAWC_Bot = 0.01;
			absContrastAirWall = 1000;
			dContrast = 0.04;
		}*/
		
		//find illumination level 			
		HistogramStuff.IlluminationInfo jII = probeIlluminationLevel(myIP, prelimCC, true);
		double expectedCVOfWallgrayValues = 0.5; //empirically.. (formerly: colCon.jCFS.convInterCVOfWallBrightness * colCon.jCFS.coeffVarOfWallBrightness;)
		int myContrast = jII.quantile99 - jII.lowerQuartile; //evaluates image containing the outside of the wall and half the wall... therefore, q99 corresponds to wall brightness2
		int minAirColGradient = (int) Math.round(relativeAirWallContrast * (double) myContrast);			
		int minWallBrightness = (int) Math.round((double)(1 - expectedCVOfWallgrayValues) * (double)jII.quantile99);
		int maxWallBrightness = (int) Math.round((double)(1 + expectedCVOfWallgrayValues) * (double)jII.quantile99);
		
		//adjust min and max wall brightnesses to predefined values
		if (jCFS.minColGV > 0) minWallBrightness = jCFS.minColGV;
		if (jCFS.maxColGV > 0) maxWallBrightness = jCFS.maxColGV;
		
		//set a lower threshold for inner edge detection
		double lowerStandardThreshold4InnerPeri = minAirColGradient / 5;
		double upperStandardThreshold4InnerPeri = minAirColGradient * 5;
		
		//xy coordinates of radial wall-interface search-lines
		double[][] x = new double[myAngle.length][(int)Math.round(startingPoint + 1)];
		double[][] y = new double[myAngle.length][(int)Math.round(startingPoint + 1)];
		
		boolean[] innerPerimeterFound = new boolean[myAngle.length];
		boolean[] outerPerimeterFound = new boolean[myAngle.length];
		
		//start the radial line search
		int trialCounter = 0;
		FitStuff.FittedEllipse fO = null;
		FitStuff.FittedEllipse fI = null;
		double R2 = 0;
		
		//remember thresholds that worked
		ArrayList<Integer> triedThreshes = new ArrayList<Integer>(); // remmber whcih thresholds have been tried out
		
		while (fO == null | (R2 < jCFS.r2Thresh & trialCounter <= jCFS.maxFittingAttempts)) {
			
			trialCounter++;
			
			for (int j = 0 ; j < myAngle.length ; j++) {			// angle from which to find outside wall
			
				double[] grayAtThisAngle = new double[(int)Math.round(startingPoint + 1)];
				locationOD[j] = 0; //nil it!
				
				//re-find outer column wall
				cc=0;
				for (double checkRadius = startingPoint ; checkRadius > 0 ; checkRadius--) {					
					x[j][cc] = Math.sin(myAngle[j]) * checkRadius + xmid;
					y[j][cc] = Math.cos(myAngle[j]) * checkRadius + ymid;				
					grayAtThisAngle[cc] = myIP.getPixelValue((int)x[j][cc], (int)y[j][cc]);					
					cc++;
				}
				
				//apply medianFilter
				double[] mgrayAtThisAngle = maths.oneDMedianFilter(grayAtThisAngle, footprintOfMedianFilter);
				
				//calculate the first derivative, dMG
				double[] dMG = maths.oneDFirstDerivative(mgrayAtThisAngle, footprintOfMedianFilter);  //footprint is needed to know the amount of data filtered away from the fringes of the vector
				
				//claculate std of gray values
				double stdOfRadial = maths.stdWithoutZeros(mgrayAtThisAngle);
				
				//set wall found trackers to false
				innerPerimeterFound[j] = false;
				outerPerimeterFound[j] = false;				
				
				for (int k = myWindowSize + 4 ; k < (int)Math.round(startingPoint) - myWindowSize - 6; k++) {
					
					double[] theNextInTheWindow = new double[myWindowSize + 1];
					double[] thePrevInTheWindow = new double[myWindowSize + 1];
					for(int l = 0 ; l < theNextInTheWindow.length ; l++) {
						theNextInTheWindow[l] = mgrayAtThisAngle[k+l+2];
						thePrevInTheWindow[l] = mgrayAtThisAngle[k-l-2];
					}
					double meanPrev = StatUtils.mean(thePrevInTheWindow);
					double meanNext = StatUtils.mean(theNextInTheWindow);
					double dNextPrev = meanNext - meanPrev;
					
					
					//check whether candidate for outer wall location has been found
					boolean checkBool = (StatUtils.min(theNextInTheWindow) >= minWallBrightness & 
										StatUtils.max(theNextInTheWindow) <= maxWallBrightness & 
										dNextPrev > stdThreshold * stdOfRadial  &  
										StatUtils.min(thePrevInTheWindow) > 0 &
										locationOD[j] == 0);
					
					if (jCFS.isAlreadyNormalized) {
						
						boolean nextIsLarger = StatUtils.min(theNextInTheWindow) >= jCFS.fixedWallGrayValue - jCFS.stdFixedWallGrayValue;
						boolean prevIsSmaller = StatUtils.max(thePrevInTheWindow) <= jCFS.fixedWallGrayValue - jCFS.stdFixedWallGrayValue;
						
						checkBool = nextIsLarger & prevIsSmaller;
								
					}
							
					
					//if column is still found, find the maximum contrast in the vicinity and save the results 
					if (checkBool) {
						
						double maxD = dNextPrev;
						int myJ = k;
						
						for (int checkMore = k + 1; checkMore < k + 5 ; checkMore++) {
							
							for(int l = 0 ; l < theNextInTheWindow.length ; l++) {
								theNextInTheWindow[l] = mgrayAtThisAngle[k+l+2];
								thePrevInTheWindow[l] = mgrayAtThisAngle[k-l-2];
							}
							meanPrev = StatUtils.mean(thePrevInTheWindow);
							meanNext = StatUtils.mean(theNextInTheWindow);
							dNextPrev = meanNext - meanPrev;
										
							if (dNextPrev > maxD) {
								maxD = dNextPrev;
								myJ = checkMore;
							}
						}
						
						xOD[j] = (float)x[j][myJ];
						yOD[j] = (float)y[j][myJ];
						locationOD[j] = myJ;
						outerPerimeterFound[j] = true;
						break;
					}
				}
			}
			
			//fit an elliptic ROI to the x and y located on the outer wall perimeter
			fO = fit.doRobustEllipseFit(sliceNum, xOD, yOD, myAngle, myIP, jCFS, relativeAirWallContrast);
			R2 = fO.R2;	
			
			//displayFoundOutlines(sliceNum, myIP, fO, null, xOD, yOD, xID, yID, i2D, jCFS, adaptedCVWallThicknessThresh);
								
			//check if sufficient wall pieces were found
		/*	int missesCounter = 0;
			for (int k = 0 ; k < outerPerimeterFound.length ; k++) {
				if (!outerPerimeterFound[k]) missesCounter++;
			}
			if (missesCounter > jCFS.maxNumberOfOutliers4PerimeterFits) R2 = 0;*/
			
			//change column seeking parameters in case that column was not found..
			if (R2 != R2 | R2 < jCFS.r2Thresh) {
				double meanPrelimRadius = (prelimCC.outerMajorRadius + prelimCC.outerMinorRadius) / 2; 
				double meanFittedRadius = (fO.majorRadius + fO.minorRadius) / 2;
				
				//if wall edge was overlooked by some of the line searches..
				if (meanPrelimRadius > meanFittedRadius) { //try to decrease the detection threshold..
					if (relativeAirWallContrast > rAWC_Bot) {
						relativeAirWallContrast *= dContrast;
						if (relativeAirWallContrast < rAWC_Bot) relativeAirWallContrast = 0.98 * rAWC_Top;
					}						
				}
				
				//if wall edge detected prematurely by some of the line searches..
				if (meanPrelimRadius < meanFittedRadius) { //try to increase the detection threshold..
					if (relativeAirWallContrast <= rAWC_Top) {
						relativeAirWallContrast /= dContrast;
						if (relativeAirWallContrast > rAWC_Top) relativeAirWallContrast = 1.02 * rAWC_Bot;
					}						
				}
				
				//id R2 is NaN
				if (R2 != R2) {
					
					relativeAirWallContrast = rAWC_Top - relativeAirWallContrast;
					
				}
				
				minAirColGradient = (int)Math.round(relativeAirWallContrast * (double)myContrast);
				
				//if wall edge was overlooked by some of the line searches..
				if (meanPrelimRadius > meanFittedRadius) { //try to decrease the detection threshold..
					if (stdThreshold > rAWC_Bot) {
						stdThreshold *= dContrast;
						if (stdThreshold < rAWC_Bot) stdThreshold = rAWC_Bot;
					}						
				}
				
				//if wall edge detected prematurely by some of the line searches..
				if (meanPrelimRadius < meanFittedRadius) { //try to increase the detection threshold..
					if (stdThreshold < rAWC_Top) {
						stdThreshold /= dContrast;
						if (stdThreshold > rAWC_Top) stdThreshold = rAWC_Top;
					}						
				}
			}
			
		}
					
		//look for inner column interface .. only in case aluminium columns are used
		if (jCFS.isAlu & fO.R2 > jCFS.r2Thresh) {					
			
			//make nice set of outer edge locations
			int[] cOD = aa.kickOutStrangeValues(locationOD); // remove strange values
			cOD = aa.fillMissingInts(cOD, outerPerimeterFound); // remove non finds by average neighbors
			
			//reset loop trackers
			trialCounter=0;			
			R2 = 0;
						
			//try to find inner edge.			
			triedThreshes.add((int)Math.round(rDStdThresh));
			int[] choiceWasMade = new int[myAngle.length];
			for (int j = 0 ; j < choiceWasMade.length ; j++) choiceWasMade[j] = 0; // code: 0 nothing found; 1 classical; 2 absolute value; 3 both classical and abs
			
			ArrayList<Integer> remOutliers = new ArrayList<Integer>();
			int stepsize = 100; //stepsize for setting threshold for wall finding..
			int innerMaxFittingAttempts = jCFS.maxFittingAttempts;
			while (fI == null | (R2 < jCFS.r2Thresh & trialCounter <= innerMaxFittingAttempts)) {
				
				trialCounter++;
				int countFinds = 0;				
				
				for (int j = 0 ; j < myAngle.length ; j++) {
			
					double[] grayAtThisAngle = new double[(int)Math.round(startingPoint + 1)];
					innerPerimeterFound[j] = false; 
				
					//re-sample radial line for wall perimeter search
					cc=0;
					for (double checkRadius = startingPoint ; checkRadius > 0 ; checkRadius--) {					
						x[j][cc] = Math.sin(myAngle[j]) * checkRadius + xmid;
						y[j][cc] = Math.cos(myAngle[j]) * checkRadius + ymid;				
						grayAtThisAngle[cc] = myIP.getPixelValue((int)x[j][cc], (int)y[j][cc]);					
						cc++;
					}
					
					//apply medianFilter
					double[] mgrayAtThisAngle = maths.oneDMedianFilter(grayAtThisAngle, footprintOfMedianFilter);
					
				
/*					DisplayThings disp = new DisplayThings();
					double[] xaxis = new double[mgrayAtThisAngle.length];
					for (int k = 0 ; k < xaxis.length ; k++) xaxis[k] = k;
					if (sliceNum > 350) disp.plotXYXY(xaxis, mgrayAtThisAngle, xaxis, mgrayAtThisAngle, "", "r", "GV");*/
					
					
					//calculate the first derivative, dMG
					double[] dMG = maths.oneDFirstDerivative(mgrayAtThisAngle, footprintOfMedianFilter);  //footprint is needed to know the amount of data filtered away from the fringes of the vector
				
					double[] thePrevInTheWindow = new double[myWindowSize];
					double[] theNextInTheWindow = new double[myWindowSize];			
					
					double[] meanPrev = new double[(int)Math.round(startingPoint) / 4];
					double[] stdPrev = new double[(int)Math.round(startingPoint) / 4];
					double[] minPrev = new double[(int)Math.round(startingPoint) / 4];
					double[] maxPrev = new double[(int)Math.round(startingPoint) / 4];
					double[] meanNext = new double[(int)Math.round(startingPoint) / 4];
					double[] stdNext = new double[(int)Math.round(startingPoint) / 4];
					double[] minNext = new double[(int)Math.round(startingPoint) / 4];
					double[] maxNext = new double[(int)Math.round(startingPoint) / 4];
					double[] dMean = new double[(int)Math.round(startingPoint) / 4];
					double[] dStd = new double[(int)Math.round(startingPoint) / 4];
					double[] dMinMax = new double[(int)Math.round(startingPoint) / 4];
										
					int startSearching4Inner = cOD[j] + myWindowSize + 1;
					int stopSearching4Inner = (int)Math.round(startingPoint) / 4;	
														
					for (int k = startSearching4Inner ; k < stopSearching4Inner ; k++) {
						
						//collect values before and after k
						for(int l = 1 ; l < myWindowSize ; l++) thePrevInTheWindow[l] = mgrayAtThisAngle[k-l];
						for(int l = 1 ; l < myWindowSize ; l++) theNextInTheWindow[l] = mgrayAtThisAngle[k+l];						
						
						meanPrev[k] = StatUtils.mean(thePrevInTheWindow);
						stdPrev[k] = Math.sqrt(StatUtils.variance(thePrevInTheWindow));
						minPrev[k] = StatUtils.min(thePrevInTheWindow);
						maxPrev[k] = StatUtils.max(thePrevInTheWindow);
						
						meanNext[k] = StatUtils.mean(theNextInTheWindow);
						stdNext[k] = Math.sqrt(StatUtils.variance(theNextInTheWindow));
						minNext[k] = StatUtils.min(theNextInTheWindow);
						maxNext[k] = StatUtils.max(theNextInTheWindow);
						
						double minmax = Math.abs(maxNext[k] - minPrev[k]);
						double maxmin = Math.abs(maxPrev[k] - minNext[k]);
						
						dMinMax[k] = maxmin;
						if (minmax > maxmin) dMinMax[k] = minmax;
						
						dMean[k] = meanPrev[k] - meanNext[k];
						dStd[k] = stdPrev[k] - stdNext[k];						
					}
					
					//calculate wall abs and rel. mean and variation..
					double[] wallDMean = new double[2 * myWindowSize];					
					double[] wallDStd = new double[2 * myWindowSize];
					double[] wallAbs = new double[2 * myWindowSize];
					double[] wallGrad = new double[2 * myWindowSize];
					for (int k = startSearching4Inner; k < startSearching4Inner + 2 * myWindowSize ; k++) {
						wallDMean[k - (startSearching4Inner)] = dMean[k];
						wallDStd[k - (startSearching4Inner)] = dStd[k];
						wallAbs[k - (startSearching4Inner)] = mgrayAtThisAngle[k];
						wallGrad[k - (startSearching4Inner)] = dMG[k];
					}
					
					//calculate relative deviations of mean and std gray values..
					double[] rDMean = new double[dMean.length];
					double[] rDStd = new double[dStd.length];
					double[] rWallAbs = new double[dMean.length];
					double[] rGradient = new double[dStd.length];
					for (int k = startSearching4Inner; k < stopSearching4Inner ; k++) {
						rDMean[k] =  Math.abs(dMean[k] - StatUtils.percentile(wallDMean,50));
						rDStd[k] =  dStd[k] - Math.abs(StatUtils.percentile(wallDStd,50));
						rWallAbs[k] = Math.abs(mgrayAtThisAngle[k] - StatUtils.percentile(wallAbs,50));
						rGradient[k] = Math.abs(dMG[k] + Math.abs(StatUtils.percentile(wallGrad,50)));
					}		
						
					//if (sliceNum > 350) disp.plotXYXY(xaxis, rDStd, xaxis, rDStd, "", "r", "GV");
					//if (sliceNum > 350) disp.plotXYXY(xaxis, rDMean, xaxis, rDStd, "", "r", "GV");
					
					//find candidates for inner edge					
					for (int k = startSearching4Inner ; k < stopSearching4Inner ; k++) {						
													
						//check candidate in case of negative std difference double peak and pos mean peak
						if (rDStd[k] > rDStdThresh) {
							
							//define range in which peaks are sought						
							int myWindowEnd = k + myWindowSize;
							
							//if range exceed data series break
							if (myWindowEnd > rDStd.length) break;
							
							double peakDMean = 0;
							double minDStd = 0;
							double maxDStd = 0;
							double maxDAbs = 0;
							int peakPos = 0;
							int minPos = 0;
							int maxPos = 0;
							for (int l = k ; l < myWindowEnd ; l++) {															
								if (rDMean[l] > peakDMean) {
									peakDMean = rDMean[l];
									peakPos = l;
								}
								if (rDStd[l] < minDStd) {
									minDStd = rDStd[l];
									minPos = l;
								}
								if (Math.abs(rDStd[l]) > maxDStd) {
									maxDStd = Math.abs(rDStd[l]);
									maxPos = l;
								}		
															
							}	
							
							//calculate mean location of edge
							int meanPos = 0;
														
							//the classical										
							/*if (((minPos <= peakPos & peakPos <= maxPos) | 
								(minPos >= peakPos & peakPos >= maxPos)) & 
								(minDStd < - rDStdThresh & maxDStd > rDStdThresh)) {
								meanPos = Math.floorDiv(minPos + maxPos + peakPos,3);		
								choiceWasMade[j] = 1;
							}*/
							
							meanPos = maxPos;
							
							//in case the std does not work
							/*if (maxDAbs > rMeanThresh) {								
								if (choiceWasMade[j] == 0) {
									choiceWasMade[j] = 2;
									meanPos = peakPos;
								}
								else choiceWasMade[j] = 3;
							}	*/					
													
							
							if (meanPos > 0) {	
								xID[j] = (float)x[j][meanPos];
								yID[j] = (float)y[j][meanPos];
								locationID[j] = meanPos;
								innerPerimeterFound[j] = true; 
								countFinds++;
								break;								
							}													
																
							//disp.plotXYXY(xaxis, rDMean, xaxis, rDStd, "findInnerEdge", "x", "deltaDiff");													
							
						}
						
					}
					
					//display radial profile
					jCFS.showRadialProfiles = false; 
					//jCFS.debug = true;
					if (jCFS.showRadialProfiles & jCFS.debug) {
						DisplayThings dT = new DisplayThings();				
						double[] X = new double[dMG.length];
						for (int kk = 0 ; kk < X.length ; kk++) X[kk] = (double)kk;
						double[] foundPoints = new double[2];
						foundPoints[0] = locationID[j];foundPoints[1] = locationOD[j];
					
						dT.plotSpecialWallFinderRadii(X, mgrayAtThisAngle, dMG, foundPoints);
					}
				}	
						
				//in case it is an alu column, fit fI
				fI = fit.doRobustEllipseFit(sliceNum, xID, yID, myAngle, myIP, jCFS, rDStdThresh);				
				R2 = fI.R2;			
				
				//check quality of fit and wall find
				double wallThickness = 0.5 * (fO.majorRadius /  - fI.majorRadius + fO.minorRadius - fI.minorRadius);		
				double wallThicknessCV = calculateCVOfWallThickness(fI, fO);
				double CVOfWallBrightness = calculateCVOfWallBrightness(myIP, fI, fO);
				// allow larger wall Brightness thresolds for small wall thicknesses
				if (wallThickness < 30) adaptedCVWallThicknessThresh = (double)30 / wallThickness * jCFS.maxCVOfWallThickness;
				
				if (R2 > jCFS.r2Thresh) {			
					if (wallThickness < myWindowSize + 3) R2 = 0.99 * jCFS.r2Thresh;
					if (CVOfWallBrightness > jCFS.CVWallBrightnessThresh) R2 = 0.99 * jCFS.r2Thresh;
					if (wallThicknessCV > adaptedCVWallThicknessThresh) R2 = 0.99 * jCFS.r2Thresh; // also kick out the ones with too many zeros			
					if (xID.length - countFinds > jCFS.maxNumberOfOutliers4PerimeterFits) R2 = 0.99 * jCFS.r2Thresh; // also kick out the ones with too many zeros
				}
				
				///check what went wrong in case that things went wrong
				int prematures = 0;
				int overlooked = 0;
				int outliers = 0;
				if (fI.relDevFromMedianRadiusOfDiscPoints != null) {
					outliers = fI.relDevFromMedianRadiusOfDiscPoints.size();
					remOutliers.add(outliers);
					for (int l = 0 ; l < fI.relDevFromMedianRadiusOfDiscPoints.size() ; l++) {
						if (fI.relDevFromMedianRadiusOfDiscPoints.get(l) < 0) overlooked++;
						else prematures++;
					}
				}
				
				//if wall edge was overlooked by some of the line searches.							
				if (remOutliers.isEmpty()) {
					stepsize /= 2;
					rDStdThresh -= stepsize;			
					if (rDStdThresh < lowerStandardThreshold4InnerPeri) rDStdThresh = 3.4 * lowerStandardThreshold4InnerPeri; 
				}				
				else {
					if (remOutliers.size() > 1) {
				
						int olast = remOutliers.size() - 1;
						if (remOutliers.get(olast) - remOutliers.get(olast - 1) > 0) {
							stepsize /= 2;
							if (triedThreshes.size() > 1) {
								int lastThreshIndex = triedThreshes.size() - 1;
								rDStdThresh -= (triedThreshes.get(lastThreshIndex) - triedThreshes.get(lastThreshIndex - 1)) / 2;
							}
							
						}
					}
					else {
						if (prematures < overlooked) {   //try to decrease the detection threshold..	
						
							if (rDStdThresh > lowerStandardThreshold4InnerPeri) {
								rDStdThresh -= stepsize;
								while (triedThreshes.contains(rDStdThresh)) rDStdThresh -= stepsize;
								if (rDStdThresh <= lowerStandardThreshold4InnerPeri) {
									rDStdThresh = (lowerStandardThreshold4InnerPeri + 2*upperStandardThreshold4InnerPeri) / 3;
									while (triedThreshes.contains((int)Math.round(rDStdThresh))) rDStdThresh -= stepsize;
								}
								triedThreshes.add((int)Math.round(rDStdThresh));
							}											
						}
					
						//if wall edge detected prematurely by some of the line searches..
						if (overlooked <= prematures ) { //try to increase the detection threshold..					
							if (rDStdThresh <= upperStandardThreshold4InnerPeri) {
								rDStdThresh += stepsize;
								while (triedThreshes.contains(rDStdThresh)) rDStdThresh += stepsize;
								if (rDStdThresh >= upperStandardThreshold4InnerPeri+1) {
									rDStdThresh = (2*lowerStandardThreshold4InnerPeri + upperStandardThreshold4InnerPeri) / 3;
									while (triedThreshes.contains((int)Math.round(rDStdThresh))) rDStdThresh += stepsize;								
								}
								triedThreshes.add((int)Math.round(rDStdThresh));
							}						
						}	
					}
				}				
				
				//simply derive inner diameter from median wall thickness derived from the fits if unsuccessful for long..				
				if (R2 < jCFS.r2Thresh & trialCounter >= innerMaxFittingAttempts - 1 & wallThickness > 30) { //estimate inner diameter from wall thickness
					
					double medianLoc = Math.round(StatUtils.percentile(locationID,50));
					
					ArrayList<Integer> goodList = new ArrayList<Integer>();
					
					for (int l = 0 ; l < cOD.length ; l++) 
						if (Math.abs(locationID[l] - medianLoc) / medianLoc < 0.15)
							goodList.add(l);
					
					double[] wThickness = new double[goodList.size()];
					for (int l = 0 ; l < wThickness.length ; l++) {
						wThickness[l] = locationID[goodList.get(l)] - cOD[goodList.get(l)];
					}
					
					int medianWallThickness = (int)Math.round(StatUtils.percentile(wThickness, 50));
					
					int[] derivedLocID = new int[locationOD.length];
					for (int l = 0 ; l < cOD.length ; l++) {						
						derivedLocID[l] = cOD[l] + medianWallThickness;
						
						//assign x and y coords
						xID[l] = (float)x[l][derivedLocID[l]];
						yID[l] = (float)y[l][derivedLocID[l]];
					}
					
					//fit again..
					fI = fit.doRobustEllipseFit(sliceNum, xID, yID, myAngle, myIP, jCFS, -1);					
					R2 = fI.R2;
					
					//check quality of fit
					wallThickness = 0.5 * (fO.majorRadius - fI.majorRadius + fO.minorRadius - fI.minorRadius);		
					double ratioOuterInnerRadius = (fI.majorRadius + fI.minorRadius) / (fO.majorRadius + fO.minorRadius);				
					wallThicknessCV = calculateCVOfWallThickness(fI, fO);
					CVOfWallBrightness = calculateCVOfWallBrightness(myIP, fI, fO);
					
					if (wallThickness < myWindowSize + 3) R2 = 0.99 * jCFS.r2Thresh;
					if (CVOfWallBrightness > adaptedCVWallThicknessThresh) R2 = 0.99 * jCFS.r2Thresh;
					if (Double.isNaN(wallThicknessCV)) R2 = 0.99 * jCFS.r2Thresh;					
					if (wallThicknessCV > jCFS.maxCVOfWallThickness) R2 = 0.99 * jCFS.r2Thresh; // also kick out the ones with too many zeros			
					if (xID.length - countFinds > jCFS.maxNumberOfOutliers4PerimeterFits) R2 = 0.99 * jCFS.r2Thresh; // also kick out the ones with too many zeros
									
					fI.R2 = R2;
				}
			}			
		}
		
		//set inner column interface for non-Aluminium columns
		if (!jCFS.isAlu & fO.R2 > jCFS.r2Thresh) {	
			fI = fO.duplicate();
			fI.minorRadius -= jCFS.wallThickness;	
			fI.majorRadius -= jCFS.wallThickness;
		}
		
		//compare new fits with formerly found properties and decide whether the column was at this depth or not..			
		boolean colDetected = true; //first set to true
		
		if (fO.R2 < jCFS.r2Thresh) colDetected = false;
		if (colDetected & fI != null) if (fO.xCenter == 0 | fI.xCenter == 0 | fO.yCenter == 0 | fI.yCenter == 0) colDetected = false;
		if (colDetected) if (Double.isNaN(fO.R2)) colDetected = false;
		if (colDetected) if (fO.R2 < jCFS.r2Thresh) colDetected = false;	
		if (colDetected) if (myContrast < jCFS.stdThreshold) colDetected = false; //if there is no column, this will be the case..
		
		double wallThickness = 0;
		if (jCFS.isAlu & colDetected) {
			wallThickness = 0.5 * (fO.majorRadius - fI.majorRadius + fO.minorRadius - fI.minorRadius);		
			double ratioOuterInnerRadius = (fI.majorRadius + fI.minorRadius) / (fO.majorRadius + fO.minorRadius);
			double wallThicknessCV = calculateCVOfWallThickness(fI, fO);
			double CVOfWallBrightness = calculateCVOfWallBrightness(myIP, fI, fO);			
			
			if (Double.isNaN(fI.R2)) colDetected = false;
			if (fI.R2 < jCFS.r2Thresh) colDetected = false;				
			if (wallThickness < myWindowSize + 3) colDetected = false;
			if (ratioOuterInnerRadius < jCFS.ratioBetweenInnerAndOuterRadius) colDetected = false;
			if (wallThicknessCV > adaptedCVWallThicknessThresh) colDetected = false;
			if (Double.isNaN(wallThicknessCV)) colDetected = false;
			if (CVOfWallBrightness > jCFS.CVWallBrightnessThresh) colDetected = false;
		}
		
		//transfer results to output vector
		i2D.columnIsAtThisDepth = colDetected;
			
		i2D.airWallContrast = (float)relativeAirWallContrast;
	
		i2D.wallSoilStdContrastThreshold = rDStdThresh;	
		
		if (fO != null){
			i2D.xCenter = fO.xCenter;	
			i2D.yCenter = fO.yCenter;
			i2D.zCenter = sliceNum;
			i2D.outerMajorRadius = fO.majorRadius;
			i2D.outerMinorRadius = fO.minorRadius;
			i2D.theta = fO.theta;
			i2D.outerR2 = fO.R2;
		}
		
		if (fI != null){
			i2D.ixCenter = fI.xCenter;
			i2D.iyCenter = fI.yCenter;		
			i2D.innerMajorRadius = fI.majorRadius;
			i2D.innerMinorRadius = fI.minorRadius;
			i2D.itheta = fI.theta;
			i2D.innerR2 = fI.R2;
			
			i2D.wallThickness = wallThickness;
			i2D.CVOfWallThickness = calculateCVOfWallThickness(fI, fO);
			i2D.ratioBetweenOuterAndInnerRadius = calculateRatioBetweenOuterAndInnerRadius(fI, fO);
			i2D.CVOfWallBrightness = calculateCVOfWallBrightness(myIP, fI, fO);
		}
				
		//have a look at the found perimeters .. if program is run in debug mode..
		if (jCFS.debug) {
						
			if (fO.R2 > 0.98 & !colDetected) {
				IJ.showStatus("There was something funny with the initial parameter guesses..\n" 
						+ "Let me try to fit find the outlines for this depth again!");
				IJ.wait(4000);
			}
			else {
				jCFS.debug = displayFoundOutlines(sliceNum, myIP, fO, fI, xOD, yOD, xID, yID, i2D, jCFS, adaptedCVWallThicknessThresh);
			}

		}
				
		return i2D;
		
	}	
	
	public boolean displayFoundOutlines(int sliceNum, ImageProcessor myIP, FitStuff.FittedEllipse fO, FitStuff.FittedEllipse fI, 
			float[] xOD, float[] yOD, float[] xID, float[] yID, ColCoords2D i2D, MenuWaiter.ColumnFinderMenuReturn jCFS, 
			double adaptedCVWallThicknessThresh) { 
		
		boolean doDebug = true; 
		
		RoiHandler roi = new RoiHandler();
		
		ImagePlus newImg = new ImagePlus("layer " + sliceNum, myIP);
		
		//add Overlays
		PolygonRoi outerPerimeter = roi.makeRoiFromFittedEllipse(fO);			
		Overlay myO = new Overlay(outerPerimeter);
		
		if(fI != null) {
			PolygonRoi innerPerimeter = roi.makeRoiFromFittedEllipse(fI);
			myO.add(innerPerimeter);
		}
		
		ContrastEnhancer myCE = new ContrastEnhancer();
		myO.setStrokeColor(Color.RED);
		
		PointRoi outerFoundEdges = new PointRoi(xOD, yOD, xOD.length);
		myO.add(outerFoundEdges);  
		PointRoi innerFoundEdges = new PointRoi(xID, yID, xID.length);
		myO.add(innerFoundEdges);  
		newImg.setOverlay(myO);
		myCE.stretchHistogram(myIP, 0.5);
		newImg.updateAndDraw();
		newImg.show();	

		//show precision of wall findings
		GenericDialog gd = new GenericDialog("ColRecogDebugger");
										
		gd.addMessage("column found: " + i2D.columnIsAtThisDepth);
		gd.addMessage("... because ...");
		gd.addMessage("outer R2: " + String.format("%1.6f\t",(float)fO.R2) + " while criterion was set to R2 > " + String.format("%1.6f\t",(float)jCFS.r2Thresh));
		if (!jCFS.isPVC & fI != null) gd.addMessage("inner R2: " + String.format("%1.6f\t",(float)fI.R2) + " while criterion was set to R2 > " + String.format("%1.6f\t",(float)jCFS.r2Thresh));		
		if (!jCFS.isPVC & fI != null) gd.addMessage("Wall thickness: " + String.format("%3.2f\t",(float)i2D.wallThickness));
		
		if (!jCFS.isPVC & fI != null) gd.addMessage("CV of wall thickness: " + String.format("%3.2f\t",(float)i2D.CVOfWallThickness) + " while criterion was set to < " + String.format("%1.4f\t",(float)adaptedCVWallThicknessThresh));
		if (!jCFS.isPVC & fI != null) gd.addMessage("ratio between outer and inner radius was " + String.format("%1.4f\t",(float)i2D.ratioBetweenOuterAndInnerRadius) + " while criterion was set to > " + String.format("%1.4f\t",(float)jCFS.ratioBetweenInnerAndOuterRadius));				
		if (!jCFS.isPVC & fI != null) gd.addMessage("CV of wall gray values was " + String.format("%1.4f\t",(float)calculateCVOfWallBrightness(myIP, fI, fO)) + " while criterion was set to > " + String.format("%1.4f\t",(float)jCFS.CVWallBrightnessThresh));				
		
		gd.addCheckbox("Switch off Visualization mode", false);
		
		gd.hideCancelButton();
		
		gd.showDialog();
	    if (gd.wasCanceled()) return false;
	    else {
	    	doDebug = !gd.getNextBoolean();
	    }
		
	    newImg.hide();
	    newImg.flush();
	    
	    return doDebug;	    
	}
	
	
	public ColCoords3D findColumnWalls3D(ImagePlus nowTiff, ColCoords3D prelimCC, MenuWaiter.ColumnFinderMenuReturn jCFS, int[] sampleSlices) {
	// ************** routine to define the column walls (sample + cylinder)
	// ****************** - 
		
		RankFilters rF = new RankFilters();					// implements Mean, Minimum, Maximum, Variance, Median, Remove Outliers, Remove NaNs, Despeckle commands
		ImageManipulator jIM = new ImageManipulator();		// contains cutImageInXYPlane, calibrate grey values, createROIFromInnerCircle, beam dehardening, biopores, fill holes, etc, stack functions
		ImageManipulator.FloatIPCalculator iC = jIM.new FloatIPCalculator();  // class containing add & substract for image processors
		HistogramStuff hS = new HistogramStuff();			// contains functions for illumination, percentiles, etc.
		
		//init column wall coordinates
		ColCoords3D preciseCC = new ColCoords3D();
		int colHeight = nowTiff.getNSlices();
				
		double[] xCenter = new double[colHeight];			// all innerCircle varibles (inner = sample, outer = cylinder) as arrays with length of image "height"
		double[] yCenter = new double[colHeight];
		double[] zCenter = new double[colHeight];
		double[] outerMinorRadius = new double[colHeight];
		double[] outerMajorRadius = new double[colHeight];
		double[] deltaMajorMinor =  new double[colHeight];
		double[] theta = new double[colHeight];		
		double[] outerR2 = new double[colHeight];
		double[] ixCenter = new double[colHeight];
		double[] iyCenter = new double[colHeight];			
		double[] innerMinorRadius = new double[colHeight];
		double[] innerMajorRadius = new double[colHeight];				
		double[] itheta = new double[colHeight];		
		double[] innerR2 = new double[colHeight];		
		boolean[] columnIsAtThisDepth = new boolean[colHeight];		
		double[] wallThickness = new double[colHeight];		
		
		double[] airWallContrast = new double[colHeight];
		double[] absoluteWallgrayValue = new double[colHeight];
	
		double[] wallSoilStdContrastThreshold = new double[colHeight];	
		
		//best ColCoords so far
		ColCoords2D bestCoords = median3DCoords(prelimCC);			// take the median coordinates of the perliminary Column outlines for a first best guess
		
		//init wall finder parameters
		bestCoords.airWallContrast = jCFS.airWallContrast;			// take the air / wall contrast that was entered by the user in the menu (jCFS)
		bestCoords.wallSoilStdContrastThreshold = jCFS.wallSoilStdContrastThreshold;	// take the contrast of soil sample vs. cylinder wall entered by user (jCFS)
					
		//probe illumination level for 3-D column
		HistogramStuff.IlluminationInfo jII = hS.new IlluminationInfo(); 	// IlluminationInfo contains ints with different quantiles, mean, median, ect. (of grey values I guess)
		for (int i = 0 ; i < colHeight ; i++) {  							// loop through the 3Dtiff (ImagePlus) nowTiff and probe Illumination Level
			nowTiff.setPosition(i+1);
			ImageProcessor nowIP = nowTiff.getProcessor();		
			ImageProcessor myIP = nowIP.duplicate();
			jII = probeIlluminationLevel(myIP, bestCoords, true);			// probeIlluminationLevel = ????????????????
		} // end colHeight
		
		//re-find the column's outer wall		
		for (int i = 0 ; i < colHeight ; i++) {  
			
			IJ.showStatus("Searching column outlines in slice " + (sampleSlices[i] + 1) + " ...");
							
			//set tiff to the correct position and get Processor etc..
			nowTiff.setPosition(i+1);
			ImageProcessor nowIP = nowTiff.getProcessor();		
			ImageProcessor myIP = nowIP.duplicate();
			int[] myHist = myIP.getHistogram();						// get the single slice (image processor), duplicate it and get histogram
			
			//ImagePlus tesImg = new ImagePlus("", myIP); tesImg.show();
			
			//fill zero values with background if image is not already calibrated...
			if (!jCFS.isAlreadyNormalized & myHist[0] > 0) {		// user can say whether image is already calibrated (user settings in jCFS)
			
				ImageProcessor zeroIP = myIP.duplicate();			// if not, duplicate and do so now
				zeroIP.scaleAndSetThreshold(0, 1, 1);
				zeroIP.threshold(2000);
				zeroIP.invert();
				//zeroIP.dilate();
				
				//ImagePlus tesImg = new ImagePlus("", zeroIP); tesImg.show();
				
				ArrayList<Integer> fringeVoxels = new ArrayList<Integer>();
				for (int x = 0 ; x < zeroIP.getWidth() ; x++) {
					for (int y = 0 ; y < zeroIP.getHeight() ; y++) {
						int oldP = (int)Math.round(myIP.getPixelValue(x, y));
						int newP = (int)Math.round(zeroIP.getPixelValue(x, y));
						if (oldP > 0 & newP > 0) fringeVoxels.add((int)Math.round(myIP.getPixelValue(x, y)));					
					}
				}
				double[] fillValue = new double[fringeVoxels.size()];
				for (int j = 0 ; j < fringeVoxels.size() ; j++) fillValue[j] = fringeVoxels.get(j);
				int meanFillValue = (int)Math.round(StatUtils.mean(fillValue));
				myIP.min(meanFillValue);
			}
			
			//tesImg = new ImagePlus("", myIP); tesImg.show();
			 			
			//apply an additional median filter		
			if (jCFS.medianFilter2D > 0) rF.rank(myIP, jCFS.medianFilter2D, RankFilters.MEDIAN);
			
			//apply an additional unsharp mask			
			if (jCFS.applyUmask) {
				double sigma = 2;
				float weight = 0.7f;				
				FloatProcessor fIP = myIP.convertToFloatProcessor();
				ImageProcessor bIP = fIP.duplicate();
				bIP.blurGaussian(sigma);				
				ImageProcessor dIP = iC.subtract(fIP, bIP);				
				dIP.multiply(weight);				
				ImageProcessor outIP = iC.add(dIP, fIP);
				myIP = outIP.convertToShortProcessor(false);
			}
			
			//find do gray-scale pre-normalization...
			
			
			//ImagePlus tesImg = new ImagePlus("", myIP); tesImg.show();
			
                    
			//find the wall;
			ColCoords2D i2D = findColumnWalls2D(sampleSlices[i], myIP, bestCoords, jCFS);	// for the sample sclies find the column wall outlines, returns 2D elliptical coordinates
			
			//catch bug that makes first fit worse than other due to bad initial parameter estimation
			if (!i2D.columnIsAtThisDepth & i2D.outerR2 > 0.98) {
				bestCoords = i2D;
				i2D = findColumnWalls2D(sampleSlices[i], myIP, bestCoords, jCFS);
			}
											
			//transfer fitting results to vectors for export
			columnIsAtThisDepth[i] = i2D.columnIsAtThisDepth;
			if (i2D.columnIsAtThisDepth) bestCoords = i2D;
			
			xCenter[i] = i2D.xCenter;
			yCenter[i] = i2D.yCenter;
			zCenter[i] = i2D.zCenter;
			outerMinorRadius[i] = i2D.outerMinorRadius;
			outerMajorRadius[i] = i2D.outerMajorRadius;	
			deltaMajorMinor[i] = outerMajorRadius[i] - outerMinorRadius[i];			
			theta[i] = i2D.theta;		
			outerR2[i] = i2D.outerR2;			
			
			//transfer results to array
			ixCenter[i] = i2D.ixCenter;
			iyCenter[i] = i2D.iyCenter;			
			innerMinorRadius[i] = i2D.innerMinorRadius;
			innerMajorRadius[i] = i2D.innerMajorRadius;				
			itheta[i] = i2D.itheta;		
			innerR2[i] = i2D.innerR2;	
			
			wallThickness[i] = i2D.wallThickness;
			
			airWallContrast[i] = i2D.airWallContrast;
			absoluteWallgrayValue[i] = i2D.absoluteWallgrayValue;
			
			wallSoilStdContrastThreshold[i] = i2D.wallSoilStdContrastThreshold;
		
		} // loop through all image slices
		
		//assign results to colCon
		preciseCC.xmid = xCenter;		
		preciseCC.ymid = yCenter;
		preciseCC.zmid = zCenter;
		preciseCC.outerMinorRadius = outerMinorRadius;
		preciseCC.outerMajorRadius = outerMajorRadius;
		preciseCC.deltaMajorMinorRadius = deltaMajorMinor;
		preciseCC.theta = theta;		
		preciseCC.outerR2 = outerR2;
		preciseCC.ixmid = ixCenter;
		preciseCC.iymid = iyCenter;			
		preciseCC.innerMinorRadius = innerMinorRadius;
		preciseCC.innerMajorRadius = innerMajorRadius;				
		preciseCC.itheta = itheta;
		preciseCC.innerR2 = innerR2;
		preciseCC.columnIsAtThisDepth = columnIsAtThisDepth;
		preciseCC.wallThickness = wallThickness;
		
		preciseCC.airWallContrast = airWallContrast;
		preciseCC.absoluteWallgrayValue = absoluteWallgrayValue;
		
		preciseCC.wallSoilStdContrastThreshold = wallSoilStdContrastThreshold;
				
		if (jCFS.putColumnStraight) {
		
			int cc = 0;
			for (int i = 0 ; i < preciseCC.outerR2.length ; i++) if (preciseCC.outerR2[i] > jCFS.r2Thresh) cc++;
			
			//filter out bad R2s
			double[] xFiltered = new double[cc];
			double[] yFiltered = new double[cc];
			double[] zFiltered = new double[cc];
			cc = 0;
			for (int i = 0 ; i < preciseCC.outerR2.length ; i++) {
				if (preciseCC.outerR2[i] > jCFS.r2Thresh) {
					xFiltered[cc]=xCenter[i];
					yFiltered[cc]=yCenter[i];
					zFiltered[cc]=zCenter[i];			
					cc++;
				}
			}
			
			//find column Tilt!		
			double[] tilt = findColumnTilt(xFiltered, yFiltered, zFiltered);
			preciseCC.tiltInXZ = tilt[0];
			preciseCC.tiltInYZ = tilt[1];
			preciseCC.tiltTotal = tilt[2];
		}
		
		return preciseCC;
	}
	
	public ColCoords2D median3DCoords(ColCoords3D thD) {
		
		ColCoords2D twD = new ColCoords2D();
				
		twD.xCenter = StatUtils.percentile(thD.xmid ,50);
		twD.yCenter = StatUtils.percentile(thD.ymid ,50);
		twD.outerMajorRadius = StatUtils.percentile(thD.outerMajorRadius ,50);
		twD.outerMinorRadius = StatUtils.percentile(thD.outerMinorRadius ,50);
		twD.theta = StatUtils.percentile(thD.theta ,50);
			
		return twD;
	}
	

	
	public ColumnContainer findTopAndBottomOfSoilColumn(ColumnContainer colCon) {
		
		colCon.preciseCC.topOfColumn = 0;		
		colCon.preciseCC.bottomOfColumn = 0;
		if (colCon.jCFS.try2FindColumnTopAndBottom) {
			for (int i = 0 ; i < colCon.nowTiff.getNSlices() ; i++) {
				if ((colCon.preciseCC.topOfColumn == 0) && colCon.preciseCC.columnIsAtThisDepth[i]) colCon.preciseCC.topOfColumn = i;
				if ((colCon.preciseCC.topOfColumn > 0) && colCon.preciseCC.columnIsAtThisDepth[i]) colCon.preciseCC.bottomOfColumn = i;
			}
		}
		else {
			colCon.preciseCC.topOfColumn =	colCon.jCFS.topOfColumn; 
			if (colCon.jCFS.bottomOfColumn == 1) colCon.preciseCC.bottomOfColumn = colCon.nowTiff.getNSlices();
			else colCon.preciseCC.bottomOfColumn = colCon.jCFS.bottomOfColumn;
		}		
		colCon.preciseCC.heightOfColumn = colCon.preciseCC.bottomOfColumn - colCon.preciseCC.topOfColumn;
		
		return colCon;
	}
	
	public ColumnContainer flagDoubtworthyResults(ColumnContainer colCon) {
		
		//calc medians of column diameters
		double outer = 0.5 * (StatUtils.percentile(colCon.preciseCC.outerMajorRadius, 50) + StatUtils.percentile(colCon.preciseCC.outerMinorRadius, 50));		
		double inner = 0.5 * (StatUtils.percentile(colCon.preciseCC.innerMajorRadius, 50) + StatUtils.percentile(colCon.preciseCC.innerMinorRadius, 50));	
		double wall = StatUtils.percentile(colCon.preciseCC.wallThickness, 50);	
				
		//set R2s for columnNotInThisDepth to 0;
		for (int z = 0 ; z < colCon.preciseCC.outerR2.length ; z++) {
			double nowOuter = 0.5 * (colCon.preciseCC.outerMajorRadius[z] + colCon.preciseCC.outerMinorRadius[z]);
			double nowInner = 0.5 * (colCon.preciseCC.innerMajorRadius[z] + colCon.preciseCC.innerMinorRadius[z]);
			double nowWall = colCon.preciseCC.wallThickness[z];
								
			if (colCon.preciseCC.outerR2[z] == 0) colCon.preciseCC.columnIsAtThisDepth[z] = false;
									
			if (!colCon.preciseCC.columnIsAtThisDepth[z] | nowWall > 1.1 * wall | nowOuter > 1.05 * outer | nowInner < 0.95 * inner) {
				colCon.preciseCC.innerR2[z] = 0;
				colCon.preciseCC.outerR2[z] = 0;
				colCon.preciseCC.columnIsAtThisDepth[z] = false;
			}
		}
		
		//kick-out layers with sudden jump in diameter at the end of the column 
		
		return colCon;
	}
		
	public ColCoords3D initColCords3D(int size) {
				
		ColCoords3D outCoords = new ColCoords3D();
		
		//init outCoords even more		
		outCoords.xmid = new double[size];			//x midpoint
		outCoords.ymid = new double[size];													//y midpoint
		outCoords.zmid = new double[size];				//y midpoint
		
		outCoords.ixmid = new double[size];				//x midpoint (inner circle)
		outCoords.iymid = new double[size];				//y midpoint (inner circle)
		
		outCoords.outerMajorRadius = new double[size];	
		outCoords.innerMajorRadius = new double[size];	
		outCoords.outerMinorRadius = new double[size];	
		outCoords.innerMinorRadius = new double[size];	
		outCoords.deltaMajorMinorRadius = new double[size];	
		outCoords.wallThickness = new double[size];	
		outCoords.theta = new double[size];	  //angle of major ellipse axis
		outCoords.itheta = new double[size];	  //angle of major ellipse axis (inner circle)
		
		outCoords.outerR2 = new double[size];	
		outCoords.innerR2 = new double[size];	
		
		outCoords.columnIsAtThisDepth = new boolean[size];	
						
		outCoords.airWallContrast = new double[size];	
		outCoords.absoluteWallgrayValue = new double[size];	
		outCoords.wallSoilStdContrastThreshold = new double[size];			
		
		//assign topSlice
		outCoords.xmid = new double[size];			//x midpoint
		outCoords.ymid = new double[size];													//y midpoint
		outCoords.zmid = new double[size];				//y midpoint
		
		outCoords.ixmid = new double[size];				//x midpoint (inner circle)
		outCoords.iymid = new double[size];				//y midpoint (inner circle)
		
		outCoords.outerMajorRadius = new double[size];	
		outCoords.innerMajorRadius = new double[size];	
		outCoords.outerMinorRadius = new double[size];	
		outCoords.innerMinorRadius = new double[size];	
		outCoords.deltaMajorMinorRadius = new double[size];	
		outCoords.wallThickness = new double[size];	
		outCoords.theta = new double[size];	  //angle of major ellipse axis
		outCoords.itheta = new double[size];	  //angle of major ellipse axis (inner circle)
		
		outCoords.outerR2 = new double[size];	
		outCoords.innerR2 = new double[size];	
		
		outCoords.columnIsAtThisDepth = new boolean[size];	
						
		outCoords.airWallContrast = new double[size];	
		outCoords.absoluteWallgrayValue = new double[size];	
		outCoords.wallSoilStdContrastThreshold = new double[size];
		
		return outCoords;
		
	}
	
	public ColCoords3D plug2DInto3DCoords(int i, ColCoords2D inCoords, ColCoords3D outCoords){
		
		outCoords.xmid[i] = inCoords.xCenter;				//x midpoint
		outCoords.ymid[i] = inCoords.yCenter;													//y midpoint
		outCoords.zmid[i] = inCoords.zCenter;				//y midpoint
		
		outCoords.ixmid[i] = inCoords.ixCenter;				//x midpoint (inner circle)
		outCoords.iymid[i] = inCoords.iyCenter;				//y midpoint (inner circle)
		
		outCoords.outerMajorRadius[i] = inCoords.outerMajorRadius;	
		outCoords.innerMajorRadius[i] = inCoords.innerMajorRadius;	
		outCoords.outerMinorRadius[i] = inCoords.outerMinorRadius;	
		outCoords.innerMinorRadius[i] = inCoords.innerMinorRadius;	
		
		outCoords.wallThickness[i] = inCoords.wallThickness;	
		outCoords.theta[i] = inCoords.theta;	  //angle of major ellipse axis
		outCoords.itheta[i] = inCoords.itheta;	  //angle of major ellipse axis (inner circle)
		
		outCoords.outerR2[i] = inCoords.outerR2;	
		outCoords.innerR2[i] = inCoords.innerR2;	
		
		outCoords.columnIsAtThisDepth[i] = inCoords.columnIsAtThisDepth;	
						
		outCoords.airWallContrast[i] = inCoords.airWallContrast;	
		outCoords.absoluteWallgrayValue[i] = inCoords.absoluteWallgrayValue;	
		outCoords.wallSoilStdContrastThreshold[i] = inCoords.wallSoilStdContrastThreshold;
		
		return outCoords;
	}
	
	public ColCoords3D plug3DInto3DCoords(int[] placings, int[] goodSlices, ColCoords3D inCoords, ColCoords3D outCoords){
					
		for (int j = 0 ; j < goodSlices.length ; j++) {
			
			//int nowZ = (int)Math.round(inCoords.zmid[goodSlices[j - i]] - outCoords.zmid[0] );
			
			outCoords.xmid[placings[j]] = inCoords.xmid[goodSlices[j]];			//x midpoint
			outCoords.ymid[placings[j]] = inCoords.ymid[goodSlices[j]];			//y midpoint
			outCoords.zmid[placings[j]] = inCoords.zmid[goodSlices[j]];				//y midpoint
			
			outCoords.ixmid[placings[j]] = inCoords.ixmid[goodSlices[j]];				//x midpoint (inner circle)
			outCoords.iymid[placings[j]] = inCoords.iymid[goodSlices[j]];				//y midpoint (inner circle)
			
			outCoords.outerMajorRadius[placings[j]] = inCoords.outerMajorRadius[goodSlices[j]];	
			outCoords.innerMajorRadius[placings[j]] = inCoords.innerMajorRadius[goodSlices[j]];	
			outCoords.outerMinorRadius[placings[j]] = inCoords.outerMinorRadius[goodSlices[j]];	
			outCoords.innerMinorRadius[placings[j]] = inCoords.innerMinorRadius[goodSlices[j]];	
			
			outCoords.wallThickness[placings[j]] = inCoords.wallThickness[goodSlices[j]];	
			outCoords.theta[placings[j]] = inCoords.theta[goodSlices[j]];	  //angle of major ellipse axis
			outCoords.itheta[placings[j]] = inCoords.itheta[goodSlices[j]];	  //angle of major ellipse axis (inner circle)
			
			outCoords.outerR2[placings[j]] = inCoords.outerR2[goodSlices[j]];	
			outCoords.innerR2[placings[j]] = inCoords.innerR2[goodSlices[j]];	
			
			outCoords.columnIsAtThisDepth[placings[j]] = inCoords.columnIsAtThisDepth[goodSlices[j]];	
							
			outCoords.airWallContrast[placings[j]] = inCoords.airWallContrast[goodSlices[j]];	
			outCoords.absoluteWallgrayValue[placings[j]] = inCoords.absoluteWallgrayValue[goodSlices[j]];	
			outCoords.wallSoilStdContrastThreshold[placings[j]] = inCoords.wallSoilStdContrastThreshold[goodSlices[j]];
		}
		
		return outCoords;
	}
	
	public ColCoords3D findColumnsBevel(InputOutput.MyFileCollection mFC, ColCoords3D samCoords, MenuWaiter.ColumnFinderMenuReturn jCFS) {
				
		FitStuff fit = new FitStuff();
		RankFilters rF=new RankFilters();
		
		//find entries with outer perimeter data 
		ArrayList<Integer> validDataAtLocationZ = new ArrayList<Integer>();
		for (int i = 0 ; i < samCoords.zmid.length ; i++) {
			if (samCoords.xmid[i] > 0 & samCoords.ymid[i] > 0 & samCoords.ixmid[i] > 0 & samCoords.iymid[i] > 0) {
				validDataAtLocationZ.add(i);
			}
		}
		
		//define the median wall thickness above the bevel..
		ArrayList<Double> outerRadiusWithoutBevel = new ArrayList<Double>();
		for (int i = 0 ; i < validDataAtLocationZ.size() ; i++) {
			int j = validDataAtLocationZ.get(i);
			if (samCoords.zmid[j] < samCoords.topOfColumn + 0.5 * samCoords.heightOfColumn) {
				outerRadiusWithoutBevel.add(0.5 * samCoords.outerMajorRadius[j] + 0.5 * samCoords.outerMinorRadius[j]);
			}
		}
		double[] outerRadius = new double[outerRadiusWithoutBevel.size()]; 
		for (int i = 0 ; i < outerRadius.length ; i++) outerRadius[i] = outerRadiusWithoutBevel.get(i);
		double medianOuterRadius = StatUtils.percentile(outerRadius, 50);
		
		//find all parts of the column with outer radii smaller than xxx * median WallThickness
		ArrayList<Double> bevelRadius = new ArrayList<Double>();
		ArrayList<Double> z = new ArrayList<Double>();
		for (int i = 0 ; i < validDataAtLocationZ.size() ; i++) {
			int j = validDataAtLocationZ.get(i);
			double outerRadiusNow = 0.5 * samCoords.outerMajorRadius[j] + 0.5 * samCoords.outerMinorRadius[j];
			if (outerRadiusNow < 0.975 * medianOuterRadius) {
				bevelRadius.add(outerRadiusNow);
				z.add(samCoords.zmid[j]);
			}
		}
		
		//fit straight to bevel 
		double[] Z = new double[z.size()];for (int i = 0 ; i < Z.length ; i++) Z[i] = z.get(i);
		double[] bw = new double[bevelRadius.size()];for (int i = 0 ; i < bw.length ; i++) bw[i] = bevelRadius.get(i);
		samCoords.bevelFunction = fit.fitLinearFunction(Z, bw);
		
		//evaluate start of bevel
		if (bevelRadius.size() > 2) {
			
			double nowColumnWall = bevelRadius.get(0);		
			double zz = z.get(0);
			while (nowColumnWall < medianOuterRadius ) {
				zz--;
				nowColumnWall = fit.evalLinearFunction(samCoords.bevelFunction, zz);
			}
			samCoords.bevelStartsAt = (int)Math.round(zz++);
			
			//kick out values that have increasing wall-thickness gradient after the bevel has been found.. no, better do not do this
			/*boolean theFormerGradientWasPositiveYuck = false;
			boolean theGradientWasPositiveYuck = false;	
			double ketchupGap = 0;
			
			for (int i = samCoords.bevelStartsAt + 1; i < samCoords.outerR2.length ; i++) {
								
				if (theGradientWasPositiveYuck) theFormerGradientWasPositiveYuck = true;		
				else theFormerGradientWasPositiveYuck = false;
				theGradientWasPositiveYuck = false;
				
				double outerRadiusNow = 0.5 * samCoords.outerMajorRadius[i] + 0.5 * samCoords.outerMinorRadius[i];
				double outerRadiusThen = 0.5 * samCoords.outerMajorRadius[i-1] + 0.5 * samCoords.outerMinorRadius[i-1];
				double outerRadiusGradient = outerRadiusNow - outerRadiusThen; //should be negative once the bevel has been reached
				
				if (outerRadiusGradient > 0) { 
					theGradientWasPositiveYuck = true;
					ketchupGap = outerRadiusNow - outerRadiusThen;
				}				
				if (theFormerGradientWasPositiveYuck) if (outerRadiusGradient < - ketchupGap) theGradientWasPositiveYuck = false;
				else theGradientWasPositiveYuck = true;
				
				if (theGradientWasPositiveYuck) {
					samCoords.innerR2[i] = 0;
					samCoords.outerR2[i] = 0;
					samCoords.columnIsAtThisDepth[i] = false;
				}
				else {
					ketchupGap = 0;
					if (jCFS.isPVC) {
						double gap = medianOuterRadius - fit.evalLinearFunction(samCoords.bevelFunction, i);				
						samCoords.innerMajorRadius[i] += gap;
						samCoords.innerMinorRadius[i] += gap;
					}
				}
				
			}*/
		}
		
		if (samCoords.bevelStartsAt > samCoords.topOfColumn + 3 & samCoords.bevelStartsAt < samCoords.bottomOfColumn - 3) {
			//sample x images from top and check when the gray values stabilize..		
			int validationRange = 5; // check the topmost 10 slices
			int[] valiSeries = new int[validationRange];
			int valiStart = (int)samCoords.bevelStartsAt - validationRange / 2; 
			for (int i = valiStart ; i < valiStart + validationRange ; i++) valiSeries[i - valiStart] = i;
			
			InputOutput jIO = new InputOutput();
			ImagePlus valiTiff = new ImagePlus();
			if (mFC.imageHasBeenLoaded) {
				ImageStack newStack = new ImageStack(mFC.nowTiff.getWidth(), mFC.nowTiff.getHeight());
		    	for (int j = 0 ; j < valiSeries.length ; j++) {
		    		
		    		mFC.nowTiff.setPosition(valiSeries[j]);	    		
		    		ImageProcessor nowIP = mFC.nowTiff.getProcessor();
		    		newStack.addSlice(nowIP);    		
		    		
		    	}
		    	ImagePlus outTiff = new ImagePlus("", newStack);
		    	valiTiff = outTiff;
			}
			else valiTiff = jIO.openTiff3DSomeSlices(mFC.nowTiffPath, mFC.nowWidth, mFC.nowHeight, valiSeries);
			
			int closestZ = 0;
			for (int i = 0 ; i < validDataAtLocationZ.size() ; i++) if (validDataAtLocationZ.get(i) + samCoords.topOfColumn < valiStart & samCoords.columnIsAtThisDepth[validDataAtLocationZ.get(i)]) {
				closestZ = validDataAtLocationZ.get(i);
				if (i > valiStart) break;		
			}
			
			ColCoords2D bevCoords = new ColCoords2D();
			
			//jCFS.debug = true;
			
			for (int i = 0 ; i < valiTiff.getNSlices() ; i++) {
				
				IJ.showStatus("Checking for bevel find in slice " + (valiSeries[i] + 1) + " ...");
				
				valiTiff.setPosition(i + 1);
				ImageProcessor myIP = valiTiff.getProcessor();
				rF.rank(myIP, jCFS.medianFilter2D, RankFilters.MEDIAN);
				ColCoords2D prelimCoords = extract2DCoordsFrom3DCoords(samCoords, closestZ);
				bevCoords = findColumnWalls2D(valiSeries[i], myIP, prelimCoords, jCFS);
				
				//if column was found, add it to the big Coords
				if (bevCoords.columnIsAtThisDepth) samCoords = plug2DInto3DCoords(valiSeries[i] - samCoords.topOfColumn, bevCoords, samCoords);
			}
			
			//init correction 4 inner diameter below bevel for PVC columns
			ArrayList<Integer> newValidDataAtLocationZ = new ArrayList<Integer>();
			for (int i = 0 ; i < samCoords.zmid.length ; i++) {
				if (samCoords.xmid[i] > 0 & samCoords.ymid[i] > 0 & samCoords.ixmid[i] > 0 & samCoords.iymid[i] > 0) {
					newValidDataAtLocationZ.add(i);
				}
			}
					
			//find the layers below the bevel		
			ArrayList<Integer> toBeCorrected = new ArrayList<Integer>();
			for (int i = 0 ; i < newValidDataAtLocationZ.size() ; i++) {
				int j = newValidDataAtLocationZ.get(i);
				if (samCoords.zmid[j] > samCoords.bevelStartsAt - 5) {  // include some layer as buffer
					toBeCorrected.add(j);
				}
			}
			
			//correct inner diameters for bevel
			if (!jCFS.isAlu) {
				for (int i = 0 ; i < toBeCorrected.size() ; i++) {		
					double outerRadiusNow = 0.5 * samCoords.outerMajorRadius[toBeCorrected.get(i)] + 0.5 * samCoords.outerMinorRadius[toBeCorrected.get(i)];
					double corrFactor = medianOuterRadius - outerRadiusNow;
					samCoords.innerMajorRadius[toBeCorrected.get(i)] += corrFactor;
					samCoords.innerMinorRadius[toBeCorrected.get(i)] += corrFactor;
				}
			}
		}
		
		return samCoords;
	}
	
	public ColCoords3D imputeMissingLayers(ColCoords3D samCoords, MenuWaiter.ColumnFinderMenuReturn jCFS) {
				
		////////////////////////////////////
		// 2016/04/05: still missing: imputation does not work for ellipsoid columns.. bevel finder would need to be adjusted for that too..
		////////////////////////////////////
		
		//init mendedCC
		ColCoords3D mendedCoords = samCoords;
		
		//find layers with complete values
		ArrayList<Integer> hasOuterPerimeter = new ArrayList<Integer>();
		ArrayList<Integer> hasInnerPerimeter = new ArrayList<Integer>();
		
		for (int z = 0 ; z < samCoords.zmid.length ; z++) {
			if (samCoords.outerR2[z] >= jCFS.r2Thresh) hasOuterPerimeter.add(z);
			if (jCFS.isAlu & samCoords.innerR2[z] >= jCFS.r2Thresh) hasInnerPerimeter.add(z);
			else if (jCFS.isPVC & samCoords.outerR2[z] >= jCFS.r2Thresh) hasInnerPerimeter.add(z);
		} // 0 to zmid.length
		
		//impute the missing values for inner perimeter...
		for (int z = 0 ; z < samCoords.zmid.length ; z++) {
			
			if (!hasInnerPerimeter.contains(z)) {			
			
				mendedCoords.numberOfImputedLayers++;
				
				//find existing entries above and below the missing v1alue
				int data4Imputation = -1;
				for (int j = 0; j < hasInnerPerimeter.size() ; j++) {
					if (hasInnerPerimeter.get(j) > z) {
						data4Imputation = j;
						break;					
					}
				}
									
				if (data4Imputation > 0){
					
					//find available data above and below
					int jabove = hasInnerPerimeter.get(data4Imputation - 1);
					int jbelow = hasInnerPerimeter.get(data4Imputation);
					
					//calculate weights according to distance to available finds
					double totalDistance = jbelow - jabove;
					double wAbove = (double)(jbelow - z) / totalDistance;
					double wBelow = (double)(z - jabove) / totalDistance;
					
					mendedCoords.ixmid[z] = wAbove * samCoords.ixmid[jabove] + wBelow * samCoords.ixmid[jbelow];
					mendedCoords.iymid[z] = wAbove * samCoords.iymid[jabove] + wBelow * samCoords.iymid[jbelow];
					mendedCoords.innerMinorRadius[z] = wAbove * samCoords.innerMinorRadius[jabove] + wBelow * samCoords.innerMinorRadius[jbelow];
					mendedCoords.innerMajorRadius[z] = wAbove * samCoords.innerMajorRadius[jabove] + wBelow * samCoords.innerMajorRadius[jbelow];
					mendedCoords.itheta[z] = wAbove * samCoords.itheta[jabove] + wBelow * samCoords.itheta[jbelow];
					
				}
				
				if (data4Imputation == 0) {
					
					mendedCoords.ixmid[z] = samCoords.ixmid[data4Imputation];
					mendedCoords.iymid[z] = samCoords.iymid[data4Imputation];						
					mendedCoords.innerMinorRadius[z] = samCoords.innerMinorRadius[data4Imputation];
					mendedCoords.innerMajorRadius[z] = samCoords.innerMajorRadius[data4Imputation];
					mendedCoords.itheta[z] = samCoords.itheta[data4Imputation];	
										
				}
				
				if (data4Imputation < 0) {
					
					mendedCoords.ixmid[z] = samCoords.ixmid[hasInnerPerimeter.get(hasInnerPerimeter.size() - 1)];
					mendedCoords.iymid[z] = samCoords.iymid[hasInnerPerimeter.get(hasInnerPerimeter.size() - 1)];
					mendedCoords.innerMinorRadius[z] = samCoords.innerMinorRadius[hasInnerPerimeter.get(hasInnerPerimeter.size() - 1)];
					mendedCoords.innerMajorRadius[z] = samCoords.innerMajorRadius[hasInnerPerimeter.get(hasInnerPerimeter.size() - 1)];
					mendedCoords.itheta[z] = samCoords.itheta[hasInnerPerimeter.get(hasInnerPerimeter.size() - 1)];
					
				}
			}	
		}
		
		//impute the missing values for outer perimeter...
		for (int z = 0 ; z < samCoords.zmid.length ; z++) {
			
			if (!hasOuterPerimeter.contains(z)) {			
			
				//find existing entries above and below the missing value
				int data4Imputation = -1;
				for (int j = 0 ; j < hasOuterPerimeter.size() ; j++) {
					if (hasOuterPerimeter.get(j) > z) {
						data4Imputation = j;
						break;					
					}
				}
				
				//find available data above	
				if (data4Imputation > 0){
					
					//find available data above and below
					int jabove = hasOuterPerimeter.get(data4Imputation - 1);
					int jbelow = hasOuterPerimeter.get(data4Imputation);
					
					//calculate weights according to distance to available finds
					double totalDistance = jbelow - jabove;
					double wAbove = (double)(jbelow - z) / totalDistance;
					double wBelow = (double)(z - jabove) / totalDistance;
					
					mendedCoords.xmid[z] = wAbove * samCoords.xmid[jabove] + wBelow * samCoords.xmid[jbelow];
					mendedCoords.ymid[z] = wAbove * samCoords.ymid[jabove] + wBelow * samCoords.ymid[jbelow];
					mendedCoords.zmid[z] = wAbove * samCoords.zmid[jabove] + wBelow * samCoords.zmid[jbelow];
					
					mendedCoords.theta[z] = wAbove * samCoords.theta[jabove] + wBelow * samCoords.theta[jbelow];
					
					//take care of bevel.. not necessary in the new version
					/*if (z < samCoords.bevelStartsAt - samCoords.topOfColumn) {
						mendedCoords.outerMinorRadius[z] = fit.evalLinearFunction(samCoords.bevelFunction, z);
						mendedCoords.outerMajorRadius[z] = fit.evalLinearFunction(samCoords.bevelFunction, z);
					}
					else {*/
						mendedCoords.outerMinorRadius[z] = wAbove * samCoords.outerMinorRadius[jabove] + wBelow * samCoords.outerMinorRadius[jbelow];
						mendedCoords.outerMajorRadius[z] = wAbove * samCoords.outerMajorRadius[jabove] + wBelow * samCoords.outerMajorRadius[jbelow];
					//}
					
				}
				if (data4Imputation == 0) {
					
					mendedCoords.xmid[z] = samCoords.xmid[data4Imputation];
					mendedCoords.ymid[z] = samCoords.ymid[data4Imputation];
					mendedCoords.zmid[z] = samCoords.zmid[data4Imputation];
					
					mendedCoords.theta[z] = samCoords.theta[data4Imputation];
					
				/*	//take care of bevel.. not necessary in the new version
					if (z > samCoords.bevelStartsAt) {
						mendedCoords.outerMinorRadius[z] = fit.evalLinearFunction(samCoords.bevelFunction, z);
						mendedCoords.outerMajorRadius[z] = fit.evalLinearFunction(samCoords.bevelFunction, z);
					}
					else {*/
						mendedCoords.outerMinorRadius[z] = samCoords.outerMinorRadius[data4Imputation];
						mendedCoords.outerMajorRadius[z] = samCoords.outerMajorRadius[data4Imputation];						
					//}
					
					
				}
				if (data4Imputation < 0) {
					
					mendedCoords.xmid[z] = samCoords.xmid[hasOuterPerimeter.get(hasOuterPerimeter.size() - 1)];
					mendedCoords.ymid[z] = samCoords.ymid[hasOuterPerimeter.get(hasOuterPerimeter.size() - 1)];	
					mendedCoords.zmid[z] = samCoords.zmid[hasOuterPerimeter.get(hasOuterPerimeter.size() - 1)];	
					
					mendedCoords.theta[z] = samCoords.theta[hasOuterPerimeter.get(hasOuterPerimeter.size() - 1)];
					
					/*//take care of bevel.. not necessary in the new version
					if (z > samCoords.bevelStartsAt) {
						mendedCoords.outerMinorRadius[z] = fit.evalLinearFunction(samCoords.bevelFunction, z);
						mendedCoords.outerMajorRadius[z] = fit.evalLinearFunction(samCoords.bevelFunction, z);
					}*/
					//else{
						mendedCoords.outerMinorRadius[z] = samCoords.outerMinorRadius[hasOuterPerimeter.get(hasOuterPerimeter.size() - 1)];
						mendedCoords.outerMajorRadius[z] = samCoords.outerMajorRadius[hasOuterPerimeter.get(hasOuterPerimeter.size() - 1)];
					//}
				
				}
			}	
		}
		
		//re-calculate wall thickness
		for (int z = 0 ; z < samCoords.zmid.length ; z++) {
			mendedCoords.wallThickness[z] = 0.5 * (mendedCoords.outerMajorRadius[z] - mendedCoords.innerMajorRadius[z] + mendedCoords.outerMinorRadius[z] - mendedCoords.innerMinorRadius[z]);
		}
		
		return mendedCoords;
	}
	
	public ColumnContainer rubberBandDetection(ColumnContainer colCon) {
		
		//check if there are layers in between where the wall was not properly detected (because of dirt or rubber band)
		double[] checker = new double[colCon.preciseCC.zmid.length];
		int diffCheckWindowSize = 1;
			
		double[] xf2 = makeChecker(colCon.preciseCC.xmid, diffCheckWindowSize);
		double[] yf2 = makeChecker(colCon.preciseCC.ymid, diffCheckWindowSize);
		double[] mif2 = makeChecker(colCon.preciseCC.outerMinorRadius, diffCheckWindowSize);
		double[] maf2 = makeChecker(colCon.preciseCC.outerMajorRadius, diffCheckWindowSize);		
		for (int i = 0 ; i < checker.length ; i++) checker[i] = Math.sqrt(xf2[i] + yf2[i] + mif2[i] + maf2[i]) / 4;				
		
		//find high and low probability outlines
		int filterSize = 16;
		int[] goodones = findSmallCheckers(checker, filterSize);
		
		//kick out values above and below column
		//for (int i = 0 ; i < goodones.length ; i++) if (i < colCon.preciseCC.topOfColumn | i > colCon.preciseCC.bottomOfColumn) goodones[i] = 0;
		
		//apply correction		
	    int MaximumGoodValueInVicinityNumber = 5;
	    colCon.preciseCC = correct4RubberBandArtifacts(colCon.preciseCC, goodones, MaximumGoodValueInVicinityNumber);		
	    
	    return colCon;
	}
	
	public double[] makeChecker(double[] values, int diffCheckWindowSize) {
		
		double[] diff = new double[values.length];
		diff[values.length - 1] = 0;
		
		for (int i = 0 ; i < values.length - 1; i++) {
			
			ArrayList<Double> diffcheck = new ArrayList<Double>();			
			
			//downstream
			for (int j = i - 1 ; j > i - 2 ; j--) if (j > 0) diffcheck.add((values[j] - values[i]) * (values[j] - values[i]));
			for (int j = i + 1 ; j < i + 2 ; j++) if (j < values.length) diffcheck.add((values[j] - values[i]) * (values[j] - values[i]));
			
			double[] diffCheckVector = new double[diffcheck.size()];
			for (int j = 0 ; j < diffcheck.size() ; j++) diffCheckVector[j] = diffcheck.get(j);
			
			diff[i] = StatUtils.max(diffCheckVector);
			
		}
		
		return diff;
		
	}
	
	public ColCoords3D correct4RubberBandArtifacts(ColCoords3D preciseCC, int[] goodone, int MaximumGoodValueInVicinityNumber) {
		
		int j,k,cc;
		ColCoords3D superPreciseCC = new ColCoords3D();
		
		double[] xmid = preciseCC.xmid;
		double[] ymid = preciseCC.ymid;				
		double[] zmid = preciseCC.zmid;
		double[] ixmid = preciseCC.ixmid; 
		double[] iymid = preciseCC.iymid;	
		double[] outerMinorRadius = preciseCC.outerMinorRadius;
		double[] outerMajorRadius = preciseCC.outerMajorRadius;
		double[] innerMinorRadius = preciseCC.innerMinorRadius; 
		double[] innerMajorRadius = preciseCC.innerMajorRadius;
		double[] wallThickness = preciseCC.wallThickness;
		double[] theta = preciseCC.theta;
		double[] itheta = preciseCC.itheta;
		double[] outerR2 = preciseCC.outerR2;
		double[] innerR2 = preciseCC.innerR2;
		boolean[] columnIsInThisDepth = preciseCC.columnIsAtThisDepth;
		int topOfColumn = preciseCC.topOfColumn;
		int heightOfColumn = preciseCC.heightOfColumn;
		int numberOfImputedLayers = preciseCC.numberOfImputedLayers;
		
		for (j = 0 ; j < goodone.length ; j++) {
			if (goodone[j] == 0 & j >= topOfColumn & j <= heightOfColumn + topOfColumn) {
				
				numberOfImputedLayers++;
				
				//search downstream
				k = j - 1; cc = 1;				
				ArrayList<Integer> downStream = new ArrayList<Integer>();
				while (k >= topOfColumn & cc <= MaximumGoodValueInVicinityNumber) {
					if (goodone[k] == 1) {
						downStream.add(k);
						cc++;
					}
					k--;
				}
				
				//search upstream
				k = j + 1; cc = 1;				
				ArrayList<Integer> upStream = new ArrayList<Integer>();
				while (k <= heightOfColumn + topOfColumn & cc <= MaximumGoodValueInVicinityNumber & k < goodone.length) {
					if (goodone[k] == 1) {
						upStream.add(k);
						cc++;
					}
					k++;
				}
				
				//calculate corrected value
				xmid[j] = plugItTogether(downStream, upStream, preciseCC.xmid);
				ymid[j] = plugItTogether(downStream, upStream, preciseCC.ymid);				
				ixmid[j] = plugItTogether(downStream, upStream, preciseCC.ixmid);
				iymid[j] = plugItTogether(downStream, upStream, preciseCC.iymid);
				outerMinorRadius[j] = plugItTogether(downStream, upStream, preciseCC.outerMinorRadius);
				outerMajorRadius[j] = plugItTogether(downStream, upStream, preciseCC.outerMajorRadius);
				innerMinorRadius[j] = plugItTogether(downStream, upStream, preciseCC.innerMinorRadius); 
				innerMajorRadius[j] = plugItTogether(downStream, upStream, preciseCC.innerMajorRadius);
				wallThickness[j] = plugItTogether(downStream, upStream, preciseCC.wallThickness);
				theta[j] = plugItTogether(downStream, upStream, preciseCC.theta);
				itheta[j] = plugItTogether(downStream, upStream, preciseCC.theta);
				outerR2[j] = plugItTogether(downStream, upStream, preciseCC.outerR2);
				innerR2[j] = plugItTogether(downStream, upStream, preciseCC.outerR2);
				
			}
		}
		
		//assign corrected values to output struct
		superPreciseCC.xmid = xmid; 
		superPreciseCC.ymid = ymid;				
		superPreciseCC.zmid = zmid;
		superPreciseCC.ixmid = ixmid; 
		superPreciseCC.iymid = iymid;	
		superPreciseCC.outerMinorRadius = outerMinorRadius;
		superPreciseCC.outerMajorRadius = outerMajorRadius;
		superPreciseCC.innerMinorRadius = innerMinorRadius; 
		superPreciseCC.innerMajorRadius = innerMajorRadius;
		superPreciseCC.wallThickness = wallThickness;
		superPreciseCC.theta = theta;
		superPreciseCC.itheta = theta;
		superPreciseCC.outerR2 = outerR2;
		superPreciseCC.innerR2 = outerR2;
		superPreciseCC.columnIsAtThisDepth = columnIsInThisDepth;
		superPreciseCC.topOfColumn = topOfColumn;
		superPreciseCC.heightOfColumn = heightOfColumn;
		superPreciseCC.numberOfImputedLayers = numberOfImputedLayers;
		
		return superPreciseCC;
		
	}
           
    public double plugItTogether(ArrayList<Integer> downStream, ArrayList<Integer> upStream, double[] ccFeature) {
    	
    	Median jMed = new Median();
    	double outValue;
    	double[] imputors =  new double[downStream.size() + upStream.size()];
    	
    	for (int k = 0 ; k < downStream.size() ; k++) imputors[k] = ccFeature[downStream.get(k)];
    	for (int k = downStream.size() ; k < imputors.length; k++) imputors[k] = ccFeature[upStream.get(k - downStream.size())];
    	
    	outValue = jMed.evaluate(imputors);
    	
    	return outValue;
    }

	public int[] findSmallCheckers(double[] checker, int filterSize) {
		
	      int[] goodone = new int[checker.length];
	   
	      for (int j = filterSize ; j < checker.length - filterSize ; j++) {
	    	  
	    	  double[] consider = new double[2 * filterSize + 1];
	    	  for (int i = 0 ; i < consider.length ; i++) consider[i] = checker[j - filterSize + i];
	    	  if (StatUtils.max(consider) > 0.1) goodone[j] = 0;
	    	  else goodone[j] = 1;
	      }

	      //also assign values to top and bottom bits of column
	      for (int j = 0 ; j < filterSize ; j++) goodone[j] = 0;
	      for (int j = checker.length - filterSize ; j < checker.length ; j++) goodone[j] = 0;
	      
	      return goodone;
	}
	
	
	public double[] findColumnTilt(double[] xCenter, double[] yCenter, double[] zCenter) {
		
		double[] xyzAngle = new double[3];
		double[] xzLine = new double[3];
		double[] yzLine = new double[3];
		double[] initialGuess = {2, 555};
			
		// first fit in the xz-plane
		CurveFitter xz = new CurveFitter(zCenter, xCenter);
		xz.setInitialParameters(initialGuess);
		xz.doFit(0);
		xzLine = xz.getParams();
				
		// first fit in the yz-plane				
		CurveFitter yz = new CurveFitter(zCenter, yCenter);
		yz.setInitialParameters(initialGuess);
		yz.doFit(0);				
		yzLine = yz.getParams();
			
		//calculate angle of tilt in xz-plane and yz-plane
		xyzAngle[0] = Math.atan(xzLine[1]);
		xyzAngle[1] = Math.atan(yzLine[1]);
		xyzAngle[2] = Math.atan(Math.sqrt(Math.tan(xyzAngle[0]) * Math.tan(xyzAngle[0]) + Math.tan(xyzAngle[1]) * Math.tan(xyzAngle[1])));
		
		return xyzAngle;
	}
	
	public ImageProcessor smoothFoundSurface(ImageProcessor prelimTopIP, PolygonRoi pRoi) {
		
		GaussianBlur gB = new GaussianBlur();
		
		//cut away things outsiude the ROI
		prelimTopIP.setRoi(pRoi);
		prelimTopIP.setBackgroundValue(0);
		prelimTopIP.setColor(0);
		prelimTopIP.fillOutside(pRoi);

		//init smoothSurface
		ImageProcessor smoothSurf = prelimTopIP.duplicate();

		//smooth surface majorly...
		//rf.rank(smoothSurf, 10, RankFilters.MEDIAN);
		double sigma = (prelimTopIP.getWidth() + prelimTopIP.getHeight()) / 100;
		gB.blurGaussian(smoothSurf, sigma, sigma, 0.01);
		
		//map smaller value of 
		ImageProcessor niceSurf = new ShortProcessor(prelimTopIP.getWidth(), prelimTopIP.getHeight());
		for (int x = 0 ; x < smoothSurf.getWidth() ; x++) {
			for (int y = 0 ; y < smoothSurf.getWidth() ; y++) {
				int prelimPix = prelimTopIP.getPixel(x, y);
				int smoothPix = smoothSurf.getPixel(x, y);
				if (prelimPix < smoothPix) niceSurf.putPixel(x, y, prelimPix);
				else niceSurf.putPixel(x, y, smoothPix);				
			}
		}
		
		return niceSurf;
		
	}
	
	public ImageProcessor findSurfaceSubroutine(ImagePlus nowTiff0, int minimumObjectThickness, String topOrBottom) {
				
		ImageManipulator jIM = new ImageManipulator();
		DisplayThings disp = new DisplayThings();
		RollerCaster rC = new RollerCaster();
		MorphologyAnalyzer mA = new MorphologyAnalyzer();
		FloodFillComponentsLabeling3D myFCCL = new FloodFillComponentsLabeling3D(26, 32);
		
		Dilate_ mD = new Dilate_();
		Erode_ mE = new Erode_();
		
		int i;
		int startSlice = 0;
		int stopSlice = 300;
		int slicesPerChunk = 2;		
		int volumeOfChunks2Exclude = 100000;
		//PolygonRoi[] pRoi = new PolygonRoi[pRoi0.length];
		double xmid = nowTiff0.getWidth() / 2;
		double ymid = nowTiff0.getHeight() / 2;		
		
		ImageStack inverseStack = new ImageStack(nowTiff0.getWidth(), nowTiff0.getHeight());
		ImagePlus inverseTop = new ImagePlus();
		ImagePlus nowTiff = new ImagePlus();
	
		if (topOrBottom.contentEquals("bottom")) {
			nowTiff = jIM.flipTiff(nowTiff0);
			//for (i = pRoi0.length - 1; i >= 0 ; i--) pRoi[pRoi0.length - 1 - i] = pRoi0[i];			
		} else {
			//pRoi = pRoi0;
			nowTiff = nowTiff0;
		}
		
		//make ROI mask of first slice: all pixels that are already 0 will not be considered 
		Selection mS = new Selection(); 
		ThresholdToSelection mT2S = new ThresholdToSelection();
		nowTiff.setSlice(1);
		ImageProcessor firstSlice = nowTiff.getProcessor().duplicate();
		firstSlice.dilate();
		firstSlice.setThreshold(0, 0, ImageProcessor.NO_LUT_UPDATE);				
		Roi fRoi = mT2S.convert(firstSlice);		
		firstSlice.invert();
		Roi iRoi = mT2S.convert(firstSlice);
					
		//find approximate location of upper soil surface
		int approximateTopSurface = 9999;		
		for (i = 0 ; i < nowTiff.getNSlices() ; i++) {
			
			IJ.showStatus("Finding approximate location of soil top surface #" + (i + 1) + "/" + nowTiff.getNSlices());
			
			nowTiff.setPosition(i + 1);			
			ImageProcessor myIP = nowTiff.getProcessor();		
			
			myIP.setRoi(fRoi);		//only consider canvas that had been airfilled in the fisrt slice
			int[] myHist = myIP.getHistogram();
			float solids = myHist[0];
			float all = myHist[0] + myHist[255];
			float myFractionOfSoil = solids / all;
			if (approximateTopSurface == 9999 & myFractionOfSoil > 0.5) {
				approximateTopSurface = i;
				break;
			}
						
		}
		if (approximateTopSurface > 400) startSlice = approximateTopSurface - 400;
		stopSlice = approximateTopSurface + 300;
		
		//create inverse of binary 
		for (i = startSlice ; i < stopSlice ; i++) {
			
			IJ.showStatus("Inverting slice #" + (i + 1 - startSlice) + "/" + (stopSlice - startSlice));
			
			nowTiff.setPosition(i + 1);
			ImageProcessor invIP = nowTiff.getProcessor().duplicate(); 
			
			invIP.invert();
			invIP.setRoi(fRoi);
			invIP.setColor(0);			
			invIP.fill(fRoi);
			
			inverseStack.addSlice(invIP);			
		}
		inverseTop.setStack(inverseStack);
		
		//inverseTop.updateAndDraw();
		//inverseTop.show();
				
		//find the real soil matrix..  (get rid of smaller aggregates on top);
		IJ.showStatus("Removing loose soil aggregates lying on top of the surface using the ParticleAnalyzer  ...");		
		ImagePlus purifiedTop = new ImagePlus();
		ImageStack myLabelStack = myFCCL.computeLabels(inverseTop.getStack());			
		purifiedTop.setStack(BinaryImages.binarize(LabelImages.volumeOpening(myLabelStack, volumeOfChunks2Exclude)));
				
		IJ.freeMemory();IJ.freeMemory();
		
		//create xy grid
		ImageProcessor outIP = new ShortProcessor(nowTiff.getWidth(), nowTiff.getHeight());
		int[][] xy = new int[outIP.getWidth()][outIP.getHeight()];
		for (int x = 0 ; x < outIP.getWidth() ; x++) {
			for (int y = 0 ; y < outIP.getHeight() ; y++) {
				xy[x][y] = 0;
			}
		}
	
		//find the soil surface		
		int latestValue = 0;
		for (int x = 0 ; x < purifiedTop.getWidth() ; x++) {
			for (int y = 0 ; y < purifiedTop.getHeight() ; y++) {
				
				IJ.showStatus("Searching in pixel (x | y) = " + x + ", " + y + " ... found surface at depth " + latestValue + " vx");
				
				pixelLoop:
				for (int z = 0 ; z < purifiedTop.getNSlices() ; z++) {
				
					purifiedTop.setPosition(z+1);
					ImageProcessor myIP = purifiedTop.getProcessor();
					
					int thisImgPixel = myIP.getPixel(x, y);
					
					if (thisImgPixel > 0) {						
						int theValue2Put = z + startSlice;
						xy[x][y] = theValue2Put;	
						latestValue = theValue2Put;
						break pixelLoop;		
					}
				}		
			}
		}

		//transfer results to image		
		for (int x = 0 ; x < outIP.getWidth() ; x++) {
			for (int y = 0 ; y < outIP.getHeight() ; y++) {
				if (xy[x][y] == 9999) xy[x][y] = 0;
				outIP.putPixelValue(x, y, xy[x][y]);
			}			
		}
		
		//calculate ruggedness of surface
		HistogramStuff mH = new HistogramStuff();
		outIP.setRoi(iRoi);
		int[] myHist = outIP.getHistogram();
		float mean = mH.findMeanFromHistogram(myHist);
		int median = mH.findMedianFromHistogram(myHist);
		int min = mH.findMinFromHistogram(myHist);
		int max = mH.findMaxFromHistogram(myHist);
		double stdev = mH.findStdFromHistogram(myHist, mean);
		
		//apply median filter to outIP to create lower cut-off		
		//RankFilters mRF = new RankFilters();
		//mRF.rank(outIP, 25.0, RankFilters.MEDIAN);		
		
		//ImagePlus result = new ImagePlus("", outIP);
		//result.show();
		
		return outIP;				
	}
	
	public ImagePlus findSoilSurface(ImagePlus nowTiff, MenuWaiter.SurfaceFinderReturn mSFR) {
				
		//init objects
		InputOutput jIO = new InputOutput();
		RoiHandler roi = new RoiHandler();
				
		//init varis
		ImagePlus mySurface = new ImagePlus();
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		ImageProcessor topIP = new ShortProcessor(nowTiff.getWidth(), nowTiff.getHeight());
		ImageProcessor bottomIP = new ShortProcessor(nowTiff.getWidth(), nowTiff.getHeight());
	
		int minimumObjectThickness = mSFR.neglectCrumbsOfLessThan;
		
		//find top and bottom surface considering the air-phase / non-air-phase boundary
		ImageProcessor prelimTopIP = findSurfaceSubroutine(nowTiff, minimumObjectThickness, "top");  
		ImageProcessor prelimBottomIP = findSurfaceSubroutine(nowTiff, minimumObjectThickness, "bottom");
		
		//transfer surface IPs to imagePlus
		outStack.addSlice(prelimTopIP);
		outStack.addSlice(prelimBottomIP);		
		mySurface.setStack(outStack);
		
		//mySurface.draw();mySurface.show();
		
		return mySurface;		
	}
	
	/*private ImageProcessor smoothAwaySmallPores(ImageProcessor origIP, double smootherRadius, PolygonRoi pRoi) {
		
		if (smootherRadius > 0) {
		
			for (int i = 0 ; i < (int)Math.round(origIP.getMax()) ; i++) {
				
				IJ.showStatus("Filling pores opening to the surface in layer " + (i + 1));
				
				ImageProcessor smoothIP = origIP.duplicate();	
				smoothIP.min(i);
				smoothIP.max(i + 1);
				smoothIP.subtract(i);
				smoothIP.multiply(255);
				
				ImageProcessor binIP = smoothIP.convertToByte(false);	
				ImagePlus binTiff = new ImagePlus("", binIP);
				
				//binTiff.draw();binTiff.show();
	
				//apply the particle counter to measure the different pore sizes;
				int imgVol = binTiff.getWidth() * binTiff.getHeight() * binTiff.getNSlices();
				Counter3D myOC = new Counter3D(binTiff, 128, 0, imgVol, false, false);
				myOC.getObjects();		
				ImagePlus myHoles = myOC.getObjMap();
				ImageProcessor holeIP = myHoles.getProcessor();		
			
				Vector allObjects = myOC.getObjectsList();
				int numOfObjects = allObjects.size();
				
				double[] volumes = new double[numOfObjects];
				int[][] bBoxXYZ = new int[numOfObjects][3];
				int[][] bBoxEXEYEZ = new int[numOfObjects][3];
				double[][] centroids = new double[numOfObjects][3];
				ArrayList<Integer> nowOut = new ArrayList<Integer>();
				for (int iii = 0; iii < numOfObjects; iii++){		         
					Object3D currObj=(Object3D) allObjects.get(iii);
					
					//volumes
					volumes[iii] = currObj.size;				
					
					//spot the holes that should be filled
					if (volumes[iii] < smootherRadius * smootherRadius * Math.PI) {
					
						int[] tmpArrayInt=currObj.bound_cube_TL;
						
						//fill holes that were spotted
						int xmin = tmpArrayInt[0] - 1;
						int xmax = xmin + currObj.bound_cube_width + 1;
						int ymin = tmpArrayInt[1] - 1;
						int ymax = ymin + currObj.bound_cube_height + 1;
						
						for (int x = xmin ; x < xmax + 1 ; x++) {
							for (int y = ymin ; y < ymax + 1 ; y++) {
								int holeCheck = (int)Math.round(holeIP.getPixelValue(x, y));
								if (holeCheck == iii) {
									origIP.putPixel(x, y, iii);	//if is hole than fill it..
								}
							}
						}						
					}

				}				
				
				IJ.freeMemory();IJ.freeMemory();				
			}
		
		}
		
		return origIP;
		
	}*/
	
	public double[] calcHeightAndDiameterOfColumn(ImagePlus soilSurface, int numberOfSlices) {
		
		double[] hod = new double[2];
		int iW = soilSurface.getWidth();
		int iH = soilSurface.getHeight();
		
		//calculate top of soil
		soilSurface.setPosition(1);
		ImageProcessor surIP = soilSurface.getProcessor();
		double depthSum = 0;
		double area = 0;
		for (int x = 0 ; x < iW ; x++) {
			for (int y = 0 ; y < iH ; y++) {
				int nowPixel = (int)surIP.getPixelValue(x, y);
				if (nowPixel > 0) {
					area++;
					depthSum += nowPixel;
				}
			}
		}
		double topMinus = depthSum / area;
		hod[1] = 2 * Math.sqrt(area / Math.PI);
		
		//calculate top of soil
		soilSurface.setPosition(2);
		surIP = soilSurface.getProcessor();
		depthSum = 0;
		area = 0;
		for (int x = 0 ; x < iW ; x++) {
			for (int y = 0 ; y < iH ; y++) {
				int nowPixel = (int)surIP.getPixelValue(x, y);
				if (nowPixel > 0) {
					area++;
					depthSum += nowPixel;
				}
			}
		}
		
		double botMinus = depthSum / area;
		
		hod[0] = numberOfSlices - topMinus - botMinus;		
		
		return hod;
	}
	
	public double[] findMedianWallgrayValues(ImagePlus nowTiff, ObjectDetector.ColCoords3D preciseCC) {

		RoiHandler roi = new RoiHandler();
		
		int i;		
		
		ImagePlus origTiff = nowTiff.duplicate();
		
		PolygonRoi[] iRoi = roi.makeMeAPolygonRoiStack4RawImage("inner", "wide", preciseCC, 0);
		PolygonRoi[] oRoi = roi.makeMeAPolygonRoiStack4RawImage("outer", "tight", preciseCC, 0);
		
		double[] wallMedian = new double[preciseCC.heightOfColumn];
		
		//sample column wall gray
		for (i = preciseCC.topOfColumn ; i < preciseCC.topOfColumn + preciseCC.heightOfColumn ; i++) { 
			
			IJ.showStatus("Sampling illumination of slice " + (i - preciseCC.topOfColumn + 1) + "/" + (preciseCC.heightOfColumn));	
			
			origTiff.setPosition(i+1);
			ImageProcessor myIP = origTiff.getProcessor();		
			
			//origTiff.draw();
			//origTiff.show();
					
			myIP.setRoi(oRoi[i - preciseCC.topOfColumn]);
			myIP.setValue(0);
			//myIP.setColor(Color.white);
			//myIP.draw(oRoi[i]);
			myIP.fillOutside(oRoi[i - preciseCC.topOfColumn]);		
			
			myIP.setRoi(iRoi[i - preciseCC.topOfColumn]);
			myIP.setValue(0);
			//myIP.setColor(Color.white);
			//myIP.draw(iRoi[i]);
			myIP.fill(iRoi[i - preciseCC.topOfColumn]);	
			
			myIP.resetRoi();
			
			//origTiff.draw();
			//origTiff.show();
				
			int[] myHist = myIP.getHistogram();
			float[] cHist = new float[myHist.length];	
			cHist[0] = 0;
			for (int j = 1 ; j < myHist.length ; j++) cHist[j] = cHist[j-1] + (float)myHist[j];
			for (int j = 1 ; j < myHist.length ; j++) cHist[j] = cHist[j] / cHist[cHist.length - 1];			
			for (int j = 1 ; j < myHist.length ; j++) if (cHist[j] > 0.5) {
				wallMedian[i - preciseCC.topOfColumn] = j;
				break;				
			}
	
		}
		
/*		double[] xax =  new double[preciseCC.heightOfColumn];
		for (i = 0 ; i < preciseCC.heightOfColumn ; i++) xax[i] = i; 		
		Plot nP = new Plot("","","",xax, wallMedian);
		nP.setLimits(0, xax.length, 6500, 7500);
		nP.draw();
		nP.show();*/
		
		return wallMedian;
		
	}	
	
	public int[] postProcessSteelgrayValues(double[] medianOfWallgrayValue) {
		
		int[] topBotHeight = new int[3];
		Median jMed = new Median();
		
		int steelMedian = (int)Math.round(jMed.evaluate(medianOfWallgrayValue));	
		
		int firstGood = -1;
		int lastGood = medianOfWallgrayValue.length - 1;
		for (int i = 0 ; i < medianOfWallgrayValue.length ; i++) {
			if (firstGood == -1) if (medianOfWallgrayValue[i] > 0.9 * steelMedian) firstGood = i;
			if (firstGood > 0 & lastGood == medianOfWallgrayValue.length - 1) if (medianOfWallgrayValue[i] < 0.9 * steelMedian) lastGood = i - 1;
		}
				
		topBotHeight[0] = firstGood;
		topBotHeight[1] = lastGood; 
		topBotHeight[2] = lastGood - firstGood + 1;
		
		return topBotHeight;
		
	}
	
	

	public MorphologyAnalyzer.SurfaceStatistics extractSurfaceStatistics(String myOutPath, String myTiff, ImagePlus nowTiff, String nowGaugePath, int numberOfSlices) {

		MorphologyAnalyzer morph = new MorphologyAnalyzer();
		MorphologyAnalyzer.SurfaceStatistics mSS = morph.new SurfaceStatistics();
		RoiHandler jRH = new RoiHandler();	
		InputOutput jIO = new InputOutput();
		HistogramStuff hist = new HistogramStuff();
		
		//read InnerCircle file		
		ColCoords3D jCO = new ColCoords3D();
		int versio = jIO.checkInnerCircleFileVersion(nowGaugePath);			
		if (versio == 0) jCO = jIO.readInnerCircleVer0(nowGaugePath);	
		else jCO = jIO.readInnerCircleVer1(nowGaugePath);	
		
		PolygonRoi[] pRoi = jRH.makeMeAPolygonRoiStack("inner", "wide", jCO, 0);
		
		//get top stats
		nowTiff.setPosition(1);
		ImageProcessor topIP = nowTiff.getProcessor();
		topIP.setRoi(pRoi[0]);		
		int[] topHist = topIP.getHistogram();
		topHist[0] = 0;
		
		mSS.highestElevation = hist.findMinFromHistogram(topHist);
		mSS.medianElevation = hist.findMedianFromHistogram(topHist);
		mSS.meanElevation = (int)Math.round(hist.findMeanFromHistogram(topHist));
		mSS.lowestElevation = hist.findMaxFromHistogram(topHist);
		
		//get bottom stats
		nowTiff.setPosition(2);
		ImageProcessor bottomIP = nowTiff.getProcessor();
		bottomIP.setRoi(pRoi[pRoi.length - 1]);
	
		int[] bottomHist = bottomIP.getHistogram();
		mSS.highestIntrusion = numberOfSlices - hist.findMaxFromHistogram(bottomHist);
		mSS.medianIntrusion = numberOfSlices - hist.findMedianFromHistogram(bottomHist);
		mSS.meanIntrusion = numberOfSlices - (int)Math.round(hist.findMeanFromHistogram(bottomHist));
		mSS.lowestIntrusion = numberOfSlices - hist.findMinFromHistogram(bottomHist);
		
		return mSS;
				
	}
	
	public MorphologyAnalyzer.SurfaceStatistics extractSurfaceStatistics(InputOutput.MyFileCollection mFC) {

		String myOutPath = mFC.myOutFolder;
		String myTiff = mFC.fileName; 
		int numberOfSlices = mFC.nOfSlices;
		
		MorphologyAnalyzer morph = new MorphologyAnalyzer();
		MorphologyAnalyzer.SurfaceStatistics mSS = morph.new SurfaceStatistics();
		RoiHandler roi = new RoiHandler();	
		InputOutput jIO = new InputOutput();
		HistogramStuff hist = new HistogramStuff();
		
		//find correct InnerCircle and Surface Files
		String[] GandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC);
		
		//read InnerCircle file and create outline ROI	
		PolygonRoi[] pRoi = new PolygonRoi[mFC.nowTiff.getNSlices()];
		if (GandS[0].contains("Gauge")) {
			
			ColCoords3D jCO = new ColCoords3D();				
			int versio = jIO.checkInnerCircleFileVersion(mFC.nowInnerCirclePath);			
			if (versio == 0) jCO = jIO.readInnerCircleVer0(mFC.nowInnerCirclePath);	
			else jCO = jIO.readInnerCircleVer1(mFC.nowInnerCirclePath);	
			
			int wallThickness = (int) (0.9 * Math.round(StatUtils.max(jCO.wallThickness)));				
		
			pRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 4);
			
		}
		
		else {

			//read column coordinates
			ObjectDetector.EggShapedColCoords3D jCO = new EggShapedColCoords3D();						
			jCO = jIO.readInnerCircleSteel(mFC);
			
			pRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 4);		
	
		}
		
		//load surface file
		String mySurfaceFilePath = mFC.mySurfaceFolder + mFC.pathSep + GandS[1];
		ImagePlus nowTiff =  jIO.openTiff3D(mySurfaceFilePath);
		
		//get top stats
		nowTiff.setPosition(1);
		ImageProcessor topIP = nowTiff.getProcessor();
		topIP.setRoi(pRoi[0]);		
		int[] topHist = topIP.getHistogram();
		topHist[0] = 0;
		
		mSS.highestElevation = hist.findMinFromHistogram(topHist);
		mSS.medianElevation = hist.findMedianFromHistogram(topHist);
		mSS.meanElevation = (int)Math.round(hist.findMeanFromHistogram(topHist));
		mSS.lowestElevation = hist.findMaxFromHistogram(topHist);
		
		//get bottom stats
		nowTiff.setPosition(2);
		ImageProcessor bottomIP = nowTiff.getProcessor();
		bottomIP.setRoi(pRoi[pRoi.length - 1]);
	
		int[] bottomHist = bottomIP.getHistogram();
		mSS.highestIntrusion = numberOfSlices - hist.findMaxFromHistogram(bottomHist);
		mSS.medianIntrusion = numberOfSlices - hist.findMedianFromHistogram(bottomHist);
		mSS.meanIntrusion = numberOfSlices - (int)Math.round(hist.findMeanFromHistogram(bottomHist));
		mSS.lowestIntrusion = numberOfSlices - hist.findMinFromHistogram(bottomHist);
		
		return mSS;
				
	}
	
	public double calculateCVOfWallBrightness(ImageProcessor myIP, FitStuff.FittedEllipse fI, FitStuff.FittedEllipse fO) {
		
		RoiHandler roi = new RoiHandler();
		HistogramStuff hist = new HistogramStuff();
		
		PolygonRoi inner = roi.makeRoiFromFittedEllipse(fI);
		PolygonRoi outer = roi.makeRoiFromFittedEllipse(fO);
		
		ImageProcessor nowIP = myIP.duplicate();
		nowIP.setValue(0);
		nowIP.setBackgroundValue(0);
		
		//get wall illumination
		nowIP.fillOutside(outer);			
		nowIP.fill(inner);
		
		int[] myHist = nowIP.getHistogram();
		myHist[0] = 0; //set zeros to zero
		
		float myMean = hist.findMeanFromHistogram(myHist);
		double myStd = hist.findStdFromHistogram(myHist, myMean);
		
		double myCVofWallgrayValues = myStd / myMean;
		
		return myCVofWallgrayValues;
		
	}
	
	public double calculateCVOfWallThickness(FitStuff.FittedEllipse fI, FitStuff.FittedEllipse fO) {
		
		RoiHandler roi = new RoiHandler();
		
		PolygonRoi iRoi = roi.makeRoiFromFittedEllipse(fI);
		PolygonRoi oRoi = roi.makeRoiFromFittedEllipse(fO);
		
		int[] iX = iRoi.getXCoordinates();
		int[] iY = iRoi.getYCoordinates();
		int[] oX = oRoi.getXCoordinates();
		int[] oY = oRoi.getYCoordinates();
		
		double ixb = iRoi.getXBase();
		double iyb = iRoi.getYBase();
		double oxb = oRoi.getXBase();
		double oyb = oRoi.getYBase();
		
		double[] wallThickness = new double[iX.length];
		
		for (int i = 0 ; i < iX.length ; i++) {
			
			double dX = ((double)oX[i] + oxb) - ((double)iX[i] + ixb);
			double dY = ((double)oY[i] + oyb) - ((double)iY[i] + iyb);
			
			wallThickness[i] = Math.sqrt(dX * dX + dY * dY);
			
		}
		
		double stdev = Math.sqrt(StatUtils.variance(wallThickness));
		double avg = StatUtils.mean(wallThickness);
		
		double CV = stdev / avg;
		
		return CV;
		
	}

	public double calculateRatioBetweenOuterAndInnerRadius(FitStuff.FittedEllipse fI, FitStuff.FittedEllipse fO) {
		
		RoiHandler roi = new RoiHandler();
		
		PolygonRoi iRoi = roi.makeRoiFromFittedEllipse(fI);
		PolygonRoi oRoi = roi.makeRoiFromFittedEllipse(fO);
		
		int[] iX = iRoi.getXCoordinates();
		int[] iY = iRoi.getYCoordinates();
		int[] oX = oRoi.getXCoordinates();
		int[] oY = oRoi.getYCoordinates();
		
		double ixb = iRoi.getXBase();
		double iyb = iRoi.getYBase();
		double oxb = oRoi.getXBase();
		double oyb = oRoi.getYBase();

		double[] ratios = new double[iX.length];
		
		for (int i = 0 ; i < iX.length ; i++) {
			
			double iLx = (double)iX[i] + ixb;
			double iLy = (double)iY[i] + iyb;
			double oLx = (double)oX[i] + oxb;
			double oLy = (double)oY[i] + oyb;
						
			double rInner = Math.sqrt(iLx * iLx + iLy * iLy);
			double rOuter = Math.sqrt(oLx * oLx + oLy * oLy);
			
			ratios[i] = rOuter / rInner;
			
		}
				
		return StatUtils.min(ratios);
		
	}
	
	public ImagePlus extractAirFilledPores(int tensionStep, ImagePlus nowTiff, ImagePlus surfTiff, MenuWaiter.WRCCalculatorMenu mWRC, InputOutput.MyFileCollection mFC) {
		
		ImageManipulator jIM = new ImageManipulator();
		MorphologyAnalyzer mA = new MorphologyAnalyzer();
		
		double columnHeight = nowTiff.getNSlices() * mWRC.voxelSizeInMicroMeter;
		double bboundary = mWRC.tensionStepsInMM[tensionStep] * 1000;  //calculate in micrometer
		double capillaryConstantInMicroMeter = 1.48e7 * Math.cos(mWRC.wettingAngle / 180 * Math.PI);
		
		ImagePlus airTiff = jIM.binarize3DImage(nowTiff, 0.5, (double)tensionStep + 1.5);
		
		IJ.showStatus("Identifying connected pore-clusters ...");
		double constant = (2 * capillaryConstantInMicroMeter) / (mWRC.connectingSubresolutionDiameter * mWRC.voxelSizeInMicroMeter);
		double tensionTerm = (constant - bboundary) / mWRC.voxelSizeInMicroMeter;
		double columnHeightInVX = columnHeight / mWRC.voxelSizeInMicroMeter;
		double drainingDepthInVX = columnHeightInVX - tensionTerm + 0.5;
		if (drainingDepthInVX < 0) drainingDepthInVX = 0;
		ImagePlus outTiff = mA.findClusterConnected2Top(airTiff, drainingDepthInVX);
		
		return outTiff;
		
	}
	
	public ImagePlus extractAirFilledPores(double tension, ImagePlus nowTiff, ImagePlus surfTiff, MenuWaiter.DrainageSimulatorOptions mDS, InputOutput.MyFileCollection mFC) {
		
		MorphologyAnalyzer mA = new MorphologyAnalyzer();
		
		double columnHeight = nowTiff.getNSlices() * mDS.voxelSizeInMicroMeter;
		double bboundary = tension * 1000;  //calculate in micrometer
		double capillaryConstantInMicroMeter = 1.48e7 * Math.cos(mDS.wettingAngle / 180 * Math.PI);
		
		ImagePlus outTiff = new ImagePlus();
		ImageStack zwischiStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImagePlus zwischiTiff = new ImagePlus();
		
		//get all pores with a diameter larger than the one corresponding to the respective matrix potential	
		for (int z = 1 ; z <= nowTiff.getNSlices() ; z++) {
			
			nowTiff.setPosition(z);
			ImageProcessor nowIP = nowTiff.getProcessor();
			double columnNow = (nowTiff.getNSlices() - (z - 0.5)) * mDS.voxelSizeInMicroMeter;  //-0.5 to get to the voxel midpoint
			double cutoffRadius = capillaryConstantInMicroMeter / (bboundary + columnNow);
			double cutoffThickness = (2 * cutoffRadius) / mDS.voxelSizeInMicroMeter;			
						
			IJ.showStatus("Extracting air-filled pores at depth " + String.format("%2.2f", (columnHeight - columnNow) / 10000) + " cm");
			
			//init out image
			ImageProcessor outIP = new ByteProcessor(nowTiff.getWidth(), nowTiff.getHeight());
			
			for (int x = 0 ; x < outIP.getWidth() ; x++) {
				for (int y = 0 ; y < outIP.getHeight(); y++) {
					double nowPix = nowIP.getPixelValue(x, y);
					if (nowPix <= cutoffThickness) outIP.putPixel(x, y, 0);
					else outIP.putPixel(x, y, 255);
				}
			}
			
			zwischiStack.addSlice(outIP);			
		}
		
		zwischiTiff.setStack(zwischiStack);		

		
		//zwischiTiff.updateAndDraw();zwischiTiff.show();
					
		IJ.showStatus("Identifying connected pore-clusters ...");
		double constant = (2 * capillaryConstantInMicroMeter) / mDS.voxelSizeInMicroMeter;
		double tensionTerm = (constant - bboundary) / mDS.voxelSizeInMicroMeter;
		double columnHeightInVX = columnHeight / mDS.voxelSizeInMicroMeter;
		double drainingDepthInVX = columnHeightInVX - tensionTerm + 0.5;
		if (drainingDepthInVX < 0) drainingDepthInVX = 0;
		outTiff = mA.findClusterConnected2Top(zwischiTiff, drainingDepthInVX);
		
		return outTiff;
		
	}
	
	public ImagePlus createDrainageMap(double[] tension, ImagePlus nowTiff, ImagePlus surfTiff, MenuWaiter.WRCCalculatorMenu mDS, InputOutput.MyFileCollection mFC) {
			
		double columnHeight = nowTiff.getNSlices() * mDS.voxelSizeInMicroMeter;
		double[] bboundary = new double[tension.length];
		for (int i = 0 ; i < tension.length ; i++) bboundary[i] = tension[i] * 1000;  //calculate in micrometer 
		double capillaryConstantInMicroMeter = 1.48e7 * Math.cos(mDS.wettingAngle / 180 * Math.PI);
		
		ImageStack zwischiStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImagePlus zwischiTiff = new ImagePlus();
		
		//get all pores with a diameter larger than the one corresponding to the respective matrix potential	
		for (int z = 1 ; z <= nowTiff.getNSlices() ; z++) {
			
			nowTiff.setPosition(z);
			ImageProcessor nowIP = nowTiff.getProcessor();
			double columnNow = (nowTiff.getNSlices() - (z - 0.5)) * mDS.voxelSizeInMicroMeter;  //-0.5 to get to the voxel midpoint
			double[] cutOffRadius = new double[tension.length]; 
			for (int i = 0 ; i < tension.length ; i++) cutOffRadius[i] = capillaryConstantInMicroMeter / (bboundary[i] + columnNow);
			double[] cutOffThickness = new double[tension.length]; 
			for (int i = 0 ; i < tension.length ; i++) cutOffThickness[i] = (2 * cutOffRadius[i]) / mDS.voxelSizeInMicroMeter;			
						
			IJ.showStatus("Creating drainage map at depth " + String.format("%2.2f", (columnHeight - columnNow) / 10000) + " cm");
			
			//init out image
			ImageProcessor outIP = new ByteProcessor(nowTiff.getWidth(), nowTiff.getHeight());
			
			for (int x = 0 ; x < outIP.getWidth() ; x++) {
				for (int y = 0 ; y < outIP.getHeight(); y++) {
					double nowPix = nowIP.getPixelValue(x, y);
					for (int i = 0 ; i < tension.length ; i++) {
						if (nowPix > 0) {
							if (nowPix <= cutOffThickness[i]) outIP.putPixel(x, y, tension.length + 1);
							else {
								outIP.putPixel(x, y, i + 1);
								break;
							}
						}
					}
				}
			}
			
			zwischiStack.addSlice(outIP);			
		}
		
		zwischiTiff.setStack(zwischiStack);		
		
		//zwischiTiff.updateAndDraw();zwischiTiff.show();
					
		return zwischiTiff;
	}	

	public class easyAccessTable_2DMV {
	// 2d matrices with easy-access functions to get 1D array (column), access multiple columns combined with row, etc. (a little like in R but not thought through as well =)
		
		public double [][] tab2D;
		
		public double[][] getColumnA (double [][] tab2D, String columnName, int rownumber){
			return tab2D;
		}
				
	}
}
	
//public RadialModes getRadialPVCIlluminationDEPRECATED(ImagePlus nowTiff, ColCoords3D jCO, MenuWaiter.BeamDeHardeningReturn mBDH) {
//
////init units
//RoiHandler roi = new RoiHandler();
//HistogramStuff hist = new HistogramStuff();
//
////init variables
//int cc=0;
//RadialModes myRM = new RadialModes(); 		
//int[] myThresh = new int[jCO.heightOfColumn];
//RankFilters myRF = new RankFilters();
//
////prepare Rois..
//PolygonRoi[] nowRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 2);
//double averageRadius = StatUtils.percentile(jCO.innerMinorRadius,50);
//for (int frisbee = 0 ; frisbee < (1 - mBDH.cutoff) * averageRadius  ; frisbee += mBDH.stepsize * averageRadius) cc++;
//double[][] maskedModeBrightness = new double[nowTiff.getNSlices()][cc + 2];
//double[][] maskedMinimumBrightness = new double[nowTiff.getNSlices()][cc + 2];
//
//PolygonRoi[][] allRoi = new PolygonRoi[jCO.heightOfColumn][cc + 2];
//for (int frisbee = 0 ; frisbee < cc ; frisbee++) {
//int nowFris = -(int)Math.round(frisbee * mBDH.stepsize * averageRadius); 
//PolygonRoi[] frisbeeRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, nowFris);
//for (int i = 0 ; i < frisbeeRoi.length ; i++) allRoi[i][frisbee] = frisbeeRoi[i];			
//}
//double[] r = new double[cc + 2];
//r[0] = 0;											// add a pseudo point at 0 
//r[cc + 1] = averageRadius;  // and as well one at the column wall...
//for (int i = 0 ; i < cc ; i++) r[cc - i] = averageRadius - (mBDH.stepsize * averageRadius / 2 + i * mBDH.stepsize * averageRadius);		
//
////sample illumination
////for (int i = 599; i < 800 ; i++) {
//for (int i = 0; i < jCO.heightOfColumn ; i++) {
//
//IJ.showStatus("Sampling radial illumination of slice #" + (i + 1) + "/" + nowTiff.getNSlices());
//
////set image to next slice
//nowTiff.setPosition(i+1);
//ImageProcessor nowIP = nowTiff.getProcessor();
//
////apply threshold to mask out the low density phase
//ImageProcessor maskedIP = nowIP.duplicate();
//maskedIP.setRoi(nowRoi[i]);	
//if (mBDH.maskingMethod != "None") {
//maskedIP.setAutoThreshold(mBDH.maskingMethod);
//myThresh[i] = maskedIP.getAutoThreshold();
//}
//else myThresh[i] = 0;
//
////get illumination
//for (int frisbee = 0 ; frisbee < cc - 1; frisbee++) {			
//ImageProcessor ceoIP = maskedIP.duplicate();
//
////apply a median filter
//if (mBDH.radiusOfMedianFilter > 1) myRF.rank(ceoIP, mBDH.radiusOfMedianFilter, RankFilters.MEDIAN);
//
////create donut shaped image segment
//ceoIP.setValue(0);	
//ceoIP.fillOutside(allRoi[i][frisbee]);
//ceoIP.fill(allRoi[i][frisbee+1]);
//ceoIP.resetRoi();
//
////get histogram 
//int[] nowHist = ceoIP.getHistogram();
//
////sample minimum gray values as a proxy for the air phase brightness
//maskedMinimumBrightness[i][cc - frisbee] = hist.findNonZeroMinimumFromHistogram(nowHist, (int)mBDH.airPhaseClassMemberMinimum);
//
////set all values below the threshold to 0.. so that only the soil matrix is left in the images 
//for (int j = 0 ; j < myThresh[i] ; j++) nowHist[j] = 0;
//
////pick out mode gray value
//maskedModeBrightness[i][cc - frisbee] = hist.findNonZeroModeFromHistogram(nowHist);
//
///*			if (i == 100) {
//ImagePlus testImg = new ImagePlus("Slice",ceoIP);
//testImg.updateAndDraw();
//testImg.show();
//}*/
//}
////and also do the same for the central area ..
//ImageProcessor ceoIP = maskedIP.duplicate();
//ceoIP.setValue(0);	
//ceoIP.fillOutside(allRoi[i][cc - 1]);
//ceoIP.resetRoi();
//int[] nowHist = ceoIP.getHistogram();	
//maskedMinimumBrightness[i][1] = hist.findNonZeroMinimumFromHistogram(nowHist, (int)mBDH.airPhaseClassMemberMinimum);
//for (int j = 0 ; j < myThresh[i] ; j++) nowHist[j] = 0;
//maskedModeBrightness[i][1] = hist.findNonZeroModeFromHistogram(nowHist);	
//
////assign the neighboring gray values to the two pseudo points
//maskedMinimumBrightness[i][0] = maskedMinimumBrightness[i][1];
//maskedMinimumBrightness[i][cc + 1] = maskedMinimumBrightness[i][cc];
//maskedModeBrightness[i][0] = maskedModeBrightness[i][1];
//maskedModeBrightness[i][cc + 1] = maskedModeBrightness[i][cc];
//
////fix skip
//double skipTarget = maskedModeBrightness[i][cc - mBDH.skip];
//for (int j = cc - mBDH.skip + 1 ; j < cc + 2 ; j++) maskedModeBrightness[i][j] = skipTarget;
//}
//
//myRM.maskedRadialModes = maskedModeBrightness;
//myRM.maskedRadialMinima = maskedMinimumBrightness;
//myRM.radius = r;
//myRM.maskingThreshold = myThresh;
//
//return myRM;
//
//}

/*public double[] findMedianSteelgrayValues(ImagePlus nowTiff, ObjectDetector.ColumnCoordinates preciseCC) {

int i;		

ImagePlus origTiff = nowTiff.duplicate();

double[] wallQ = new double[nowTiff.getNSlices()];
double[] q1 = new double[nowTiff.getNSlices()];		

//sample column wall gray
for (i = 0 ; i < nowTiff.getNSlices() ; i++) { 
	
	IJ.showStatus("Sampling illumination of slice " + (i + 1) + "/" + nowTiff.getNSlices());	
	
	origTiff.setPosition(i+1);
	ImageProcessor myIP = origTiff.getProcessor();		
	
	float[] myX = new float[preciseCC.xID[i].length];for (int j = 0 ; j < myX.length ; j++) myX[j] = preciseCC.xID[i][j];
	float[] myY = new float[preciseCC.yID[i].length];for (int j = 0 ; j < myY.length ; j++) myY[j] = preciseCC.yID[i][j];
	PolygonRoi iRoi = new PolygonRoi(myX,myY,Roi.POLYLINE);
	float[] moX = new float[preciseCC.xOD[i].length];for (int j = 0 ; j < moX.length ; j++) moX[j] = preciseCC.xOD[i][j];
	float[] moY = new float[preciseCC.yOD[i].length];for (int j = 0 ; j < moY.length ; j++) moY[j] = preciseCC.yOD[i][j];			
	PolygonRoi oRoi = new PolygonRoi(moX,moY,Roi.POLYLINE);
			
	//fit a spline to it..
	//iRoi.fitSpline();
	//oRoi.fitSpline();
	
	//trick imagej to fill the outside of oRoi
	ImageProcessor copyIP = myIP.duplicate();
	copyIP.setValue(0);
	copyIP.fill();
	
	copyIP.setRoi(oRoi);
	copyIP.setValue(1);
	copyIP.fill(oRoi);
	
	copyIP.setRoi(iRoi);
	copyIP.setValue(0);
	copyIP.fill(iRoi);
	
	//multiply both image processors
	for (int xx = 0 ; xx < myIP.getWidth() ; xx++) {
		for (int yy = 0 ; yy < myIP.getHeight() ; yy++) {
			if (copyIP.getPixel(xx, yy) == 0) myIP.putPixel(xx, yy, 0);
		}
	}
	
	//origTiff.updateAndDraw();			
	//origTiff.show();
		
	int[] myHist = myIP.getHistogram();
	float[] cHist = new float[myHist.length];	
	cHist[0] = 0;
	for (int j = 1 ; j < myHist.length ; j++) cHist[j] = cHist[j-1] + (float)myHist[j];
	for (int j = 1 ; j < myHist.length ; j++) cHist[j] = cHist[j] / cHist[cHist.length - 1];			
	
	//evaluate thickness of wall
	//double wallThickness = Math.sqrt((preciseCC.xOD[i][0] - preciseCC.xID[i][0]) * (preciseCC.xOD[i][0] -preciseCC.xID[i][0]) + (preciseCC.yOD[i][0] -preciseCC.yID[i][0]) * (preciseCC.yOD[i][0] -preciseCC.yID[i][0])); 			 
	
	for (int j = 1 ; j < myHist.length ; j++) if (cHist[j] > 0.5) {
		q1[i] = j;
		break;				
	}
	
	for (int j = 1 ; j < myHist.length ; j++) if (cHist[j] > 0.9) {
		wallQ[i] = j;
		break;			
	}
	
}

//filter samples with dark q1
//double mq1 = jMed.evaluate(q1);
//for (int j = 0 ; j < q1.length ; j++) if (q1[j] < 0.5 * mq1) wallQ[j] = 0;

//origTiff.updateAndDraw();			
//origTiff.show();

return wallQ;

}	*/


/*	public ColumnContainer findColumnsWalls3D(ColumnContainer colCon) {

//init column wall coordinates
int colHeight = colCon.nowTiff.getNSlices();
double[] xCenter = new double[colHeight];
double[] yCenter = new double[colHeight];
double[] zCenter = new double[colHeight];
double[] outerMinorRadius = new double[colHeight];
double[] outerMajorRadius = new double[colHeight];
double[] deltaMajorMinor =  new double[colHeight];
double[] theta = new double[colHeight];		
double[] outerR2 = new double[colHeight];
double[] ixCenter = new double[colHeight];
double[] iyCenter = new double[colHeight];			
double[] innerMinorRadius = new double[colHeight];
double[] innerMajorRadius = new double[colHeight];				
double[] itheta = new double[colHeight];		
double[] innerR2 = new double[colHeight];		
boolean[] columnIsAtThisDepth = new boolean[colHeight];		
double[] wallThickness = new double[colHeight];	
double[] CVofWallthickness = new double[colHeight];	
double[] ratioOuterInnerRadius = new double[colHeight];	*/

/*public ColumnContainer findPVCOrAluWalls(ImagePlus nowTiff, ColCoords3D prelimCC, MenuWaiter.ColumnFinderMenuReturn jCFS, int[] startIncrStop) {


	
//flag doubtworthy results
colCon = flagDoubtworthyResults(colCon);

//look if the detected top and bottom of the columns are likely correct		
colCon = findTopAndBottomOfSoilColumn(colCon);

//find bevel of column
if (jCFS.hasBevel == true) colCon = findColumnsBevel(colCon);
else colCon.preciseCC.bevelStartsAt = 99999;  // a very large number to prevent the program to search for the bevel

//check if the column had been spuriously detected below the bevel		
colCon = findTopAndBottomOfSoilColumn(colCon);

//impute layer outline on top if R2 < chosen confidence criterion
colCon = imputeMissingLayers(colCon);

//detect rubber band
//colCon = rubberBandDetection(colCon);

return colCon;

}*/
