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

import java.io.File;

import org.apache.commons.math3.stat.StatUtils;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

/** 
 * RoiHandler is a SoilJ class containing all subroutines for the creation of 
 * regions of interests (ROI) overlays or overlay stacks.
 * 
 * @author John Koestel
 *
 */

public class RoiHandler implements PlugIn {

	
	public void run(String arg) {
				//ok, this is not needed..
	}
	
	public class ColumnRoi {
		
		public double area;
		public ObjectDetector.ColCoords3D jCO;
		public PolygonRoi[] pRoi;
		public PolygonRoi[] iRoi;
		public ImagePlus nowTiff;
		public ImagePlus surfaceNotCut;
		
	}
	
	public ObjectDetector.ColCoords3D scaleColumnCoordinates(double[] iDims, ObjectDetector.ColCoords3D inCo, MenuWaiter.PoreSpaceAnalyzerOptions mPSA, double scalingFactor, double[] oddVoxelContribution) {
		
		ObjectDetector jOD = new ObjectDetector();
		ObjectDetector.ColCoords3D outCo = jOD.new ColCoords3D();
		
		int cFT = mPSA.mRSO.cutAwayFromTop;
		int iHeight = mPSA.mRSO.heightOfRoi;
		
		//check if height of ROI is definedin percent of original height
		if (iHeight == 0) {
			if (mPSA.mRSO.cutZPercent == true) {
				if ((mPSA.mRSO.cutAwayFromTop == 50 & mPSA.mRSO.cutAwayFromBottom == 0) | (mPSA.mRSO.cutAwayFromTop == 0 & mPSA.mRSO.cutAwayFromBottom == 50)) {							
					iHeight = inCo.xmid.length / 2;
					if (mPSA.mRSO.cutAwayFromTop == 50) cFT = inCo.xmid.length / 2 - 1; 
				}			
				else {
					iHeight = inCo.xmid.length - (mPSA.mRSO.cutAwayFromTop + mPSA.mRSO.cutAwayFromBottom) * inCo.xmid.length / 100;
					cFT = inCo.xmid.length - mPSA.mRSO.cutAwayFromTop * inCo.xmid.length / 100;
				}
			}
			else {
				iHeight = inCo.xmid.length - mPSA.mRSO.cutAwayFromTop - mPSA.mRSO.cutAwayFromBottom;
			}
		}
		
		//check if radius of ROI is defined in percent of original radius
		double cutAwayFromXY = mPSA.mRSO.cutAwayFromWall;
		if (mPSA.mRSO.cutXYPercent == true) {
			double myRadius = (StatUtils.mean(inCo.innerMajorRadius) + StatUtils.mean(inCo.innerMinorRadius)) / 2;
			cutAwayFromXY = mPSA.mRSO.cutAwayFromWall * myRadius / 100; 
		}
		
		double[] innerMajorRadius = new double[iHeight];
		double[] innerMinorRadius = new double[iHeight];
		double[] xmid = new double[iHeight]; 
	    double[] ymid = new double[iHeight];
		double[] ixmid = new double[iHeight];
		double[] iymid = new double[iHeight];
		double[] theta = new double[iHeight];	
		double[] itheta = new double[iHeight]; 
		
		for (int i = 0 ; i < iHeight ; i ++) {
			
			innerMajorRadius[i] = scalingFactor * (inCo.innerMajorRadius[cFT + i] - cutAwayFromXY) - oddVoxelContribution[0];
			innerMinorRadius[i] = scalingFactor * (inCo.innerMinorRadius[cFT + i] - cutAwayFromXY) - oddVoxelContribution[1];
			xmid[i] = iDims[0] / 2 - oddVoxelContribution[0];
			ymid[i] = iDims[1] / 2 - oddVoxelContribution[0];
			ixmid[i] = iDims[0] / 2 - oddVoxelContribution[0];
			iymid[i] = iDims[1] / 2 - oddVoxelContribution[0];
			theta[i] = inCo.theta[cFT + i]; 	
			itheta[i] = inCo.itheta[cFT + i]; 	
			
		}
		
		outCo.innerMajorRadius = innerMajorRadius;
		outCo.innerMinorRadius = innerMinorRadius;
		outCo.xmid = xmid;
		outCo.ymid = ymid;
		outCo.ixmid = ixmid;	
		outCo.iymid = iymid;
		outCo.theta = theta;
		outCo.itheta = itheta; 
		outCo.topOfColumn = 0;
		outCo.heightOfColumn = iHeight;
		
		return outCo;
		
	}
	
	public PolygonRoi makeRoiFromFittedEllipse(FitStuff.FittedEllipse fE) {
		
		TailoredMaths math = new TailoredMaths();
		
		int j, cc;
		int maxAlpha = 360;
		int dAlpha = 5;	
		double[] myAngles = new double[maxAlpha/dAlpha];
		float[] x = new float[maxAlpha/dAlpha];
		float[] y = new float[maxAlpha/dAlpha];
		double[][] xy = new double[maxAlpha/dAlpha][maxAlpha/dAlpha];
						
		cc = 0;
		for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {myAngles[cc] = angle; cc++;}
			
		//create xy coordinates of ellipse..
		xy = math.getXYOfEllipseFromAngle(myAngles, fE.xCenter, fE.yCenter, fE.majorRadius, fE.minorRadius, fE.theta);
		
//		if (dx >= 0) alpha = Math.acos(dy / nowRadius);				
//		if (dx < 0) alpha = 2 * Math.PI - Math.acos(dy / nowRadius);	
//		
//		//reconstruct coordinates of column wall at this angle..			
//		float dx = x - (float)jCO.xmid[i];
//		float dy = y - (float)jCO.ymid[i];
//		float nowRadius = (float)Math.sqrt((double)(dx*dx) + (double)(dy*dy));
		
		//cast it into two arrays..
		for (j = 0 ; j < myAngles.length ; j++) x[j] = (float)xy[j][0];
		for (j = 0 ; j < myAngles.length ; j++) y[j] = (float)xy[j][1]; 
					
		PolygonRoi pRoi = new PolygonRoi(x, y, Roi.POLYGON);		
		
		return pRoi;		
		
	}
	/**
	 * 
	 * @param innerOrOuter		[String]
	 * @param wideTightOrExact	[String]
	 * @param jCO				[Object.Detector.EggShapedColCoord3D]
	 * @param cutAwayFromWalls	[int]
	 * @return
	 */
	public PolygonRoi[] makeMeAPolygonRoiStack(String innerOrOuter, String wideTightOrExact, ObjectDetector.EggShapedColCoords3D jCO, int cutAwayFromWalls) {

		TailoredMaths math = new TailoredMaths();
		
		int maxAlpha = 360;
		int dAlpha = maxAlpha / jCO.anglesChecked;	
		double[] myAngles = new double[maxAlpha/dAlpha];
		
		int cc = 0;
		for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {myAngles[cc] = angle; cc++;}
		
		float[] x = new float[jCO.anglesChecked];
		float[] y = new float[jCO.anglesChecked];
		double clipper0 = 0;		
		double clipper = clipper0;
		
		PolygonRoi[] nRoi = new PolygonRoi[jCO.heightOfColumn];
	
		if (wideTightOrExact.contentEquals("wide")) clipper0 = 1;
		if (wideTightOrExact.contentEquals("tight")) clipper0 = -1;
		if (wideTightOrExact.contentEquals("exact")) clipper0 = 0;
		if (wideTightOrExact.contentEquals("manual")) clipper0 = cutAwayFromWalls;
		
		for (int i = 0 ; i < jCO.heightOfColumn ; i++) { 
			
			if (innerOrOuter.contentEquals("inner")) {
				clipper = clipper0;
				for (int j = 0 ; j < jCO.anglesChecked ; j++) x[j] = (float)jCO.xID[i][j];
				for (int j = 0 ; j < jCO.anglesChecked ; j++) y[j] = (float)jCO.yID[i][j];				
			} 
			if (innerOrOuter.contentEquals("outer")) {				
				clipper = -clipper0;
				for (int j = 0 ; j < jCO.anglesChecked ; j++) x[j] = (float)jCO.xOD[i][j];
				for (int j = 0 ; j < jCO.anglesChecked ; j++) y[j] = (float)jCO.yOD[i][j];
			}
			if (innerOrOuter.contentEquals("outerFromInner")) {			
				if (jCO.wallThickness[i] - 10 > 3) clipper = jCO.wallThickness[i] - clipper0 ;
				else clipper = 3;
				for (int j = 0 ; j < jCO.anglesChecked ; j++) x[j] = (float)jCO.xID[i][j];
				for (int j = 0 ; j < jCO.anglesChecked ; j++) y[j] = (float)jCO.yID[i][j];
			}
			
			//apply the clipper
			for (int j = 0 ; j < jCO.anglesChecked ; j++) {
				
				double nowAngle = myAngles[j];
				double xnow = x[j];
				double ynow = y[j];
				double xmid = jCO.xmid[i];
				double ymid = jCO.ymid[i];
				
				double dx = xnow - xmid;
				double dy = ynow - ymid;
				double r = Math.sqrt(dx*dx + dy*dy);
				
				double clipRel = (r + clipper) / r;
				
				//clip
				double xnew = xnow + (Math.cos(nowAngle) * clipRel);	
				double ynew = ynow + (Math.sin(nowAngle) * clipRel);
				
				x[j] = (float)xnew;
				y[j] = (float)ynew;
				
			}
			
			nRoi[i - jCO.topOfColumn] = new PolygonRoi(x, y, Roi.POLYGON);
		}
		
		return nRoi;
	}
	
	
	/**
	 * 
	 * @param innerOrOuter			[String]	cut away from inner or from outer wall?
	 * @param wideTightOrExact		[String]	
	 * @param preciseCC				[ObjectDetector.ColCoords3D]	The column outlines as in innerCircle files
	 * @param cutAwayFromWalls		[int]		voxels to cut away
	 * @return PolygonRoi[heigthOfColumn] 
	 */
	public PolygonRoi[] makeMeAPolygonRoiStack(String innerOrOuter, String wideTightOrExact, ObjectDetector.ColCoords3D preciseCC, int cutAwayFromWalls) {

		TailoredMaths math = new TailoredMaths();
		
		int cc;
		int maxAlpha = 360;		// one full round
		int dAlpha = 5;			//steps
		double[] myAngles = new double[maxAlpha/dAlpha];	// to store all angles
		float[] x = new float[maxAlpha/dAlpha];				// to store x from each angle
		float[] y = new float[maxAlpha/dAlpha];
		double[][] xy = new double[maxAlpha/dAlpha][2];		// to store x and y from each angle
		double majRad = 0;
		double minRad = 0;
		double clipper0 = -100;								// how much to cut away? depending on wide/tight/exact >> set on some random number(?) for now 
		double clipper = clipper0;
		
		PolygonRoi[] nRoi = new PolygonRoi[preciseCC.heightOfColumn];	// to store polygon with entry for each slice
		
		cc = 0;
		for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {myAngles[cc] = angle; cc++;}
		
		if (wideTightOrExact.contentEquals("wide")) clipper0 = 1;
		if (wideTightOrExact.contentEquals("tight")) clipper0 = -1;
		if (wideTightOrExact.contentEquals("exact")) clipper0 = 0;
		if (wideTightOrExact.contentEquals("manual")) clipper0 = cutAwayFromWalls;
		
		for (int i = 0 ; i < preciseCC.heightOfColumn ; i++) { // loop through column
			
			if (innerOrOuter.contentEquals("inner")) {
				majRad = preciseCC.innerMajorRadius[i];
				minRad = preciseCC.innerMinorRadius[i];
				clipper = clipper0;
			} // inner
			if (innerOrOuter.contentEquals("outer")) {
				majRad = preciseCC.outerMajorRadius[i];
				minRad = preciseCC.outerMinorRadius[i];
				clipper = -clipper0;
			} // outer
			if (innerOrOuter.contentEquals("outerFromInner")) {  // what is outerFromInner?
				majRad = preciseCC.innerMajorRadius[i];
				minRad = preciseCC.innerMinorRadius[i];
				if (preciseCC.wallThickness[i] - 10 > 3) clipper = preciseCC.wallThickness[i] - clipper0 ;
				else clipper = 3;
			} // outerFromInner (???)
			
			xy = math.getXYOfEllipseFromAngle(					// calculate the ROI outlines
					myAngles, preciseCC.xmid[i], 
					preciseCC.ymid[i], 
					majRad + clipper, 
					minRad + clipper, 
					preciseCC.theta[i]
			);
		
			for (int j = 0 ; j < myAngles.length ; j++) x[j] = (float)xy[j][0];		// in each slice the x coordinates of all angles
			for (int j = 0 ; j < myAngles.length ; j++) y[j] = (float)xy[j][1]; 
					
			nRoi[i - preciseCC.topOfColumn] = new PolygonRoi(x, y, Roi.POLYGON);	//  at height - top of column
		} // loop through column
		
		return nRoi;
	}
	
	
	/**
	 * 
	 * @param imageDimensions
	 * @param mRSO
	 * @return
	 */
	public PolygonRoi makeMeAnIndependentRoi(int[] imageDimensions, MenuWaiter.ROISelectionOptions mRSO) {

		TailoredMaths math = new TailoredMaths();
		
		int maxAlpha = 360;
		int dAlpha = 1;	
		double[] myAngles = new double[maxAlpha/dAlpha];		
		double[][] xy = new double[maxAlpha/dAlpha][2];
		PolygonRoi pRoi = null;
		
		//if maximal cylinder is selected, assign midpoint and radius
		if (mRSO.maximalCylinder) {
			
			mRSO.cylX = imageDimensions[0] / 2;
			mRSO.cylY = imageDimensions[1] / 2;
			
			mRSO.cylRadius = (imageDimensions[0] + imageDimensions[1]) / 4;
			
		}
		
		if (mRSO.choiceOfRoi.equals("Cuboid")) {	
			
			float x1 = mRSO.cubeX1;
			float y1 = mRSO.cubeY1;
			float x2 = mRSO.cubeX2;
			float y2 = mRSO.cubeY2;
			
			if (mRSO.cubeX2 == 0) {
				x2 = imageDimensions[0];
			}
			
			if (mRSO.cubeY2 == 0) {
				y2 = imageDimensions[1];
			}
			
			float[] x = new float[]{x1, x2, x2, x1};
			float[] y = new float[]{y1, y1, y2, y2};	
			
			pRoi = new PolygonRoi(x, y, Roi.POLYGON);
			
		} else if (mRSO.choiceOfRoi.equals("Cylinder")) {
			
			int cc = 0;
			for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {
				myAngles[cc] = angle; 
				cc++;
			}		
			xy = math.getXYOfEllipseFromAngle(myAngles, mRSO.cylX, mRSO.cylY, mRSO.cylRadius, mRSO.cylRadius, 0);
			
			float[] x = new float[maxAlpha/dAlpha];
			float[] y = new float[maxAlpha/dAlpha];
			
			for (int i = 0 ; i < myAngles.length ; i++) {
				x[i] = (float)xy[i][0];
				y[i] = (float)xy[i][1];
			}	
			
			pRoi = new PolygonRoi(x, y, Roi.POLYGON);
		}
		
		return pRoi;
	}
	
	public PolygonRoi[] makeMeASplineRoiStack(String innerOrOuter, String wideTightOrExact, ObjectDetector.ColCoords3D preciseCC, double cutAway) {

		TailoredMaths math = new TailoredMaths();
		
		int i, j, cc;
		int maxAlpha = 360;
		int dAlpha = 5;	
		double[] myAngles = new double[maxAlpha/dAlpha];
		float[] x = new float[maxAlpha/dAlpha];
		float[] y = new float[maxAlpha/dAlpha];
		double[][] xy = new double[maxAlpha/dAlpha][maxAlpha/dAlpha];	
		double majRad, minRad;
		double clipper = -100;
		
		PolygonRoi[] pRoi = new PolygonRoi[preciseCC.heightOfColumn];
		
		cc = 0;
		for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {myAngles[cc] = angle; cc++;}
		
		if (wideTightOrExact.contentEquals("wide")) clipper = 5;
		if (wideTightOrExact.contentEquals("tight")) clipper = -3;
		if (wideTightOrExact.contentEquals("exact")) clipper = 0;
		if (wideTightOrExact.contentEquals("manual")) clipper = -cutAway;
		
		for (i = 0 ; i < preciseCC.heightOfColumn ; i++) { 
			
			if (innerOrOuter.contentEquals("inner")) {
				majRad = preciseCC.innerMajorRadius[i];
				minRad = preciseCC.innerMinorRadius[i];
			} 
			else {
				majRad = preciseCC.outerMajorRadius[i];
				minRad = preciseCC.outerMinorRadius[i];
			}
			
			xy = math.getXYOfEllipseFromAngle(myAngles, preciseCC.xmid[i], preciseCC.ymid[i], majRad + clipper, minRad + clipper, preciseCC.theta[i]);
		
			for (j = 0 ; j < myAngles.length ; j++) x[j] = (float)xy[j][0];
			for (j = 0 ; j < myAngles.length ; j++) y[j] = (float)xy[j][1]; 
					
			pRoi[i] = new PolygonRoi(x, y, Roi.POLYGON);
		}
		
		return pRoi;
	}

	public PolygonRoi[] makeMeAPolygonRoiStack4RawImage(String innerOrOuter, String wideTightOrExact, ObjectDetector.ColCoords3D preciseCC, double cutAway) {

		//init units
		TailoredMaths math = new TailoredMaths();
		
		//init variables
		int i, j, cc;
		int maxAlpha = 360;
		int dAlpha = 5;	
		double[] myAngles = new double[maxAlpha/dAlpha];
		float[] x = new float[maxAlpha/dAlpha];
		float[] y = new float[maxAlpha/dAlpha];
		double[][] xy = new double[maxAlpha/dAlpha][maxAlpha/dAlpha];	
		double majRad, minRad;
		double clipper = -100;
		
		PolygonRoi[] pRoi = new PolygonRoi[preciseCC.heightOfColumn];
		
		cc = 0;
		for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {myAngles[cc] = angle; cc++;}
		
		if (wideTightOrExact.contentEquals("wide")) clipper = 1;
		if (wideTightOrExact.contentEquals("tight")) clipper = -1;
		if (wideTightOrExact.contentEquals("exact")) clipper = 0;
		if (wideTightOrExact.contentEquals("manual")) clipper = -cutAway;
		
		for (i = preciseCC.topOfColumn ; i < preciseCC.topOfColumn + preciseCC.heightOfColumn ; i++) { 
			
			if (innerOrOuter.contentEquals("inner")) {
				majRad = preciseCC.innerMajorRadius[i];
				minRad = preciseCC.innerMinorRadius[i];
			} 
			else {
				majRad = preciseCC.outerMajorRadius[i];
				minRad = preciseCC.outerMinorRadius[i];
			}
			
			xy = math.getXYOfEllipseFromAngle(myAngles, preciseCC.xmid[i], preciseCC.ymid[i], majRad + clipper, minRad + clipper, preciseCC.theta[i]);
		
			for (j = 0 ; j < myAngles.length ; j++) x[j] = (float)xy[j][0];
			for (j = 0 ; j < myAngles.length ; j++) y[j] = (float)xy[j][1]; 
					
			pRoi[i - preciseCC.topOfColumn] = new PolygonRoi(x, y, Roi.POLYGON);
		}
		
		return pRoi;
	}

	
	/**
	 * 
	 * @param mFC  					[InputOutput.MyFileCollection]
	 * @param nowTiff				[ImagePlus]
	 * @param mRSO					[MenuWater.ROISelectionOptions]
	 * @param imagePhase2BeAnalyzed	[String]
	 * @param nameOfAnalyzedPhase	[String]
	 * @return ColumnRoi() colRoi
	 */
	public ColumnRoi prepareDesiredRoi(InputOutput.MyFileCollection mFC, ImagePlus nowTiff, MenuWaiter.ROISelectionOptions mRSO, 
			String imagePhase2BeAnalyzed, String nameOfAnalyzedPhase) {
		
		String pathSep = "/";
		
		//init units
		InputOutput jIO = new InputOutput();
		ObjectDetector jOD = new ObjectDetector();
		ImageManipulator jIM = new ImageManipulator();
		HistogramStuff hist = new HistogramStuff();
		RoiHandler roi = new RoiHandler();
		RollerCaster rC = new RollerCaster();
		ColumnRoi colRoi = new ColumnRoi();				// contains ColCoords3D, PolygonRoi pRoi & iRoi, nowTiff, surfaceNotCut (ImagePlus)
		
		//init some important variables
		int startSlice = 0;
		int stopSlice = nowTiff.getNSlices();		
		ImagePlus[] outTiffs = new ImagePlus[2];		
		String imgName = mFC.colName;
		double area = 0;
		boolean hasSurface = mRSO.includeSurfaceTopography;
		
		//check if mFC and mRSO values match
		if (mRSO.cutAwayFromTop != mFC.startSlice) mFC.startSlice = mRSO.cutAwayFromTop;
		if (mRSO.cutAwayFromBottom > 0 & mRSO.cutAwayFromBottom != nowTiff.getNSlices() - mFC.stopSlice) mFC.stopSlice = nowTiff.getNSlices() - mRSO.cutAwayFromBottom;
		
		//binarize4 if necessary
		nowTiff = jIM.extractPhaseOfInterest(nowTiff, imagePhase2BeAnalyzed, nameOfAnalyzedPhase);	
			
		//and here we go..
		if (mRSO.choiceOfRoi.equals("RealSample")) {
			
			ImagePlus soilSurface = new ImagePlus();
						
			//find correct InnerCircle and Surface Files
			String[] GandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC);
			mFC.nowInnerCirclePath = mFC.myInnerCircleFolder + mFC.pathSep + GandS[0];
			
			//read InnerCircle file and create outline ROI
			if (GandS[0].contains("Gauge")) {
				
				ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();						// initiate new ColCoords3D for ROI	
				int versio = jIO.checkInnerCircleFileVersion(mFC.nowInnerCirclePath);		// 
				if (versio == 0) jCO = jIO.readInnerCircleVer0(mFC.nowInnerCirclePath);	
				else jCO = jIO.readInnerCircleVer1(mFC.nowInnerCirclePath);					// read inner circle file
				colRoi.jCO = jCO;															// 
				
				//check if cut away from wall is given in percent
				int cutAwayFromXY = mRSO.cutAwayFromWall;									// as selected by user in menu
				if (mRSO.cutXYPercent) {
					double myRadius = (StatUtils.mean(jCO.innerMajorRadius) + StatUtils.mean(jCO.innerMinorRadius)) / 2; // mean radius of major and minor inner circle
					cutAwayFromXY = (int) Math.floor(mRSO.cutAwayFromWall * myRadius / 100); // caluculate voxels to cut away >> why not inner & outer je zB 10% wegschneiden?
				}

				//calculate average area
				double[] avgRadius = new double[jCO.innerMajorRadius.length];				// array for storing avg radius in each slice
				for (int radiusFinder = 0 ; radiusFinder < jCO.innerMajorRadius.length ; radiusFinder++) {
					avgRadius[radiusFinder] = (jCO.innerMajorRadius[radiusFinder] + jCO.innerMinorRadius[radiusFinder]) / 2; // avg of major and minor inner radius
				}
				double myRadius = StatUtils.mean(avgRadius) - cutAwayFromXY;				// mean of the slice-averages
				area = myRadius * myRadius * Math.PI;										// > area
				
				int averageRadius = (int)Math.round(myRadius);								// int
				PolygonRoi[] pRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, -cutAwayFromXY);	// returns coordinates from all (chosen) angles for each slice
				PolygonRoi[] iRoi = null;
				if (cutAwayFromXY == 0) iRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 10);	// if cutAwayFromY is 0 then add 10???
				if (mRSO.cutAwayFromCenter > 0) {											// cutAwayFromCenter: number of voxels or percent
					int cutAwayFromCenter = mRSO.cutAwayFromCenter;
					if (mRSO.cutXYPercent) {
						cutAwayFromCenter = (int)Math.floor(averageRadius * mRSO.cutAwayFromCenter / 100);	// adapt to voxels if in percent
					}
					iRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, -averageRadius + cutAwayFromCenter); // -Radius+cutAwayFromCenter
				}
				
				colRoi.pRoi = pRoi;
				colRoi.iRoi = iRoi;
				
			}	// if (GandS[0].contains("Gauge"))		
			else {

				//read column coordinates
				ObjectDetector.EggShapedColCoords3D jCO = jOD.new EggShapedColCoords3D();						
				jCO = jIO.readInnerCircleSteel(mFC);
				colRoi.jCO = null;
				
				double myRadius = StatUtils.mean(jCO.innerRadius) - mRSO.cutAwayFromWall;			
				area = myRadius * myRadius * Math.PI;
				
				int averageRadius = (int)Math.round(myRadius);
				PolygonRoi[] pRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 4);		
				PolygonRoi[] iRoi = null;
				if (mRSO.cutAwayFromWall == 0) iRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 10);
				if (mRSO.cutAwayFromCenter > 0) iRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, -averageRadius);
				
				colRoi.pRoi = pRoi;
				colRoi.iRoi = iRoi;
			
			}	// else if (GandS[0].contains("Gauge"))		
					
			//load surfaces
			if (mRSO.includeSurfaceTopography) {
				String surfPath = mFC.mySurfaceFolder + pathSep + GandS[1];
				soilSurface = jIO.openTiff3D(surfPath);
			}
			
	// + + + + + + + + + + if the ROI is defined by the choice of layer below the surface
	// + + + + + + + + + + and the height of the ROI
			if ((mRSO.cutAwayFromTop > 0 & mRSO.cutAwayFromBottom == 0 ) & 
					mRSO.includeSurfaceTopography & mRSO.heightOfRoi > 0 & mRSO.heightOfRoi < 4000) {  // ROI below top and height of column
				
				//assign soil top and bottom surfaces		
				ImageStack surStack = soilSurface.getStack();						
				ImageProcessor botIP = surStack.getProcessor(2);
							
				int[] botSurfHist = botIP.getHistogram();
				int maxBotDepression = mFC.nOfSlices - hist.findMaxFromHistogram(botSurfHist);
				
				ImageStack[] tempStack = new ImageStack[2];
				ImageStack dimStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
				tempStack[0] = dimStack;tempStack[1] = dimStack;
			
				//assemble
				ImagePlus[] topHalf = jIM.addColumnRegion2Stack(nowTiff, 0, maxBotDepression, colRoi.pRoi, colRoi.iRoi, hasSurface);
							
				outTiffs = topHalf;
				//outTiffs[1] = null;
				
				//outTiffs[0].updateAndDraw();
				//outTiffs[0].show();				
				
				//also cut ROI
				for (int i = mFC.startSlice ; i < mFC.stopSlice - mFC.startSlice ; i++) {
					colRoi.pRoi[i - mFC.startSlice] = colRoi.pRoi[i];
					if (colRoi.iRoi != null) colRoi.iRoi[i - mFC.startSlice] = colRoi.iRoi[i];
				}				
			} // if ROI below top and height of column
			
			//if only bottom voxels need to be removed below surface
			if ((mRSO.cutAwayFromTop > 0 & mRSO.cutAwayFromBottom == 0 ) & mRSO.includeSurfaceTopography & (mRSO.heightOfRoi == 0 | mRSO.heightOfRoi > 4000)) {
				
				//assign soil top and bottom surfaces		
				ImageStack surStack = soilSurface.getStack();						
				ImageProcessor botIP = surStack.getProcessor(2);
							
				int[] botSurfHist = botIP.getHistogram();
				int maxBotDepression = mFC.nOfSlices - hist.findMaxFromHistogram(botSurfHist);
				
				ImageStack[] tempStack = new ImageStack[2];
				ImageStack dimStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
				tempStack[0] = dimStack;tempStack[1] = dimStack;
							
				//assemble
				ImagePlus[] topHalf = jIM.addColumnRegion2Stack(nowTiff, 0, maxBotDepression, colRoi.pRoi, colRoi.iRoi, hasSurface);
				ImagePlus[] botHalf = jIM.removePoresBelowSurface(nowTiff,maxBotDepression, botIP, colRoi.pRoi, colRoi.iRoi);
				
				outTiffs = jIM.stackImagePlusArrays(topHalf, botHalf);
				
				//also cut ROI
				for (int i = mFC.startSlice ; i < mFC.stopSlice - mFC.startSlice ; i++) {
					colRoi.pRoi[i - mFC.startSlice] = colRoi.pRoi[i];
					if (colRoi.iRoi != null) colRoi.iRoi[i - mFC.startSlice] = colRoi.iRoi[i];
				}					
			}
			
			//if only top voxels need to be removed below surface
			if ((mRSO.cutAwayFromTop == 0 & mRSO.cutAwayFromBottom > 0 ) & mRSO.includeSurfaceTopography & (mRSO.heightOfRoi == 0 | mRSO.heightOfRoi > 4000)) {
				
				//assign soil top and bottom surfaces		
				ImageStack surStack = soilSurface.getStack();
				ImageProcessor topIP = surStack.getProcessor(1);		
				
				//get max surface depressions				
				int[] topSurfHist = topIP.getHistogram();
				int maxTopDepression = mFC.nOfSlices - hist.findMaxFromHistogram(topSurfHist);			
				
				//extract image phase that is to be analyzed		
				
				
				//assemble
				ImagePlus[] topHalf = jIM.removePoresAboveSurface(nowTiff, maxTopDepression, topIP, colRoi.pRoi, colRoi.iRoi);
				ImagePlus[] botHalf = jIM.addColumnRegion2Stack(nowTiff,maxTopDepression, nowTiff.getNSlices(), colRoi.pRoi, colRoi.iRoi, hasSurface);
				
				outTiffs = jIM.stackImagePlusArrays(topHalf, botHalf);
				
				//also cut ROI
				for (int i = mFC.startSlice ; i < mFC.stopSlice - mFC.startSlice ; i++) {
					colRoi.pRoi[i - mFC.startSlice] = colRoi.pRoi[i];
					if (colRoi.iRoi != null) colRoi.iRoi[i - mFC.startSlice] = colRoi.iRoi[i];
				}		
				
			}
			
			//if both soil surfaces should be included
			if (mRSO.cutAwayFromTop == 0 & mRSO.includeSurfaceTopography & mRSO.cutAwayFromBottom == 0 & (mRSO.heightOfRoi == 0 | mRSO.heightOfRoi > 4000)) {
				
				//assign soil top and bottom surfaces		
				ImageStack surStack = soilSurface.getStack();
				ImageProcessor topIP = surStack.getProcessor(1);		
				ImageProcessor botIP = surStack.getProcessor(2);
				
				//get max surface depressions				
				int[] topSurfHist = topIP.getHistogram();
				int maxTopDepression = hist.findMaxFromHistogram(topSurfHist);			
				int[] botSurfHist = botIP.getHistogram();
				int maxBotDepression = mFC.nOfSlices - hist.findMaxFromHistogram(botSurfHist);
				
				//assemble
				ImagePlus[] topOfC = jIM.removePoresAboveSurface(nowTiff, maxTopDepression, topIP, colRoi.pRoi, colRoi.iRoi);
				ImagePlus[] middleOfC= jIM.addColumnRegion2Stack(nowTiff,maxTopDepression, maxBotDepression, colRoi.pRoi, colRoi.iRoi, hasSurface);
				ImagePlus[] botOfC = jIM.removePoresBelowSurface(nowTiff,maxBotDepression, botIP, colRoi.pRoi, colRoi.iRoi);
			
				outTiffs = jIM.stackImagePlusArrays(topOfC, middleOfC);
				outTiffs = jIM.stackImagePlusArrays(outTiffs, botOfC);		
			
			}		
			
			//if no surface topographie should be included
			if (!mRSO.includeSurfaceTopography) { 
					
				outTiffs[0] = nowTiff;			
				outTiffs[1] = null;
				
				//also cut ROI
				for (int i = mFC.startSlice ; i < mFC.stopSlice - mFC.startSlice ; i++) {
					colRoi.pRoi[i - mFC.startSlice] = colRoi.pRoi[i];
					if (colRoi.iRoi != null) colRoi.iRoi[i - mFC.startSlice] = colRoi.iRoi[i];
				}
				
			}
		
			//cut image in XY-plane ...  this is done in a rather complicated fashion because of the surface topography image. 
			int maxX = 0;	//initialize maxX with the minimally possible value so that the true value is larger
			int maxY = 0;
			int minX = nowTiff.getWidth();  //initialize minX with the maximally possible value so that the true value is larger
			int minY = nowTiff.getHeight();
			for (int i = 0 ; i < colRoi.pRoi.length ; i++) {
				int[] nowX = colRoi.pRoi[i].getXCoordinates();
				int[] nowY = colRoi.pRoi[i].getXCoordinates();
				
				double nmiX = StatUtils.min(rC.castInt2Double(nowX)) + colRoi.pRoi[i].getXBase();
				double nmiY = StatUtils.min(rC.castInt2Double(nowY)) + colRoi.pRoi[i].getYBase();
				double nmaX = StatUtils.max(rC.castInt2Double(nowX)) + colRoi.pRoi[i].getXBase();
				double nmaY = StatUtils.max(rC.castInt2Double(nowY)) + colRoi.pRoi[i].getYBase();
				
				if (nmiX < minX) minX = (int)Math.round(nmiX);
				if (nmiY < minY) minY = (int)Math.round(nmiY);
				if (nmaX > maxX) maxX = (int)Math.round(nmaX);
				if (nmaY > maxY) maxY = (int)Math.round(nmaY);
				
			}
			
			//catch impossible values..
			if (minX < 0) minX = 0;
			if (minY < 0) minY = 0;
			if (maxX >= nowTiff.getWidth()) maxX = nowTiff.getWidth() - 1;
			if (maxY >= nowTiff.getHeight()) maxY = nowTiff.getHeight() - 1;
		
			//do the cut.. this is done in a rather complicated fashion because of the surface topography image. 
			int[] imageDimensions = {maxX-minX, maxY-minY}; 
			mRSO.choiceOfRoi = "Cuboid";
			mRSO.cubeX1 = minX;
			mRSO.cubeY1 = minY;
			mRSO.cubeX2 = maxX;
			mRSO.cubeY2 = maxY;
			PolygonRoi cutRoi = makeMeAnIndependentRoi(imageDimensions, mRSO);			
			ImageStack cutStack0 = new ImageStack(maxX-minX, maxY-minY);
			ImageStack cutStack1 = new ImageStack(maxX-minX, maxY-minY);
			
			//nowTiff.updateAndDraw();nowTiff.show();			
		
			for (int i = 0 ; i < outTiffs[0].getNSlices() ; i++) {
				
				outTiffs[0].setPosition(i+1);
				ImageProcessor nowIP = outTiffs[0].getProcessor().duplicate();
				
				//remove voxels outside ROI
				nowIP.setRoi(colRoi.pRoi[i]);
				nowIP.setColor(0);
				nowIP.fillOutside(colRoi.pRoi[i]);
				
				//also remove voxels inside inner Roi, is selected
				if (mRSO.cutAwayFromCenter > 0) {
					nowIP.setRoi(colRoi.iRoi[i]);
					nowIP.setColor(0);
					nowIP.fill(colRoi.iRoi[i]);
				}
				
				nowIP.setRoi(cutRoi);
				ImageProcessor cutIP = nowIP.crop();
				
				cutStack0.addSlice(cutIP);
				
				if (outTiffs[1] != null) {
					
					outTiffs[1].setPosition(i+1);
					ImageProcessor now1IP = outTiffs[1].getProcessor();
					
					now1IP.setRoi(cutRoi);
					ImageProcessor cutIP1 = now1IP.crop();
					
					cutStack1.addSlice(cutIP1);					
				}
			}
			
			outTiffs[0].setStack(cutStack0);
			if (outTiffs[1] != null) outTiffs[1].setStack(cutStack1);
			
			//outTiffs[0].updateAndDraw();outTiffs[0].show();
			
			mRSO.choiceOfRoi = "RealSample";
			
			//if there is a surface File, also cut this one..
			if (mRSO.includeSurfaceTopography & (mRSO.cutAwayFromTop > 0 | mRSO.cutAwayFromBottom > 0)) {
				
				ImageStack cutStackS = new ImageStack(maxX-minX, maxY-minY);
					
				int roiEnds = 0;
				for (int i = 0 ; i < 2 ; i++) {
					
					soilSurface.setPosition(i + 1);					
					ImageProcessor nowIP = soilSurface.getProcessor().duplicate();
					
					int[] myHist = nowIP.getHistogram();					
					int surfMedian = hist.findMedianFromHistogram(myHist);									
					int roiStartsAt = surfMedian + mRSO.cutAwayFromTop;		
					
					if (i == 0)	roiEnds = roiStartsAt + mRSO.heightOfRoi;
					
					if (i == 1) {
						roiStartsAt = mFC.nOfSlices - roiEnds;
					}
					
					nowIP.subtract(roiStartsAt);
					nowIP.add(1);  //make sure that the new surface is not 0. 
					
					nowIP.setRoi(cutRoi);
					ImageProcessor cutIP = nowIP.crop();
					
					//ImagePlus test = new ImagePlus("",cutIP);
					//test.updateAndDraw();
					//test.show();
					
					cutStackS.addSlice(cutIP);
				}
				
				ImagePlus cutSurfacefile = new ImagePlus();
				cutSurfacefile.setStack(cutStackS);
						
				//save it
				jIO.tiffSaver(mFC.myCutSurfaceFolder, mFC.fileName, cutSurfacefile);				
			}
			
		}
		
		if (!mRSO.choiceOfRoi.equals("RealSample")) {
			
			startSlice = 0;
			stopSlice = nowTiff.getNSlices();
			
			ImagePlus outTiff = new ImagePlus();
			
			//if maximal cylinder is selected, assign midpoint and radius
			if (mRSO.maximalCylinder) {
				
				mRSO.cylX = nowTiff.getWidth() / 2;
				mRSO.cylY = nowTiff.getHeight() / 2;
				
				mRSO.cylRadius = (nowTiff.getWidth() + nowTiff.getHeight()) / 4;
				
			}
			
			if (mRSO.choiceOfRoi.equals("Cuboid")) {
				area = (mRSO.cubeX2 - mRSO.cubeX1) *  (mRSO.cubeY2 - mRSO.cubeY1);
			}
			
			if (mRSO.choiceOfRoi.equals("Cylinder")) {
				area = mRSO.cylRadius * mRSO.cylRadius * Math.PI;
			}
			
			if (mRSO.choiceOfRoi.equals("TopOfEveryThing") | mRSO.choiceOfRoi.equals("BottomOfEveryThing")) {
				
				if (mRSO.areaOfInterest == 0) area = nowTiff.getWidth() * nowTiff.getHeight();
				else area = mRSO.areaOfInterest;
				
				startSlice = 0;
				stopSlice = nowTiff.getNSlices();
				
				if (mRSO.choiceOfRoi.equals("TopOfEveryThing")) stopSlice = (int)Math.round(stopSlice / 2);
				if (mRSO.choiceOfRoi.equals("BottomOfEveryThing")) startSlice = (int)Math.round(stopSlice / 2);
				
			}
			
			if (mRSO.choiceOfRoi.equals("Everything!")) {
				
				if (mRSO.areaOfInterest == 0) area = nowTiff.getWidth() * nowTiff.getHeight();
				else area = mRSO.areaOfInterest;
				
				outTiff = nowTiff.duplicate();
				
			}
			
			int[] imageDimensions = {nowTiff.getWidth(), nowTiff.getHeight()};
			
			PolygonRoi pRoi = makeMeAnIndependentRoi(imageDimensions, mRSO);
			PolygonRoi[] pRoiOut = new PolygonRoi[stopSlice - startSlice];
			
			ImageStack roiStack = new ImageStack(imageDimensions[0], imageDimensions[1]);
			if (mRSO.cutCanvas) roiStack = new ImageStack(pRoi.getBounds().width, pRoi.getBounds().height);
			
			//IJ.showMessage(imageDimensions[0] + " --- " + imageDimensions[1]);
			//IJ.showMessage(mPSA.cubeX1 + " --- " + mPSA.cubeX2);
			//IJ.showMessage(mPSA.cubeY1 + " --- " + mPSA.cubeY2);
			
			for (int i = startSlice ; i < stopSlice ; i++) {
				
				nowTiff.setPosition(i + 1);
				IJ.showStatus("Compiling cut-out slice #" + (i + 1 - startSlice));
				
				ImageProcessor modIP = nowTiff.getProcessor().duplicate();
				modIP.setRoi(pRoi);
			
				if (!mRSO.choiceOfRoi.equals("TopOfEveryThing") | !mRSO.choiceOfRoi.equals("BottomOfEveryThing") | !mRSO.choiceOfRoi.equals("EveryThing!")) {					
					modIP.setColor(0);
					modIP.fillOutside(pRoi);
				}
				
				if (mRSO.cutCanvas) {
					ImageProcessor cutIP = modIP.crop();
					roiStack.addSlice(cutIP);
				}
				else roiStack.addSlice(modIP);
				
				pRoiOut[i] = pRoi;
			}
			
			if (!mRSO.choiceOfRoi.equals("Everything!")) {
				outTiff.setStack(roiStack);
			}
			
			outTiffs[0] = outTiff;
			outTiffs[1] = null;
			
			//save ROI 
			colRoi.pRoi = pRoiOut;
		}		
		
		IJ.showStatus("Desired ROI has been compiled ...");

		//save the cut-out samples		
		if (mRSO.saveROI) {
			if (mRSO.cutAwayFromTop > 0 | mRSO.cutAwayFromWall > 0 | mRSO.choiceOfRoi.equals("RealSample")) {			
				jIO.tiffSaver(mFC.myPreOutFolder, imgName + ".tif", outTiffs[0]);						
			}
			if (mRSO.cutAwayFromTop == 0 & mRSO.choiceOfRoi.equals("RealSample") & !mRSO.includeSurfaceTopography) {
				jIO.tiffSaver(mFC.myPreOutFolder, imgName + ".tif", outTiffs[0]);						
			}
			if (mRSO.cutAwayFromTop == 0 & mRSO.choiceOfRoi.equals("RealSample") & mRSO.includeSurfaceTopography) {
				String thicknessTiffLocation = mFC.myPreOutFolder + pathSep + "Tiff4ThicknessAnalyses";
				new File(thicknessTiffLocation).mkdir();
				jIO.tiffSaver(thicknessTiffLocation, imgName + ".tif", outTiffs[1]);
			} 
			
			if (!mRSO.choiceOfRoi.equals("RealSample")) {
				jIO.tiffSaver(mFC.myPreOutFolder, imgName + ".tif", outTiffs[0]);
			}
		}
		
		//save also the cross-sectional area to allow for macroporosity calculations.. 
		mRSO.areaOfInterest = area;
		colRoi.area = area;
		String path4Area = mFC.myPreOutFolder + pathSep + "area.asc";
		jIO.writeStringIntoAsciiFile(path4Area, "" +  String.format("%8.0f\t",area));
	
		colRoi.nowTiff = outTiffs[0].duplicate();
		if (outTiffs[1] != null) {
			colRoi.surfaceNotCut = outTiffs[1].duplicate();
		}
		
		//outTiffs[0].updateAndDraw();
		//outTiffs[0].show();		
	
		return colRoi;
	}
	
	public ColumnRoi assembleInnerCircleROIs(InputOutput.MyFileCollection mFC, ImagePlus nowTiff, boolean useInnerCircle, int extraInnerCircle, boolean cutCanvas) {
		
		//init units
		InputOutput jIO = new InputOutput();
		ObjectDetector jOD = new ObjectDetector();
		ImageManipulator jIM = new ImageManipulator();
		RoiHandler roi = new RoiHandler();		
		ColumnRoi colRoi = new ColumnRoi();
		
		if (useInnerCircle) {		
			
			//find correct InnerCircle and Surface Files
			String[] GandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC);
			
			//read InnerCircle file and create outline ROI
			if (GandS[0].contains("Gauge")) {
				
				ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();				
				int versio = jIO.checkInnerCircleFileVersion(mFC.nowInnerCirclePath);			
				if (versio == 0) jCO = jIO.readInnerCircleVer0(mFC.nowInnerCirclePath);	
				else jCO = jIO.readInnerCircleVer1(mFC.nowInnerCirclePath);	
				colRoi.jCO = jCO;

				//calculate average area
				double[] avgRadius = new double[jCO.innerMajorRadius.length];
				for (int radiusFinder = 0 ; radiusFinder < jCO.innerMajorRadius.length ; radiusFinder++) {
					avgRadius[radiusFinder] = (jCO.innerMajorRadius[radiusFinder] + jCO.innerMinorRadius[radiusFinder]) / 2;
				}
				double myRadius = StatUtils.mean(avgRadius) - extraInnerCircle;			
				colRoi.area = myRadius * myRadius * Math.PI;
				
				int averageRadius = (int)Math.round(myRadius);
				PolygonRoi[] pRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, -extraInnerCircle);
			
				colRoi.pRoi = pRoi;
				
			}			
			else {

				//read column coordinates
				ObjectDetector.EggShapedColCoords3D jCO = jOD.new EggShapedColCoords3D();						
				jCO = jIO.readInnerCircleSteel(mFC);
				colRoi.jCO = null;
				
				double myRadius = StatUtils.mean(jCO.innerRadius) - extraInnerCircle;			
				colRoi.area = myRadius * myRadius * Math.PI;
				
				int averageRadius = (int)Math.round(myRadius);
				PolygonRoi[] pRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, -extraInnerCircle);		
				
				colRoi.pRoi = pRoi;
			
			}
					
			//Cut out image
			ImagePlus outTiff = jIM.cutImageInXYPlane(nowTiff, colRoi.pRoi, cutCanvas);	
			colRoi.nowTiff = outTiff;	
			
		}		
		else {
			
			colRoi.nowTiff = nowTiff;
			colRoi.area = nowTiff.getWidth() * nowTiff.getHeight();
			colRoi.pRoi = null;
			colRoi.jCO = null;
			
		}
		
		return colRoi;
	}
}