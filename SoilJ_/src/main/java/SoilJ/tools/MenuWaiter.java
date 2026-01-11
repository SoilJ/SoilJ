package SoilJ.tools;

import java.awt.Font;
import java.io.File;
import java.io.Serializable;

import org.apache.commons.math3.stat.StatUtils;

import SoilJ.tools.ObjectDetector.ApproximateColumnIllumination;

/**
 *SoilJ.tools is a collection of classes for SoilJ,
 *a collection of ImageJ plugins for the semi-automatized processing of 3-D X-ray images of soil columns
 *Copyright 2014 2015 2016 2017 2018 John Koestel
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
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.AutoThresholder;
import ij.process.AutoThresholder.Method;
import inra.ijpb.binary.ChamferWeights3D;

/**
 * MenuWaiter is a SoilJ class in which all pop-up menus are collected that are featured in SoilJ
 *
 * @author John Koestel
 *
 */

public class MenuWaiter implements PlugIn {


	public class BeamDeHardeningReturn {

		public boolean isSteelColumn = true;
		public int approximateRadius = 700;
		public int anglesChecked = 72;
		public double upperQ = 0.6;
		public double lowerQ = 0.8;
		public double correctionFunctionMedianFilter = 20;
		public String maskingMethod;

		double radiusOfMedianFilter = 2;
		double stepsize = 0.05;
		double cutoff = 0.2;
		int skip = 2;
		int maxEval = 200;
		double goodnessCriterion = 0.9;
		double severeGoodnessCriterion = 0.98;
		double gammaGoodnessCriterion = 0.7;
		double airPhaseClassMemberMinimum = 2;

	}
	
	public class LorenzosOptions {
		
		public String[] myTiffs;
		public String myInFolder;
		public String myAboveFolder;
		public String myOutFolder;
		public String myName;
		public int stepsize;
		public int average;
		public boolean runInBatchMode; 
		
	}
	
	public class ElsasOptions {
		
		public String[] myTiffs;
		public String myInFolder;
		public String myOutFolder;
		public String myTargetImage;
		public String myTransformationFile;
		
	}
	
	public class BioPoreExtractionOptions {
		
		public boolean doNotProcessOriginalResolution = false;
		public double thresholdVesselness; 
		public double smallesAllowedElongation;
		public int maximumBlurring;
		
	}
	
	public class RootExtractionOptions {
		
		public int minThreshold = 0;
		public int maxThreshold = 0;		
		public boolean doNotProcessOriginalResolution = false;
		public double thresholdVesselness; 
		public double smallesAllowedElongation;
		public int maximumBlurring;
		
	}
	
	public class CrackExtractionOptions {
		
		public boolean doNotProcessOriginalResolution = false;
		public double thresholdPlanarity;
		public double smallesAllowedElongation;
		public int maximumBlurring;
				
	}
		
	public class GravelExtractionOptions {
		
		public double voxelsize;		
		public int threshold;		
		public boolean filterOutSandFraction;
		public double medianFilter; 
		
	}

	public class ChiCalculatorReturn {

		public int wholeColumn = 0;
		public Boolean allClusters = false;
		public int clustersTakenIntoAccount = 101;
		public int minClusterSize = 5000;

	}

	public class MedianFilterAndUnsharpMaskReturn {

		public float medianFilterSizeXDir = 2;
		public float medianFilterSizeYDir = 2;
		public float medianFilterSizeZDir = 2;
		public double uMaskStandardDeviationXDir = 2;
		public double uMaskStandardDeviationYDir = 2;
		public double uMaskStandardDeviationZDir = 2;
		public float uMaskSharpeningWeight = 0.6f;

	}

	public class SelectFiles {
		
		public boolean selectChosen;
		
	}
	
	public class REVAnalyzerOptions {

		public String choiceOfRoi;
		public String choiceOfMethod;

		public int cubeX1;
		public int cubeX2;
		public int cubeY1;
		public int cubeY2;
		public int cubeZ1;
		public int cubeZ2;

		public int edgeX;
		public int edgeY;
		public int edgeZ;

		public int startDivNumber;
		public int stopDivNumber;
		
		public int stepNumber;
		public int stepLength;
		
		public int moveEdgeX;
		public int moveEdgeY;
		public int moveEdgeZ;
		public int movesX;
		public int movesY;
		public int movesZ;

		public int cylX;
		public int cylY;
		public int cylZ1;
		public int cylZ2;
		public int cylRadius;

		public boolean globVolume;
		public boolean globThickness;
		public boolean calcFractal;
		public boolean globAnisotropy;

		public boolean performParticleAnalyses;

		public boolean calcVolume;
		public boolean calcEuler;
		public boolean calcThickness;
		public boolean calcCriticalPoreDiameter;
		public boolean calcAnisotropy;
		public boolean calcPercolation;

		public boolean plotBinary;
		public boolean plotLabels;
		public boolean plotVolume;
		public boolean plotThickness;
		public boolean plotPercolation;

	}

	public class Extract2DHistogramOptions {
		
		public boolean calcGradientImage;
		
		public ROISelectionOptions mRSO;
		
	}
	
	public class ROISelectionOptions {
		
		public String choiceOfRoi;
		public String choiceOfZRoi;
		public String choiceOfXYRoi;

		public boolean cutZPercent; 
		public boolean cutXYPercent;		
		
		public int heightOfRoi;
		public int cutAwayFromTop;
		public int cutAwayFromBottom;
		
		public int cutAwayFromWall;
		public int cutAwayFromCenter;		
		public int keepRadiusFromCenter;	
		
		public boolean includeSurfaceTopography;
		public boolean useInnerCircleFiles;
		public boolean cutCanvas = false;
		 
		public String imagePhase2BeAnalyzed; 
		
		public int cubeX1;
		public int cubeX2;
		public int cubeY1;
		public int cubeY2;
		public int cubeZ1;
		public int cubeZ2;

		public int cylX;
		public int cylY;
		public int cylZ1;
		public int cylZ2;
		public int cylRadius;
		
		public double areaOfInterest;
		
		public boolean saveROI;
		
		public boolean maximalCylinder;
		
	}
	
	public class PoreSpaceAnalyzerOptions {
		
		public ROISelectionOptions mRSO;

		public int includeBreaks;
		
		public String imagePhase2BeAnalyzed = "255";
		public String nameOfAnalyzedPhase;		
		
		public boolean removeHoles;
		
		public boolean globVolume;
		public boolean globSurface;
		public boolean globThickness;
		public boolean calcCriticalPoreDiameter;
		public boolean calcChi;
		public boolean calcFractal;
		public boolean globAnisotropy;
		public boolean calcDistance;
		public boolean calcPercolatingVolume;
		public boolean calcVolCon2Top;

		public boolean performParticleAnalyses;

		public boolean calcVolume;
		public boolean calcSurface;
		public boolean calcMoments;
		public boolean calcUnitVectors;
		public boolean calcEuler;
		public boolean calcThickness;
		public boolean calcCorrLength;
		public boolean calcAnisotropy;
		public boolean calcSkeleton;
		public boolean calcInclination;
		public boolean calcPercolation;

		public boolean plotLabels;				
		public boolean plotThickness;
		public boolean plotPercolation;		
		public boolean plotDistanceMap;
		public boolean plotPoresConnected2Top;
		
		//skeletonization
		public boolean includeLoopGeneratorsOnTopAndBottom;
		
		public int numberOfSlicesAdded;
		public int numberOfBlankSlices2Add;
		public int numberOfCopiedSlices2Add;	

	}

	
	public class SurfaceFinderReturn {

		int neglectCrumbsOfLessThan;

	}

	public class MultiRegionSegmentation {		
		
		public int[] thresholds;
		public int numberOfThresholds;
		
		public boolean convert2EightBit;
		public int lowerBoundary;
		public int upperBoundary;
	
		public boolean calcHistograms;
		public boolean useInnerCircle;
		public boolean useJointHistogram;
		public boolean useAverageThresholds;
		
		public boolean calcOtsu;
		public boolean calcMaxEntro;
		public boolean calcMinErr;
		
		public int numberOfMethods;
		
	}
	
	public class ThresholderMenuReturn {

		public boolean useInnerCircle;

		public boolean filterImages = false;
		public String filterTag;
		public boolean setMaxgray2Wallgray;

		public Method myPrimaryMethod = AutoThresholder.Method.Default;
		public Method mySecondaryMethod = AutoThresholder.Method.Default;

		public int minThreshold = 0;
		public int maxThreshold = 0;

		public boolean useConstantThreshold;

		public boolean save3DImage;
		public boolean save4Evaluation;
		public boolean save4GeoDict;

	}
	
	public class POMThresholderMenuReturn {

		public int minThreshold = 0;
		public int maxThreshold = 0;
		public int minPOMVoxelNumber;
		public int openingSteps = 0;
		
		public int blackTopHatThreshold;
		public int gradientThreshold;
		
		public double windowSize;
		public double overlap;

		public boolean save3DImage;
		public boolean save4Evaluation;

	}

	public class OMFinderSettingsDEPRECATED {

		int mingrayValue;
		int maxgrayValue;
		double overlap;
		double windowSize;

	}

	public class OMFinderSettings {
	
		public int mingrayValue;
		public int maxgrayValue;
		public double overlap;
		public double windowSize;

	}
	
	public class RandomClusterGenerator {

		public final double pc =  0.0976d;

		public String shape;
		public int domainX;
		public int domainY;
		public int domainZ;
		public double[] porosityBounds;
		public double standardDeviation;
		public int numOfCopies;
		public String mode;
		public double[] porosityList;

	}
	
	public class AggregateMaskOptionMenu {
		
		public boolean useInnerCircleFiles;
		public int extraInnerCircle;
		public boolean cutCanvas;
		public int closingVoxels;
		public double erosionOvershoot;
		public double filterSize;
		
		//MorphoLibJ dist trans Watershed parameters
		public boolean normalize;
		public int dynamic;
		public ChamferWeights3D weights;
		
	}
	
	public class TobiasWatershedOptionMenu {
		
		public boolean useInnerCircleFiles;
		public int extraInnerCircle;
		public boolean cutCanvas;
		public boolean invertImage;
		public double filterSize;
		
		//MorphoLibJ dist trans Watershed parameters
		public boolean normalize;
		public int dynamic;
		public ChamferWeights3D weights;
		
	}
	
	public class GrayValues2Rescale {
		
		public int oldLow;
		public int oldHigh;
		
		public int newLow;
		public int newHigh;
		
	}
	
	public class ClosingMaskOptionMenu {
		
		public boolean useInnerCircleFiles;
		public int extraInnerCircle;
		public boolean cutCanvas;
		public int closingVoxels;
		public boolean noErosion;
		public double erosionOvershoot;
		public double filterSize;
		
	}
	

	public void run(String arg) {
				//ok, this is not needed..
	}
	
	public AggregateMaskOptionMenu showAggregateMaskOptionMenu() {
		
		GenericDialog gd = new GenericDialog("Aggregate mask generator menu");
		gd.addMessage("Aggregate mask generator created by Tobias Boelscher and John Koestel, Januar 2018\n");
		AggregateMaskOptionMenu aMO = new AggregateMaskOptionMenu();
		
		gd.addCheckbox("Do you want to use InnerCircle files?", false);
		gd.addNumericField("How many voxels should I cut away in addition to the inner circle? (irrelevant if above is unchecked)", 3, 0, 3, "vxs");
		gd.addCheckbox("Do you want to cut the canvas in the XY-plane accordingly?", false);
		
		gd.addNumericField("Morphological closing step number ", 20, 0, 3, "voxels");
		gd.addNumericField("Erosion overshoot factor ", 1.5, 1, 3, "");
		gd.addNumericField("Filter out particles smaller than ", 5, 2, 6, "% of the ROI size");
		
		Font newFont = new Font("SansSerif", Font.BOLD, 12);
		gd.addMessage("MorphoLibJ - Distance map options:", newFont);
		
		String standardWeight = ChamferWeights3D.WEIGHTS_3_4_5_7.toString();
		gd.addChoice("Distances", ChamferWeights3D.getAllLabels(), standardWeight);

		gd.addCheckbox("Normalize weights", false);
		
		gd.addMessage("Watershed options:", newFont);
		gd.addNumericField("Dynamic", 4, 2);
		
		gd.addHelp( "http://imagej.net/MorphoLibJ#Utilities_for_binary_images" );
		
		String myReference = "\nIf you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.\n\n";
		gd.addMessage(myReference);
		myReference = "Legland, D.; Arganda-Carreras, I. & Andrey, P. (2016), MorphoLibJ: integrated library and plugins for mathematical \n";
		myReference += "morphology with ImageJ, Bioinformatics (Oxford Univ Press) 32(22): 3532-3534, PMID 27412086, doi:10.1093/bioinformatics/btw413";
		gd.addMessage(myReference);

		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	aMO.useInnerCircleFiles = gd.getNextBoolean();	    	
	    	aMO.extraInnerCircle = (int)gd.getNextNumber();
	    	aMO.cutCanvas = gd.getNextBoolean();
	    	aMO.closingVoxels = (int)gd.getNextNumber();
	    	aMO.erosionOvershoot = gd.getNextNumber();
	    	aMO.filterSize = gd.getNextNumber();
	    	
			String weightLabel = gd.getNextChoice();			
			aMO.normalize = gd.getNextBoolean();
			aMO.dynamic = (int)gd.getNextNumber();
			
			// identify which weights should be used
			aMO.weights = ChamferWeights3D.fromLabel( weightLabel );
	    }
		
		return aMO;
	}
	
	public TobiasWatershedOptionMenu showTobiasWatershedOptionMenu() {
		
		GenericDialog gd = new GenericDialog("Tobias watershed menu");
		TobiasWatershedOptionMenu aMO = new TobiasWatershedOptionMenu();
		
		gd.addCheckbox("Do you want to use InnerCircle files?", false);
		gd.addNumericField("How many voxels should I cut away in addition to the inner circle? (irrelevant if above is unchecked)", 3, 0, 3, "vxs");
		gd.addCheckbox("Do you want to cut the canvas in the XY-plane accordingly?", false);
		
		gd.addCheckbox("Invert original Image? (check if air has the value 255)", false);		
		gd.addNumericField("Filter out particles smaller than ", 5, 2, 6, "% of the ROI size");
		
		Font newFont = new Font("SansSerif", Font.BOLD, 12);
		gd.addMessage("MorphoLibJ - Distance map options:", newFont);
		
		String standardWeight = ChamferWeights3D.WEIGHTS_3_4_5_7.toString();
		gd.addChoice("Distances", ChamferWeights3D.getAllLabels(), standardWeight);
		
		gd.addCheckbox("Normalize weights", false);
		
		gd.addMessage("Watershed options:", newFont);
		gd.addNumericField("Dynamic", 4, 2);
		
		gd.addHelp( "http://imagej.net/MorphoLibJ#Utilities_for_binary_images" );
		
		String myReference = "\nIf you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.\n\n";
		gd.addMessage(myReference);
		myReference = "Legland, D.; Arganda-Carreras, I. & Andrey, P. (2016), MorphoLibJ: integrated library and plugins for mathematical \n";
		myReference += "morphology with ImageJ, Bioinformatics (Oxford Univ Press) 32(22): 3532-3534, PMID 27412086, doi:10.1093/bioinformatics/btw413";
		gd.addMessage(myReference);

		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    
	    	aMO.useInnerCircleFiles = gd.getNextBoolean();
	    	aMO.extraInnerCircle = (int)gd.getNextNumber();
	    	aMO.cutCanvas = gd.getNextBoolean();
	    	
	    	aMO.invertImage = gd.getNextBoolean();
	    	aMO.filterSize = gd.getNextNumber();
	    	
			String weightLabel = gd.getNextChoice();			
			aMO.normalize = gd.getNextBoolean();
			aMO.dynamic = (int)gd.getNextNumber();
			
			// identify which weights should be used
			aMO.weights = ChamferWeights3D.fromLabel( weightLabel );
	    }
		
		return aMO;
	}
	
	public ClosingMaskOptionMenu showClosingMaskOptionMenu() {
		
		GenericDialog gd = new GenericDialog("Closing mask generator menu");
		gd.addMessage("Closing mask generator created by Tobias Boelscher and John Koestel, Januar 2018\n");
		ClosingMaskOptionMenu aMO = new ClosingMaskOptionMenu();
		
		gd.addCheckbox("Do you want to use InnerCircle files?", false);
		gd.addNumericField("How many voxels should I cut away in addition to the inner circle? (irrelevant if above is unchecked)", 3, 0, 3, "vxs");
		gd.addCheckbox("Do you want to cut the canvas in the XY-plane accordingly?", false);
		
		gd.addNumericField("Morphological closing step number ", 20, 0, 3, "voxels");
		gd.addCheckbox("No erosion?", false);
		gd.addNumericField("Erosion overshoot factor (irrelevant if 'no erosion' is checked)", 1.5, 1, 3, "");
		gd.addNumericField("Filter out particles smaller than ", 5, 2, 6, "% of the ROI size");
		
		String myReference = "\nIf you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.\n\n";		
		gd.addMessage(myReference);

		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	aMO.useInnerCircleFiles = gd.getNextBoolean();
	    	aMO.extraInnerCircle = (int)gd.getNextNumber();
	    	aMO.cutCanvas = gd.getNextBoolean();
	    	
	    	aMO.closingVoxels = (int)gd.getNextNumber();
	    	aMO.noErosion = gd.getNextBoolean();
	    	aMO.erosionOvershoot = gd.getNextNumber();
	    	aMO.filterSize = gd.getNextNumber();	    
	    }
		
		return aMO;
	}
	
	public class ClipXYMenu {
		
		public int clipAway;		
		
	}
	
	public class Divide2ComplexMRIImages {
		
		public ImagePlus w;
		public ImagePlus z;
		
		public int width;
		public int height;
		public int slices;
		
		boolean littleEndian;
		

	}
		
	/*public ImagePlus[] showDivide2ComplexMRIImages() {
		
			
		GenericDialog gd = new GenericDialog("Divide two complex MRI images");
			
		gd.add
		gd.addNumericField("Minimal considered POM volume per individual particle (set to 0 if watershed option is selected)", 10, 2, 6, "");
		gd.addCheckbox("Do you want to add the small POM to the pore space? (otherwise it is added to the matrix and solid phase)", false);
		
		gd.addNumericField("Specify the lower threshold in original image (obsolete if a POM region mask image region is used)", 9500, 0, 5, "");
		gd.addNumericField("Specify the upper threshold in original image (obsolete if a POM region mask image region is used)", 15000, 0, 5, "");
		gd.addNumericField("Specify threshold in gradient image (obsolete if a POM region mask image region is used)", 500, 0, 5, "");
		

		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	
	    	tS.usePOMRegionMask = gd.getNextBoolean();
	    	tS.useWatershed = gd.getNextBoolean();
	    	tS.minimalPOMVolume = gd.getNextNumber();
	    	tS.addSmallPOM2AirPhase = gd.getNextBoolean();
	    	tS.lowerThresh = (int)Math.round(gd.getNextNumber());
	    	tS.upperThresh = (int)Math.round(gd.getNextNumber());
	    	tS.gradientThresh = (int)Math.round(gd.getNextNumber());
	    	
	    }
		
		return tS;		
		
	}*/
	
	
	
	public class TernarySegmentation {
		
		public int lowerThresh;
		public int upperThresh;
		public int gradientThresh;
		public double minimalPOMVolume;
		public boolean usePOMRegionMask;
		public boolean useWatershed;
		public boolean addSmallPOM2AirPhase;
		
	}
	
	public MultiRegionSegmentation showMultiRegionSegmentationMenu() {
		
		MultiRegionSegmentation mRS = new MultiRegionSegmentation();
		
		GenericDialog gd = new GenericDialog("Multi region segmentation");
	
		gd.addCheckbox("I am using 16-bit images", true);
		gd.addNumericField("Please enter the lower boundary of the relevant grayscale range.", 5000, 0, 6, "");
		gd.addNumericField("Please enter the upper boundary of the relevant grayscale range.", 20000, 0, 6, "");
		
		gd.addCheckbox("I want to apply the multi region segmentation to the joint histogram of all my images!", true);
		gd.addMessage("Do not check if you have images with uncalibrated gray-scales!");
		
		gd.addNumericField("Please enter the number of regions you want to extract from the image", 3, 0, 3, "");
		
		String[] myAlgorithms = {"Otsu", "MaxEntropy", "MinError"};
		boolean[] myChoices = {true, true, true};
	
		gd.addMessage("Which thresholding algorithms do you want to apply?");
		gd.addCheckboxGroup(3, 1, myAlgorithms, myChoices);
		
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
 
	    	mRS.convert2EightBit = gd.getNextBoolean();
	    	mRS.lowerBoundary = (int)Math.round(gd.getNextNumber());
	    	mRS.upperBoundary = (int)Math.round(gd.getNextNumber());
	    	
	    	mRS.useJointHistogram = gd.getNextBoolean();	    	
	    	
	    	mRS.numberOfThresholds = (int)Math.round(gd.getNextNumber()) - 1;
	    	
	    	mRS.calcOtsu = gd.getNextBoolean();
	    	mRS.calcMaxEntro = gd.getNextBoolean();
	    	mRS.calcMinErr = gd.getNextBoolean();
	    }
		
		return mRS;
	}
	
	public ClipXYMenu showClipXYMenu() {
		
		GenericDialog gd = new GenericDialog("Clip image canvas in XY-plane relative to InnerCircle");
		gd.addMessage("This plugin requires InnerCircle files!!\n");
		
		ClipXYMenu myClip = new ClipXYMenu();
		
		gd.addNumericField("Specify margin around the inner perimeter that should be left (also negative numbers possible)", -3, 0, 8, "");

		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	myClip.clipAway = -(int)Math.round(gd.getNextNumber());  //get negative of typed in value
	    }
	    
	    return myClip;
	}
	
	public BeamDeHardeningReturn showBeamDeHardeningMenu() {

		BeamDeHardeningReturn mBDH = new BeamDeHardeningReturn();

		GenericDialog preg = new GenericDialog("Beam de-hardening menu");
		
		String[] typeOfColumn = new String[2];
		typeOfColumn[0] = "steel";
		typeOfColumn[1] = "aluminium or PVC";
		preg.addRadioButtonGroup("Please choose a column material!", typeOfColumn, 1, 2, "aluminium or PVC");
		
		preg.addNumericField("Enter the quantile for the lower reference value ", 0.6, 2, 5, "");
		preg.addNumericField("Enter the quantile for the upper reference value ", 0.8, 2, 5, "");
		
		preg.addMessage("Please remember that this routine requires InnerCircle files!!");

		preg.showDialog();
		if (preg.wasCanceled()) return null;
		else {
			String myChoice = preg.getNextRadioButton();
			
	    	int myChoiceIndex = 0;
	    	for (int i = 0 ; i < typeOfColumn.length ; i++) if (myChoice.equalsIgnoreCase(typeOfColumn[i])) {
	    		myChoiceIndex = i;
	    		break;
	    	}
	    	
	    	if (myChoiceIndex == 0) mBDH.isSteelColumn = true;
	    	else mBDH.isSteelColumn = false;	    	
	    	
	    	mBDH.lowerQ = preg.getNextNumber();	    	
	    	mBDH.upperQ = preg.getNextNumber();

		}
			
		return mBDH;

	}

	public ChiCalculatorReturn showChiCalculatorMenu() {

		GenericDialog gd = new GenericDialog("chi calculator menu");

		ChiCalculatorReturn mCCR = new ChiCalculatorReturn();

		String[] choiceOfROI = new String[2];
		choiceOfROI[0] = "sub-sample";
		choiceOfROI[1] = "whole column";
		choiceOfROI[2] = "random field";
		gd.addChoice("Please choose a region of interest!", choiceOfROI, "sub-sample");

		gd.addCheckbox("Do you want to consider all pore-clusters in the sample?", false);

		String[] choiceOfCutoff = new String[3];
		choiceOfCutoff[0] = "11 largest pore-clusters";
		choiceOfCutoff[1] = "101 largest pore-clusters";
		choiceOfCutoff[2] = "1001 largest pore-clusters";
		gd.addChoice("If not, how many do you want to consider?", choiceOfCutoff, "101 largest pore-clusters");

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	mCCR.wholeColumn = gd.getNextChoiceIndex();
	    	mCCR.allClusters = gd.getNextBoolean();
	    	int howManyIndex = gd.getNextChoiceIndex();
	    	switch (howManyIndex) {
	    		case 0: mCCR.clustersTakenIntoAccount = 11;break;
	    		case 1: mCCR.clustersTakenIntoAccount = 101;break;
	    		case 2: mCCR.clustersTakenIntoAccount = 1001;break;
	    	}
	    	if (mCCR.allClusters == true) mCCR.clustersTakenIntoAccount = 0;

	    }

		return mCCR;

	}

	public class ClipperMenuReturn {

		public boolean isSoilColumn = false;
		public boolean isCylinder = false;
		public boolean isRectangular = false;
		public boolean preserveOriginialCanvasSize = false;
		public boolean referenceIsSoilSurface = false;
		public boolean referenceIsTopMostSlice = false;
		public boolean addCanvasExccedance2Bottom = false;

		public int clipFromInnerPerimeter = 0;
		public int clipFromCanvasEdge = 0;
		public int canvasExceedsBy = 0;

		public int startAtSlice = 0;
		public int stopAtSlice = 0;
		public int heightOfROI = 0;

	}

	/*public ClipperMenuReturn showClipperDialog() {

		//construct objects
		GenericDialog gd = new GenericDialog("Clipper dialog");

		ClipperMenuReturn cMR = new ClipperMenuReturn();

		//build menu
		gd.addCheckbox("Do you want to use soil column outline coordinates saved as the Inner Circle?", true);

		String[] myChoices = new String[2];
		myChoices[0]="Cylindrical";myChoices[1]="Rectangular";
		gd.addRadioButtonGroup("What kind of shape do you want to cut out?", myChoices, 1, 2, myChoices[0]);

		gd.addCheckbox("Do you want to preserve the original canvas size?", false);

		gd.addNumericField("How many voxel do you want to cut away from the inner perimeter (canvas edge)?", 0, 0, 5, "");
		gd.addNumericField("By how many voxels should the canvas exceed the clipped image (if canvas size is not to be preserved)?", 2, 0, 3, "");
		gd.addCheckbox("Shall the blank canvas be also exceeded at the bottom of the column?", false);

		String[] myRefs = new String[2];
		myRefs[0]="The median soil surface as detected in SurfaceOfColumn";myRefs[1]="The topmost slice of the 3-D image";
		gd.addRadioButtonGroup("Choose a reference slice:", myRefs, 1, 2, myRefs[0]);

		gd.addNumericField("How many voxel below the reference depth is the upper boundary of your ROI?", 0, 0, 5, "");
		gd.addNumericField("How tall is your ROI in voxels?", 0, 0, 3, "");


		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	//get segmentation options
	    	cMR.isSoilColumn = gd.getNextBoolean();

	    	int shapeChoice = 0;
	    	String choice = gd.getNextRadioButton();
	    	for (int i = 0 ; i < myChoices.length ; i++) if (myChoices[i].equalsIgnoreCase(choice)) {
	    		shapeChoice = i;
	    		break;
	    	}
	    	switch (shapeChoice) {
				case 0 : {cMR.isCylinder = true; cMR.isRectangular = false; break;}
				case 1 : {cMR.isCylinder = false; cMR.isRectangular = true; break;}
	    	}

		    cMR.preserveOriginialCanvasSize = gd.getNextBoolean();

		    cMR.clipFromInnerPerimeter = (int)Math.round(gd.getNextNumber());
		    cMR.clipFromCanvasEdge = cMR.clipFromInnerPerimeter;
		    cMR.canvasExceedsBy = (int)Math.round(gd.getNextNumber());
		    cMR.addCanvasExccedance2Bottom = gd.getNextBoolean();

	    	int refChoice = 0;
	    	String rchoice = gd.getNextRadioButton();
	    	for (int i = 0 ; i < myRefs.length ; i++) if (myRefs[i].equalsIgnoreCase(rchoice)) {
	    		refChoice = i;
	    		break;
	    	}
	    	switch (refChoice) {
				case 0 : {cMR.referenceIsSoilSurface = true; cMR.referenceIsTopMostSlice = false; break;}
				case 1 : {cMR.referenceIsSoilSurface = false; cMR.referenceIsTopMostSlice = true; break;}
	    	}

		    cMR.startAtSlice = (int)Math.round(gd.getNextNumber());
		    cMR.heightOfROI = (int)Math.round(gd.getNextNumber());
		    cMR.stopAtSlice = cMR.startAtSlice + cMR.heightOfROI;

	    	return cMR;
	    }

	}*/

	public ThresholderMenuReturn showManualThresholdingDialog(ThresholderMenuReturn mTMR) {

		//construct objects
		GenericDialog gd = new GenericDialog("Thank you for choosing a user defined (constant) threshold for all images!");

		gd.addMessage("Make sure that you have run 'CalibrateGrayValues' when you apply this option!!!");

		gd.addNumericField("Enter the lower threshold for segmenting the image ", 0, 0, 5, "");

		gd.addNumericField("Enter the upper threshold for segmenting the image ", 0, 0, 5, "");
	
		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	mTMR.minThreshold = (int)Math.round(gd.getNextNumber());

	    	mTMR.maxThreshold = (int)Math.round(gd.getNextNumber());
	    	
	    	return mTMR;
	    }
	}

	public Boolean[] showFinalizerDialog() {

		//construct objects
		GenericDialog gd = new GenericDialog("");
		gd.addCheckbox("I want to merge two binaries", false);
		gd.addCheckbox("I want to cut out my sample using a binary mask", true);
		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	Boolean mergeOrNot = gd.getNextBoolean();
	    	Boolean cutOrNot = gd.getNextBoolean();
	    	Boolean[] finalizer = {mergeOrNot, cutOrNot};
	    	return finalizer;
	    }
	}

	public class Calc3DMenuReturn {

		public String operation;
		public String operationTag;
		public String filterTag;

		public boolean useInnerCircle;
		public boolean filterImages;

	}

	public Calc3DMenuReturn show3DCalcDialog(String A, String B) {

		//construct objects
		GenericDialog gd = new GenericDialog("3D Calculator dialog");
		Calc3DMenuReturn m3D = new Calc3DMenuReturn();

		//construct dialog window
		gd.addMessage(A);
		String[] items = {"Plus", "Minus"};

		gd.addChoice("", items, items[1]);

		gd.addMessage(B);

		gd.addMessage("");

		//add choices to menu
		gd.addCheckbox("Do you want to use the column outlines defined in 'inner circle' as a ROI?", false);

		gd.addMessage("Search-string for only segmenting some of the files in this folder. Leave blank for segmenting them all.");
		gd.addStringField("", "", 25);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	String myChoice = gd.getNextChoice();

	    	if (myChoice.equalsIgnoreCase(items[0])) {
	    		m3D.operation = "+";
	    		m3D.operationTag = "Plus";
	    	}

	    	if (myChoice.equalsIgnoreCase(items[1])) {
	    		m3D.operation = "-";
	    		m3D.operationTag = "Minus";
	    	}

	    	m3D.useInnerCircle = gd.getNextBoolean();

	    	m3D.filterTag = gd.getNextString();
	    	if (m3D.filterTag.hashCode() == 0) m3D.filterImages = false;
	    	else m3D.filterImages = true;

	      	return m3D;
	    }

	}
	
	

	public class HistogramMenuReturn {

		public boolean useInnerCircle;
		public boolean useSoilSurfaceFile;

		public boolean filterImages = false;
		public String filterTag;

		public boolean saveAsciis = true;
		public boolean saveImages = false;
		public boolean calcLumpedHistogram = false;

	}

	public HistogramMenuReturn showHistogramDialog() {

		//construct objects
		GenericDialog gd = new GenericDialog("HistoGrammar Dialog");

		HistogramMenuReturn hMR = new HistogramMenuReturn();

		//add choices to menu
		gd.addCheckbox("Do you want to use the column outlines defined in 'inner circle' as a ROI?", false);
		gd.addCheckbox("Do you want to use the soil surface topography?", false);

		/*gd.addMessage("Search-string for only segmenting some of the files in this folder. Leave blank for segmenting them all.");
		gd.addStringField("", "", 25);

		gd.addCheckbox("Save histograms as ASCII files?", true);

		gd.addCheckbox("Save histograms as images?", false);

		gd.addCheckbox("Calculate lumped histogram of all images in this folder?", true);
		gd.addMessage("(This option only makes sense if your images have normalized gray values, e.g. Hounsfield units)");
*/
		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	hMR.useInnerCircle = gd.getNextBoolean();
	    	hMR.useSoilSurfaceFile = gd.getNextBoolean();

	    	hMR.filterTag = gd.getNextString();
	    	if (hMR.filterTag.hashCode() == 0) hMR.filterImages = false;
	    	else hMR.filterImages = true;

	    	hMR.saveAsciis = true;//gd.getNextBoolean();
	    	hMR.saveImages = false;//gd.getNextBoolean();
	    	hMR.calcLumpedHistogram = true;//gd.getNextBoolean();

	      	return hMR;
	    }
	}

	public MedianFilterAndUnsharpMaskReturn showMedianAndUsharpMaskMenu() {

		GenericDialog gd = new GenericDialog("3-D median filter and 3-D unsharp mask");

		MedianFilterAndUnsharpMaskReturn mMU = new MedianFilterAndUnsharpMaskReturn();

		//set up menu
		gd.addNumericField("Median filter radius in X-direction ", mMU.medianFilterSizeXDir, 1, 6, "");
		gd.addNumericField("Median filter radius in Y-direction ", mMU.medianFilterSizeYDir, 1, 6, "");
		gd.addNumericField("Median filter radius in Z-direction ", mMU.medianFilterSizeZDir, 1, 6, "");
		gd.addNumericField("Standard deviation of the Gaussian kernel for the unsharp mask in X-direction ", mMU.uMaskStandardDeviationXDir, 1, 6, "");
		gd.addNumericField("Standard deviation of the Gaussian kernel for the unsharp mask in Y-direction ", mMU.uMaskStandardDeviationYDir, 1, 6, "");
		gd.addNumericField("Standard deviation of the Gaussian kernel for the unsharp mask in Z-direction ", mMU.uMaskStandardDeviationZDir, 1, 6, "");
		gd.addNumericField("Weight for unsharp mask application ", mMU.uMaskSharpeningWeight, 2, 6, "");

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	mMU.medianFilterSizeXDir = (float)gd.getNextNumber();
	    	mMU.medianFilterSizeYDir = (float)gd.getNextNumber();
	    	mMU.medianFilterSizeZDir = (float)gd.getNextNumber();

	    	mMU.uMaskStandardDeviationXDir = gd.getNextNumber();
	    	mMU.uMaskStandardDeviationYDir = gd.getNextNumber();
	    	mMU.uMaskStandardDeviationZDir = gd.getNextNumber();
	    	mMU.uMaskSharpeningWeight = (float)gd.getNextNumber();
	    }

		return mMU;

	}

	public class RingArtifactRemoverOptions {

		public boolean useInnerCircle;
		public int incrementBetweenSlices2BeCheckedForRings;

	}

	public RingArtifactRemoverOptions returnRingArtifactRemoverOptions() {

		GenericDialog gd = new GenericDialog("3-D median filter and 3-D unsharp mask");

		RingArtifactRemoverOptions rARO = new RingArtifactRemoverOptions();

		// set up dialog box
		gd.addCheckbox("Do you want to use the InnerCircle files?", false);

		gd.addNumericField("Increment between horizontal to be checked for ring artifacts", 10, 0);


		//show dialog
		gd.showDialog();
		if (gd.wasCanceled()) return null;
		else {

			rARO.useInnerCircle = gd.getNextBoolean();

			rARO.incrementBetweenSlices2BeCheckedForRings = (int)Math.round(gd.getNextNumber());
		}

		return rARO;

	}

	public class CalibrationReferences {

		public boolean useInnerCircle;
		public boolean useSoilSurfaceFile;

		public String material;
		public String lowRef;
		public String hiRef;

		public double lowerTarget;
		public double upperTarget;

		public double lowerReference;
		public double upperReference;

		public double divisorMinimumMode;

		public boolean sampleLowerWithinSoil;
		public boolean sampleUpperWithinSoil;

		public double[] originalLower;
		public double[] originalUpper;

		public String lowerTag;
		public String upperTag;
		
		public int windowSize;
		public String lowerSmoothingFilter;
		public String upperSmoothingFilter;		

	}

	public CalibrationReferences showCalibrationMenu() {

		MenuWaiter.CalibrationReferences mNR = new CalibrationReferences();

		//construct objects
		GenericDialog gd = new GenericDialog("Choose references for calibrating the gray-values!");

		// set up dialog box
		gd.addCheckbox("Do you want to use the InnerCircle files?", false);
	
		String[] columnMaterials = {"steel", "not-steel"};
		gd.addRadioButtonGroup("What kind of material is your column made of?", columnMaterials, 2, 1, columnMaterials[1]);

		gd.addMessage("");

		gd.addNumericField("lower normalization target value", 5000, 0, 6, "");
		gd.addNumericField("upper normalization target value", 20000, 0, 6, "");

		gd.addMessage("");

		String[] lowRef = {"quantile", "wall"};
		gd.addRadioButtonGroup("How do you want to define you lower reference gray value?", lowRef, 2, 1, lowRef[0]);

		String[] hiRef = {"quantile", "wall"};
		gd.addRadioButtonGroup("How do you want to define you upper reference gray value?", hiRef, 2, 1, hiRef[1]);

		String[] samplingLocation = {"inside the column", "outside the column, close to the wall", "outside the column, with distance to the wall"};

		String[] smoothingFilter = {"Median", "Mean", "Minimum", "Maximum"};
		
		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	mNR.useInnerCircle = gd.getNextBoolean();	    	

	    	mNR.material = gd.getNextRadioButton();
	    	if (mNR.material.equalsIgnoreCase("not-steel")) {
	    		mNR.material = "aluminium";
	    	}
	    	else mNR.material = "steel";

	    	mNR.lowerTarget = gd.getNextNumber();
	    	mNR.upperTarget = gd.getNextNumber();

	    	mNR.lowRef = gd.getNextRadioButton();
	    	mNR.hiRef = gd.getNextRadioButton();

	    }

	    if (mNR.lowRef.equalsIgnoreCase("quantile") | (mNR.hiRef.equalsIgnoreCase("quantile"))) {

	    	GenericDialog newGD = new GenericDialog("Please fill in missing input parameters!");

		    if (mNR.lowRef.equalsIgnoreCase("quantile")) {
		    	newGD.addNumericField("quantile for lower normalization reference value", 0.001, 5, 7, "");
		    	if (mNR.useInnerCircle) {
		    		newGD.addRadioButtonGroup("where do you want to sample the lower reference value?", samplingLocation, 3, 1, samplingLocation[0]);
		    	}
			}

		    if (mNR.hiRef.equalsIgnoreCase("quantile")) {
		    	newGD.addNumericField("quantile for upper normalization reference value", 0.001, 5, 7, "");
		    	if (mNR.useInnerCircle) {
		    		newGD.addRadioButtonGroup("where do you want to sample the upper reference value?", samplingLocation, 3, 1, samplingLocation[0]);
		    	}
			}

		    //add option for smoothing the reference values
		    newGD.addMessage("");
		    newGD.addNumericField("Window size for smoothing the reference value profiles", 15, 0);
		    
		    newGD.addRadioButtonGroup("Choose smoothing method for lower reference profile", smoothingFilter, 2, 2, smoothingFilter[0]);
		    newGD.addRadioButtonGroup("Choose smoothing method for upper reference profile", smoothingFilter, 2, 2, smoothingFilter[0]);
		      
		    //show dialog
		    String sampleLowerWithinSoil = "";
		    String sampleUpperWithinSoil = "";

			myReference = "If you are using this plugin please cite the following references: \n\n";
			newGD.setInsets(40, 0, 0);newGD.addMessage(myReference);
			myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
			myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
			newGD.setInsets(0, 0, 0);newGD.addMessage(myReference);

		    newGD.showDialog();

	  	    if (newGD.wasCanceled()) return null;
	  	    else {
	  	    	if (mNR.lowRef.equalsIgnoreCase("quantile")) {
	  	    		mNR.lowerReference = newGD.getNextNumber();
	  	    		
	  	    		String mylo = String.format("%1.4f", mNR.lowerReference);
	  		    	mNR.lowerTag = "Quantile" + mylo.substring(2, 5);
	  		    	
	  		    	mNR.sampleLowerWithinSoil = true;

	  		    	if (mNR.useInnerCircle) {
		  	    		
	  		    		sampleLowerWithinSoil = newGD.getNextRadioButton();
		  	    		
		  		    	if (sampleLowerWithinSoil == samplingLocation[0]) {
		  		    		mNR.sampleLowerWithinSoil = true;
		  		    		mNR.lowerTag += "Inside";
		  		    	}
		  		    	else {
		  		    		mNR.sampleLowerWithinSoil = false;
		  		    		if (sampleLowerWithinSoil == samplingLocation[1]) {	  		    	   		
		  		    	   		mNR.lowerTag += "Outside";
		  		    		}
		  		    		else {
		  		    			mNR.lowerTag += "FarOutside";
		  		    		}
		  		    	}	
	  	    		}
	  	    	}
  	    	
	  	    	if (mNR.hiRef.equalsIgnoreCase("quantile")) {
	  	    		mNR.upperReference = newGD.getNextNumber();
	  	    		
	  	    		String myup = String.format("%1.4f", mNR.upperReference);
	  		    	mNR.upperTag = "Quantile" + myup.substring(2, 5);
	  		    	
	  		    	mNR.sampleUpperWithinSoil = true;

	  		    	if (mNR.useInnerCircle) {
	  		    		sampleUpperWithinSoil = newGD.getNextRadioButton();
	  		    		  		    	
		  		    	if (sampleUpperWithinSoil == samplingLocation[0]) {
		  		    		mNR.sampleUpperWithinSoil = true;
		  		    		mNR.upperTag += "Inside";
		  		    	}
		  		    	else {
		  		    		mNR.sampleUpperWithinSoil = false;
		  		    		if (sampleUpperWithinSoil == samplingLocation[1]) {	  		    	   		
		  		    	   		mNR.upperTag += "Outside";
		  		    		}
		  		    		else {
		  		    			mNR.upperTag += "FarOutside";
		  		    		}
		  		    	}
	  		    	}
	  	    	}
	  	    	
	  	    	//read smoothing method
	  	    	int windowSize = (int)Math.round(newGD.getNextNumber() / 2);
	  	    	mNR.windowSize = 2 * windowSize - 1;
	  	    	
	  	    	mNR.lowerSmoothingFilter = newGD.getNextRadioButton();
	  	    	mNR.upperSmoothingFilter = newGD.getNextRadioButton();	  	    	
	  	    }
  	    }

  	    if (mNR.lowRef.equalsIgnoreCase("wall")) mNR.lowerTag = "Wall";

  	    if (mNR.hiRef.equalsIgnoreCase("wall")) mNR.upperTag = "Wall";
  
	    return mNR;

	}

	public OMFinderSettingsDEPRECATED showOMFinderMenuDEPRECATED() {

		OMFinderSettingsDEPRECATED oMF = new OMFinderSettingsDEPRECATED();

		GenericDialog gd = new GenericDialog("Fresh organic matter extractor dialog");

		gd.addNumericField("What is the minimum gray value of the organic matter? ", 9500, 0, 6, "");
		gd.addNumericField("What is the maximum gray value of the organic matter?", 15000, 0, 6, "");
		//gd.addNumericField("Window size relative to span between min and max gray values?", 0.33, 2, 4, "");
		//gd.addNumericField("Please choose an overlap (%)", 0.5, 2, 4, "");

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	oMF.mingrayValue = (int)Math.round(gd.getNextNumber());
	    	oMF.maxgrayValue = (int)Math.round(gd.getNextNumber());
	    	oMF.windowSize = 0.33f;//gd.getNextNumber();
	    	oMF.overlap = 0.5f;//gd.getNextNumber();

			return oMF;
	    }

	}
	
	public OMFinderSettings showOMFinderMenu() {

		OMFinderSettings oMF = new OMFinderSettings();

		GenericDialog gd = new GenericDialog("Menu for extracting a 2D histogram");

		gd.addCheckbox("I still need to calculate a gradient image", true);

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	return oMF;
	    }

	}
	
	public ROISelectionOptions regionOfInterestSelection() {

		//construct objects
		GenericDialog gd = new GenericDialog("Please tell me what you want to get analysed!");

		ROISelectionOptions mRSO = new ROISelectionOptions();
	    mRSO.useInnerCircleFiles = false;
	    mRSO.includeSurfaceTopography = false;
	    mRSO.saveROI = false;

		String[] choiceOfRoi = new String[7];
		choiceOfRoi[0] = "No ROI but the entire 3D image!";
		choiceOfRoi[1] = "The entire XY canvas, but cut away some slices in the vertical direction";
		choiceOfRoi[2] = "Cylindrical column with outlines defined under 'Inner Circle'";
		choiceOfRoi[3] = "Cylindrical ROI";
		choiceOfRoi[4] = "Cuboid ROI";
		choiceOfRoi[5] = "No ROI but the top half of the entire image!";
		choiceOfRoi[6] = "No ROI but the bottom half of the entire image!";
		gd.addRadioButtonGroup("Please choose a region of interest!", choiceOfRoi, 4, 1, choiceOfRoi[0]);

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);		
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	String myChoice = gd.getNextRadioButton();
	    	int myChoiceIndex = 0;
	    	for (int i = 0 ; i < choiceOfRoi.length ; i++) if (myChoice.equalsIgnoreCase(choiceOfRoi[i])) {
	    		myChoiceIndex = i;
	    		break;
	    	}	    	
	    	switch (myChoiceIndex) {
	    		case 0: mRSO.choiceOfRoi = "Everything!"; break;
	    		case 1: mRSO.choiceOfRoi = "EverythingInXY"; break;
	    		case 2: mRSO.choiceOfRoi = "RealSample"; break;
	    		case 3: mRSO.choiceOfRoi = "Cylinder"; break;
	    		case 4: mRSO.choiceOfRoi = "Cuboid"; break;
	    		case 5: mRSO.choiceOfRoi = "TopOfEveryThing"; break;
	    		case 6: mRSO.choiceOfRoi = "BottomOfEveryThing"; break;
	    	}
	    				
	    }

	    gd.removeAll();	 
	    
	    if (mRSO.choiceOfRoi.equalsIgnoreCase("EverythingInXY")) {
	    	
	    	gd.addNumericField("How many slices do you want to cut away from the top? ", 0, 0, 6, "");
			gd.addNumericField("How many slices do you want to cut away from the bottom? ", 0, 0, 6, "");

			myReference = "If you are using this plugin please cite the following references: \n\n";
			gd.setInsets(50, 0, 0);gd.addMessage(myReference);
			myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
			myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
			gd.setInsets(0, 0, 0);gd.addMessage(myReference);

	    	//show dialog
			gd.showDialog();
		    if (gd.wasCanceled()) return null;
		    else {
		    	//domain definitions		    	
		    	mRSO.cutAwayFromTop = (int)Math.round(gd.getNextNumber());
		    	mRSO.cutAwayFromBottom = (int)Math.round(gd.getNextNumber());
		    }
	    }

		if (mRSO.choiceOfRoi.equalsIgnoreCase("RealSample")) {
			
			GenericDialog gdR = new GenericDialog("Please tell me more");
			
			mRSO.useInnerCircleFiles = true;
			
			String[] choiceOfSubRoi = new String[5];
			choiceOfSubRoi[0] = "The entire soil column";
			choiceOfSubRoi[1] = "The central part of the soil column (define depth under soil surface and column height in voxel)";
			choiceOfSubRoi[2] = "The central part of the soil column (define number of voxels to be cut away from top and/or bottom)";
			choiceOfSubRoi[3] = "The central part of the soil column (define percentage to be cut away from top and/or bottom)";
			choiceOfSubRoi[4] = "The top or bottom half of the soil column";
			gdR.addRadioButtonGroup("Please choose a region of interest along the Z-axis", choiceOfSubRoi, 5, 1, "The entire soil column");

			String[] choiceOfSubRoi2 = new String[6];
			choiceOfSubRoi2[0] = "The entire soil column";
			choiceOfSubRoi2[1] = "The central part of the soil column (define radius around central axis in voxels)";
			choiceOfSubRoi2[2] = "The central part of the soil column (define number of voxels to be cut away from the wall)";
			choiceOfSubRoi2[3] = "The central part of the soil column (define percentage to be cut away from the wall)";
			choiceOfSubRoi2[4] = "A hollow cylinder cut out from the soil column (define number of voxels to be cut away from the center and/or the wall)";
			choiceOfSubRoi2[5] = "A hollow cylinder cut out from the soil column (define percentage to be cut away from the center and/or the wall)";
			gdR.addRadioButtonGroup("Please choose a region of interest in the XY-plane", choiceOfSubRoi2, 6, 1, "The entire soil column");
			
			String[] referencePoint = new String[2];
			referencePoint[0] = "The top and bottom soil surfaces as found by the SoilJ SoilSurfaceFinder";
			referencePoint[1] = "The top and bottom of the 3-D image canvas";
			gdR.addRadioButtonGroup("Please choose reference points for the top and bottom of your column", referencePoint, 2, 1, "The top and bottom of the 3-D image canvas");
		
			myReference = "If you are using this plugin please cite the following references: \n\n";
			gdR.setInsets(50, 0, 0);gdR.addMessage(myReference);
			myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
			myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
			gdR.setInsets(0, 0, 0);gdR.addMessage(myReference);

			//show dialog
			gdR.showDialog();
		    if (gdR.wasCanceled()) return null;
		    else {
		    	
		    	mRSO.cutZPercent = false;
    			mRSO.cutXYPercent = false;
    			
    			mRSO.cutAwayFromBottom = 0;
    			mRSO.cutAwayFromTop = 0;
    			mRSO.cutAwayFromCenter = 0;
    			mRSO.cutAwayFromWall = 0;
    			mRSO.heightOfRoi = 0;
		    	
		    	//domain definitions Z
		    	String myChoice = gdR.getNextRadioButton();
		    	int myChoiceIndex = 0;
		    	for (int i = 0 ; i < choiceOfSubRoi.length ; i++) if (myChoice.equalsIgnoreCase(choiceOfSubRoi[i])) {
		    		myChoiceIndex = i;
		    		break;
		    	}
		    
		    	switch (myChoiceIndex) {
		    		case 0: {
		    			mRSO.choiceOfZRoi = "everything";		    			
		    			break;
		    		}
		    		case 1: {
		    			
		    			mRSO.choiceOfZRoi = "fixedHeight";		    			
		    			GenericDialog gd3 = new GenericDialog("Please tell me more!");
		    			gd3.addNumericField("Cut away from top", 0, 0, 5, "voxel");
		    			gd3.addNumericField("Column height", 0, 0, 5, "voxel");
		    			
		    			gd3.showDialog();
		    		    if (gd3.wasCanceled()) return null;
		    		    else {
		    		    	mRSO.cutAwayFromTop = (int) Math.round(gd3.getNextNumber());
		    		    	mRSO.heightOfRoi = (int) Math.round(gd3.getNextNumber());
		    		    }
		    		
		    			break;
		    			
		    		}
		    		case 2: {
		    		
		    			mRSO.choiceOfZRoi= "relativeHeightVoxels"; 
		    		    			
		    			GenericDialog gd3 = new GenericDialog("Please tell me more!");
		    			gd3.addNumericField("Cut away from top", 0, 0, 5, "voxel");
		    			gd3.addNumericField("Cut away from bottom", 0, 0, 5, "voxel");
		    			
		    			gd3.showDialog();
		    		    if (gd3.wasCanceled()) return null;
		    		    else {
		    		    	mRSO.cutAwayFromTop = (int) Math.round(gd3.getNextNumber());
		    		    	mRSO.cutAwayFromBottom = (int) Math.round(gd3.getNextNumber());
		    		    	mRSO.heightOfRoi = 0;
		    		    }	    		    
		    			
		    			break;
		    		}
		    		case 3: {
		    			mRSO.choiceOfZRoi = "relativeHeightPercent";
			    	
		    			GenericDialog gd3 = new GenericDialog("Please tell me more!");
		    			gd3.addNumericField("Cut away from top", 0, 0, 5, "%");
		    			gd3.addNumericField("Cut away from bottom", 0, 0, 5, "%");
		    			
		    			gd3.showDialog();
		    		    if (gd3.wasCanceled()) return null;
		    		    else {
		    		    	mRSO.cutZPercent = true;
		    		    	mRSO.cutAwayFromTop = (int) Math.round(gd3.getNextNumber());
		    		    	mRSO.cutAwayFromBottom = (int) Math.round(gd3.getNextNumber());
		    		    	mRSO.heightOfRoi = 0;
		    		    }
		    			break;
		    		}
		    		case 4: {
		    			mRSO.choiceOfZRoi = "TopOrBottom";
		    			
		    			GenericDialog gd3 = new GenericDialog("Please tell me more!");
		    			
		    			String[] topBot= new String[2];
		    			topBot[0] = "Top half of column";
		    			topBot[1] = "Bottom half of column";
		    			gd3.addRadioButtonGroup("Please choose your preference", topBot, 2, 1, "Top half of column");
		    					    			
		    			gd3.showDialog();
		    		    if (gd3.wasCanceled()) return null;
		    		    else {
		    		    	myChoice = gd3.getNextRadioButton();
		    		    	myChoiceIndex = 0;
		    		    	for (int i = 0 ; i < topBot.length ; i++) if (myChoice.equalsIgnoreCase(topBot[i])) {
		    		    		myChoiceIndex = i;
		    		    		break;
		    		    	}
		    		    
		    		    	switch (myChoiceIndex) {
		    		    		case 0: {
		    		    			mRSO.cutZPercent = true;
				    		    	mRSO.cutAwayFromTop = 0;
				    		    	mRSO.cutAwayFromBottom = 50;		
				    		    	mRSO.heightOfRoi = 0;
		    		    			break;
		    		    		}
		    		    		case 1: {
		    		    			mRSO.cutZPercent = true;
				    		    	mRSO.cutAwayFromTop = 50;
				    		    	mRSO.cutAwayFromBottom = 0;	
				    		    	mRSO.heightOfRoi = 0;
		    		    			break;
		    		    		}
		    		    	}
		    		    }
		    			break;
		    		}
		    	}
		    	
		    	myChoice = gdR.getNextRadioButton();
		    	myChoiceIndex = 0;
		    	for (int i = 0 ; i < choiceOfSubRoi2.length ; i++) if (myChoice.equalsIgnoreCase(choiceOfSubRoi2[i])) {
		    		myChoiceIndex = i;
		    		break;
		    	}
		    
		    	switch (myChoiceIndex) {
		    		case 0: {
		    			mRSO.choiceOfXYRoi = "everything"; break;
		    		}
		    		case 1: {
		    			mRSO.choiceOfXYRoi = "centralRadiusVoxels"; 
		    			
		    			GenericDialog gd3 = new GenericDialog("Please tell me more!");
		    			gd3.addNumericField("Radius around the central axis", 0, 0, 5, "voxels");
		    			
		    			gd3.showDialog();
		    		    if (gd3.wasCanceled()) return null;
		    		    else {
		    		    	mRSO.keepRadiusFromCenter = (int) Math.round(gd3.getNextNumber());		    		    
		    		    }
		    			
		    			break;
		    		}
		    		case 2: {
		    			mRSO.choiceOfXYRoi = "centralVoxels"; 
		    			
		    			GenericDialog gd3 = new GenericDialog("Please tell me more!");
		    			gd3.addNumericField("Cut away from wall", 0, 0, 5, "voxels");
		    			
		    			gd3.showDialog();
		    		    if (gd3.wasCanceled()) return null;
		    		    else {
		    		    	mRSO.cutAwayFromWall = (int) Math.round(gd3.getNextNumber());		    		    
		    		    }
		    			
		    			break;
		    		}
		    		case 3: {		    			
		    			
		    			mRSO.choiceOfXYRoi= "centralPercent";	
		    			mRSO.cutXYPercent = true;
		    			
		    			GenericDialog gd3 = new GenericDialog("Please tell me more!");
		    			gd3.addNumericField("Cut away from wall", 0, 0, 5, "%");
		    			
		    			gd3.showDialog();
		    		    if (gd3.wasCanceled()) return null;
		    		    else {
		    		    	mRSO.cutXYPercent = true;
		    		    	mRSO.cutAwayFromWall = (int) Math.round(gd3.getNextNumber());		   
		    		    }
		    			
		    			break;
		    		}
		    		case 4: {
		    			mRSO.choiceOfXYRoi = "donutVoxels"; 
		    			
		    			GenericDialog gd3 = new GenericDialog("Please tell me more!");
		    			gd3.addNumericField("Cut away from wall", 0, 0, 5, "voxels");
		    			gd3.addNumericField("Cut away from center", 0, 0, 5, "voxels");	
		    			
		    			gd3.showDialog();
		    		    if (gd3.wasCanceled()) return null;
		    		    else {
		    		    	mRSO.cutAwayFromWall = (int) Math.round(gd3.getNextNumber());
		    		    	mRSO.cutAwayFromCenter = (int) Math.round(gd3.getNextNumber());
		    		    }
		    			
		    			break;
		    		}
		    		case 5: {
		    			mRSO.choiceOfXYRoi = "donutPercent"; 
		    			
		    			mRSO.cutXYPercent = true;
		    			
		    			GenericDialog gd3 = new GenericDialog("Please tell me more!");
		    			gd3.addNumericField("Cut away from wall", 0, 0, 5, "%");
		    			gd3.addNumericField("Cut away from center", 0, 0, 5, "%");			    		
		    			
		    			gd3.showDialog();
		    		    if (gd3.wasCanceled()) return null;
		    		    else {
		    		    	mRSO.cutXYPercent = true;
		    		    	mRSO.cutAwayFromWall = (int) Math.round(gd3.getNextNumber());
		    		    	mRSO.cutAwayFromCenter = (int) Math.round(gd3.getNextNumber());
		    		    }
		    			
		    			break;		    		
		    		}
		    	}
		   
		    	myChoice = gdR.getNextRadioButton();
		    	myChoiceIndex = 0;
		    	for (int i = 0 ; i < referencePoint.length ; i++) if (myChoice.equalsIgnoreCase(referencePoint[i])) {
		    		myChoiceIndex = i;
		    		break;
		    	}
		    
		    	switch (myChoiceIndex) {
		    		case 0: mRSO.includeSurfaceTopography = true;break;
		    		case 1: mRSO.includeSurfaceTopography = false;break;		    	    		
		    	}		
		    }
		    
		}

	    if (mRSO.choiceOfRoi.equalsIgnoreCase("Cylinder")) {
				    	
			gd.addNumericField("Z-coordinate of cylinder top ", 0, 0, 6, "");
			gd.addNumericField("Z-coordinate of cylinder bottom ", 0, 0, 6, "(leave 0 if you want to analyze all slices or if you want to cut away slices from the bottom)");
			gd.addMessage(" ");
			gd.addNumericField("Alternatively, specify how many slices you want to cut away from the bottom ", 0, 0, 6, "(leave 0 if you want to analyze all slices or specify a Z-coordinate)");
	    	
			gd.addMessage("      ");	    	
			
	    	gd.addCheckbox("Do you want to analyze the ROI corresponding to the maximally inscribable circle in the XY plane?", true);
	    	
	    	gd.addMessage("      ");	    	
	    	gd.addMessage("... or do you want to specify the coordinates of the ROI manually?");	    	
	    	gd.addMessage("(the coordinates below will be ignored if the checkbox is checked)");
	    	gd.addMessage("      ");
	    	
	    	gd.addNumericField("X-coordinate of cylinder center ", 0, 0, 6, "");
			gd.addNumericField("Y-coordinate of cylinder center ", 0, 0, 6, "");
			gd.addNumericField("Radius of cylinder ", 0, 0, 6, "");

			myReference = "If you are using this plugin please cite the following references: \n\n";
			gd.setInsets(50, 0, 0);gd.addMessage(myReference);
			myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
			myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
			gd.setInsets(0, 0, 0);gd.addMessage(myReference);

			//show dialog
			gd.showDialog();
		    if (gd.wasCanceled()) return null;
		    else {

		    	//domain definitions
		    	mRSO.cylZ1 = (int)Math.round(gd.getNextNumber());
		    	mRSO.cylZ2 = (int)Math.round(gd.getNextNumber());
		    	mRSO.cutAwayFromBottom = (int)Math.round(gd.getNextNumber());
		    	
		    	mRSO.maximalCylinder = gd.getNextBoolean();
		    	
		    	mRSO.cylX = (int)Math.round(gd.getNextNumber());
		    	mRSO.cylY = (int)Math.round(gd.getNextNumber());		    	
		    	mRSO.cylRadius = (int)Math.round(gd.getNextNumber());
		    }
		}

		if (mRSO.choiceOfRoi.equalsIgnoreCase("Cuboid")) {
		
			gd.addNumericField("Left X-coordinate of cuboid ", 0, 0, 6, "");
			gd.addNumericField("Right X-coordinate of cuboid (0 corresponds to the maximum x)", 0, 0, 6, "");
			gd.addNumericField("Topmost Y-coordinate of cuboid ", 0, 0, 6, "");
			gd.addNumericField("Bottommost Y-coordinate of cuboid (0 corresponds to the maximum y)", 0, 0, 6, "");
			gd.addNumericField("Z-coordinate of cuboid top ", 0, 0, 6, "");
			gd.addNumericField("Z-coordinate of cuboid bottom ", 0, 0, 6, "");

			//show dialog
			gd.showDialog();
		    if (gd.wasCanceled()) return null;
		    else {
	
		    	//domain definitions
		    	mRSO.cubeX1 = (int)Math.round(gd.getNextNumber());
		    	mRSO.cubeX2 = (int)Math.round(gd.getNextNumber());
		    	mRSO.cubeY1 = (int)Math.round(gd.getNextNumber());
		    	mRSO.cubeY2 = (int)Math.round(gd.getNextNumber());
		    	mRSO.cubeZ1 = (int)Math.round(gd.getNextNumber());
		    	mRSO.cubeZ2 = (int)Math.round(gd.getNextNumber());
		    }
		}

		gd.removeAll();
		
		return mRSO;
		
	}
	
	public ROISelectionOptions showHistogramExtractionMenu() {

		ROISelectionOptions mRSO = regionOfInterestSelection();
				
	    return mRSO;
	}	

	public PoreSpaceAnalyzerOptions showNetworkSpaceAnalyzerMenu() {
		
		PoreSpaceAnalyzerOptions mPSAO = new PoreSpaceAnalyzerOptions(); 
		mPSAO.mRSO = regionOfInterestSelection();  //select Region of interest
		
		GenericDialog gd2 = new GenericDialog("Please tell me what you want to get analysed!");
		
		gd2.addStringField("Image phase that is analyzed", "255", 5);
		gd2.addStringField("Name of the image phase that is analyzed", "Pores", 40);
		
		gd2.addMessage("");
		gd2.addNumericField("Number of copies of first and last slices to be added on top and bottom?", 10, 0);
		
		String myReference = "This plugin implements Ignacio Arganda-Carreras Skeletonize3D and AnalyzeSkleleton plugins.\n";
		myReference += "The skeletonization procedure was adapted to yield percolation theory related measures using an idea of Diego Soto.\n\n";
		myReference += "If you are using this plugin please cite the following reference: \n\n";
		gd2.setInsets(0, 0, 0);gd2.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd2.setInsets(0, 0, 0);gd2.addMessage(myReference);
		
		//show dialog
		gd2.showDialog();
	    if (gd2.wasCanceled()) return null;	    
	    else {

	    	mPSAO.calcCriticalPoreDiameter = false;
	    	
	    	//image phase 2 be analyzed
	    	mPSAO.imagePhase2BeAnalyzed = gd2.getNextString();
	    	mPSAO.nameOfAnalyzedPhase = gd2.getNextString();
	    	
	    	//number of slices to be added
	    	mPSAO.numberOfCopiedSlices2Add = (int)Math.round(gd2.getNextNumber());
	    	
	    }
		
		return mPSAO;
		
	}


	public ElsasOptions selectElsaOptions() {
		
		ElsasOptions mEO = new ElsasOptions();
		InputOutput jIO = new InputOutput();
		
		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";	
		
		InputOutput.MyFileCollection mFC = jIO.fileSelector("Please choose a file or folder with your earthorm data");
		
		mEO.myInFolder = mFC.myBaseFolder;
		mEO.myTiffs = mFC.myTiffs;
		
		mEO.myOutFolder = mFC.myBaseFolder + pathSep + "Unwarped";
		new File(mEO.myOutFolder).mkdir();
		
		mEO.myTargetImage = jIO.chooseAFile("Please choose a target image", mEO.myInFolder, "target.jpg");
		
		mEO.myTransformationFile = jIO.chooseAFile("Please choose the transformation file", mEO.myInFolder, "target.jpg");
	
		return mEO;
		
	}
	
	public LorenzosOptions selectLorenzosOptions() {
		
		LorenzosOptions mLO = new LorenzosOptions();
		InputOutput jIO = new InputOutput();
		
		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";	
					
		GenericDialog gd = new GenericDialog("Please tell me what time-step you want to extract");
		
		gd.addCheckbox("Run in batch mode", true);
		
		gd.addNumericField("Timestep in seconds", 5, 1);
		gd.addNumericField("Number of consecutive images to average ", 5, 0);
		
		//show dialog				
		gd.showDialog();			    
		if (gd.wasCanceled()) return null;	    
		else {
			
			mLO.runInBatchMode = gd.getNextBoolean();
			
			double myStep = gd.getNextNumber();
			mLO.stepsize = (int)Math.round(myStep / 0.133);
			if (mLO.stepsize < 1) mLO.stepsize = 1;
			
			mLO.average = (int)gd.getNextNumber();
			if (mLO.average < 1) mLO.average = 1;
			
		}
		
		if (mLO.runInBatchMode) {
						
			mLO.myInFolder = jIO.chooseAFolder("Please select a folder with your infiltration data");
			mLO.myAboveFolder = jIO.getTheFolderAbove(mLO.myInFolder, pathSep);
			
			mLO.myOutFolder = mLO.myAboveFolder + pathSep + "InfMovies_Step" + mLO.stepsize + "frames_Avg" + mLO.average;
			
			new File(mLO.myOutFolder).mkdir();
			
		}
		else {
			
			InputOutput.MyFileCollection mFC = jIO.fileSelectorAll("Please choose a file or folder with your infiltration data");
		
			mLO.myInFolder = mFC.myBaseFolder;			
			mLO.myAboveFolder = jIO.getTheFolderAbove(mFC.myBaseFolder, pathSep);
			
			mLO.myOutFolder = mLO.myAboveFolder + pathSep + "InfMovies_Step" + mLO.stepsize + "frames_Avg" + mLO.average;
		
			new File(mLO.myOutFolder).mkdir();
			
			
		}
		
		return mLO;
		
	}
	
	public PoreSpaceAnalyzerOptions showLoopDiameterDistMenu() {
		
		PoreSpaceAnalyzerOptions mLDD = new PoreSpaceAnalyzerOptions(); 
		mLDD.mRSO = regionOfInterestSelection();  //select Region of interest

		GenericDialog gd2 = new GenericDialog("Please tell me what you want to get analysed!");
		
		gd2.addStringField("Image phase that is analyzed", "255", 5);
		gd2.addStringField("Name of the image phase that is analyzed", "Pores", 40);
		gd2.addCheckbox("Do you want to fill holes in the analyzed image phase?", true);
		
		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd2.setInsets(50, 0, 0);gd2.addMessage(myReference);
		myReference = "Legland, D.; Arganda-Carreras, I. & Andrey, P. 2016. MorphoLibJ: integrated library and plugins for mathematical morphology with ImageJ.\n, ";
		myReference += "Bioinformatics (Oxford Univ Press) 32(22): 3532-3534, PMID 27412086, doi:10.1093/bioinformatics/btw413\n\n";
		gd2.setInsets(0, 0, 0);gd2.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd2.setInsets(0, 0, 0);gd2.addMessage(myReference);

		//show dialog
		gd2.showDialog();
	    if (gd2.wasCanceled()) return null;	    
	    else {

	    	mLDD.calcCriticalPoreDiameter = false;
	    
	    	//image phase 2 be analyzed
	    	mLDD.imagePhase2BeAnalyzed = gd2.getNextString();
	    	mLDD.nameOfAnalyzedPhase = gd2.getNextString();
	    	mLDD.removeHoles = gd2.getNextBoolean();
	    	
	    	//global measures
	    	mLDD.calcAnisotropy = false;
	    	mLDD.calcFractal = false;
	    	mLDD.calcThickness = false;
	    	mLDD.calcDistance = false;
	    	mLDD.calcCriticalPoreDiameter = false;
			
			//local measures
	    	// maybe re-include later..

			//plotting options
			mLDD.plotLabels = false;			
			mLDD.plotThickness = false;
			mLDD.plotPercolation =  false;
			mLDD.plotDistanceMap =  false;
			mLDD.plotPoresConnected2Top =  false;
			
			mLDD.performParticleAnalyses = true;
			
			//and the area of the horizontal cross section
			if (mLDD.mRSO.choiceOfRoi.equalsIgnoreCase("Everything!") | mLDD.mRSO.choiceOfRoi.equalsIgnoreCase("TopOfEveryThing") | mLDD.mRSO.choiceOfRoi.equalsIgnoreCase("BottomOfEveryThing")) {
				  
				mLDD.mRSO.areaOfInterest = gd2.getNextNumber();				
				
			}
			
			mLDD.mRSO.saveROI = false;

	    }

	    return mLDD;
		
	}
	
	public PoreSpaceAnalyzerOptions showPoreSpaceAnalyzerMenu() {

		PoreSpaceAnalyzerOptions mPSAO = new PoreSpaceAnalyzerOptions(); 
		mPSAO.mRSO = regionOfInterestSelection();  //select Region of interest

		GenericDialog gd2 = new GenericDialog("Please tell me what you want to get analysed!");
		
		gd2.addStringField("Image phase that is analyzed", "255", 5);
		gd2.addStringField("Name of the image phase that is analyzed", "Pores", 40);
		gd2.addCheckbox("Do you want to fill holes in the analyzed image phase?", true);

		String[] whatIntegralMeasuresShallIInvestigate = new String[]{"Anisotropy","Fractal Dimension","Thickness","Distance","Critical Diameter"};
		boolean[] myIntChoices = new boolean[]{true,true,true,true,true};
		gd2.setInsets(20, 100, 0);gd2.addMessage("\nWhich global morphological measures should I calculate in addition to the Minkowski functionals (volume, surface area, curvature, Euler number)?");
		gd2.setInsets(0, 100, 0);gd2.addCheckboxGroup(2, 3, whatIntegralMeasuresShallIInvestigate, myIntChoices);

		String[] whichImagesShallIPlot = new String[]{"Cluster Label","Thickness","Percolating Clusters","Distance Map","VolCon2Top"};
		boolean[] myPlotChoices = new boolean[]{false,false,false,false, true};
		gd2.setInsets(20, 200, 0);gd2.addMessage("\nWhich images shall I save?");
		gd2.setInsets(0, 200, 0);gd2.addCheckboxGroup(3, 3, whichImagesShallIPlot, myPlotChoices);
		
		if (mPSAO.mRSO.choiceOfRoi.equalsIgnoreCase("Everything!") | mPSAO.mRSO.choiceOfRoi.equalsIgnoreCase("TopOfEveryThing") | mPSAO.mRSO.choiceOfRoi.equalsIgnoreCase("BottomOfEveryThing")) {
		    	
	    	gd2.addNumericField("What is the cross-sectional area of the object of interest in pixels? ", 0, 0, 8, "");
	    	gd2.addMessage("(This is needed to calculate the porosity (or volume fraction of the investigated phase)");
	  
		}
		
		gd2.addMessage("");
		gd2.addCheckbox("Do you want to save the selected ROI?", false);
    	
		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd2.setInsets(50, 0, 0);gd2.addMessage(myReference);
		myReference = "Legland, D.; Arganda-Carreras, I. & Andrey, P. 2016. MorphoLibJ: integrated library and plugins for mathematical morphology with ImageJ.\n, ";
		myReference += "Bioinformatics (Oxford Univ Press) 32(22): 3532-3534, PMID 27412086, doi:10.1093/bioinformatics/btw413\n\n";
		gd2.setInsets(0, 0, 0);gd2.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd2.setInsets(0, 0, 0);gd2.addMessage(myReference);

		//show dialog
		gd2.showDialog();
	    if (gd2.wasCanceled()) return null;	    
	    else {

	    	mPSAO.calcCriticalPoreDiameter = false;
	    
	    	//image phase 2 be analyzed
	    	mPSAO.imagePhase2BeAnalyzed = gd2.getNextString();
	    	mPSAO.nameOfAnalyzedPhase = gd2.getNextString();
	    	mPSAO.removeHoles = gd2.getNextBoolean();
	    	
	    	//global measures
	    	mPSAO.calcAnisotropy = gd2.getNextBoolean();
	    	mPSAO.calcFractal = gd2.getNextBoolean();
	    	mPSAO.calcThickness = gd2.getNextBoolean();
	    	mPSAO.calcDistance = gd2.getNextBoolean();
	    	mPSAO.calcCriticalPoreDiameter = gd2.getNextBoolean();
			
			//local measures
	    	// maybe re-include later..

			//plotting options
			mPSAO.plotLabels = gd2.getNextBoolean();			
			mPSAO.plotThickness = gd2.getNextBoolean();
			mPSAO.plotPercolation = gd2.getNextBoolean();
			mPSAO.plotDistanceMap = gd2.getNextBoolean();
			mPSAO.plotPoresConnected2Top = gd2.getNextBoolean();
			
			mPSAO.performParticleAnalyses = true;
			
			//and the area of the horizontal cross section
			if (mPSAO.mRSO.choiceOfRoi.equalsIgnoreCase("Everything!") | mPSAO.mRSO.choiceOfRoi.equalsIgnoreCase("TopOfEveryThing") | mPSAO.mRSO.choiceOfRoi.equalsIgnoreCase("BottomOfEveryThing")) {
				  
				mPSAO.mRSO.areaOfInterest = gd2.getNextNumber();				
				
			}
			
			mPSAO.mRSO.saveROI = gd2.getNextBoolean();

	    }

	    return mPSAO;

	}
	
	public Extract2DHistogramOptions show2DHistogramExtractionMenu() {
		
		Extract2DHistogramOptions e2DH = new Extract2DHistogramOptions();		
		
		//construct objects
		GenericDialog gd = new GenericDialog("");
		
		gd.addMessage("The 2D histogram also analyses the gradient of your images!\n Please be aware that the plugin does not work for image filenames starting with 'Join' or 'POMR'!");
		
		gd.addCheckbox("I want to calculate the gradient images on the fly!", false);		
			
		//show dialog
		gd.showDialog();
		if (gd.wasCanceled()) return null;
		else {
			
			e2DH.calcGradientImage = gd.getNextBoolean();
			
		}				

		ROISelectionOptions mRSO = showHistogramExtractionMenu();
		e2DH.mRSO = mRSO;

		
		return e2DH;
		
	}
	

	
	public REVAnalyzerOptions showREVAnalyzerMenu() {

		///////////////////////////////////////
		// type of ROI and Method
		///////////////////////////////////////

		//construct objects
		GenericDialog gd = new GenericDialog("Please tell me what you want to get analysed!");

		REVAnalyzerOptions mRA = new REVAnalyzerOptions();

		String[] choiceOfRoi = new String[2];
		choiceOfRoi[0] = "Cuboid ROI";
		choiceOfRoi[1] = "Cylindrical ROI (defunct)";
		gd.addRadioButtonGroup("Please choose a region of interest!", choiceOfRoi, 1, 2, "Cuboid ROI");

		String[] choiceOfMethod = new String[3];
		choiceOfMethod[0] = "Sub-ROIs by division";
		choiceOfMethod[1] = "Sub-ROIs by shrinkage";
		//choiceOfMethod[2] = "Sub-ROIs by division and shrinkage";
		choiceOfMethod[2] = "Moving sub-ROIs";
		gd.addRadioButtonGroup("Please choose a region of interest!", choiceOfMethod, 3, 1, "Sub-ROIs by division");

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	//get ROI typ
	    	String myChoice = gd.getNextRadioButton();
	    	int myChoiceIndex = 0;
	    	for (int i = 0 ; i < choiceOfRoi.length ; i++) if (myChoice.equalsIgnoreCase(choiceOfRoi[i])) {
	    		myChoiceIndex = i;
	    		break;
	    	}
	    	switch (myChoiceIndex) {
	    		case 0: mRA.choiceOfRoi = "Cuboid"; break;
	    		case 1: IJ.error("This option is not yet implemented! Sorry!"); return null;  //mRA.choiceOfRoi = "Cylinder"; break;
	    	}

	    	//get Method
	    	myChoice = gd.getNextRadioButton();
	    	myChoiceIndex = 0;
	    	for (int i = 0 ; i < choiceOfMethod.length ; i++) if (myChoice.equalsIgnoreCase(choiceOfMethod[i])) {
	    		myChoiceIndex = i;
	    		break;
	    	}
	    	switch (myChoiceIndex) {
	    		case 0: mRA.choiceOfMethod = "Sub-ROIs by division"; break;
	    		case 1: mRA.choiceOfMethod = "Sub-ROIs by shrinkage"; break;
	    		//case 2: mRA.choiceOfMethod = "Sub-ROIs by division and shrinkage"; break;
	    		case 2: mRA.choiceOfMethod = "Moving sub-ROIs";; break;
	    	}
	    }

	    gd.removeAll();

		///////////////////////////////////////
		// ROI coordinates
		///////////////////////////////////////

	    if (mRA.choiceOfRoi.equalsIgnoreCase("Cuboid")) {

	    	GenericDialog cubeGD = new GenericDialog("Please input cube coordinates for the largest region of interest!");

	    	cubeGD.addNumericField("Left X-coordinate of cuboid ", 175, 0, 6, "");
	    	cubeGD.addNumericField("Right X-coordinate of cuboid (0 corresponds to the maximum x)", 575, 0, 6, "");
	    	cubeGD.addNumericField("Topmost Y-coordinate of cuboid ", 175, 0, 6, "");
	    	cubeGD.addNumericField("Bottommost Y-coordinate of cuboid (0 corresponds to the maximum y)", 575, 0, 6, "");
	    	cubeGD.addNumericField("Z-coordinate of cuboid top ", 300, 0, 6, "");
	    	cubeGD.addNumericField("Z-coordinate of cuboid bottom ", 700, 0, 6, "");

			//show dialog
	    	cubeGD.showDialog();
		    if (cubeGD.wasCanceled()) return null;
		    else {
		    	//domain definitions
		    	mRA.cubeX1 = (int)Math.round(cubeGD.getNextNumber());
		    	mRA.cubeX2 = (int)Math.round(cubeGD.getNextNumber());
		    	mRA.cubeY1 = (int)Math.round(cubeGD.getNextNumber());
		    	mRA.cubeY2 = (int)Math.round(cubeGD.getNextNumber());
		    	mRA.cubeZ1 = (int)Math.round(cubeGD.getNextNumber());
		    	mRA.cubeZ2 = (int)Math.round(cubeGD.getNextNumber());

		    	mRA.edgeX = mRA.cubeX2 - mRA.cubeX1;
		    	mRA.edgeY = mRA.cubeY2 - mRA.cubeY1;
		    	mRA.edgeZ = mRA.cubeZ2 - mRA.cubeZ1;
		    }
		    cubeGD.removeAll();
		}

	    if (mRA.choiceOfRoi.equalsIgnoreCase("Cylinder")) {
			gd.addNumericField("X-coordinate of cylinder center ", 0, 0, 6, "");
			gd.addNumericField("Y-coordinate of cylinder center ", 0, 0, 6, "");
			gd.addNumericField("Z-coordinate of cylinder top ", 0, 0, 6, "");
			gd.addNumericField("Z-coordinate of cylinder bottom ", 0, 0, 6, "");
			gd.addNumericField("Radius of cylinder ", 0, 0, 6, "");

			//show dialog
			gd.showDialog();
		    if (gd.wasCanceled()) return null;
		    else {
		    	//domain definitions
		    	mRA.cylX = (int)Math.round(gd.getNextNumber());
		    	mRA.cylY = (int)Math.round(gd.getNextNumber());
		    	mRA.cylZ1 = (int)Math.round(gd.getNextNumber());
		    	mRA.cylZ2 = (int)Math.round(gd.getNextNumber());
		    	mRA.cylRadius = (int)Math.round(gd.getNextNumber());
		    }

		    //check for evenness
	    	mRA.edgeX = 2 * mRA.cylRadius;
	    	mRA.edgeY = 2 * mRA.cylRadius;
	    	mRA.edgeZ = mRA.cubeZ2 - mRA.cubeZ1;
	    	if (Math.floorMod(mRA.edgeZ, 2) > 0) {
	    		mRA.edgeZ--;
	    		mRA.cubeZ2--;
	    	}
		}

		///////////////////////////////////////
		// the number of steps to be investigated for the ROI
		///////////////////////////////////////

		GenericDialog stepDG = new GenericDialog("Please tell me what you want to get analysed!");

		if (mRA.choiceOfMethod.equalsIgnoreCase("Sub-ROIs by division")) {
			
			stepDG.addNumericField("Please enter the starting division factor (1 stands for the whole domain) ", 1, 0);
			stepDG.addNumericField("Please enter the final division factor (must be smaller or equal to 8)", 8, 0);

			//show dialog
			stepDG.showDialog();
		    if (stepDG.wasCanceled()) return null;
		    else {
		    	mRA.startDivNumber = (int)Math.round(stepDG.getNextNumber());
		    	mRA.stopDivNumber = (int)Math.round(stepDG.getNextNumber());
		    }

		    //determine number of ROIS to analyze
		    String myMesg = "";
		    
		    //show message recapitulating the choices
		    myMesg += "You are about to run a ROI analysis using divisons.\n\n";
		
		    int roiNumber = 0;

		    for (int i=mRA.startDivNumber ; i <= mRA.stopDivNumber ; i++) {
		    	
		    	roiNumber += (int)Math.round(Math.pow(i,3));
		    	
			    myMesg += "This corresponds to "+ String.format("%3.0f",Math.pow(i, 3)) + " ROIs with edge lengths of " + mRA.edgeX/i +
			    		" vx, " + mRA.edgeY/i  + " vx and " + mRA.edgeZ/i +
			    		" vx in X, Y, and Z directions.\n\n";
		    }		    


		    myMesg += "You will analyze " + roiNumber + " different ROIs.";

		    if (IJ.showMessageWithCancel("Are you sure?", myMesg));
		    else return null;

		}

		if (mRA.choiceOfMethod.equalsIgnoreCase("Sub-ROIs by shrinkage")) {
			
			stepDG.addMessage("The shrinkage step-size s will be calculated as follows:\n" +
						"s = 0.66 L / (n + 1)\n" + 
						"where L is the edge length of the investigated domain and n is the number of shrinkage steps.\n" + 
						"The investigated scales will be L-s, L-2s, L-3s,... L-ns\n");
			stepDG.addNumericField("Please enter the number shrinkage steps n you want to investigate ", 10, 0);			

			//show dialog
			stepDG.showDialog();
		    if (stepDG.wasCanceled()) return null;
		    else {
		    	mRA.stepNumber = (int)Math.round(stepDG.getNextNumber());
		    }
		
		    //calculate shrinkage step-sizes
		    float sx = 0.666667f * (float)mRA.edgeX / ((float)mRA.stepNumber + 1);
		    float sy = 0.666667f * (float)mRA.edgeY / ((float)mRA.stepNumber + 1);
		    float sz = 0.666667f * (float)mRA.edgeZ / ((float)mRA.stepNumber + 1);
		    
		    //show message recapitulating the choices
		    String myMesg = "You are about to run a ROI analysis using " + mRA.stepNumber + " ROI shrinkages\n\n" +		    		

		    		"This corresponds to each 27 ROIs with edge lengths of ...\n";

		    for (int i = 1 ; i <= mRA.stepNumber ; i++) {
		    	myMesg += (int)Math.round((mRA.edgeX - i*sx)) + " vx, " + 
		    			  (int)Math.round((mRA.edgeY - i*sy)) + " vx and " + 
		    			  (int)Math.round((mRA.edgeZ - i*sz)) + " vx in X, Y, and Z directions.\n";
		    }

		    myMesg += "\nYou will analyze " + (27 * mRA.stepNumber) + " different ROIs.";

		    if (IJ.showMessageWithCancel("Are you sure?", myMesg));
		    else return null;

		}

		/*if (mRA.choiceOfMethod.equalsIgnoreCase("Sub-ROI by division and shrinkage")) {
			stepDG.addNumericField("Please enter the number division steps you want to investigate ", 1, 0);
			stepDG.addNumericField("Please enter the number shrinkage steps you want to investigate ", 15, 0);
			stepDG.addNumericField("Please enter the number size of the shrinkage steps ", 20, 0);

			//show dialog
			stepDG.showDialog();
		    if (stepDG.wasCanceled()) return null;
		    else {
		    	mRA.divNumber = (int)Math.round(stepDG.getNextNumber());
		    	mRA.stepNumber = (int)Math.round(stepDG.getNextNumber());
		    	mRA.stepLength = (int)Math.round(stepDG.getNextNumber());
		    }

		    //check if ROI is divisible by 2^stepnumber
		    String myMesg = "";

		    int maxdivisor = (int)Math.round(Math.pow(8, mRA.divNumber));
		    int restX = mRA.edgeX % maxdivisor; //Math.floorMod(mRA.edgeX,maxdivisor);
		    int restY = mRA.edgeY % maxdivisor; //Math.floorMod(mRA.edgeY,maxdivisor);
		    int restZ = mRA.edgeY % maxdivisor; //Math.floorMod(mRA.edgeZ,maxdivisor);

		    if (restX > 0) {
		    	myMesg += "X-edge is not divisible by " + maxdivisor + ". Changing rightmost X to " + (mRA.cubeX2 - restX) + ".\n\n";
		    	mRA.cubeX2 -= restX; mRA.edgeX -= restX;
		    }

		    if (restX > 0) {
		    	myMesg += "Y-edge is not divisible by " + maxdivisor + ". Changing bottommost Y to " + (mRA.cubeY2 - restY) + ".\n\n";
		    	mRA.cubeY2 -= restY; mRA.edgeY -= restY;
		    }

		    if (restX > 0) {
		    	myMesg += "Z-edge is not divisible by " + maxdivisor + ". Changing bottommost Z to " + (mRA.cubeZ2 - restZ) + ".\n\n";
		    	mRA.cubeZ2 -= restZ; mRA.edgeZ -= restZ;
		    }

		    //check whether step length is even
		    boolean hasBeenDecreased = false;
		    if (mRA.stepLength % 2 > 0) {
		    	mRA.stepLength++;
		    	hasBeenDecreased = true;
		    }

		    //check whether starting volume is enough to conduct all steps..
		    double[] edges = {mRA.edgeX, mRA.edgeY, mRA.edgeZ};
		    if (mRA.stepNumber * mRA.stepLength >= StatUtils.min(edges) / Math.pow(2, mRA.divNumber) - 1) {
		    	double smallestEdge = StatUtils.min(edges) / Math.pow(2, mRA.divNumber);
		    	mRA.stepLength = (int)Math.floor(smallestEdge / mRA.stepNumber) - 1;
		    	if (mRA.stepLength % 2 > 0) mRA.stepLength--;		//make sure that step-length is even.
		    	myMesg += "Step-number times step-length exceeds minimum edge length!\n";
		    	myMesg += "Step-length has therefore been reduced to " + mRA.stepLength + "!\n";
		    	myMesg += " \n";
		    }
		    else if (hasBeenDecreased) {
		    	myMesg += "Step-length must be even!\n";
		    	myMesg += "Step-length has therefore been increased to " + mRA.stepLength + "!\n";
		    	myMesg += " \n";
		    }

		    //check if step length is still 2 or larger
		    if (mRA.stepLength < 2) {
		    	IJ.error("This is silly! Please choose more reasonable parameters!\nBailing now...!");
		    	return null;
		    }

		    //recapitulate
		    myMesg += "You are about to run a ROI analysis using " + mRA.divNumber + " ROI divisons followed by " + mRA.stepNumber + " ROI shrinkages\n\n";
		    myMesg += "with an increment of " + mRA.stepLength + " vx.\n\n";
		    myMesg += 		"This corresponds to 1 ROI with edge lengths of " + mRA.edgeX + " vx, " +
		    		mRA.edgeY  + " vx and " + mRA.edgeZ + " vx in X, Y, and Z directions.\n\n";

		    int roiNumber = 1;

		    for (int i = 1 ; i <= mRA.divNumber ; i++) {
		    	myMesg += "This corresponds to "+ String.format("%3.0f",Math.pow(8, i)) + " ROIs with edge lengths of " + mRA.edgeX/(i*2) +
		    			" vx, " + mRA.edgeY/(i*2)  + " vx and " + mRA.edgeZ/(i*2) +
		    			" vx in X, Y, and Z directions.\n\n";
		    	roiNumber += Math.pow(8, i);
		    }

		    for (int i = 1 ; i <= mRA.stepNumber ; i++) {
		    	myMesg += (int)Math.round(Math.pow(8, mRA.divNumber)) + " ROIs with "+ (mRA.edgeX/(mRA.divNumber*2) - i*mRA.stepLength) +
		    			" vx, " + (mRA.edgeY/(mRA.divNumber*2) - i*mRA.stepLength)  + " vx and " +
		    			(mRA.edgeZ/(mRA.divNumber*2) - i*mRA.stepLength)+ " vx in X, Y, and Z directions.\n";
		    }

		    myMesg += "\nYou will analyze " + (roiNumber + mRA.stepNumber*(int)Math.round(Math.pow(8, mRA.divNumber))) + " different ROIs.";

		    if (IJ.showMessageWithCancel("Are you sure?", myMesg));
		    else return null;

		}*/
		
		if (mRA.choiceOfMethod.equalsIgnoreCase("Moving sub-ROIs")) {
	
			stepDG.addNumericField("Please enter the edge length of the moving ROI in X-direction", 250, 0);
			stepDG.addNumericField("Please enter the edge length of the moving ROI in Y-direction", 250, 0);
			stepDG.addNumericField("Please enter the edge length of the moving ROI in Z-direction", 250, 0);
			
			stepDG.addNumericField("Please enter the number of ROI moves in X-direction", 10, 0);
			stepDG.addNumericField("Please enter the number of ROI moves in Y-direction", 10, 0);
			stepDG.addNumericField("Please enter the number of ROI moves in Z-direction", 10, 0);

			//show dialog
			stepDG.showDialog();
		    if (stepDG.wasCanceled()) return null;
		    else {
		    	mRA.moveEdgeX = (int)Math.round(stepDG.getNextNumber());
		    	mRA.moveEdgeY = (int)Math.round(stepDG.getNextNumber());
		    	mRA.moveEdgeZ = (int)Math.round(stepDG.getNextNumber());
		    	
		    	mRA.movesX = (int)Math.round(stepDG.getNextNumber());
		    	mRA.movesY = (int)Math.round(stepDG.getNextNumber());
		    	mRA.movesZ = (int)Math.round(stepDG.getNextNumber());
		    }

		    //check whether step length is even
		    String myMesg = "";
		    boolean hasBeenDecreased = false;
		    if (mRA.stepLength % 2 > 0) {
		    	mRA.stepLength++;
		    	hasBeenDecreased = true;
		    }

		    //check whether starting volume is enough to conduct all steps..
		    double[] edges = {mRA.edgeX, mRA.edgeY, mRA.edgeZ};
		    if (mRA.stepNumber * mRA.stepLength >= StatUtils.min(edges)) {
		    	mRA.stepLength = (int)Math.round(Math.floor(StatUtils.min(edges) / mRA.stepNumber)) - 1;
		    	if (mRA.stepLength % 2 > 0) mRA.stepLength--;		//make sure that step-length is even.
		    	myMesg += "Step-number times step-length exceeds minimum edge length!\n";
		    	myMesg += "Step-size has therefore been reduced to " + mRA.stepLength + "!\n";
		    	myMesg += " \n";
		    }
		    else if (hasBeenDecreased) {
		    	myMesg += "Step-size must be even!\n";
		    	myMesg += "Step-size has therefore been increased to " + mRA.stepLength + "!\n";
		    	myMesg += " \n";
		    }

		    //show message recapitulating the choices
		    myMesg += "You are about to run a ROI analysis using " + (mRA.movesX * mRA.movesY * mRA.movesZ) + " moving ROIs.\n" +
		    
		    		"Each ROI will have an edge lengths of ...\n";
		    
		    myMesg += mRA.moveEdgeX + " vx, " + mRA.moveEdgeY  + " vx and " + mRA.moveEdgeZ + " vx in X, Y, and Z directions.\n";
		
		    if (IJ.showMessageWithCancel("Are you sure you want this?", myMesg));
		    else return null;

		}

		stepDG.removeAll();


		///////////////////////////////////////
		// the morphological measures needed
		///////////////////////////////////////

		GenericDialog gd2 = new GenericDialog("Please tell me what you want to get analysed!");

		String[] whatIntegralMeasuresShallIInvestigate = new String[]{"Macropore Volume","Average Thickness","Fractal Dimension","Anisotropy"};
		boolean[] myIntChoices = new boolean[]{true,true,true,true};
		gd2.setInsets(20, 200, 0);gd2.addMessage("\nWhich global morphological measures should I calculate?");
		gd2.setInsets(0, 200, 0);gd2.addCheckboxGroup(3, 3, whatIntegralMeasuresShallIInvestigate, myIntChoices);

		String[] whatShallIInvestigate = new String[]{"Volume", "Euler Number","Thickness",
				"Anisotropy","Percolation","Crit. Pore Diameter"};
		boolean[] myChoices = new boolean[]{true,true,true,true,true,true};
		gd2.setInsets(20, 200, 0);gd2.addMessage("\nWhich morphological measures should I calculate for each cluster?");
		gd2.setInsets(0, 200, 0);gd2.addCheckboxGroup(3, 3, whatShallIInvestigate, myChoices);

		String[] whichImagesShallIPlot = new String[]{"Binary ROI","Cluster Label","Volume","Thickness","Percolating Clusters","Anisotropy"};
		boolean[] myPlotChoices = new boolean[]{true, true,false,true,false,false};
		gd2.setInsets(20, 200, 0);gd2.addMessage("\nWhich images shall I save?");
		gd2.setInsets(0, 200, 0);gd2.addCheckboxGroup(3, 3, whichImagesShallIPlot, myPlotChoices);

		//show dialog
		gd2.showDialog();
	    if (gd2.wasCanceled()) return null;
	    else {

	    	//global measures
	    	mRA.globVolume = gd2.getNextBoolean();
			mRA.globThickness = gd2.getNextBoolean();
			mRA.calcFractal = gd2.getNextBoolean();
			mRA.globAnisotropy = gd2.getNextBoolean();

			//local measures
	    	mRA.calcVolume = gd2.getNextBoolean();
	    	mRA.calcEuler = gd2.getNextBoolean();
			mRA.calcThickness = gd2.getNextBoolean();
			mRA.calcAnisotropy = gd2.getNextBoolean();
			mRA.calcPercolation = gd2.getNextBoolean();
			mRA.calcCriticalPoreDiameter = gd2.getNextBoolean();

			//plotting options
			mRA.plotBinary = gd2.getNextBoolean();
			mRA.plotLabels = gd2.getNextBoolean();
			mRA.plotVolume = gd2.getNextBoolean();
			mRA.plotThickness = gd2.getNextBoolean();
			mRA.plotPercolation = gd2.getNextBoolean();

			mRA.performParticleAnalyses = false;
			if (mRA.calcVolume == true ||
					mRA.calcEuler == true ||
					mRA.calcThickness == true ||
					mRA.calcAnisotropy == true ||
					mRA.calcPercolation == true ||
					mRA.calcCriticalPoreDiameter == true ||

					//plotting options
					mRA.plotLabels == true ||
					mRA.plotVolume == true ||
					mRA.plotThickness == true ||
					mRA.plotPercolation == true) {

				mRA.performParticleAnalyses = true;

			}

	    }

	    //IJ.error("Plot thickness is " + mRA.plotThickness);

	    return mRA;

	}

	public class AnalyzeMeanGrayValueMenu {
		
		public int lowerThreshold;
		public int upperThreshold;
		
	}
	
	public AnalyzeMeanGrayValueMenu showAnalyzeMeanGrayValueMenu() {
		
		GenericDialog gd = new GenericDialog("Mean gray value calculator");

		AnalyzeMeanGrayValueMenu mMGV = new AnalyzeMeanGrayValueMenu();
		
		gd.addNumericField("Lower Threshold", 0, 0);
		gd.addNumericField("Upper Threshold", 65536, 0);
		
		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	
	    	mMGV.lowerThreshold = (int)gd.getNextNumber();
	    	mMGV.upperThreshold = (int)gd.getNextNumber();
	    	
	    }
	    
	    return mMGV;
	}	
	
	public class WRCCalculatorMenu {

		public ROISelectionOptions mRSO;
		
		public String choiceOfWRCType;
		
		public double voxelSizeInMicroMeter;		
		public double columnHeightInMM;
		public double wettingAngle;
		public double connectingSubresolutionDiameter;
		
		public double[] tensionStepsInMM;		
		public double[] tensionAtTop;		
		public double[] tensionAtCenter;
		public double[] tensionAtBottom;
		
		public boolean hasSurfaceFiles;
		public boolean saveAirWaterImages;	
		public boolean enforceIntegerTensions;
				
		public String imagePhase2BeAnalyzed = "255";
		
		public boolean calcAir;
		public boolean calcMatrixDistance;
		public boolean calcWater;
		
	}
	
	public class DrainageSimulatorOptions {

		public ROISelectionOptions mRSO;
				
		public double voxelSizeInMicroMeter;		
		public double columnHeightInMM;
		public double wettingAngle;
					
		public double[] pressureStepsInMM;		
		public double[] pressureAtTopInMM;		
		public double[] pressureAtCenterInMM;
		public double[] pressureAtBottomInMM;
		
		public boolean hasSurfaceFiles;
		public boolean saveAirWaterImages;	
		public boolean enforceIntegerPressures;
				
		public String imagePhase2BeAnalyzed = "255";
		
	}

	public WRCCalculatorMenu showWRCMenu() {
		
		//construct objects
		GenericDialog gd = new GenericDialog("Water retention curve calculator");
		WRCCalculatorMenu mWRC = new WRCCalculatorMenu();
		
		String myText = "This plugin simulates the water retention curve on a pore-network assuming the bottommost layer of the ROI as reference height.\n ";
		myText += "Isolated visible pore clusters are assumed to be connected by pores with diameters below the image resolution.\n";
		myText += "The connecting sub-resolution pore diameter is specified by a fraction of the voxel edgle length\n";
		myText += "(a value of 1 assumes a diameter of the voxel edge length; a value of 0 assumes the absence of any subresoluton pores)\n\n";
		gd.setInsets(0, 0, 0);gd.addMessage(myText);

		//define volume to be analyzed		
		mWRC.mRSO = regionOfInterestSelection();  //select Region of interest
		
		//choose draining or wetting curves..
		String[] choiceOfWRC = new String[1];
		choiceOfWRC[0] = "drainage";
		//choiceOfWRC[1] = "wetting from bottom";
		//choiceOfWRC[2] = "gravity flow";
		//choiceOfWRC[3] = "gravity flow with seepage face";
		gd.addRadioButtonGroup("Please choose a type of water retention curve", choiceOfWRC, 1, 1, "drainage");
		
		gd.addNumericField("Please enter the voxel edge length in micrometer", 43, 0, 6, "");
		gd.addNumericField("Please specify the diameter of the sub-resolution connecting pore diameter in voxels (value must be smaller than 2)", 1, 3, 6, "");
		gd.addNumericField("Please enter wetting angle (0 - 90 degree, 0 is perfect wetting)", 0, 1, 6, "");
		gd.addNumericField("Please enter approximate height of ROI in mm", 48, 2, 6, "");
		gd.addNumericField("How many tension steps do you want to calculate?", 5, 0, 3, "");
		
		//gd.addCheckbox("Do you want to take the soils top and bottom topographies into account?", false);

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
	    int tensionStepNumber = 0;
		gd.showDialog();
		int myChoiceIndex = 0;
	    if (gd.wasCanceled()) return null;
	    else {	    	
	    	//drainage or wetting
	    	String myChoice = gd.getNextRadioButton();
	    	
	    	for (int i = 0 ; i < choiceOfWRC.length ; i++) if (myChoice.equalsIgnoreCase(choiceOfWRC[i])) {
	    		myChoiceIndex = i;
	    		break;
	    	}	    	
	    	switch (myChoiceIndex) {
	    		case 0: mWRC.choiceOfWRCType = "drainage"; break;
	    		case 1: mWRC.choiceOfWRCType = "wetting from bottom"; break;
	    		case 2: mWRC.choiceOfWRCType = "gravity flow"; break;
	    		case 3: mWRC.choiceOfWRCType = "gravity flow with seepage face"; break;
	    	}

	    	mWRC.voxelSizeInMicroMeter = gd.getNextNumber();
	    	mWRC.connectingSubresolutionDiameter = gd.getNextNumber();
	    	mWRC.wettingAngle = gd.getNextNumber();
	    	mWRC.columnHeightInMM = gd.getNextNumber();
	    	tensionStepNumber = (int)Math.round(gd.getNextNumber());
	    	
	    	//mWRC.hasSurfaceFiles = gd.getNextBoolean();

	    }
	    
	    //calculate the highest reasonable tension to consider
		double capillaryConstant4MicroMeter = 1.48e7 * Math.cos(mWRC.wettingAngle / 180 * Math.PI);
		double maximalTensionInMicrometer = capillaryConstant4MicroMeter / mWRC.voxelSizeInMicroMeter;
		double maxTensionInMM = maximalTensionInMicrometer / 1000;
		double reasonableTension = maxTensionInMM;
	    if (mWRC.choiceOfWRCType.equalsIgnoreCase("drainage")) {
	    	reasonableTension = maxTensionInMM;
	    }
	    double[] tensionSteps = new double[tensionStepNumber]; 	
	    
	    //construct objects
	    GenericDialog gd2 = new GenericDialog("Water retention curve calculator - choose tension steps");
	    
		//choose draining or wetting curves..
		String[] wetEndTension = new String[3];
		wetEndTension[0] = "zero tension at top of ROI";
		wetEndTension[1] = "zero tension at vertical center of ROI";
		wetEndTension[2] = "zero tension at bottom of ROI";
		gd2.addRadioButtonGroup("Please choose a wet end reference point", wetEndTension, 3, 1, "zero tension at vertical center of ROI");
	    
	    gd2.addMessage("The tension at which water is drained from all visible pores connected to top surface is ca. " + String.format("%4.1f\t",reasonableTension)  + " mm.");
	    
		//choose draining or wetting curves..
		String[] dryEndTension = new String[3];
		dryEndTension[0] = "all visible pores at top surface of ROI are drained";
		dryEndTension[1] = "all visible connected pores at vertical center of ROI are drained";
		dryEndTension[2] = "all visible connected pores at bottom surface of ROI are drained";
		gd2.addRadioButtonGroup("Please choose a dry end reference point", dryEndTension, 3, 1, "all visible connected pores at bottom surface of ROI are drained");
	    	    
	    gd2.addCheckbox("Do you want to enforce integer values at the tension steps?", false);
	    
		myReference = "If you are using this plugin please cite the following references: \n\n";
		gd2.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd2.setInsets(0, 0, 0);gd.addMessage(myReference);
	    
		//show dialog
		boolean roundValues = false;
		gd2.showDialog();
	    if (gd2.wasCanceled()) return null;	    
	    else {
	    	
	    	//wet end
	    	String myChoice = gd2.getNextRadioButton();
	    	
	    	for (int i = 0 ; i < wetEndTension.length ; i++) {
	    		boolean itsTrueNow = myChoice.equalsIgnoreCase(wetEndTension[i]);
	    		if (itsTrueNow) {	    		
	    			myChoiceIndex = i;
	    			break;
	    		}
	    	}	    	
	    	switch (myChoiceIndex) {
	    		case 0: tensionSteps[0] = -mWRC.columnHeightInMM; break;
	    		case 1: tensionSteps[0] = -mWRC.columnHeightInMM / 2; break;
	    		case 2: tensionSteps[0] = 0; break;	    		
	    	}

	    	//dry end
	    	myChoice = gd2.getNextRadioButton();
	    	
	    	for (int i = 0 ; i < dryEndTension.length ; i++) if (myChoice.equalsIgnoreCase(dryEndTension[i])) {
	    		myChoiceIndex = i;
	    		break;
	    	}	    	
	    	switch (myChoiceIndex) {
	    		case 0: tensionSteps[tensionStepNumber - 1] = reasonableTension - mWRC.columnHeightInMM; break;
	    		case 1: tensionSteps[tensionStepNumber - 1] = reasonableTension - mWRC.columnHeightInMM / 2; break;
	    		case 2: tensionSteps[tensionStepNumber - 1] = reasonableTension; break;	    		
	    	}
	    	
	    	roundValues = gd2.getNextBoolean();    	
	    	mWRC.enforceIntegerTensions = roundValues;
	    }
	    
	    //assign tension steps
	    double stepSize = (reasonableTension - tensionSteps[0]) / (tensionStepNumber - 1);
	    if (tensionStepNumber > 2) for (int i = 1 ; i < tensionStepNumber - 1; i++) {
	    	tensionSteps[i] = i * stepSize + tensionSteps[0];
	    }
	    
	    //round values in case option was checked
	    if (roundValues) {
	    	tensionSteps[0] = Math.round(tensionSteps[0]); 
	    	tensionSteps[tensionStepNumber - 1] =  Math.round(tensionSteps[tensionStepNumber - 1]);	    	
		    if (tensionStepNumber > 2) for (int i = 1 ; i < tensionStepNumber - 1; i++) {
		    	tensionSteps[i] = Math.round(i * stepSize + tensionSteps[0]);
		    }	    	
	    }
	    
	    //construct objects
	    GenericDialog gd3 = new GenericDialog("Water retention curve calculator - modify tension steps (reference depth is the bottom surface of te ROI)");
	    
	    if (roundValues) {
		    for (int i = 0 ; i < tensionStepNumber ; i++) gd3.addNumericField("Modify tension step #" + i  + " at the bottom of the ROI?", tensionSteps[i], 0, 6, " mm");
	    }
		else {	    
			for (int i = 0 ; i < tensionStepNumber ; i++) gd3.addNumericField("Modify tension step #" + i  + " at the bottom of the ROI?", tensionSteps[i], 2, 6, " mm");			
		}
	    
	    gd3.addCheckbox("Do you want to save the images of air and water?", false);
	    
		myReference = "If you are using this plugin please cite the following references: \n\n";
		gd3.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd3.setInsets(0, 0, 0);gd3.addMessage(myReference);
	    
	    gd3.showDialog();
	    if (gd3.wasCanceled()) return null;	    
	    else {
	    	for (int i = 0 ; i < tensionStepNumber ; i++) {
	    		tensionSteps[i] = gd3.getNextNumber();
	    	}
	    	
	    	mWRC.saveAirWaterImages = gd3.getNextBoolean();
	    }
	  
	    mWRC.tensionStepsInMM = tensionSteps;	    
	    
//	    //construct objects
//	    GenericDialog gd4 = new GenericDialog("Water retention curve calculator - reduce memory requirements?");
//	    
//	    myText = "Calculating all implemented properties of air, water and matrix phases is very memory intensive.\n ";
//		myText += "In case you run out of memory, you may want to deselect some of the following calculations.\n";
//		gd4.setInsets(0, 0, 0);gd4.addMessage(myText);
//	    
//	    gd4.addCheckbox("Yes, please calculate morphological properties of the air phase", true);
//	    gd4.addCheckbox("Yes, please calculate morphological properties of the water phase", true);
//	    gd4.addCheckbox("Yes, please calculate distances to the next air-filled pore", true);
//	    
//	    gd4.showDialog();
//	    if (gd4.wasCanceled()) return null;	    
//	    else {	    	
//	    	mWRC.calcAir = gd4.getNextBoolean();
//	    	mWRC.calcWater = gd4.getNextBoolean();
//	    	mWRC.calcMatrixDistance = gd4.getNextBoolean();
//	    }
//	    
//		myReference = "If you are using this plugin please cite the following references: \n\n";
//		gd4.setInsets(40, 0, 0);gd.addMessage(myReference);
//		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
//		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
//		gd4.setInsets(0, 0, 0);gd4.addMessage(myReference);
	    
		return mWRC;
	}
	
	public DrainageSimulatorOptions showDrainageSimulatorMenu() {
		
		//construct objects
		GenericDialog gd = new GenericDialog("Drainage simulation");
		DrainageSimulatorOptions mDS = new DrainageSimulatorOptions();
		
		String myText = "This plugin simulates water and air phases in the pore network  under drainage. Hydraulic equilibrium is assumed\n ";
		gd.setInsets(0, 0, 0);gd.addMessage(myText);

		
		//define volume to be analyzed		
		mDS.mRSO = regionOfInterestSelection();  //select Region of interest
		
		//choose draining or wetting curves..
		String[] choiceOfWRC = new String[1];
		choiceOfWRC[0] = "drainage";
		//choiceOfWRC[1] = "wetting from bottom";
		//choiceOfWRC[2] = "gravity flow";
		//choiceOfWRC[3] = "gravity flow with seepage face";
		gd.addRadioButtonGroup("Please choose a type of water retention curve", choiceOfWRC, 1, 1, "drainage");
		
		gd.addNumericField("Please enter the voxel edge length in micrometer", 43, 0, 6, "");
		gd.addNumericField("Please enter wetting angle (0 - 90 degree, 0 is perfect wetting)", 0, 1, 6, "");
		gd.addNumericField("Please enter approximate height of ROI in mm", 48, 2, 6, "");
		
		gd.setInsets(0, 0, 0);gd.addMessage("");
		
		gd.addNumericField("How many pressure steps do you want to calculate?", 10, 0, 3, "");
		gd.addNumericField("Initial pressure at upper boundary in mm", 0, 2, 3, "");
		gd.addNumericField("Constant pressure at lower boundary in mm", 0, 2, 3, "");
		
		gd.addCheckbox("Do you want to enforce integer values at the tension steps?", false);
		
		//gd.addCheckbox("Do you want to take the soils top and bottom topographies into account?", false);
		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
	    int pressureStepNumber = 0;
		boolean roundValues = false;
		double initialPressureAtTop = 0;
		double pressureAtBottom = 0;
		
		gd.showDialog();
		int myChoiceIndex = 0;
	    if (gd.wasCanceled()) return null;
	    else {	    	
	    		    	
	    	mDS.voxelSizeInMicroMeter = gd.getNextNumber();
	    	mDS.wettingAngle = gd.getNextNumber();
	    	mDS.columnHeightInMM = gd.getNextNumber();
	    	
	    	pressureStepNumber = (int)Math.round(gd.getNextNumber());
	    	initialPressureAtTop = gd.getNextNumber();
	    	pressureAtBottom = gd.getNextNumber();
	    	
	    	roundValues = gd.getNextBoolean();    	
	    	mDS.enforceIntegerPressures = roundValues;
	    	
	    	//mDS.hasSurfaceFiles = gd.getNextBoolean();

	    }
	    
	    //calculate the highest reasonable pressure to consider
	    double[] pressureSteps = new double[pressureStepNumber]; 	
	    
	    //assign pressure steps
	    pressureSteps[0] = initialPressureAtTop + mDS.columnHeightInMM;
	    pressureSteps[pressureStepNumber - 1] = pressureAtBottom;
	    double stepSize = (pressureSteps[0] - pressureAtBottom) / (pressureStepNumber - 1);
	    if (pressureStepNumber > 2) for (int i = 1 ; i < pressureStepNumber - 1; i++) {
	    	pressureSteps[i] = pressureSteps[0] - i * stepSize;
	    }
	    
	    //round values in case option was checked
	    if (roundValues) {
	    	pressureSteps[0] = Math.round(pressureSteps[0]); 
	    	pressureSteps[pressureStepNumber - 1] =  Math.round(pressureSteps[pressureStepNumber - 1]);	    	
		    if (pressureStepNumber > 2) for (int i = 1 ; i < pressureStepNumber - 1; i++) {
		    	pressureSteps[i] = Math.round(pressureSteps[0] - i * stepSize);
		    }	    	
	    }
	    
	    //construct objects
	    GenericDialog gd3 = new GenericDialog("Drainage simulator - modify pressure steps (reference depth is the bottom surface of te ROI)");
	    
	    if (roundValues) {
		    for (int i = 0 ; i < pressureStepNumber ; i++) gd3.addNumericField("Modify pressure step #" + i  + " at the bottom of the ROI?", pressureSteps[i], 0, 6, " mm");
	    }
		else {	    
			for (int i = 0 ; i < pressureStepNumber ; i++) gd3.addNumericField("Modify pressure step #" + i  + " at the bottom of the ROI?", pressureSteps[i], 2, 6, " mm");			
		}	    
	    gd3.addCheckbox("Do you want to save the images of air and water?", false);
	    
		myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);
	    
	    gd3.showDialog();
	    if (gd3.wasCanceled()) return null;	    
	    else {
	    	for (int i = 0 ; i < pressureStepNumber ; i++) {
	    		pressureSteps[i] = gd3.getNextNumber();
	    	}
	    	
	    	mDS.saveAirWaterImages = gd3.getNextBoolean();
	    }
	  
	    mDS.pressureStepsInMM = pressureSteps;	    
	    
		return mDS;
	}

	public class ColumnFinderMenuReturn implements Serializable {

		public boolean isAlreadyNormalized;
		public boolean isSteel;
		public boolean isPVC;
		public boolean isAlu;
		public boolean try2FindColumnTopAndBottom;
		public boolean hasBevel;
		public boolean grayValueOfWallIsKnown;
		public boolean putColumnStraight;
		public boolean doAResliceFirst;
		
		//cut off of top and bottom
		public double topCutOff;
		public double botCutOff;

		//debug settings
		public boolean debug;
		public boolean showRadialProfiles;
		public boolean showFit;

		//???
		public int segmentLength;

		//a priori fixed column properties..
		public double outerDiameter;
		public int minColGV;
		public int maxColGV;
		public int wallThickness;
		public int topOfColumn;
		public int bottomOfColumn;
		public int fixedWallGrayValue;
		public int stdFixedWallGrayValue;		

		//fitting parameters needed for finding outer wall
		public double airWallContrast;
		public double airDividedByWall;

		//fitting parameter for finding inner wall
		public double wallSoilStdContrastThreshold;

		//parameters needed to decide whether the column is at this depth
		public double CVWallBrightnessThresh;
		public double maxCVOfWallThickness;
		public double stdThreshold;
		public double ratioBetweenInnerAndOuterRadius;
		public double percentToleratedDifference;

		//if column has already been normalized
		public double absoluteWallgrayValue;

		//optimaization parameters
		public boolean applyUmask;
		public int medianFilter2D;
		public int medianFilter;
		public double r2Thresh;
		public int maxFittingAttempts;		//attempts allowed to fiddle with "airWallContrast" and "wallSoilStdContrastThreshold" to find the edge
		public int maxNumberOfOutliers4PerimeterFits;

	}
	
	public ColumnFinderMenuReturn deepCopyCFS(ColumnFinderMenuReturn virginCFS, ColumnFinderMenuReturn jCFS) {
		
		virginCFS.isAlreadyNormalized = jCFS.isAlreadyNormalized;
		virginCFS.isSteel = jCFS.isSteel;
		virginCFS.isPVC = jCFS.isPVC;
		virginCFS.isAlu = jCFS.isAlu;
		virginCFS.try2FindColumnTopAndBottom = jCFS.try2FindColumnTopAndBottom;
		virginCFS.hasBevel = jCFS.hasBevel;
		virginCFS.grayValueOfWallIsKnown = jCFS.grayValueOfWallIsKnown;
		virginCFS.putColumnStraight = jCFS.putColumnStraight;
		virginCFS.doAResliceFirst = jCFS.doAResliceFirst;
		
		//cut off of top and bottom
		virginCFS.topCutOff = jCFS.topCutOff;
		virginCFS.botCutOff = jCFS.botCutOff;

		//debug settings
		virginCFS.debug = jCFS.debug;
		virginCFS.showRadialProfiles = jCFS.showRadialProfiles;
		virginCFS.showFit = jCFS.showFit;

		//???
		virginCFS.segmentLength = jCFS.segmentLength;

		//a priori fixed column properties..
		virginCFS.outerDiameter = jCFS.outerDiameter;
		virginCFS.minColGV = jCFS.minColGV;
		virginCFS.maxColGV = jCFS.maxColGV;
		virginCFS.wallThickness = jCFS.wallThickness;
		virginCFS.topOfColumn = jCFS.topOfColumn;
		virginCFS.bottomOfColumn = jCFS.bottomOfColumn;
		virginCFS.fixedWallGrayValue = jCFS.fixedWallGrayValue;
		virginCFS.stdFixedWallGrayValue = jCFS.stdFixedWallGrayValue;		

		//fitting parameters needed for finding outer wall
		virginCFS.airWallContrast = jCFS.airWallContrast;
		virginCFS.airDividedByWall = jCFS.airDividedByWall;

		//fitting parameter for finding inner wall
		virginCFS.wallSoilStdContrastThreshold = jCFS.wallSoilStdContrastThreshold;

		//parameters needed to decide whether the column is at this depth
		virginCFS.CVWallBrightnessThresh = jCFS.CVWallBrightnessThresh;
		virginCFS.maxCVOfWallThickness = jCFS.maxCVOfWallThickness;
		virginCFS.stdThreshold = jCFS.stdThreshold;
		virginCFS.ratioBetweenInnerAndOuterRadius = jCFS.ratioBetweenInnerAndOuterRadius;
		virginCFS.percentToleratedDifference = jCFS.percentToleratedDifference;

		//if column has already been normalized
		virginCFS.absoluteWallgrayValue = jCFS.absoluteWallgrayValue;

		//optimaization parameters
		virginCFS.applyUmask = jCFS.applyUmask;
		virginCFS.medianFilter2D = jCFS.medianFilter2D;
		virginCFS.medianFilter = jCFS.medianFilter;
		virginCFS.r2Thresh = jCFS.r2Thresh;
		virginCFS.maxFittingAttempts = jCFS.maxFittingAttempts;		//attempts allowed to fiddle with "airWallContrast" and "wallSoilStdContrastThreshold" to find the edge
		virginCFS.maxNumberOfOutliers4PerimeterFits = jCFS.maxNumberOfOutliers4PerimeterFits;
		
		return virginCFS;
		
	}

	public ColumnFinderMenuReturn showColumnStraightenerMenu() {

		//construct objects
		GenericDialog gd = new GenericDialog("PVC Column Finder Dialog");

		ColumnFinderMenuReturn mCFS = new ColumnFinderMenuReturn();

		//for now (Feb 2017), do not show anything..
		
		//add choices to menu
	/*	String[] columnMaterials = {"steel", "PVC", "aluminium"};
		gd.addRadioButtonGroup("What kind of material is your column made of?", columnMaterials, 3, 1, columnMaterials[2]);

		gd.addNumericField("Minimum contrast between air and wall relative to image contrast (0..1)", 0.25, 3, 6, "");
		gd.addNumericField("Footprint of median filter for radial wall position search", 5, 0, 3, "");

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2017. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	String myChoice = gd.getNextRadioButton();
	    	int myNumericChoice = 0;
	    	for (int i = 0 ; i < columnMaterials.length ; i++) if (myChoice.equalsIgnoreCase(columnMaterials[i])) myNumericChoice = i;
	    	switch (myNumericChoice) {
	    		case 0:	{
	    			mCFS.isSteel = true;
	    			mCFS.isPVC = false;
	    			mCFS.isAlu = false;
	    			break;
	    		}
	    		case 1: {
	    			mCFS.isSteel = false;
	    			mCFS.isPVC = true;
	    			mCFS.isAlu = false;
	    			break;
	    		}
	    		case 2: {
	    			mCFS.isSteel = false;
	    			mCFS.isPVC = false;
	    			mCFS.isAlu = true;
	    			break;
	    		}
	    	}

	    	mCFS.airWallContrast = (float)gd.getNextNumber();
	    	mCFS.medianFilter = (int)Math.round(gd.getNextNumber());
	    }*/

		//for now, assign the properties manually 
		mCFS.isSteel = false;
		mCFS.isPVC = false;
		mCFS.isAlu = true;
		
		mCFS.airWallContrast = 0.25f;
    	mCFS.medianFilter = 5;
		
	    return mCFS;

	}

	public ColumnFinderMenuReturn showColumnFinderDialog() {

		//construct objects
		GenericDialog gd = new GenericDialog("Column Finder Dialog");

		ColumnFinderMenuReturn mCFS = new ColumnFinderMenuReturn();
		
		//add choices to menu
		String[] viewPoints = {"perpendicular to Z-axis","parallel to Z-axis"};		
		gd.addRadioButtonGroup("How are the 2-D slices in your 3-D image aligned with the Z-axis?", viewPoints, 2, 1, viewPoints[0]);
		
		gd.addCheckbox("Do you want the program to move your column upright into the center of the canvas?", false);
		gd.addCheckbox("This column has already calibrated gray-values", false);
		
		gd.addMessage("");
		
		gd.addCheckbox("These columns are bevelled!", true);
		gd.addCheckbox("Do you want the program to find top and bottom of your column?", true);
		//gd.addNumericField("Height at which the top of the 3-D canvas should be cut off above the found top (% of found column height)", 10, 0);
		//gd.addNumericField("Depth at which the bottom of the 3-D canvas should be cut off below the found bottom (% of found column height)", 10, 0);
		gd.addMessage("");

		String[] columnMaterials = {"aluminium or PVC column and a BAD contrast between soil matrix and wall",
				"aluminium or PVC column and a GOOD contrast between soil matrix and wall", "steel"};
		gd.addRadioButtonGroup("What kind of material is your column made of?", columnMaterials, 3, 1, columnMaterials[0]);

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	String myView = gd.getNextRadioButton();
	    	if (myView.equalsIgnoreCase("parallel to Z-axis")) mCFS.doAResliceFirst = true;
	    	else mCFS.doAResliceFirst = false;
	    	
	    	mCFS.putColumnStraight = gd.getNextBoolean();
	    	mCFS.isAlreadyNormalized = gd.getNextBoolean();
	    	
			mCFS.hasBevel = gd.getNextBoolean();
	    	mCFS.try2FindColumnTopAndBottom = gd.getNextBoolean();
	    	//mCFS.topCutOff = gd.getNextNumber();
	    	//mCFS.botCutOff = gd.getNextNumber();
	    	
	    	String myChoice = gd.getNextRadioButton();
	    	int myNumericChoice = 0;
	    	for (int i = 0 ; i < columnMaterials.length ; i++) if (myChoice.equalsIgnoreCase(columnMaterials[i])) myNumericChoice = i;
	    	switch (myNumericChoice) {	    		
	    		case 0: {
	    			mCFS.isSteel = false;
	    			mCFS.isPVC = true;
	    			mCFS.isAlu = false;
	    			break;
	    		}
	    		case 1: {
	    			mCFS.isSteel = false;
	    			mCFS.isPVC = false;
	    			mCFS.isAlu = true;
	    			break;
	    		}
	    		case 2:	{  
	    			mCFS.isSteel = true;
	    			mCFS.isPVC = false;
	    			mCFS.isAlu = false;
	    			break;
	    		}
	    	}
	    }

	    gd.removeAll();

		/*gd.addMessage("");
		gd.addMessage("Advanced Options");
		gd.addMessage("");*/
	    
	    GenericDialog gd2 = new GenericDialog("Column Finder Dialog");
	    String[] noiseLevel = {"insanely high (R2 = 0.95)","incredibly high (R2 = 0.99)","very high (R2 = 0.995)", 
	    		"high (R2 = 0.999)", "medium (R2 = 0.9995)", "low (R2 = 0.9999)", "very low (R2 = 0.99995)"};
	    			   
	    //gd2.addNumericField("optional: what is the outer diameter of the column in pixels?", 0, 0, 6, "");	
	    gd2.addNumericField("optional: what is the minimum possible gray value of the column material?", 0, 0, 6, "");	
	    gd2.addNumericField("optional: what is the maximum possible gray value of the column material?", 0, 0, 6, "");	
	    
		if (mCFS.isPVC == true | mCFS.isSteel == true) gd2.addNumericField("How thick is the column wall in pixels? ", 20, 0, 6, "");
		if (!mCFS.try2FindColumnTopAndBottom) {
			gd2.addNumericField("In which layer is the top of the column? (if you do not want to search for it)", 1, 0, 4, "");
			gd2.addNumericField("In which layer is the bottom of the column? (enter '1' for the bottommost layer)", 1, 0, 4, "");
		}

		//gd2.addMessage("Criteria for finding the outer wall perimeter");
		if (mCFS.isAlreadyNormalized) {
			gd2.addNumericField("Mean gray-value of column wall ", 1000, 0, 6, "");
			gd2.addNumericField("Standard deviation gray-value of column wall ", 150, 0, 6, "");
		}
		
		String insertSpace = "\n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(insertSpace);		
		
	/*	else {
			gd2.addNumericField("Minimum contrast between air and wall relative to image contrast (0..1)", 0.3, 3, 6, "");
			gd2.addNumericField("Absolute minimum contrast between air and wall", 100, 0, 6, "");
		}*/

		/*if (mCFS.isAlu == true) {
			gd2.addMessage("Criteria for finding the inner wall perimeter");
			gd2.addNumericField("StDev. contrast threshold between wall and soil (50 - 2000) ", 500, 0, 6, "");
		}*/
		
		//criteria to check whether there is a wall
		if (!mCFS.isSteel) {
			//gd2.addMessage("Criteria to check whether there is a wall");		
			//gd2.addNumericField("Ratio between inner and outer radius must be larger than ...", 0.75, 3, 6, "");
			//gd2.addNumericField("Maximum allowed coefficient of variation of wall-thickness", 0.14, 4, 6, "");
			//gd2.addNumericField("Maximal CV of wall gray values", 0.14, 4, 6, "");
			//gd2.addNumericField("Neighboring column outline coordinates may be different by max", 0.5, 4, 6, "%");
	
			//gd2.addMessage("Filtering options");
			//gd2.addNumericField("Footprint of 2-D median filter image cross-section", 0, 0, 3, "");
			//gd2.addCheckbox("Use an unsharp mask", true);
	
			//gd2.addMessage("Criteria for goodnes of fit");
			//gd2.addNumericField("Maximal number of edge finding trials", 10, 0, 6, "");
			//gd2.addNumericField("Maximal number of outliers during ellipse-fit", 10, 0, 6, "");
			gd2.addNumericField("Gray value outside the wall divided by gray value of wall? ", 0.2, 1, 3, "");						
			gd.setInsets(40, 0, 0);gd.addMessage("(increase this value if the outer perimeter is found inside the column;\n "
					+ "decrease this value if the perimeter is found outside the column)");
		    gd.setInsets(40, 0, 0);gd.addMessage(insertSpace);					
			
			gd2.addChoice("Select an R2 criterion according to the noise level in image and your image size", noiseLevel, noiseLevel[5]);		
			gd.setInsets(40, 0, 0);gd.addMessage("The R2 corresponds to the coefcicient of determination of the ellipse fit to the found positions of the column wall. "
					+ "The larger the R2 criteria is put, the more precise column location will be found;/n"
					+ "However, the more noise inherent to the image, the less accurate the colum location can be determined.\n"
					+ "Moreover, the less pixels there are in an image in the XY-plane, the lower the R2 of the ellipse fit will be by trend.\n"
					+ "In such cases, the column may have been detected but discarded because of a too large R2 criterion.\n"
					+ "Check whether this is the case in the pop-up in the column finder visualization adn reduce the R2 criterion if needed.");
		}
			
		gd2.addCheckbox("Visualize column wall finds", false);

		myReference = "If you are using this plugin please cite the following references: \n\n";
		gd2.setInsets(40, 0, 0);gd2.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd2.setInsets(0, 0, 0);gd2.addMessage(myReference);
		

		//show dialog
		gd2.showDialog();

		if (gd2.wasCanceled()) return null;
		else {
			
			//mCFS.outerDiameter = gd2.getNextNumber();
			mCFS.minColGV = (int)Math.round(gd2.getNextNumber());
			mCFS.maxColGV = (int)Math.round(gd2.getNextNumber());
			
			if (mCFS.isPVC | mCFS.isSteel) mCFS.wallThickness = (int)Math.round(gd2.getNextNumber());
			if (!mCFS.try2FindColumnTopAndBottom) {
				mCFS.topOfColumn = (int)Math.round(gd2.getNextNumber());
				mCFS.bottomOfColumn = (int)Math.round(gd2.getNextNumber());				
			}

			//get parameters to get outer edge
			if (mCFS.isAlreadyNormalized) {
				mCFS.fixedWallGrayValue = (int)gd2.getNextNumber();
				mCFS.stdFixedWallGrayValue = (int)gd2.getNextNumber();
			}
			else {
				mCFS.airWallContrast = 0.3;//(float)gd2.getNextNumber();
				mCFS.stdThreshold = 1;//(float)gd2.getNextNumber();
			}

			// get parameters to find inner edge
			if (mCFS.isAlu) mCFS.wallSoilStdContrastThreshold = 500;//gd2.getNextNumber();

			//get parameters to check whether the wall is there
			if (!mCFS.isSteel) {
				mCFS.ratioBetweenInnerAndOuterRadius = 0.75;//gd2.getNextNumber();			
				mCFS.maxCVOfWallThickness = 0.14;//gd2.getNextNumber();
				mCFS.CVWallBrightnessThresh = 0.14;//gd2.getNextNumber();
				mCFS.percentToleratedDifference = 1;//gd2.getNextNumber();
	
				//get optimization parameters;				
				mCFS.medianFilter2D = 0;//(int)Math.round(gd2.getNextNumber());
				mCFS.applyUmask = true;//gd2.getNextBoolean();
				mCFS.maxFittingAttempts = 10;//(int)Math.round(gd2.getNextNumber());
				mCFS.maxNumberOfOutliers4PerimeterFits = 10;//(int)Math.round(gd2.getNextNumber());
				
				mCFS.airDividedByWall = gd2.getNextNumber();
				
				int myNoiseLevel = gd2.getNextChoiceIndex();
				if (myNoiseLevel == 0) mCFS.r2Thresh = 0.95;
				if (myNoiseLevel == 1) mCFS.r2Thresh = 0.99;
				if (myNoiseLevel == 2) mCFS.r2Thresh = 0.995;
				if (myNoiseLevel == 3) mCFS.r2Thresh = 0.999;
				if (myNoiseLevel == 4) mCFS.r2Thresh = 0.9995;
				if (myNoiseLevel == 5) mCFS.r2Thresh = 0.9999;
				if (myNoiseLevel == 6) mCFS.r2Thresh = 0.99995;
			}

			mCFS.debug = gd2.getNextBoolean();
			
			return mCFS;
		}
	}

	public ColumnFinderMenuReturn showTuneThreshold4ColumnDetectionMenu(ApproximateColumnIllumination aCI) {

		//construct objects
		GenericDialog gd = new GenericDialog("PVC Column Finder Dialog");

		ColumnFinderMenuReturn mCFS = new ColumnFinderMenuReturn();

		gd.addMessage("The likely median gray value outside the column is " + aCI.outside);
		gd.addMessage("The likely median gray value of the column wall is " + aCI.wall);
		gd.addMessage("The likely median gray value inside the column is " + aCI.inside);
		gd.addMessage("The global contrast for this image is " + aCI.globalContrast);
		gd.addMessage("");

		double suggestedAirWallThreshold = 0.7 * ((double)(aCI.wall - aCI.outside) / (double)aCI.globalContrast);
		double suggestedWallSoilThreshold = 0.7 * ((double)(aCI.wall - aCI.inside) / (double)aCI.globalContrast);

		gd.addNumericField("Minimum contrast between air and wall relative to average gray level (0..1)",suggestedAirWallThreshold, 3, 6, "");
		gd.addNumericField("Minimum contrast between wall and soil relative to average gray level (0..1)", suggestedWallSoilThreshold, 3, 6, "");

		//show dialog
		gd.showDialog();

		if (gd.wasCanceled()) return null;
		else {

			mCFS.airWallContrast = (float)gd.getNextNumber();

			return mCFS;
		}

	}

	public RandomClusterGenerator showRandomClusterGeneratorMenu() {

		GenericDialog gd = new GenericDialog("random cluster generator menu");

		RandomClusterGenerator mRCG = new RandomClusterGenerator();

		gd.addNumericField("What is the width (x) of your domain (vx)? ", 200, 0, 4, "");
		gd.addNumericField("What is the length (y) of your domain (vx)? ", 200, 0, 4, "");
		gd.addNumericField("What is the height (z) of your domain (vx)? ", 200, 0, 4, "");

		String[] choiceOfPoroClasses = new String[3];
		choiceOfPoroClasses[0] = "Range between two bounds";
		choiceOfPoroClasses[1] = "Gaussian around pc";
		choiceOfPoroClasses[2] = "List of porosities in a file";
		gd.addChoice("Please choose porosity classes!", choiceOfPoroClasses, "Gaussian around pc");

		String[] choiceOfShape = new String[2];
		choiceOfShape[0] = "Cubic";
		choiceOfShape[1] = "Cylindric";
		gd.addChoice("Please choose a shape of your domain!", choiceOfShape, "cylinder");

		gd.addNumericField("In case you chose 'Range between two bounds', please give me a lower bound for the random porosities! ", 0.0f, 3, 5, "");
		gd.addNumericField("In case you chose 'Range between two bounds', please give me a upper bound for the random porosities! ", 0.2f, 3, 5, "");

		gd.addNumericField("In case you chose 'Gaussian around pc', please give me a relative standard deviation! ", 0.1f, 3, 5, "");

		gd.addNumericField("How many copies of this random field do you want to create? ", 1, 0, 4, "");

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	mRCG.domainX = (int)Math.round(gd.getNextNumber());
	    	mRCG.domainY = (int)Math.round(gd.getNextNumber());
	    	mRCG.domainZ = (int)Math.round(gd.getNextNumber());
	    	double[] porosityRange = new double[]{gd.getNextNumber(), gd.getNextNumber()};
	    	double standardDeviation = gd.getNextNumber();
	    	int poroClassIndex = gd.getNextChoiceIndex();
	    	switch (poroClassIndex) {
				case 0: {
					mRCG.mode = "range";
					mRCG.porosityBounds = porosityRange;
					break;
					}
				case 1: {
					mRCG.mode = "Gaussian";
					mRCG.standardDeviation = standardDeviation;
					break;
				}
				case 2: {
					mRCG.mode = "predefinedList";
					break;
				}
	    	}
	    	int typeChoiceIndex = gd.getNextChoiceIndex();
		    switch (typeChoiceIndex) {
		    	case 0: {
					mRCG.shape = "Cubic";
					break;
					}
				case 1: {
					mRCG.shape = "Cylindric";
					break;
				}
	    	}
	    	mRCG.numOfCopies = (int)Math.round(gd.getNextNumber());
	    }

		return mRCG;

	}

	public SurfaceFinderReturn showSurfaceFinderMenu() {

		//construct objects
		GenericDialog gd = new GenericDialog("Soil surface finder");

		//init SurfaceFinderReturn
		SurfaceFinderReturn mSFR = new SurfaceFinderReturn();
	
		//assign pores that are open to the surface to the soil bulk volume?
		gd.addMessage("Do you want filter out small soil crumbs loosely on the soil surface?");

		gd.addNumericField("Diameter of soil crumbs to filter out (vx) ", 8, 0, 3, "");

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show the dialog and harvest the information
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	mSFR.neglectCrumbsOfLessThan = (int)Math.round(gd.getNextNumber());

	    	return mSFR;
	    }

	}

	public ThresholderMenuReturn showThresholdingDialog() {

		//construct objects
		GenericDialog gd = new GenericDialog("Thresholding dialog");

		ThresholderMenuReturn mTMR = new ThresholderMenuReturn();

		mTMR.useConstantThreshold = false; //set constant threshold option to a default false

		int filterChoices = 0;

		String[] myChoices = new String[12];
		myChoices[0] = "Otsu";
		myChoices[1] = "Renyi Entropy";
		myChoices[2] = "Maximum Entropy";
		myChoices[3] = "Mean";
		myChoices[4] = "Moments";
		myChoices[5] = "Huang";
		myChoices[6] = "IsoData";
		myChoices[7] = "DefaultIsoData";
		myChoices[8] = "IJ_IsoData";
		myChoices[9] = "Minimum";
		myChoices[10] = "Triangle";
		myChoices[11] = "UserDefinedThreshold";

		String[] mySecondaryChoices = new String[myChoices.length];
		mySecondaryChoices[0]="no two-step segmentation";
		for (int i = 1 ; i < myChoices.length ; i++) mySecondaryChoices[i] = myChoices[i - 1];

		//add threshold method choices
		gd.addMessage("Segmentation Options");

		gd.addCheckbox("Do you want to use the column outlines defined in 'inner circle' as a ROI?", false);

		//gd.addMessage("Search-string for only segmenting some of the files in this folder. Leave blank for segmenting them all.");
		//gd.addStringField("", "", 25);

		gd.addCheckbox("Do you want to cut away all gray values brighter than the wall?", false);

		gd.addMessage("Please choose a primary thresholding algorithm!");
		gd.addRadioButtonGroup("", myChoices, 5, 2, myChoices[11]);

		gd.addMessage("In case you want to perform a two-step segmentation approach, choose a secondary thresholding algorithm!");
		gd.addRadioButtonGroup("", mySecondaryChoices, 5, 2, mySecondaryChoices[0]);

		//saving options
		gd.addMessage("");
		gd.addMessage("Saving options");

		gd.addCheckbox("Save the segmented, 3-D binary images", false);

		gd.addCheckbox("Save segmentation results as overlays in some sample slices", true);

		gd.addCheckbox("Export TIFF stacks for GeoDict", false);

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	//get segmentation options
	    	mTMR.useInnerCircle = gd.getNextBoolean();

	    /*	mTMR.filterTag = gd.getNextString();
	    	if(mTMR.filterTag.isEmpty()) mTMR.filterImages = false;
	    	else mTMR.filterImages = true;*/

	    	mTMR.setMaxgray2Wallgray = gd.getNextBoolean();

	    	String choice = gd.getNextRadioButton();
	    	for (int i = 0 ; i < myChoices.length ; i++) if (myChoices[i].equalsIgnoreCase(choice)) {
	    		filterChoices = i;
	    		break;
	    	}
	    	switch (filterChoices) {
				case 0 : {mTMR.myPrimaryMethod = AutoThresholder.Method.Otsu; break;}
				case 1 : {mTMR.myPrimaryMethod = AutoThresholder.Method.RenyiEntropy; break;}
				case 2 : {mTMR.myPrimaryMethod = AutoThresholder.Method.MaxEntropy; break;}
				case 3 : {mTMR.myPrimaryMethod = AutoThresholder.Method.Mean; break;}
				case 4 : {mTMR.myPrimaryMethod = AutoThresholder.Method.Moments; break;}
				case 5 : {mTMR.myPrimaryMethod = AutoThresholder.Method.Huang; break;}
				case 6 : {mTMR.myPrimaryMethod = AutoThresholder.Method.IsoData; break;}
				case 7 : {mTMR.myPrimaryMethod = AutoThresholder.Method.Default; break;}
				case 8 : {mTMR.myPrimaryMethod = AutoThresholder.Method.IJ_IsoData; break;}
				case 9 : {mTMR.myPrimaryMethod = AutoThresholder.Method.Minimum; break;}
				case 10 : {mTMR.myPrimaryMethod = AutoThresholder.Method.Triangle; break;}
				case 11 : {mTMR.myPrimaryMethod = null;
							mTMR.useConstantThreshold = true;
							break;}
	    	}

	    	choice = gd.getNextRadioButton();
	    	for (int i = 0 ; i < myChoices.length ; i++) if (mySecondaryChoices[i].equalsIgnoreCase(choice)) {
	    		filterChoices = i;
	    		break;
	    	}
	    	switch (filterChoices) {
	    		case 0 : {mTMR.mySecondaryMethod = null; break;}
	    		case 1 : {mTMR.mySecondaryMethod = AutoThresholder.Method.Otsu; break;}
				case 2 : {mTMR.mySecondaryMethod = AutoThresholder.Method.RenyiEntropy; break;}
				case 3 : {mTMR.mySecondaryMethod = AutoThresholder.Method.MaxEntropy; break;}
				case 4 : {mTMR.mySecondaryMethod = AutoThresholder.Method.Mean; break;}
				case 5 : {mTMR.mySecondaryMethod = AutoThresholder.Method.Moments; break;}
				case 6 : {mTMR.mySecondaryMethod = AutoThresholder.Method.Huang; break;}
				case 7 : {mTMR.mySecondaryMethod = AutoThresholder.Method.IsoData; break;}
				case 8 : {mTMR.mySecondaryMethod = AutoThresholder.Method.Default; break;}
				case 9 : {mTMR.mySecondaryMethod = AutoThresholder.Method.IJ_IsoData; break;}
				case 10 : {mTMR.myPrimaryMethod = AutoThresholder.Method.Minimum; break;}
				case 11 : {mTMR.myPrimaryMethod = AutoThresholder.Method.Triangle; break;}
	    	}

	    	//get saving options
	    	mTMR.save3DImage = gd.getNextBoolean();
	    	mTMR.save4Evaluation = gd.getNextBoolean();
	    	mTMR.save4GeoDict = gd.getNextBoolean();

	      	return mTMR;
	    }
	}
	
	public POMThresholderMenuReturn showPOMThresholdingDialog() {

		//construct objects
		GenericDialog gd = new GenericDialog("Thresholding dialog");

		POMThresholderMenuReturn mTMR = new POMThresholderMenuReturn();
		
		//add threshold method choices
		gd.addMessage("Organic Material Segmentation Options");
	
		gd.addMessage("");
		gd.addMessage("Make sure that you have run 'CalibrateGrayValues' when you use this plugin!!!");

		gd.addNumericField("Enter the lower threshold for segmenting the image ", 10000, 0, 5, "");

		gd.addNumericField("Enter the upper threshold for segmenting the image ", 16500, 0, 5, "");
		
		gd.addNumericField("Enter black top hat image threshold ", 900, 0, 5, "");
		
		gd.addNumericField("Enter gradient image threshold ", 9000, 0, 5, "");
		
		//gd.addNumericField("Enter gray-value range of search window ", 3500, 0, 5, "");
		
		//gd.addNumericField("How many erosions steps do you want to run on each range ", 2, 0, 1, "");

		//gd.addNumericField("Enter search window overlap in percent ", 80, 0, 2, "");
		
		gd.addNumericField("Filter away all POM or root clusters with less than ", 50, 0, 5, " voxels (0 will skip this step)");

		//saving options
		gd.addMessage("");
		gd.addMessage("Saving options");

		gd.addCheckbox("Save the segmented, 3-D binary images", true);

		gd.addCheckbox("Save segmentation results as overlays in some sample slices", false);		
		
		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. 2018. SoilJ: An ImageJ plugin for the semiautomatic processing of three-dimensional X-ray images of soils.\n ";
		myReference += "Vadose Zone Journal, doi:10.2136/vzj2017.03.0062.";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	mTMR.minThreshold = (int)Math.round(gd.getNextNumber());
	    	mTMR.maxThreshold = (int)Math.round(gd.getNextNumber());

	    	//mTMR.windowSize = gd.getNextNumber();
	    	//mTMR.openingSteps = (int)Math.round(gd.getNextNumber());
	    	//mTMR.overlap = gd.getNextNumber();

	    	mTMR.blackTopHatThreshold = (int)Math.round(gd.getNextNumber());
	    	mTMR.gradientThreshold = (int)Math.round(gd.getNextNumber());
	    	
	    	mTMR.minPOMVoxelNumber = (int)Math.round(gd.getNextNumber());
	    	
	    	//get saving options
	    	mTMR.save3DImage = gd.getNextBoolean();
	    	mTMR.save4Evaluation = gd.getNextBoolean();	

	      	return mTMR;
	    }
	}
	
	public BioPoreExtractionOptions showBioPoreExtractionMenu() {
		
		BioPoreExtractionOptions mBEO = new BioPoreExtractionOptions();
		
		GenericDialog gd = new GenericDialog("Give me some biopore extraction parameters, please.");

		//prepare radio box with the choices
		gd.addCheckbox("Do you want to skip calculating on the original image resoluttion (check if your image sizes are large)", false);
		gd.addNumericField("Please enter the minimum vesselness a biopore must have", 0.6, 1);
		gd.addNumericField("Please enter the minimum length in voxels a biopore must have", 5, 1);
		//gd.addNumericField("Please enter maximum footprint of the Gaussian Blur", 100, 0);
		
		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	//get ROI typ
	    	mBEO.doNotProcessOriginalResolution = gd.getNextBoolean();
	    	mBEO.thresholdVesselness = gd.getNextNumber();
	    	mBEO.smallesAllowedElongation = gd.getNextNumber();
	    	//mBEO.maximumBlurring = (int)gd.getNextNumber();
	    }
		
		return mBEO;
	}
	
	public RootExtractionOptions showRootsExtractionMenu() {
		
		RootExtractionOptions mREO = new RootExtractionOptions();
		
		GenericDialog gd = new GenericDialog("Give me some root extraction parameters, please.");
		
		//gd.addMessage("");
		//gd.addMessage("Make sure that you have run 'CalibrateGrayValues' when you use this plugin!!!");

		//gd.addNumericField("Enter the lower threshold for segmenting the image ", 10000, 0, 5, "");

		//gd.addNumericField("Enter the upper threshold for segmenting the image ", 16500, 0, 5, "");

		gd.addMessage("");
		
		//prepare radio box with the choices
		gd.addCheckbox("Do you want to skip calculating on the original image resoluttion (check if your image sizes are large)", false);
		gd.addNumericField("Please enter the minimum vesselness a root must have", 0, 1, 3, " (a value of 0 switches of this filter criterion)");
		gd.addNumericField("Please enter the minimum length in voxels a root must have", 30, 1);
		//gd.addNumericField("Please enter maximum footprint of the Gaussian Blur", 100, 0);
		
		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	
	    	//mREO.minThreshold = (int)Math.round(gd.getNextNumber());
	    	//mREO.maxThreshold = (int)Math.round(gd.getNextNumber());
	    	
	    	mREO.doNotProcessOriginalResolution = gd.getNextBoolean();
	    	mREO.thresholdVesselness = gd.getNextNumber();
	    	mREO.smallesAllowedElongation = gd.getNextNumber();
	    	
	    }
		
		return mREO;
	}
	
	public CrackExtractionOptions showCrackExtractionMenu() {
		
		CrackExtractionOptions mCEO = new CrackExtractionOptions();
		
		GenericDialog gd = new GenericDialog("Give me some biopore extraction parameters, please.");

		//prepare radio box with the choices
		gd.addCheckbox("Do you want to skip calculating on the original image resoluttion (check if your image sizes are large)", false);
		gd.addNumericField("Please enter the minimum platyness a crack must have", 0.3, 1);
		gd.addNumericField("Please enter the minimum extension in voxels a crack must have", 20, 1);
		//gd.addNumericField("Please enter maximum footprint of the Gaussian Blur", 100, 0);
		
		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	//get ROI typ
	    	mCEO.doNotProcessOriginalResolution = gd.getNextBoolean();
	    	mCEO.thresholdPlanarity = gd.getNextNumber();
	    	mCEO.smallesAllowedElongation = gd.getNextNumber();
	    	//mCEO.maximumBlurring = (int)gd.getNextNumber();
	    }
		
		return mCEO;
	}
	
	public GravelExtractionOptions showGravelExtractionMenu() {
		
		GravelExtractionOptions mGEO = new GravelExtractionOptions();
		
		GenericDialog gd = new GenericDialog("Give me some gravel and stones extraction parameters, please.");
		
		gd.addNumericField("Please let me know the voxel size of the image in millimeter", 0.04, 3);
		gd.addNumericField("Please give me a threshold for the gravel fraction", 17000, 0);
		gd.addNumericField("Please tell me the footprint of the median filter I should apply", 5, 0);
		gd.addCheckbox("Shall I filter out the sand fraction (diameter < 2 mm)?", true);
		
		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	
	    	mGEO.voxelsize = gd.getNextNumber();
	    	mGEO.threshold = (int)gd.getNextNumber();	   
	    	mGEO.medianFilter = gd.getNextNumber();
	    	mGEO.filterOutSandFraction = gd.getNextBoolean();
	    	
	    }
		
		return mGEO;
	}	
	
	public GrayValues2Rescale showGrayValueRescaleMenu() {
		
		GrayValues2Rescale myGrays = new GrayValues2Rescale();
		
		GenericDialog gd = new GenericDialog("Let me know how you want to rescale your gray values, please.");
		
		gd.addNumericField("Old lower reference value", 5000, 0);
		gd.addNumericField("Old higher reference value", 20000, 0);
		
		gd.addMessage("");
		gd.addMessage("These values will be scaled to ...");
		gd.addMessage("");
		
		gd.addNumericField("New lower reference value", 7000, 0);
		gd.addNumericField("New higher reference value", 20000, 0);
		
		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	
	    	myGrays.oldLow = (int)gd.getNextNumber();
	    	myGrays.oldHigh = (int)gd.getNextNumber();
	    	myGrays.newLow = (int)gd.getNextNumber();
	    	myGrays.newHigh = (int)gd.getNextNumber();
	    	
	    }
		
		return myGrays;
		
	}
	
	public SelectFiles showFileSelectionMenu(File file) {

		//construct objects
		InputOutput jIO = new InputOutput();
		GenericDialog gd = new GenericDialog("Do you want to process the selected file or all files in the selected folder?");

		SelectFiles mSF = new SelectFiles();
		
		//list files in directory
		File myPath = file.getParentFile();
		String[] myTiffs = jIO.listTiffsInFolder(myPath);
		
		if (myTiffs.length > 1) {
			
			String altChoice = ""; 
			if (myTiffs.length == 2) altChoice = myTiffs[0] + "     " + myTiffs[1] + "   ";
			if (myTiffs.length == 3) altChoice = myTiffs[0] + "     " + myTiffs[1] + "     " + myTiffs[2]+ "   ";
			if (myTiffs.length == 4) altChoice = myTiffs[0] + "     " + myTiffs[1] + "     " + myTiffs[2] + "     " + myTiffs[3]+ "   ";
			if (myTiffs.length > 4) altChoice = myTiffs[0] + "     " + myTiffs[1] + "     " + myTiffs[2] + "     " + myTiffs[3] + " ...   ";

			//prepare radio box with the choices
			String[] choiceOfRoi = new String[2];
			choiceOfRoi[0] = file.getName();
			choiceOfRoi[1] = altChoice;
			gd.addRadioButtonGroup("Please enter your choice", choiceOfRoi, 2, 1, file.getName());			

			//show dialog
			gd.showDialog();
		    if (gd.wasCanceled()) return null;
		    else {
		    	//get ROI typ
		    	String myChoice = gd.getNextRadioButton();
		    	int myChoiceIndex = 0;
		    	for (int i = 0 ; i < choiceOfRoi.length ; i++) if (myChoice.equalsIgnoreCase(choiceOfRoi[i])) {
		    		myChoiceIndex = i;
		    		break;
		    	}
		    	switch (myChoiceIndex) {
		    		case 0: mSF.selectChosen = true; break;
		    		case 1: mSF.selectChosen = false; 
		    	}
	
		    } 
		}
		
		else mSF.selectChosen = true;
		
		return mSF;
	}
	
	public SelectFiles showHistogramSelectionMenu(File file, int bitDepth) {

		//construct objects
		InputOutput jIO = new InputOutput();
		GenericDialog gd = new GenericDialog("Do you want to process the selected file or all files in the selected folder?");

		SelectFiles mSF = new SelectFiles();
		
		//list files in directory
		File myPath = file.getParentFile();
		String[] myHists = jIO.listHistsInFolder(myPath);
		if (bitDepth == 8) myHists = jIO.listHistsInFolder8(myPath);
		
		if (myHists.length > 1) {
			
			String altChoice = ""; 
			if (myHists.length == 2) altChoice = myHists[0] + "     " + myHists[1] + "   ";
			if (myHists.length == 3) altChoice = myHists[0] + "     " + myHists[1] + "     " + myHists[2]+ "   ";
			if (myHists.length == 4) altChoice = myHists[0] + "     " + myHists[1] + "     " + myHists[2] + "     " + myHists[3]+ "   ";
			if (myHists.length > 4) altChoice = myHists[0] + "     " + myHists[1] + "     " + myHists[2] + "     " + myHists[3] + " ...   ";

			//prepare radio box with the choices
			String[] choiceOfRoi = new String[2];
			choiceOfRoi[0] = file.getName();
			choiceOfRoi[1] = altChoice;
			gd.addRadioButtonGroup("Please enter your choice", choiceOfRoi, 2, 1, file.getName());			

			//show dialog
			gd.showDialog();
		    if (gd.wasCanceled()) return null;
		    else {
		    	//get ROI typ
		    	String myChoice = gd.getNextRadioButton();
		    	int myChoiceIndex = 0;
		    	for (int i = 0 ; i < choiceOfRoi.length ; i++) if (myChoice.equalsIgnoreCase(choiceOfRoi[i])) {
		    		myChoiceIndex = i;
		    		break;
		    	}
		    	switch (myChoiceIndex) {
		    		case 0: mSF.selectChosen = true; break;
		    		case 1: mSF.selectChosen = false; 
		    	}
	
		    } 
		}
		
		else mSF.selectChosen = true;
		
		return mSF;
	}

}



