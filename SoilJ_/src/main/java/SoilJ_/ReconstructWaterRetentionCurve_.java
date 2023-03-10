package SoilJ_;

import java.io.File;

import SoilJ.tools.DisplayThings;
import SoilJ.tools.ImageManipulator;
import SoilJ.tools.ImageManipulator.StackCalculator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.InputOutput.MyFileCollection;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.MorphologyAnalyzer;
import SoilJ.tools.MorphologyAnalyzer.WRCPhaseProps;
import SoilJ.tools.ObjectDetector;
import SoilJ.tools.RoiHandler;
import ij.IJ;

/**
 *SoilJ is a collection of ImageJ plugins for the semi-automatized processing of 3-D X-ray images of soil columns
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

import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

/** 
 * ExtractAirFilledMacropores is a SoilJ plugin for extracting the air-filled pore space 
 * from a pore-diameter image (thickness image) for soil columns with a surface exposed to te atmosphere.
 * The plugin is based in the capillary equation. The soil is assumed to be in equilibrium with a specific 
 * matrix potential defined at the lower boundary of the column.
 * The plugin does not take hysteretic effects into account.
 * The binary images of the air filled pore space are saved in a newly created folder.
 * 
 * @author John Koestel
 *
 */

public class ReconstructWaterRetentionCurve_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {
		
		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";

		//set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = ReconstructWaterRetentionCurve_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ObjectDetector jOD = new ObjectDetector();
		ImageManipulator jIM = new ImageManipulator();
		MenuWaiter menu = new MenuWaiter();
		StackCalculator sC= jIM.new StackCalculator();
		RoiHandler roi = new RoiHandler();
		MorphologyAnalyzer morph = new MorphologyAnalyzer();
		DisplayThings disp = new DisplayThings();
		
		// init variables
		int i;
		MenuWaiter.WRCCalculatorMenu mWRC = menu.new WRCCalculatorMenu();
		mWRC = menu.showWRCMenu();
		if (mWRC == null) return;
		
		//construct image related objects
		ImagePlus airwaterTiff;  //output
		
		//read base folder and number of 3D images
		MyFileCollection mFC = jIO.fileSelector("Please choose a file or folder with your image data");				
						
		//create output folder
		String myOutFolder = "WaterRetentionData";	
		String myOutPath = mFC.myBaseFolder + pathSep + myOutFolder;
		new File(myOutPath).mkdir();		
		mFC.myOutFolder = myOutPath;
		
	    //monitor memory
	    Runtime rt = Runtime.getRuntime();
	    long currInUse = (rt.totalMemory() - rt.freeMemory()) / 1000000000;  //in GB
						
		//loop over 3D images
		for (i = 0 ; i < mFC.myTiffs.length ; i++) {  

			//try to free up some memory
			System.gc();System.gc();
			
			//assemble image for analyses
			mFC.fileName = mFC.myTiffs[i];			
			mFC = jIO.addCurrentFileInfo(mFC);
			
			//create new output folder
			String nowOutPath = myOutPath + mFC.pathSep + mFC.colName;
			new File(nowOutPath).mkdir();	
			mFC.myOutFolder = nowOutPath;
			
			//load file			
			int[] startStopSlices = jIO.findStartAndStopSlices(mFC, mWRC.mRSO);
			int[] colSlices = new int[startStopSlices[1] - startStopSlices[0] + 1];
			for (int j = 0 ; j < colSlices.length ; j++) colSlices[j] = startStopSlices[0] + j;
		
			//cut roi		
			mFC.startSlice = startStopSlices[0];
			mFC.stopSlice = startStopSlices[1];
			RoiHandler.ColumnRoi colRoi = roi.prepareDesiredRoi(mFC, jIO.openTiff3DSomeSlices(mFC, colSlices), mWRC.mRSO, "null", "null");
			System.gc();System.gc();
			currInUse = (rt.totalMemory() - rt.freeMemory()) / 1000000000; 	
			
			//calculate tension at top, center and bottom of column
			mWRC.columnHeightInMM = mWRC.voxelSizeInMicroMeter * colRoi.nowTiff.getNSlices() / 1000;
			
			//colRoi.nowTiff.updateAndDraw(); colRoi.nowTiff.show();
					
			//make a binarized version of the current Tiff..
			ImageStack binStack = new ImageStack(colRoi.nowTiff.getWidth(), colRoi.nowTiff.getHeight(), colRoi.nowTiff.getNSlices());			
			for (int z = 0 ; z < colRoi.nowTiff.getNSlices() ; z++) {
				colRoi.nowTiff.setPosition(z + 1);
				ImageProcessor thIP = colRoi.nowTiff.getProcessor();
				ImageProcessor cpIP = thIP.duplicate();
				cpIP.max(2);
				ImageProcessor binIP = cpIP.convertToByte(false);
				binIP.multiply(200);	//so that all non-zero values are 255  (400 is clipped at 255 for 8-bit images)
				binStack.setProcessor(binIP, z+1);		
			}
			ImagePlus binTiff = new ImagePlus("bin", binStack);
			System.gc();System.gc();
			currInUse = (rt.totalMemory() - rt.freeMemory()) / 1000000000; 
			
			//create drainage map...
			IJ.showStatus("Creating drainage map..");
			ImagePlus drainageTiff = jOD.createDrainageMap(mWRC.tensionStepsInMM, colRoi.nowTiff, null, mWRC, mFC);
			colRoi.nowTiff.unlock();colRoi.nowTiff.flush();
			System.gc();System.gc();
			currInUse = (rt.totalMemory() - rt.freeMemory()) / 1000000000; 
			
			//binTiff.updateAndDraw();binTiff.show();
			
			//init output structure
			WRCPhaseProps myWRCProperties = morph.new WRCPhaseProps();
			
			//init output arrays
			double[] tensionAtBottomInMM = new double[mWRC.tensionStepsInMM.length];
			double[] tensionAtCenterInMM = new double[mWRC.tensionStepsInMM.length];
			double[] tensionAtTopInMM = new double[mWRC.tensionStepsInMM.length];

			double[] areaInMM2 = new double[mWRC.tensionStepsInMM.length];
			double[] heightInMM = new double[mWRC.tensionStepsInMM.length];
			double[] bulkVolumeInMM3 = new double[mWRC.tensionStepsInMM.length];
			
			double[] theta_w = new double[mWRC.tensionStepsInMM.length];			
			double[] sigma_w = new double[mWRC.tensionStepsInMM.length];
			double[] chi_w = new double[mWRC.tensionStepsInMM.length];			
			double[] gamma_w = new double[mWRC.tensionStepsInMM.length];
			double[] fractalDim_w = new double[mWRC.tensionStepsInMM.length];
			double[] percolates_w = new double[mWRC.tensionStepsInMM.length]; 				
			double[] thetaLC_w = new double[mWRC.tensionStepsInMM.length];
			double[] thetaPerc_w = new double[mWRC.tensionStepsInMM.length];
			
			double[] theta_a = new double[mWRC.tensionStepsInMM.length];
			double[] sigma_a = new double[mWRC.tensionStepsInMM.length];
			double[] chi_a = new double[mWRC.tensionStepsInMM.length];
			double[] gamma_a = new double[mWRC.tensionStepsInMM.length];
			double[] fractalDim_a = new double[mWRC.tensionStepsInMM.length];
			double[] percolates_a = new double[mWRC.tensionStepsInMM.length];
			double[] thetaLC_a = new double[mWRC.tensionStepsInMM.length];
			double[] thetaPerc_a = new double[mWRC.tensionStepsInMM.length];
			double[] depthOfPenetrationInMM_a = new double[mWRC.tensionStepsInMM.length];
			
		    double[] averageDistanceFromAeratedPoreInMM = new double[mWRC.tensionStepsInMM.length]; 
		    double[] fractionLessThan3MMFromPhaseBoundary = new double[mWRC.tensionStepsInMM.length]; 
		    double[] fractionMoreThan3MMFromPhaseBoundary = new double[mWRC.tensionStepsInMM.length]; 
		    
			//double[] fractionLessThan5VXFromAeratedPore = new double[mWRC.tensionStepsInMM.length];
			//double[] fraction5to10VXFromAeratedPore = new double[mWRC.tensionStepsInMM.length];
			//double[] fraction10to20VXFromAeratedPore = new double[mWRC.tensionStepsInMM.length];
			//double[] fractionMoreThan10to20VXFromAeratedPore = new double[mWRC.tensionStepsInMM.length];
			
			//find air and water-filled pores	
			for (int j = 0 ; j < mWRC.tensionStepsInMM.length ; j++) {
				
				//set up calculate properties for air phase		
				IJ.showStatus("Calc air phase props at h = " + String.format("%03.0f", mWRC.tensionStepsInMM[j]) + " mm");
				RoiHandler.ColumnRoi airRoi = roi.new ColumnRoi();				
				airRoi.area = colRoi.area;
				airRoi.pRoi = colRoi.pRoi;
				airRoi.voxelSizeInMicron = mWRC.voxelSizeInMicroMeter;	
				airRoi.nowTiff = jOD.extractAirFilledPores(j, drainageTiff, null, mWRC, mFC);
				System.gc();System.gc();
				currInUse = (rt.totalMemory() - rt.freeMemory()) / 1000000000; 
				
				//drainageTiff.updateAndDraw();drainageTiff.show();
				//airRoi.nowTiff.updateAndDraw();airRoi.nowTiff.show();
			    
				//save image
			    if (mWRC.saveAirWaterImages) {			    		    
			    	airwaterTiff = sC.subtract(binTiff, sC.subtractNumber(airRoi.nowTiff, 128));				
			    	if (mWRC.tensionStepsInMM[j] < 0) mFC.fileName = "AnW_At_neg" + String.format("%04.0f",-1 * mWRC.tensionStepsInMM[j]) + "mm.tif";
			    	else mFC.fileName = "AnW_At_" + String.format("%04.0f",mWRC.tensionStepsInMM[j]) + "mm.tif";
			    	jIO.tiffSaver(mFC, airwaterTiff);		
			    	airwaterTiff.unlock();airwaterTiff.flush();
			    	System.gc();System.gc();
			    }
				currInUse = (rt.totalMemory() - rt.freeMemory()) / 1000000000; 
			    
			    //calculate water phase properties
			    IJ.showStatus("Calc water morphs at h = " + String.format("%03.0f", mWRC.tensionStepsInMM[j]) + " mm");
			    RoiHandler.ColumnRoi waterRoi = roi.new ColumnRoi(); 
			    waterRoi.nowTiff = sC.subtract(binTiff, airRoi.nowTiff);
			    waterRoi.area = colRoi.area;
			    waterRoi.pRoi = colRoi.pRoi; 
			    MorphologyAnalyzer.ROIMorphoProps myWaterProps = morph.getSomeSimpleMorphoProps(mFC, waterRoi, mWRC.mRSO);
			    waterRoi.nowTiff.unlock();waterRoi.nowTiff.flush();	
				System.gc();System.gc();
				currInUse = (rt.totalMemory() - rt.freeMemory()) / 1000000000; 
	    
			    //calculate distances to next air-filled pore.
				IJ.showStatus("Calc distances to air at h = " + String.format("%03.0f", mWRC.tensionStepsInMM[j]) + " mm");
			    MorphologyAnalyzer.ROIMorphoProps myDistProps = morph.getDistances2AirFilledPores(mFC, airRoi, mWRC.mRSO);
				System.gc();System.gc();
				currInUse = (rt.totalMemory() - rt.freeMemory()) / 1000000000; 
			    
			    //calculate properties for air phase
			    IJ.showStatus("Calc air morphs at h = " + String.format("%03.0f", mWRC.tensionStepsInMM[j]) + " mm");
				MorphologyAnalyzer.ROIMorphoProps myAirProps = morph.getSomeSimpleMorphoProps(mFC, airRoi, mWRC.mRSO);					
				System.gc();System.gc();
				currInUse = (rt.totalMemory() - rt.freeMemory()) / 1000000000; 
				
 				//populate output structure
 				tensionAtBottomInMM[j] = mWRC.tensionStepsInMM[j];
 				tensionAtCenterInMM[j] = mWRC.tensionStepsInMM[j] + mWRC.columnHeightInMM / 2;
 				tensionAtTopInMM[j] = mWRC.tensionStepsInMM[j] + mWRC.columnHeightInMM;

 				areaInMM2[j] = colRoi.area * mWRC.voxelSizeInMicroMeter * mWRC.voxelSizeInMicroMeter / 1000000;
 				heightInMM[j] = mWRC.columnHeightInMM;
 				bulkVolumeInMM3[j] = myAirProps.roiBulkVolume * (Math.pow(mWRC.voxelSizeInMicroMeter / 1000, 2));
 				
 				theta_w[j] = myWaterProps.phaseVolumeFraction; 				
				sigma_w[j] = myWaterProps.surfaceArea / bulkVolumeInMM3[j];
 				chi_w[j] = myWaterProps.eulerNumber; 				
 				gamma_w[j] = myWaterProps.gamma;
 				fractalDim_w[j] = myWaterProps.surfaceFractalDimension;
 				percolates_w[j] = myWaterProps.phasePercolates; 				
 				thetaLC_w[j] = myWaterProps.largesPhaseClusterVolumeFraction;
 				thetaPerc_w[j] = myWaterProps.percolatingVolumeFraction;
 				
 				theta_a[j] = myAirProps.phaseVolumeFraction;
 				sigma_a[j] = myAirProps.surfaceArea / bulkVolumeInMM3[j];
 				chi_a[j] = myAirProps.eulerNumber;
 				gamma_a[j] = myAirProps.gamma;
 				fractalDim_a[j] = myAirProps.surfaceFractalDimension;
 				percolates_a[j] = myAirProps.phasePercolates;
 				thetaLC_a[j] = myAirProps.largesPhaseClusterVolumeFraction;
 				thetaPerc_a[j] = myAirProps.percolatingVolumeFraction;
 				depthOfPenetrationInMM_a[j] = myAirProps.depthOfPhasePenetration * mWRC.voxelSizeInMicroMeter / 1000;
 				
 				averageDistanceFromAeratedPoreInMM[j] = myDistProps.averageDistance2PhaseBoundary * mWRC.voxelSizeInMicroMeter / 1000; 
 				fractionLessThan3MMFromPhaseBoundary[j] = myDistProps.fractionLessThan3MMFromPhaseBoundary;
 				fractionMoreThan3MMFromPhaseBoundary[j] = myDistProps.fractionMoreThan3MMFromPhaseBoundary;
 				
			}
			
			//populate output structure
			myWRCProperties.tensionAtBottomInMM = tensionAtBottomInMM;
			myWRCProperties.tensionAtCenterInMM = tensionAtCenterInMM;
			myWRCProperties.tensionAtTopInMM = tensionAtTopInMM;

			myWRCProperties.areaInMM2 = areaInMM2;
			myWRCProperties.heightInMM = heightInMM;
			myWRCProperties.bulkVolumeInMM3 = bulkVolumeInMM3;
			
			myWRCProperties.theta_w = theta_w; 				
			myWRCProperties.sigma_w = sigma_w;
			myWRCProperties.chi_w = chi_w; 				
			myWRCProperties.gamma_w = gamma_w;
			myWRCProperties.fractalDim_w = fractalDim_w;
			myWRCProperties.percolates_w = percolates_w; 				
			myWRCProperties.thetaLC_w = thetaLC_w;
			myWRCProperties.thetaPerc_w = thetaPerc_w;
			
			myWRCProperties.theta_a = theta_a;
			myWRCProperties.sigma_a = sigma_a;
			myWRCProperties.chi_a = chi_a;
			myWRCProperties.gamma_a = gamma_a;
			myWRCProperties.fractalDim_a = fractalDim_a;
			myWRCProperties.percolates_a = percolates_a;
			myWRCProperties.thetaLC_a = thetaLC_a;
			myWRCProperties.thetaPerc_a = thetaPerc_a;
			myWRCProperties.depthOfPenetrationInMM_a = depthOfPenetrationInMM_a;
			
			myWRCProperties.averageDistanceFromAeratedPoreInMM = averageDistanceFromAeratedPoreInMM; 
			myWRCProperties.fractionLessThan3MMFromPhaseBoundary = fractionLessThan3MMFromPhaseBoundary;
			myWRCProperties.fractionMoreThan3MMFromPhaseBoundary = fractionMoreThan3MMFromPhaseBoundary;
			
			myWRCProperties.enforceIntegerTensions = mWRC.enforceIntegerTensions;
			
			//save results as ascii
			jIO.writeWRCResultsInASCII(mFC, myWRCProperties);
						
		}
		
		IJ.showMessage("Water retention curves have been reconstructed for all samples!!");
			
	}
}