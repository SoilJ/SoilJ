package SoilJ_;

import java.io.File;

import SoilJ.tools.ImageManipulator;
import SoilJ.tools.ImageManipulator.StackCalculator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.InputOutput.MyFileCollection;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.MorphologyAnalyzer;
import SoilJ.tools.MorphologyAnalyzer.WRCPhaseProps;
import SoilJ.tools.ObjectDetector;
import SoilJ.tools.RoiHandler;

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

public class SimulateDrainage_ extends ImagePlus implements PlugIn  {

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

		
		// init variables
		int i;
		MenuWaiter.DrainageSimulatorOptions mDS = menu.new DrainageSimulatorOptions();
		mDS = menu.showDrainageSimulatorMenu();
		if (mDS == null) return;
		
		//construct image related objects
		ImagePlus airTiff;  
		ImagePlus waterTiff;
		ImagePlus airwaterTiff;  //output
		
		//read base folder and number of 3D images
		MyFileCollection mFC = jIO.fileSelector("Please choose a file or folder with your image data");				
						
		//create output folder
		String myOutFolder = "DrainageSimulation";
		String mySubBaseFolder = jIO.getTheFolderAbove(mFC.myBaseFolder, pathSep);
		String myOutPath = mySubBaseFolder + pathSep + myOutFolder;
		new File(myOutPath).mkdir();		
		mFC.myOutFolder = myOutPath;
						
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
			int[] startStopSlices = jIO.findStartAndStopSlices(mFC, mDS.mRSO);
			int[] colSlices = new int[startStopSlices[1] - startStopSlices[0] + 1];
			for (int j = 0 ; j < colSlices.length ; j++) colSlices[j] = startStopSlices[0] + j;
		
			//cut roi		
			mFC.startSlice = startStopSlices[0];
			mFC.stopSlice = startStopSlices[1];
			RoiHandler.ColumnRoi colRoi = roi.prepareDesiredRoi(mFC, jIO.openTiff3DSomeSlices(mFC, colSlices), mDS.mRSO, "null", "null");
						
			//calculate  at top, center and bottom of column
			mDS.columnHeightInMM = mDS.voxelSizeInMicroMeter * colRoi.nowTiff.getNSlices() / 1000;
			
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
			
			//binTiff.updateAndDraw();binTiff.show();
			
			//init output structure
			WRCPhaseProps myDrainageProps = morph.new WRCPhaseProps();
			
			//init output arrays
			double[] pressureAtBottomInMM = new double[mDS.pressureStepsInMM.length];
			double[] pressureAtCenterInMM = new double[mDS.pressureStepsInMM.length];
			double[] pressureAtTopInMM = new double[mDS.pressureStepsInMM.length];

			double[] areaInMM2 = new double[mDS.pressureStepsInMM.length];
			double[] heightInMM = new double[mDS.pressureStepsInMM.length];
			double[] bulkVolumeInMM3 = new double[mDS.pressureStepsInMM.length];
			
			double[] theta_w = new double[mDS.pressureStepsInMM.length];			
			double[] sigma_w = new double[mDS.pressureStepsInMM.length];
			double[] chi_w = new double[mDS.pressureStepsInMM.length];			
			double[] gamma_w = new double[mDS.pressureStepsInMM.length];
			double[] fractalDim_w = new double[mDS.pressureStepsInMM.length];
			double[] percolates_w = new double[mDS.pressureStepsInMM.length]; 				
			double[] thetaLC_w = new double[mDS.pressureStepsInMM.length];
			double[] thetaPerc_w = new double[mDS.pressureStepsInMM.length];
			
			double[] theta_a = new double[mDS.pressureStepsInMM.length];
			double[] sigma_a = new double[mDS.pressureStepsInMM.length];
			double[] chi_a = new double[mDS.pressureStepsInMM.length];
			double[] gamma_a = new double[mDS.pressureStepsInMM.length];
			double[] fractalDim_a = new double[mDS.pressureStepsInMM.length];
			double[] percolates_a = new double[mDS.pressureStepsInMM.length];
			double[] thetaLC_a = new double[mDS.pressureStepsInMM.length];
			double[] thetaPerc_a = new double[mDS.pressureStepsInMM.length];
			double[] depthOfPenetrationInMM_a = new double[mDS.pressureStepsInMM.length];
						
			//find air and water-filled pores
			RoiHandler.ColumnRoi airRoi = roi.new ColumnRoi();
			RoiHandler.ColumnRoi waterRoi = roi.new ColumnRoi();
			airRoi.area = colRoi.area;
			airRoi.pRoi = colRoi.pRoi;
			waterRoi.area = colRoi.area;
			waterRoi.pRoi = colRoi.pRoi;
			for (int j = 0 ; j < mDS.pressureStepsInMM.length ; j++) {
				
				//extract air-filled pores 
				if (-mDS.pressureStepsInMM[j] <= -mDS.columnHeightInMM) airRoi.nowTiff = sC.subtractNumber(binTiff, 255);  //if all pores in the image will be water-filled, don't bother to calculate it explicitly.
				else airRoi.nowTiff = jOD.extractAirFilledPores(-mDS.pressureStepsInMM[j], colRoi.nowTiff, null, mDS, mFC);
				
				//airTiff.updateAndDraw();airTiff.show();
				
				//extract water-filled pores
				waterRoi.nowTiff = sC.subtract(binTiff, airRoi.nowTiff);
				
				//calculate properties for air phase
			    MorphologyAnalyzer.ROIMorphoProps myAirProps = morph.getSomeSimpleMorphoProps(mFC, airRoi, mDS.mRSO);
				
				//calculate properties for air phase
				MorphologyAnalyzer.ROIMorphoProps myWaterProps = morph.getSomeSimpleMorphoProps(mFC, waterRoi, mDS.mRSO);
		
				//try to free up some memory
 				System.gc();System.gc();
 				
 				//populate output structure
 				pressureAtBottomInMM[j] = mDS.pressureStepsInMM[j];
 				pressureAtCenterInMM[j] = mDS.pressureStepsInMM[j] - mDS.columnHeightInMM / 2;
 				pressureAtTopInMM[j] = mDS.pressureStepsInMM[j] - mDS.columnHeightInMM;

 				areaInMM2[j] = colRoi.area * mDS.voxelSizeInMicroMeter * mDS.voxelSizeInMicroMeter / 1000000;
 				heightInMM[j] = mDS.columnHeightInMM;
 				bulkVolumeInMM3[j] = myAirProps.roiBulkVolume * (Math.pow(mDS.voxelSizeInMicroMeter / 1000, 2));
 				
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
 				depthOfPenetrationInMM_a[j] = myAirProps.depthOfPhasePenetration * mDS.voxelSizeInMicroMeter / 1000;
 				
				//save image
 				if (mDS.saveAirWaterImages) {
 					airwaterTiff = sC.subtract(binTiff, sC.subtractNumber(airRoi.nowTiff, 128));
 					if (mDS.pressureStepsInMM[j] < 0) mFC.fileName = "AnW_At_neg" + String.format("%04.0f",-1 * mDS.pressureStepsInMM[j]) + "mm.tif";
 					else mFC.fileName = "AnW_At_" + String.format("%04.0f",mDS.pressureStepsInMM[j]) + "mm.tif";
 					jIO.tiffSaver(mFC, airwaterTiff);
 				}
 				
			}
			
			//populate output structure
			myDrainageProps.pressureAtBottomInMM = pressureAtBottomInMM;
			myDrainageProps.pressureAtCenterInMM = pressureAtCenterInMM;
			myDrainageProps.pressureAtTopInMM = pressureAtTopInMM;

			myDrainageProps.areaInMM2 = areaInMM2;
			myDrainageProps.heightInMM = heightInMM;
			myDrainageProps.bulkVolumeInMM3 = bulkVolumeInMM3;
			
			myDrainageProps.theta_w = theta_w; 				
			myDrainageProps.sigma_w = sigma_w;
			myDrainageProps.chi_w = chi_w; 				
			myDrainageProps.gamma_w = gamma_w;
			myDrainageProps.fractalDim_w = fractalDim_w;
			myDrainageProps.percolates_w = percolates_w; 				
			myDrainageProps.thetaLC_w = thetaLC_w;
			myDrainageProps.thetaPerc_w = thetaPerc_w;
			
			myDrainageProps.theta_a = theta_a;
			myDrainageProps.sigma_a = sigma_a;
			myDrainageProps.chi_a = chi_a;
			myDrainageProps.gamma_a = gamma_a;
			myDrainageProps.fractalDim_a = fractalDim_a;
			myDrainageProps.percolates_a = percolates_a;
			myDrainageProps.thetaLC_a = thetaLC_a;
			myDrainageProps.thetaPerc_a = thetaPerc_a;
			myDrainageProps.depthOfPenetrationInMM_a = depthOfPenetrationInMM_a;
			
			myDrainageProps.enforceIntegerTensions = mDS.enforceIntegerPressures;
			
			//save results as ascii
			jIO.writeDrainageSimulationResultsInASCII(mFC, myDrainageProps);
						
		}
			
	}
}
