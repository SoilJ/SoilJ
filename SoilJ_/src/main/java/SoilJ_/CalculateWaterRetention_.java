package SoilJ_;

import java.io.File;

import SoilJ.tools.HistogramStuff;
import SoilJ.tools.ImageManipulator;
import SoilJ.tools.ImageManipulator.StackCalculator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.ObjectDetector;
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

public class CalculateWaterRetention_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {
		
		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";

		//set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = CalculateWaterRetention_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ObjectDetector jOD = new ObjectDetector();
		ImageManipulator jIM = new ImageManipulator();
		MenuWaiter menu = new MenuWaiter();
		HistogramStuff hist = new HistogramStuff();
		
		// init variables
		int i;
		MenuWaiter.WRCCalculatorMenu mWRC = menu.new WRCCalculatorMenu();
		mWRC = menu.showWRCMenu();
		if (mWRC == null) return;
		
		//construct image related objects
		ImagePlus nowTiff;  //input
		ImagePlus airTiff;  //zwischiput
		ImagePlus waterTiff;  //output
		
		//read base folder and number of 3D images
	    String[] myTiffs = null;String myBaseFolder = null;
	    while (myTiffs == null) {
	    	myBaseFolder = jIO.chooseAFolder("Please choose the folder with thickness images.");
			if (myBaseFolder == null) return;
			myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));
			if (myTiffs == null) {
				IJ.error("Please choose a folder with TIFF images or cancel. Thank you.");	
			}
		}
		String myTiffName;
		String myOutFolder = "WaterRetentionData";
						
		//create output folder
		String mySubBaseFolder = jIO.getTheFolderAbove(myBaseFolder, pathSep);
		String myOutPath = mySubBaseFolder + pathSep + myOutFolder;
		new File(myOutPath).mkdir();
		
		//make sub directories for the tension steps
		String[] dirNames = new String[mWRC.tensionSteps.length];
		for (i = 0 ; i < mWRC.tensionSteps.length ; i++) {
			dirNames[i] = myOutPath + pathSep + String.format("%2.1f",mWRC.tensionSteps[i]) + "cm";
			dirNames[i].replace(',', '.');
			new File(dirNames[i]).mkdir();
		}
		
		//init information storage files 
		InputOutput.MyFileCollection mFC = jIO.new MyFileCollection();
		mFC.myOutFolder = myOutPath;
		mFC.myBaseFolder = myBaseFolder;
		
		//calculate tension at top, center and bottom of column
		if (mWRC.choiceOfWRCType.equalsIgnoreCase("draining")) {
			for (i = 0 ; i < mWRC.tensionSteps.length ; i++) {
				mWRC.tensionAtBottom[i] = mWRC.tensionSteps[i];
				mWRC.tensionAtCenter[i] = mWRC.tensionSteps[i] + mWRC.columnHeight / 2;
				mWRC.tensionAtTop[i] = mWRC.tensionSteps[i] + mWRC.columnHeight;
			}
		}
				
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  

			//try to free up some memory
			System.gc();System.gc();
			
			//assign current TIFF
			myTiffName = myTiffs[i];			
			
			//assign tiff file
			mFC.fileName = myTiffs[i];
			mFC = jIO.addCurrentFileInfo(mFC);
			
			nowTiff = jIO.openTiff3D(myBaseFolder + pathSep + myTiffName);
			
			//make a binarized version of the current Tiff..
			ImageStack binStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
			ImagePlus binTiff = new ImagePlus();
			for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
				nowTiff.setPosition(z + 1);
				ImageProcessor thIP = nowTiff.getProcessor();
				ImageProcessor cpIP = thIP.duplicate();
				cpIP.max(2);
				ImageProcessor binIP = cpIP.convertToByte(true);
				binIP.threshold(1);
				binStack.addSlice(binIP);
			}
			binTiff.setStack(binStack);
			
			//calculate pore volume in voxel
			int[] phiHist = hist.sampleHistogram(binTiff);
			int phi = phiHist[255];
			
			int[] air = new int[mWRC.tensionSteps.length];
			int[] water = new int[mWRC.tensionSteps.length];
			double[] wSM = new double[mWRC.tensionSteps.length];
			double[] aSM = new double[mWRC.tensionSteps.length];
			
			//binTiff.updateAndDraw();
			//binTiff.show();
			
			IJ.error("Does not take sirface topography into account for now. Sorry.");
			
			//find water-filled pores
			for (int j = 0 ; j < mWRC.tensionSteps.length ; j++) {
				
				//extract air-filled pores
				airTiff = jOD.extractAirFilledPores(mWRC.tensionSteps[j], nowTiff, null, mWRC, mFC);
				
				//airTiff.updateAndDraw();
				//airTiff.show();
				
				//extract water-filled pores
				StackCalculator sC= jIM.new StackCalculator();
				waterTiff = sC.subtract(binTiff, airTiff);
				
				//waterTiff.updateAndDraw();
				//waterTiff.show();
				
				//calculate water content etc..	
				int[] airHist = hist.sampleHistogram(airTiff);
				int[] waterHist = hist.sampleHistogram(waterTiff);
								
				air[j] = airHist[255];
				water[j] = waterHist[255];
				wSM[j] = (float)water[j] / (float)phi;
				aSM[j] = (float)air[j] / (float)phi;
				
				//save images
				if (mWRC.saveAirImages) jIO.tiffSaver(dirNames[j], "AirAt" + String.format("%2.1f\t",mWRC.tensionSteps[j]) + "cm_" + myTiffName, airTiff);
				if (mWRC.saveWaterImages) jIO.tiffSaver(dirNames[j], "WaterAt" + String.format("%2.1f\t",mWRC.tensionSteps[j]) + "cm_" + myTiffName, waterTiff);
				
				//try to free up some memory
 				System.gc();System.gc();
			}
			
			//save results as ascii
			jIO.writeWRCResultsInASCII(mFC, mWRC, air, water, wSM, aSM);
						
		}
			
	}
}