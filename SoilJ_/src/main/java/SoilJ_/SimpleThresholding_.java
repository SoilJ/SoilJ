package SoilJ_;

import java.io.File;

import SoilJ.tools.ImageManipulator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;

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
import ij.plugin.PlugIn;

/** 
 * SegmentThis is a SoilJ plugin offering several choices for image segmentation.
 * 
 * @author John Koestel
 *
 */

public class SimpleThresholding_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";
		
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = SimpleThresholding_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ImageManipulator jIM = new ImageManipulator();
		MenuWaiter menu = new MenuWaiter();
		
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();		
		ImagePlus[] outTiff = {null, null, null};

		//ask for threshold choice
		MenuWaiter.ThresholderMenuReturn mTMR = menu.showThresholdingDialog();
		if (mTMR == null) return;		
	
		
		String myOutFolder;
		String myOutPath = null;		
		if (mTMR.useConstantThreshold == true) {
			
			mTMR = menu.showManualThresholdingDialog(mTMR);
			String thresholds = "";
			if (mTMR.minThreshold > 0) thresholds = "" + mTMR.minThreshold;
			if (mTMR.maxThreshold > 0 & thresholds.equalsIgnoreCase("")) thresholds = "" + mTMR.maxThreshold;
			
			myOutFolder = "ConstantThreshold" + thresholds;
			
		}
		else {myOutFolder = mTMR.myPrimaryMethod.toString();		
			if (mTMR.filterImages) myOutFolder = mTMR.filterTag + "_" + myOutFolder;
			if (mTMR.mySecondaryMethod != null) myOutFolder = mTMR.myPrimaryMethod.toString() + "_" + mTMR.mySecondaryMethod.toString();
			if (mTMR.setMaxgray2Wallgray) myOutFolder = myOutFolder + "_WallIsMaxgray";
		}		
		
		//select file or files
	  	InputOutput.MyFileCollection mFC = jIO.fileSelector("Please choose a file or folder with your image data");
	 
		//add inner circle folder location if necessary
		if (mTMR.useInnerCircle) {
			mFC = jIO.addInnerCircleFolder(mFC);
		}
		
		myOutPath = mFC.myBaseFolder + pathSep + myOutFolder;						
		new File(myOutPath).mkdir();
		
		mFC.myOutFolder = myOutPath;
				
		//create folder for each horizontal cross-section for checking
		int[] eDC = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100};  //eDC stands for evaluation depth choices 
			
		//loop over 3D images
		//int errorCounts = 0;
		for (int i = 0 ; i < mFC.myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();
			
			//assign tiff file
			mFC.fileName = mFC.myTiffs[i];
			mFC = jIO.addCurrentFileInfo(mFC);
			
			//find correct innerCircle file..
			if (mTMR.useInnerCircle) {
				String[] GandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC);
				mFC.nowInnerCirclePath = mFC.myInnerCircleFolder + mFC.pathSep + GandS[0];
			}
			
			//define evaluation depths
			int[] myZ = new int[eDC.length];
			for (int j = 0 ; j < eDC.length ; j++) {
				double checkZ = (double)eDC[j] * (double)mFC.nOfSlices / 100d;
				myZ[j] = (int)Math.round(checkZ);
			}
			if (myZ[0] == 0) myZ[0] = 1;
			if (myZ[eDC.length - 1] != mFC.nOfSlices) myZ[eDC.length - 1] = mFC.nOfSlices;
				
			//load image
			if (mTMR.save3DImage | mTMR.save4GeoDict) nowTiff = jIO.openTiff3D(mFC.nowTiffPath);
			else nowTiff = jIO.openTiff3DSomeSlices(mFC, myZ);
			
			//apply segmentation			
			outTiff = jIM.applyAChosenThresholdingMethod(nowTiff, mFC, mTMR);			 
		
			//save result
			if (mTMR.save3DImage == true) jIO.tiffSaver(myOutPath, mFC.myTiffs[i], outTiff[0]);
			if (mTMR.save4GeoDict == true) jIO.saveAsTiffStack4GeoDict(myOutPath, mFC.myTiffs[i], outTiff[0]);
			
			if (i == 0) jIO.writeStringIntoAsciiFile(myOutPath + pathSep + "Path2Gauge.txt", mFC.myInnerCircleFolder + "\n" + mTMR.filterTag); 	//also save the gauge path..			
		 			
		}
		
	}
	
}