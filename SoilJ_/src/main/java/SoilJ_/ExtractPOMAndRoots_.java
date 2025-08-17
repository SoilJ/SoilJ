package SoilJ_;

import java.io.File;

import SoilJ.tools.ImageManipulator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.TailoredMaths;

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
 * ExtractFreshOM is a SoilJ plugin aiming at segmenting the fresh soil organic matter and/or water phase in the column
 * and saves the respective binary images in a newly created folder.  
 *  
 * @author John Koestel
 */

public class ExtractPOMAndRoots_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";
		
		//set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = ExtractPOMAndRoots_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ImageManipulator jIM = new ImageManipulator();
		MenuWaiter menu = new MenuWaiter();
		TailoredMaths jTM = new TailoredMaths();

		//init images
		ImagePlus nowTiff = new ImagePlus();		
		ImagePlus outTiff = new ImagePlus();
		
		//ask for threshold choice
		String myOutFolder;
		String myOutPath = null;		
		MenuWaiter.POMThresholderMenuReturn mTMR = menu.showPOMThresholdingDialog();
		if (mTMR == null) return;
		String thresholds = "";
		thresholds = "" + mTMR.minThreshold / 100 + "h";
		thresholds += "_" + mTMR.maxThreshold / 100 + "h";
		
		myOutFolder = "POM" + thresholds + "_BTHTh" + (int)Math.round(mTMR.blackTopHatThreshold / 100) + "h_GradTh" + (int)Math.round(mTMR.gradientThreshold / 100) + "h";
	
		//select file or files
	  	InputOutput.MyFileCollection mFC = jIO.fileSelector("Please choose a file or folder with your image data");
		
		myOutPath = mFC.myBaseFolder + pathSep + myOutFolder;						
		new File(myOutPath).mkdir();
		
		mFC.myOutFolder = myOutPath;

				
		//loop over 3D images
		for (int i = 0 ; i < mFC.myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();
			
			//assign tiff file
			mFC.fileName = mFC.myTiffs[i];
			mFC = jIO.addCurrentFileInfo(mFC);					
				
			//load image
			nowTiff = jIO.openTiff3D(mFC.nowTiffPath);			
			
			//do the cutting
			outTiff = jIM.extractPOMandRoots(nowTiff, mFC, mTMR);	
			
			//save the results
			if (mTMR.save3DImage == true) jIO.tiffSaver(myOutPath, mFC.myTiffs[i], outTiff);
			
			//also save threshold comparison images			
			if (mTMR.save4Evaluation) {
				int[] myZ = jTM.convert2Array(mFC.nOfSlices);
				ImagePlus graySnapTiff = jIM.prepareSnapshots4ThresholdComparison(nowTiff, myZ);
				ImagePlus binSnapTiff = jIM.prepareSnapshots4ThresholdComparison(outTiff, myZ);
				jIO.writeSnapshots4ThresholdComparison(mFC, graySnapTiff, binSnapTiff, myZ);	
			}
						
			if (i == 0) jIO.writeStringIntoAsciiFile(myOutPath + pathSep + "Path2Gauge.txt", mFC.myInnerCircleFolder + "\n"); 	//also save the gauge path..
					
		}
	}
}