package SoilJ_;

import java.io.File;

import SoilJ.tools.ImageManipulator;
import SoilJ.tools.InputOutput;
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
import ij.plugin.PlugIn;

/** 
 * ExtractFreshOM is a SoilJ plugin aiming at segmenting the fresh soil organic matter and/or water phase in the column
 * and saves the respective binary images in a newly created folder.  
 *  
 * @author John Koestel
 */

public class TernarySegmentation_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";
		
		//set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = Extract2DHistograms_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();
		//MenuWaiter menu = new MenuWaiter();
		ImageManipulator jIM = new ImageManipulator();
		
		// init variables
		int i;
		
		//init folder collection
		InputOutput.MyFileCollection mFC = jIO.new MyFileCollection();
		
		//read base folder and number of 3D images
	    String[] myTiffs = null;String myBaseFolder = null;
	    while (myTiffs == null) {
	    	myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
			if (myBaseFolder == null) return;
			myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));
			if (myTiffs == null) {
				IJ.error("Please choose a folder with TIFF images or cancel. Thank you.");	
			}
		}
	    if (myBaseFolder.substring(myBaseFolder.length() - 1).equalsIgnoreCase(pathSep)) {
	    	myBaseFolder = myBaseFolder.substring(0, myBaseFolder.length() - 1);
	    }
	    
	    String myGradFolder = "";
    	if (jIO.testIfFolderIsPresent(new File(myBaseFolder), "Gradients")) {
    		myGradFolder = myBaseFolder + pathSep + "Gradients";
    	}
    	else {	    	
    		myTiffs = null;
			while (myTiffs == null) {
				myGradFolder = jIO.chooseAFolder("Please choose the folder with your gradient images");
				if (myGradFolder == null) return;
				myTiffs = jIO.listTiffsInFolder(new File(myGradFolder));
				if (myTiffs == null) {				
					IJ.error("Please choose a folder with TIFF images or cancel. Thank you.");
				}
			}		
    	}
    	    			    
	    //create output paths	
  		mFC.myBaseFolder = myBaseFolder;
  		
  		String myOutPath =  myBaseFolder + pathSep + "TernarySegmentation";
  		new File(myOutPath).mkdir();		
  		mFC.myOutFolder = myOutPath; 	//also add it to the folder collection
  		  		  		
  		//save pathSep
  		mFC.pathSep = pathSep;

		//remember gradient folder
		mFC.myGradientFolder = myGradFolder;
		
		//load pomRegion
		ImagePlus pomRegion = new ImagePlus();
		String nowDir = IJ.getDirectory("current");
		mFC.nowTiffPath = jIO.chooseAFile("Please choose the file with your POM seed region mask", nowDir, "POMRegionMask.tif");
		pomRegion = jIO.openTiff3D(mFC);
		
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  

			//assign current TIFF
			mFC.fileName = myTiffs[i];
			mFC = jIO.addCurrentFileInfo(mFC);
			mFC.nowTiffPath = mFC.myBaseFolder + pathSep + mFC.fileName;
			ImagePlus nowTiff = jIO.openTiff3D(mFC.nowTiffPath);
						
			//load or create gradient Image		
			mFC.nowTiffPath = mFC.myGradientFolder + pathSep + mFC.fileName;
			ImagePlus gradTiff = jIO.openTiff3D(mFC.nowTiffPath);	
			
			//segment ternary and save
			jIM.segmentTernaryAndSave(mFC, nowTiff, gradTiff, pomRegion);
		
		}
		
		IJ.showStatus("Done! Flushing memory.. ");
		
		//try to clear memory
		IJ.freeMemory();IJ.freeMemory();
		
		IJ.showStatus("Memory should be free again!");
		
	}
}