package SoilJ_;

import java.io.File;

import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.MorphologyAnalyzer;
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

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;

/** 
 * CreateAnAggregateMask is a SoilJ plugin that creates a mask defining outlines of soil aggregates.
 * It was developed by Tobias Bölshcwer and John Koestel. It makes use of the ImageJ plugins BoneJ (fast particle analyzer) and
 * MorphoLibJ (3D watershed distance transform).
 * 
 * The results of the analyses are written into image and ASCII files. 
 * 
 * @author John Koestel
 *
 */

public class CreateClosingMask_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";
		
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = PoreSpaceAnalyzer_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();
		ObjectDetector jOD = new ObjectDetector();
		MenuWaiter menu = new MenuWaiter();
		RoiHandler roi = new RoiHandler();
		MorphologyAnalyzer morph = new MorphologyAnalyzer();
		
		// init variables
		int i;
		
		//tell me what I should do!
		MenuWaiter.ClosingMaskOptionMenu aMO = menu.showClosingMaskOptionMenu();
		if (aMO == null) return;
		
		InputOutput.MyFileCollection mFC = jIO.new MyFileCollection();
		
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();
		ImagePlus closingTiff = new ImagePlus();
		
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
					
		//create output paths
		String myPreOutFolder = "";		
		myPreOutFolder = "ClosMask_Clos" + aMO.closingVoxels + "vxs";
		
		//populate data information file	
		mFC = jIO.getAllMyNeededFolders(myBaseFolder, myTiffs, "", aMO.useInnerCircleFiles, false);   //"" put because no possibility for a filterTag implemented yet

		//add basic folders to folder collection
		mFC.myBaseFolder = myBaseFolder;
		mFC.myTiffs = myTiffs;
		
		String myPreOutPath = myBaseFolder + myPreOutFolder;
		new File(myPreOutPath).mkdir();		
		mFC.myPreOutFolder = myPreOutPath; 	//also add it to the folder collection
			
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory			
			IJ.freeMemory();IJ.freeMemory();		
			
			//assemble image for analyses
			mFC.fileName = myTiffs[i];
			mFC = jIO.addCurrentFileInfo8Bit(mFC);
						
			//load file			
			nowTiff = jIO.openTiff3D(mFC);
			
			//cut image
			RoiHandler.ColumnRoi colRoi = roi.assembleInnerCircleROIs(mFC, nowTiff, aMO.useInnerCircleFiles, aMO.extraInnerCircle, aMO.cutCanvas);					
						
			//try to free up some memory			
			IJ.freeMemory();IJ.freeMemory();		
			
			//apply analyzes
			closingTiff = jOD.makeAggregateMask(mFC, colRoi, aMO);		
			
			//save file
			jIO.tiffSaver(mFC.myPreOutFolder, mFC.myTiffs[i], closingTiff);
			
		}
				
		//try to free up some memory			
		IJ.freeMemory();IJ.freeMemory();	
		
	}
}