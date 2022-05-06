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
 * BeamDeHardening is a SoilJ plugin aiming at the removal of concentric imaging artifacts
 * as they are occurring due to beam hardening. 
 * 
 * @author John Koestel
 *
 */

public class ExtractVoxelWiseWRC_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {
		
		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";

		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = BeamDeHardening_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ImageManipulator jIM = new ImageManipulator();
								
		// init variables
		int i;
		
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();
			
		//ask for threshold choice
	    //MenuWaiter.BeamDeHardeningReturn mBDH = menu.showBeamDeHardeningMenu();
	    //if (mBDH == null) return;
			    
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
		
		//read gauges in the gauge folder
		boolean useInnerCircle = true;
		InputOutput.MyFileCollection mFC = jIO.getAllMyNeededFolders(myBaseFolder, myTiffs, "", false, false);

		//save and create output folder		
		String myOutFolder;
		String myOutPath = null;
		double scalingFactor = 0.05;
		myOutFolder = "ScaledBy1over" + (int)(1 / scalingFactor);			
		myOutPath = myBaseFolder + myOutFolder;
		new File(myOutPath).mkdir();
		String myThetaFolder = myOutPath + pathSep + "theta";
		new File(myThetaFolder).mkdir();
				
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();		
			
			//load image
			nowTiff = jIO.openTiff3D(myBaseFolder + pathSep + myTiffs[i]);	
			
			//scale by factor 10			
			IJ.run(nowTiff , "Scale...", "x=0.05 y=0.05 z=0.05 width=26 height=26 depth=17 interpolation=None average process create");
			ImagePlus smallTiff = IJ.getImage();
			smallTiff.hide();
			
			//nowTiff = jIM.scaleIsotropicly3D(nowTiff, scalingFactor);			
									
			//save image result
			jIO.tiffSaver(myOutPath, myTiffs[i], smallTiff);
			
			//save WRC data
		}
		
		System.gc();
		
		//subtract ddry
		String dryPath = myOutPath + pathSep + "KB1_pdry.tif";
		ImagePlus dryTiff = jIO.openTiff3D(dryPath);	
		
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  //myTiffs.length
	
			//try to free up some memory
			System.gc();		
			
			//load image
			if (!myTiffs[i].equalsIgnoreCase("KB1_pdry.tif")) {
				
				nowTiff = jIO.openTiff3D(myOutPath + pathSep + myTiffs[i]);
			
				//SUBTRACT  ddry
				ImagePlus diffTiff = jIM.calc3D(nowTiff, "-", dryTiff);
				
				//divide by 5000... the value for water!
				IJ.run(diffTiff, "32-bit", "false");
				IJ.run(diffTiff, "Divide...", "value=5000 stack");		
				//diffTiff = IJ.getImage();
				//diffTiff.hide();
				
				//save image result
				jIO.tiffSaver(myThetaFolder, myTiffs[i], diffTiff);
				
			}
		}
	}
}