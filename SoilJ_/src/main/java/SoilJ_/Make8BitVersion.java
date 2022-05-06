package SoilJ_;

import java.io.File;

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
import ij.ImageStack;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

/** 
 * Make8BitVersion is a SoilJ plugin that converts a 16 bit grey scale image into an 8 bit one.
 * 
 * @author John Koestel
 *
 */

public class Make8BitVersion extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";
		
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = Make8BitVersion.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
	
		// init variables
		int i;
	
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();
		
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
		String myOutFolder = "8bit";		
		String myOutPath = myBaseFolder + myOutFolder;
		new File(myOutPath).mkdir();
		
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();
			
			//load image
			nowTiff = jIO.openTiff3D(myBaseFolder + pathSep + myTiffs[i]);	
			ImagePlus outTiff = new ImagePlus();
			ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
			
			for (int j = 1 ; j < nowTiff.getNSlices() + 1 ; j++) {
				
				nowTiff.setSlice(j);
				ImageProcessor nowIP = nowTiff.getProcessor();
				nowIP.max(1);
				nowIP.multiply(255);
				
				ImageProcessor binIP = nowIP.convertToByteProcessor(false);
				
				outStack.addSlice(binIP);
				
			}
			
			outTiff.setStack(outStack);
										
			//apply segmentation	
			jIO.tiffSaver(myOutPath, myTiffs[i], outTiff);
			
		}
	}
}