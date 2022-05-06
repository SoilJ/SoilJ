package SoilJ_;

import java.io.File;

import SoilJ.tools.InputOutput;
import SoilJ.tools.MorphologyAnalyzer;
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
 * ExtractPercolatingPorosity is a SoilJ plugin that extracts the percolating phase from a binary image.
 * 
 * @author John Koestel
 *
 */

public class ExtractPercolatingPorosity extends ImagePlus implements PlugIn  {

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
		//MenuWaiter menu = new MenuWaiter();
		//RoiHandler roi = new RoiHandler();
		MorphologyAnalyzer morph = new MorphologyAnalyzer();
		
		// init variables
		int i;
		
		InputOutput.MyFileCollection mFC = jIO.new MyFileCollection();
		
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();
		
		//read base folder and number of 3D images
	    String[] myTiffs = null;String myBaseFolder = null;
	    while (myTiffs == null) {
	    	myBaseFolder = jIO.chooseAFolder("Please choose the folder with your cluster label images");
			if (myBaseFolder == null) return;
			myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));
			if (myTiffs == null) {
				IJ.error("Please choose a folder with TIFF images or cancel. Thank you.");	
			}
		}

     	mFC.myBaseFolder = myBaseFolder;
		mFC.myTiffs = myTiffs;
		
		String myPreOutPath = myBaseFolder;	
		mFC.myPreOutFolder = myPreOutPath; 	//also add it to the folder collection
		
		double[] classBounds = {2,2.8,3.4,4,4.4,4.8,5.6,6,6.3,6.6,6.9,7.2,7.4,8,8.2,8.4,8.7,8.9,9.1,9.3,
				9.7,10,10.1,10.3,10.7,10.9,11.3,11.4,11.6,12,12.5,13,13.5,14,14.5,15,15.5,16,
				16.5,17,17.5,18,19,20,21,22,23,24,25,26,27,28,29,30,32,34,36,38,40,
				42,44,46,48,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,160,170,180,190,200,
				210,220,230,240,250,10000};
		
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();
			
			//load image
			nowTiff = jIO.openTiff3D(myBaseFolder + pathSep + myTiffs[i]);	
									
			String myName = myTiffs[i];
			myName = myName.substring(0,myName.length()-4) + ".psd";
												
			//apply segmentation							}
			morph.extractPoresizeDistro(myPreOutPath + pathSep + myName, nowTiff, classBounds);
			
		}
		
	}
}