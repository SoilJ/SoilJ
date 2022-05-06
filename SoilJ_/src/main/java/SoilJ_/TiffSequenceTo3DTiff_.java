package SoilJ_;

import java.io.File;
import java.io.FilenameFilter;

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
import ij.io.FileInfo;
import ij.io.Opener;
import ij.plugin.PlugIn;

/** 
 * CombineTiffStack2Tiff is a SoilJ plugin that fuses all TIFFs contained in subfolders into
 * 3-D TIFFs and saves the files into a newly created folder named "3-D"
 * 
 * @author John Koestel
 *
 */

public class TiffSequenceTo3DTiff_ extends ImagePlus implements PlugIn {

	// loads a stack of tiff images from the Offer sampling campaign (June 2013)
	// and aligns them as such that they stand 'perfectly' straight

	public void run(String arg) {
		
		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";
		
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = TiffSequenceTo3DTiff_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//init objects		
		InputOutput jIO = new InputOutput();
	
		ImagePlus imgStack;
		
		String myBaseFolder;
		String mySubBaseFolder;
		String[] myFolders;
		String myOutFolder = "Tiff3D";
		String myOutPath;
		int i;
		
		//read base folder and number of 3D images
		myFolders = null; myBaseFolder = null;
		
		while (myFolders == null) {
			myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
			if (myBaseFolder == null) return;
			myFolders = jIO.listFoldersInFolder(new File(myBaseFolder));
			if (myFolders == null) {
				IJ.error("Please choose a folder containing folders with TIFF image-sequences or cancel. Sorry if I have confused you.");	
			}
			else { 		//check if this folder contains folders with single layer tiffs
				
				File myFolder = new File(myBaseFolder + pathSep + myFolders[0]);
				
				FilenameFilter tiffFilter = new FilenameFilter() {
					public boolean accept(File dir, String name) {
						String lowercaseName = name.toLowerCase();
						if (lowercaseName.endsWith(".tif") || lowercaseName.endsWith(".tiff")) {
							return true;
						} else {
							return false;
						}
					}
				};
				
				File[] myTiffs = myFolder.listFiles(tiffFilter);
				int numOTiffs = myTiffs.length;
						
				//sample first and last image
				FileInfo[] fI1 = Opener.getTiffFileInfo(myTiffs[0].getPath());		 
				FileInfo[] fI2 = Opener.getTiffFileInfo(myTiffs[numOTiffs - 1].getPath());		
				
				if (fI1[0].nImages > 1) {					
					IJ.error("Computer says that file " + fI1[0].fileName + " is a 3-D TIFF containing " + fI1[0].nImages + " 2-D slices.\n"
							+ "Please select a folder containing exclusively folders with 2-D TIFFs.\nThank you.");
					myFolders = null;
				}
				
				else if (fI2[0].nImages > 1) {					
					IJ.error("Computer says that file " + fI2[0].fileName + " is a 3-D TIFF containing " + fI2[0].nImages + " 2-D slices.\n"
							+ "Please select a folder containing exclusively folders with 2-D TIFFs.\nThank you.");
					myFolders = null;
				}
				
				if (myFolders != null) {
					if (fI1[0].width != fI2[0].width || fI1[0].height != fI2[0].height) {
						IJ.error("The canvas size must be identical or all images. Thank you.");
						myFolders = null;break;
					}
				}	
				
			}
		}
		
		//create output folder		
		mySubBaseFolder = jIO.getTheFolderAbove(myBaseFolder, pathSep);
		myOutPath = mySubBaseFolder + pathSep + myOutFolder;
		new File(myOutPath).mkdir(); 
		
		//go through each tiffstack and save it as one 3D tiff in the output folder
		for (i = 0 ; i < myFolders.length ; i++) {
			
			File myFolder = new File(myBaseFolder + pathSep + myFolders[i]);
			imgStack = jIO.openStack(myFolder.getAbsolutePath());
			jIO.tiffSaver(myOutPath, myFolders[i] + ".tif", imgStack);						
					
		}	
	}
}
	