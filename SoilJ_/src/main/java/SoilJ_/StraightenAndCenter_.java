package SoilJ_;

import java.io.File;

import SoilJ.tools.ImageManipulator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.ObjectDetector;

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
import ij.Prefs;
import ij.plugin.PlugIn;

/** 
 * StraightenAndCenter is a SoilJ class aiming at detecting the location adn diameter of a cylindrical soil column
 * as well as its inclination. The column is thereupon rotated into an upright position and unsused canvas in the XY-planes
 * are clipped away from the image. 
 * 
 * @author John Koestel
 *
 */

public class StraightenAndCenter_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";
		
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = StraightenAndCenter_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//fix memory settings		
		Prefs.keepUndoBuffers=false;		
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ObjectDetector jOD = new ObjectDetector();
		ImageManipulator jIM = new ImageManipulator();
		MenuWaiter menu = new MenuWaiter();
						
		// init variables
		int i;
		
		//construct image related objects
		ImagePlus upTiff = new ImagePlus();  //perfectly upright columns
					
		//get details on PVC column finding
		MenuWaiter.ColumnFinderMenuReturn jCFS = menu.new ColumnFinderMenuReturn(); 
		jCFS = menu.showColumnStraightenerMenu();
		if (jCFS == null) return;
		
		//select file or files
	  	InputOutput.MyFileCollection mFC = jIO.fileSelector("Please choose a file or folder with your image data");
	 
	  	//Create output folder
		String myOutFolder = "StraightAndCentered";
		
		//create output folder		
		String myOutPath = mFC.myBaseFolder + pathSep + myOutFolder;
		new File(myOutPath).mkdir();
		
		//loop over 3D images
		int errorCounts = 0;
		for (i = 0 ; i < mFC.myTiffs.length ; i++) {  

			//assign current TIFF	
			mFC.fileName = mFC.myTiffs[i];
			jIO.addCurrentFileInfo(mFC);
			
			//check if everything is in order
			if (mFC.columnHasNotBeenFound) {
				IJ.error(mFC.eMsg);
				return;
			}
			
			//load file
			ImagePlus nowTiff = jIO.openTiff3D(mFC.myBaseFolder + pathSep + mFC.fileName);
			
			IJ.freeMemory();IJ.freeMemory();
			
			try {
				//find the inner edge of the PVC column
				boolean look4PreciseCoords = false;
				ObjectDetector.ColCoords3D prelimCC = jOD.findOrientationOfPVCOrAluColumn(nowTiff, jCFS, look4PreciseCoords);
				
				IJ.freeMemory();IJ.freeMemory();	
						
				//check if column is already upright					
				upTiff = jIM.putColumnUprightInCenter(nowTiff, prelimCC, jCFS);
				
				nowTiff.flush();
				IJ.freeMemory();IJ.freeMemory();	
			
				//save the results
				jIO.tiffSaver(myOutPath, mFC.fileName, upTiff);				
				
			}
			catch(Exception e){
				
				String eMsg = "Something went wrong when trying to find column '" + mFC.fileName + "'!\n";
				
				errorCounts++;
				
				IJ.log(eMsg);				
							
			}
			
			upTiff.flush();
			IJ.freeMemory();IJ.freeMemory();IJ.freeMemory();IJ.freeMemory();					
						
		}
		
		if (errorCounts > 0) {
			
			String eMsg = "\n\nThis plugin requires images on which the entire outer perimeter of the column wall\n";		
			eMsg += "is visible. If the wall is only partially visible, the plugin will probably not work.\n\n";
			eMsg += "In the case that the unused canvas in horizontal cross-sections is very large,\n"; 
			eMsg += "try cutting away parts of it before trying this plugin again. \n";
			eMsg += "Make sure that the column is more or less in the center of the XY-plane after cutting the canvas\n";
			eMsg += "of the XY-plane.\n\n";
			eMsg += "In the case that the image is very noisy, use a filter to reduce the image noise\n";
			eMsg += "Suitable filters are a non-local means filter or, alternatively, a median filter \n";			
			eMsg += "followed by an unsharp mask.\n\n";
			eMsg += "Thank you for using SoilJ!";			
			
			IJ.log(eMsg);
		}
	}
}