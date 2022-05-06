package SoilJ_;

import java.io.File;

import SoilJ.tools.DisplayThings;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.MorphologyAnalyzer;
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

import ij.ImagePlus;
import ij.plugin.PlugIn;

/** 
 * SoilSurfaceFinder is a SoilJ plugin aiming at detecting the top and bottom topographies of soil columns. 
 * The results are saved in two-layer TIFFs depicting the distance of the topmost (bottommost) soil solid phase voxel 
 * from the top (bottom) of the 3-D image canvas. The topography of the soil top surface is sved as the 
 * first of the two images.
 * 
 * @author John Koestel
 *
 */

public class SurfaceDetection_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";
		
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = SurfaceDetection_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ObjectDetector jOD = new ObjectDetector();
		MenuWaiter menu = new MenuWaiter();		
		DisplayThings dT = new DisplayThings();			
		
		// init variables
		int i;
	
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();	
		ImagePlus outTiff = new ImagePlus();
		
		//show menu for SoilSurfaceFinder
		MenuWaiter.SurfaceFinderReturn mSFR = menu.new SurfaceFinderReturn();
		mSFR = menu.showSurfaceFinderMenu();
		if (mSFR == null) return;
		
		//select file or files
	  	InputOutput.MyFileCollection mFC = jIO.fileSelector("Please choose a file or folder with your image data");
		
		//create output path for top surface
		String myOutFolder = "SurfaceOfColumn";
		String myOutPath = mFC.myBaseFolder + pathSep + myOutFolder;
		new File(myOutPath).mkdir();
		mFC.myOutFolder = myOutPath;
		
		//init surface statistics output
		MorphologyAnalyzer.SurfaceStatistics[] mSS = new MorphologyAnalyzer.SurfaceStatistics[mFC.myTiffs.length]; 
		
		//loop over 3D images
		for (i = 0 ; i < mFC.myTiffs.length ; i++) {  //myTiffs.length		

			//try to free up some memory
			System.gc();
			
			//load image
			nowTiff = jIO.openTiff3D(mFC.myBaseFolder + pathSep + mFC.myTiffs[i]);	
			mFC.fileName = mFC.myTiffs[i];
			mFC.mySurfaceFolder = null;
			mFC = jIO.addCurrentFileInfo8Bit(mFC);			
					
			//apply segmentation
			outTiff = jOD.findSoilSurface(nowTiff, mSFR);			
						
			//analyze surfaces
			//mSS[i] = jOD.extractSurfaceStatistics(myOutPath, mFC.myTiffs[i], outTiff, mFC.myInnerCircleFiles[myGandS[0]], nowTiff.getNSlices());
			
			//save file to 
			dT.displayColumnSurfacesByZ(nowTiff, mFC, outTiff);				
			
			//save result
			jIO.tiffSaver(myOutPath, mFC.myTiffs[i], outTiff);			
		}
		
		//save statistics
		jIO.saveSurfaceStatistics(myOutPath, mFC.myTiffs, mSS);
		
	}
}