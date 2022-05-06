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
 * NormalizeTheseImages is a SoilJ plugin that standardizes the greyscale of an image with respect to 
 * grey value quantiles of horizontal slices and/or the grey value of the columns wall.
 * 
 * @author John Koestel
 *
 */

public class CalibrateGrayValuesManually_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {
		
		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";
		
		
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = CalibrateGrayValues_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ImageManipulator jIM = new ImageManipulator();	
		MenuWaiter menu = new MenuWaiter();

		//select file or files
	  	InputOutput.MyFileCollection mFC = jIO.fileSelector("Please choose a file or folder with your image data");
			
		//set folder structure
		String myOutFolder;
		String myOutPath = null;		
		myOutFolder = "ManualStandard_5000_20000";		
		myOutFolder.replace(',', '.');
		myOutPath = mFC.myBaseFolder + pathSep + myOutFolder;
		new File(myOutPath).mkdir();
		mFC.myOutFolder = myOutPath;
		
		//loop over 3D images
		for (int i = 0 ; i < mFC.myTiffs.length ; i++) {  //myTiffs.length
		
			//assemble image for analyses
			mFC.fileName = mFC.myTiffs[i];
			jIO.addCurrentFileInfo(mFC);
	
			InputOutput.CaliRoiLocations myCRL = jIO.readGrayValueCalibrationFile(mFC.myBaseFolder + pathSep + mFC.colName + ".txt");
						
			ImagePlus nowTiff = jIO.openTiff3D(mFC.nowTiffPath);
			
			ImagePlus outTiff = jIM.calibrateGrayValuesFromList(nowTiff, myCRL);
		
			jIO.tiffSaver(myOutPath, mFC.fileName, outTiff);
	
		}
	
	}
	
}