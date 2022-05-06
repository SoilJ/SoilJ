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
 * BeamDeHardening is a SoilJ plugin aiming at the removal of concentric imaging artifacts
 * as they are occurring due to beam hardening. 
 * 
 * @author John Koestel
 *
 */

public class BeamDeHardening_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {
		
		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";

		/*
		 * // set the plugins.dir property to make the plugin appear in the Plugins menu
		 * Class<?> clazz = BeamDeHardening_.class; String url = clazz.getResource("/" +
		 * clazz.getName().replace('.', '/') + ".class").toString(); String pluginsDir =
		 * url.substring(5, url.length() - clazz.getName().length() - 6);
		 * System.setProperty("plugins.dir", pluginsDir);
		 */
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ImageManipulator jIM = new ImageManipulator();
		MenuWaiter menu = new MenuWaiter();
						
		// init variables
		int i;
		
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();
		ImagePlus outTiff = new ImagePlus();
		
		//ask for threshold choice
	    MenuWaiter.BeamDeHardeningReturn mBDH = menu.showBeamDeHardeningMenu();
	    if (mBDH == null) return;

	    //read file or files
	  	InputOutput.MyFileCollection mFC = jIO.fileSelector("Please choose a file or folder with your image data");
	  	
	  	//also have a look for innerCircle folder
	  	mFC = jIO.addInnerCircleFolder(mFC);
	  	
	  	//save and create output folder		
		String myOutFolder;
		String myOutPath = null;
		myOutFolder = "BeamDeHardened";			
		myOutPath = mFC.myBaseFolder + pathSep + myOutFolder;
		new File(myOutPath).mkdir();		
				
		//loop over 3D images
		for (i = 0 ; i < mFC.myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();		
			
			//echelon info on current file
			mFC.fileName = mFC.myTiffs[i];
			mFC = jIO.addCurrentFileInfo(mFC);
			
			//select the correct gauge and surface files			
			String[] myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC);
			mFC.nowInnerCirclePath = mFC.myInnerCircleFolder + pathSep + myGandS[0];
			
			//load image
			nowTiff = jIO.openTiff3D(mFC.myBaseFolder + pathSep + mFC.myTiffs[i]);	
						
			outTiff = jIM.beamDeHardenThis(nowTiff, mFC, mBDH);
			
			//save result
			jIO.tiffSaver(myOutPath, mFC.myTiffs[i], outTiff);
		}
		
	}
	
}