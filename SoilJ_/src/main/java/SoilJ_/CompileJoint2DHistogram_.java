package SoilJ_;

import java.io.File;

import SoilJ.tools.HistogramStuff;
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

public class CompileJoint2DHistogram_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";
		
		/*
		 * //set the plugins.dir property to make the plugin appear in the Plugins menu
		 * Class<?> clazz = Extract2DHistograms_.class; String url =
		 * clazz.getResource("/" + clazz.getName().replace('.', '/') +
		 * ".class").toString(); String pluginsDir = url.substring(5, url.length() -
		 * clazz.getName().length() - 6); System.setProperty("plugins.dir", pluginsDir);
		 */
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();				
		HistogramStuff hist = new HistogramStuff();		

		//read file or files
		InputOutput.MyFileCollection mFC = jIO.fileSelector("Please choose the folder with your 2D histogram data");
    	    
	    //set output folder  		 		
  		String myOutPath = mFC.myBaseFolder + pathSep + "Statistics";
  		new File(myOutPath).mkdir();	
  		mFC.myOutFolder = myOutPath;
  		
  		//save pathSep
  		mFC.pathSep = pathSep;
		
		//loop over 3D images
		hist.analyze2DHistograms(mFC);
		
		//try to clear memory
		IJ.freeMemory();IJ.freeMemory();
		
		IJ.showStatus("The 2D histograms have now been analyzed!");
		
	}
}