package SoilJ_;

import java.io.File;

import org.apache.commons.math3.stat.StatUtils;

import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.MorphologyAnalyzer;
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

import ij.ImagePlus;
import ij.plugin.PlugIn;

/** 
 * ExtractPoresizeDistributio is a SoilJ plugin that extracts the pore size distribution from an
 * image depicting pore diameters (thicknesses) and writes it into an ASCII file.
 * 
 * @author John Koestel
 *
 */

public class Extract32BitHistogram_ extends ImagePlus implements PlugIn  {

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
		MenuWaiter menu = new MenuWaiter();
		RoiHandler roi = new RoiHandler();
		MorphologyAnalyzer morph = new MorphologyAnalyzer();
		
		//select ROI
		MenuWaiter.ROISelectionOptions mRSO = menu.regionOfInterestSelection();
		
		//choose files
		InputOutput.MyFileCollection mFC = jIO.fileSelector("Please choose a file or folder with your image data");

		mFC.myOutFolder = mFC.myBaseFolder + pathSep + "Histograms"; 
		new File(mFC.myOutFolder).mkdir();
		
		//loop over 3D images
		for (int i = 0 ; i < mFC.myTiffs.length ; i++) {  //myTiffs.length

			//assemble image for analyses
			mFC.fileName = mFC.myTiffs[i];
			mFC = jIO.addCurrentFileInfo(mFC);
			
			//load file			
			int[] startStopSlices = jIO.findStartAndStopSlices(mFC, mRSO);
			int[] colSlices = new int[startStopSlices[1]  - startStopSlices[0]];
			for (int j = 0 ; j < colSlices.length ; j++) colSlices[j] = startStopSlices[0] + j;
		
			//cut roi		
			mFC.startSlice = startStopSlices[0];
			mFC.stopSlice = startStopSlices[1];
			RoiHandler.ColumnRoi colRoi = roi.prepareDesiredRoi(mFC, jIO.openTiff3DSomeSlices(mFC, colSlices), mRSO, "null", "all");
			
			//load image
			String myName = mFC.fileName.substring(0,mFC.fileName.length()-4) + ".hist";
		
			//define class bounds.. needed due to sphere fitting artefacts..
			double[] maximum = {colRoi.nowTiff.getWidth(), colRoi.nowTiff.getHeight()};
			double[] classBounds = new double[(int)StatUtils.max(maximum)];
			for (int j = 0 ; j < classBounds.length ; j++) classBounds[j] = j;				
			
			//apply segmentation							}
			morph.extract32BitHisto(mFC.myOutFolder + pathSep + myName, colRoi.nowTiff, classBounds);			
			
		}
		
	}
}