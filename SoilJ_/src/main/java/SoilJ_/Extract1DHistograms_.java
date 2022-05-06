package SoilJ_;

import java.io.File;

import SoilJ.tools.HistogramStuff;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.ObjectDetector;
import SoilJ.tools.RoiHandler;
import SoilJ.tools.RollerCaster;
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
 * ExtractPoresizeDistributio is a SoilJ plugin that extracts the pore size distribution from an
 * image depicting pore diameters (thicknesses) and writes it into an ASCII file.
 * 
 * @author John Koestel
 *
 */

public class Extract1DHistograms_ extends ImagePlus implements PlugIn  {

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
		ObjectDetector jOD = new ObjectDetector();
		MenuWaiter menu = new MenuWaiter();
		RoiHandler roi = new RoiHandler();
		HistogramStuff hist = new HistogramStuff();
		RollerCaster rC = new RollerCaster();
		
		// init variables
		int i;
		
		//tell me what I should do!
		MenuWaiter.ROISelectionOptions mRSO = menu.showHistogramExtractionMenu();
		if (mRSO == null) return;
		
		//create Folder structure
		InputOutput.MyFileCollection mFC = jIO.createFolders4SubROIData(mRSO);
		
		//also create a folder for the cut surfaces
		if (mRSO.includeSurfaceTopography) {
			String myCutSurfacePath = mFC.myPreOutFolder + pathSep + "CutSurfaceFiles";
			new File(myCutSurfacePath).mkdir();
			mFC.myCutSurfaceFolder = myCutSurfacePath;
		}
		
		//probe type of images that are going to be processed..
		mFC.fileName = mFC.myTiffs[0];
		mFC = jIO.addCurrentFileInfo(mFC);
		
		float[][] allHists = null;
		if (mFC.bitDepth == 16) allHists = new float[mFC.myTiffs.length][(int)Math.round(Math.pow(2, 16))];
		else allHists = new float[mFC.myTiffs.length][(int)Math.round(Math.pow(2, 8))];
		
		//loop over 3D images
		for (i = 0 ; i < mFC.myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();
	
			//write imageProperties in folder collection
			mFC.fileName = mFC.myTiffs[i];
			if (mFC.bitDepth == 16) mFC = jIO.addCurrentFileInfo(mFC);
			else mFC = jIO.addCurrentFileInfo8Bit(mFC);
			
			//load file                                                                                                               
			int[] startStopSlices = jIO.findStartAndStopSlices(mFC, mRSO);			
			int[] colSlices = new int[startStopSlices[1]  - startStopSlices[0]];
			for (int j = 0 ; j < colSlices.length ; j++) colSlices[j] = startStopSlices[0] + j;
			ImagePlus nowTiff = jIO.openTiff3DSomeSlices(mFC, colSlices);
			
			//cut roi		
			mFC.startSlice = startStopSlices[0];
			mFC.stopSlice = startStopSlices[1];
			
			//cut image... the ROI handler needs the PoreSpaceAnalyzerOptions as input.. should be changed in the future..			
			MenuWaiter.PoreSpaceAnalyzerOptions mPSA = menu.new PoreSpaceAnalyzerOptions(); 
			mPSA.mRSO = mRSO;
			mPSA.imagePhase2BeAnalyzed = "null";
			RoiHandler.ColumnRoi colRoi = roi.prepareDesiredRoi(mFC, nowTiff, mRSO, mPSA.imagePhase2BeAnalyzed, "");					
						
			//try to free up some memory			
			IJ.freeMemory();IJ.freeMemory();		
			
			//apply segmentation	
			int[] histo = new int[(int)Math.pow(2, mFC.bitDepth)];
			if (mFC.bitDepth == 16) histo = hist.extractHistograms16(mFC, colRoi.nowTiff);
			else histo = hist.extractHistograms8(mFC, colRoi.nowTiff);

			//save it
			jIO.writeHistogram(mFC, histo);				
		}
		
	}
}