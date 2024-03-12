package SoilJ_;

import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.MorphologyAnalyzer;
import SoilJ.tools.RoiHandler;
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
 * ExtractLoopDiameterDistribution is a SoilJ plugin for extracting the diameter of loops (redundant connections) in 
 * the analyzed image phase
 * 
 * @author John Koestel
 *
 */

public class ExtractLoopDiameterDistribution_ extends ImagePlus implements PlugIn  {

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
		
		// init variables
		int i;
		
		//tell me what I should do!
		MenuWaiter.PoreSpaceAnalyzerOptions mLDD = menu.showLoopDiameterDistMenu();
		if (mLDD== null) return;
		
		//create Folder structure
		InputOutput.MyFileCollection mFC = jIO.createFolders4SubROIData(mLDD);
					
		//loop over 3D images
		for (i = 0 ; i < mFC.myTiffs.length ; i++) {  //myTiffs.length
			
			//assemble image for analyses
			mFC.fileName = mFC.myTiffs[i];
			mFC = jIO.addCurrentFileInfo8Bit(mFC);
			
			//load file			
			int[] startStopSlices = jIO.findStartAndStopSlices(mFC, mLDD.mRSO);
			int[] colSlices = new int[startStopSlices[1] - startStopSlices[0]];
			for (int j = 0 ; j < colSlices.length ; j++) colSlices[j] = startStopSlices[0] + j;
		
			//cut roi		
			mFC.startSlice = startStopSlices[0];
			mFC.stopSlice = startStopSlices[1];
			RoiHandler.ColumnRoi colRoi = roi.prepareDesiredRoi(mFC, jIO.openTiff3DSomeSlices(mFC, colSlices), mLDD.mRSO, mLDD.imagePhase2BeAnalyzed, mLDD.nameOfAnalyzedPhase);
			
			//apply analyzes
			morph.getLoopDiameterDistribution(i, mFC, colRoi, mLDD);		
			
			//try to free up some memory
			colRoi.nowTiff.unlock();			
			colRoi.nowTiff.flush();

			System.gc();System.gc();
			IJ.showStatus("Trying to clean up memory ... ");
			//IJ.wait(mLDD.includeBreaks);
			IJ.wait(1000);
			IJ.freeMemory();IJ.freeMemory();
			System.gc();System.gc();
			
			IJ.showStatus("Pore space analyses has been completed. ");
			
		}		
	}
			
	
}