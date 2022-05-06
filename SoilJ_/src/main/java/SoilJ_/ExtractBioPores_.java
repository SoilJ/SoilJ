package SoilJ_;

import java.io.File;

import SoilJ.tools.ImageManipulator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
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

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;

/** 
 * PoreSpaceAnalyzer is a SoilJ plugin that calculates several morphologic properties from 
 * binary images. To a large part it makes use of plugins collected in the BoneJ library:
 * Doube M, Kłosowski MM, Arganda-Carreras I, Cordeliéres F, Dougherty RP, Jackson J, Schmid B, Hutchinson JR, Shefelbine SJ. (2010) BoneJ: free and extensible bone image analysis in ImageJ. Bone 47:1076-9. doi: 10.1016/j.bone.2010.08.023 
 * www.bonej.org
 * 
 * The results of the analyses are written into image and ASCII files. 
 * 
 * @author John Koestel
 *
 */

public class ExtractBioPores_ extends ImagePlus implements PlugIn  {

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
		ImageManipulator jIM = new ImageManipulator(); 
		
		// init variables
		int i;
		
		//tell me what I should do!
		MenuWaiter.ROISelectionOptions mRSO = menu.regionOfInterestSelection();
		if (mRSO == null) return;
		
		//get some bio pore extraction parameters
		MenuWaiter.BioPoreExtractionOptions mBEO = menu.showBioPoreExtractionMenu();
		
		//create Folder structure
		InputOutput.MyFileCollection mFC = jIO.createFolders4SubROIData(mRSO);
		
		String myOutFolder = mFC.myOutFolder;
		mFC.myOutFolder = myOutFolder + pathSep + "BioPores" + (int)mBEO.numberOfSigmas + "Sigmas"; 	 	
		new File(mFC.myOutFolder).mkdir();
		mFC.myOutFolder2 = myOutFolder + pathSep + "NonBioPores" + (int)mBEO.numberOfSigmas + "Sigmas"; 	 	
		new File(mFC.myOutFolder2).mkdir();
					
		//loop over 3D images
		for (i = 0 ; i < mFC.myTiffs.length ; i++) {  //myTiffs.length
			
			//assemble image for analyses
			mFC.fileName = mFC.myTiffs[i];
			mFC = jIO.addCurrentFileInfo8Bit(mFC);
			
			//load file			
			int[] startStopSlices = jIO.findStartAndStopSlices(mFC, mRSO);
			int[] colSlices = new int[startStopSlices[1]  - startStopSlices[0]];
			for (int j = 0 ; j < colSlices.length ; j++) colSlices[j] = startStopSlices[0] + j;
		
			//cut roi		
			mFC.startSlice = startStopSlices[0];
			mFC.stopSlice = startStopSlices[1];
			RoiHandler.ColumnRoi colRoi = roi.prepareDesiredRoi(mFC, jIO.openTiff3DSomeSlices(mFC, colSlices), mRSO, "255", "");
			
			System.gc();System.gc();
			IJ.showStatus("Trying to clean up memory ... ");
			IJ.freeMemory();			
			IJ.wait(1000);
			IJ.freeMemory();IJ.freeMemory();
			System.gc();System.gc();			
			
			//apply analyzes			
			ImagePlus outTiff = jIM.extractBioPores(mFC, colRoi.nowTiff, mBEO);		
						
			//save result
			jIO.tiffSaver(mFC, outTiff);
			
			//also save non-bio pores
			jIO.tiffSaver2(mFC, jIM.extractNonBioPores(mFC, colRoi.nowTiff, outTiff, mBEO));	
			
			//try to free up some memory
			colRoi.nowTiff.unlock();			
			colRoi.nowTiff.flush();
			outTiff.unlock();
			outTiff.flush();

			System.gc();System.gc();
			IJ.showStatus("Trying to clean up memory ... ");
			IJ.freeMemory();IJ.freeMemory();
			System.gc();System.gc();
			
			IJ.showStatus("Biopores have been extracted.");
			
		}		
	}
}