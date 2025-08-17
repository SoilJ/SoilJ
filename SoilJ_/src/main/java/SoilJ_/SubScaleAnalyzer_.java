package SoilJ_;

import java.io.File;

import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.MorphologyAnalyzer;
import SoilJ.tools.REVAnalyses;
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
 * REVAnalyzer is a SoilJ plugin for investigating the existence and size of REV for soil 
 * pore network morphological properties
 * 
 * @author John Koestel
 *
 */

public class SubScaleAnalyzer_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		String pathSep = "\\";
		
		//probe operating system and adjust pathSep if necessary
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
		
		//shall I cut away something?
		MenuWaiter.REVAnalyzerOptions mRA = menu.showREVAnalyzerMenu();
		if (mRA == null) return;
		
		if (mRA != null) {						//check whether mRA has returned reasonable values
		
			//construct image related objects
			ImagePlus nowTiff = new ImagePlus();
			
			//select file or files
		  	InputOutput.MyFileCollection mFC = jIO.fileSelector("Please choose a file or folder with your image data");
		  						
			//create output paths
			String myPreOutFolder = "";
			if (mRA.choiceOfRoi.equals("Cuboid")) {
				myPreOutFolder += "Cub_XL" + mRA.cubeX1 + "XR" + mRA.cubeX2 + "YL" + mRA.cubeY1 + "YR" + mRA.cubeY2 + "ZT" + mRA.cubeZ1 + "ZB" + mRA.cubeZ2;
			}
			//if (mRA.choiceOfRoi.equals("Cylinder")) myPreOutFolder = "REVCyl_X" + mRA.cylX + "Y" + mRA.cylY+ "ZT" + mRA.cylZ1 + "ZB" + mRA.cylZ2 + "R" + mRA.cylRadius;	
			if (mRA.choiceOfMethod.equals("Sub-ROIs by division")) {
				myPreOutFolder = "Div" + myPreOutFolder;
			}
			if (mRA.choiceOfMethod.equals("Sub-ROIs by shrinkage")) {
				myPreOutFolder = "Shr" + myPreOutFolder;
			}
	/*		if (mRA.choiceOfMethod.equals("Sub-ROI by division and shrinkage")) {
				myPreOutFolder = "DiSh" + myPreOutFolder;
			}*/
			if (mRA.choiceOfMethod.equals("Moving sub-ROIs")) {
				myPreOutFolder = "Mov_EX" + mRA.moveEdgeX + "EY" + mRA.moveEdgeY + "EZ" + mRA.moveEdgeZ + "in" + myPreOutFolder;
			}
			
			myPreOutFolder = "REV" + myPreOutFolder;
					
			String myPreOutPath = mFC.myBaseFolder + myPreOutFolder;
			new File(myPreOutPath).mkdir();		
			mFC.myPreOutFolder = myPreOutPath; 	//also add it to the folder collection
			
			//loop over 3D images			
			for (i = 0 ; i < mFC.myTiffs.length ; i++) {  //myTiffs.length
		
				//try to free up some memory
				System.gc();
				
				//assign tiff file
				mFC.fileName = mFC.myTiffs[i];
				mFC = jIO.addCurrentFileInfo(mFC);
				
				//create a new folder for saving the results
				mFC.myPreOutFolder = myPreOutPath + pathSep + mFC.myTiffs[i].substring(0, mFC.myTiffs[i].length() - 4);
				new File(mFC.myPreOutFolder).mkdir();
				
				//load image
				nowTiff = jIO.openTiff3D(mFC.myBaseFolder + pathSep + mFC.myTiffs[i]);	
				
				//prepare coordinates for each sub-ROI
				REVAnalyses reAn = new REVAnalyses();
				REVAnalyses.REVAnalysesPack revPack = reAn.createROIs4REVAnalyses(mRA);
				
				//plug in REV option into poresSpace Analyzer options
				MenuWaiter.PoreSpaceAnalyzerOptions mPSA = reAn.castRA2PSA(mRA);
				
				mFC.myTiffs = revPack.roiName;
				
				//loop over all ROIs
				for (int j = 0 ; j < revPack.numberOfROIs ; j++) {
					
					//prepare ROI				
					mPSA.mRSO.cubeX1 = revPack.x1[j];mPSA.mRSO.cubeX2 = revPack.x2[j];
					mPSA.mRSO.cubeY1 = revPack.y1[j];mPSA.mRSO.cubeY2 = revPack.y2[j];
					mPSA.mRSO.cubeZ1 = revPack.z1[j];mPSA.mRSO.cubeZ2 = revPack.z2[j];
					mPSA.mRSO.choiceOfRoi = mRA.choiceOfRoi;
					
					mFC.colName = revPack.roiName[j];
					RoiHandler.ColumnRoi colRoi = roi.prepareDesiredRoi(mFC, nowTiff, mPSA.mRSO, mPSA.imagePhase2BeAnalyzed, mPSA.nameOfAnalyzedPhase);
													
					//apply analyzes
					morph.aTailoredPoreSpaceAnalyzer(j, mFC, colRoi, mPSA);
				}
			}
		}
	}
}