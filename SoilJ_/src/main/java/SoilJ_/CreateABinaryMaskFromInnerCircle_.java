package SoilJ_;

import org.apache.commons.math3.stat.StatUtils;

import SoilJ.tools.ImageManipulator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.ObjectDetector;
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
import ij.gui.PolygonRoi;
import ij.plugin.PlugIn;

/** 
 * CreateABinaryMaskFromInnerCircle a SoilJ plugin that creates a binary mask from InnerCircle coordinates for the use in elastix 
 * 
 * The results of the analyses are written into image and ASCII files. 
 * 
 * @author John Koestel
 *
 */

public class CreateABinaryMaskFromInnerCircle_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";
		
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = CreateABinaryMaskFromInnerCircle_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();	
		MenuWaiter menu = new MenuWaiter();
		RoiHandler roi = new RoiHandler();
		ImageManipulator jIM = new ImageManipulator();	
		ObjectDetector jOD = new ObjectDetector();
		
		// init variables
		int i;
		
		//tell me what I should do!
		MenuWaiter.ROISelectionOptions mRSO = menu.regionOfInterestSelection();
		if (mRSO == null) return;
		
		//create Folder structure
		InputOutput.MyFileCollection mFC = jIO.createFolders4SubROIData(mRSO);		
					
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

			String[] GandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC);
			mFC.nowInnerCirclePath = mFC.myInnerCircleFolder + mFC.pathSep + GandS[0];
			
			//read InnerCircle file and create outline ROI
			if (GandS[0].contains("Gauge")) {
				
				ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();				
				int versio = jIO.checkInnerCircleFileVersion(mFC.nowInnerCirclePath);			
				if (versio == 0) jCO = jIO.readInnerCircleVer0(mFC.nowInnerCirclePath);	
				else jCO = jIO.readInnerCircleVer1(mFC.nowInnerCirclePath);	
				

				//calculate average area
				double[] avgRadius = new double[jCO.innerMajorRadius.length];
				for (int radiusFinder = 0 ; radiusFinder < jCO.innerMajorRadius.length ; radiusFinder++) {
					avgRadius[radiusFinder] = (jCO.innerMajorRadius[radiusFinder] + jCO.innerMinorRadius[radiusFinder]) / 2;
				}
				double myRadius = StatUtils.mean(avgRadius) - mRSO.cutAwayFromWall;
				
				int averageRadius = (int)Math.round(myRadius);
				PolygonRoi[] pRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, -mRSO.cutAwayFromWall);
				PolygonRoi[] iRoi = null;
				if (mRSO.cutAwayFromWall == 0) iRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 10);
				if (mRSO.cutAwayFromCenter > 0) iRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, -averageRadius);
				
				//create and save mask
				jIM.createROIFromInnerCircle(mFC, pRoi);
			}			
		}
		
		//try to free up some memory			
		IJ.freeMemory();IJ.freeMemory();	
		
	}
}