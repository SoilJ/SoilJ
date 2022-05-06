package SoilJ_;

import java.io.File;
import java.util.ArrayList;

import SoilJ.tools.DisplayThings;
import SoilJ.tools.ImageManipulator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
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

import ij.IJ;
import ij.ImagePlus;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.plugin.Slicer;

/** 
 * FindColumnOutLines is a SoilJ plugin for the automatized soil column outline detection
 * 
 * @author John Koestel
 *
 */

public class FindColumnOutlines_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ObjectDetector jOD = new ObjectDetector();
		MenuWaiter menu = new MenuWaiter();
		DisplayThings disp = new DisplayThings();
		ImageManipulator jIM = new ImageManipulator();
				
		// init variables
		int i;
					
		//get details on PVC column finding
		MenuWaiter.ColumnFinderMenuReturn jCFS = menu.new ColumnFinderMenuReturn(); 
		jCFS = menu.showColumnFinderDialog();
		if (jCFS == null) return;
		
		//read file or files
		InputOutput.MyFileCollection mFC = jIO.fileSelector("Please choose a file or folder with your image data");
	 		
		//create output folder 
		String myOutFolderName = "WallCoordinateIdentified";
		String myOutFolder = mFC.myBaseFolder + pathSep + myOutFolderName;
		new File(myOutFolder).mkdir();
		mFC.myOutFolder = myOutFolder;
		
		//create innerCircle Folder
		String myInnerCircleFolder = myOutFolder + pathSep + "InnerCircle";
		new File(myInnerCircleFolder).mkdir();
		mFC.myInnerCircleFolder = myInnerCircleFolder;
		
		//remember failed column detections
		ArrayList<String> failedColumnDetection = new ArrayList<String>();
		
		//loop over 3D images
		int errorCounts = 0;
		for (i = 0 ; i < mFC.myTiffs.length ; i++) {  				
	
			//assign tiff file
			mFC.fileName = mFC.myTiffs[i];
			mFC = jIO.addCurrentFileInfo(mFC);
			mFC.startSlice = 1;
			mFC.stopSlice = mFC.nOfSlices + 1;
			if (!jCFS.try2FindColumnTopAndBottom) {
				mFC.startSlice = jCFS.topOfColumn;			
				if (jCFS.bottomOfColumn > 1) mFC.stopSlice = jCFS.bottomOfColumn;
			}
			mFC.pathSep = pathSep;
			mFC = jIO.addCurrentFileInfo(mFC);
			
			//check if everything is in order
			if (mFC.somethingIsWrong) {
				IJ.error(mFC.eMsg);
				return;
			}
			
			//reslice in case that the column is looked at from the wrong angle
			mFC.imageHasBeenLoaded = false;
			ImagePlus nowTiff = new ImagePlus();
			if (jCFS.doAResliceFirst) {
				
				Slicer sly = new Slicer();
				
				mFC.imageHasBeenLoaded = true;				
				
				nowTiff = jIO.openTiff3D(mFC);	
				
				Calibration cal = nowTiff.getCalibration();				
				mFC.caliZ = cal.pixelDepth;
				mFC.caliX = cal.pixelWidth;
				mFC.caliY = cal.pixelHeight;
				cal.disableDensityCalibration();
				cal.pixelDepth = 1;
				cal.pixelWidth = 1;
				cal.pixelHeight = 1;
				nowTiff.setCalibration(cal);
				
				ImagePlus reslicedTiff = sly.reslice(nowTiff);		
						
				mFC.nowTiff = reslicedTiff;
				
			}
			
			//try to find column outline but catch in case they were not found			
			try {
				
				//if column is not a steel column
				if (!jCFS.isSteel) {
		
					//load a stack of sample images			
					InputOutput.SampleTiffWrapper sTW = jIO.assembleRepresentativeSample(mFC);
									
					//re-determine the column outlines			
					boolean look4PreciseCoords = true;
					ObjectDetector.ColCoords3D prelimCC = jOD.findOrientationOfPVCOrAluColumn(sTW.samTiff, jCFS, look4PreciseCoords);
					
					//find some outlines seriously
					ObjectDetector.ColCoords3D samCoords = jOD.findColumnWalls3D(sTW.samTiff, prelimCC, jCFS, sTW.samSlices);
					
					//check if column stands straight and if it should be rotated upright into the center of the canvas
					if (jCFS.putColumnStraight) {
						
						if (!mFC.imageHasBeenLoaded) {
							mFC.nowTiff = jIO.openTiff3D(mFC.nowTiffPath);
							mFC.imageHasBeenLoaded = true;
						}
		
						mFC.nowTiff = jIM.putColumnUprightInCenter(mFC.nowTiff, samCoords, jCFS);				
						
						//store new image Facts
						mFC.nowWidth = mFC.nowTiff.getWidth();
						mFC.nowHeight = mFC.nowTiff.getHeight();
						mFC.nOfSlices = mFC.nowTiff.getNSlices();				
						
						//find column outlines once more...	
						sTW = jIO.assembleRepresentativeSample(mFC);				
						prelimCC = jOD.findOrientationOfPVCOrAluColumn(sTW.samTiff, jCFS, look4PreciseCoords);
						samCoords = jOD.findColumnWalls3D(sTW.samTiff, prelimCC, jCFS, sTW.samSlices);
					}
					
					//find upper end of column				
					ObjectDetector.ColCoords3D topCoords = null;
					if (jCFS.try2FindColumnTopAndBottom) topCoords = jOD.findColumnsTop(mFC, samCoords, jCFS);
					else topCoords = jOD.findClosestXYSlice2Top(mFC, samCoords, jCFS);
					
					//find lower end of column and create complete set of slice coordinates		
					ObjectDetector.ColCoords3D roughCoords = null;
					if (jCFS.try2FindColumnTopAndBottom) roughCoords = jOD.findColumnsBottom(mFC, topCoords, jCFS);
					else roughCoords = jOD.findClosestXYSlice2Bottom(mFC, topCoords, jCFS);				
			
					//find bevel
					ObjectDetector.ColCoords3D bevCoords = jOD.new ColCoords3D();
					if (jCFS.hasBevel) bevCoords = jOD.findColumnsBevel(mFC, roughCoords, jCFS);
					else bevCoords = roughCoords;
					
					//interpolate column outlines for non-samples slices
					ObjectDetector.ColCoords3D colCoords = jOD.imputeMissingLayers(bevCoords, jCFS);
					
					//define vertical extent of column that should be saved	
					int nowStart = colCoords.topOfColumn;
					int nowStop = colCoords.bottomOfColumn;
					colCoords.heightOfColumn = nowStop - nowStart;
					/* .. choice of cutting of top and bottom of 3-D canvas.. will be implemented later maybe..
					 * nowStart = (int)Math.round((double)nowStart - jCFS.topCutOff / 100 *
					 * (double)colCoords.heightOfColumn); nowStop = (int)Math.round((double)nowStop
					 * + jCFS.botCutOff / 100 * (double)colCoords.heightOfColumn);
					 * colCoords.topCutOff = nowStart; colCoords.botCutOff = nowStop;
					 */
					if (nowStart < 0) nowStart = 0;
					if (nowStop > mFC.nOfSlices) nowStop = mFC.nOfSlices;
					int[] loadSlices = new int[nowStop - nowStart];				
					for (int j = nowStart ; j < nowStop ; j++) loadSlices[j - nowStart] = j;
					
					//load image and cut out vertical ROI...
					if (!mFC.imageHasBeenLoaded) mFC.nowTiff = jIO.openTiff3D(mFC);						
					ImagePlus myTiff = jIM.selectSlicesFromTiff(mFC.nowTiff, loadSlices);
					
					//create and save Z-projections depicting the outlines of the detected wall
					disp.displayColumnOutlinesByZ(mFC.nowTiff, colCoords, mFC, null);				
					
					//save inner circle file
					String myInnerCirclePath = mFC.myInnerCircleFolder + pathSep + "Gauge" + mFC.colName + ".txt";
					jIO.writeInnerCircleVer1(myInnerCirclePath, colCoords);	
					
					//apply calibration function again
					Calibration cal = myTiff.getCalibration();
					cal.pixelDepth = mFC.caliZ;
					cal.pixelWidth = mFC.caliX;
					cal.pixelHeight = mFC.caliY;
					myTiff.setCalibration(cal);
					
					//and save... oufff!!
					jIO.tiffSaver(mFC.myOutFolder, mFC.fileName, myTiff);
					
					//free memory
					myTiff.flush();
					IJ.freeMemory();IJ.freeMemory();		
				}
				
				else {			// if it is a steel column				
					
					if (!mFC.imageHasBeenLoaded) {
						nowTiff = jIO.openTiff3D(mFC);
					}
					ObjectDetector.ColCoords3D prelimCC = jOD.findOrientationOfPVCOrAluColumn(nowTiff, jCFS, true);
					ObjectDetector.EggShapedColCoords3D eggshapedCC = jOD.findEggShapedWalls(nowTiff, prelimCC, jCFS.wallThickness);
					
					//create and save Z-projections depicting the outlines of the detected wall
					disp.displayColumnOutlinesByZ(nowTiff, eggshapedCC, mFC, null);
					
					//cut off top and bottom of image...
					ImagePlus outTiff = jIM.cutOffTopAndBottomSteel(nowTiff, eggshapedCC.topOfColumn, eggshapedCC.bottomOfColumn); 
					
					//and save... oufff!!
					jIO.tiffSaver(mFC.myOutFolder, mFC.fileName, outTiff);
					
					//save column coordinates and cut image...
					jIO.writeInnerCircleSteel(eggshapedCC, mFC);
					
				}
			}
			catch(Exception e) {
				failedColumnDetection.add(mFC.myTiffs[i]);
			}
		}
		
		if (failedColumnDetection.size() == 0) {
			IJ.log("It appears that the column detection went smoothly for all samples.\n "
					+ "Please check the column detection performance in WallCoordinateidentified/ReviewFoundOutlines "
					+ "if you are happy with the found column outlines.");			
		}
		else {
			String eMsg = "No column could be detected for the following samples:\n\n";
			for (i = 0 ; i < failedColumnDetection.size() ; i++) eMsg += failedColumnDetection.get(i) + "\n";
			eMsg += "\n Please re-run the column detection for these columns using a lower R2 value "
					+ "or using less strict criteria for a successful column find. ";					
		}
		
		if (errorCounts > 0) {
			
			String eMsg = "\n\nThis plugin requires images on which the entire outer perimeter of the column wall\n";		
			eMsg += "is visible. If the wall is only partially visible, the plugin will probably not work.\n\n";		
			eMsg += "In the case that the contrast between wall and canvas gray-values is very small, the plugin is\n";
			eMsg += "doomed to fail.\n\n";
			eMsg += "In the case that the image is very noisy, use a filter to reduce the image noise\n";
			eMsg += "Suitable filters are a non-local means filter or, alternatively, a median filter \n";			
			eMsg += "followed by an unsharp mask.\n\n";
			if (jCFS.isAlu) {
				eMsg += "In the case that the contrast between column wall and soil is not very strong, \n";
				eMsg += "consider hacking in the wall thickness manually\n.";
				if(jCFS.hasBevel) {
					eMsg += "The location of the bevel should be found nevertheless.\n\n";
				}
				else eMsg += "\n";		
			}
			eMsg += "Thank you for using SoilJ!";			
			
			IJ.log(eMsg);
		}		
		
		IJ.freeMemory();IJ.freeMemory();
	}
		
}

