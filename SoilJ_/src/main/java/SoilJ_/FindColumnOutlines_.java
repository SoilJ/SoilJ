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

		String pathSep = System.getProperty("file.separator");
//		String pathSep = "\\";
//		//probe operating system and adjust PathSep if necessary
//		String myOS = System.getProperty("os.name");
//		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ObjectDetector jOD = new ObjectDetector();		// contains ColCoords3D, but also illumination, find orientation, etc. 
		MenuWaiter menu = new MenuWaiter();
		DisplayThings disp = new DisplayThings();
		ImageManipulator jIM = new ImageManipulator();
					
		//get details on PVC column finding
		MenuWaiter.ColumnFinderMenuReturn jCFS = menu.new ColumnFinderMenuReturn(); 	// get "user settings" selected in the menu
		jCFS = menu.showColumnFinderDialog();											// menu where user chooses whether to find top & bottom, whether z-axis, good or bad contrast, steel cylinder, etc.
		if (jCFS == null) return;
		
		//remember virgin jCFS values
	    double virCVWallBrightnessThresh = jCFS.CVWallBrightnessThresh;
	    int virFixedWallGrayValue = jCFS.fixedWallGrayValue;
	    boolean virGrayValueIsKnown = jCFS.grayValueOfWallIsKnown;
	    boolean virIsAlreadyNormalized = jCFS.isAlreadyNormalized;
	    int virMaxColGV = jCFS.maxColGV;
	    int virMinColGV = jCFS.minColGV;
	    int stdFixedWallGrayValue = jCFS.stdFixedWallGrayValue;
	    
		//read file or files
		InputOutput.MyFileCollection mFC = jIO.fileSelector("Please choose a file or folder with your image data"); 	// chose the images for which to find outline
	 		
		//create output folder 
		String myOutFolderName = "WallCoordinateIdentified";
		String myOutFolder = mFC.myBaseFolder + pathSep + myOutFolderName;
		new File(myOutFolder).mkdir();	// create the folder
		mFC.myOutFolder = myOutFolder;	// store the information in the mFC (mFC is kind of info storage for everything about current image and total of images
		
		//create innerCircle Folder
		String myInnerCircleFolder = myOutFolder + pathSep + "InnerCircle";
		new File(myInnerCircleFolder).mkdir();			// create the folder
		mFC.myInnerCircleFolder = myInnerCircleFolder;	// store name in mFC info collection
		
		//remember failed column detections
		ArrayList<String> failedColumnDetection = new ArrayList<String>();	// create array for storing "failures"
		
		//loop over 3D images >> look for column outlines in each image separately
		int errorCounts = 0;
		for (int i = 0 ; i < mFC.myTiffs.length ; i++) {
	
			//try to purge unnecessary memory
			System.gc();
			IJ.wait(200);
			System.gc();
			
			//restart with the virgin jCFS file values
		    jCFS.CVWallBrightnessThresh = virCVWallBrightnessThresh;
		    jCFS.fixedWallGrayValue = virFixedWallGrayValue;
		    jCFS.grayValueOfWallIsKnown = virGrayValueIsKnown;
		    jCFS.isAlreadyNormalized = virIsAlreadyNormalized;
		    jCFS.maxColGV = virMaxColGV;
		    jCFS.minColGV = virMinColGV;
		    jCFS.stdFixedWallGrayValue = stdFixedWallGrayValue;
			
			//assign tiff file
			mFC.fileName = mFC.myTiffs[i];				// fileName in mFC conatins current information that is then overwritten for next image
			mFC = jIO.addCurrentFileInfo(mFC);	
			
			mFC.startSlice = 1;							// start at top
			mFC.stopSlice = mFC.nOfSlices;			   // to bottom
			if (!jCFS.try2FindColumnTopAndBottom) {		// promting settings (chosen by user in menu) >> should program try to find the top and the bottom ? only either currently not possible
				mFC.startSlice = jCFS.topOfColumn;			// if it should not look for bottom >> take top as given in menu
				if (jCFS.bottomOfColumn > 1) mFC.stopSlice = jCFS.bottomOfColumn; // if entered bottom > 1 >> take bottom as indicated by user in menu
			}
			mFC.pathSep = pathSep; // (why?)
			mFC = jIO.addCurrentFileInfo(mFC); 	//again???		
			
			//check if everything is in order
			if (mFC.columnHasNotBeenFound) {
				IJ.error(mFC.eMsg);
				return;
			}
			
			//save some feedback to log
			IJ.log("//////////////////////////////////////////////////////////////////////////////////////////////////////");
			IJ.log("Processing " + mFC.fileName + " ...");
			IJ.log("//////////////////////////////////////////////////////////////////////////////////////////////////////");
			IJ.log("");
			
			//reslice in case that the column is looked at from the wrong angle
			mFC.imageHasBeenLoaded = false;
			ImagePlus nowTiff = new ImagePlus();			// nowTiff is current tiff to look at
			if (jCFS.doAResliceFirst) {
				
				Slicer sly = new Slicer();					// slicer:  image is sampled at constant _pixel_ increments (distance 1)
				
				mFC.imageHasBeenLoaded = true;				
				nowTiff = jIO.openTiff3D(mFC);				// open the current image (in image loop)  >> probably uses the variable fileName in mFC which is current image
				
				Calibration cal = nowTiff.getCalibration();	// returns the imagecalibration (if e.g. mm voxel size)  >> store that info in mFC
				mFC.caliZ = cal.pixelDepth;
				mFC.caliX = cal.pixelWidth;
				mFC.caliY = cal.pixelHeight;
				cal.disableDensityCalibration();
				cal.pixelDepth = 1;							// set the calibration to 1 (voxel size is always 1 voxel) instead of anything else
				cal.pixelWidth = 1;
				cal.pixelHeight = 1;
				nowTiff.setCalibration(cal);
				
				ImagePlus reslicedTiff = sly.reslice(nowTiff);	// then reslice the current tiff
						
				mFC.nowTiff = reslicedTiff;					// overwrite the current tiff to analyze with the resliced version of itself
				
			}
			
			//try to find column outline but catch in case they were not found			
			//try {
				
				//if column is not a steel column
				if (!jCFS.isSteel) {					// menu selection whether steel cylinder
		
					//load a stack of sample images			
					InputOutput.SampleTiffWrapper sTW = jIO.assembleRepresentativeSample(mFC);		// returns class sampleTiffWrapper (class of public ImagePlus, int[] samTiffSliceNumbers, int[] samSlices, boolean hasConverged)
						
					//count fingding attempts
					int findingAttemptNumber = 0;
					
					//try all over again once the algorithm has hit the wall..
					while (findingAttemptNumber < 5) {
						
						//try to purge unneccessary memory
						System.gc();
						IJ.wait(200);
						System.gc();
						
						//count attempt
						findingAttemptNumber++;
						
						//break loop in case the attempts did not make things better. 
						if (findingAttemptNumber > 4) {
							if (findingAttemptNumber < 100) {
								IJ.log("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
								IJ.log("Column in " + mFC.fileName + " could definitely not be found.. \nskipping this sample.. please have a closer look at this sample");
								IJ.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
								failedColumnDetection.add(mFC.fileName);
							}
							break;
						}
						
						//change slices that are to be checked
						if (findingAttemptNumber > 1) {
							sTW = jIO.assembleRepresentativeSample(mFC);
						    //jCFS.CVWallBrightnessThresh = virCVWallBrightnessThresh;
						    //jCFS.fixedWallGrayValue = virFixedWallGrayValue;
						    //jCFS.grayValueOfWallIsKnown = virGrayValueIsKnown;
						    //jCFS.isAlreadyNormalized = virIsAlreadyNormalized;						    
						    //jCFS.maxColGV = virMaxColGV;
						    //jCFS.minColGV = virMinColGV;
						    //jCFS.stdFixedWallGrayValue = stdFixedWallGrayValue;
						}
						
						IJ.log("Attempt #" + findingAttemptNumber + " to find column for " + mFC.fileName + " ...");
												
						//re-determine the column outlines
						boolean look4PreciseCoords = true;
						ObjectDetector.ColCoords3D prelimCC = jOD.findOrientationOfPVCOrAluColumn(sTW.samTiff, jCFS, look4PreciseCoords);  // returns ColCoords3D object
						
						//find some outlines seriously
						ObjectDetector.ColCoords3D samCoords = jOD.findColumnWalls3D(sTW.samTiff, prelimCC, jCFS, sTW.samSlices);			// returns ColCoords3D object with the "definitive" ROIs
						mFC.columnHasNotBeenFound = false;						
						
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
						
						//try {
						
							//find upper end of column				
							ObjectDetector.ColCoords3D topCoords = null;
							if (jCFS.try2FindColumnTopAndBottom) {
								topCoords = jOD.findColumnsTop(mFC, samCoords, jCFS);
								if (topCoords == null) {
									IJ.log("Top of " + mFC.fileName + " could not be found.. let's try all over again\n");
									mFC.columnHasNotBeenFound = true;
								}
							}
							else topCoords = jOD.findClosestXYSlice2Top(mFC, samCoords, jCFS);					
			
							//find lower end of column and create complete set of slice coordinates							
							if (!mFC.columnHasNotBeenFound) {
								ObjectDetector.ColCoords3D roughCoords = null;						
								if (jCFS.try2FindColumnTopAndBottom) {						
									roughCoords = jOD.findColumnsBottom(mFC, topCoords, jCFS);				
									if (roughCoords == null) {
										IJ.log("Bottom of " + mFC.fileName + " could not be found.. let's try all over again\n");
										mFC.columnHasNotBeenFound = true;
									}						
								}
								else roughCoords = jOD.findClosestXYSlice2Bottom(mFC, topCoords, jCFS);
						
								//find bevel (abgeschrägt)						
								ObjectDetector.ColCoords3D bevCoords = jOD.new ColCoords3D();
								if (jCFS.hasBevel && ! mFC.columnHasNotBeenFound) bevCoords = jOD.findColumnsBevel(mFC, roughCoords, jCFS);
								else bevCoords = roughCoords;
								
								//if everything went well so far, flag column as found and break while loop
								if (!mFC.columnHasNotBeenFound) {
									IJ.log("Success!! Column has been found!!"); //add a new line in the log
									findingAttemptNumber = 100;  // break while loop..
																
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
							}
						//}
						//catch(Exception e) {
						//	IJ.log("Algorithm has crashed.. trying to restart...");
						//}
					}
				} // if not steel
				
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
				
				IJ.log(""); //add a new line in the log
				
				//switch warning lamp off again
				mFC.columnHasNotBeenFound = false;
 			}
//			catch(Exception e) {
//				failedColumnDetection.add(mFC.myTiffs[i]);
//			}
//		}
		
		if (failedColumnDetection.size() == 0) {
			IJ.log("\nIt appears that the column detection went smoothly for all samples.\n "
					+ "Please check the column detection performance in WallCoordinateidentified/ReviewFoundOutlines "
					+ "if you are happy with the found column outlines.");			
		}
		else {
			String eMsg = "No column could be detected for the following samples:\n\n";
			for (int i = 0 ; i < failedColumnDetection.size() ; i++) eMsg += failedColumnDetection.get(i) + "\n";
			eMsg += "\n Please re-run the column detection for these columns using a lower R2 value "
					+ "or using less strict criteria for a successful column find. ";		
			IJ.log(eMsg);
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

