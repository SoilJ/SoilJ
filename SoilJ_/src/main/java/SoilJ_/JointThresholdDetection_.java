package SoilJ_;

import java.io.File;

import SoilJ.tools.DisplayThings;
import SoilJ.tools.HistogramStuff;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
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
import ij.process.AutoThresholder;

/** 
 * HistoGrammar is a SoilJ plugin that extracts the histograms of an ensemble of 3-D X-ray images and writes them into 
 * one ASCII file. The joint histogram may be optionally calculated.
 * 
 * @author John Koestel
 *
 */

public class JointThresholdDetection_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";
		
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = JointThresholdDetection_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();
		MenuWaiter menu = new MenuWaiter();
		InputOutput.MyFileCollection mFC = jIO.new MyFileCollection();
		HistogramStuff hist = new HistogramStuff();
		AutoThresholder thresh = new AutoThresholder();
								
		// init variables
		int i;
		MenuWaiter.HistogramMenuReturn hMR; 
		
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();
		
		//ask for threshold choice
		hMR = menu.showHistogramDialog();
		if (hMR == null) return;
		
		//read base folder and number of 3D images
	    String[] myTiffs0 = null;String myBaseFolder = null;
	    while (myTiffs0 == null) {
	    	myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
			if (myBaseFolder == null) return;
			myTiffs0 = jIO.listTiffsInFolder(new File(myBaseFolder));
			if (myTiffs0 == null) {
				IJ.error("Please choose a folder with TIFF images or cancel. Thank you.");	
			}
		}
		
		//if not all tiffs shall be thresholded
		String[] myTiffs;
		if (hMR.filterImages) {myTiffs = jIO.filterFolderList(myTiffs0, hMR.filterTag);} else myTiffs = myTiffs0;
		
		//read gauges in the gauge folder
		mFC = jIO.getAllMyNeededFolders(myBaseFolder, myTiffs, hMR.filterTag, hMR.useInnerCircle, hMR.useSoilSurfaceFile);
		
		//set output path	
		String outPath = mFC.myBaseFolder + pathSep + "Histograms.16b";
		if (hMR.filterImages) outPath = mFC.myBaseFolder  + pathSep + hMR.filterTag + "Histograms.16b";
		mFC.pathSep = pathSep;
		mFC.myOutFolder = outPath;		
				
		//check if there is already an old version of the histogram file and flag it for overwriting in if it exists 
		String[] myFiles = jIO.listFilesInFolder(new File(mFC.myBaseFolder));
		
		//init histogram vectors
		int[] myBulkHistogram = new int[256 * 256];
		int[] nowHist = new int[256 * 256];
		int[][] allHists = new int [mFC.myTiffs.length][256 *256];
						
		//loop over 3D images
		for (i = 0 ; i < mFC.myTiffs.length ; i++) {  //myTiffs.length 

			//try to free up some memory			
			System.gc();
			
			//load image
			nowTiff = jIO.openTiff3D(myBaseFolder + pathSep + myTiffs[i]);	
						
			//select the correct gauge and surface files
			if (hMR.useInnerCircle) {
				String[] myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC);
				mFC.nowInnerCirclePath = mFC.myInnerCircleFolder + pathSep + myGandS[0];
			}
			
			//get 16-bit Histogram			
			nowHist = hist.sampleThisHistogram16(nowTiff, myTiffs[i], mFC, hMR);
						
			//add to bulk histogram of all images
			for (int j = 0 ; j < nowHist.length ; j++) myBulkHistogram[j] += nowHist[j];
			
			//save individual histograms			
			jIO.writeHistogram(mFC, nowHist);
		
		}
		
		
		
		//find thresholds.. still under construction
		RollerCaster cast = new RollerCaster();
		double[] medianFilteredBulkHistogram = cast.castInt2Double(myBulkHistogram);//, 11);  //maths.oneDMedianFilter
		double[] mFBH256 = hist.convert16to8BitHistogram(medianFilteredBulkHistogram);
		double[] greyValues = new double[medianFilteredBulkHistogram.length];
		for (i = 0 ; i < greyValues.length ; i++) greyValues[i] = i;		
		
		double[] myThreshes = new double[9];
		myThreshes[0] = 256 * thresh.getThreshold(AutoThresholder.Method.Default, cast.castDouble2Int(mFBH256));
		myThreshes[1] = 256 * thresh.getThreshold(AutoThresholder.Method.Otsu, cast.castDouble2Int(mFBH256));		
		myThreshes[2] = 256 * thresh.getThreshold(AutoThresholder.Method.Huang, cast.castDouble2Int(mFBH256));
		myThreshes[3] = 256 * thresh.getThreshold(AutoThresholder.Method.MaxEntropy, cast.castDouble2Int(mFBH256));
		myThreshes[4] = 256 * thresh.getThreshold(AutoThresholder.Method.Minimum, cast.castDouble2Int(mFBH256));
		myThreshes[5] = 256 * thresh.getThreshold(AutoThresholder.Method.MinError, cast.castDouble2Int(mFBH256));
		myThreshes[6] = 256 * thresh.getThreshold(AutoThresholder.Method.RenyiEntropy, cast.castDouble2Int(mFBH256));
		myThreshes[7] = 256 * thresh.getThreshold(AutoThresholder.Method.Triangle, cast.castDouble2Int(mFBH256));
		myThreshes[8] = 256 * thresh.getThreshold(AutoThresholder.Method.IsoData, cast.castDouble2Int(mFBH256));
		
		DisplayThings disp = new DisplayThings();	
		disp.plotHistogramWithThresholds(greyValues, medianFilteredBulkHistogram, myThreshes, "HistogramAndThresholds", "grey value", "frequency");
		
	}
	
}