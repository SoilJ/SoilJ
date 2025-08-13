package SoilJ_;


import java.io.File;

import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.RoiHandler;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import ij.plugin.Slicer;
import ij.plugin.frame.ContrastAdjuster;
import ij.process.ImageProcessor;

/** 
 * ExtractVerticalCrossSections is a SoilJ plugin.  
 * 
 * @author John Koestel
 *
 */

public class ExtractVerticalCrossSections_ extends ImagePlus implements PlugIn  {

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
		Slicer slicer = new Slicer();
		ContrastAdjuster jCA = new ContrastAdjuster();
		
		// init variables
		int i;
		
		//tell me what I should do!
		MenuWaiter.ROISelectionOptions mRSO = menu.regionOfInterestSelection();
		if (mRSO == null) return;
		
		//create Folder structure
		InputOutput.MyFileCollection mFC = jIO.createFolders4SubROIData(mRSO, true);				
		String myInspectionFolder =  mFC.myOutFolder + pathSep + "2DCrossSections";
		new File(myInspectionFolder).mkdir();			
		
		//loop over 3D images
		for (i = 0 ; i < mFC.myTiffs.length ; i++) {  //myTiffs.length
			
			//assemble image for analyses
			mFC.fileName = mFC.myTiffs[i];
			mFC = jIO.addCurrentFileInfo8Bit(mFC);
			
			//load file			
			int[] startStopSlices = jIO.findStartAndStopSlices(mFC, mRSO);
			int[] colSlices = new int[startStopSlices[1] - startStopSlices[0]];
			for (int j = 0 ; j < colSlices.length ; j++) colSlices[j] = startStopSlices[0] + j;
		
			//cut roi		
			mFC.startSlice = startStopSlices[0];
			mFC.stopSlice = startStopSlices[1];
			RoiHandler.ColumnRoi colRoi = roi.prepareDesiredRoi(mFC, jIO.openTiff3DSomeSlices(mFC, colSlices), mRSO, "null", "null");
			
			//reslice
			ImagePlus sliceTiff = slicer.reslice(colRoi.nowTiff);
			
			//get central slice
			sliceTiff.setSlice(sliceTiff.getNSlices() / 2);
			ImageProcessor outIP = sliceTiff.getProcessor();
			
			//adjust brightness
			outIP.setMinAndMax(5000, 25000);
			
			//convert to RGB
			outIP = outIP.convertToRGB();			
			
			//convert to Tiff
			ImagePlus outTiff = new ImagePlus("", outIP);			

			//outTiff.updateAndDraw();
			//outTiff.show();
			
			//wait a bit..
			IJ.wait(1000);

			//and save... oufff!!
			jIO.save2DTiff(myInspectionFolder, mFC.fileName, outTiff);
			
			IJ.showStatus("A vertical column cross-section has been cut out and saved..");
			
		}		
	}
}
		
