package SoilJ_;

import SoilJ.tools.HistogramStuff;
import SoilJ.tools.InputOutput;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

/** 
 * PlotVerticalProfile is a SoilJ plugin that extracts vertical profiles of the greyscale
 * statistics within horizontal slices of 3-D images. 
 * 
 * @author John Koestel
 *
 */

public class CorrectSand_ extends ImagePlus implements PlugIn  {

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
		
		InputOutput jIO = new InputOutput();
		HistogramStuff hist = new HistogramStuff();
		
		InputOutput.MyFileCollection mFC = jIO.new MyFileCollection();
		mFC.nowTiffPath = "Z:\\FreezingMRI\\MRI\\Quantification\\SA1_warped16.tif";
		mFC.nowWidth = 256;
		mFC.nowHeight = 417;
		int times = 50;
		
		ImagePlus[] nowTiff = new ImagePlus[times];
		
		int[] nowRange = new int[248];
		
		for (int j = 0 ; j < 248 ; j++) nowRange[j] = j;
		
		for (int i = 0 ; i < times ; i++) {
			
			for (int j = 0 ; j < 248 ; j++) nowRange[j] = nowRange[j] + 248;
		
			nowTiff[i] = jIO.openTiff3DSomeSlices(mFC, nowRange);
		}
		
		//find histogram peak		
		int[] lowerRef = new int[times];
		int[] upperRef = new int[times];
		for (int i = 0 ; i < times ; i++) {
			
			int[] hist3D = new int[256*256];
						
			//lowerRef
			for (int z = 0 ; z < 248 ; z++) {
				
				nowTiff[i].setPosition(z+1);
				ImageProcessor nowIP = nowTiff[i].getProcessor();
				nowIP.resetRoi();
				
				int[] histo = nowIP.getHistogram();
				
				for (int j = 0 ; j < histo.length ; j++) hist3D[j] += histo[j];
				
			}
			
			lowerRef[i] = hist.findModeFromHistogram(hist3D);
			
			//upperRef
			nowTiff[i].setPosition(50);
			ImageProcessor nowIP = nowTiff[i].getProcessor();
			
			nowIP.setRoi(213, 97, 7, 13);
			
			upperRef[i] = hist.findMaxFromHistogram(nowIP.getHistogram());
			
		}
		
		//do brightness correction
		int newLow = 600;
		int newHigh = 8000;
			
		for (int i = 0 ; i < times ; i++) {
			
			//lowerRef
			for (int z = 0 ; z < 248 ; z++) {
				
				nowTiff[i].setPosition(z+1);
				ImageProcessor nowIP = nowTiff[i].getProcessor();
									
				for (int x = 0 ; x < nowTiff[0].getWidth() ; x++) {					
					for (int y = 0 ; y < nowTiff[0].getHeight(); y++) {
						
						int nowPix = nowIP.getPixel(x, y);
						double nowCorrect = ((double)nowPix - (double)lowerRef[i]) / ((double)upperRef[i] - (double)lowerRef[i]);
						
						int newPix = (int)Math.round(nowCorrect * (double)(newHigh - newLow) + (double)newLow);
							
						nowIP.putPixel(x, y, newPix);
						
					}
				}
			}
			
			nowTiff[i].updateAndDraw();
			
			//save it
			mFC.myOutFolder = "Z:\\FreezingMRI\\MRI\\Quantification";
			mFC.fileName = "SA" + (1000 + i) + ".tif";
			jIO.tiffSaver(mFC, nowTiff[i]);
			
		}
		
	
		
		
		
		
		
	}
	
}