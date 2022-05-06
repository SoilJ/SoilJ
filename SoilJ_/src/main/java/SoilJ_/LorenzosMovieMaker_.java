package SoilJ_;

import java.io.File;

import org.apache.commons.math3.stat.StatUtils;

import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

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

public class LorenzosMovieMaker_ extends ImagePlus implements PlugIn  {

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
		MenuWaiter menu = new MenuWaiter();
		InputOutput jIO = new InputOutput();
		
		//tell me what I should do!
		MenuWaiter.LorenzosOptions mLO = menu.selectLorenzosOptions();
		if (mLO == null) return;
		
		int lastRun = 1;
		String[] myRuns = {"none"};
		if (mLO.runInBatchMode) {
			myRuns = jIO.listFoldersInFolder(new File(mLO.myAboveFolder));
			
			lastRun = myRuns.length;
		}
		
		for (int k = 0 ; k < lastRun ; k++) {
			
			if (!mLO.runInBatchMode | (!myRuns[0].equalsIgnoreCase("none") & !myRuns[k].contains("hide") & !myRuns[k].contains("Inf"))) {
				
				if (myRuns[0].equalsIgnoreCase("none")) {
	
					mLO.myTiffs = jIO.listTiffsInFolder(new File(mLO.myInFolder));
					mLO.myName = mLO.myInFolder.replace(mLO.myAboveFolder + pathSep,"");				
									
				}
				
				else {
					
					mLO.myInFolder = mLO.myAboveFolder + pathSep + myRuns[k];
					mLO.myTiffs = jIO.listTiffsInFolder(new File(mLO.myInFolder));
					mLO.myName = mLO.myInFolder.replace(mLO.myAboveFolder + pathSep,"");
					
				}				
							
				//loop over dry images
				ImagePlus testImage = jIO.openTiff2D(mLO.myInFolder + pathSep + mLO.myTiffs[0]);
				ImageStack dryStack = new ImageStack(testImage.getWidth(), testImage.getHeight());
				IJ.showStatus("Checking for dry images for " + mLO.myName + " ...");
				int dryTiffNumber = 0;
				for (int i = 0 ; i < mLO.myTiffs.length ; i++) {  //myTiffs.length
					
					String testString = mLO.myTiffs[i].substring(0, 1);
					
					if (testString.equalsIgnoreCase("1")) {
						ImagePlus nowTiff = jIO.openTiff2D(mLO.myInFolder + pathSep + mLO.myTiffs[i]);
						ImageProcessor nowIP = nowTiff.getProcessor();
						dryStack.addSlice(nowIP);
						
						dryTiffNumber++;
					}
				}
				ImagePlus dryTiffs = new ImagePlus();
				dryTiffs.setStack(dryStack);
				
				//calculate mean dry image
				IJ.showStatus("Calculating average dry image for " + mLO.myName + " ...");
				ImageProcessor dryIP = new ShortProcessor(dryTiffs.getWidth(), dryTiffs.getHeight());
				for (int x = 0 ; x < dryTiffs.getWidth() ; x++) {
					for (int y = 0 ; y < dryTiffs.getHeight() ; y++) {
						
						double[] allDryAtThisPixel = new double[dryTiffs.getNSlices()];
						
						for (int i = 0 ; i < dryTiffs.getNSlices() ; i++) {
							
							dryTiffs.setSlice(i);
							ImageProcessor nowIP = dryTiffs.getProcessor();
							allDryAtThisPixel[i] = nowIP.getPixel(x, y);					
						}
						
						dryIP.putPixel(x, y, (int)StatUtils.mean(allDryAtThisPixel));
					}
				}	
				ImagePlus dryTiff = new ImagePlus("DryTiff", dryIP);
				
				//check if average is smaller than stepsize
				int average = mLO.average;
				if (average > mLO.stepsize) average = mLO.stepsize;		
				
				//loop over wet images		
				ImageStack wetStack = new ImageStack(testImage.getWidth(), testImage.getHeight());
				int cc = 0;		
				for (int i = 0 ; i < mLO.myTiffs.length ; i += mLO.stepsize) {  //myTiffs.length
			
					if (!mLO.myTiffs[i].substring(0).equalsIgnoreCase("1")) {				
						
						cc++;
						IJ.showStatus("Creating infiltration image of " + mLO.myName + ", frame " + cc + " of " + (int)Math.floor((double)(mLO.myTiffs.length - dryTiffNumber) / mLO.stepsize + 1) + " ...");
						
						ImageProcessor dry = dryTiff.getProcessor();
						
						ImageProcessor water = new FloatProcessor(dry.getWidth(), dry.getHeight());
						
						for (int j = 0 ; j < average ; j++) {
							
							if (i + j < mLO.myTiffs.length) {
							
								ImagePlus nowTiff = jIO.openTiff2D(mLO.myInFolder + pathSep + mLO.myTiffs[i + j]);
											
								ImageProcessor wet = nowTiff.getProcessor();
												
								for (int x = 0 ; x < dryTiffs.getWidth() ; x++) {
									for (int y = 0 ; y < dryTiffs.getHeight() ; y++) {
										
										float nowWater = water.getPixelValue(x, y);
										
										float dryPix = dry.getPixel(x, y);
										float wetPix = wet.getPixel(x, y);
										
										float newValue = (1f / (float)average) * (float)(wetPix - dryPix) + nowWater;
										
										water.putPixelValue(x,y,newValue);						
									}
								}
							}
						}
						
						wetStack.addSlice(water);
					}		
				}		
				
				ImagePlus waterTiff = new ImagePlus();
				waterTiff.setStack(wetStack);
			
				//save infiltration movie..
				jIO.tiffSaver(mLO.myOutFolder, mLO.myName + ".tif", waterTiff);
				
				IJ.showStatus("Infiltration image created!");
			}
		}
		
		IJ.showStatus("All infiltration images created!");
	}
}