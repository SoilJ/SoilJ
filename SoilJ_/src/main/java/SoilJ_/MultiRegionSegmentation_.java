package SoilJ_;

import java.io.File;

import SoilJ.tools.ImageManipulator;
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
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

/** 
 *  
 * @author John Koestel
 */

public class MultiRegionSegmentation_ extends ImagePlus implements PlugIn {
	
	public void run(String arg) {
		
		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";
		
		//This is how command line parameter transfer needs to be done:		
		//run("MultiRegionSegmentation", "nomenu joint otsu maxentro minerror classes=3 path=Y:/SoilJ/8Bit/Histo_EntireImage/Histograms.8b")
		GenericDialog gd = new GenericDialog("");
		gd.addCheckbox("nomenu", false);
		gd.addCheckbox("joint", true);
		gd.addCheckbox("otsu", true);
		gd.addCheckbox("maxentro", true);
		gd.addCheckbox("minerror", true);
		gd.addNumericField("classes", 3, 0);
		gd.addStringField("path", "myPath");
		
		boolean noMenu = gd.getNextBoolean();
		boolean joint = gd.getNextBoolean();
		boolean otsu = gd.getNextBoolean();
		boolean maxEntro = gd.getNextBoolean();
		boolean minError = gd.getNextBoolean();
		int classes = (int)gd.getNextNumber();
		String path = gd.getNextString();
		
				
		//adapt path to pathSep
		String newPath = "";
		if (noMenu) {			
			for (int i = 0 ; i < path.length() ; i++) {
			
				if (path.substring(i,i+1).equalsIgnoreCase("/")) newPath += pathSep;
				else newPath += path.substring(i,i+1);
			
			}			
		}
		
		//set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = Extract2DHistograms_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();
		MenuWaiter menu = new MenuWaiter();
		ImageManipulator jIM = new ImageManipulator();
		RollerCaster rC = new RollerCaster();
		
		//tell me what to do..
		MenuWaiter.MultiRegionSegmentation mRS = menu.new MultiRegionSegmentation();
		
		//catch command line start..
		String myHistFile = "";
		String nowDir = IJ.getDirectory("current");
		if (noMenu) {
			mRS.useJointHistogram = joint;
			mRS.calcOtsu = otsu;
			mRS.calcMaxEntro = maxEntro;
			mRS.calcMinErr = minError;
			mRS.numberOfThresholds = classes - 1;
			myHistFile = newPath;
		}
		else {
			mRS = menu.showMultiRegionSegmentationMenu();	
			if (mRS == null) return;
			myHistFile = jIO.chooseAFile("Please choose the file with your 1D histograms", nowDir, "Histograms.16b");			
		}
		
		//init folder collection
		InputOutput.MyFileCollection mFC = jIO.new MyFileCollection();
		    	 
		File nowFile = new File(myHistFile);
		String myOutPath = nowFile.getParent();		
		
		//create output folder  				
  		mFC.myOutFolder = myOutPath; 	//also add it to the folder collection
  		
 		//save pathSep
  		mFC.pathSep = pathSep;
  		  			
		// load histograms	
  		InputOutput.Histograms hists = jIO.readHistogram(myHistFile);
  		
  		//define bitDepth
  		if (hists.histogramClasses > 256) mFC.bitDepth = 16;
		
		if (!mRS.useJointHistogram) {
			for (int i = 0 ; i < hists.numberOfHistograms ; i++) {  
				
				mFC.colName = "MultiThreshs_" + hists.sampleNames.get(i) + "_" + (mRS.numberOfThresholds + 1) + "Phases.asc";
				
				int[] nowHist = new int[hists.histogramClasses];
				for (int j = 0 ; j < hists.histogramClasses ; j++) nowHist[j] = hists.histograms[i][j];

				int[][] myThresholds = new int[3][mRS.numberOfThresholds];
				
				if (mRS.calcOtsu) myThresholds[0] = jIM.multiOtsu(mFC, nowHist, mRS);
				if (mRS.calcMaxEntro) myThresholds[1] = jIM.multiMaxEntropy(mFC, nowHist, mRS);
				if (mRS.calcMinErr) myThresholds[2] = jIM.multiMinError(mFC, nowHist, mRS);
				
				jIO.writeMultiThresholds(mFC, myThresholds, mRS);
				
			}
		}
		
		else {  //if joint histogram shall be used..
			
			double[] jointHist = new double[hists.histogramClasses];
			
			for (int i = 0 ; i < hists.numberOfHistograms ; i++) {
				for (int j = 0 ; j < hists.histogramClasses ; j++) {
					jointHist[j] += (double)hists.histograms[i][j] / hists.numberOfHistograms; 
				}
			}
			
			mFC.colName = "JointHistMultiThreshs_" + (mRS.numberOfThresholds + 1) + "Phases.asc";
			
			int[][] myThresholds = new int[3][mRS.numberOfThresholds];
	
			if (mRS.calcOtsu) myThresholds[0] = jIM.multiOtsu(mFC, rC.castDouble2Int(jointHist), mRS);
			if (mRS.calcMaxEntro) myThresholds[1] = jIM.multiMaxEntropy(mFC, rC.castDouble2Int(jointHist), mRS);
			if (mRS.calcMinErr) myThresholds[2] = jIM.multiMinError(mFC, rC.castDouble2Int(jointHist), mRS);
			
			jIO.writeMultiThresholds(mFC, myThresholds, mRS);

		}  		
		
		//try to clear memory
		IJ.freeMemory();IJ.freeMemory();	
		
	}
}