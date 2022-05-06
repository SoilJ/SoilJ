package SoilJ_;


	import java.io.File;

import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.MorphologyAnalyzer;
import SoilJ.tools.ObjectDetector;
import SoilJ.tools.RoiHandler;
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

public class ExtractVerticalProfileStatistics_ extends ImagePlus implements PlugIn  {

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
		ObjectDetector jOD = new ObjectDetector();	
		MorphologyAnalyzer jMA = new MorphologyAnalyzer();
		
		MorphologyAnalyzer.ProfileStatistics mPS = jMA.new ProfileStatistics();
				
		//tell me what I should do!
		MenuWaiter.ROISelectionOptions mRSO = menu.regionOfInterestSelection();
		if (mRSO == null) return;
		
		//get some bio pore extraction parameters
		//MenuWaiter.AnalyzeMeanGrayValueMenu mAGV = menu.showAnalyzeMeanGrayValueMenu();
		
		//create Folder structure
		InputOutput.MyFileCollection mFC = jIO.fileSelector("Please choose a file or folder with your image data");
		mFC = jIO.getAllMyNeededFolders(mFC.myBaseFolder, mFC.myTiffs, "", mRSO.useInnerCircleFiles, mRSO.includeSurfaceTopography);   //"" put because no possibility for a filterTag implemented yet
				
		mFC.myOutFolder = mFC.myBaseFolder + pathSep + "VerticalProfileStats"; 	 	
		new File(mFC.myOutFolder).mkdir();
					
		//loop over 3D images
		for (int i = 0 ; i < mFC.myTiffs.length ; i++) {  //myTiffs.length
			
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
			RoiHandler.ColumnRoi colRoi = roi.prepareDesiredRoi(mFC, jIO.openTiff3DSomeSlices(mFC, colSlices), mRSO, "null", "");
					
			//extract mean gray value
			mPS = jMA.findProfileStatistics(colRoi.nowTiff, colRoi.pRoi);
			
			//save it
			String fname = mFC.colName + ".mgv";
			jIO.writeVerticalProfileStats(mFC.myOutFolder + pathSep + fname, mPS);
			
		}
			
		IJ.showStatus("Mean gray values are calculated and saved!");
	}
}
		
