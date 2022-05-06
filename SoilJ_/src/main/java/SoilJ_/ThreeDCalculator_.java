package SoilJ_;

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
 * is a SoilJ class
 * 
 * @author John Koestel
 *
 */

public class ThreeDCalculator_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

	/*	String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";
		
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = SimpleThresholding_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ImageManipulator jIM = new ImageManipulator();
		MenuWaiter menu = new MenuWaiter();
						
		// init variables
		MenuWaiter.Calc3DMenuReturn m3D;
		int i;
		
		//construct image related objects
		ImagePlus nowTiffA = new ImagePlus();
		ImagePlus nowTiffB = new ImagePlus();
		ImagePlus outTiff = new ImagePlus();;
		
		//read base folder and number of 3D images
		String myBaseFolderA = jIO.chooseAFolder("Please choose the folder with image A ..");
		String[] myTiffs0 = jIO.listTiffsInFolder(new File(myBaseFolderA));
		
		String myBaseFolderB = jIO.chooseAFolder("Please choose the folder with image B ..");
		
		//ask for threshold choice
		m3D = menu.show3DCalcDialog(myBaseFolderA, myBaseFolderB);
	
		//if not all tiffs shall be thresholded
		String[] myTiffs = null;
		if (m3D.filterImages) {
			myTiffs = jIO.filterFolderList(myTiffs0, m3D.filterTag);
		} 
		else myTiffs = myTiffs0;
		
		//read gauges in the gauge folder
		String myGaugeFolder = myBaseFolderA;
		String[] myGauges = null;
		if (m3D.useInnerCircle) {
			myGaugeFolder = jIO.findTheInnerCircle(myBaseFolderA);		
			myGauges = jIO.listInnerCircleFiles(myGaugeFolder, m3D.filterTag);
		}
		
		String myOutFolder = "Calc" + m3D.operationTag;
		String myOutPath = myBaseFolderA + myOutFolder;						
		new File(myOutPath).mkdir();	
		
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();
			
			//load image
			nowTiffA = jIO.openTiff3D(myBaseFolderA + pathSep + myTiffs[i]);
			nowTiffB = jIO.openTiff3D(myBaseFolderB + pathSep + myTiffs[i]);
			
			//select the correct gauge and surface files	
			if (m3D.useInnerCircle){
				String[] myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC);
				nowGauge = myGauges[myGandS[0]];
			}
			
			//apply segmentation			
			outTiff = jIM.calc3D(nowTiffA, nowGauge, m3D, nowTiffB);			 
		
			//save result
			jIO.tiffSaver(myOutPath, myTiffs[i], outTiff);
			
		}*/
		
	}
	

	
}