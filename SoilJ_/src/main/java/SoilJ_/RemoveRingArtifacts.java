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
 * RemoveRingArtifact is a SoilJ plugin that may one day be written to remove ring-shaped image artifacts
 * as they e.g. occur because of splashes of target material on a beryllium window of the X-ray tube.
 * 
 * @author John Koestel
 *
 */

public class RemoveRingArtifacts extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		/*
		 * String pathSep = "\\";
		 * 
		 * //probe operating system and adjust PathSep if necessary String myOS =
		 * System.getProperty("os.name"); if (myOS.equalsIgnoreCase("Linux")) pathSep =
		 * "/";
		 * 
		 * // set the plugins.dir property to make the plugin appear in the Plugins menu
		 * Class<?> clazz = RemoveRingArtifacts.class; String url =
		 * clazz.getResource("/" + clazz.getName().replace('.', '/') +
		 * ".class").toString(); String pluginsDir = url.substring(5, url.length() -
		 * clazz.getName().length() - 6); System.setProperty("plugins.dir", pluginsDir);
		 * 
		 * //construct biggish objects InputOutput jIO = new InputOutput();
		 * ImageManipulator jIM = new ImageManipulator(); MenuWaiter menu = new
		 * MenuWaiter();
		 * 
		 * // init variables int i;
		 * 
		 * //construct image related objects ImagePlus nowTiff = new ImagePlus();
		 * ImagePlus outTiff = new ImagePlus();
		 * 
		 * //ask for threshold choice MenuWaiter.RingArtifactRemoverOptions rARO =
		 * menu.returnRingArtifactRemoverOptions();
		 * 
		 * //read base folder and number of 3D images String[] myTiffs = null;String
		 * myBaseFolder = null; while (myTiffs == null) { myBaseFolder =
		 * jIO.chooseAFolder("Please choose the folder with your image data"); if
		 * (myBaseFolder == null) return; myTiffs = jIO.listTiffsInFolder(new
		 * File(myBaseFolder)); if (myTiffs == null) {
		 * IJ.error("Please choose a folder with TIFF images or cancel. Thank you."); }
		 * }
		 * 
		 * //read gauges in the gauge folder boolean useInnerCircle =
		 * rARO.useInnerCircle; InputOutput.MyFileCollection mFC =
		 * jIO.getAllMyNeededFolders(myBaseFolder, myTiffs, "", useInnerCircle, true);
		 * 
		 * //save and create output folder String myOutFolder; String myOutPath = null;
		 * myOutFolder = "RingRemoved"; myOutPath = myBaseFolder + myOutFolder; new
		 * File(myOutPath).mkdir();
		 * 
		 * //loop over 3D images for (i = 0 ; i < myTiffs.length ; i++) {
		 * //myTiffs.length
		 * 
		 * //try to free up some memory System.gc();
		 * 
		 * //load image nowTiff = jIO.openTiff3D(myBaseFolder + pathSep + myTiffs[i]);
		 * 
		 * //select the correct inner circle if this option has been selected.. String
		 * innerCirclePath = ""; if (rARO.useInnerCircle) { int[] myGandS =
		 * jIO.getTheCorrectGaugeNSurfaceFiles(myTiffs[i], mFC.myInnerCircleFiles,
		 * null); innerCirclePath = mFC.myInnerCircleFiles[myGandS[0]]; }
		 * 
		 * //outTiff = jIM.removeRingArtifacts(nowTiff, innerCirclePath, rARO);
		 * 
		 * //save result jIO.tiffSaver(myOutPath, myTiffs[i], outTiff); }
		 */
		
	}
	
}