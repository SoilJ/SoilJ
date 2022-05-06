package SoilJ_;

import java.io.File;

import SoilJ.tools.ArtificialPoreNetworkCreator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;

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
import ij.plugin.PlugIn;

/** 
 * GenerateRandomPoreClusters is a SoilJ plugin that creates random pore networks.
 * 
 * @author John Koestel
 *
 */

public class GenerateRandomPoreClusters_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		String pathSep = "\\";
		
		//probe operating system and adjust PathSep if necessary
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux")) pathSep = "/";
		
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = GenerateRandomPoreClusters_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();	
		MenuWaiter menu = new MenuWaiter();		
		ArtificialPoreNetworkCreator creator = new ArtificialPoreNetworkCreator();
		
		MenuWaiter.RandomClusterGenerator mRCG = menu.new RandomClusterGenerator();
				
		// init variables
		int i;
		
		//get information on what shall be done..
		mRCG = menu.showRandomClusterGeneratorMenu();
		if (mRCG == null) return;
		
		//read base folder and number of 3D images
		String myBaseFolder = jIO.chooseAFolder("Please choose the folder into which you want to save your images");		
		String myOutFolder = myBaseFolder + mRCG.shape + mRCG.domainX + "x" + mRCG.domainY + "x" + mRCG.domainZ;
		new File(myOutFolder).mkdir();		
		
		//in case that the porosities should be taken from a list..
		if (mRCG.mode.equalsIgnoreCase("predefinedList")) {
			String myPorosityListFile = jIO.chooseAFile("Please choose the file that contains the porosities..",myBaseFolder,"*.*");
			double[] myPorosities = jIO.readSingleColumnDoublesFromAscii(myPorosityListFile);	
			mRCG.porosityList = myPorosities;
		}
		
		//loop over 3D images
		for (i = 0 ; i < mRCG.numOfCopies ; i++) {  //myTiffs.length
			
			//apply segmentation			
			creator.createASetOfRandomFields(myOutFolder, i, mRCG);
			
			IJ.freeMemory();IJ.freeMemory();
			
		}
		
	}
}