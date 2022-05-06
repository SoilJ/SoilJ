/*
 * To the extent possible under law, the Fiji developers have waived
 * all copyright and related or neighboring rights to this tutorial code.
 *
 * See the CC0 1.0 Universal license for details:
 *     http://creativecommons.org/publicdomain/zero/1.0/
*/ 

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
 
/**
 *SoilJ is a plugin for the semi-automatized processing of 3-D X-ray images of soil columns
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
 *
 *
 * @author The Fiji Team (template programming..)
 * @John Koestel 
 */

public class SoilJ_ implements PlugInFilter {


	@Override
	public int setup(String arg, ImagePlus imp) {
		return 1;
	}

	@Override
	public void run(ImageProcessor ip) {
		IJ.error("AAAAAAAAAAAARRRRGGGHHHHHH!!!!!!!!!");
	}

	/**
	 * Start SoilJ!
	 * @param args
	 */
	
	public static void main(String[] args) throws Exception {

		// set the plugins.dir property to make the plugin appear in the Plugins menu
		// see: https://stackoverflow.com/a/7060464/1207769
		//Class<?> clazz = SoilJ_.class;	
		
		//System.setProperty("plugins.dir", "D:\\Eclipse\\Maven");
		System.setProperty("plugins.dir", "D:\\GitSoilJ");
		//Debug.run("Soil_", "");
		
		// start ImageJ
		new ImageJ();

		// run the plugin 
		//IJ.runPlugIn(clazz.getName(), "");
		 
	}
}












