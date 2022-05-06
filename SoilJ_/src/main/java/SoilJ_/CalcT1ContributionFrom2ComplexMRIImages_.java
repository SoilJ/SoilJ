package SoilJ_;

import ij.ImagePlus;
import ij.plugin.PlugIn;

/** 
 * PlotVerticalProfile is a SoilJ plugin that extracts vertical profiles of the greyscale
 * statistics within horizontal slices of 3-D images. 
 * 
 * @author John Koestel
 *
 */

public class CalcT1ContributionFrom2ComplexMRIImages_ extends ImagePlus implements PlugIn  {

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
		
		
		
		
	}
	
}