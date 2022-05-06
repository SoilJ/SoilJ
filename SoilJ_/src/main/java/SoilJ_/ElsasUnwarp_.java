package SoilJ_;

import SoilJ.tools.MenuWaiter;
import bunwarpj.bUnwarpJ_;
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

public class ElsasUnwarp_ extends ImagePlus implements PlugIn  {

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
		
		// init variables
		int i;
		
		//tell me what I should do!
		MenuWaiter.ElsasOptions mEO = menu.selectElsaOptions();
		if (mEO == null) return;
		
		//loop over 3D images
		for (i = 0 ; i < mEO.myTiffs.length ; i++) {  //myTiffs.length
	
			if (!mEO.myTiffs[i].equalsIgnoreCase(mEO.myTargetImage)) {
				
				String nowPath = mEO.myInFolder + pathSep + mEO.myTiffs[i];			
				String outPath = mEO.myOutFolder + pathSep + mEO.myTiffs[i];
				bUnwarpJ_.elasticTransformImageMacro(mEO.myTargetImage, nowPath, mEO.myTransformationFile, outPath);
				
			}
			
		}		
	}
}