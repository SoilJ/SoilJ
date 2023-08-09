package SoilJ_;

import SoilJ.tools.DisplayThings;
import SoilJ.tools.HistogramStuff;
import SoilJ.tools.InputOutput;
import SoilJ.tools.InputOutput.Histograms;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.ObjectDetector;
import SoilJ.tools.RoiHandler;
import SoilJ.tools.RollerCaster;
import SoilJ.tools.TailoredMaths;

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
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;

/** 
 * ExtractPoresizeDistributio is a SoilJ plugin that extracts the pore size distribution from an
 * image depicting pore diameters (thicknesses) and writes it into an ASCII file.
 * 
 * @author John Koestel
 *
 */

public class Analyze1DHistograms_ extends ImagePlus implements PlugIn  {

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
		ObjectDetector jOD = new ObjectDetector();
		MenuWaiter menu = new MenuWaiter();
		RoiHandler roi = new RoiHandler();
		HistogramStuff hist = new HistogramStuff();
		RollerCaster rC = new RollerCaster();
		TailoredMaths math = new TailoredMaths();
		DisplayThings disp = new DisplayThings();	
		ResultsTable myTable = new ResultsTable();
		
		//select a histogram file
		InputOutput.MyFileCollection mFC = jIO.selectAHistogramFile("Please select a histogram file (.8b or .16b)");
		
		//load histogram
		Histograms myHist = jIO.readHistogram(mFC);
		
		//calculate joint histogram.
		if (myHist.numberOfHistograms == 1) myHist.jHisto = rC.castInt2Float(myHist.histograms[0]); 
		else myHist = hist.calcJointHistogram(myHist.floatHistograms);		
		
		//calculate thresholds
		HistogramStuff.Thresholds myThreshes = hist.calcThreshes(myHist.jHisto);
			
		//plot joint histogram		
		disp.plotHistogramAdvanced(myHist.jHisto, 0); 
		
		//save histogram
		if (mFC.bitDepth == 8) jIO.writeJointHistogram8(mFC, rC.castFloat2Int(myHist.jHisto));
		if (mFC.bitDepth == 16) jIO.writeJointHistogram(mFC, rC.castFloat2Int(myHist.jHisto));
		
		//plot table with thresholds
		ResultsTable tab = new ResultsTable();
		
		tab.incrementCounter();
		tab.addValue("Otsu", myThreshes.otsu);
		tab.addValue("Minimum", myThreshes.minimumJohns);
		tab.addValue("Isodata", myThreshes.isodata);
		tab.addValue("Renyi Entropy", myThreshes.renyiEntropy);		
		tab.addValue("Maximum Entropy", myThreshes.maxEntropy);
		tab.addValue("Triangle", myThreshes.triangle);
		tab.addValue("Huang", myThreshes.huang);
		tab.addValue("Li", myThreshes.li);
		tab.addValue("Yen", myThreshes.yen);
		
		tab.show("Thresholds");
			
	}
}