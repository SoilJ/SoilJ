package SoilJ.tools;

/**
 *SoilJ.tools is a collection of classes for SoilJ, 
 *a collection of ImageJ plugins for the semi-automatized processing of 3-D X-ray images of soil columns
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
import ij.ImageStack;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.PolygonRoi;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;

/** 
 * HistogramStuff is a SoilJ class contains all subroutines that deal with histogram information.
 * 
 * @author John Koestel
 *
 */

public class HistogramStuff implements PlugIn {

	
	public void run(String arg) {
				//ok, this is not needed..
	}

	public class IlluminationInfo {
		
		public int min;
		public int quantile01;
		public int quantile05;
		public int quantile10;
		public int lowerQuartile;
		public int median;
		public int upperQuartile;
		public int quantile90;
		public int quantile95;
		public int quantile99;
		public int max;
		
		public int specialQuantile;
		
		public int mean;
		
	}
	
	public int findQuantileFromHistogram(int[] myHist, double quantile) {
		
		double[] cumHist = calcCumulativeHistogram(myHist);
		int myQuantile = findPercentileFromCumHist(cumHist, quantile);
		
		return myQuantile;
	}

	public int findPercentileFromCumHist(double[] cumHist, double threshold) {		
		
		int myThresh = 0;		
		double maxCumHist = cumHist[cumHist.length - 1]; 
		
		for (int i = 1 ; i < cumHist.length ; i++) {
			if (cumHist[i] / maxCumHist >= threshold) {
				myThresh = i;
				break;
			}
		}
		
		return myThresh;		
	}
	
	public ImagePlus extractJoint2DHistogramFrom2DHistograms(InputOutput.MyFileCollection mFC) {
		
		InputOutput jIO = new InputOutput();
		
		ImageProcessor meanIP = new FloatProcessor(6554, 6554);
		
		for (int i = 0 ; i < mFC.myTiffs.length ; i++) {
			
			String startOfTiff = mFC.myTiffs[i].substring(0, 4);
			
			if (!startOfTiff.equalsIgnoreCase("Join") & !startOfTiff.equalsIgnoreCase("POMR")) {
				
				IJ.showStatus("Processing 2D histogram of " + mFC.myTiffs[i]);
				
				String nowTiffPath = mFC.myBaseFolder + mFC.pathSep + mFC.myTiffs[i];
				ImagePlus nowTiff = jIO.openTiff2D(nowTiffPath);

				ImageProcessor nowIP = nowTiff.getProcessor();
				
				
				for (int x = 0 ; x < 6554 ; x++) {
					for (int y = 0 ; y < 6554 ; y++) {
				
						double newPixelValue = meanIP.getPixelValue(x, y) + nowIP.getPixelValue(x, y);
						meanIP.putPixelValue(x, y, newPixelValue);
						
					}
				}
			}
			
		}
		
		ImagePlus outTiff = new ImagePlus("Joint2DHistogram", meanIP);
		
		return outTiff;
		
	}
	
	public int[] extractHistograms16(InputOutput.MyFileCollection mFC, ImagePlus nowTiff) {
		
		RollerCaster rC = new RollerCaster();
		double[] histo = new double[(int)Math.round(Math.pow(2, 16))];
		
		//get histogram
		for (int j = 0 ; j < nowTiff.getNSlices() ; j++) { 
		
			nowTiff.setPosition(j+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			int[] nowHist = nowIP.getHistogram();
			nowHist[0] = 0;  //set zero entries to 0
			for (int i = 0 ; i < histo.length ; i++) histo[i] += (double)nowHist[i]*100/nowTiff.getNSlices();
			
		}		
		
		return rC.castDouble2Int(histo);
		
	}
	
	public int[] extractHistograms8(InputOutput.MyFileCollection mFC, ImagePlus nowTiff) {
		
		RollerCaster rC = new RollerCaster();
		double[] histo = new double[(int)Math.round(Math.pow(2, 8))];
		
		//get histogram
		for (int j = 0 ; j < nowTiff.getNSlices() ; j++) { 
		
			nowTiff.setPosition(j+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
		
			int[]  nowHist = nowIP.getHistogram();
			nowHist[0] = 0;  //set zero entries to 0
			for (int i = 0 ; i < histo.length ; i++) histo[i] += (double)nowHist[i]*100/nowTiff.getNSlices();
			
		}		
		
		return rC.castDouble2Int(histo);
		
	}
	
	public void analyze2DHistograms(InputOutput.MyFileCollection mFC) {
				
		InputOutput jIO = new InputOutput();
		ImageManipulator jIM = new ImageManipulator();
		
		ImageProcessor sumIP = new FloatProcessor(6554, 6554);
		
		int numOfGoodTiffs = 0;		

		for (int i = 0 ; i < mFC.myTiffs.length ; i++) {
			
			numOfGoodTiffs++;
				
			IJ.showStatus("Processing 2D histogram of " + mFC.myTiffs[i]);
			
			String nowTiffPath = mFC.myBaseFolder + mFC.pathSep + mFC.myTiffs[i];
			ImagePlus nowTiff = jIO.openTiff2D(nowTiffPath);

			ImageProcessor nowIP = nowTiff.getProcessor();
			
			for (int x = 0 ; x < 6554 ; x++) {
				for (int y = 0 ; y < 6554 ; y++) {
			
					double newPixelValue = sumIP.getPixelValue(x, y) + nowIP.getPixelValue(x, y) / mFC.myTiffs.length;
					sumIP.putPixelValue(x, y, newPixelValue);
					
				}
			}
						
		}
		
		//flip result in the vertical direction
		sumIP.flipVertical();

		//save mean Tiff
		ImagePlus outTiff = new ImagePlus("Mean2DHistogram", sumIP);
		
		//scale the result in the vertical by factor 0.5		
		outTiff = jIM.scaleWithoutMenu(outTiff, outTiff.getWidth(), outTiff.getHeight() / 2, 1, ImageProcessor.BILINEAR);
		
		mFC.nowTiffPath = mFC.myOutFolder + mFC.pathSep + "MeanHistogram2D.tif";
		jIO.save2DTiff(mFC, outTiff);
		
	}
	
	public int[][] extract2DHistogram(ImagePlus nowTiff, ImagePlus gradTiff) {
		
		int[][] hist2D = new int[256][256];
		
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			
			IJ.showStatus("Extracting 2D histogram at layer " + (z+1) + "/" + nowTiff.getNSlices());
			
			nowTiff.setSlice(z + 1);
			ImageProcessor mP = nowTiff.getProcessor();
			
			gradTiff.setSlice(z + 1);
			ImageProcessor gP = gradTiff.getProcessor().convertToShort(false);
			
			for (int x = 0 ; x < mP.getWidth(); x++) {
				
				for (int y = 0 ; y < mP.getHeight(); y++) {
					
					int mPix = (int)Math.floor((float)mP.getPixel(x, y) / 256);
					int gPix = (int)Math.floor((float)gP.getPixel(x, y) / 25.6);
					
					if (gPix > 255) gPix = 255;
					
					hist2D[mPix][gPix]++;
					
				}
			}			
		}		
		
		return hist2D;
		
	}
	
	public int findModeFromHistogram(int[] myHist) {
		
		int max = 0;
		int maxVal = 0;
		
		for (int i = 1 ; i < myHist.length ; i++) {
			if (myHist[i] > maxVal) {
				maxVal = myHist[i];
				max = i;
			}
		}		
		
		return max;
		
	}

	public int findPercentileFromTiff(ImageStack nowStack, double threshold) {
		
		ImagePlus nowTiff = new ImagePlus();
		nowTiff.setStack(nowStack);
		
		int[] myHist = sampleHistogram(nowTiff);
		double[] cumHist = calcCumulativeHistogram(myHist);		
		int myP = findPercentileFromCumHist(cumHist,threshold);
		
		return myP;
	}

	public int findNonZeroMinimumFromHistogram(int[] myHist, int largerThan) {
		
		int myMinimum = 0;
		int cc = 1;
		
		while (myMinimum == 0) {
			if (myHist[cc] > largerThan) {	
				myMinimum = cc;				
			}
			cc++;
		}
		
		return myMinimum;
	}

	public int findNonZeroModeFromHistogram(int[] myHist) {
		
		int myModeNumber = 0;
		int myMode = 0;
		for (int i = 1 ; i < myHist.length ; i++) if (myHist[i] > myModeNumber) {
			myModeNumber = myHist[i];
			myMode = i;
		}
		
		return myMode;
	}

	public int findMinFromHistogram(int[] myHist) {
		
		int myMinimum = 0;
		for (int i = 0 ; i < myHist.length ; i++) if (myHist[i] > 0) {
			myMinimum = i;
			break;
		}
		
		return myMinimum;
	}

	public float findMeanFromHistogram(int[] myHist) {
		
		float myMean = 0;
		float myWeightedSum = 0;
		float myNumberOfEntries = 0;
		for (int i = 0 ; i < myHist.length ; i++) if (myHist[i] > 0) {
			myNumberOfEntries += myHist[i];
			myWeightedSum += i * myHist[i];
		}
		
		myMean = myWeightedSum / myNumberOfEntries;
		
		return myMean;
	}

	public double findStdFromHistogram(int[] myHist, float myMean) {
				
		float myWeightedSDSum = 0;
		float myNumberOfEntries = 0;
		for (int i = 0 ; i < myHist.length ; i++) if (myHist[i] > 0) {
			myNumberOfEntries += myHist[i];
			myWeightedSDSum += (i - myMean) * (i - myMean) * myHist[i];
		}
		
		double myStd = Math.sqrt(myWeightedSDSum / (myNumberOfEntries - 1));
		
		return myStd;
	}
	
	public int findMaxFromHistogram(int[] myHist) {
				
		int myMaximum = 0;
		for (int i = 0 ; i < myHist.length ; i++) if (myHist[i] > 0) myMaximum = i;
		
		return myMaximum;
	}

	public int findMedianFromHistogram(int[] myHist) {
		
		double[] cumHist = calcCumulativeHistogram(myHist);
		int myMedian = findPercentileFromCumHist(cumHist, 0.5);
		
		return myMedian;
	}

	public double[] calcCumulativeHistogram(int[] myHist) {
		
		double[] cumHist =  new double[myHist.length];
	
		cumHist[0]=0;
		for (int i = 1 ; i < myHist.length ; i++) {
			cumHist[i] = cumHist[i-1] + myHist[i];
		}	
		
		return cumHist;
	}

	public int[] sampleHistogram(ImagePlus nowTiff) {
	
		int numberOfGreyValues = 0;
		switch (nowTiff.getBitDepth()){
			case 8 : numberOfGreyValues = 256;break;
			case 16 : numberOfGreyValues = 65536;break;
			case 32 : {
				numberOfGreyValues = 65536;
				ImageConverter.setDoScaling(true);
				ImageConverter sIC = new ImageConverter(nowTiff);		
				sIC.convertToGray16();
				break;
			}
		}	
		
		//nowTiff.show();
		
		int[] myHist = new int[numberOfGreyValues];
		ImageProcessor outIP = null;
		
		//sample histogram
		for (int i = 1 ; i < nowTiff.getStackSize() + 1 ; i++) {  
			nowTiff.setPosition(i);			
			outIP = nowTiff.getProcessor();
			int[] newHist=outIP.getHistogram();	
			for (int j = 1 ; j < myHist.length ; j++) {
				myHist[j] = myHist[j] + newHist[j];
			}			
		}
		
		return myHist;
	}
	
	public double[] convert16to8BitHistogram(double[] my16Hist) {
		
		double[] my8Hist = new double[256];
		
		int cc = 0;
		for (int i = 0 ; i < my16Hist.length ; i += 256) {			
			
			double mysum = 0;
			
			for (int j = 0 ; j < 256 ; j++) {
				mysum += my16Hist[i + j];
			}	
			
			my8Hist[cc] = mysum / 256; 
			cc++;
		}
		
		return my8Hist;
		
	}
	
	public int[] convert16to8BitHistogram(int[] my16Hist, int lowTh, int highTh) {
		
		int[] my8Hist = new int[256];
		
		double reCalcFac = (double)(highTh - lowTh) / 256;
		
		for (int i = 0 ; i < 256 ; i++) {
			
			double newValue = 0;
			
			double start = i * reCalcFac;
			double stop = (i + 1) * reCalcFac;
							
			int startj = (int)Math.floor(start);
			double startFrac = startj + 1 - start;
			
			int stopj = (int)Math.floor(stop);
			double stopFrac = stop - stopj;
			
			//add fractions of cut classes
			newValue += my16Hist[startj] * startFrac;
			newValue += my16Hist[stopj] * stopFrac;
			
			//add full classes
			for (int j = startj + 1 ; j < stopj ; j++) newValue += my16Hist[j];
			
			my8Hist[i] = (int)Math.round(newValue);
			
		}
		
		return my8Hist;
		
	}

	public int[] sampleThisHistogram16(ImagePlus nowTiff, String myName, InputOutput.MyFileCollection mFC, MenuWaiter.HistogramMenuReturn hMR) {
	
		InputOutput jIO = new InputOutput();	
		RoiHandler roi = new RoiHandler();
							
		int[] myHist = new int[256 * 256];
		int[] newHist;	
		int i, j;
				
		//read gauge file if desired
		PolygonRoi[] pRoi = new PolygonRoi[nowTiff.getNSlices()];		
		if (hMR.useInnerCircle) {
			
			//read inner circle
			ObjectDetector jOD = new ObjectDetector();
			if (mFC.nowInnerCirclePath.contains("Steel")) {
				
				ObjectDetector.EggShapedColCoords3D	jCO = jOD.new EggShapedColCoords3D();
				jCO = jIO.readInnerCircleSteel(mFC);
				pRoi = roi.makeMeAPolygonRoiStack("inner", "exact", jCO, 0);	
			}
			else {
				
				ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();				
				jCO = jIO.readInnerCircleVer1(mFC.nowInnerCirclePath);				
				pRoi = roi.makeMeAPolygonRoiStack("inner", "exact", jCO, 0);
			}
			
			
		}
		
		//get the stacked histogram		
		for (i = 1 ; i < nowTiff.getNSlices() + 1 ; i++) {
			
			IJ.showStatus("Getting 16-bit histogram of slice #" + i + "/" + (nowTiff.getNSlices()));
			
			nowTiff.setPosition(i);		
			
			ImageProcessor myIP = nowTiff.getProcessor();
			ImageProcessor modIP = myIP.duplicate();
			
			//cut out everything outside column
			if (hMR.useInnerCircle) {
				modIP.setRoi(pRoi[i - 1]);
				modIP.setColor(0);
				modIP.fillOutside(pRoi[i - 1]);				
			}
			
			newHist=modIP.getHistogram();
			newHist[0] = 0; //set zero GV to zero
			for (j = 0 ; j < myHist.length ; j++) {
				myHist[j] = myHist[j] + newHist[j];
			}			
			
		}		
		
		return myHist;
		
	}

	public int findTheKnee(int[] myHist) {
		
		int myKnee = 0;
		
		//find maximum and last entry
		double maxVal = 0;
		int myMax = 0;		
		int myLast = 0;
		for (int i = 0 ; i < myHist.length ; i++) {
			if (myHist[i] > maxVal) {
				maxVal = myHist[i];
				myMax = i;
			}
			if (myHist[i] > 0) myLast = i;
		}
		
		//find the coefficients of the spanned string
		double ssA = -maxVal / (myLast - myMax);
		double ssAlpha = - Math.atan(ssA);
		double ssB = Math.tan(ssAlpha) * myLast;
	
		//prepare plot
		double[] xAxis = new double[myHist.length];
		double[] doubleHist = new double[myHist.length];
		double[] spannedString = new double[myHist.length];
		for (int x = 0 ; x < myLast ; x++) {
			xAxis[x] = x;
			doubleHist[x] = myHist[x];
			spannedString[x] = ssA * x + ssB; 
		}
		
		//find line perpendicular to spanned string
		double pBeta = Math.PI / 2 - ssAlpha;
		double pA = Math.tan(pBeta);
		
		//find the longest distance between myHist and spanned string
		double theLongest = 0;		
		for (int x = myMax + 1; x < myLast ; x++){
			double s0 = myHist[x];
			double xNew = (ssB + pA * x - s0) / (pA - ssA);
			double x1 = xNew - x;
			double y = ssA * xNew + ssB;
			double y1 = y - s0; 
			double l =  Math.sqrt(x1*x1 + y1*y1);
			if (l > theLongest) {
				theLongest = l;
				myKnee = x;
			}
			
			//plot it for verification
			double[] lx = new double[(int)Math.round(x1)];
			double[] ly = new double[(int)Math.round(x1)];
			for (int xx = x ; xx < (int)Math.round(xNew) ; xx++) {
				lx[xx - x] = xx;
				ly[xx - x] = s0 + pA * (xx - x); 
			}
			Plot myPlot = new Plot("verification", "greyscale","frequency", xAxis, doubleHist, 4);			
			myPlot.draw();
			myPlot.addPoints(xAxis, spannedString,2);
			myPlot.draw();
			myPlot.addPoints(lx, ly,1);
			//myPlot.setLimits(x - 5, x1 + 5, 0, maxVal);
			myPlot.draw();			
			PlotWindow nowWindow = myPlot.show();
			nowWindow.close();
	
		}
		
		return myKnee;
	}
	
}