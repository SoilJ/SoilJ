package SoilJ.tools;

import org.apache.commons.math3.stat.StatUtils;

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
	
	public class Thresholds {
		
		public int otsu;
		public int minimumJohns;
		public int isodata;
		public int renyiEntropy;
		public int huang;
		public int triangle;
		public int maxEntropy;
		public int minimumIJ;
		public int minError;
		public int li;
		public int yen;		
		
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
		int myMinValue = 100000000;
		for (int i = 0 ; i < myHist.length ; i++) if (myHist[i] > 0 & myHist[i] < myMinValue) {
			myMinimum = i;
			myMinValue = myHist[i];
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
	
	public InputOutput.Histograms calcJointHistogram(float[][] histo) {
		
		InputOutput jIO = new InputOutput();
		
		float[] jHisto = new float[histo[0].length];
		float[] jHistoSD = new float[histo[0].length];

		for (int j = 0 ; j < jHisto.length ; j++) {
			double[] column = new double[histo.length];
			for (int i = 0 ; i < histo.length ; i++) {
				column[i] = histo[i][j];
			}
			jHisto[j] = (float)StatUtils.mean(column);					
			jHistoSD[j] = (float)Math.sqrt(StatUtils.variance(column));  //calculate standard deviation for a later version
		}
		
		InputOutput.Histograms outHist = jIO.new Histograms();
		outHist.floatHistograms = histo;
		outHist.jHisto = jHisto;
		outHist.jHistoSD = jHistoSD;
		outHist.histogramClasses = histo[0].length;
		outHist.numberOfHistograms = histo.length;
		
		return outHist;
		
	}
	
	public int findHistogramMinimum(float[] histo) {
		
		//finds the minimum in between two histogram peaks.. only works for bimodal soil histograms 
		
		RollerCaster rC = new RollerCaster();
		
		//find location of data on histogram
		double xLeft = findQuantileFromHistogram(rC.castFloat2Int(histo), 0.001);

		//find "matrix" peak
		int myMode = findModeFromHistogram(rC.castFloat2Int(histo));
		
		//find secondary peak (mostly corresponding to air)
		float[] leftHisto = new float[(myMode - (int)xLeft) / 2 + (int)xLeft];
		for (int i = 0 ; i < leftHisto.length ; i++) leftHisto[i] = histo[i];
		int myOtherMode = findModeFromHistogram(rC.castFloat2Int(leftHisto));
		
		//check whether secondary mode is not part of matrix peak.. 
		if (leftHisto.length - myOtherMode > 16) {
			leftHisto = new float[(myMode - (int)xLeft) / 3 + (int)xLeft];
			for (int i = 0 ; i < leftHisto.length ; i++) leftHisto[i] = histo[i];
			myOtherMode = findModeFromHistogram(rC.castFloat2Int(leftHisto));
		}

		//find minimum
		float[] myMini = new float[myMode - myOtherMode];
		for (int i = 0 ; i < myMini.length ; i++) myMini[i] = histo[i + myOtherMode];
		int nowMini = findMinFromHistogram(rC.castFloat2Int(myMini));	
		
		return nowMini + myOtherMode;
		
	}
	
	public HistogramStuff.Thresholds calcThreshes(float[] histo) {
		
		RollerCaster rC = new RollerCaster();
		
		Thresholds myThreshes = new Thresholds();
		
		int[] hist = rC.castFloat2Int(histo);

		//find location of data on histogram
		double xLeft = findQuantileFromHistogram(hist, 0.0005);
		double xRight = findQuantileFromHistogram(hist, 0.9995);
		
		//convert to 8 bit
		int[] hist8 = convert16to8BitHistogram(hist, 0, 65335);
		
		//find Minimum using John's algorithm
		myThreshes.minimumJohns = findHistogramMinimum(histo);
		
		//find IJ thresholds
		int bitMultiplier = 1;
		if (histo.length > 1000) bitMultiplier = 256;  //check if histogram is 16 bit 
		myThreshes.otsu = Otsu(hist8) * bitMultiplier;
		myThreshes.maxEntropy = MaxEntropy(hist8) * bitMultiplier;
		myThreshes.isodata = IsoData(hist8) * bitMultiplier;
		myThreshes.renyiEntropy = RenyiEntropy(hist8) * bitMultiplier;
		myThreshes.minimumIJ = Minimum(hist8) * bitMultiplier;
		myThreshes.minError = MinError(hist8) * bitMultiplier;
		myThreshes.li = Li(hist8) * bitMultiplier;
		myThreshes.yen = Yen(hist8) * bitMultiplier;
		myThreshes.huang = Huang(hist8) * bitMultiplier;
		myThreshes.triangle = Triangle(hist8) * bitMultiplier;
				
		return myThreshes;
		
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// the following thresholding routines are copied from AutoThresholder.java in ImageJ
	// input histograms need to be 16-bit
	///////////////////////////////////////////////////////////////////////////////
	
	public int Huang(int[] data) {
	        // Implements Huang's fuzzy thresholding method 
	        // Uses Shannon's entropy function (one can also use Yager's entropy function) 
	        // Huang L.-K. and Wang M.-J.J. (1995) "Image Thresholding by Minimizing  
	        // the Measures of Fuzziness" Pattern Recognition, 28(1): 41-51
	        // M. Emre Celebi  06.15.2007
	        // Ported to ImageJ plugin by G. Landini from E Celebi's fourier_0.8 routines
	        int threshold=-1;
	        int ih, it;
	        int first_bin;
	        int last_bin;
	        double sum_pix;
	        double num_pix;
	        double term;
	        double ent;  // entropy 
	        double min_ent; // min entropy 
	        double mu_x;

	        /* Determine the first non-zero bin */
	        first_bin=0;
	        for (ih = 0; ih < 256; ih++ ) {
	            if ( data[ih] != 0 ) {
	                first_bin = ih;
	                break;
	            }
	        }

	        /* Determine the last non-zero bin */
	        last_bin=255;
	        for (ih = 255; ih >= first_bin; ih-- ) {
	            if ( data[ih] != 0 ) {
	                last_bin = ih;
	                break;
	            }
	        }
	        term = 1.0 / ( double ) ( last_bin - first_bin );
	        double [] mu_0 = new double[256];
	        sum_pix = num_pix = 0;
	        for ( ih = first_bin; ih < 256; ih++ ){
	            sum_pix += (double)ih * data[ih];
	            num_pix += data[ih];
	            /* NUM_PIX cannot be zero ! */
	            mu_0[ih] = sum_pix / num_pix;
	        }

	        double [] mu_1 = new double[256];
	        sum_pix = num_pix = 0;
	        for ( ih = last_bin; ih > 0; ih-- ){
	            sum_pix += (double)ih * data[ih];
	            num_pix += data[ih];
	            /* NUM_PIX cannot be zero ! */
	            mu_1[ih - 1] = sum_pix / ( double ) num_pix;
	        }

	        /* Determine the threshold that minimizes the fuzzy entropy */
	        threshold = -1;
	        min_ent = Double.MAX_VALUE;
	        for ( it = 0; it < 256; it++ ){
	            ent = 0.0;
	            for ( ih = 0; ih <= it; ih++ ) {
	                /* Equation (4) in Ref. 1 */
	                mu_x = 1.0 / ( 1.0 + term * Math.abs ( ih - mu_0[it] ) );
	                if ( !((mu_x  < 1e-06 ) || ( mu_x > 0.999999))) {
	                    /* Equation (6) & (8) in Ref. 1 */
	                    ent += data[ih] * ( -mu_x * Math.log ( mu_x ) - ( 1.0 - mu_x ) * Math.log ( 1.0 - mu_x ) );
	                }
	            }

	            for ( ih = it + 1; ih < 256; ih++ ) {
	                /* Equation (4) in Ref. 1 */
	                mu_x = 1.0 / ( 1.0 + term * Math.abs ( ih - mu_1[it] ) );
	                if ( !((mu_x  < 1e-06 ) || ( mu_x > 0.999999))) {
	                    /* Equation (6) & (8) in Ref. 1 */
	                    ent += data[ih] * ( -mu_x * Math.log ( mu_x ) - ( 1.0 - mu_x ) * Math.log ( 1.0 - mu_x ) );
	                }
	            }
	            /* No need to divide by NUM_ROWS * NUM_COLS * LOG(2) ! */
	            if ( ent < min_ent ) {
	                min_ent = ent;
	                threshold = it;
	            }
	        }
	        return threshold;
	    }

	    public boolean bimodalTest(double [] y) {
	        int len=y.length;
	        boolean b = false;
	        int modes = 0;
	 
	        for (int k=1;k<len-1;k++){
	            if (y[k-1] < y[k] && y[k+1] < y[k]) {
	                modes++;
	                if (modes>2)
	                    return false;
	            }
	        }
	        if (modes == 2)
	            b = true;
	        return b;
	    }

	    public int Intermodes(int[] data) {
	        // J. M. S. Prewitt and M. L. Mendelsohn, "The analysis of cell images," in
	        // Annals of the New York Academy of Sciences, vol. 128, pp. 1035-1053, 1966.
	        // ported to ImageJ plugin by G.Landini from Antti Niemisto's Matlab code (GPL)
	        // Original Matlab code Copyright (C) 2004 Antti Niemisto
	        // See http://www.cs.tut.fi/~ant/histthresh/ for an excellent slide presentation
	        // and the original Matlab code.
	        //
	        // Assumes a bimodal histogram. The histogram needs is smoothed (using a
	        // running average of size 3, iteratively) until there are only two local maxima.
	        // j and k
	        // Threshold t is (j+k)/2.
	        // Images with histograms having extremely unequal peaks or a broad and
	        // flat valleys are unsuitable for this method.
	        
	        int minbin=-1, maxbin=-1;
	        for (int i=0; i<data.length; i++)
	            if (data[i]>0) maxbin = i;
	        for (int i=data.length-1; i>=0; i--)
	            if (data[i]>0) minbin = i;
	        int length = (maxbin-minbin)+1;
	        double [] hist = new double[length];
	        for (int i=minbin; i<=maxbin; i++)
	            hist[i-minbin] = data[i];
	            
	        int iter = 0;
	        int threshold=-1;
	        while (!bimodalTest(hist) ) {
	             //smooth with a 3 point running mean filter
	            double previous=0, current=0, next=hist[0];
	            for (int i=0; i<length-1; i++) {
	                previous = current;
	                current = next;
	                next = hist[i + 1];
	                hist[i] = (previous+current+next)/3;
	            }
	            hist[length-1] = (current+next)/3;
	            iter++;
	            if (iter>10000) {
	                threshold = -1;
	                IJ.log("Intermodes Threshold not found after 10000 iterations.");
	                return threshold;
	            }
	        }

	        // The threshold is the mean between the two peaks.
	        int tt=0;
	        for (int i=1; i<length - 1; i++) {
	            if (hist[i-1] < hist[i] && hist[i+1] < hist[i]){
	                tt += i;
	                //IJ.log("mode:" +i);
	            }
	        }
	        threshold = (int) Math.floor(tt/2.0);
	        return threshold+minbin;
	    }

	    public int IsoData(int[] data) {
	        // Also called intermeans
	        // Iterative procedure based on the isodata algorithm [T.W. Ridler, S. Calvard, Picture 
	        // thresholding using an iterative selection method, IEEE Trans. System, Man and 
	        // Cybernetics, SMC-8 (1978) 630-632.] 
	        // The procedure divides the image into objects and background by taking an initial threshold,
	        // then the averages of the pixels at or below the threshold and pixels above are computed. 
	        // The averages of those two values are computed, the threshold is incremented and the 
	        // process is repeated until the threshold is larger than the composite average. That is,
	        //  threshold = (average background + average objects)/2
	        // The code in ImageJ that implements this function is the getAutoThreshold() method in the ImageProcessor class. 
	        //
	        // From: Tim Morris (dtm@ap.co.umist.ac.uk)
	        // Subject: Re: Thresholding method?
	        // posted to sci.image.processing on 1996/06/24
	        // The algorithm implemented in NIH Image sets the threshold as that grey
	        // value, G, for which the average of the averages of the grey values
	        // below and above G is equal to G. It does this by initialising G to the
	        // lowest sensible value and iterating:

	        // L = the average grey value of pixels with intensities < G
	        // H = the average grey value of pixels with intensities > G
	        // is G = (L + H)/2?
	        // yes => exit
	        // no => increment G and repeat
	        //
	        int i, l, totl, g=0;
	        double toth, h;
	        for (i = 1; i < 256; i++) {
	            if (data[i] > 0){
	                g = i + 1;
	                break;
	            }
	        }
	        while (true){
	            l = 0;
	            totl = 0;
	            for (i = 0; i < g; i++) {
	                 totl = totl + data[i];
	                 l = l + (data[i] * i);
	            }
	            h = 0;
	            toth = 0;
	            for (i = g + 1; i < 256; i++){
	                toth += data[i];
	                h += ((double)data[i]*i);
	            }
	            if (totl > 0 && toth > 0){
	                l /= totl;
	                h /= toth;
	                if (g == (int) Math.round((l + h) / 2.0))
	                    break;
	            }
	            g++;
	            if (g > 254)
	                return -1;
	        }
	        return g;
	    }
	    
	    public int defaultIsoData(int[] data) {
	        // This is the modified IsoData method used by the "Threshold" widget in "Default" mode
	        int n = data.length;
	        int[] data2 = new int[n];
	        int mode=0, maxCount=0;
	        for (int i=0; i<n; i++) {
	            int count = data[i];
	            data2[i] = data[i];
	            if (data2[i]>maxCount) {
	                maxCount = data2[i];
	                mode = i;
	            }
	        }
	        int maxCount2 = 0;
	        for (int i = 0; i<n; i++) {
	            if ((data2[i]>maxCount2) && (i!=mode))
	                maxCount2 = data2[i];
	        }
	        int hmax = maxCount;
	        if ((hmax>(maxCount2*2)) && (maxCount2!=0)) {
	            hmax = (int)(maxCount2 * 1.5);
	            data2[mode] = hmax;
	        }
	        return IJIsoData(data2);
	    }

	    public int IJIsoData(int[] data) {
	        // This is the original ImageJ IsoData implementation, here for backward compatibility.
	        int level;
	        int maxValue = data.length - 1;
	        double result, sum1, sum2, sum3, sum4;
	        int count0 = data[0];
	        data[0] = 0; //set to zero so erased areas aren't included
	        int countMax = data[maxValue];
	        data[maxValue] = 0;
	        int min = 0;
	        while ((data[min]==0) && (min<maxValue))
	            min++;
	        int max = maxValue;
	        while ((data[max]==0) && (max>0))
	            max--;
	        if (min>=max) {
	            data[0]= count0; data[maxValue]=countMax;
	            level = data.length/2;
	            return level;
	        }
	        int movingIndex = min;
	        int inc = Math.max(max/40, 1);
	        do {
	            sum1=sum2=sum3=sum4=0.0;
	            for (int i=min; i<=movingIndex; i++) {
	                sum1 += (double)i*data[i];
	                sum2 += data[i];
	            }
	            for (int i=(movingIndex+1); i<=max; i++) {
	                sum3 += (double)i*data[i];
	                sum4 += data[i];
	            }           
	            result = (sum1/sum2 + sum3/sum4)/2.0;
	            movingIndex++;
	        } while ((movingIndex+1)<=result && movingIndex<max-1);
	        data[0]= count0; data[maxValue]=countMax;
	        level = (int)Math.round(result);
	        return level;
	    }


	    public int Li(int[] data) {
	        // Implements Li's Minimum Cross Entropy thresholding method
	        // This implementation is based on the iterative version (Ref. 2) of the algorithm.
	        // 1) Li C.H. and Lee C.K. (1993) "Minimum Cross Entropy Thresholding" 
	        //    Pattern Recognition, 26(4): 617-625
	        // 2) Li C.H. and Tam P.K.S. (1998) "An Iterative Algorithm for Minimum 
	        //    Cross Entropy Thresholding"Pattern Recognition Letters, 18(8): 771-776
	        // 3) Sezgin M. and Sankur B. (2004) "Survey over Image Thresholding 
	        //    Techniques and Quantitative Performance Evaluation" Journal of 
	        //    Electronic Imaging, 13(1): 146-165 
	        //    http://citeseer.ist.psu.edu/sezgin04survey.html
	        // Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines
	        int threshold;
	        double num_pixels;
	        double sum_back; /* sum of the background pixels at a given threshold */
	        double sum_obj;  /* sum of the object pixels at a given threshold */
	        double num_back; /* number of background pixels at a given threshold */
	        double num_obj;  /* number of object pixels at a given threshold */
	        double old_thresh;
	        double new_thresh;
	        double mean_back; /* mean of the background pixels at a given threshold */
	        double mean_obj;  /* mean of the object pixels at a given threshold */
	        double mean;  /* mean gray-level in the image */
	        double tolerance; /* threshold tolerance */
	        double temp;

	        tolerance=0.5;
	        num_pixels = 0;
	        for (int ih = 0; ih < 256; ih++ ) 
	            num_pixels += data[ih];

	        /* Calculate the mean gray-level */
	        mean = 0.0;
	        for (int ih = 0 + 1; ih < 256; ih++ ) //0 + 1?
	            mean += (double)ih * data[ih];
	        mean /= num_pixels;
	        /* Initial estimate */
	        new_thresh = mean;

	        do {
	            old_thresh = new_thresh;
	            threshold = (int) (old_thresh + 0.5);   /* range */
	            /* Calculate the means of background and object pixels */
	            /* Background */
	            sum_back = 0;
	            num_back = 0;
	            for (int ih = 0; ih <= threshold; ih++ ) {
	                sum_back += (double)ih * data[ih];
	                num_back += data[ih];
	            }
	            mean_back = ( num_back == 0 ? 0.0 : ( sum_back / ( double ) num_back ) );
	            /* Object */
	            sum_obj = 0;
	            num_obj = 0;
	            for (int ih = threshold + 1; ih < 256; ih++ ) {
	                sum_obj += (double)ih * data[ih];
	                num_obj += data[ih];
	            }
	            mean_obj = ( num_obj == 0 ? 0.0 : ( sum_obj / ( double ) num_obj ) );

	            /* Calculate the new threshold: Equation (7) in Ref. 2 */
	            //new_thresh = simple_round ( ( mean_back - mean_obj ) / ( Math.log ( mean_back ) - Math.log ( mean_obj ) ) );
	            //simple_round ( double x ) {
	            // return ( int ) ( IS_NEG ( x ) ? x - .5 : x + .5 );
	            //}
	            //
	            //#define IS_NEG( x ) ( ( x ) < -DBL_EPSILON ) 
	            //DBL_EPSILON = 2.220446049250313E-16
	            temp = ( mean_back - mean_obj ) / ( Math.log ( mean_back ) - Math.log ( mean_obj ) );

	            if (temp < -2.220446049250313E-16)
	                new_thresh = (int) (temp - 0.5);
	            else
	                new_thresh = (int) (temp + 0.5);
	            /*  Stop the iterations when the difference between the
	            new and old threshold values is less than the tolerance */
	        }
	        while ( Math.abs ( new_thresh - old_thresh ) > tolerance );
	        return threshold;
	    }

	    public int MaxEntropy(int[] data) {
	        // Implements Kapur-Sahoo-Wong (Maximum Entropy) thresholding method
	        // Kapur J.N., Sahoo P.K., and Wong A.K.C. (1985) "A New Method for
	        // Gray-Level Picture Thresholding Using the Entropy of the Histogram"
	        // Graphical Models and Image Processing, 29(3): 273-285
	        // M. Emre Celebi
	        // 06.15.2007
	        // Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines
	        int threshold=-1;
	        int ih, it;
	        int first_bin;
	        int last_bin;
	        double tot_ent;  /* total entropy */
	        double max_ent;  /* max entropy */
	        double ent_back; /* entropy of the background pixels at a given threshold */
	        double ent_obj;  /* entropy of the object pixels at a given threshold */
	        double [] norm_histo = new double[256]; /* normalized histogram */
	        double [] P1 = new double[256]; /* cumulative normalized histogram */
	        double [] P2 = new double[256]; 

	        double total =0;
	        for (ih = 0; ih < 256; ih++ ) 
	            total+=data[ih];

	        for (ih = 0; ih < 256; ih++ )
	            norm_histo[ih] = data[ih]/total;

	        P1[0]=norm_histo[0];
	        P2[0]=1.0-P1[0];
	        for (ih = 1; ih < 256; ih++ ){
	            P1[ih]= P1[ih-1] + norm_histo[ih];
	            P2[ih]= 1.0 - P1[ih];
	        }

	        /* Determine the first non-zero bin */
	        first_bin=0;
	        for (ih = 0; ih < 256; ih++ ) {
	            if ( !(Math.abs(P1[ih])<2.220446049250313E-16)) {
	                first_bin = ih;
	                break;
	            }
	        }

	        /* Determine the last non-zero bin */
	        last_bin=255;
	        for (ih = 255; ih >= first_bin; ih-- ) {
	            if ( !(Math.abs(P2[ih])<2.220446049250313E-16)) {
	                last_bin = ih;
	                break;
	            }
	        }

	        // Calculate the total entropy each gray-level
	        // and find the threshold that maximizes it 
	        max_ent = Double.MIN_VALUE;

	        for ( it = first_bin; it <= last_bin; it++ ) {
	            /* Entropy of the background pixels */
	            ent_back = 0.0;
	            for ( ih = 0; ih <= it; ih++ )  {
	                if ( data[ih] !=0 ) {
	                    ent_back -= ( norm_histo[ih] / P1[it] ) * Math.log ( norm_histo[ih] / P1[it] );
	                }
	            }

	            /* Entropy of the object pixels */
	            ent_obj = 0.0;
	            for ( ih = it + 1; ih < 256; ih++ ){
	                if (data[ih]!=0){
	                ent_obj -= ( norm_histo[ih] / P2[it] ) * Math.log ( norm_histo[ih] / P2[it] );
	                }
	            }

	            /* Total entropy */
	            tot_ent = ent_back + ent_obj;

	            // IJ.log(""+max_ent+"  "+tot_ent);
	            if ( max_ent < tot_ent ) {
	                max_ent = tot_ent;
	                threshold = it;
	            }
	        }
	        return threshold;
	    }

	    public int Mean(int[] data) {
	        // C. A. Glasbey, "An analysis of histogram-based thresholding algorithms,"
	        // CVGIP: Graphical Models and Image Processing, vol. 55, pp. 532-537, 1993.
	        //
	        // The threshold is the mean of the greyscale data
	        int threshold = -1;
	        double tot=0, sum=0;
	        for (int i=0; i<256; i++){
	            tot+= data[i];
	            sum+=((double)i*data[i]);
	        }
	        threshold =(int) Math.floor(sum/tot);
	        return threshold;
	    }

	    public int MinError(int[] data) {
	        // Kittler and J. Illingworth, "Minimum error thresholding," Pattern Recognition, vol. 19, pp. 41-47, 1986.
	        // C. A. Glasbey, "An analysis of histogram-based thresholding algorithms," CVGIP: Graphical Models and Image Processing, vol. 55, pp. 532-537, 1993.
	        // Ported to ImageJ plugin by G.Landini from Antti Niemisto's Matlab code (GPL)
	        // Original Matlab code Copyright (C) 2004 Antti Niemisto
	        // See http://www.cs.tut.fi/~ant/histthresh/ for an excellent slide presentation
	        // and the original Matlab code.

	        int threshold = Mean(data); //Initial estimate for the threshold is found with the MEAN algorithm.
	        int Tprev =-2;
	        double mu, nu, p, q, sigma2, tau2, w0, w1, w2, sqterm, temp;
	        //int counter=1;
	        while (threshold!=Tprev){
	            //Calculate some statistics.
	            mu = B(data, threshold)/A(data, threshold);
	            nu = (B(data, data.length - 1)-B(data, threshold))/(A(data, data.length - 1)-A(data, threshold));
	            p = A(data, threshold)/A(data, data.length - 1);
	            q = (A(data, data.length - 1)-A(data, threshold)) / A(data, data.length - 1);
	            sigma2 = C(data, threshold)/A(data, threshold)-(mu*mu);
	            tau2 = (C(data, data.length - 1)-C(data, threshold)) / (A(data, data.length - 1)-A(data, threshold)) - (nu*nu);

	            //The terms of the quadratic equation to be solved.
	            w0 = 1.0/sigma2-1.0/tau2;
	            w1 = mu/sigma2-nu/tau2;
	            w2 = (mu*mu)/sigma2 - (nu*nu)/tau2 + Math.log10((sigma2*(q*q))/(tau2*(p*p)));

	            //If the next threshold would be imaginary, return with the current one.
	            sqterm = (w1*w1)-w0*w2;
	            if (sqterm < 0) {
	                IJ.log("MinError(I): not converging.");
	                return threshold;
	            }

	            //The updated threshold is the integer part of the solution of the quadratic equation.
	            Tprev = threshold;
	            temp = (w1+Math.sqrt(sqterm))/w0;

	            if (Double.isNaN(temp))
	                threshold = Tprev;
	            else
	                threshold =(int) Math.floor(temp);
	        }
	        return threshold;
	    }

	    private double A(int[] y, int j) {
	        if (j>=y.length) j=y.length-1;
	        double x = 0;
	        for (int i=0;i<=j;i++)
	            x+=y[i];
	        return x;
	    }

	    private double B(int[] y, int j) {
	        if (j>=y.length) j=y.length-1;
	        double x = 0;
	        for (int i=0;i<=j;i++)
	            x+=i*y[i];
	        return x;
	    }

	    private double C(int[] y, int j) {
	        if (j>=y.length) j=y.length-1;
	        double x = 0;
	        for (int i=0;i<=j;i++)
	            x+=i*i*y[i];
	        return x;
	    }
	    
	    public int Minimum(int[] data) {
	        // J. M. S. Prewitt and M. L. Mendelsohn, "The analysis of cell images," in
	        // Annals of the New York Academy of Sciences, vol. 128, pp. 1035-1053, 1966.
	        // ported to ImageJ plugin by G.Landini from Antti Niemisto's Matlab code (GPL)
	        // Original Matlab code Copyright (C) 2004 Antti Niemisto
	        // See http://www.cs.tut.fi/~ant/histthresh/ for an excellent slide presentation
	        // and the original Matlab code.
	        //
	        // Assumes a bimodal histogram. The histogram needs is smoothed (using a
	        // running average of size 3, iteratively) until there are only two local maxima.
	        // Threshold t is such that yt-1 > yt <= yt+1.
	        // Images with histograms having extremely unequal peaks or a broad and
	        // flat valleys are unsuitable for this method.
	        int iter =0;
	        int threshold = -1;
	        double [] iHisto = new double [256];
	        for (int i=0; i<256; i++)
	            iHisto[i]=(double) data[i];
	        double [] tHisto = new double[iHisto.length] ;

	        while (!bimodalTest(iHisto) ) {
	             //smooth with a 3 point running mean filter
	            for (int i=1; i<255; i++)
	                tHisto[i]= (iHisto[i-1] + iHisto[i] +iHisto[i+1])/3;
	            tHisto[0] = (iHisto[0]+iHisto[1])/3; //0 outside
	            tHisto[255] = (iHisto[254]+iHisto[255])/3; //0 outside
	            System.arraycopy(tHisto, 0, iHisto, 0, iHisto.length) ;
	            iter++;
	            if (iter>10000) {
	                threshold = -1;
	                IJ.log("Minimum: threshold not found after 10000 iterations.");
	                return threshold;
	            }
	        }
	        // The threshold is the minimum between the two peaks.
	        for (int i=1; i<255; i++) {
	            if (iHisto[i-1] > iHisto[i] && iHisto[i+1] >= iHisto[i]) {
	                threshold = i;
	                break;
	            }
	        }
	        return threshold;
	    }

	    public int Moments(int[] data) {
	        //  W. Tsai, "Moment-preserving thresholding: a new approach," Computer Vision,
	        // Graphics, and Image Processing, vol. 29, pp. 377-393, 1985.
	        // Ported to ImageJ plugin by G.Landini from the the open source project FOURIER 0.8
	        // by  M. Emre Celebi , Department of Computer Science,  Louisiana State University in Shreveport
	        // Shreveport, LA 71115, USA
	        //  http://sourceforge.net/projects/fourier-ipal
	        //  http://www.lsus.edu/faculty/~ecelebi/fourier.htm
	        double total =0;
	        double m0=1.0, m1=0.0, m2 =0.0, m3 =0.0, sum =0.0, p0=0.0;
	        double cd, c0, c1, z0, z1;  /* auxiliary variables */
	        int threshold = -1;

	        double [] histo = new  double [256];

	        for (int i=0; i<256; i++)
	            total+=data[i];

	        for (int i=0; i<256; i++)
	            histo[i]=(double)(data[i]/total); //normalised histogram

	        /* Calculate the first, second, and third order moments */
	        for ( int i = 0; i < 256; i++ ) {
	            double di = i;
	            m1 += di * histo[i];
	            m2 += di * di * histo[i];
	            m3 += di * di * di * histo[i];
	        }
	        /* 
	        First 4 moments of the gray-level image should match the first 4 moments
	        of the target binary image. This leads to 4 equalities whose solutions 
	        are given in the Appendix of Ref. 1 
	        */
	        cd = m0 * m2 - m1 * m1;
	        c0 = ( -m2 * m2 + m1 * m3 ) / cd;
	        c1 = ( m0 * -m3 + m2 * m1 ) / cd;
	        z0 = 0.5 * ( -c1 - Math.sqrt ( c1 * c1 - 4.0 * c0 ) );
	        z1 = 0.5 * ( -c1 + Math.sqrt ( c1 * c1 - 4.0 * c0 ) );
	        p0 = ( z1 - m1 ) / ( z1 - z0 );  /* Fraction of the object pixels in the target binary image */

	        // The threshold is the gray-level closest  
	        // to the p0-tile of the normalized histogram 
	        sum=0;
	        for (int i=0; i<256; i++){
	            sum+=histo[i];
	            if (sum>p0) {
	                threshold = i;
	                break;
	            }
	        }
	        return threshold;
	    }

	    public int Otsu(int[] data) {
	        // Otsu's threshold algorithm
	        // C++ code by Jordan Bevik <Jordan.Bevic@qtiworld.com>
	        // ported to ImageJ plugin by G.Landini
	        int k,kStar;  // k = the current threshold; kStar = optimal threshold
	        double N1, N;    // N1 = # points with intensity <=k; N = total number of points
	        double BCV, BCVmax; // The current Between Class Variance and maximum BCV
	        double num, denom;  // temporary bookeeping
	        double Sk;  // The total intensity for all histogram points <=k
	        double S, L=256; // The total intensity of the image

	        // Initialize values:
	        S = N = 0;
	        for (k=0; k<L; k++){
	            S += (double)k * data[k];   // Total histogram intensity
	            N += data[k];       // Total number of data points
	        }

	        Sk = 0;
	        N1 = data[0]; // The entry for zero intensity
	        BCV = 0;
	        BCVmax=0;
	        kStar = 0;

	        // Look at each possible threshold value,
	        // calculate the between-class variance, and decide if it's a max
	        for (k=1; k<L-1; k++) { // No need to check endpoints k = 0 or k = L-1
	            Sk += (double)k * data[k];
	            N1 += data[k];

	            // The float casting here is to avoid compiler warning about loss of precision and
	            // will prevent overflow in the case of large saturated images
	            denom = (double)( N1) * (N - N1); // Maximum value of denom is (N^2)/4 =  approx. 3E10

	            if (denom != 0 ){
	                // Float here is to avoid loss of precision when dividing
	                num = ( (double)N1 / N ) * S - Sk;  // Maximum value of num =  255*N = approx 8E7
	                BCV = (num * num) / denom;
	            }
	            else
	                BCV = 0;

	            if (BCV >= BCVmax){ // Assign the best threshold found so far
	                BCVmax = BCV;
	                kStar = k;
	            }
	        }
	        // kStar += 1;  // Use QTI convention that intensity -> 1 if intensity >= k
	        // (the algorithm was developed for I-> 1 if I <= k.)
	        return kStar;
	    }


	    public int Percentile(int[] data) {
	        // W. Doyle, "Operation useful for similarity-invariant pattern recognition,"
	        // Journal of the Association for Computing Machinery, vol. 9,pp. 259-267, 1962.
	        // ported to ImageJ plugin by G.Landini from Antti Niemisto's Matlab code (GPL)
	        // Original Matlab code Copyright (C) 2004 Antti Niemisto
	        // See http://www.cs.tut.fi/~ant/histthresh/ for an excellent slide presentation
	        // and the original Matlab code.

	        int iter =0;
	        int threshold = -1;
	        double ptile= 0.5; // default fraction of foreground pixels
	        double [] avec = new double [256];

	        for (int i=0; i<256; i++)
	            avec[i]=0.0;

	        double total =partialSum(data, 255);
	        double temp = 1.0;
	        for (int i=0; i<256; i++){
	            avec[i]=Math.abs((partialSum(data, i)/total)-ptile);
	            //IJ.log("Ptile["+i+"]:"+ avec[i]);
	            if (avec[i]<temp) {
	                temp = avec[i];
	                threshold = i;
	            }
	        }
	        return threshold;
	    }


	    double partialSum(int [] y, int j) {
	        double x = 0;
	        for (int i=0;i<=j;i++)
	            x+=y[i];
	        return x;
	    }


	    public int RenyiEntropy(int[] data) {
	        // Kapur J.N., Sahoo P.K., and Wong A.K.C. (1985) "A New Method for
	        // Gray-Level Picture Thresholding Using the Entropy of the Histogram"
	        // Graphical Models and Image Processing, 29(3): 273-285
	        // M. Emre Celebi
	        // 06.15.2007
	        // Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines

	        int threshold; 
	        int opt_threshold;

	        int ih, it;
	        int first_bin;
	        int last_bin;
	        int tmp_var;
	        int t_star1, t_star2, t_star3;
	        int beta1, beta2, beta3;
	        double alpha;/* alpha parameter of the method */
	        double term;
	        double tot_ent;  /* total entropy */
	        double max_ent;  /* max entropy */
	        double ent_back; /* entropy of the background pixels at a given threshold */
	        double ent_obj;  /* entropy of the object pixels at a given threshold */
	        double omega;
	        double [] norm_histo = new double[256]; /* normalized histogram */
	        double [] P1 = new double[256]; /* cumulative normalized histogram */
	        double [] P2 = new double[256]; 

	        double total =0;
	        for (ih = 0; ih < 256; ih++ ) 
	            total+=data[ih];

	        for (ih = 0; ih < 256; ih++ )
	            norm_histo[ih] = data[ih]/total;

	        P1[0]=norm_histo[0];
	        P2[0]=1.0-P1[0];
	        for (ih = 1; ih < 256; ih++ ){
	            P1[ih]= P1[ih-1] + norm_histo[ih];
	            P2[ih]= 1.0 - P1[ih];
	        }

	        /* Determine the first non-zero bin */
	        first_bin=0;
	        for (ih = 0; ih < 256; ih++ ) {
	            if ( !(Math.abs(P1[ih])<2.220446049250313E-16)) {
	                first_bin = ih;
	                break;
	            }
	        }

	        /* Determine the last non-zero bin */
	        last_bin=255;
	        for (ih = 255; ih >= first_bin; ih-- ) {
	            if ( !(Math.abs(P2[ih])<2.220446049250313E-16)) {
	                last_bin = ih;
	                break;
	            }
	        }

	        /* Maximum Entropy Thresholding - BEGIN */
	        /* ALPHA = 1.0 */
	        /* Calculate the total entropy each gray-level
	        and find the threshold that maximizes it 
	        */
	        threshold =0; // was MIN_INT in original code, but if an empty image is processed it gives an error later on.
	        max_ent = 0.0;

	        for ( it = first_bin; it <= last_bin; it++ ) {
	            /* Entropy of the background pixels */
	            ent_back = 0.0;
	            for ( ih = 0; ih <= it; ih++ )  {
	                if ( data[ih] !=0 ) {
	                    ent_back -= ( norm_histo[ih] / P1[it] ) * Math.log ( norm_histo[ih] / P1[it] );
	                }
	            }

	            /* Entropy of the object pixels */
	            ent_obj = 0.0;
	            for ( ih = it + 1; ih < 256; ih++ ){
	                if (data[ih]!=0){
	                ent_obj -= ( norm_histo[ih] / P2[it] ) * Math.log ( norm_histo[ih] / P2[it] );
	                }
	            }

	            /* Total entropy */
	            tot_ent = ent_back + ent_obj;

	            // IJ.log(""+max_ent+"  "+tot_ent);

	            if ( max_ent < tot_ent ) {
	                max_ent = tot_ent;
	                threshold = it;
	            }
	        }
	        t_star2 = threshold;

	        /* Maximum Entropy Thresholding - END */
	        threshold =0; //was MIN_INT in original code, but if an empty image is processed it gives an error later on.
	        max_ent = 0.0;
	        alpha = 0.5;
	        term = 1.0 / ( 1.0 - alpha );
	        for ( it = first_bin; it <= last_bin; it++ ) {
	            /* Entropy of the background pixels */
	            ent_back = 0.0;
	            for ( ih = 0; ih <= it; ih++ )
	                ent_back += Math.sqrt ( norm_histo[ih] / P1[it] );

	            /* Entropy of the object pixels */
	            ent_obj = 0.0;
	            for ( ih = it + 1; ih < 256; ih++ )
	                ent_obj += Math.sqrt ( norm_histo[ih] / P2[it] );

	            /* Total entropy */
	            tot_ent = term * ( ( ent_back * ent_obj ) > 0.0 ? Math.log ( ent_back * ent_obj ) : 0.0);

	            if ( tot_ent > max_ent ){
	                max_ent = tot_ent;
	                threshold = it;
	            }
	        }

	        t_star1 = threshold;

	        threshold = 0; //was MIN_INT in original code, but if an empty image is processed it gives an error later on.
	        max_ent = 0.0;
	        alpha = 2.0;
	        term = 1.0 / ( 1.0 - alpha );
	        for ( it = first_bin; it <= last_bin; it++ ) {
	            /* Entropy of the background pixels */
	            ent_back = 0.0;
	            for ( ih = 0; ih <= it; ih++ )
	                ent_back += ( norm_histo[ih] * norm_histo[ih] ) / ( P1[it] * P1[it] );

	            /* Entropy of the object pixels */
	            ent_obj = 0.0;
	            for ( ih = it + 1; ih < 256; ih++ )
	                ent_obj += ( norm_histo[ih] * norm_histo[ih] ) / ( P2[it] * P2[it] );

	            /* Total entropy */
	            tot_ent = term *( ( ent_back * ent_obj ) > 0.0 ? Math.log(ent_back * ent_obj ): 0.0 );

	            if ( tot_ent > max_ent ){
	                max_ent = tot_ent;
	                threshold = it;
	            }
	        }

	        t_star3 = threshold;

	        /* Sort t_star values */
	        if ( t_star2 < t_star1 ){
	            tmp_var = t_star1;
	            t_star1 = t_star2;
	            t_star2 = tmp_var;
	        }
	        if ( t_star3 < t_star2 ){
	            tmp_var = t_star2;
	            t_star2 = t_star3;
	            t_star3 = tmp_var;
	        }
	        if ( t_star2 < t_star1 ) {
	            tmp_var = t_star1;
	            t_star1 = t_star2;
	            t_star2 = tmp_var;
	        }

	        /* Adjust beta values */
	        if ( Math.abs ( t_star1 - t_star2 ) <= 5 )  {
	            if ( Math.abs ( t_star2 - t_star3 ) <= 5 ) {
	                beta1 = 1;
	                beta2 = 2;
	                beta3 = 1;
	            }
	            else {
	                beta1 = 0;
	                beta2 = 1;
	                beta3 = 3;
	            }
	        }
	        else {
	            if ( Math.abs ( t_star2 - t_star3 ) <= 5 ) {
	                beta1 = 3;
	                beta2 = 1;
	                beta3 = 0;
	            }
	            else {
	                beta1 = 1;
	                beta2 = 2;
	                beta3 = 1;
	            }
	        }
	        //IJ.log(""+t_star1+" "+t_star2+" "+t_star3);
	        /* Determine the optimal threshold value */
	        omega = P1[t_star3] - P1[t_star1];
	        opt_threshold = (int) (t_star1 * ( P1[t_star1] + 0.25 * omega * beta1 ) + 0.25 * t_star2 * omega * beta2  + t_star3 * ( P2[t_star3] + 0.25 * omega * beta3 ));

	        return opt_threshold;
	    }


	    public int Shanbhag(int[] data) {
	        // Shanhbag A.G. (1994) "Utilization of Information Measure as a Means of
	        //  Image Thresholding" Graphical Models and Image Processing, 56(5): 414-419
	        // Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines
	        int threshold;
	        int ih, it;
	        int first_bin;
	        int last_bin;
	        double term;
	        double tot_ent;  /* total entropy */
	        double min_ent;  /* max entropy */
	        double ent_back; /* entropy of the background pixels at a given threshold */
	        double ent_obj;  /* entropy of the object pixels at a given threshold */
	        double [] norm_histo = new double[256]; /* normalized histogram */
	        double [] P1 = new double[256]; /* cumulative normalized histogram */
	        double [] P2 = new double[256]; 

	        double total =0;
	        for (ih = 0; ih < 256; ih++ ) 
	            total+=data[ih];

	        for (ih = 0; ih < 256; ih++ )
	            norm_histo[ih] = data[ih]/total;

	        P1[0]=norm_histo[0];
	        P2[0]=1.0-P1[0];
	        for (ih = 1; ih < 256; ih++ ){
	            P1[ih]= P1[ih-1] + norm_histo[ih];
	            P2[ih]= 1.0 - P1[ih];
	        }

	        /* Determine the first non-zero bin */
	        first_bin=0;
	        for (ih = 0; ih < 256; ih++ ) {
	            if ( !(Math.abs(P1[ih])<2.220446049250313E-16)) {
	                first_bin = ih;
	                break;
	            }
	        }

	        /* Determine the last non-zero bin */
	        last_bin=255;
	        for (ih = 255; ih >= first_bin; ih-- ) {
	            if ( !(Math.abs(P2[ih])<2.220446049250313E-16)) {
	                last_bin = ih;
	                break;
	            }
	        }

	        // Calculate the total entropy each gray-level
	        // and find the threshold that maximizes it 
	        threshold =-1;
	        min_ent = Double.MAX_VALUE;

	        for ( it = first_bin; it <= last_bin; it++ ) {
	            /* Entropy of the background pixels */
	            ent_back = 0.0;
	            term = 0.5 / P1[it];
	            for ( ih = 1; ih <= it; ih++ )  { //0+1?
	                ent_back -= norm_histo[ih] * Math.log ( 1.0 - term * P1[ih - 1] );
	            }
	            ent_back *= term;

	            /* Entropy of the object pixels */
	            ent_obj = 0.0;
	            term = 0.5 / P2[it];
	            for ( ih = it + 1; ih < 256; ih++ ){
	                ent_obj -= norm_histo[ih] * Math.log ( 1.0 - term * P2[ih] );
	            }
	            ent_obj *= term;

	            /* Total entropy */
	            tot_ent = Math.abs ( ent_back - ent_obj );

	            if ( tot_ent < min_ent ) {
	                min_ent = tot_ent;
	                threshold = it;
	            }
	        }
	        return threshold;
	    }


	    public int Triangle(int[] data) {
	        //  Zack, G. W., Rogers, W. E. and Latt, S. A., 1977,
	        //  Automatic Measurement of Sister Chromatid Exchange Frequency,
	        // Journal of Histochemistry and Cytochemistry 25 (7), pp. 741-753
	        //
	        //  modified from Johannes Schindelin plugin
	        // 
	        // find min and max
	        int min = 0, dmax=0, max = 0, min2=0;
	        for (int i = 0; i < data.length; i++) {
	            if (data[i]>0){
	                min=i;
	                break;
	            }
	        }
	        if (min>0) min--; // line to the (p==0) point, not to data[min]

	        // The Triangle algorithm cannot tell whether the data is skewed to one side or another.
	        // This causes a problem as there are 2 possible thresholds between the max and the 2 extremes
	        // of the histogram.
	        // Here I propose to find out to which side of the max point the data is furthest, and use that as
	        //  the other extreme.
	        for (int i = 255; i >0; i-- ) {
	            if (data[i]>0){
	                min2=i;
	                break;
	            }
	        }
	        if (min2<255) min2++; // line to the (p==0) point, not to data[min]

	        for (int i =0; i < 256; i++) {
	            if (data[i] >dmax) {
	                max=i;
	                dmax=data[i];
	            }
	        }
	        // find which is the furthest side
	        //IJ.log(""+min+" "+max+" "+min2);
	        boolean inverted = false;
	        if ((max-min)<(min2-max)){
	            // reverse the histogram
	            //IJ.log("Reversing histogram.");
	            inverted = true;
	            int left  = 0;          // index of leftmost element
	            int right = 255; // index of rightmost element
	            while (left < right) {
	                // exchange the left and right elements
	                int temp = data[left]; 
	                data[left]  = data[right]; 
	                data[right] = temp;
	                // move the bounds toward the center
	                left++;
	                right--;
	            }
	            min=255-min2;
	            max=255-max;
	        }

	        if (min == max){
	            //IJ.log("Triangle:  min == max.");
	            return min;
	        }

	        // describe line by nx * x + ny * y - d = 0
	        double nx, ny, d;
	        // nx is just the max frequency as the other point has freq=0
	        nx = data[max];   //-min; // data[min]; //  lowest value bmin = (p=0)% in the image
	        ny = min - max;
	        d = Math.sqrt(nx * nx + ny * ny);
	        nx /= d;
	        ny /= d;
	        d = nx * min + ny * data[min];

	        // find split point
	        int split = min;
	        double splitDistance = 0;
	        for (int i = min + 1; i <= max; i++) {
	            double newDistance = nx * i + ny * data[i] - d;
	            if (newDistance > splitDistance) {
	                split = i;
	                splitDistance = newDistance;
	            }
	        }
	        split--;

	        if (inverted) {
	            // The histogram might be used for something else, so let's reverse it back
	            int left  = 0; 
	            int right = 255;
	            while (left < right) {
	                int temp = data[left]; 
	                data[left]  = data[right]; 
	                data[right] = temp;
	                left++;
	                right--;
	            }
	            return (255-split);
	        }
	        else
	            return split;
	    }


	    public int Yen(int[] data) {
	        // Implements Yen  thresholding method
	        // 1) Yen J.C., Chang F.J., and Chang S. (1995) "A New Criterion 
	        //    for Automatic Multilevel Thresholding" IEEE Trans. on Image 
	        //    Processing, 4(3): 370-378
	        // 2) Sezgin M. and Sankur B. (2004) "Survey over Image Thresholding 
	        //    Techniques and Quantitative Performance Evaluation" Journal of 
	        //    Electronic Imaging, 13(1): 146-165
	        //    http://citeseer.ist.psu.edu/sezgin04survey.html
	        //
	        // M. Emre Celebi
	        // 06.15.2007
	        // Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines
	        int threshold;
	        int ih, it;
	        double crit;
	        double max_crit;
	        double [] norm_histo = new double[256]; /* normalized histogram */
	        double [] P1 = new double[256]; /* cumulative normalized histogram */
	        double [] P1_sq = new double[256]; 
	        double [] P2_sq = new double[256]; 

	        double total =0;
	        for (ih = 0; ih < 256; ih++ ) 
	            total+=data[ih];

	        for (ih = 0; ih < 256; ih++ )
	            norm_histo[ih] = data[ih]/total;

	        P1[0]=norm_histo[0];
	        for (ih = 1; ih < 256; ih++ )
	            P1[ih]= P1[ih-1] + norm_histo[ih];

	        P1_sq[0]=norm_histo[0]*norm_histo[0];
	        for (ih = 1; ih < 256; ih++ )
	            P1_sq[ih]= P1_sq[ih-1] + norm_histo[ih] * norm_histo[ih];

	        P2_sq[255] = 0.0;
	        for ( ih = 254; ih >= 0; ih-- )
	            P2_sq[ih] = P2_sq[ih + 1] + norm_histo[ih + 1] * norm_histo[ih + 1];

	        /* Find the threshold that maximizes the criterion */
	        threshold = -1;
	        max_crit = Double.MIN_VALUE;
	        for ( it = 0; it < 256; it++ ) {
	            crit = -1.0 * (( P1_sq[it] * P2_sq[it] )> 0.0? Math.log( P1_sq[it] * P2_sq[it]):0.0) +  2 * ( ( P1[it] * ( 1.0 - P1[it] ) )>0.0? Math.log(  P1[it] * ( 1.0 - P1[it] ) ): 0.0);
	            if ( crit > max_crit ) {
	                max_crit = crit;
	                threshold = it;
	            }
	        }
	        return threshold;
	    }
	
}