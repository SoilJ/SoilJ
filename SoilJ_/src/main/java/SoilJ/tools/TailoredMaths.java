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

import org.apache.commons.math3.stat.StatUtils;

import ij.plugin.PlugIn;

/** 
 * TailoredMaths is a SoilJ class collecting mathematical subroutines.
 * 
 * @author John Koestel
 *
 */

public class TailoredMaths implements PlugIn {
	
	public void run(String arg) {
				//ok, this is not needed..
	}

	public double[] oneDMedianFilter(double[] grayValue, int footprint) {
		
		int j, k;
		
		//catch even footprints
		if (Math.floorMod(footprint,2) == 0) footprint--;
		
		//calculate first and last entries with values
		int firstEntry = Math.floorDiv(footprint, 2) + 1;
		int lastEntry = grayValue.length - 1 - firstEntry;
		
		//apply a median filter to greyAtThisAngle
		double[] mGreyAtThisAngle = new double[grayValue.length];
		double[] withinTheFilter = new double[footprint];		
		for (j = (footprint - 1) / 2 ; j < grayValue.length - (footprint - 1) / 2; j++) {
			for (k = 0 ; k < footprint ; k++) withinTheFilter[k] = grayValue[j - (footprint - 1) / 2 + k]; 
			mGreyAtThisAngle[j] = StatUtils.percentile(withinTheFilter, 50);
		}
		for (j = 0 ; j < firstEntry ; j++) mGreyAtThisAngle[j] = mGreyAtThisAngle[firstEntry];
		for (j = lastEntry + 1 ; j < grayValue.length ; j++) mGreyAtThisAngle[j] = mGreyAtThisAngle[lastEntry];
		
		return mGreyAtThisAngle;
	}
	
	public int[][] applyCorrFactor(int[][] interCount, int xl, int yl, double correctionFactor) {
		
		int[][] cInterCount = interCount.clone();
		
		for (int i = 0 ; i < xl ; i++) {
			for (int j = 0 ; j < yl ; j++) {
				cInterCount[i][j] =  (int)Math.round((double)interCount[i][j] / correctionFactor);
			}
		}
				
		return cInterCount;
	}
	
	public double[] compileStatistics(int[][] interCount, int xl, int yl, double correctionFactor) {
		
		RollerCaster rC = new RollerCaster();
		AndAllTheRest aatr = new AndAllTheRest();
		
		int[][] cInterCount = applyCorrFactor(interCount, xl, yl, correctionFactor);
		
		double[] aRe = new double[4]; 
				
		aRe[0] = StatUtils.sum(rC.castInt2Double(aatr.make1DIntFrom2DInt(cInterCount, xl, yl)));
		aRe[1]= StatUtils.mean(rC.castInt2Double(aatr.make1DIntFrom2DInt(cInterCount, xl, yl)));
		aRe[2] = Math.sqrt(StatUtils.variance(rC.castInt2Double(aatr.make1DIntFrom2DInt(cInterCount, xl, yl))));
		
		aRe[3] = StatUtils.sum(rC.castInt2Double(aatr.make1DIntFrom2DInt(cInterCount, xl, yl)));
 		
		return aRe;
		
	}
	
	public static double calculateWeightedAverage(double[] target, int[] weight) {
		
		double weightsum = 0;
		double targetsum = 0;	
		
		for (int i = 0 ; i < target.length ; i++) {
			weightsum += (double)weight[i];
			targetsum += target[i] * (double)weight[i];
		}
		
		if (weightsum > 0) return targetsum / weightsum;
		else return 0;
	}
	
	public static double meanWithoutZeros(double[] target) {
		
		double cc = 0;
		double targetsum = 0;
		
		for (int i = 0 ; i < target.length ; i++) {
			
			targetsum += target[i];
			
			if (target[i] > 0) {
				cc++;
			}
			
		}
		if (cc == 0 | targetsum == 0) return 1;
		else return targetsum / cc;

	}
	
	public double[] oneDMeanFilter(double[] grayValue, int footprint) {
		
		int j, k;
		
		//catch even footprints
		if (Math.floorMod(footprint,2) == 0) footprint--;
		
		//calculate first and last entries with values
		int firstEntry = Math.floorDiv(footprint, 2) + 1;
		int lastEntry = grayValue.length - 1 - firstEntry;
		
		//apply a median filter to greyAtThisAngle
		double[] newGrayValue = new double[grayValue.length];
		double[] withinTheFilter = new double[footprint];
		for (j = 0 ; j < grayValue.length ; j++) newGrayValue[j] = 0; //nullning av minimum filtered values
		for (j = (footprint - 1) / 2 ; j < grayValue.length - (footprint - 1) / 2; j++) {
			for (k = 0 ; k < footprint ; k++) withinTheFilter[k] = grayValue[j - (footprint - 1) / 2 + k]; 
			newGrayValue[j] = StatUtils.mean(withinTheFilter);
		}
		for (j = 0 ; j < firstEntry ; j++) newGrayValue[j] = newGrayValue[firstEntry];
		for (j = lastEntry + 1 ; j < grayValue.length ; j++) newGrayValue[j] = newGrayValue[lastEntry];
		
		return newGrayValue;
	}	
	
	public double[] oneDMinimumFilter(double[] grayValue, int footprint) {
		
		int j, k;
		
		//apply a median filter to greyAtThisAngle
		double[] newGrayValue = new double[grayValue.length];
		double[] withinTheFilter = new double[footprint];
		for (j = 0 ; j < grayValue.length ; j++) newGrayValue[j] = 0; //nullning av minimum filtered values
		for (j = (footprint - 1) / 2 ; j < grayValue.length - (footprint - 1) / 2; j++) {
			for (k = 0 ; k < footprint ; k++) withinTheFilter[k] = grayValue[j - (footprint - 1) / 2 + k]; 
			newGrayValue[j] = StatUtils.min(withinTheFilter);
		}
		
		return newGrayValue;
	}

	public double[] oneDMaximumFilter(double[] grayValue, int footprint) {
		
		int j, k;
		
		//apply a median filter to greyAtThisAngle
		double[] newGrayValue = new double[grayValue.length];
		double[] withinTheFilter = new double[footprint];
		for (j = 0 ; j < grayValue.length ; j++) newGrayValue[j] = 0; //nullning av minimum filtered values
		for (j = (footprint - 1) / 2 ; j < grayValue.length - (footprint - 1) / 2; j++) {
			for (k = 0 ; k < footprint ; k++) withinTheFilter[k] = grayValue[j - (footprint - 1) / 2 + k]; 
			newGrayValue[j] = StatUtils.max(withinTheFilter);
		}
		
		return newGrayValue;
	}
	
	public double[] filterAnArray(double[] myArray, int windowSize, String myFilter) {
		
		double[] myOutRay = new double[myArray.length];
		
		if (myFilter.compareToIgnoreCase("Minimum") == 0) myOutRay = oneDMinimumFilter(myArray, windowSize);
		if (myFilter.compareToIgnoreCase("Mean") == 0) myOutRay = oneDMeanFilter(myArray, windowSize);
		if (myFilter.compareToIgnoreCase("Median") == 0) myOutRay = oneDMedianFilter(myArray, windowSize);
		if (myFilter.compareToIgnoreCase("Maximum") == 0) myOutRay = oneDMaximumFilter(myArray, windowSize);
	
		return myOutRay;
	}
	
	public double[] fillStartAndEndZeros(double[] myArray) {
		
		int first = 0;
		int last = 0;
		double[] outArray = myArray;
		
		for (int i = 0 ; i < myArray.length ; i++) {
			if (first == 0) if (myArray[i] != 0) first = i;
			if (first > 0) if (myArray[i] > 0) last = i;
		}
		
		for (int i = 0 ; i < first ; i++) {
			outArray[i] = myArray[first];
		}
		for (int i = last + 1 ; i < myArray.length ; i++) {
			outArray[i] = myArray[last];
		}
		
		return outArray;
		
	}
	
	public double[] oneDFirstDerivative(double[] mGreyAtThisAngle, int footprintOfMedianFilter) {
		
		int j;		
		double[] dMG = new double[mGreyAtThisAngle.length];
		for (j = 0 ; j < mGreyAtThisAngle.length ; j++) dMG[j] = 0; //nullning av first derivative
		for (j = (footprintOfMedianFilter - 1) / 2 ; j < mGreyAtThisAngle.length - (footprintOfMedianFilter - 1) / 2 - 1; j++) {
			dMG[j] = mGreyAtThisAngle[j + 1] - mGreyAtThisAngle[j];
		}
		
		return dMG;
	}
	
	public double stdWithoutZeros(double[] vector) {
		
		int cc=0;
		for (int i = 0 ; i < vector.length ; i++) if (vector[i] == 0) cc++;
		
		double[] newVec = new double[vector.length - cc];
		int ii = 0;
		for (int i = 0 ; i < vector.length ; i++) {
			if (vector[i] > 0) {
				newVec[ii] = vector[i];
				ii++;
			}
		}
		
		return Math.sqrt(StatUtils.variance(newVec));
		
	}

	public double[] LinearLOESSFilter(double[] data, int windowHalfSize) {
				
		FitStuff fs =  new FitStuff();
		
		double[] smoothie = new double[data.length];
		
		//create prolonged data series for smoothing
		double[] xData = new double[data.length + 2 * windowHalfSize];
		
		for (int i = 0 ; i < windowHalfSize ; i++) {
			xData[i] = data[0];
			xData[windowHalfSize + data.length + i] = data[data.length - 1];
		}
		
		for (int i = 0 ; i < data.length ; i++) xData[windowHalfSize + i] = data[i];
		
		//go and smooth!
		for (int i = 0 ; i < data.length ; i++) {
					
			double[] nowWindow = new double[2 * windowHalfSize];
			double[] xWindow = new double[2 * windowHalfSize];
			
			for (int j = 0 ; j < 2 * windowHalfSize ; j++) nowWindow[j] = xData[i + j];
			
			for (int j = 0 ; j < nowWindow.length ; j++) xWindow[j] = j;
			
			smoothie[i] = fs.evalLinearFunction(fs.fitLinearFunction(xWindow, nowWindow), windowHalfSize);				
					
		}
		
		return smoothie;
	}

	public double[][] getXYOfEllipseFromAngle(double[] angleAtThisAngle, double xCenter, double yCenter, double major, double minor, double theta) {
		
		int j;
	
		double[][] xy = new double[angleAtThisAngle.length][2];
	
		for (j = 0 ; j < angleAtThisAngle.length ; j++) {			
			double alpha = angleAtThisAngle[j] - theta + Math.PI / 2;
			double a = major;
			double b = minor;			
			xy[j][0] = xCenter - a * Math.cos(alpha) * Math.cos(theta) + b * Math.sin(alpha) * Math.sin(theta);
			xy[j][1] = yCenter + a * Math.cos(alpha) * Math.sin(theta) + b * Math.sin(alpha) * Math.cos(theta);			
		}
	
		return xy;
	}

	public float[] logitCorrectionFunctionCalculator(int i, double[][][] bhc, int bhcEntry, float[] sR) {
		
		float[] ftot = new float[sR.length];
		int standardRadius = sR.length;
		
		double my = bhc[i][0][bhcEntry];
		double sd = bhc[i][1][bhcEntry];
		double thresh = bhc[i][2][bhcEntry];
		double skipslast = bhc[i][3][bhcEntry];
		double r = bhc[i][4][bhcEntry];
		double dy = bhc[i][5][bhcEntry];		
		double[] f = new double[(int)Math.round(standardRadius - skipslast - thresh)];
		
		for (int j = 0 ; j < (int)Math.round(standardRadius - skipslast - thresh); j++) {
			f[(int)Math.round(standardRadius - skipslast - thresh) - j - 1] = dy * (1 - 1 / (1 + Math.exp(-(j - my) / sd))) + r;	 //switch vector around!!!! important!!
		}
		
		int cc=0;
		for (int j = 0 ; j < skipslast ; j++) {
			ftot[cc] = (float)StatUtils.min(f);
			cc++;
		}
		for (int j = 0 ; j < f.length; j++) {
			ftot[cc] = (float)f[j];
			cc++;
		}
		for (int j = 0 ; j < thresh ; j++) {
			ftot[cc] = (float)StatUtils.max(f);
			cc++;
		}
			
		return ftot;
	}

	public double min(double[] myArray) {		
		double myMin = StatUtils.min(myArray);		
		return myMin;
	}
	
	public double minNoZero(double[] myArray) {
		double myMax = StatUtils.max(myArray);
		for (int i = 0 ; i < myArray.length ; i++) if (myArray[i] == 0) myArray[i] = myMax;
		double myMin = StatUtils.min(myArray);		
		return myMin;
	}

	public double max(double[] myArray) {		
		double myMax = StatUtils.max(myArray);		
		return myMax;
	}

	public static int factorial(int n) {
	    int fact = 1; // this  will be the result
	    for (int i = 1; i <= n; i++) {
	        fact *= i;
	    }
	    return fact;
	}
	
	public static int intSum(int[] myInt) {
		
		int mySum = 0;
		
		for (int i = 0 ; i < myInt.length ; i++) mySum += myInt[i];
		
		return mySum;
		
	}
	
}
	