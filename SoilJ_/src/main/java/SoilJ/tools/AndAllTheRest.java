package SoilJ.tools;

/**
 *SoilJ.tools is a collection of classes for SoilJ, 
 *an ImageJ plugin for the semi-automatized processing of 3-D X-ray images of soil columns
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

import java.awt.Color;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

/** 
 * AndAllTheRest is a SoilJ class containing some silly minor functions
 * 
 * @author John Koestel
 *
 */

public class AndAllTheRest implements PlugIn {

	
	public void run(String arg) {
				//ok, this is not needed..
	}

	public ArrayList<Integer> findFirstPositionInArray(double[] myArray, double value) {
		
		ArrayList<Integer> myPosition = new ArrayList<Integer>();
		for (int i = 0 ; i < myArray.length ; i++) if (myArray[i] == value) myPosition.add(i);
		
		return myPosition;
				
	}
	
	public int[] make1DIntFrom2DInt(int[][] someInt, int x, int y) {
		
		int[] outInt = new int[x * y];
		int cc = 0;
		
		for (int i = 0 ; i < x ; i++) {
			for (int j = 0 ; j < y ; j++) {				
				outInt[cc] = someInt[i][j];
				cc++;
			}
		}
		
		return outInt;
		
	}
	
	public boolean isContainedIn(int j, int[] array) {
	
		boolean isItTrue = false;
		
		for (int i = 0 ; i < array.length ; i++) {
			if (j == array[i]) {
				isItTrue = true;
			}
		}
		
		return isItTrue;
	}

	public static int[] getIndicesInOrder(double[] array) { //stolen from a JAVA student in the WWW.. god bless him!!
	    Map<Integer, Double> map = new HashMap<Integer, Double>(array.length);
	    for (int i = 0; i < array.length; i++)
	        map.put(i, array[i]);
	
	    List<Entry<Integer, Double>> l = new ArrayList<Entry<Integer, Double>>(map.entrySet());
	
	    Collections.sort(l, new Comparator<Entry<?, Double>>() {
	            @Override
	            public int compare(Entry<?, Double> e1, Entry<?, Double> e2) {
	                return e2.getValue().compareTo(e1.getValue());
	            }
	        });
	
	    int[] result = new int[array.length];
	    for (int i = 0; i < result.length; i++)
	        result[i] = l.get(i).getKey();
	
	    return result;
	}
	
	public static int[] getIndicesInOrder(int[] array) { //stolen from a JAVA student in the WWW.. god bless him!!
	    Map<Integer, Integer> map = new HashMap<Integer, Integer>(array.length);
	    for (int i = 0; i < array.length; i++)
	        map.put(i, array[i]);
	
	    List<Entry<Integer, Integer>> l = new ArrayList<Entry<Integer, Integer>>(map.entrySet());
	
	    Collections.sort(l, new Comparator<Entry<?, Integer>>() {
	            @Override
	            public int compare(Entry<?, Integer> e1, Entry<?, Integer> e2) {
	                return e2.getValue().compareTo(e1.getValue());
	            }
	        });
	
	    int[] result = new int[array.length];
	    for (int i = 0; i < result.length; i++)
	        result[i] = l.get(i).getKey();
	
	    return result;
	}
	
	public int[] fillMissingInts(int[] inInt, boolean[] isFilled) {
				
		int[] outInt = new int[inInt.length];
		
		ArrayList<Integer> isThere = new ArrayList<Integer>();
		ArrayList<Integer> isNotThere = new ArrayList<Integer>();
				
		for (int i = 0 ; i < inInt.length ; i++) {
			if (isFilled[i]) {
				isThere.add(i);
				outInt[i]=inInt[i];
			}
			else isNotThere.add(i);
		}
		
		if (isNotThere.isEmpty()) return inInt;
				
		for (int i = 0 ; i < isNotThere.size() ; i++) {
			
			int nows = isNotThere.get(i);
			int lower = -1;
			int upper = -1;
			if (nows < isThere.get(0)) {
				lower = isThere.get(isThere.size() - 1);
				upper = isThere.get(0);
			}
			if (nows > isThere.get(isThere.size() - 1)) {
				lower = isThere.get(isThere.size() - 1);
				upper = isThere.get(0);
			}
			if (nows > isThere.get(0) & nows < isThere.get(isThere.size() - 1)) {
				for (int k = nows - 1 ; k >= 0 ; k--) if(isThere.contains(k)) lower = k;
				int k = nows + 1; 
				while (upper < 0) {
					if (isThere.contains(k)) upper = k;
					else k++;
				}
			}
			
			outInt[isNotThere.get(i)] = (int)Math.round((inInt[lower] + inInt[upper]) / 2);
			
		}
		
		return outInt;
	}
	
	public int[] kickOutStrangeValues(int[] inInt) {
		
		Median jMed = new Median();
		RollerCaster rC = new RollerCaster();
		double myMedian = jMed.evaluate(rC.castInt2Double(inInt));
		
		for (int i = 0 ; i < inInt.length ; i++) 
			if (Math.abs((double)inInt[i] - myMedian) / myMedian > 0.25) 
				inInt[i] = 0; 
			
		return inInt;
		
	}
	
	
}