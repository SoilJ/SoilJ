package SoilJ.tools;

import java.util.ArrayList;

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

import java.util.Date;

import SoilJ.tools.MorphologyAnalyzer.PlusEdge;
import ij.plugin.PlugIn;
import sc.fiji.analyzeSkeleton.Edge;
import sc.fiji.analyzeSkeleton.Graph;
import sc.fiji.analyzeSkeleton.Point;
import sc.fiji.analyzeSkeleton.Vertex;

/** 
 * RollerCaster is a SoilJ class with subroutines for casting arrays.
 * 
 * @author John Koestel
 *
 */

public class RollerCaster implements PlugIn {
	
	final double MilliSecondsOfOneDay = 24*60*60*1000;

	
	public MorphologyAnalyzer.PlusGraph cast2PlusGraph(Graph nowGraph, int numberOfLastSlice, int numberOfSlicesAdded) {
		
		MorphologyAnalyzer mA = new MorphologyAnalyzer();		
		
		MorphologyAnalyzer.PlusGraph outGraph = mA.new PlusGraph();
		
		ArrayList<Edge> nowEdges = nowGraph.getEdges();
		
		for (int i = 0 ; i < nowEdges.size() ; i++)	{
			
			//check if vertices are in added slices
			Edge nowEdge = nowEdges.get(i);			
						
			//only add edge if it is not completely outside the added slices otherwise modify or kick out
			PlusEdge outEdge = checkIfEdgeIsToBeKept(nowEdge, numberOfLastSlice, numberOfSlicesAdded);
					
			if (outEdge != null) {				
				outGraph.plusEdge.add(outEdge);
			}
			
		}
		
		ArrayList<Vertex> nowVerts = nowGraph.getVertices();
		for (int i = 0 ; i < nowVerts.size() ; i++) {
			
			Vertex nowVert = nowVerts.get(i);
			ArrayList<Point> nowPoints = nowVert.getPoints();
			ArrayList<Edge> nowEdgesV = nowVert.getBranches();
			Edge nowPrede = nowVert.getPredecessor();
			MorphologyAnalyzer.PlusVertex newVert = mA.new PlusVertex();
			
			//check if vertex is not within added slices
			double zSum = 0;		
			for (int j = 0 ; j < nowPoints.size() ; j++) {
				zSum += nowPoints.get(j).z;		
			}
			
			//only add edge if it is not completely outside the added slices otherwise modify or kick out
			boolean keepVertex = true;					
	
			if (zSum / (double)nowPoints.size() > numberOfSlicesAdded & zSum / (double)nowPoints.size() < numberOfLastSlice / 2) keepVertex = true;			
			if (zSum / (double)nowPoints.size() <= numberOfLastSlice + numberOfSlicesAdded & zSum / (double)nowPoints.size() > numberOfLastSlice / 2) keepVertex = true;			
			if (zSum / (double)nowPoints.size() <= numberOfSlicesAdded & zSum / (double)nowPoints.size() < numberOfLastSlice / 2) keepVertex = false;
			if (zSum / (double)nowPoints.size() >= numberOfLastSlice + numberOfSlicesAdded & zSum / (double)nowPoints.size() > numberOfLastSlice / 2) keepVertex = false;

			if (keepVertex) {
				
				for (int j = 0 ; j < nowPoints.size() ; j++) newVert.addPoint(nowPoints.get(j));
			
				newVert.setPredecessor(nowPrede);
				
				for (int j = 0 ; j < nowEdgesV.size() ; j++) {
					
					Edge nowEdge = nowEdgesV.get(j);				
	
					//check if vertices are in added slices	
					PlusEdge outEdge = checkIfEdgeIsToBeKept(nowEdge, numberOfLastSlice, numberOfSlicesAdded);		
								
					if (outEdge != null) {					
						newVert.plusBranches.add(outEdge);
					}
				}
			}
			outGraph.addVertex(newVert);
		}
		
		return outGraph;
	}
	
	public MorphologyAnalyzer.PlusGraph cast2PlusGraph(Graph nowGraph) {
		
		MorphologyAnalyzer mA = new MorphologyAnalyzer();		
		
		MorphologyAnalyzer.PlusGraph outGraph = mA.new PlusGraph();
		
		ArrayList<Edge> nowEdges = nowGraph.getEdges();
		
		for (int i = 0 ; i < nowEdges.size() ; i++)	{
			
			//check if vertices are in added slices
			Edge nowEdge = nowEdges.get(i);		
			
			Vertex v1 = nowEdge.getV1();
			Vertex v2 = nowEdge.getV2();
			ArrayList<Point> slabs = nowEdge.getSlabs();
			double edgeLength = nowEdge.getLength();
						
			//only add edge if it is not completely outside the added slices otherwise modify or kick out
			PlusEdge outEdge = mA.new PlusEdge(v1, v2, slabs, edgeLength);
			
			outGraph.plusEdge.add(outEdge);						
		}
		
		ArrayList<Vertex> nowVerts = nowGraph.getVertices();
		for (int i = 0 ; i < nowVerts.size() ; i++) {			
			Vertex nowVert = nowVerts.get(i);
			outGraph.addVertex(nowVert);
		}
		
		return outGraph;
	}
	
	public MorphologyAnalyzer.PlusVertex cast2PlusVertex(Vertex nowVert) {
		
		MorphologyAnalyzer mA = new MorphologyAnalyzer();		
			
		ArrayList<Point> nowPoints = nowVert.getPoints();
		ArrayList<Edge> nowEdgesV = nowVert.getBranches();
		Edge nowPrede = nowVert.getPredecessor();
		MorphologyAnalyzer.PlusVertex newVert = mA.new PlusVertex();
			
			
		for (int j = 0 ; j < nowPoints.size() ; j++) newVert.addPoint(nowPoints.get(j));
			
		newVert.setPredecessor(nowPrede);
				
		for (int j = 0 ; j < nowEdgesV.size() ; j++) {
					
			Edge nowEdge = nowEdgesV.get(j);		
			Vertex v1 = nowEdge.getV1();
			Vertex v2 = nowEdge.getV2();
			ArrayList<Point> slabs = nowEdge.getSlabs();
			double edgeLength = nowEdge.getLength();			
	
			PlusEdge outEdge = mA.new PlusEdge(v1, v2, slabs, edgeLength);			
								
			newVert.plusBranches.add(outEdge);		
		}
		
		return newVert;
	}
	
	public MorphologyAnalyzer.PlusEdge checkIfEdgeIsToBeKept(Edge nowEdge, int numberOfLastSlice, int numberOfSlicesAdded) {
		
		MorphologyAnalyzer mA = new MorphologyAnalyzer();
		
		Vertex v1 = nowEdge.getV1();
		Vertex v2 = nowEdge.getV2();
		ArrayList<Point> slabs = nowEdge.getSlabs();
		double edgeLength = nowEdge.getLength();
		ArrayList<Point> v1Points = nowEdge.getV1().getPoints();
		ArrayList<Point> v2Points = nowEdge.getV2().getPoints();
		
		boolean v1out = false;	
		boolean v2out = false;
		if ((v1Points.size() == 1) & (v1Points.get(0).z < numberOfSlicesAdded)) v1out = true;
		if ((v1Points.size() == 1) & (v1Points.get(0).z >= numberOfLastSlice + numberOfSlicesAdded)) v1out = true;
		if ((v2Points.size() == 1) & (v2Points.get(0).z < numberOfSlicesAdded)) v2out = true;
		if ((v2Points.size() == 1) & (v2Points.get(0).z >= numberOfLastSlice + numberOfSlicesAdded)) v2out = true;
		
		boolean keepEdge = true;
	
		if (v1out & v2out) keepEdge = false;
		
		if (!v1out & v2out) {
			
			double zSum = 0;
			int zMax = 0;
			int zMin = Integer.MAX_VALUE;
			for (int j = 0 ; j < v1Points.size() ; j++) {
				zSum += v1Points.get(j).z;
				if (v1Points.get(j).z > zMax) zMax = v1Points.get(j).z;
				if (v1Points.get(j).z < zMin) zMin = v1Points.get(j).z;
			}
							
			//in case the vertex it at the top of the column
			if (zSum / (double)v1Points.size() > numberOfSlicesAdded & zSum / (double)v1Points.size() < numberOfLastSlice / 2) {					
				int newV2Z = zMax - 1;
				ArrayList<Integer> toBeRemoved = new ArrayList<Integer>();
				if (zMax > numberOfSlicesAdded) {
					newV2Z = numberOfSlicesAdded;
					for (int k = 0 ; k < slabs.size() ; k++) if (slabs.get(k).z < newV2Z) toBeRemoved.add(k);
					slabs.removeAll(toBeRemoved);
				}
				else {
					while (slabs.size() > 0) slabs.remove(0);		
				}
				
				Point newPoint = new Point(v1Points.get(0).y, v1Points.get(0).y, newV2Z);					
							
				v2Points.remove(0);
				v2Points.add(newPoint);			
				
				edgeLength = 1;
			}
			
			//in case the vertex it at the bottom of the column
			if (zSum / (double)v1Points.size() < numberOfLastSlice + numberOfSlicesAdded & zSum / (double)v1Points.size() > numberOfLastSlice / 2) {					
				
				int newV2Z = zMax + 1;
				ArrayList<Integer> toBeRemoved = new ArrayList<Integer>();
				if (zMax < numberOfLastSlice + numberOfSlicesAdded) {
					newV2Z = numberOfLastSlice + numberOfSlicesAdded;
					for (int k = 0 ; k < slabs.size() ; k++) if (slabs.get(k).z > newV2Z) toBeRemoved.add(k);
					slabs.removeAll(toBeRemoved);
				}
				else {
					while (slabs.size() > 0) slabs.remove(0);
				}
				
				Point newPoint = new Point(v1Points.get(0).y, v1Points.get(0).y, newV2Z);					
				v2Points.remove(0);					
				v2Points.add(newPoint);	
				
				edgeLength = 1;
			}
			
			//also check if the second vertex contains of more than one point but is within the added slices
			if (zSum / (double)v1Points.size() <= numberOfSlicesAdded & zSum / (double)v1Points.size() < numberOfLastSlice / 2) {
				keepEdge = false;
			}
			
			if (zSum / (double)v1Points.size() >= numberOfLastSlice + numberOfSlicesAdded & zSum / (double)v1Points.size() > numberOfLastSlice / 2) {
				keepEdge = false;
			}				
		}
		
		if (v1out & !v2out) {
			
			double zSum = 0;
			int zMax = 0;
			int zMin = Integer.MAX_VALUE;
			for (int j = 0 ; j < v2Points.size() ; j++) {
				zSum += v2Points.get(j).z;
				if (v2Points.get(j).z > zMax) zMax = v2Points.get(j).z;
				if (v2Points.get(j).z < zMin) zMin = v2Points.get(j).z;
			}
							
			//in case the vertex it at the top of the column
			//in case the vertex it at the top of the column
			if (zSum / (double)v2Points.size() > numberOfSlicesAdded & zSum / (double)v2Points.size() < numberOfLastSlice / 2) {					
				int newV1Z = zMax - 1;
				ArrayList<Integer> toBeRemoved = new ArrayList<Integer>();
				if (zMax > numberOfSlicesAdded) {
					newV1Z = numberOfSlicesAdded;
					for (int k = 0 ; k < slabs.size() ; k++) if (slabs.get(k).z < newV1Z) toBeRemoved.add(k);
					slabs.removeAll(toBeRemoved);
				}
				else {
					while (slabs.size() > 0) slabs.remove(0);		
				}
				
				Point newPoint = new Point(v2Points.get(0).y, v2Points.get(0).y, newV1Z);					
							
				v1Points.remove(0);
				v1Points.add(newPoint);			
				
				edgeLength = 1;
			}
			
			//in case the vertex it at the bottom of the column
			if (zSum / (double)v2Points.size() < numberOfLastSlice + numberOfSlicesAdded & zSum / (double)v2Points.size() > numberOfLastSlice / 2) {					
				
				int newV2Z = zMax + 1;
				ArrayList<Integer> toBeRemoved = new ArrayList<Integer>();
				if (zMax < numberOfLastSlice + numberOfSlicesAdded) {
					newV2Z = numberOfLastSlice + numberOfSlicesAdded;
					for (int k = 0 ; k < slabs.size() ; k++) if (slabs.get(k).z > newV2Z) toBeRemoved.add(k);
					slabs.removeAll(toBeRemoved);
				}
				else {
					while (slabs.size() > 0) slabs.remove(0);
				}
				
				Point newPoint = new Point(v2Points.get(0).y, v2Points.get(0).y, newV2Z);					
				v1Points.remove(0);					
				v1Points.add(newPoint);	
				
				edgeLength = 1;
			}
			
			//also check if the second vertex contains of more than one point but is within the added slices
			if (zSum / (double)v2Points.size() <= numberOfSlicesAdded & zSum / (double)v2Points.size() < numberOfLastSlice / 2) {
				keepEdge = false;
			}
			
			if (zSum / (double)v2Points.size() >= numberOfLastSlice + numberOfSlicesAdded & zSum / (double)v2Points.size() > numberOfLastSlice / 2) {
				keepEdge = false;
			}			
		}
		
		MorphologyAnalyzer.PlusEdge outEdge = null;
		if (keepEdge) {
			outEdge = mA.new PlusEdge(v1, v2, slabs, edgeLength);
			outEdge.plusV1 = cast2PlusVertex(v1);
			outEdge.plusV2 = cast2PlusVertex(v2);
		}
				
		return outEdge;
	}
	
	public void run(String arg) {
				//ok, this is not needed..
	}

	public float[] intVector2Float(int[] myVector){
		
		float[] outVector =  new float[myVector.length];
		
		for (int i = 0 ; i < myVector.length ; i++) {
			outVector[i] = (float)myVector[i];
		}
		
		return outVector;
	}
	
	public double[] intVector2Double(int[] myVector){
		
		double[] outVector =  new double[myVector.length];
		
		for (int i = 0 ; i < myVector.length ; i++) {
			outVector[i] = (double)myVector[i];
		}
		
		return outVector;
	}

	public double[] castInt2Double(int[] myArray) {
		
		double[] out = new double[myArray.length];
		
		for (int i = 0 ; i < out.length ; i++) out[i] = (double)myArray[i];
		
		return out;
	}

	public float[] castInt2Float(int[] myArray) {
		
		float[] out = new float[myArray.length];
		
		for (int i = 0 ; i < out.length ; i++) out[i] = (float)myArray[i];
		
		return out;
	}

	public float[] castDouble2Float(double[] myArray) {
		
		float[] out = new float[myArray.length];
		
		for (int i = 0 ; i < out.length ; i++) out[i] = (float)myArray[i];
		
		return out;
	}
	
	public int[] castDouble2Int(double[] myArray) {
		
		int[] out = new int[myArray.length];
		
		for (int i = 0 ; i < out.length ; i++) out[i] = (int)Math.round(myArray[i]);
		
		return out;
	}
	
	public ArrayList<Integer> castIntArray2List(int[] myArray) {
		
		ArrayList<Integer> myList = new ArrayList<Integer>();
		
		for (int i = 0 ; i < myArray.length ; i++) myList.add((int)Math.round(myArray[i]));
		
		return myList;
	}
	
	public ArrayList<Integer> castDoubleArray2IntegerList(double[] myArray) {
		
		ArrayList<Integer> myList = new ArrayList<Integer>();
		
		for (int i = 0 ; i < myArray.length ; i++) myList.add((int)Math.round(myArray[i]));
		
		return myList;
	}
	
	public int[] castFloat2Int(float[] myArray) {
		
		int[] out = new int[myArray.length];
		
		for (int i = 0 ; i < out.length ; i++) out[i] = (int)Math.round(myArray[i]);
		
		return out;
	}
	
	public double[] castFloat2Double(float[] myArray) {
		
		double[] out = new double[myArray.length];
		
		for (int i = 0 ; i < out.length ; i++) out[i] = (double)myArray[i];
		
		return out;
	}
	
	public int[][] castFloatFloat2IntInt(float[][] myArray) {
		
		int[][] out = new int[myArray.length][myArray[0].length];
		
		for (int i = 0 ; i < out.length ; i++) {
			for (int j = 0 ; j < out[0].length ; j++) {
				out[i][j] = (int)Math.round(myArray[i][j]);
			}
		}
		
		return out;
	}
	
	public double date2Days(Date date){
	    //  convert a date to an integer and back again
	    double currentTime=date.getTime();
	    currentTime = currentTime/MilliSecondsOfOneDay;
	    return currentTime; 
	}

	public double[] reverseDoubleArray(double[] myArray) {
		
		double[] oldArray = myArray;
		for (int i = 0 ; i < myArray.length ; i++) myArray[i] = oldArray[myArray.length - i - 1];
		
		return myArray;
		
	}
	
}