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


import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.stat.StatUtils;

import SoilJ.tools.RoiHandler.ColumnRoi;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import inra.ijpb.binary.ChamferWeights3D;
import inra.ijpb.binary.conncomp.FloodFillComponentsLabeling3D;
import inra.ijpb.binary.distmap.DistanceTransform3D;
import inra.ijpb.binary.distmap.DistanceTransform3DFloat;
import inra.ijpb.binary.distmap.DistanceTransform3DShort;
import inra.ijpb.geometry.Box3D;
import inra.ijpb.label.LabelImages;
import inra.ijpb.measure.region3d.BoundingBox3D;
import inra.ijpb.measure.region3d.IntrinsicVolumesAnalyzer3D;
import inra.ijpb.measure.region3d.IntrinsicVolumesAnalyzer3D.Result;
import process3d.Dilate_;
import sc.fiji.analyzeSkeleton.AnalyzeSkeleton_;
import sc.fiji.analyzeSkeleton.Edge;
import sc.fiji.analyzeSkeleton.Graph;
import sc.fiji.analyzeSkeleton.Point;
import sc.fiji.analyzeSkeleton.SkeletonResult;
import sc.fiji.analyzeSkeleton.Vertex;
import sc.fiji.localThickness.LocalThicknessWrapper;
import sc.fiji.skeletonize3D.Skeletonize3D_;

/** 
 * is a SoilJ class
 * 
 * @author John Koestel
 *
 */

public class MorphologyAnalyzer implements PlugIn {

	public void run(String arg) {
		//ok, this is not needed..
	}

	public class SurfaceStatistics {
		
		public int highestElevation;
		public int medianElevation;
		public int meanElevation;
		public int lowestElevation;
		
		public int highestIntrusion;
		public int medianIntrusion;
		public int meanIntrusion;
		public int lowestIntrusion;		
	}		
	
	public class ProfileStatistics {
		
		public int[] numberOfNonZeros;
		public double[] mean;
		public double[] geomean;
		public double[] median;
		public double[] mode;
		public double[] std;
		public double[] mini;
		public double[] maxi;		
		
	}
	
	public class FractalProperties {
		
		public double volumeFractalDim;
		public double surfaceFractalDim;
		
	}
	
	public class MLJLabellingResults {
		
		public double[] volume;
		public double[] surface;
		public double[] euler;
		
	}
	
	public class FloydWarshallReturn {
		
		public double[][] adjacencyMatrix;
		public double[][] distanceMatrix;
		public int[][] predecessorMatrix;	
		
		ArrayList<Vertex> topVertices;
		ArrayList<Vertex> bottomVertices;
		ArrayList<Integer> touchesTop;
		ArrayList<Integer> touchesBottom;
		
		double[][] top2BottomConnection;
		double[][] top2BottomEuklideanConnection;
		
		ArrayList<PlusEdge> edgeList;
		ArrayList<Vertex> vertexList;
		
		double shortestPath;
		
		ReconstructedPath shortest;
		
	}
	
	public class ReconstructedPath {
		
		public ArrayList<Point> pathPoints;
		public ArrayList<Double> pathDistances;
		
		public double edgeLength;
		public double edgeEuclideanDistance;
		public double edgeTortuosity;
		public double globalTortuosity;
		public double bottleNeck;
		
		public int numberOfVertices;
		public double sumOfCN;
		
	}
	
	public class ConnectingGraphs {
		
		ArrayList<Integer> percolating = new ArrayList<Integer>();
		ArrayList<Integer> con2Top = new ArrayList<Integer>();
		ArrayList<Integer> con2Bot = new ArrayList<Integer>();
		
	}
	
	public class NetworkProps {
		
		FloydWarshallReturn[] mFWR;
		
		double fastestPathResistance;
		double fastestPathBottleneck;
		double fastestPathEdgeTortuosity;
		double fastestPathGlobalTortuosity;
		int fastestPathVertices;
		int fastestPathSlabs;
		int fastestPathVerticesPerSlab;
		double meanFastestPathCN;   //CN == coordination number
		
		int globalEdges; 
		int globalSlabs;
		int globalVertices;
		int globalDeadEnds;
		int numberOfSingleVertexGraphs;
		int numberOfPercolatingTopPoints;
		int numberOfPercolatingBotPoints;
		int numberOfCon2TopPoints;
		int numberOfCon2BotPoints;
				
		double globalEdgeLength;
		double globalEuclidEdgeLength;
		double globalTortuosity;
		double globalMeanTortuosity;
		double globalMeanCN; 		//CN == coordination number
		double globalMeanAngleFromZ;
		
		double percolatingTortuosity;
		int percolatingSlabs;
		int percolatingVertices;
		int percolatingEdges;
		int percolatingDeadEnds;
		double percolatingMeanCN; 		//CN == coordination number		
		double percolatingMeanAngleFromZ;
		
		double con2TopTortuosity;
		int con2TopSlabs;
		int con2TopVertices;
		int con2TopEdges;
		int con2TopDeadEnds;
		double con2TopMeanCN; 		//CN == coordination number		
		double con2TopMeanAngleFromZ;
		
		double con2BotTortuosity;
		int con2BotSlabs;
		int con2BotVertices;
		int con2BotEdges;
		int con2BotDeadEnds;
		double con2BotMeanCN; 		//CN == coordination number		
		double con2BotMeanAngleFromZ;
		
	}
	
	public class WRCPhaseProps {
		
		public double[] tensionAtTopInMM;
		public double[] tensionAtCenterInMM;
		public double[] tensionAtBottomInMM;
		
		public double[] pressureAtTopInMM;
		public double[] pressureAtCenterInMM;
		public double[] pressureAtBottomInMM;
		
		public double[] areaInMM2;
		public double[] heightInMM;
		public double[] bulkVolumeInMM3;

		public double[] theta_w;
		public double[] sigma_w;
		public double[] chi_w;
		public double[] gamma_w;
		public double[] fractalDim_w;
		public double[] percolates_w;
		public double[] thetaLC_w;
		public double[] thetaPerc_w;
		public double[] dc_w;
		
		public double[] theta_a;
		public double[] sigma_a;
		public double[] chi_a;
		public double[] gamma_a;
		public double[] fractalDim_a;
		public double[] percolates_a;
		public double[] thetaLC_a;
		public double[] thetaPerc_a;
		public double[] dc_a;
		public double[] depthOfPenetrationInMM_a;
		
		public double[] averageDistanceFromAeratedPoreInMM; 
		public double[] fractionLessThan3MMFromPhaseBoundary;
		public double[] fractionMoreThan3MMFromPhaseBoundary;
		
		public boolean enforceIntegerTensions;
		
	}
	
	public void extractPoresizeDistro(String savePath, ImagePlus nowTiff, double[] classBounds) {
		
			
		InputOutput jIO = new InputOutput();
		
		int[] psd = new int[classBounds.length];
		
		//get histogram
		for (int j = 0 ; j < nowTiff.getNSlices() ; j++) { 
		
			IJ.showStatus("Adding to histogram in layer " + (j + 1));
			
			nowTiff.setPosition(j+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
						 
			for (int x = 0 ; x < nowTiff.getWidth() ; x++) {
				for (int y = 0 ; y < nowTiff.getHeight() ; y++) {
					double nowPix = nowIP.getPixelValue(x, y);
					
					//collect according for class boundaries	
					if (nowPix > 0) {						
						if (nowPix >= 0) {							
							//catch columns with very large pores
							int enterHere = (int)Math.floor(nowPix*2);
							if (enterHere > classBounds.length - 1) enterHere = classBounds.length - 1;
							psd[enterHere]++;
						}
					}
				}
			}				
		}		
		
		//save histogram
		jIO.writePoreSizeDistribution(savePath, psd, classBounds);
	}	
	
	public void extract32BitHisto(String savePath, ImagePlus nowTiff, double[] classBounds) {
		
		
		InputOutput jIO = new InputOutput();
		
		int[] psd = new int[classBounds.length];
		
		//get histogram
		for (int j = 0 ; j < nowTiff.getNSlices() ; j++) { 
		
			IJ.showStatus("Adding to histogram in layer " + (j + 1));
			
			nowTiff.setPosition(j+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
						 
			for (int x = 0 ; x < nowTiff.getWidth() ; x++) {
				for (int y = 0 ; y < nowTiff.getHeight() ; y++) {
					double nowPix = nowIP.getPixelValue(x, y);
					
					//collect according for class boundaries	
					if (nowPix > 0) {
//						for (int i = 0 ; i < classBounds.length - 1 ; i++) { 
//							if (nowPix >= classBounds[i] & nowPix < classBounds[i+1]) value2PutIn = i; 
//						}
						int myEntry = (int)Math.round(nowPix);
						if (myEntry >= classBounds.length) myEntry = classBounds.length - 1;
						if (nowPix >= 0) psd[myEntry]++;
					}
				}
			}				
		}		
		
		//save histogram
		jIO.writePoreSizeDistribution(savePath, psd, classBounds);
	}

	public ImagePlus extractSpecificPoreLabels(ImagePlus idImp, ArrayList<Integer> percolatingClusters) {
		
		ImageStack outStack = new ImageStack(idImp.getWidth(), idImp.getHeight());
		
		//idImp.updateAndDraw();idImp.show();		
		
		//transform List to Array for better preformance
		IJ.showStatus("Sorting out detected air clusters..");
		Collections.sort(percolatingClusters);		
	
		for (int z = 1 ; z <= idImp.getNSlices() ; z++) {
			
			IJ.showStatus("Extracting respective clusters in slice " + z + "/" + idImp.getNSlices());
		
			ImageProcessor binIP = new ByteProcessor(idImp.getWidth(), idImp.getHeight());
			idImp.setSlice(z);			
			ImageProcessor idIP = idImp.getProcessor();
			
			for (int x = 0 ; x < binIP.getWidth() ; x++) {
				
				for (int y = 0 ; y < binIP.getWidth() ; y++) {
					
					int id = (int)Math.round(idIP.getPixelValue(x, y));
					if (percolatingClusters.contains(id)) binIP.putPixel(x, y, 255);
					else binIP.putPixel(x, y, 0);
					
				}
				
			}
			
			outStack.addSlice(binIP); 
			
		}
		
		ImagePlus outTiff = new ImagePlus();
		outTiff.setStack(outStack);
		
		return outTiff;
		
	}
	
	public PoreClusterProps calculateChiOfPercolation(ImagePlus nowTiff, InputOutput.MyFileCollection mFC, PoreClusterProps mCP) {
		
			int i,j;
			int x,y,z;
			
			double dSum = 0;
			double gSumSum = 0;
			
			double[] connectedCorrelationLength = new double[mCP.id.length];
			double numerator = 0;
			double gSum = 0;
			double[] chi = new double[mCP.id.length];
			
			IJ.showStatus("Calculating Chi of percolation ...");
			IJ.freeMemory();IJ.freeMemory();
			
			//sweep through the clusters
			int numberOfConsideredClusters = mCP.id.length;  //all clusters would be  mCP.id.length			
			for (i = 0 ; i < numberOfConsideredClusters ; i++) {			 
				
				int nowCluster = mCP.id[i];
					
				if (mCP.volume[i] == 1) {
					
					connectedCorrelationLength[i] = 0;
					numerator = 0;
					gSum = 0;
					
				} else if (!mCP.isPercolating[i]) {
					
					//long startTime = System.currentTimeMillis();		
					
					double rmax = 0;		//maximally possible distance between two voxels		
					
					//find all voxels belonging to the cluster	
					List<Integer> nX = new ArrayList<Integer>();
					List<Integer> nY = new ArrayList<Integer>();
					List<Integer> nZ = new ArrayList<Integer>();
					
					//define bounding box around the cluster 
					int x0 = mCP.minX[i];
					int y0 = mCP.minY[i];
					int z0 = mCP.minZ[i];
					int xe = mCP.maxX[i];
					int ye = mCP.maxY[i];
					int ze = mCP.maxZ[i];
					
					double xy = Math.sqrt((xe - x0) * (xe - x0) + (ye - y0) * (ye - y0));
					rmax = Math.sqrt((ze - z0) * (ze - z0) + xy * xy);
											
					for (z = z0 ; z <= ze ; z++) { 
						nowTiff.setPosition(z + 1);
						ImageProcessor nowIP = nowTiff.getProcessor();
						
						for (x = x0 ; x <= xe ; x++) {
							for (y = y0 ; y <= ye ; y++) {
								int nowPixel = (int)nowIP.getPixelValue(x, y);
								if (nowPixel == nowCluster) { 
									nX.add(x);
									nY.add(y);
									nZ.add(z);
								}
							}
						}					
					}				
							
					//give a warning if the cluster size is too large
					if (nX.size() > Integer.MAX_VALUE) IJ.error("cluster size exceeds Java Integer range.. arrgghhh...");
					
					//rewrite listed coordinates as array
					int nowX[] = new int[nX.size()];
					int nowY[] = new int[nX.size()];
					int nowZ[] = new int[nX.size()];
					for(j = 0 ; j < nX.size() ; j++) {
						nowX[j] = nX.get(j);
						nowY[j] = nY.get(j);
						nowZ[j] = nZ.get(j);
					}
					
					//if cluster size is larger than w * rmax take a random subset of all voxels below this number
					int numOfConsideredVoxels = nowX.length;
					int requiredMultipleOfRMax = 100;
					ArrayList<Integer> retain = new ArrayList<Integer>();
					if (numOfConsideredVoxels > requiredMultipleOfRMax * rmax) {
						Random randomGenerator = new Random();					
						while (retain.size() < (int)Math.floor(requiredMultipleOfRMax * rmax) + 1) {						
							int newNum = randomGenerator.nextInt(nowX.length);
							if (!retain.contains(newNum)) retain.add(newNum);
						}
						numOfConsideredVoxels = (int)Math.floor(requiredMultipleOfRMax * rmax);										
					}
					
					//calculate distances	
					int iRmax = (int)Math.ceil(rmax) + 1;
					int[] myHist = new int[iRmax];
					for (j = 0 ; j < myHist.length ; j++) myHist[j] = 0;						
					for (int a = 0 ; a < numOfConsideredVoxels ; a++) {
						for (int b = 0 ; b < numOfConsideredVoxels ; b++) {
							int ia = a; int ib = b;
							if (!retain.isEmpty()) {
								ia = retain.get(a);
								ib = retain.get(b);
							}
							int dx = nowX[ia] - nowX[ib];
							int dy = nowY[ia] - nowY[ib];
							int dz = nowZ[ia] - nowZ[ib];	
							double dxy = Math.sqrt(dx * dx + dy *dy);
							double dist = Math.sqrt(dxy * dxy + dz * dz);
							myHist[(int)Math.round(dist)]++;											
						}					
					}
									
					//normalize everything
					double[] r = new double[iRmax];
					for (j = 0 ; j < iRmax ; j++) r[j] = (j + 0.5); 
					double[] g = new double[iRmax];
					for (j = 0 ; j < iRmax; j++) g[j] = myHist[j] / (mCP.volume[i] - 1);
					
					//calculate pseudo chi
					double numer = 0;
					double denom = 0;
					for (j = 0 ; j < myHist.length ; j++) numer += r[j] * r[j] * g[j];
					for (j = 0 ; j < myHist.length ; j++) denom += g[j];
					connectedCorrelationLength[i] = Math.sqrt(numer / denom);	
					numerator = numer;
					gSum = denom;
											
	/*				float interval = (System.currentTimeMillis() - startTime) / 1000;		
					DecimalFormat df = new DecimalFormat();
					df.setMaximumFractionDigits(2);		
					IJ.error("Calculation time for " + numOfConsideredVoxels + " voxels was " + df.format(interval) + " seconds ..");*/
					
				} else {
					
					connectedCorrelationLength[i] = -1;
					numerator = 0;
					gSum = 0;
					
				}
				
				//calculate chi			
				dSum += numerator;
				gSumSum += gSum;			
				chi[i] = Math.sqrt(dSum / gSumSum);		
			}	
			
			IJ.freeMemory();IJ.freeMemory();
			
			//assign values to PoreClusterProperties
			//mCP.connectedCorrelationLength = connectedCorrelationLength;
			//mCP.chiOfPercolation = chi;
			
			return mCP;		
		}

	public double[] calculateTheBulkSoilVolume(String nowGaugePath, ImagePlus soilSurface, MenuWaiter.ROISelectionOptions mRSO, InputOutput.MyFileCollection mFC) {
		
		InputOutput jIO = new InputOutput();
		RoiHandler rH  = new RoiHandler();
		
		double[] bulkVolume = new double[4];
		
		//read InnerCircle file
		ObjectDetector jOD = new ObjectDetector();
		ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();
		int versio = jIO.checkInnerCircleFileVersion(nowGaugePath);			
		if (versio == 0) jCO = jIO.readInnerCircleVer0(nowGaugePath);	
		else jCO = jIO.readInnerCircleVer1(nowGaugePath);
		
		PolygonRoi[] pRoi = rH.makeMeAPolygonRoiStack("inner", "tight", jCO, -mRSO.cutAwayFromWall);
		
		//assign soil top and bottom surfaces		
		ImageStack surStack = soilSurface.getStack();
		ImageProcessor topIP = surStack.getProcessor(1);		
		ImageProcessor botIP = surStack.getProcessor(2);
		
		//cut away outside of column of top ...
		topIP.setRoi(pRoi[0]);
		topIP.setColor(0);
		topIP.fillOutside(pRoi[0]);
		
		// ... and bottom
		botIP.setRoi(pRoi[0]);
		botIP.setColor(0);
		botIP.fillOutside(pRoi[pRoi.length - 1]);
				
		// count volume above the soil surface ...		
		double volumeAboveColumn = 0;
		for (int x = 0 ; x < topIP.getWidth() ; x++) {
			for (int y = 0 ; y < topIP.getHeight() ; y++) {
				int surPix = topIP.getPixel(x, y);
				volumeAboveColumn += surPix;
			}
		}
		
		// .. and below..
		double volumeBelowColumn = 0;
		for (int x = 0 ; x < botIP.getWidth() ; x++) {
			for (int y = 0 ; y < botIP.getHeight() ; y++) {
				int surPix = botIP.getPixel(x, y);
				volumeBelowColumn += surPix;
			}
		}
		
		//calculate volume of cut-out column (in case the parts close to the walls have been cut away)
		double columnVolume = 0;
		if (mRSO.heightOfRoi == 0) {			
			
			//calculate volume of complete column
			double volume2Walls = 0;
			double originalVolume = 0;
			for (int i = 0 ; i < jCO.heightOfColumn ; i++) {
				originalVolume += Math.round(Math.PI * jCO.innerMajorRadius[i] * jCO.innerMinorRadius[i]);						
			}
			volume2Walls = originalVolume - columnVolume;
			
			for (int i = 0 ; i < jCO.heightOfColumn ; i++) {
				columnVolume += Math.round(Math.PI * (jCO.innerMajorRadius[i] - mRSO.cutAwayFromWall) * (jCO.innerMinorRadius[i] - mRSO.cutAwayFromWall));
			}			
			bulkVolume[0] = columnVolume - volumeAboveColumn - volumeBelowColumn;
			bulkVolume[1] = volumeAboveColumn;
			bulkVolume[2] = volumeBelowColumn;
			bulkVolume[3] = volume2Walls;
		}
		else {
			for (int i = mFC.startSlice ; i < mFC.startSlice + mRSO.heightOfRoi ; i++) {
				columnVolume += Math.round(Math.PI * (jCO.innerMajorRadius[i] - mRSO.cutAwayFromWall) * (jCO.innerMinorRadius[i] - mRSO.cutAwayFromWall));
			}	
			bulkVolume[0] = columnVolume;
		}
		
		return bulkVolume;
		
	}

	public ArrayList<Integer> check4TouchingTheBottom(ImagePlus nowTiff, ImagePlus origTiff, ImagePlus surfTiff) {
	
		ArrayList<Integer> conBot = new ArrayList<Integer>();
		
		String pathSep = "/";
	
		InputOutput jIO = new InputOutput(); 
		
		//get cluster image			
		int i,x,y;
		int iW = nowTiff.getWidth();
		int iH = nowTiff.getHeight();		
		int stackHeight = nowTiff.getNSlices();
		
		//define neighborhood variables
		int n;
		
		if (surfTiff != null) {
			
			surfTiff.setPosition(2);
			ImageProcessor botIP = surfTiff.getProcessor();
			int cc = 0;
			
			for (x = 0 ; x < iW ; x++) {
				for (y = 0 ; y < iH ; y++) {
					
					cc++;
					//IJ.showStatus("Finding pore clusters connected to the bottom surface in kilo-pixel " + (cc / 1000) + "/" + (iW * iH / 1000));
					
					int myVox = botIP.getPixel(x, y);
					if (myVox > 0) {  //check if pixel is within the soil column outlines					
						for (int ix = -1 ; ix < 2 ; ix++) {  
							for (int iy = -1 ; iy < 2 ; iy++) {
								if (ix != 0 & iy != 0) {
									n = botIP.getPixel(x + ix, y + iy);  //get values of neighborhood 
									if (n > myVox) {
										for (i = 0 ; i < n - myVox ; i++) {
											nowTiff.setPosition(stackHeight - (myVox + i));
											ImageProcessor nowIP = nowTiff.getProcessor();
											int nPixel = (int)nowIP.getPixelValue(x + ix, y + iy) - 1;										
											if (nPixel > -1) if (!conBot.contains(nPixel)) conBot.add(nPixel);										
										}
									} 
									if (n > 0 & n < myVox) {
										for (i = 0 ; i < myVox - n ; i++) {
											nowTiff.setPosition(stackHeight - (myVox - i));
											ImageProcessor nowIP = nowTiff.getProcessor();
											int nPixel = (int)nowIP.getPixelValue(x + ix, y + iy) - 1;								
											if (nPixel > -1) if (!conBot.contains(nPixel)) conBot.add(nPixel);										
										}
									}
									if (n == myVox) {										
										nowTiff.setPosition(stackHeight - myVox);
										ImageProcessor nowIP = nowTiff.getProcessor();
										int nPixel = (int)nowIP.getPixelValue(x + ix, y + iy) - 1;								
										if (nPixel > -1) if (!conBot.contains(nPixel)) conBot.add(nPixel);
									}											
								}
							}
						}
					}
				}
			}	
		}
		
		if (surfTiff == null) {
			
			//check bottom of sub-sample
			nowTiff.setPosition(stackHeight);
			origTiff.setPosition(stackHeight);
			ImageProcessor nowIP = nowTiff.getProcessor();
			ImageProcessor oriIP = origTiff.getProcessor();
			int cc = 0;		
			for (x = 0 ; x < iW ; x++) {
				for (y = 0 ; y < iH ; y++) {				
					cc++;
					IJ.showStatus("Finding pore clusters connected to the bottom surface in kilo-pixel " + (cc / 1000) + "/" + (iW * iH / 1000));
					
					int nPixel = (int)nowIP.getPixelValue(x, y);
					int oPixel = (int)oriIP.getPixelValue(x, y);		
					if (oPixel != 0) if (!conBot.contains(nPixel)) conBot.add(nPixel);	
				}
			}	
		}
		
		return conBot;
	
	}

	public ArrayList<Integer> check4TouchingTheTop(ImagePlus nowTiff, ImagePlus origTiff, ImagePlus surfTiff) {
	
		ArrayList<Integer> conTop = new ArrayList<Integer>();
		
		String pathSep = "/";
		
		//get cluster image			
		int i,x,y;
		int iW = nowTiff.getWidth();
		int iH = nowTiff.getHeight();		
		
		//define neighborhood variables
		int n;
		
		if (surfTiff != null) {
			
			surfTiff.setPosition(1);
			ImageProcessor surIP = surfTiff.getProcessor();
			int cc = 0;
			
			for (x = 0 ; x < iW ; x++) {
				for (y = 0 ; y < iH ; y++) {
				
					cc++;
					//IJ.showStatus("Finding pore clusters connected to the top surface in kilo-pixel " + (cc / 1000) + "/" + (iW * iH / 1000));
					
					int myVox = surIP.getPixel(x, y);
					if (myVox > 0) {  //check if pixel is within the soil column outlines					
						for (int ix = -1 ; ix < 2 ; ix++) {  
							for (int iy = -1 ; iy < 2 ; iy++) {
								if (ix != 0 & iy != 0) {
									n = surIP.getPixel(x + ix, y + iy);  //get values of neighborhood 
									if (n > myVox) {
										for (i = 0 ; i < n - myVox ; i++) {
											nowTiff.setPosition(myVox + i);
											ImageProcessor nowIP = nowTiff.getProcessor();
											int nPixel = (int)nowIP.getPixelValue(x + ix, y + iy) - 1;							
											if (nPixel > -1) if (!conTop.contains(nPixel)) conTop.add(nPixel);										
										}
									} 
									if (n > 0 & n < myVox) {
										for (i = 0 ; i < myVox - n ; i++) {
											nowTiff.setPosition(myVox - i);
											ImageProcessor nowIP = nowTiff.getProcessor();
											int nPixel = (int)nowIP.getPixelValue(x + ix, y + iy) - 1;								
											if (nPixel > -1) if (!conTop.contains(nPixel)) conTop.add(nPixel);										
										}
									}									
									if (n == myVox) {										
										nowTiff.setPosition(myVox);
										ImageProcessor nowIP = nowTiff.getProcessor();
										int nPixel = (int)nowIP.getPixelValue(x + ix, y + iy) - 1;								
										if (nPixel > -1) if (!conTop.contains(nPixel)) conTop.add(nPixel);
									}									
								}
							}
						}
					}
				}
			}
		}
		
		if (surfTiff == null) {
			
			int cc = 0;
			
			//check top of sub sample
			nowTiff.setPosition(1);
			origTiff.setPosition(1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			ImageProcessor oriIP = origTiff.getProcessor();
			for (x = 0 ; x < iW ; x++) {
				for (y = 0 ; y < iH ; y++) {
					cc++;
					IJ.showStatus("Finding pore clusters connected to the top surface in kilo-pixel " + (cc / 1000) + "/" + (iW * iH / 1000));
				
					int nPixel = (int)nowIP.getPixelValue(x, y);							
					int oPixel = (int)oriIP.getPixelValue(x, y);		
					if (oPixel != 0) if (!conTop.contains(nPixel)) conTop.add(nPixel);	
				}
			}
		}
		
		return conTop;
		
	}

	public ArrayList<Integer> checkForPercolatingClusters(int nParticles, ArrayList<Integer> conTop, ArrayList<Integer> conBot) {
			
			int i;		
			ArrayList<Integer> conPerk= new ArrayList<Integer>();
			
			//check if 
			if (nParticles < 1) {
				return conPerk; 
			}
							
			//sort Lists
			Collections.sort(conTop);
			Collections.sort(conBot);
			
			//checking if a cluster is percolating..
			//IJ.showStatus("Pinning down percolating clusters ...");
					
			//find percolating clusters			
	        for (int temp : conTop) {
	        	if (conBot.contains(temp)) conPerk.add(temp);
	        }
		
			return conPerk;	
		}
	
	public FractalProperties calculateFractalProperties(RoiHandler.ColumnRoi colRoi, MenuWaiter.ROISelectionOptions mRSO) {
		
		FractalProperties myFracs = new FractalProperties();
		HistogramStuff hist = new HistogramStuff();
		Dilate_ dil = new Dilate_();
		RoiHandler roi = new RoiHandler();
		ImageManipulator jIM = new ImageManipulator();
		DisplayThings disp = new DisplayThings();
		ImageManipulator.StackCalculator sC = jIM.new StackCalculator();
		FitStuff fit = new FitStuff();
		
		double[] l = {1, 2, 4, 8, 16};				//box edge length		
		//double[] v = {1, 8, 64, 512, 4096};			//box edge volumes
		double[] volBC = new double[l.length];		//volume box count
		double[] surfBC = new double[l.length];	//surface box count
		
		//make a copy of the fist image....
		ImagePlus copyTiff = colRoi.nowTiff.duplicate();
		
		for (int i = 0 ; i < l.length ; i++) {
	
			//scale if necessary
			if (l[i] > 1) {				
				copyTiff = jIM.binaryScale2HalfSize(copyTiff);
			}
			
			//copyTiff.updateAndDraw();
			//copyTiff.show();
			
			//pore volume		
			int[] myHist = hist.sampleHistogram(copyTiff);
			volBC[i] = myHist[255];
			
			//interface (surfaces)
			ImagePlus dilTiff = dil.dilate(copyTiff, 255, true);
			dilTiff = sC.subtract(dilTiff, copyTiff);
			
			//dilTiff.updateAndDraw();
			//dilTiff.show();
			
			if (colRoi.jCO != null) {				
				double[] oddVoxelContribution = {0, 0};
				if (Math.floorMod(copyTiff.getWidth(),2) > 0) oddVoxelContribution[0] = 1 / copyTiff.getWidth() / 2;
				if (Math.floorMod(copyTiff.getHeight(),2) > 0) oddVoxelContribution[1] = 1 / copyTiff.getHeight() / 2;				
								
				double[] iDims = {dilTiff.getWidth(), dilTiff.getHeight()};
				colRoi.jCO.heightOfColumn = dilTiff.getNSlices();
				ObjectDetector.ColCoords3D scaledCO = roi.scaleColumnCoordinates(iDims, colRoi.jCO, mRSO, 1 / l[i], oddVoxelContribution);
				PolygonRoi[] pRoi = roi.makeMeAPolygonRoiStack("inner", "tight", scaledCO, -1);					
				dilTiff = jIM.clearOutside(dilTiff, pRoi);
			}
			
			//dilTiff.updateAndDraw();
			//dilTiff.show();
			
			myHist = hist.sampleHistogram(dilTiff);
			surfBC[i] = myHist[255];
		}
		
		//logarithmize variables..
		double[] x = new double[l.length];
		double[] y1 = new double[l.length];
		double[] y2 = new double[l.length];
		
		for (int i = 0 ; i < l.length ; i++) {
			x[i] = Math.log10(l[i]);
			y1[i] = Math.log10(volBC[i]);
			y2[i] = Math.log10(surfBC[i]);
		}			
			
		//fit a line and calculate the fractal dimension
		//FitStuff.LinearFunction volFracFit = fit.fitLinearFunction(x, y1);
		//myFracs.volumeFractalDim = -volFracFit.slope;	
		
		FitStuff.LinearFunction surfFracFit = fit.fitLinearFunction(x, y2);
		myFracs.surfaceFractalDim = -surfFracFit.slope;
		
		//disp.plotXYXY(x, y1, x, y2, "fractal props", "log10 scale (vx^2)", "log10 box count");		
		
		copyTiff.unlock();copyTiff.flush();		
		IJ.freeMemory();
		
		return myFracs;
		
	}
	
	public AnisotropyResults calculateAnisotropy(boolean roiIsCylinder, ImagePlus nowTiff, double area) {
		 
		/////////////////////////////////////////////////////////////////////
		// OUTSIDE ROI VOXEL STATISTICS WORK ONLY FOR FAIRLY CIRCULAR COLUMNS //
		/////////////////////////////////////////////////////////////////////		

		TailoredMaths maths = new TailoredMaths();
		
		AnisotropyResults aRe = new AnisotropyResults();
		
		//init counts
		int xl = nowTiff.getWidth();
		int yl = nowTiff.getHeight();
		int zl = nowTiff.getNSlices();
			
		int[][] interCountZ = new int[xl][yl];
		int[][] interCountX = new int[zl][yl];
		int[][] interCountY = new int[xl][zl];
		int[][] interCountXY = new int[xl][yl];
		int[][] interCountYX = new int[xl][yl];

		int[][] interCountDXZ = new int[xl][yl];
		int[][] interCountDYZ = new int[xl][yl];
		int[][] interCountDXYZ = new int[xl][yl];
		int[][] interCountDYXZ = new int[xl][yl];
		
		int[][] interCountUXZ = new int[xl][yl];
		int[][] interCountUYZ = new int[xl][yl];
		int[][] interCountUXYZ = new int[xl][yl];
		int[][] interCountUYXZ = new int[xl][yl];
		
		int oldPixel = 0;
		int uldPixel = 0;
		
		ImageProcessor nowIP = null;
		ImageProcessor oldIP = null;
		ImageProcessor uldIP = null;
		
/*		double xCent = 0;
		double yCent = 0;
		double radius = 0;
		int outSideCount = 0;*/
			
		//count intercepts in the vertical and horizontal directions
		int voxelCounter = 0;
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {		
			
			IJ.showStatus("Analyzing anisotropies in layer " + (z + 1) + " of " + nowTiff.getNSlices());
			
			nowTiff.setSlice(z+1);
			ImageProcessor newIP = nowTiff.getProcessor().duplicate();
			
			if (z > 0) {
				nowTiff.setSlice(z);
				oldIP = nowTiff.getProcessor().duplicate();
			}
			
			if (z < nowTiff.getNSlices() - 1) {
				nowTiff.setSlice(z + 2);
				uldIP = nowTiff.getProcessor().duplicate();
			}
							
			//calculate center of ROI	
			double radius = 0;
			double xCent = 0;
			double yCent = 0;
	
			if (roiIsCylinder) {			
				xCent = nowTiff.getWidth() / 2;
				yCent = nowTiff.getHeight() / 2;				
				radius = (xCent+yCent) / 2;
			}
						
			boolean isWithinRoi = true;						
			for (int x = 0 ; x < xl ; x++) {
				for (int y = 0 ; y < yl ; y++) {
					
					//check if pixel is within ROI... 
					if (radius > 0) {
						double nowxDist = x - xCent;						
						double nowyDist = y - yCent;
						
						double radialCoord = Math.sqrt(nowxDist * nowxDist + nowyDist * nowyDist);
						if (radialCoord < radius - 1) isWithinRoi = true;
						else isWithinRoi = false;
					}
					
					//sample neighborhood and count
					if (isWithinRoi) {
						
						voxelCounter++;
						int nowPixel = newIP.getPixel(x, y);
						
						//vertical
						if (z > 0) oldPixel = oldIP.getPixel(x, y);
						else oldPixel = nowPixel;
						if (nowPixel != oldPixel) {
							interCountZ[x][y]++;
						}
	
						//horizontal X
						if (x > 0) oldPixel = newIP.getPixel(x - 1, y);
						else  oldPixel = nowPixel;
						if (nowPixel != oldPixel) interCountX[z][y]++;
						
						//horizontal Y
						if (y > 0) oldPixel = newIP.getPixel(x, y - 1);
						else  oldPixel = nowPixel;
						if (nowPixel != oldPixel) interCountY[x][z]++;
						
						//horizontal XY
						if (x > 0 & y > 0) oldPixel = newIP.getPixel(x - 1, y - 1);
						else  oldPixel = nowPixel;
						if (nowPixel != oldPixel) interCountXY[x][y]++;
						
						//horizontal YX
						if (x > 0 & y < yl - 1) oldPixel = newIP.getPixel(x - 1, y + 1);
						else  oldPixel = nowPixel;
						if (nowPixel != oldPixel) interCountYX[x][y]++;
						
						//skew 4 ones downwards  
						//XZ
						if (z > 0 & x > 0) oldPixel = oldIP.getPixel(x - 1, y);
						else  oldPixel = nowPixel;
						if (nowPixel != oldPixel) interCountDXZ[x][y]++;
						
						//YZ
						if (z > 0 & y > 0) oldPixel = oldIP.getPixel(x, y - 1);
						else  oldPixel = nowPixel;
						if (nowPixel != oldPixel) interCountDYZ[x][y]++;
						
						//XYZ
						if (z > 0 & x > 0 & y > 0) oldPixel = oldIP.getPixel(x - 1, y - 1);
						else  oldPixel = nowPixel;
						if (nowPixel != oldPixel) interCountDXYZ[x][y]++;
						
						//YXZ
						if (z > 0 & x > 0 & y < yl - 1) oldPixel = oldIP.getPixel(x - 1, y + 1);
						else  oldPixel = nowPixel;
						if (nowPixel != oldPixel) interCountDYXZ[x][y]++;
						
						//skew 4 ones upwards  
						//XZ
						if (z < nowTiff.getNSlices() - 1 & x > 0) oldPixel = uldIP.getPixel(x - 1, y);
						else  oldPixel = nowPixel;
						if (nowPixel != oldPixel) interCountUXZ[x][y]++;
						
						//YZ
						if (z < nowTiff.getNSlices() - 1 & y > 0) oldPixel = uldIP.getPixel(x, y - 1);
						else  oldPixel = nowPixel;
						if (nowPixel != oldPixel) interCountUYZ[x][y]++;
						
						//XYZ
						if (z < nowTiff.getNSlices() - 1 & x > 0 & y > 0) oldPixel = uldIP.getPixel(x - 1, y - 1);
						else  oldPixel = nowPixel;
						if (nowPixel != oldPixel) interCountUXYZ[x][y]++;
						
						//YXZ
						if (z < nowTiff.getNSlices() - 1 & x > 0 & y < yl - 1) oldPixel = uldIP.getPixel(x - 1, y + 1);
						else  oldPixel = nowPixel;
						if (nowPixel != oldPixel) interCountUYXZ[x][y]++;
					}
					
				}				
			}	
		}
		
		//correction function for anisotropy:		
			
		
		//calculate anisotropy indices
		aRe.z = maths.compileStatistics(interCountZ, xl, yl, 1);
		aRe.x = maths.compileStatistics(interCountX, zl, yl, 1);
		aRe.y = maths.compileStatistics(interCountY, xl, zl, 1);
		aRe.xy = maths.compileStatistics(interCountXY, xl, yl, Math.sqrt(2));
		aRe.yx = maths.compileStatistics(interCountYX, xl, yl, Math.sqrt(2));
		
		aRe.dxz = maths.compileStatistics(interCountDXZ, xl, yl, Math.sqrt(2));
		aRe.dyz = maths.compileStatistics(interCountDYZ, xl, yl, Math.sqrt(2));
		aRe.dxyz = maths.compileStatistics(interCountDXYZ, xl, yl, Math.sqrt(3));
		aRe.dyxz = maths.compileStatistics(interCountDYXZ, xl, yl, Math.sqrt(3));
		
		aRe.uxz = maths.compileStatistics(interCountUXZ, xl, yl, Math.sqrt(2));
		aRe.uyz = maths.compileStatistics(interCountUYZ, xl, yl, Math.sqrt(2));
		aRe.uxyz = maths.compileStatistics(interCountUXYZ, xl, yl, Math.sqrt(3));
		aRe.uyxz = maths.compileStatistics(interCountUYXZ, xl, yl, Math.sqrt(3));
		
		//calculate intermediate statistics..
		double[] aniSums = {aRe.x[0], aRe.y[0], aRe.z[0], aRe.xy[0], aRe.yx[0], aRe.dxz[0] ,aRe.dyz[0],aRe.dxyz[0],aRe.dyxz[0],aRe.uxz[0],aRe.uyz[0],aRe.uxyz[0],aRe.uyxz[0]};
				
		double aniMin = StatUtils.min(aniSums);
		double aniMax = StatUtils.max(aniSums);

		//calc anisotropy
		aRe.anisotropy = (aniMax - aniMin) / (aniMax + aniMin);
			
		//get direction of main alignment
		int[] mainAli = {0, 0, 0};		
		for (int i = 0 ; i < aniSums.length ; i++) {
			
			if (aRe.x[0] == aniMin) {mainAli[0] = 1;mainAli[1] = 0;mainAli[2] = 0;}
			if (aRe.y[0] == aniMin) {mainAli[0] = 0;mainAli[1] = 1;mainAli[2] = 0;}
			if (aRe.z[0] == aniMin) {mainAli[0] = 0;mainAli[1] = 0;mainAli[2] = 1;}
			if (aRe.xy[0] == aniMin) {mainAli[0] = 1;mainAli[1] = 1;mainAli[2] = 0;}
			if (aRe.yx[0] == aniMin) {mainAli[0] = 1;mainAli[1] = -1;mainAli[2] = 0;}
			if (aRe.dxz[0] == aniMin) {mainAli[0] = 1;mainAli[1] = 0;mainAli[2] = -1;}
			if (aRe.dyz[0] == aniMin) {mainAli[0] = 0;mainAli[1] = 1;mainAli[2] = -1;}
			if (aRe.dxyz[0] == aniMin) {mainAli[0] = 1;mainAli[1] = 1;mainAli[2] = -1;}
			if (aRe.dyxz[0] == aniMin) {mainAli[0] = 1;mainAli[1] = -1;mainAli[2] = -1;}
			if (aRe.uxz[0] == aniMin) {mainAli[0] = 1;mainAli[1] = 0;mainAli[2] = 1;}
			if (aRe.uyz[0] == aniMin) {mainAli[0] = 0;mainAli[1] = 1;mainAli[2] = 1;}
			if (aRe.uxyz[0] == aniMin) {mainAli[0] = 1;mainAli[1] = 1;mainAli[2] = 1;}
			if (aRe.uyxz[0] == aniMin) {mainAli[0] = 1;mainAli[1] = -1;mainAli[2] = 1;}
			
		}		
		
		aRe.mainAlignment = mainAli;
		
		return aRe;
	}
	
	public double calculatePercolationThreshold(ImagePlus nowTiff, int maxVol, ImagePlus surfTiff, InputOutput.MyFileCollection mFC) {
		
		Dilate_ mD = new Dilate_();
		
		double percolationThreshold = 0;
		
		boolean laPerco = false;
		int cc = 0;
		
		while (!laPerco) {
		
			nowTiff = mD.dilate(nowTiff, 255, false);		
		
			//set percolatioin threshold..
			percolationThreshold = macroPoreVolume(nowTiff);
			
			//find pore clusters in laPores and check if they form a percolating path		
			FloodFillComponentsLabeling3D myFCCL = new FloodFillComponentsLabeling3D(26, 32);
			ImageStack myLabelStack = myFCCL.computeLabels(nowTiff.getStack());						
			ImagePlus laParLaTiff = new ImagePlus();
			laParLaTiff.setStack(myLabelStack);
			
			int[] labels = LabelImages.findAllLabels(myLabelStack);
			int numOfObjects = labels.length;
			
			cc++;
			IJ.showStatus("Finding critical pore diameter, iteration " + cc + " ...");
			
			ArrayList<Integer> laConTop = check4TouchingTheTop(laParLaTiff, nowTiff, surfTiff);	
			ArrayList<Integer> laConBot = check4TouchingTheBottom(laParLaTiff, nowTiff, surfTiff);
			ArrayList<Integer> percolatingClusters = checkForPercolatingClusters(numOfObjects, laConTop, laConBot);
			
			if (!percolatingClusters.isEmpty()) {
				laPerco = true;
			}			
		}
		
		return percolationThreshold;
		
	}
	
	public double[] calculateCriticalPoreDiameter(ImagePlus distTiff, int maxVol, ImagePlus surfTiff, InputOutput.MyFileCollection mFC) {
	
		double criticalPoreRadius = 0;
		
		//find the smallest of the maximal pore diameters of all cross-sections and use it as starting point
		ProfileStatistics pS = findProfileStatistics(distTiff, null);				
		double minVertThick = StatUtils.min(pS.maxi);
		
		//set threshold to pore diameters == minVertThick 
		double upperthreshold = minVertThick;
		double threshold = upperthreshold;
		double lowerthreshold = 1;
		boolean laPerco = false;
		boolean notTheFirstTrialAnymore = false;
		
		//also remember porosity at percolation
		double percolationThreshold = 0;
				
		int cc = 0;
		//loop while laPerco is false
		while (!laPerco) {
			
			cc++;
			
			IJ.showStatus("Finding critical pore diameter, iteration " + cc + " ...");
			
			//set all values smaller than the threshold to 0
			ImageStack largerPores = new ImageStack(distTiff.getWidth(), distTiff.getHeight());
			ImagePlus laPores = new ImagePlus();
			for (int z = 1 ; z <= distTiff.getNSlices() ; z++) {	
				
				distTiff.setPosition(z);
				ImageProcessor nowIP = distTiff.getProcessor();
				ImageProcessor laIP = new ByteProcessor(distTiff.getWidth(), distTiff.getHeight());
				
				//loop through all pixel within cross-section
				for (int x = 0 ; x < distTiff.getWidth() ; x++) {
					for (int y = 0 ; y < distTiff.getHeight() ; y++) {							
						double nowPixel = nowIP.getPixelValue(x, y);
						if (nowPixel >= threshold) laIP.putPixel(x, y, 255);
						else laIP.putPixel(x, y, 0);							
					}
				}					
				largerPores.addSlice(laIP);
			}		
			laPores.setStack(largerPores);
			
			//set percolatioin threshold..
			percolationThreshold = macroPoreVolume(laPores);
			
			//laPores.updateAndDraw();
			//laPores.show();
							
			IJ.freeMemory();IJ.freeMemory();
					
			//find pore clusters in laPores and check if they form a percolating path		
			FloodFillComponentsLabeling3D myFCCL = new FloodFillComponentsLabeling3D(26, 32);
			ImageStack myLabelStack = myFCCL.computeLabels(laPores.getStack());						
			ImagePlus laParLaTiff = new ImagePlus();
			laParLaTiff.setStack(myLabelStack);
			
			int[] labels = LabelImages.findAllLabels(myLabelStack);
			int numOfObjects = labels.length;
			
			IJ.showStatus("Finding critical pore diameter, iteration " + cc + " ...");
			
			ArrayList<Integer> laConTop = check4TouchingTheTop(laParLaTiff, laPores, surfTiff);	
			ArrayList<Integer> laConBot = check4TouchingTheBottom(laParLaTiff, laPores, surfTiff);
			ArrayList<Integer> percolatingClusters = checkForPercolatingClusters(numOfObjects, laConTop, laConBot);
			
			//laParLaTiff.updateAndDraw();
			//laParLaTiff.show();
			
			//check if there is a percolating cluster and try again with a smaller threshold in case of not
			if (percolatingClusters.isEmpty()) {
				if (upperthreshold - lowerthreshold < 0.2) {
					criticalPoreRadius = threshold;
					laPerco = true;  // this will break the loop..
				}
				else {
					if (!notTheFirstTrialAnymore) notTheFirstTrialAnymore = true;
					else upperthreshold = threshold;
					threshold = lowerthreshold + (upperthreshold - lowerthreshold) / 2;					 
				}
			}
			else {
				if (upperthreshold - lowerthreshold < 0.2 | !notTheFirstTrialAnymore) { 
					laPerco = true;  // this will break the loop..
					criticalPoreRadius = threshold;
				}
				else {
					lowerthreshold = threshold;
					threshold = lowerthreshold + (upperthreshold - lowerthreshold) / 2;	;
				}
			}
		}
		
		double[] results = {2 * criticalPoreRadius, percolationThreshold};
		
		return results;
				
	}
	
	public ImagePlus extractPercolatingThicknesses(ImagePlus thickTiff, ImagePlus labelTiff, ArrayList<Integer> percolatingClusters, int startZ) {
			
		ArrayList<Double> poreDiameters = new ArrayList<Double>();		//unique pore diameters
	
		//init image of just the percolating clusters..
		ImagePlus thickPC = new ImagePlus();
		ImageStack thickPCStack = new ImageStack(thickTiff.getWidth(), thickTiff.getHeight());
			
		//loop through horizontal cross-sections
		for (int z = startZ ; z <= thickTiff.getNSlices() ; z++) {				
			
			thickTiff.setPosition(z);
			ImageProcessor nowIP = thickTiff.getProcessor();
			
			//thickImp.updateAndDraw();
			//thickImp.show();
			
			labelTiff.setPosition(z);
			ImageProcessor clustIP = labelTiff.getProcessor();
			ImageProcessor outIP = nowIP.duplicate();
			
			//loop through all pixel within cross-section
			for (int x = 0 ; x < thickTiff.getWidth() ; x++) {
				for (int y = 0 ; y < thickTiff.getHeight() ; y++) {
					
					//check if pixel is contained in a percolating cluster
					int nowCluster = (int)clustIP.getPixelValue(x, y);
					if (percolatingClusters.contains(nowCluster)) {																
						
						double nowPixel = nowIP.getPixelValue(x, y);
						outIP.putPixelValue(x, y, nowPixel); 
															
						//make an inventory of existing pore diameters..
						boolean isContained = false;
						double myRoundPix = Math.round(nowPixel);
						if (!poreDiameters.isEmpty()) {
							for (int j = 0 ; j < poreDiameters.size() ; j++) if (poreDiameters.get(j) == myRoundPix) isContained = true;										
						}
						if (!isContained) poreDiameters.add(myRoundPix);									
					}
					else outIP.putPixelValue(x, y, 0);   //cut away if not within percolating cluster							
				}							
			}
			
			thickPCStack.addSlice(outIP);
		}
		
		//sort the list and create a new image of just the thicknesses of the percolating pore clusters.. 
		Collections.sort(poreDiameters);
		thickPC.setStack(thickPCStack);		
		
		//thickPC.updateAndDraw();
		//thickPC.show();
		
		IJ.freeMemory();IJ.freeMemory();

		return thickPC;
		
	}
		
	
	public ImagePlus extractPercolatingDistances(ImagePlus distTiff, ImagePlus labelTiff, ArrayList<Integer> percolatingClusters, int startZ) {
		
		ArrayList<Double> poreDiameters = new ArrayList<Double>();		//unique pore diameters
	
		//init image of just the percolating clusters..
		ImagePlus thickPC = new ImagePlus();
		ImageStack thickPCStack = new ImageStack(distTiff.getWidth(), distTiff.getHeight());
			
		//loop through horizontal cross-sections
		for (int z = startZ ; z <= distTiff.getNSlices() ; z++) {				
			
			distTiff.setPosition(z);
			ImageProcessor nowIP = distTiff.getProcessor();
			
			//thickImp.updateAndDraw();
			//thickImp.show();
			
			labelTiff.setPosition(z);
			ImageProcessor clustIP = labelTiff.getProcessor();
			ImageProcessor outIP = nowIP.duplicate();
			
			//loop through all pixel within cross-section
			for (int x = 0 ; x < distTiff.getWidth() ; x++) {
				for (int y = 0 ; y < distTiff.getHeight() ; y++) {
					
					//check if pixel is contained in a percolating cluster
					int nowCluster = (int)clustIP.getPixelValue(x, y);
					if (percolatingClusters.contains(nowCluster)) {																
						
						double nowPixel = nowIP.getPixelValue(x, y);
						outIP.putPixelValue(x, y, nowPixel); 
															
						//make an inventory of existing pore diameters..
						boolean isContained = false;
						double myRoundPix = Math.round(nowPixel);
						if (!poreDiameters.isEmpty()) {
							for (int j = 0 ; j < poreDiameters.size() ; j++) if (poreDiameters.get(j) == myRoundPix) isContained = true;										
						}
						if (!isContained) poreDiameters.add(myRoundPix);									
					}
					else outIP.putPixelValue(x, y, 0);   //cut away if not within percolating cluster							
				}							
			}
			
			thickPCStack.addSlice(outIP);
		}
		
		//sort the list and create a new image of just the thicknesses of the percolating pore clusters.. 
		Collections.sort(poreDiameters);
		thickPC.setStack(thickPCStack);		
		
		//thickPC.updateAndDraw();
		//thickPC.show();
		
		IJ.freeMemory();IJ.freeMemory();

		return thickPC;
		
	}
	
	public class MLJParticles {
		
		public ImagePlus myPartyPic;
		public double nParticles;
		public double[] volumes;	
		public int[] xMin;
		public int[] yMin;
		public int[] zMin;
		public int[] xMax;
		public int[] yMax;
		public int[] zMax;
		double[][] centroids;		
		public double[] loops;
		public double[] cavities;
		
		
	}
	
	
	public ProfileStatistics findProfileStatistics(ImagePlus nowTiff, PolygonRoi[] pRoi) {
		
		ProfileStatistics pS = new ProfileStatistics();		
		int[] numberOfNonZeros = new int[nowTiff.getNSlices()];
		double[] mean = new double[nowTiff.getNSlices()];	
		double[] geomean = new double[nowTiff.getNSlices()];	
		double[] logMean = new double[nowTiff.getNSlices()];
		double[] median = new double[nowTiff.getNSlices()];
		double[] mode = new double[nowTiff.getNSlices()];
		double[] std = new double[nowTiff.getNSlices()];
		double[] minVal = new double[nowTiff.getNSlices()];
		double[] maxVal = new double[nowTiff.getNSlices()];	
		
		for (int i = 1 ; i <= nowTiff.getNSlices() ; i++) {
			
			IJ.showStatus("Finding statistics of along vertical profile " + i + "/" + nowTiff.getNSlices() + "...");
		
			//move to next slice
			nowTiff.setPosition(i);
			ImageProcessor nowIP = nowTiff.getProcessor();
		
			//sample all grey values in a float array list
			ArrayList<Float> myValues = new ArrayList<Float>();			
			for (int x = 0 ; x < nowIP.getWidth() ; x++) {
				for (int y = 0 ; y < nowIP.getHeight() ; y++) {
					if (nowIP.getPixelValue(x, y) > 0) myValues.add(nowIP.getPixelValue(x, y));										
				}
			}
		
			//convert the float array list to an array
			double[] greyValues = new double[myValues.size()];
			for (int j = 0 ; j < myValues.size(); j++) greyValues[j] = myValues.get(j); 
			
			//calculate statistical values from the float array
			numberOfNonZeros[i-1] = myValues.size();
			mean[i-1] = StatUtils.mean(greyValues);
			logMean[i-1] = StatUtils.sumLog(greyValues)/numberOfNonZeros[i-1];
			median[i-1] = StatUtils.percentile(greyValues,50);			
			geomean[i-1] = StatUtils.geometricMean(greyValues);			
			std[i-1] = Math.sqrt(StatUtils.variance(greyValues));
			maxVal[i-1] = StatUtils.max(greyValues);	
			minVal[i-1] = StatUtils.min(greyValues);	
		}	
		
		pS.numberOfNonZeros = numberOfNonZeros;
		pS.mean = mean;
		pS.median = median;
		pS.geomean = geomean;
		pS.mode = mode;
		pS.std = std;
		pS.maxi = maxVal;
		pS.mini = minVal;
		
		return pS;
	}
	
	public void tailoredPoreNetworkAnalyzer(int imageNumber, InputOutput.MyFileCollection mFC, RoiHandler.ColumnRoi colRoi, MenuWaiter.PoreSpaceAnalyzerOptions mPSA) {
		
		//try to free up some memory			
		IJ.freeMemory();IJ.freeMemory();
		
		String pathSep = "/";
		
		InputOutput jIO = new InputOutput();	
		ImageManipulator jIM = new ImageManipulator();
		Skeletonize3D_ jSK = new Skeletonize3D_();
		AnalyzeSkeleton_ jAS = new AnalyzeSkeleton_();
		
		//prepare image for skeletonization
		//remove holes from pore system;
		if (mPSA.removeHoles) {
			colRoi.nowTiff = jIM.removeHoles(colRoi.nowTiff);
			//colRoi.nowTiff.updateAndDraw();
			//colRoi.nowTiff.show();
		}
		
		//init some very basic variables..		
		String nowImageName = mFC.colName;
		
		//create distance map to superimpose skeleton
		float[] floatWeights = ChamferWeights3D.BORGEFORS.getFloatWeights();
		boolean normalize = true;
		DistanceTransform3D dt = new DistanceTransform3DFloat(floatWeights, normalize);
		ImageStack result = dt.distanceMap(colRoi.nowTiff.getStack());
		
		ImagePlus distTiff = new ImagePlus("dist",result);
		
		//add slices on top and bottom of column to be able to extract the backbone
		ImageManipulator.SkeletonizerOptions mSO = jIM.addSlices4Skeletonization(colRoi.nowTiff, mPSA, colRoi);
		ImagePlus nowTiff = mSO.nowTiff;
		
		//do skeletonization				
		jSK.setup("",nowTiff);
		nowTiff.setPosition((int)Math.round(nowTiff.getNSlices() / 2));
		ImageProcessor nowIP = nowTiff.getProcessor();
		jSK.run(nowIP);
				
		//fuse skeleton with distance Tiff
		distTiff = jIM.fuseSkeletonAndDistance(nowTiff, distTiff, mSO.numberOfSlicesAdded);
		
		///////////////////////////////////
		//analyze skeleton
		///////////////////////////////////
		
		//get basics
		int pruneIndex = 0;
		boolean pruneEnds = false;
		boolean shortPath = false;		
		boolean silent = true;
		boolean verbose = false;
		Roi roi = null;
		jAS.setup("", nowTiff.duplicate());
		SkeletonResult mSP = jAS.run(pruneIndex, pruneEnds, shortPath, nowTiff.duplicate(), silent, verbose, roi); 
		
		//calculate constriction factors
		PlusGraph[] mPG = getConstrictionFactors(mSP, distTiff, mSO.numberOfSlicesAdded);
	
		//find connections to top and bottom..
		ConnectingGraphs mCG = findConnectingGraphs(mPG, distTiff.getNSlices(), mSO.numberOfSlicesAdded);		
		
		//calculate shortest connection
		FloydWarshallReturn[] mFWR = new FloydWarshallReturn[mCG.percolating.size()];  
		for (int i = 0 ; i < mCG.percolating.size() ; i++) {
			try {
				mFWR[i] = runFloydWarshallAlgorithm(mPG[mCG.percolating.get(i)], distTiff.getNSlices(), mSO.numberOfSlicesAdded);
			}
			catch (Exception e) {
				mFWR[i] = null;
			}
		}
		
		//calculate tortuosity
		NetworkProps mNP = calculateTortuosityAndStuff (mPG, mFWR, mCG, distTiff.getNSlices(), mSO.numberOfSlicesAdded);
		
		//save results
		jIO.writeNetworkAnalysesResults(mFC, mNP);
		
	}
	
	public NetworkProps calculateTortuosityAndStuff (PlusGraph[] mPG, FloydWarshallReturn[] mFWR, ConnectingGraphs mCG, int numberOfLastSlice, int numberOfAddedSlices) {
			
		NetworkProps mNP = new NetworkProps();
		
		mNP.mFWR = mFWR; 
		
		int[] numberOfEdges = new int[mPG.length];
		int[] numberOfVertices = new int[mPG.length];
		int[] numberOfSlabs = new int[mPG.length];
		int[] numberOfEndPoints = new int[mPG.length];
		int[] numberOfVertexPoints = new int[mPG.length];
		int[] sumOfCoordinationNumbers = new int[mPG.length];	
		int[] numberOfTopVertices = new int[mPG.length];
		int[] numberOfBotVertices = new int[mPG.length];
		int singleVertexGraphs = 0;
		
		double[] sumOfEdgeLength = new double[mPG.length];
		double[] sumOfEuclidEdgeLength = new double[mPG.length];
		double[] meanTortuosity = new double[mPG.length];
		double[] meanAngle2Z = new double[mPG.length];
				
		for (int i = 0 ; i < mPG.length ; i++) {			
			
			ArrayList<PlusEdge> myEdges = mPG[i].plusEdge;
			ArrayList<Vertex> myVerts = mPG[i].getVertices();
						
			numberOfEdges[i] = myEdges.size(); 
			numberOfVertices[i] = myVerts.size();
			
			double[] consFacs = new double[numberOfEdges[i]];
			double[] edgeLengths = new double[numberOfEdges[i]];
			double[] euclidEdgeLengths = new double[numberOfEdges[i]];
			double[] tortuosity = new double[numberOfEdges[i]];
			double[] angle = new double[numberOfEdges[i]];
			int[] nowNumberOfSlabs = new int[numberOfEdges[i]];
			
			//go through all edges
			if (numberOfEdges[i] == 0) {
				
				numberOfSlabs[i] = 0;
				sumOfEdgeLength[i] = 0;
				sumOfEuclidEdgeLength[i] = 0;
				meanTortuosity[i] = 0;			
				meanAngle2Z[i] = 0;
				
				singleVertexGraphs++;
				
			}
			else {
				for (int j = 0 ; j < numberOfEdges[i] ; j++) {
			
					PlusEdge nowEdge = myEdges.get(j);
					
					consFacs[j] = nowEdge.constrictionFactor;
					nowNumberOfSlabs[j] = nowEdge.getSlabs().size();
					
					Vertex v1 = nowEdge.getV1();
					Vertex v2 = nowEdge.getV2();
					
					edgeLengths[j] = nowEdge.getLength() - nowEdge.lengthReduction;
					euclidEdgeLengths[j] = getEuclideanDistance(v1, v2, numberOfLastSlice, numberOfAddedSlices);
					
					if (euclidEdgeLengths[j] > edgeLengths[j] | euclidEdgeLengths[j] == 0) {
						tortuosity[j] = 0;
					}					
					else tortuosity[j] = edgeLengths[j] / euclidEdgeLengths[j];					
					
					if (euclidEdgeLengths[j] > edgeLengths[j]) euclidEdgeLengths[j] = edgeLengths[j]; 
					
					angle[j] = calculateAngleToZAxis(v1, v2); 
					
				}					
				
				numberOfSlabs[i] = TailoredMaths.intSum(nowNumberOfSlabs);
				sumOfEdgeLength[i] = StatUtils.sum(edgeLengths);
				sumOfEuclidEdgeLength[i] = StatUtils.sum(euclidEdgeLengths);
				meanTortuosity[i] = TailoredMaths.meanWithoutZeros(tortuosity);			
				meanAngle2Z[i] = StatUtils.mean(angle);
			}
			

			
			//go through all vertices
			int[] nowCoordNumber = new int[numberOfVertices[i]];
			int[] nowVoxelNumber = new int[numberOfVertices[i]];
			int nowEndPointnumber = 0;
			int nowTopNumber = 0;
			int nowBotNumber = 0;
			
			for (int j = 0 ; j < numberOfVertices[i] ; j++) {
			
				Vertex nowVert = myVerts.get(j);
				
				nowCoordNumber[j] = nowVert.getBranches().size();
				nowVoxelNumber[j] = nowVert.getPoints().size(); 			
				
				if (vertexIsConnected2TopAndBottom(nowVert, numberOfLastSlice, numberOfAddedSlices) == 1) nowTopNumber++;
				if (vertexIsConnected2TopAndBottom(nowVert, numberOfLastSlice, numberOfAddedSlices) == 2) nowBotNumber++;
				
				if (nowCoordNumber[j] == 1) nowEndPointnumber++;
				
			}
			
			numberOfEndPoints[i] = nowEndPointnumber;
			numberOfVertexPoints[i] = TailoredMaths.intSum(nowVoxelNumber);			
			sumOfCoordinationNumbers[i] = TailoredMaths.intSum(nowCoordNumber);
			numberOfTopVertices[i] = nowTopNumber; 
			numberOfBotVertices[i] = nowBotNumber;
			
		}
		
		//get numbers for percolating Graphs
		mNP.numberOfPercolatingTopPoints = 0;
		mNP.numberOfPercolatingBotPoints = 0;
		mNP.percolatingSlabs = 0;
		mNP.percolatingDeadEnds = 0;
		mNP.percolatingVertices = 0;
		mNP.percolatingEdges = 0;
		double percolatingEdgeLength = 0;
		double percolatingEuclidEdgeLength = 0;
		double sumOfPercolatingAngles = 0;
		int percolatingCN = 0;		
		if (mCG.percolating.size() > 0) {
			for (int i = 0 ; i < mCG.percolating.size() ; i++) {
			
				int perco = mCG.percolating.get(i);
			
				mNP.numberOfPercolatingTopPoints += numberOfTopVertices[perco];
				mNP.numberOfPercolatingBotPoints += numberOfBotVertices[perco];
				mNP.percolatingSlabs += numberOfSlabs[perco];				
				mNP.percolatingDeadEnds += numberOfEndPoints[perco];
				percolatingCN += sumOfCoordinationNumbers[perco];
				mNP.percolatingEdges += numberOfEdges[perco];
				mNP.percolatingVertices += numberOfVertices[perco];
				
				percolatingEdgeLength += sumOfEdgeLength[perco];
				percolatingEuclidEdgeLength += sumOfEuclidEdgeLength[perco];
				
				sumOfPercolatingAngles += meanAngle2Z[perco] * (double)numberOfEdges[perco];
				
			}
		}
				
		//get numbers for all graphs connected to the top
		mNP.numberOfCon2TopPoints = 0;
		mNP.con2TopSlabs = 0;
		mNP.con2TopDeadEnds = 0;
		mNP.con2TopVertices = 0;
		mNP.con2TopEdges = 0;
		double con2TopEdgeLength = 0;
		double con2TopEuclidEdgeLength = 0;
		double sumOfcon2TopAngles = 0;
		int con2TopCN = 0;		
		if (mCG.con2Top.size() > 0) {
			for (int i = 0 ; i < mCG.con2Top.size() ; i++) {
			
				int con2Top = mCG.con2Top.get(i);
			
				mNP.numberOfCon2TopPoints += numberOfTopVertices[con2Top];
				mNP.con2TopSlabs += numberOfSlabs[con2Top];				
				mNP.con2TopDeadEnds += numberOfEndPoints[con2Top];
				con2TopCN += sumOfCoordinationNumbers[con2Top];
				mNP.con2TopEdges += numberOfEdges[con2Top];
				mNP.con2TopVertices += numberOfVertices[con2Top];
				
				con2TopEdgeLength += sumOfEdgeLength[con2Top];
				con2TopEuclidEdgeLength += sumOfEuclidEdgeLength[con2Top];	
				
				sumOfcon2TopAngles += meanAngle2Z[con2Top] * (double)numberOfEdges[con2Top];
			}
		}
		
		//get numbers for all graphs connected to the bottom
		mNP.numberOfCon2BotPoints = 0;
		mNP.con2BotSlabs = 0;
		mNP.con2BotDeadEnds = 0;
		mNP.con2BotVertices = 0;
		mNP.con2BotEdges = 0;
		double con2BotEdgeLength = 0;
		double con2BotEuclidEdgeLength = 0;
		double sumOfcon2BotAngles = 0;
		int con2BotCN = 0;		
		if (mCG.con2Bot.size() > 0) {
			for (int i = 0 ; i < mCG.con2Bot.size() ; i++) {
			
				int con2Bot = mCG.con2Bot.get(i);
			
				mNP.numberOfCon2BotPoints += numberOfTopVertices[con2Bot];
				mNP.con2BotSlabs += numberOfSlabs[con2Bot];				
				mNP.con2BotDeadEnds += numberOfEndPoints[con2Bot];
				con2BotCN += sumOfCoordinationNumbers[con2Bot];
				mNP.con2BotEdges += numberOfEdges[con2Bot];
				mNP.con2BotVertices += numberOfVertices[con2Bot];
				
				con2BotEdgeLength += sumOfEdgeLength[con2Bot];
				con2BotEuclidEdgeLength += sumOfEuclidEdgeLength[con2Bot];				
				
				sumOfcon2BotAngles += meanAngle2Z[con2Bot] * (double)numberOfEdges[con2Bot];
			}
		}
		
		//also refine Floyd-Warshall results
		mNP.fastestPathResistance = Double.MAX_VALUE;		
		mNP.fastestPathBottleneck = 0;
		mNP.fastestPathEdgeTortuosity = 0;
		mNP.fastestPathGlobalTortuosity = 0;
		mNP.fastestPathVertices = 0;
		mNP.fastestPathSlabs = 0;
		mNP.fastestPathVerticesPerSlab = 0;;
		mNP.meanFastestPathCN = 0;   //CN == coordination number*/
		if (mFWR.length > 0) {
			for (int i = 0 ; i < mFWR.length ; i++) {
				if (mFWR[i] != null) {
					if (mFWR[i].shortestPath < mNP.fastestPathResistance) {
						mNP.fastestPathResistance = mFWR[i].shortestPath;
						mNP.fastestPathBottleneck = mFWR[i].shortest.bottleNeck;
						mNP.fastestPathEdgeTortuosity = mFWR[i].shortest.edgeTortuosity;
						mNP.fastestPathGlobalTortuosity = mFWR[i].shortest.globalTortuosity;
						mNP.fastestPathVertices = mFWR[i].shortest.numberOfVertices;
						mNP.fastestPathSlabs = mFWR[i].shortest.pathPoints.size() - mNP.fastestPathVertices;
						mNP.fastestPathVerticesPerSlab = mNP.fastestPathVertices / mNP.fastestPathSlabs;
						mNP.meanFastestPathCN = (double)mFWR[i].shortest.sumOfCN / (double)mNP.fastestPathVertices;
					}
				}
			}			
		}
		if (mNP.fastestPathResistance == Double.MAX_VALUE) mNP.fastestPathResistance = 0;
		
		//calculate final measures!
		mNP.globalEdges = TailoredMaths.intSum(numberOfEdges); 
		mNP.globalSlabs = TailoredMaths.intSum(numberOfSlabs);
		mNP.globalVertices = TailoredMaths.intSum(numberOfVertices);
		mNP.globalDeadEnds = TailoredMaths.intSum(numberOfEndPoints);
		mNP.numberOfSingleVertexGraphs = singleVertexGraphs;
		mNP.numberOfCon2TopPoints = TailoredMaths.intSum(numberOfTopVertices);
		mNP.numberOfCon2BotPoints = TailoredMaths.intSum(numberOfBotVertices);
		
		mNP.globalEdgeLength = StatUtils.sum(sumOfEdgeLength);
		mNP.globalEuclidEdgeLength = StatUtils.sum(sumOfEuclidEdgeLength);
		mNP.globalTortuosity = mNP.globalEdgeLength / mNP.globalEuclidEdgeLength;		
		mNP.globalMeanTortuosity = TailoredMaths.calculateWeightedAverage(meanTortuosity, numberOfEdges);
		mNP.globalMeanCN = (double)TailoredMaths.intSum(sumOfCoordinationNumbers) / (double)mNP.globalVertices; 		//CN == coordination number	
		mNP.globalMeanAngleFromZ = TailoredMaths.calculateWeightedAverage(meanAngle2Z, numberOfEdges);
		
		mNP.percolatingTortuosity = percolatingEdgeLength / percolatingEuclidEdgeLength;
		if (mNP.percolatingVertices == 0) mNP.percolatingMeanCN = 0;
		else mNP.percolatingMeanCN = (double)percolatingCN / (double)mNP.percolatingVertices; 		//CN == coordination number		
		mNP.percolatingMeanAngleFromZ = sumOfPercolatingAngles / (double)mNP.percolatingEdges;
		
		mNP.con2TopTortuosity = con2TopEdgeLength / con2TopEuclidEdgeLength;	
		if (mNP.con2TopVertices == 0) mNP.con2TopMeanCN = 0;
		else mNP.con2TopMeanCN = (double)con2TopCN / (double)mNP.con2TopVertices; 		//CN == coordination number		
		mNP.con2TopMeanAngleFromZ = sumOfcon2TopAngles / (double)mNP.con2TopEdges;
		
		mNP.con2BotTortuosity = con2BotEdgeLength / con2BotEuclidEdgeLength;
		if (mNP.con2BotVertices == 0) mNP.con2BotMeanCN = 0;
		else mNP.con2BotMeanCN = (double)con2BotCN / (double)mNP.con2BotVertices; 
		mNP.con2BotMeanAngleFromZ = sumOfcon2BotAngles / (double)mNP.con2BotEdges;
				
		return mNP;
		
	}
	
	public double getEuclideanDistance(Vertex v1, Vertex v2, int numberOfLastSlice, int numberOfAddedSlices) {
		
		double euclidDist = 0;
		
		Point v1c = meanVertexLocation(v1);
		Point v2c = meanVertexLocation(v2);		
		
		//check if one or both of the points is outside the proper soil volume
		boolean v1IsOutsideSoil = false;
		boolean v2IsOutsideSoil = false;
		if (v1c.z < numberOfAddedSlices + 1 | v1c.z >= numberOfLastSlice - numberOfAddedSlices) v1IsOutsideSoil = true;
		if (v2c.z < numberOfAddedSlices + 1 | v2c.z >= numberOfLastSlice - numberOfAddedSlices) v2IsOutsideSoil = true;
		
		if (v1IsOutsideSoil & v2IsOutsideSoil) {
			euclidDist = 0;
		}
		else {
			
			double dx = v1c.x - v2c.x;
			double dy = v1c.y - v2c.y;
			double dz = v1c.z - v2c.z;
			
			double lengthReduction = 0;
			
			if (v1IsOutsideSoil | v2IsOutsideSoil) {				
							
				double dxy = Math.sqrt(dx * dx + dy * dy);
				double alpha = Math.atan2(dxy, dz);
				
				if (v1IsOutsideSoil) {
					
					if (v1c.z < numberOfAddedSlices + 1) lengthReduction = (numberOfAddedSlices - v1c.z) / Math.cos(alpha);
					else lengthReduction = (numberOfAddedSlices - (numberOfLastSlice - v1c.z)) / Math.cos(alpha);
					
				}
				if (v2IsOutsideSoil) {
					
					if (v2c.z < numberOfAddedSlices + 1) lengthReduction = (numberOfAddedSlices - v2c.z) / Math.cos(alpha);
					else lengthReduction = (numberOfAddedSlices - (numberOfLastSlice - v2c.z)) / Math.cos(alpha);
					
				}
			}
			
			euclidDist = calcDistanceBetweenTwoPoints(dx, dy, dz);  
			euclidDist -= Math.abs(lengthReduction);
		}
		
		//kick out negative values..
		if (euclidDist < 0) euclidDist = 0.1; 
		
		return euclidDist;
				
	}
	
	public double calcDistanceBetweenTwoPoints(double dx, double dy, double dz) {
		
		double distance = 0;
		
		double dxy = Math.sqrt(dx * dx + dy * dy);
		
		distance = Math.sqrt(dxy * dxy + dz * dz);
		
		return distance;
		
	}
	
	public double calculateAngleToZAxis(Vertex v1, Vertex v2) {
		
		Point v1c = meanVertexLocation(v1);
		Point v2c = meanVertexLocation(v2);		
		
		double dx = v1c.x - v2c.x;
		double dy = v1c.y - v2c.y;
		double dz = v1c.z - v2c.z;
		
		double angle = 0;
		
		double dxy = Math.sqrt(dx * dx + dy * dy);
		
		angle = Math.atan2(dxy, dz);
				
		if (angle < 0) angle = -angle;
		if (angle > Math.PI / 2) angle = Math.PI - angle;
		
		angle = angle / (Math.PI / 2) * 90;
		
		return angle;
		
	}
	
	public Point meanVertexLocation(Vertex v) {
				
		ArrayList<Point> vv =  v.getPoints();
		Point center = null;
		
		if (vv.size() == 1) {			
			center = new Point(vv.get(0).x, vv.get(0).y, vv.get(0).z);			
		}
		else {
			double[] x = new double[vv.size()];
			double[] y = new double[vv.size()];
			double[] z = new double[vv.size()];
			
			for (int i = 0 ; i < x.length ; i++) {				
				x[i] = vv.get(i).x;
				y[i] = vv.get(i).y;
				z[i] = vv.get(i).z;				
			}
			
			center = new Point((int)Math.round(StatUtils.mean(x)), (int)Math.round(StatUtils.mean(y)), (int)Math.round(StatUtils.mean(z)));			
		}
		
		return center;		
	}
	
	public int vertexIsConnected2TopAndBottom(Vertex nowVert, int numberOfLastSlice, int numberOfAddedSlices) {
		
		int connection = 0;
			
		ArrayList<Point> nowPoints = nowVert.getPoints();
			
		for (int k = 0 ; k < nowPoints.size() ; k++) {
	
			if (nowPoints.get(k).z < numberOfAddedSlices + 1) {
				connection = 1;return connection;
			}
			
			if (nowPoints.get(k).z >= numberOfLastSlice - numberOfAddedSlices - 1) {
				connection = 2;return connection;
			}
		}		
		
		return connection;
	}
	
	public FloydWarshallReturn runFloydWarshallAlgorithm(PlusGraph graph, int numberOfLastSlice, int numberOfSlicesAdded) {
		
		//This function had originally been copied from Ignacio Arganda-Carreras "AnalyzeSkeleton_" algorithm for ImageJ, 
		//but it has been adapted to calculate paths of least hydraulic resistance traversing a pore network by John Koestel. 
		/**
		 * The original Disclaimer in AnalyzeSkeleton_.java:
		 * 
		 * Determine the longest shortest path using the APSP (all pairs shortest path) 
		 * warshall algorithm
		 * 
		 * @param graph the graph of a tree
		 * @param shortestPathPoints list to store the longest shortest path points
		 * @return longest shortest path length
		 * @author Huub Hovens
		*/ 
		
		//output structure
		FloydWarshallReturn mFWR = new FloydWarshallReturn();
		
	
		// local fields
		/** vertex 1 of an edge */
		Vertex v1 = null;
		/** vertex 2 of an edge */
		Vertex v2 = null;
		/** the equivalent of row in a matrix */
		int row = 0;
		/** the equivalent of column in a matrix */
		int column = 0;
		/** the value of the longest shortest path */
		
		ArrayList<PlusEdge> edgeList = graph.plusEdge;
		ArrayList<Vertex> vertexList = graph.getVertices();

		//create empty adjacency and predecessor matrix
		/** the matrix that contains the length of the shortest path from vertex a to vertex b */
		double[][] adjacencyMatrix = new double[vertexList.size()][vertexList.size()];
		double[][] distanceMatrix = new double[vertexList.size()][vertexList.size()];
		/** the matrix that contains the predecessor vertex of vertex b in the shortest path from vertex a to b */
		int[][] predecessorMatrix = new int[vertexList.size()][vertexList.size()];

		// applying initial conditions
		for(int i = 0 ; i < vertexList.size(); i++)
		{
			for(int j = 0 ; j < vertexList.size(); j++)
			{
				adjacencyMatrix[i][j]= Double.POSITIVE_INFINITY;
				distanceMatrix[i][j]= Double.POSITIVE_INFINITY;
				predecessorMatrix[i][j]= -1;
			}
		}

		//also remember if a Vertex touches the top..
		ArrayList<Integer> touchesTop = new ArrayList<Integer>();
		ArrayList<Integer> touchesBottom = new ArrayList<Integer>();
		ArrayList<Vertex> topVertices = new ArrayList<Vertex>();
		ArrayList<Vertex> bottomVertices = new ArrayList<Vertex>();
		
		int cc=0;
		for (PlusEdge edge : edgeList) {
			
			cc++;
			
			v1 = edge.getV1();
			v2 = edge.getV2();
			
			if (v1 == null | v2 == null) IJ.error("Damn!!!");
			
			// use the index of the vertices as the index in the matrix			
			row = vertexList.indexOf(v1);		
			
			if (vertexIsConnected2TopAndBottom(v1, numberOfLastSlice, numberOfSlicesAdded) == 1) {			
				if(!touchesTop.contains(row)) {
					touchesTop.add(row);
					topVertices.add(v1);
				}
			}			
			if (vertexIsConnected2TopAndBottom(v1, numberOfLastSlice, numberOfSlicesAdded) == 2) {
				if (!touchesBottom.contains(row)) {
					touchesBottom.add(row);
					bottomVertices.add(v1);
				}
			}
			
			column = vertexList.indexOf(v2);			
			if (vertexIsConnected2TopAndBottom(v2, numberOfLastSlice, numberOfSlicesAdded) == 1) {
				if (!touchesTop.contains(column)) {
					touchesTop.add(column);
					topVertices.add(v2);
				}
			}
			if (vertexIsConnected2TopAndBottom(v2, numberOfLastSlice, numberOfSlicesAdded) == 2) {
				if (!touchesBottom.contains(column)) {
					touchesBottom.add(column);
					bottomVertices.add(v2);
				}
			}
			
			//if(column == -1)
			//{
			//	IJ.log("Vertex " + v2.getPoints().get(0) + " not found in the list of vertices!");
			//	continue;
			//}*/

			/* 
			 * the diagonal is 0. 
			 * 
			 * Because not every vertex is a 'v1 vertex' 
			 * the [column][column] statement is needed as well.
			 *
			 * in an undirected graph the adjacencyMatrix is symmetric 
			 * thus A = Transpose(A)
			 */
			adjacencyMatrix[row][row] = 0;
			adjacencyMatrix[column][column] = 0;
			adjacencyMatrix[row][column] = edge.constrictionFactor;  //replace length by a resistance weighted length
			adjacencyMatrix[column][row] = edge.constrictionFactor;  //replace length by a resistance weighted length
			
			//also calculate the distancee of the paths..
			distanceMatrix[row][row] = 0;
			distanceMatrix[column][column] = 0;
			distanceMatrix[row][column] = edge.getLength() - edge.lengthReduction;  //get length but correct for the added slices  
			distanceMatrix[column][row] = edge.getLength() - edge.lengthReduction;  //get length but correct for the added slices

			/* 
			 * the diagonal remains -1.
			 * for the rest I use the index of the vertex so I can later refer to the vertexList
			 * for the correct information
			 * 
			 * Determining what belongs where requires careful consideration of the definition
			 * 
			 * the array contains the predecessor of "column" in a path from "row" to "column"
			 * therefore in the other statement it is the other way around.   
			 */
			predecessorMatrix[row][row] = -1;
			predecessorMatrix[column][column] = -1;
			predecessorMatrix[row][column] = row;
			predecessorMatrix[column][row] = column;
		}
		// matrices now have their initial conditions

		mFWR.topVertices = topVertices;
		mFWR.bottomVertices = bottomVertices;
		mFWR.touchesTop = touchesTop;
		mFWR.touchesBottom = touchesBottom;
		
		// the warshall algorithm with k as candidate vertex and i and j walk through the adjacencyMatrix
		// the predecessor matrix is updated at the same time. 

		//setup parallel Floyd-Warshall algorithm..
		int cores = Runtime.getRuntime().availableProcessors();
		int numberOfThreads = cores - 1;
		ExecutorService exec = Executors.newFixedThreadPool(numberOfThreads);
		ParallelFloydWarshall pFW = new ParallelFloydWarshall(vertexList.size(), adjacencyMatrix, distanceMatrix, predecessorMatrix,
				exec, numberOfThreads);
				
		//run it
		IJ.showStatus("Running multi-threaded Floyd-Warshall algorithm ...");
		pFW.solve();
		
		//get results
		adjacencyMatrix = pFW.getAdjacencyMatrix();
		distanceMatrix = pFW.getDistanceMatrix();
		predecessorMatrix = pFW.getPredecessorMatrix();
		
		//old, single threaded version..
		/*for (int k = 0 ; k < vertexList.size(); k++)	{						
			IJ.showStatus("Running Floyd-Warshall algorithm for vertex " + (k + 1) + " of " + vertexList.size() + "...");			
			for (int i = 0 ; i < vertexList.size(); i++) {
				for (int j = 0 ; j < vertexList.size(); j++) {
					if (adjacencyMatrix[i][k] + adjacencyMatrix[k][j] < adjacencyMatrix[i][j]) {
						adjacencyMatrix[i][j] = adjacencyMatrix[i][k] + adjacencyMatrix[k][j];
						distanceMatrix[i][j] = distanceMatrix[i][k] + distanceMatrix[k][j];
						predecessorMatrix[i][j] = predecessorMatrix[k][j];

					}
				}
			}
		}*/
		
		mFWR.adjacencyMatrix = adjacencyMatrix;
		mFWR.distanceMatrix = distanceMatrix;
		mFWR.predecessorMatrix = predecessorMatrix;

		// find the shorted connection from top to bottom
		double[][] top2BottomConnection = new double[touchesTop.size()][touchesBottom.size()];
		double[][] top2BottomDistance = new double[touchesTop.size()][touchesBottom.size()];
		
		int ro = 0;
		int co = 0;
		double shortestPath = Double.POSITIVE_INFINITY;
		for(int edget : touchesTop) {
			for(int edgeb : touchesBottom) {			
				// sometimes infinities still remain
				top2BottomConnection[ro][co] = adjacencyMatrix[edget][edgeb];
				top2BottomDistance[ro][co] = distanceMatrix[edget][edgeb];
				
				//test if this path was the shortest one...
				if (adjacencyMatrix[edget][edgeb] < shortestPath) {
					shortestPath = adjacencyMatrix[edget][edgeb];					
				}
				
				co++;				
			}
			ro++;
			co=0;
		}

		//trace back the longest shortest path
		ReconstructedPath mRP = new ReconstructedPath();
		cc=0;
		for (int i = 0 ; i < top2BottomConnection.length ; i++) {
			for (int j = 0 ; j < top2BottomConnection[i].length ; j++) {
				if (top2BottomConnection[i][j] == shortestPath) {
					mRP = reconstructPathFromFloydWarshallResults(predecessorMatrix, touchesTop.get(i), touchesBottom.get(j), edgeList, vertexList, numberOfLastSlice, numberOfSlicesAdded);
				}
			}
		}
		
		// and deliver it..
		mFWR.top2BottomConnection = top2BottomConnection;
		mFWR.top2BottomEuklideanConnection = top2BottomDistance;
		
		mFWR.edgeList = edgeList;
		mFWR.vertexList = vertexList;
				
		mFWR.shortestPath = shortestPath;
		mFWR.shortest = mRP;
		
		return mFWR;

	}
	
	public ReconstructedPath reconstructPathFromFloydWarshallResults(int[][] predecessorMatrix, int startIndex, int endIndex, ArrayList<PlusEdge> edgeList, ArrayList<Vertex> vertexList, int numberOfLastSlice, int numberOfSlicesAdded) {
		
		/**
		 * Reconstruction and visualisation of the longest shortest path found by the APSP warshall algorithm
		 *  
		 * @param predecessorMatrix the Matrix which contains the predecessor of vertex b in the shortest path from a to b
		 * @param startIndex the index of the row which contains the longest shortest path
		 * @param endIndex the index of the column which contains the longest shortest path
		 * @param edgeList the list of edges 
		 * @param vertexList the list of vertices
		 * @param shortestPathPoints contains points of the longest shortest path for each graph
		 * @author Huub Hovens
		 */
		
		// We know the first and last vertex of the longest shortest path, namely a and b
		// using the predecessor matrix we can now determine the path that is taken from a to b
		// remember a and b are indices and not the actual vertices.

		ReconstructedPath mRP = new ReconstructedPath();
		ArrayList<Point> pathPoints = new ArrayList<Point>();
		ArrayList<Double> pathDistances = new ArrayList<Double>();
		
		double edgeLength = 0;
		double edgeEuclideanDistance = 0;
		int numberOfVertices = 1;
		double sumOfCN = 1;
		double bottleNeck = Double.MAX_VALUE;
		
		int b = endIndex;
		final int a = startIndex;

		while (b != a) {
			Vertex predecessor = vertexList.get(predecessorMatrix[a][b]);
			Vertex endvertex = vertexList.get(b);
			ArrayList<PlusEdge> sp_edgeslist = new ArrayList<PlusEdge>();
			Double lengthtest = Double.POSITIVE_INFINITY;
			PlusEdge shortestedge = null;

			// search all edges for a combination of the two vertices
			for (PlusEdge edge : edgeList) {

				if ((edge.getV1()==predecessor && edge.getV2()==endvertex) || (edge.getV1()==endvertex && edge.getV2()==predecessor))
				{
					// sometimes there are multiple edges between two vertices so add them to a list
					// for a second test
					sp_edgeslist.add(edge);
				}

			}
			
			// the second test
			// this test looks which edge has the shortest length in sp_edgeslist
			for (PlusEdge edge : sp_edgeslist) {
				if (edge.getLength() < lengthtest)
				{
					shortestedge = edge;
					lengthtest = edge.getLength();
				}
			}
			
			//get edge length and calculate distance between v1 and v2
			edgeLength += shortestedge.getLength() - shortestedge.lengthReduction;
			edgeEuclideanDistance += getEuclideanDistance(shortestedge.getV1(), shortestedge.getV2(), numberOfLastSlice, numberOfSlicesAdded);
			numberOfVertices++;
			
			// add vertex 1 points
			Vertex v1 = shortestedge.getV2() != predecessor ?
						shortestedge.getV2() : shortestedge.getV1();
			double v1Dist = shortestedge.getV2() != predecessor ?
							shortestedge.v2Dist : shortestedge.v1Dist;
						
			for (Point p : v1.getPoints()) {
				if( !pathPoints.contains( p ))	{
					
					//also check also if voxel is within proper domain 
					if (p.z >= numberOfSlicesAdded & p.z <= numberOfLastSlice - numberOfSlicesAdded) {
						pathPoints.add(p);
						pathDistances.add(v1Dist);
						sumOfCN += (double)v1.getBranches().size() / (double)v1.getPoints().size();;
					}

				}
			}

			// add slab points of the shortest edge to the list of points
			ArrayList<Point> slabs = shortestedge.getSlabs();
			ArrayList<Double> distance = shortestedge.distance;
			// reverse order if needed
			if (shortestedge.getV2() != predecessor)
				Collections.reverse(slabs);
			if (shortestedge.getV2() != predecessor)
				Collections.reverse(distance);
			int cc = 0;
			for (Point p : slabs) {
				//also check also if voxel is within proper domain 
				if (p.z >= numberOfSlicesAdded & p.z <= numberOfLastSlice - numberOfSlicesAdded) {
					pathPoints.add(p);
					pathDistances.add(distance.get(cc));
				}
				cc++;
			}

			// add vertex 2 points too
			Vertex v2 = shortestedge.getV2() != predecessor ?
					shortestedge.getV1() : shortestedge.getV2();
			double v2Dist = shortestedge.getV2() != predecessor ?
							shortestedge.v1Dist : shortestedge.v2Dist;
			for (Point p : v2.getPoints()) {
				if( !pathPoints.contains( p )) {
					//also check also if voxel is within proper domain 
					if (p.z >= numberOfSlicesAdded & p.z <= numberOfLastSlice - numberOfSlicesAdded) {						
						pathPoints.add(p);
						pathDistances.add(v2Dist);
						sumOfCN += (double)v2.getBranches().size() / (double)v2.getPoints().size();
					}
				}
			}

			// now make the index of the endvertex the index of the predecessor so that the path now goes from
			// a to predecessor and repeat cycle
			b = predecessorMatrix[a][b];
		}
/*		if (pathPoints.size() != 0)
		{
		//	this.spx = shortestPathPoints.get(0).x;
		//	this.spy = shortestPathPoints.get(0).y;
		//	this.spz = shortestPathPoints.get(0).z;
		}*/

		mRP.pathPoints = pathPoints;
		mRP.pathDistances = pathDistances;
		
		mRP.edgeLength = edgeLength;
		mRP.edgeEuclideanDistance = edgeEuclideanDistance;
		mRP.edgeTortuosity = edgeLength / edgeEuclideanDistance;
		mRP.globalTortuosity = edgeLength / (numberOfLastSlice - 2 * numberOfSlicesAdded - 1);
		mRP.numberOfVertices = numberOfVertices;
		mRP.sumOfCN = sumOfCN;
		
		//get bottleneck
		for (int i = 2 ; i < pathDistances.size() - 2 ; i++) {
			
			double nowDistance = pathDistances.get(i);
			if (nowDistance < bottleNeck) bottleNeck = nowDistance;
			
		}
		
		mRP.bottleNeck = bottleNeck;
		
		return mRP;		
	}
	
	public PlusGraph[] getConstrictionFactors(SkeletonResult mSP, ImagePlus distTiff, int numberOfLayersAdded) {
		
		Graph[] myGraphs = mSP.getGraph();
		RollerCaster rC = new RollerCaster();
		PlusGraph[] newGraphs = new PlusGraph[myGraphs.length];
					
		int[] numberOfEdges = new int[myGraphs.length];
		for (int i = 0 ; i < myGraphs.length ; i++) {			
			ArrayList<Edge> myEdges = myGraphs[i].getEdges();
			numberOfEdges[i] = myEdges.size();
		}
		
		//cast graphs to PlusGraphs..	
		for (int i = 0 ; i < myGraphs.length ; i++) {			
			newGraphs[i] = rC.cast2PlusGraph(myGraphs[i]);
		}
		
		//init output graph
		PlusGraph[] outGraphs = new PlusGraph[newGraphs.length];		
		
		//to add in case that vertex or slabs are outside proper ROI
		double bigDistance = 1000000;	
		double bigArea = bigDistance * bigDistance * Math.PI;		
								
		//get largest graph and calculate coordination number
		for (int i = 0 ; i < newGraphs.length ; i++) {
			
			IJ.showStatus("Calculating constriction factors for network " + (i + 1) + " of " + newGraphs.length);
	
			PlusGraph nowGraph = newGraphs[i];		
			ArrayList<PlusEdge> nowEdges =  nowGraph.plusEdge;			
			ArrayList<PlusEdge> outEdges =  new ArrayList<PlusEdge>();
						
			for (int j = 0 ; j < nowEdges.size() ; j++) {
			
				double lengthReduction = 0;				
				PlusEdge nowEdge = (PlusEdge)nowEdges.get(j);
				
				//get vertex 1
				Vertex nowVertex = nowEdge.getV1();
				ArrayList<Point> nowPoints = nowVertex.getPoints();
				double[] allV1 = new double[nowPoints.size()];
				for (int k = 0 ; k < nowPoints.size() ; k++) {				
					
					Point nowPoint = nowPoints.get(k);					
					if (nowPoints.size() == 1 & (nowPoint.z < numberOfLayersAdded) | (nowPoint.z >= distTiff.getNSlices() - numberOfLayersAdded)) {
						allV1[k] = bigDistance;
						lengthReduction++;
					}
					else {					
						distTiff.setPosition(nowPoint.z + 1);
						ImageProcessor distIP = distTiff.getProcessor();
					
						allV1[k] = distIP.getPixelValue(nowPoint.x, nowPoint.y);
					}
				}				
				nowEdge.v1Dist = StatUtils.mean(allV1);								
				 
				//get vertex 2
				nowVertex = nowEdge.getV2();
				nowPoints = nowVertex.getPoints();
				double[] allV2 = new double[nowPoints.size()];
				for (int k = 0 ; k < nowPoints.size() ; k++) {					
					
					Point nowPoint = nowPoints.get(k);					
					if (nowPoints.size() == 1 & (nowPoint.z < numberOfLayersAdded) | (nowPoint.z >= distTiff.getNSlices() - numberOfLayersAdded)) {
						allV2[k] = bigDistance;
						lengthReduction++;
					}
					else {		
						distTiff.setPosition(nowPoint.z + 1);					
						ImageProcessor distIP = distTiff.getProcessor();
					
						allV2[k] = distIP.getPixelValue(nowPoint.x, nowPoint.y);
					}
				}				
				nowEdge.v2Dist = StatUtils.mean(allV2);
				
				//get distances for slab voxels
				ArrayList<Point> slabVoxels = nowEdge.getSlabs();
				double[] area = new double[slabVoxels.size()];
				
				//check whether this slab voxel is in the added stack and remove it if true
				ArrayList<Integer> toRemove = new ArrayList<Integer>();
				for (int k = 0 ; k < slabVoxels.size() ; k++) {
					if (slabVoxels.get(k).z < numberOfLayersAdded) toRemove.add(k);
					if (slabVoxels.get(k).z >= distTiff.getNSlices() - numberOfLayersAdded) toRemove.add(k); 
				}

				for (int k = 0 ; k < slabVoxels.size() ; k++) {
			
					if (toRemove.contains(k)) {
						nowEdge.distance.add(bigDistance);	
						area[k] = bigArea;
						lengthReduction++;
					}
					else {			
						Point nowPoint = slabVoxels.get(k);
					
						distTiff.setPosition(nowPoint.z + 1);
						ImageProcessor distIP = distTiff.getProcessor();
					
						double dist = (double)distIP.getPixelValue(nowPoint.x, nowPoint.y);
						nowEdge.distance.add(dist);			
					
						area[k] = dist * dist * Math.PI;						
					}
				}
				
				//calculate constriction factor
				double constrictionFactor = 0;
				constrictionFactor += 1 / (nowEdge.v1Dist * nowEdge.v1Dist * Math.PI);
				constrictionFactor += 1 / (nowEdge.v2Dist * nowEdge.v2Dist * Math.PI);
				for (int k = 0  ; k < slabVoxels.size() ; k++) constrictionFactor += 1 / (area[k]);
				
				nowEdge.constrictionFactor = constrictionFactor;
				
				//harmonize length Reduction with real length
				nowEdge.lengthReduction = lengthReduction;
				if (lengthReduction >= nowEdge.getLength() & nowEdge.getLength() > 0.1) nowEdge.lengthReduction = nowEdge.getLength() - 0.1;
				
				//get bottleNeck
				double bottleNeck = Double.POSITIVE_INFINITY;
				if (nowEdge.v1Dist < bottleNeck) bottleNeck = nowEdge.v1Dist;
				if (nowEdge.v2Dist < bottleNeck) bottleNeck = nowEdge.v2Dist;
				for (int k = 0  ; k < slabVoxels.size() ; k++) if (nowEdge.distance.get(k) < bottleNeck) bottleNeck = nowEdge.distance.get(k); 
				
				nowEdge.bottleNeck = bottleNeck;
				
				//add results to the outGraph
				outEdges.add(nowEdge);
				
			}
			
			nowGraph.plusEdge = outEdges;
			outGraphs[i] = nowGraph;
		}

		return outGraphs;
		
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	public void tailoredPoreSpaceAnalyzer(int imageNumber, InputOutput.MyFileCollection mFC, RoiHandler.ColumnRoi colRoi, MenuWaiter.PoreSpaceAnalyzerOptions mPSA) {
				
		//try to free up some memory			
		IJ.freeMemory();IJ.freeMemory();
		
		//check whether this is a windows OS
		String pathSep = "/";
		String checkname = mFC.myPreOutFolder;
		String testString = checkname.substring(2, 3);
		if (testString.equalsIgnoreCase("\\")) pathSep = "\\";
		
		InputOutput jIO = new InputOutput();
		ROIMorphoProps myP = new ROIMorphoProps();
		PoreClusterProps mPCP = new PoreClusterProps();		
		ImageManipulator jIM = new ImageManipulator();
		RollerCaster rC = new RollerCaster();
		
		//colRoi.nowTiff.updateAndDraw();
		//colRoi.nowTiff.show();
		
		//remove holes from pore system;
		if (mPSA.removeHoles) {
			colRoi.nowTiff = jIM.removeHoles(colRoi.nowTiff);
			colRoi.nowTiff.updateAndDraw();
			//colRoi.nowTiff.show();
		}
		
		//init some very basic variables..		
		String nowImageName = mFC.colName;
		
		//create output folders and save statistics		
		IJ.showStatus("Creating directories for the pore-cluster properties ...");
		String myOutFolder = "Stats";
		String myOutPath = mFC.myPreOutFolder + pathSep + myOutFolder;
		new File(myOutPath).mkdir();
		if (mPSA.performParticleAnalyses == true) new File(myOutPath + pathSep + "Clusters").mkdir();
		if (mPSA.calcAnisotropy == true) new File(myOutPath + pathSep + "Anisotropy").mkdir();
		if (mPSA.calcThickness == true) new File(myOutPath + pathSep + "Thickness").mkdir();		
		String outROIPath = mFC.myPreOutFolder + pathSep + "Stats" + pathSep + nowImageName + ".roi";
		String outClustPath = mFC.myPreOutFolder + pathSep + "Stats" + pathSep + "Clusters" + pathSep + nowImageName + ".clust";		
		String outAnisotPath = mFC.myPreOutFolder + pathSep + "Stats" + pathSep + "Anisotropy" + pathSep + nowImageName + ".aniso";
		String outThickPath = mFC.myPreOutFolder + pathSep + "Stats" + pathSep + "Thickness" + pathSep + nowImageName + ".psd";
		
		/////////////////////////////////////////////////////////////////
		//calculate fractal dimension
		/////////////////////////////////////////////////////////////////
		
		IJ.showStatus("Calculating fractal properties ...");
		
		if (mPSA.calcFractal == true) {
			
			FractalProperties myFP = calculateFractalProperties(colRoi, mPSA.mRSO);		

			myP.surfaceFractalDimension = myFP.surfaceFractalDim;
			
		}
		
		IJ.freeMemory();IJ.freeMemory();
		
		//colRoi.nowTiff.show();
				
		/////////////////////////////////////////////////////////////////
		//calculate anisotropy
		/////////////////////////////////////////////////////////////////
		
		IJ.showStatus("Investigating anisotropies ...");
		
		if (mPSA.calcAnisotropy == true) {

			AnisotropyResults arreArre = new AnisotropyResults(); 
			
			boolean roiIsCylindrical = false;
			if (mPSA.mRSO.choiceOfRoi.equalsIgnoreCase("Cylinder")) roiIsCylindrical = true;
			if (mPSA.mRSO.choiceOfRoi.equalsIgnoreCase("RealSample")) roiIsCylindrical = true;
			
			arreArre = calculateAnisotropy(roiIsCylindrical, colRoi.nowTiff, mPSA.mRSO.areaOfInterest);
			myP.anisotropy = arreArre.anisotropy;
			myP.mainAlignment = arreArre.mainAlignment;			
			
			//save anisotropy statistics in file
			jIO.writeROIAnisoResults(nowImageName, outAnisotPath, arreArre);

		}
		
		//colRoi.nowTiff.show();
				
		IJ.freeMemory();IJ.freeMemory();
		
		/////////////////////////////////////////////////////////////////		
		//calculate macroporosity
		/////////////////////////////////////////////////////////////////
		
		IJ.showStatus("Calculating macroporosities ...");
		
		if (myP.phaseVolume == 0) myP.phaseVolume = macroPoreVolume(colRoi.nowTiff);		
		
		double[] bulkSoilVolume = {0, 0, 0, 0};
		ImagePlus surfTiff = new ImagePlus();
		if (mPSA.mRSO.choiceOfRoi.equalsIgnoreCase("RealSample") & mPSA.mRSO.includeSurfaceTopography) {
			String[] myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC);
			surfTiff = jIO.openTiff3D(mFC.mySurfaceFolder + pathSep + myGandS[1]);
			if (mFC.myCutSurfaceFolder != null) {
				surfTiff = jIO.openTiff3D(mFC.myCutSurfaceFolder + pathSep + myGandS[1]);
			}
			bulkSoilVolume = calculateTheBulkSoilVolume(myGandS[0], surfTiff, mPSA. mRSO, mFC);
			myP.roiBulkVolume = bulkSoilVolume[0];			
		}
		else { 
			bulkSoilVolume[0] = Math.round(mPSA.mRSO.areaOfInterest) * (double)colRoi.nowTiff.getNSlices();
			myP.roiBulkVolume = bulkSoilVolume[0];
			surfTiff = null;
		}
		
		myP.phaseVolumeFraction = myP.phaseVolume / myP.roiBulkVolume;	
				
		IJ.freeMemory();IJ.freeMemory();

		/////////////////////////////////////////////////////////////////
		// do particle analyses
		/////////////////////////////////////////////////////////////////
		
		//check if the clusters need to be identified
		mPCP.containsAParticleAnalysis = false;
		if (mPSA.performParticleAnalyses == true) {
		
			mPCP.containsAParticleAnalysis = true;
			
			////////////////////////////////////////////////////////
			// get clusters	
			////////////////////////////////////////////////////////
			
			IJ.showStatus("Identifying connected pore-clusters ...");

			FloodFillComponentsLabeling3D myFCCL = new FloodFillComponentsLabeling3D(26, 32);
			ImageStack myLabelStack = myFCCL.computeLabels(colRoi.nowTiff.getStack());
											
			IJ.freeMemory();IJ.freeMemory();
			
			//calculate surface and Euler number.. using MorphoLibJ
			ImagePlus labelTiff = new ImagePlus();
			labelTiff.setStack(myLabelStack);				
			
			int[] labels = LabelImages.findAllLabels(myLabelStack);
			int numOfObjects = labels.length;
			Calibration calib = labelTiff.getCalibration();
			
			//labelTiff.updateAndDraw();
			//labelTiff.show();
			
			//surface and curvature
			IntrinsicVolumesAnalyzer3D mIRA = new IntrinsicVolumesAnalyzer3D();
			mIRA.setConnectivity(26);
			Result[] myMorphos = mIRA.analyzeRegions(myLabelStack, labels, calib);
			
			double[] morphoVols = new double[myMorphos.length];
			double[] morphoSurf = new double[myMorphos.length];
			double[] morphoCurv = new double[myMorphos.length];
			double[] morphoEuler = new double[myMorphos.length];
			for (int i = 0 ; i < myMorphos.length ; i++) {
				morphoVols[i] = myMorphos[i].volume;
				morphoSurf[i] = myMorphos[i].surfaceArea;
				morphoCurv[i] = 2 * Math.PI * myMorphos[i].meanBreadth;
				morphoEuler[i] = myMorphos[i].eulerNumber;
			}
						
			//myP.phaseVolume = StatUtils.sum(morphoVols);
			myP.surfaceArea = StatUtils.sum(morphoSurf);
			myP.meanCurvature = StatUtils.sum(morphoCurv);			
			myP.eulerNumber = StatUtils.sum(morphoEuler);	
			
			//calculate percolating clusters ... needs the other parameters.. and is needed for following parameters	
			ArrayList<Integer> conTop = new ArrayList<Integer>();
			ArrayList<Integer> conBot = new ArrayList<Integer>();
			ArrayList<Integer> isPercolating = new ArrayList<Integer>();
			if (numOfObjects > 0) {
				conTop = check4TouchingTheTop(labelTiff, colRoi.nowTiff, surfTiff);
				conBot = check4TouchingTheBottom(labelTiff, colRoi.nowTiff, surfTiff);
				isPercolating = checkForPercolatingClusters(numOfObjects, conTop, conBot);				
			}
			
			//free memory	
			IJ.freeMemory();IJ.freeMemory();
			
			//catch in case that numOfObjects < 0
			boolean hasItObjects = true;
			if (numOfObjects < 1) {
				hasItObjects = false;
				numOfObjects = 1;
			}
				
			// init variables for export of results			
			int[] id = new int[numOfObjects];	
			
			double[] volume = new double[numOfObjects];						//Vol: particle volume
			double[] surfaceArea = new double[numOfObjects];				//SA: surface area (0 if too small for mesh to be produced; see warning log)
			
			double[] xCenter = new double[numOfObjects];					//x Cent: x-coordinate of particle centroid
			double[] yCenter = new double[numOfObjects];					//y Cent: y-coordinate of particle centroid
			double[] zCenter = new double[numOfObjects];					//z Cent: z-coordinate of particle centroid
			int[] xMin = new int[numOfObjects];								//x Min: min x-coordinate of particle 
			int[] yMin = new int[numOfObjects];								//y Min: min y-coordinate of particle 
			int[] zMin = new int[numOfObjects];								//z Min: min z-coordinate of particle 
			int[] xMax = new int[numOfObjects];								//x Max: max x-coordinate of particle 
			int[] yMax = new int[numOfObjects];								//y Max: max y-coordinate of particle 
			int[] zMax = new int[numOfObjects];								//z Max: max z-coordinate of particle 
			
			boolean[] tTop = new boolean[numOfObjects];						//cluster is touching top..		
			boolean[] tBot = new boolean[numOfObjects];						//cluster is touching bot..		
			boolean[] percolating = new boolean[numOfObjects];				//cluster is percolating or not..			
		
			double[] euler = new double[numOfObjects];						//Euler (xi): Euler characteristic of the particle
			double[] meanCurvature = new double[numOfObjects];						//Holes (beta1): number of topological holes (handles) in the particle
						
			//sort results 
			int[] sortedIndices = AndAllTheRest.getIndicesInOrder(morphoVols);
			
			//write results into a biiiiig table		
			if (hasItObjects) {
				for (int i = 0; i < morphoVols.length; i++) {			
					if (morphoVols[sortedIndices[i]] > 0) {			
						
						IJ.showStatus("Sorting out properties of the individual pore-clusters ...");
						
						id[i] = sortedIndices[i] + 1;				
						volume[i] = morphoVols[sortedIndices[i]];
		
						surfaceArea[i] = morphoSurf[sortedIndices[i]];		
						meanCurvature[i] = morphoCurv[sortedIndices[i]];
						euler[i] = morphoEuler[sortedIndices[i]];
	
					}
				}
			}

			//also find percolating clusters
			for (int i = 0 ; i < id.length ; i++) {
				if (conTop.contains(id[i])) tTop[i] = true;
				else  tTop[i] = false;
				
				if (conBot.contains(id[i])) tBot[i] = true;
				else  tBot[i] = false;
				
				if (isPercolating.contains(id[i])) percolating[i] = true;
				else  percolating[i] = false;
			}
			
			//assign results to output structure		
			mPCP.hasItObjects = hasItObjects;
			mPCP.id = id;
			
			mPCP.volume = volume;
			mPCP.surfaceArea = surfaceArea;
			mPCP.meanCurvature = meanCurvature;		
			mPCP.euler = euler;
			
			mPCP.xCenter = xCenter;
			mPCP.yCenter = yCenter;		
			mPCP.zCenter = zCenter;
			mPCP.minX = xMin;
			mPCP.minY = yMin;
			mPCP.minZ = zMin;
			mPCP.maxX = xMax;
			mPCP.maxY = yMax;
			mPCP.maxZ = zMax;

			mPCP.touchesTop = tTop;
			mPCP.touchesBot = tBot;
			mPCP.isPercolating = percolating;	
							
			/*//calculate connected correlation length and chi
			if (mPSA.calcChi == true) {						
				myP = calculateChiOfPercolation(labelTiff, mFC, myP);
			}
			else {
				double[] dummy = new double[myP.id.length];
				for (int i = 0 ; i < myP.id.length ; i++) dummy[i] = -1;				
				myP.connectedCorrelationLength = dummy;
				myP.chiOfPercolation = dummy;				
			}*/			
			
			//calculate global connection probability (gamma)
			double sumOfSquaredClusterSizes = 0;
			double rF = 1000000;
			for (int i = 0 ; i < mPCP.id.length ; i++) {
				double squaredClusterSize = mPCP.volume[i] / rF * mPCP.volume[i] / rF;
				sumOfSquaredClusterSizes += squaredClusterSize;
			}
			myP.gamma = sumOfSquaredClusterSizes / (myP.phaseVolume / rF * myP.phaseVolume / rF);			
						
			//save desired images
			if (mPSA.plotLabels == true) {
				myOutFolder = "ClusterLabels";
				myOutPath = mFC.myPreOutFolder + pathSep + myOutFolder;
				new File(myOutPath).mkdir();				
				String clusterLabelImageName = nowImageName + "_ClusterLables.tif";
				jIO.tiffSaver(myOutPath, clusterLabelImageName, labelTiff);
				
				if (!mPSA.plotPercolation & !mPSA.plotPoresConnected2Top) {
					labelTiff.unlock();
					labelTiff.flush();
				}
			}
			
			//top-connected porosity calculation..
			double topConVolume = 0; 
			ArrayList<Integer> touchesTopList = new ArrayList<Integer>();
			for (int i = 0 ; i < percolating.length ; i++) if (mPCP.touchesTop[i]) {				
				touchesTopList.add(id[i]);	
				topConVolume += volume[i];
			}
			myP.volumeFractionConnected2Top = topConVolume / myP.roiBulkVolume;
			if (mPSA.plotPoresConnected2Top) {
				
				ImagePlus touchesTopImp = extractSpecificPoreLabels(labelTiff, touchesTopList);
								
				myOutFolder = "VolumeConnected2Top";
				myOutPath = mFC.myPreOutFolder + pathSep + myOutFolder;
				new File(myOutPath).mkdir();
				String nowDir = mFC.myPreOutFolder + pathSep + myOutFolder;		
				String volumeImageName = nowImageName + "_VolCon2Top.tif";	
				jIO.tiffSaver(nowDir, volumeImageName, touchesTopImp);
				
				touchesTopImp.unlock();touchesTopImp.flush();				
			}			
		
			//percolating porosity calculation..
			double percVolume = 0; 
			ArrayList<Integer> percolatingClusters = new ArrayList<Integer>();
			for (int i = 0 ; i < percolating.length ; i++) if (percolating[i]) {
				percolatingClusters.add(id[i]);
				percVolume += volume[i];
			}
			myP.percolatingVolumeFraction = percVolume / myP.roiBulkVolume;
			if (mPSA.plotPercolation & percVolume > 0) {
				
				ImagePlus percPorosityImp = extractSpecificPoreLabels(labelTiff, percolatingClusters);
				myOutFolder = "PercolatingVolume";
				myOutPath = mFC.myPreOutFolder + pathSep + myOutFolder;
				new File(myOutPath).mkdir();
				String nowDir = mFC.myPreOutFolder + pathSep + myOutFolder;		
				String volumeImageName = nowImageName + "_PercVol.tif";	
				jIO.tiffSaver(nowDir, volumeImageName, percPorosityImp);
				
				percPorosityImp.unlock();percPorosityImp.flush();				
			}
			
			labelTiff.unlock();labelTiff.flush();
			
			System.gc();System.gc();		
			IJ.freeMemory();IJ.freeMemory();
			System.gc();System.gc();
			IJ.wait(2000);
			IJ.freeMemory();
						
			////////////////////////////////////////////////////////
			//calculate thicknesses
			////////////////////////////////////////////////////////
			
			LocalThicknessWrapper myLTW = new LocalThicknessWrapper();
			myLTW.calibratePixels = false;
			myLTW.maskThicknessMap = true;
			myLTW.threshold = 128;
			myLTW.setSilence(true);
			ImagePlus thickTiff = null;
			if (mPSA.plotThickness == true | mPSA.calcThickness == true) {					
				if (mPSA.mRSO.includeSurfaceTopography) myLTW.processImage(colRoi.surfaceNotCut);
				else myLTW.processImage(colRoi.nowTiff);				
				thickTiff = myLTW.getResultImage();
				myP.averagePhaseDiameter = calculateAverageValue(thickTiff);
				
				//extract pore size distribution				
				double[] classBounds = new double[200];
				for (int i = 0 ; i < classBounds.length ; i++) classBounds[i] = 0.5 * i;
				classBounds[classBounds.length - 1] = 100000;
				
				extractPoresizeDistro(outThickPath, thickTiff, classBounds);
				
				//delete image from heap in case it is not going to be saved..
				if (mPSA.plotThickness == false) {
					thickTiff.unlock();thickTiff.flush();
				}
			}
			
			if (mPSA.plotThickness == true) {			
				myOutFolder = "Thickness";
				myOutPath = mFC.myPreOutFolder + pathSep + myOutFolder;
				new File(myOutPath).mkdir();
				String nowDir = myOutPath;		
				String thicknessImageName = nowImageName + "_Thickness.tif";
				jIO.tiffSaver(nowDir, thicknessImageName, thickTiff);			
				
				thickTiff.unlock();thickTiff.flush();
			}			
			
			System.gc();System.gc();		
			IJ.freeMemory();IJ.freeMemory();
			System.gc();System.gc();
			
			////////////////////////////////////////////////////////
			// calculate critical pore diameter
			////////////////////////////////////////////////////////			
			
			if (mPSA.calcCriticalPoreDiameter == true | mPSA.plotDistanceMap) {
				
				//check which clusters are percolating and write them into a list..
				//ArrayList<Integer> percolatingClusters = new ArrayList<Integer>();
				//for (int i = 0 ; i < isPercolating.length ; i++) if (isPercolating[i] == true) percolatingClusters.add(i + 1);
				
				//if there is not even one cluster, the critical pore diameters is <= image resolution
				myP.criticalPhaseDiameter = -1 ;
				
				ImagePlus distTiff = colRoi.nowTiff.duplicate();
				
				//create a distance map
				if ((mPSA.calcCriticalPoreDiameter & !percolatingClusters.isEmpty()) | mPSA.plotDistanceMap | mPSA.calcAverageDistance) {
					
					IJ.showStatus("Performing Euklidean distance transform ...");
					
					float[] floatWeights = ChamferWeights3D.BORGEFORS.getFloatWeights();
					boolean normalize = true;
					DistanceTransform3D dt = new DistanceTransform3DFloat(floatWeights, normalize);
					ImageStack result = dt.distanceMap(colRoi.nowTiff.getStack());
					
					distTiff = new ImagePlus("dist",result);
					
					myP.averageDistance2PhaseBoundary = calculateAverageValue(distTiff);

				}
				
				if (mPSA.plotDistanceMap) {
					myOutFolder = "DistanceMap";
					myOutPath = mFC.myPreOutFolder + pathSep + myOutFolder;
					new File(myOutPath).mkdir();
					String nowDir = mFC.myPreOutFolder + pathSep + "DistanceMap";		
					String distImageName = nowImageName + "_DistanceMap.tif";
					jIO.tiffSaver(nowDir, distImageName, distTiff);
				}
					
				//calculate critical pore diameter	
				if (mPSA.calcCriticalPoreDiameter) {
					if (!percolatingClusters.isEmpty()) {
				
						myP.phasePercolates = 1;
									
						int maxVolume = (int)Math.round(volume[0]) + 1;
						double[] results = calculateCriticalPoreDiameter(distTiff, maxVolume, surfTiff, mFC);
						myP.criticalPhaseDiameter = results[0]; 
						myP.percolationThreshold = results[1] / myP.roiBulkVolume;
						
					}
					else {
						
						myP.phasePercolates = 0;
						
						int maxVolume = (int)Math.round(volume[0]) + 1;
						
						if (myP.phaseVolumeFraction > 0) {			//catch in case there are no pores..
							myP.percolationThreshold = calculatePercolationThreshold(colRoi.nowTiff.duplicate(), maxVolume, surfTiff, mFC) / myP.roiBulkVolume;
						}
						else myP.percolationThreshold = -1;
							
					}
				}
				
				//calculate KT87 volume
				if (mPSA.calcKT87Volume) {
					
					ImagePlus kt87Volume = new ImagePlus();
					
					if (myP.criticalPhaseDiameter > 6) {
					
						double myThreshold = myP.criticalPhaseDiameter / 2 / 3; // divided by 2 because the distanceTiff shows the radii, not the diameters
						
						kt87Volume = extractKT87volume(distTiff, myThreshold, surfTiff, mFC);
						
						//calculate KT87 porosity
						//myP.KT87VolumeFraction = calcVolume(kt87Volume) / myP.roiBulkVolume;				
					
						//plot KT87volume
						if (mPSA.plotKTVolume) {
							myOutFolder = "KT87Volume";
							myOutPath = mFC.myPreOutFolder + pathSep + myOutFolder;
							new File(myOutPath).mkdir();	
							String nowDir = mFC.myPreOutFolder + pathSep + "KT87Volume";	
							String KT87ImageName = nowImageName + "_KT87Volume.tif";
							jIO.tiffSaver(nowDir, KT87ImageName, kt87Volume);									
						}						
					}
					else {						
						
						//myP.KT87VolumeFraction = myP.percolatingVolumeFraction;
						
					}			
					
					kt87Volume.unlock();kt87Volume.flush();
					
				}
				
				distTiff.unlock();distTiff.flush();
				
			}
		}		//closes if-tag checking whether the particle analyses tag
		
		System.gc();System.gc();
		IJ.showStatus("Trying to clean up memory ... ");
		IJ.freeMemory();IJ.freeMemory();
		System.gc();System.gc();
		
		//do thickness calculation without particle analyzer
		if (!mPSA.performParticleAnalyses & mPSA.plotThickness) {
			LocalThicknessWrapper myLTW = new LocalThicknessWrapper();
			myLTW.calibratePixels = false;
			myLTW.maskThicknessMap = true;
			myLTW.threshold = 128;
			myLTW.setSilence(true);
			ImagePlus thickTiff = null;								
			if (mPSA.mRSO.includeSurfaceTopography) myLTW.processImage(colRoi.surfaceNotCut);
			else myLTW.processImage(colRoi.nowTiff);				
			thickTiff = myLTW.getResultImage();						
			IJ.freeMemory();IJ.freeMemory();
			myP.averagePhaseDiameter = calculateAverageValue(thickTiff);
			
			myOutFolder = "PoreThick";
			myOutPath = mFC.myPreOutFolder + pathSep + myOutFolder;
			new File(myOutPath).mkdir();
			String nowDir = mFC.myPreOutFolder + pathSep + "PoreThick";		
			String thicknessImageName = nowImageName + "_Thickness.tif";
			jIO.tiffSaver(nowDir, thicknessImageName, thickTiff);
		}
		
		System.gc();System.gc();
		IJ.showStatus("Trying to clean up memory ... ");
		IJ.freeMemory();IJ.freeMemory();
		System.gc();System.gc();
		
		//write Results
		IJ.showStatus("Writing results ... ");
		jIO.writeROIMorphoResults(nowImageName, outROIPath, myP);
		
		if (mPSA.performParticleAnalyses == true) jIO.writeClusterMorphoResults(nowImageName, outClustPath, mPCP);
		
	}
	
	public double calcVolume(ImagePlus nowTiff) {
		
		double myVolume = 0;
		
		for (int z = 1 ; z <= nowTiff.getNSlices() ; z++) {	
			
			nowTiff.setPosition(z);
			ImageProcessor nowIP = nowTiff.getProcessor();			
			
			//loop through all pixel within cross-section
			for (int x = 0 ; x < nowTiff.getWidth() ; x++) {
				for (int y = 0 ; y < nowTiff.getHeight() ; y++) {							
					int nowPixel = nowIP.getPixel(x, y);
					if (nowPixel > 1) myVolume++;
				}
			}				
		}		
		
		return myVolume;
		
	}
	
	public ImagePlus extractKT87volume(ImagePlus distTiff, double threshold, ImagePlus surfTiff, InputOutput.MyFileCollection mFC) {
		
		Dilate_ mD = new Dilate_();
		
		//set all values smaller than the threshold to 0
		ImageStack largerPores = new ImageStack(distTiff.getWidth(), distTiff.getHeight());
		ImagePlus laPores = new ImagePlus();
		for (int z = 1 ; z <= distTiff.getNSlices() ; z++) {	
			
			distTiff.setPosition(z);
			ImageProcessor nowIP = distTiff.getProcessor();
			ImageProcessor laIP = new ByteProcessor(distTiff.getWidth(), distTiff.getHeight());
			
			//loop through all pixel within cross-section
			for (int x = 0 ; x < distTiff.getWidth() ; x++) {
				for (int y = 0 ; y < distTiff.getHeight() ; y++) {							
					double nowPixel = nowIP.getPixelValue(x, y);
					if (nowPixel >= threshold) laIP.putPixel(x, y, 255);
					else laIP.putPixel(x, y, 0);							
				}
			}					
			largerPores.addSlice(laIP);
		}		
		laPores.setStack(largerPores);
		
		//laPores.updateAndDraw();
		//laPores.show();
						
		IJ.freeMemory();IJ.freeMemory();
		
		//find pore clusters in laPores and check if they form a percolating path		
		FloodFillComponentsLabeling3D myFCCL = new FloodFillComponentsLabeling3D(26, 32);
		ImageStack myLabelStack = myFCCL.computeLabels(laPores.getStack());						
		ImagePlus laParLaTiff = new ImagePlus();
		laParLaTiff.setStack(myLabelStack);
		
		int[] labels = LabelImages.findAllLabels(myLabelStack);
		int numOfObjects = labels.length;
		
		ArrayList<Integer> laConTop = check4TouchingTheTop(laParLaTiff, laPores, surfTiff);	
		ArrayList<Integer> laConBot = check4TouchingTheBottom(laParLaTiff, laPores, surfTiff);
		ArrayList<Integer> percolatingClusters = checkForPercolatingClusters(numOfObjects, laConTop, laConBot);

		//remove all isolated clusters
		ImageStack outStack = new ImageStack(laParLaTiff.getWidth(), laParLaTiff.getHeight());
		ImagePlus outTiff = new ImagePlus();
		for(int z = 0 ; z < laParLaTiff.getNSlices() ; z++) {
			
			laParLaTiff.setSlice(z + 1);
			ImageProcessor nowIP = laParLaTiff.getProcessor();
			ImageProcessor outIP = new ByteProcessor(nowIP.getWidth(), nowIP.getHeight());
			
			for (int x = 0 ; x < nowIP.getWidth() ; x++) {
				for (int y = 0 ; y < nowIP.getHeight() ; y++) {
					
					int nowPix = (int)Math.round(nowIP.getPixelValue(x, y));
					if (percolatingClusters.contains(nowPix)) outIP.putPixel(x, y, 255);
					else outIP.putPixel(x, y, 0);
					
				}
			}
			
			outStack.addSlice(outIP);
			
		}
		
		outTiff.setStack(outStack);
		
		//dilate until resolution is reached again
		
		for (int i = 0 ; i < (int)Math.round(threshold) ; i++) outTiff = mD.dilate(outTiff, 255, false);
		 
		return outTiff;
				
	}
	
	public double calculateAverageValue(ImagePlus thickTiff) {
		
		ArrayList<Float> listOfThicknesses = new ArrayList<Float>();
		
		for (int z = 0 ; z < thickTiff.getNSlices() ; z++) {
			
			thickTiff.setSlice(z + 1);
			ImageProcessor nowIP = thickTiff.getProcessor();
			
			for (int x = 0 ; x < thickTiff.getWidth() ; x++) {
				for (int y = 0 ; y < thickTiff.getHeight() ; y++) {
					
					float nowPix = nowIP.getPixelValue(x, y);
					if (nowPix > 0) listOfThicknesses.add(nowPix);
					
				}
			}
		}
		
		double[] thickArray = new double[listOfThicknesses.size()]; 
		for (int i = 0 ; i < thickArray.length ; i++) thickArray[i] = listOfThicknesses.get(i);
		
		double avgThick = StatUtils.mean(thickArray);
		
		return avgThick;
	}
	
	public DistanceProps calculatePhaseDistanceProperties(ImagePlus thickTiff, RoiHandler.ColumnRoi pRoi) {
		
		DistanceProps mDP = new DistanceProps();
		HistogramStuff hist = new HistogramStuff();
			
		//calcluate 3 mm cutoff
		double my3mm = 3f / (pRoi.voxelSizeInMicron / 1000);
		
		//extract histogram
		int[] myHist = hist.extractHistograms16(thickTiff);
		myHist[0] = 0;
		
		//caclulate aerated voxels
		float ovox = 0;
		float avox = 0;
		for (int i = 0 ; i < my3mm ; i++) ovox += myHist[i];
		for (int i = (int)Math.ceil(my3mm) ; i < myHist.length ; i++) avox += myHist[i];
		
		//assign results
		mDP.averageDistance2PhaseBoundaryInMM = hist.findMeanFromHistogram(myHist);
		mDP.medianDistance2PhaseBoundaryInMM = hist.findMedianFromHistogram(myHist);
		mDP.fractionLessThan3MMFromPhaseBoundary = ovox / (pRoi.area * thickTiff.getNSlices());
		mDP.fractionMoreThan3MMFromPhaseBoundary = avox / (pRoi.area * thickTiff.getNSlices());;
		
		return mDP;
	}
	
	public ImagePlus findClusterConnected2Top(ImagePlus nowTiff, double drainingDepthInVx) {
		
		BoundingBox3D mBB = new BoundingBox3D();
		
		//nowTiff.updateAndDraw();nowTiff.show();
		
		//find connected pore clusters
		IJ.showStatus("Find air phases connected to top..");
		FloodFillComponentsLabeling3D myFCCL = new FloodFillComponentsLabeling3D(26, 32);
		ImageStack myLabelStack = myFCCL.computeLabels(nowTiff.getStack());
		
		//calculate surface and Euler number.. using MorphoLibJ
		ImagePlus labelTiff = new ImagePlus();
		labelTiff.setStack(myLabelStack);				
		
		//get Topmost coordinate of each cluster
		int[] labels = LabelImages.findAllLabels(myLabelStack);
		int numOfObjects = labels.length;
		Calibration calib = labelTiff.getCalibration();		
		Box3D[] myBox = mBB.analyzeRegions(myLabelStack, labels, calib);
		
		//labelTiff.updateAndDraw();labelTiff.show();
			
		// check if the cluster is draining
		ArrayList<Integer> isDraining = new ArrayList<Integer>();
		if (numOfObjects > 0 ) {	
			for (int i = 0 ; i < myBox.length ; i++) {
				if (myBox[i].getZMin() <= drainingDepthInVx) isDraining.add(i + 1);
			}
		}		
		
		//extract Image
		ImagePlus drainingClusters = extractSpecificPoreLabels(labelTiff, isDraining);
		
		//touchesTopImp.updateAndDraw();touchesTopImp.show();
							
		return drainingClusters;
				
	}
	
	public ImagePlus findClusterConnected2Bottom(ImagePlus nowTiff, double wettingHeightInVX) {
		
		BoundingBox3D mBB = new BoundingBox3D();
		
		//find connected pore clusters
		FloodFillComponentsLabeling3D myFCCL = new FloodFillComponentsLabeling3D(26, 32);
		ImageStack myLabelStack = myFCCL.computeLabels(nowTiff.getStack());
		
		//calculate surface and Euler number.. using MorphoLibJ
		ImagePlus labelTiff = new ImagePlus();
		labelTiff.setStack(myLabelStack);
		int[] labels = LabelImages.findAllLabels(myLabelStack);
		int numOfObjects = labels.length;
		Calibration calib = labelTiff.getCalibration();			
		Box3D[] myBox = mBB.analyzeRegions(myLabelStack, labels, calib);
		
		//calculate percolating clusters ... needs the other parameters.. and is needed for following parameters	
		ArrayList<Integer> isWetting = new ArrayList<Integer>();
		if (numOfObjects > 0) {	
			for (int i = 0 ; i < myBox.length ; i++) {
				if (myBox[i].getZMin() >= wettingHeightInVX) isWetting.add(i);
			}
		}		
		
		//extract Image
		ImagePlus touchesTopImp = extractSpecificPoreLabels(labelTiff, isWetting);
							
		return touchesTopImp;
				
	}

	public double tMaxLaplace(ImageStack nowStack, ImageStack lapStack, int myMode) {
		
		HistogramStuff hist = new HistogramStuff();
		ImageManipulator jIM = new ImageManipulator();
		
		ImagePlus lapTiff = new ImagePlus();
		lapTiff.setStack(lapStack);
		
		ImagePlus origTiff = new ImagePlus();
		origTiff.setStack(nowStack);
	
		//get histogram of the gradient image
		int[] myHistGradient = hist.sampleHistogram(lapTiff);				
		myHistGradient[0] = 0;
		int myCutoffGradient = hist.findTheKnee(myHistGradient);
		
		//segment gradient image..
		ImagePlus myGradientMask = lapTiff.duplicate();
		myGradientMask = jIM.binarizeGradientMask(myGradientMask, myCutoffGradient);
		//myGradientMask.show();
		
		//sample upper threshold for this segment
		double weightSum = 0;
		double weightAndSampleSum = 0;
		for (int i = 0 ; i < lapStack.getSize() ; i++) {
			
			IJ.showStatus("Calculating threshold from Laplace mask slice " + (i+1) + "/" + lapTiff.getNSlices());
			
			//set all image to the correct position
			origTiff.setPosition(i+1); //remember that this image has one slice more on top and bottom, respectively..
			lapTiff.setPosition(i+1);
			myGradientMask.setPosition(i+1);
			
			//get each images IPs
			ImageProcessor origIP = origTiff.getProcessor();
			ImageProcessor gradIP = lapTiff.getProcessor();
			ImageProcessor maskIP = myGradientMask.getProcessor();
	
			//because imageCalculator does not work... do it by hand.. again..
			for (int x = 0 ; x < origIP.getWidth(); x++) {
				for (int y = 0 ; y < origIP.getHeight(); y++) {
					if (maskIP.get(x, y) > 0 && origIP.get(x, y) > 0 && origIP.get(x, y) < myMode) {
						double nowWeight = gradIP.getPixelValue(x, y);
						double nowSample = origIP.getPixelValue(x, y);
						weightSum = weightSum + nowWeight;
						weightAndSampleSum = weightAndSampleSum + nowWeight * nowSample;
					}							
				}
			}
		}		
		
		double tMax = weightAndSampleSum / weightSum;
		
		return tMax;
	}

	public double tMaxSobel(ImageStack nowStack, ImageStack sobStack, int myMode) {
	
		//init units
		HistogramStuff hist = new HistogramStuff();
		ImageManipulator jIM = new ImageManipulator();
		
		//init variables
		ImagePlus gradTiff = new ImagePlus();
		gradTiff.setStack(sobStack);
		
		ImagePlus torigTiff = new ImagePlus();
		torigTiff.setStack(nowStack);
		
		//get histogram of the gradient image
		int[] myHistGradient = hist.sampleHistogram(gradTiff);				
		myHistGradient[0] = 0;
		int myCutoffGradient = hist.findTheKnee(myHistGradient);
		
		//segment gradient image..
		ImagePlus myGradientMask = gradTiff.duplicate();
		myGradientMask = jIM.binarizeGradientMask(myGradientMask, myCutoffGradient);
		//myGradientMask.show();
		
		//sample upper threshold for this segment
		double weightSum = 0;
		double weightAndSampleSum = 0;
		for (int i = 0 ; i < sobStack.getSize() ; i++) {
			
			IJ.showStatus("Calculating threshold from Sobel mask slice " + (i+1) + "/" + gradTiff.getNSlices());
			
			//set all image to the correct position
			torigTiff.setPosition(i+1); //remember that this image has one slice more on top and bottom, respectively..
			gradTiff.setPosition(i+1);
			myGradientMask.setPosition(i+1);
			
			//get each images IPs
			ImageProcessor origIP = torigTiff.getProcessor();
			ImageProcessor gradIP = gradTiff.getProcessor();
			ImageProcessor maskIP = myGradientMask.getProcessor();		
	
			//because imageCalculator does not work... do it by hand.. again..
			for (int x = 0 ; x < origIP.getWidth(); x++) {
				for (int y = 0 ; y < origIP.getHeight(); y++) {
					if (maskIP.get(x, y) > 0 && origIP.get(x, y) > 0 && origIP.get(x, y) < myMode) {
						double nowWeight = gradIP.getPixelValue(x, y);
						double nowSample = origIP.getPixelValue(x, y);
						weightSum = weightSum + nowWeight;
						weightAndSampleSum = weightAndSampleSum + nowWeight * nowSample;
					}							
				}
			}
		}		
		
		double tMax = weightAndSampleSum / weightSum;
		
		return tMax;
	}

	public double macroPoreVolume(ImagePlus nowTiff) {
		
		double macroPoreVolume = 0;
		
		int nslices = nowTiff.getNSlices();
		
		for (int i = 0 ; i < nslices ; i++) {
			
			nowTiff.setPosition(i+1);
			ImageProcessor myIP = nowTiff.getProcessor();
			int[] nowHist = myIP.getHistogram();
			
			macroPoreVolume += nowHist[255];
						
		}
		
		return macroPoreVolume;
	}

	public ImagePlus findTheLargestCluster(ImagePlus nowTiff, int indexOfLargestCluster) {
		
		int i;
		ImageStack lStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImagePlus lTiff = new ImagePlus();
		
		//nowTiff.draw();
		//nowTiff.show();
		
		for (i = 0 ; i < nowTiff.getNSlices() ; i++) {
			
			nowTiff.setPosition(i + 1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			IJ.showStatus("Isolating largest cluster in slice #" + (i + 1) + "/" + (nowTiff.getNSlices()+1));				
			
			ImageProcessor modIP = new ByteProcessor(nowTiff.getWidth(), nowTiff.getHeight());
			
			for (int x = 0 ; x < nowIP.getWidth() ; x++) {
				for (int y = 0 ; y < nowIP.getHeight() ; y++) {
					int nowPix = Math.round(nowIP.getPixelValue(x, y));					
					if (nowPix == indexOfLargestCluster) modIP.putPixelValue(x, y, 255);
					else modIP.putPixelValue(x, y, 0);					
				}			
			}								
			lStack.addSlice(modIP);
		}
		lTiff.setStack(lStack);
		
		return lTiff;
	}
	
	public class AnisotropyResults {
		
		public double[] x;
		public double[] y;
		public double[] z;
		
		public double[] xy;		
		public double[] yx;
		
		public double[] dxz;
		public double[] dyz;
		public double[] dxyz;
		public double[] dyxz;
		
		public double[] uxz;
		public double[] uyz;
		public double[] uxyz;
		public double[] uyxz;
		
		public double anisotropy;
		public int[] mainAlignment;			
		
	}
	
	public class PlusVertex extends Vertex {
		
		public ArrayList<PlusEdge> plusBranches = null;
		
		public PlusVertex() {
			
			super();
			
			this.plusBranches = new ArrayList<PlusEdge>();			
			
		}
		
		public ArrayList<PlusEdge> getPlusBranches() {
			return this.plusBranches;
		}
		
		public void setPlusBranches(ArrayList<PlusEdge> nowBranches) {
			this.plusBranches = nowBranches;
		}
		
	}
	
	public class PlusEdge extends Edge {
		
		private RollerCaster rC = new RollerCaster();
		public double v1Dist;
		public double v2Dist;
		public ArrayList<Double> distance = null;
		
		public double constrictionFactor; 
		public double bottleNeck;
		public double lengthReduction;
		
		public PlusVertex plusV1 = null;
		public PlusVertex plusV2 = null;
		
		public PlusEdge(Vertex v1, Vertex v2, ArrayList<Point> slabs, double length) {
			
			super(v1, v2, slabs, length);
			
			this.distance = new ArrayList<Double>();
			this.constrictionFactor = 0; 
			this.v1Dist = 0;
			this.v2Dist = 0;
			this.bottleNeck = 0;
			this.lengthReduction = 0;
			
		}
		
	}
	
	public class PlusGraph extends Graph {
		
		public ArrayList<PlusEdge> plusEdge = null;	
		public ArrayList<PlusVertex> plusVertex = null;
		public ArrayList<Vertex> topVerts = null;
		public ArrayList<Vertex> botVerts = null;
		
		public PlusGraph() {
			
			super();
			
			this.plusEdge = new ArrayList <PlusEdge>();
			this.plusVertex = new ArrayList <PlusVertex>();
			
			this.topVerts = new ArrayList<Vertex>();
			this.botVerts = new ArrayList<Vertex>();			
			
		}
		
		public ArrayList<PlusVertex> getPlusVertices() {
			return this.plusVertex;			
		}
		
	}
	
	public class ROIMorphoProps {
		
		public double roiBulkVolume;		
		public double phaseVolume;
		public double phaseVolumeFraction;
		public double surfaceArea;		
		public double meanCurvature;
		public double eulerNumber;
		
		public double gamma;
		
		public double averagePhaseDiameter;
		public double averageDistance2PhaseBoundary;
		public double fractionLessThan3MMFromPhaseBoundary;
		public double fractionMoreThan3MMFromPhaseBoundary;
		
		public double surfaceFractalDimension;					
		
		public double anisotropy;		
		public int[] mainAlignment;
			
		public int phasePercolates;	
		public double criticalPhaseDiameter;		
		public double percolatingVolumeFraction;
		public double percolationThreshold;
		public double largesPhaseClusterVolumeFraction;
		
		public double volumeFractionConnected2Top;
		public double volumeFractionConnected2Bottom;
		
		public double depthOfPhasePenetration;

		public int[] distanceHistogram;		
		
	}
	
	public class DistanceProps {
		
		public double averageDistance2PhaseBoundaryInMM;
		public double medianDistance2PhaseBoundaryInMM;
		public double maxDistance2PhaseBoundaryInMM;
		public double modeDistance2PhaseBoundaryInMM;
		public double fractionLessThan3MMFromPhaseBoundary;
		public double fractionMoreThan3MMFromPhaseBoundary;
		
	}
	
	public class PoreClusterProps {
			
		public boolean hasItObjects;
		public boolean containsAParticleAnalysis;
		public int[] id;								//ID: unique particle identifier; this number is the label used for the particle in all calculations and output
		public double[] volume;							//Vol: particle volume
		public double[] surfaceArea;					//SA: surface area (0 if too small for mesh to be produced; see warning log)
		public double[] meanCurvature;
		public double[] euler;	

		public boolean[] isPercolating;
		public boolean[] touchesTop;
		public boolean[] touchesBot;
		
		public double[] xCenter;						//x Cent: x-coordinate of particle centroid
		public double[] yCenter;						//y Cent: y-coordinate of particle centroid
		public double[] zCenter;						//z Cent: z-coordinate of particle centroid
		
		public int[] minZ;
		public int[] maxZ;
		public int[] minX; 
		public int[] maxX;
		public int[] minY;
		public int[] maxY;
	}

	public class ParallelFloydWarshall {
	  
		//Ross Anderson's parallel Floyd-Warshall algorithm
		
		private ExecutorService exec;		  
		private int numThreads;
		private double[] current;
		private double[] next;
		
		private double[] currentD;
		private double[] nextD;
		
		private int[] currentP;
		private int[] nextP;
		  
		private int[] maxIndex;
		private int numNodes;
		private boolean solved;
		
		private double[][] outAdjacency;
		private double[][] outDistance;
		private int[][] outPrede;
		
		private int getIndex(int i, int j){
			return i*numNodes+j;
	  	}
	  
	  	private int getI(int index){
	  		return index / numNodes;
	  	}
	  
	  	private int getJ(int index){
	  		return index % numNodes;
	  	}
	  
	  	/**
	   	* @param numNodes the number of nodes in the graph
	   	* @param distances the matrix of distances between nodes, indexed from 0 to
	   	*                  numNodes-1.  distances[i][j] cost of a directed edge from
	   	*                  i to j.  Must be Double.POSITIVE_INFINITY if the edge is
	   	*                  not present.  distance[i][i] is a self arc (allowed).
	   	*/
	  	public ParallelFloydWarshall(int numNodes, double[][] adjacencyMatrix, double[][] distanceMatrix, int[][] predecessorMatrix,
	  								ExecutorService exec, int numThreads) {
	  		
		  	this.exec = exec;
		  	this.numThreads = numThreads;
		  	this.numNodes = numNodes;
		  	this.current = new double[numNodes*numNodes];
		  	this.next = new double[numNodes*numNodes];
		  	this.currentD = new double[numNodes*numNodes];
		  	this.nextD = new double[numNodes*numNodes];
		  	this.currentP = new int[numNodes*numNodes];
		  	this.nextP = new int[numNodes*numNodes];
		  	this.maxIndex = new int[numNodes*numNodes];
		  	Arrays.fill(maxIndex, -1);
		  	for(int i = 0; i < numNodes; i++){
			  	for(int j = 0; j < numNodes; j++){
				  	current[getIndex(i,j)] = adjacencyMatrix[i][j];
				  	currentD[getIndex(i,j)] = distanceMatrix[i][j];
				  	currentP[getIndex(i,j)] = predecessorMatrix[i][j];
			  	}
		  	}
		  	this.solved = false;		  	
		  	
		  	this.outAdjacency = adjacencyMatrix;
		  	this.outDistance = distanceMatrix;
		  	this.outPrede = predecessorMatrix;
	  	}
	  
	  	public double[][] getAdjacencyMatrix() {
	  		
	  		for (int k = 0 ; k < numNodes*numNodes ; k++) {
	  			
			  	int i = getI(k);
			  	int j = getJ(k);
	  			
			  	outAdjacency[i][j] = current[k];
	  		}
	  		
	  		return outAdjacency;
	  		
	  	}
	  	
	  	public double[][] getDistanceMatrix() {
	  		
	  		for (int k = 0 ; k < numNodes*numNodes ; k++) {
	  			
			  	int i = getI(k);
			  	int j = getJ(k);
	  			
			  	outDistance[i][j] = currentD[k];
	  		}
	  		
	  		return outDistance;
	  		
	  	}
	  	
	  	public int[][] getPredecessorMatrix() {
	  		
	  		for (int k = 0 ; k < numNodes*numNodes ; k++) {
	  			
			  	int i = getI(k);
			  	int j = getJ(k);
	  			
			  	outPrede[i][j] = currentP[k];
	  		}
	  		
	  		return outPrede;	  		
	  	}
	  	
	  	public void solve(){
		  	if(solved){
			  	throw new RuntimeException("Already solved");
		  	}
		  	for(int k = 0; k < numNodes; k++){
			  	List<Callable<Boolean>> tasks = new ArrayList<Callable<Boolean>>();
			  	if(current.length < numThreads){
				  	for(int i = 0; i < current.length; i++){
					  	tasks.add(new FloydJob(i,i+1,k));
				  	}
			  	}
			  	else{
				  	for(int t = 0; t < numThreads; t++){
				  		int calcsPerThread = current.length/numThreads;
					  	int lo = t * calcsPerThread;
					  	int hi = (t+1) * calcsPerThread;
					  	tasks.add(new FloydJob(lo,hi,k));
				  	}
			  	}
			  	try {
				  	List<Future<Boolean>> results = this.exec.invokeAll(tasks);
				  	for(Future<Boolean> result : results){
					  	if(!result.get().booleanValue()){
						  	throw new RuntimeException();
					  	}
				  	}
			  	} catch (InterruptedException e) {
				  	throw new RuntimeException(e);
			  	} catch (ExecutionException e) {
				  	throw new RuntimeException(e);
			  	}
			  	double[] temp = current;
			  	double[] tempD = currentD;
			  	int[] tempP = currentP;
			  	current = next;
			  	currentD = nextD;
			  	currentP = nextP;
			  	next = temp;
			  	nextD = tempD;  
			  	nextP = tempP;  
		  	}
		  	next = null;
		  	nextD = null;
		  	nextP = null;
		  	solved = true;
	  	}
	  
	  	/**
	  	 * 	
	  	 * @param i must lie in in [0,numNodes)
	  	 * @param j must lie in in [0,numNodes)
	  	 * @return the length of the shortest directed path from node i to node j.
	  	 *         If i == j, gives the shortest directed cycle starting at node i
	  	 *          (note that the graph may contain nodes with self loops).  Returns
	  	 *         Double.POSITIVE_INFINITY if there is no path from i to j.
	  	 */
	  	public double shorestPathLength(int i, int j){
	  		if(!solved){
	  			throw new RuntimeException("Must solve first");
	  		}
	  		return this.current[getIndex(i,j)];
	  	}
	  	/**
	  	 * Example: If the path from node 2 to node 5 is an edge from 2 to 3 and then
	  	 * an edge from 3 to 5, the return value will be
	  	 * Arrays.asList(Integer.valueOf(2),Integer.valueOf(3),Integer.valueOf(5));
	  	 * 
	  	 * @param i the start of the directed path
	  	 * @param j the end of the directed path
	  	 * @return The shortest path starting at node i and ending at node j, or null
	  	 *         if no such path exists.
	  	 */
	  	public List<Integer> shortestPath(int i, int j){    
		  	if(current[getIndex(i,j)] == Double.POSITIVE_INFINITY){
			  	return null;
		  	}
		  	else{
			  	List<Integer> ans = new ArrayList<Integer>();      
			  	ans.add(Integer.valueOf(i));
			  	shortestPathHelper(i,j,ans);
			  	return ans;
		  	}
	  	}
	  
	  	public void shortestPathHelper(int i, int j, List<Integer> partialPath){
		  	int index = getIndex(i,j);
		  	if(this.maxIndex[index] < 0){
			  	partialPath.add(Integer.valueOf(j));
		  	}
		  	else{
			  	shortestPathHelper(i,this.maxIndex[index],partialPath);
			  	shortestPathHelper(this.maxIndex[index],j,partialPath);
		  	}
	  	}
	  
	  	private class FloydJob implements Callable<Boolean>{
		  
		  	private final int lo;
		  	private final int hi;
		  	private final int k;
		  
		  	public FloydJob(int lo, int hi, int k){
			  	this.lo = lo;
			  	this.hi = hi;
			  	this.k = k;
		  	}
		  
		  	@Override
		  	public Boolean call() throws Exception {
			  	for(int index = lo; index < hi; index++){
				  	int i = getI(index);
				  	int j = getJ(index);
				  	double alternatePathValue = current[getIndex(i,k)]
	                                   		+ current[getIndex(k,j)];
				  	
				  	double alternatePathValueD = currentD[getIndex(i,k)]
                       		+ currentD[getIndex(k,j)];
				  	
				  	if(alternatePathValue < current[index]){
					 	next[index] = alternatePathValue;
					 	nextD[index] = alternatePathValueD;
					 	nextP[index] = currentP[getIndex(k,j)];
					 	maxIndex[index] = k;
				  	}
				  	else{
					  	next[index] = current[index];
					  	nextD[index] = currentD[index];
					  	nextP[index] = currentP[index];
				  	}
			  	}
			  	return true;
		  	}
	  	}
	}
	
	public ConnectingGraphs findConnectingGraphs(PlusGraph[] mPG, int numberOfLastSlice, int numberOfAddedSlices) {
		
		ConnectingGraphs mCG = new ConnectingGraphs();
		
		for (int i = 0 ; i < mPG.length ; i++) {
			
			PlusGraph nowGraph = mPG[i];
			
			ArrayList<PlusEdge> edgeList = nowGraph.plusEdge;
			ArrayList<Vertex> vertexList = nowGraph.getVertices();
			
			//also remember if a Vertex touches the top..
			ArrayList<Integer> touchesTop = new ArrayList<Integer>();
			ArrayList<Integer> touchesBottom = new ArrayList<Integer>();
			
			int cc=0;
			for (PlusEdge edge : edgeList) {
				
				cc++;
				
				Vertex v1 = edge.getV1();
				Vertex v2 = edge.getV2();
				
				if (v1 == null | v2 == null) IJ.error("Damn!!!");
				
				// use the index of the vertices as the index in the matrix			
				int row = vertexList.indexOf(v1);		
				
				if (vertexIsConnected2TopAndBottom(v1, numberOfLastSlice, numberOfAddedSlices) == 1) {			
					if(!touchesTop.contains(row)) {
						touchesTop.add(row);	
					}
				}			
				if (vertexIsConnected2TopAndBottom(v1, numberOfLastSlice, numberOfAddedSlices) == 2) {
					if (!touchesBottom.contains(row)) {
						touchesBottom.add(row);						
					}
				}
				
				int column = vertexList.indexOf(v2);			
				if (vertexIsConnected2TopAndBottom(v2, numberOfLastSlice, numberOfAddedSlices) == 1) {
					if (!touchesTop.contains(column)) {
						touchesTop.add(column);						
					}
				}
				if (vertexIsConnected2TopAndBottom(v2, numberOfLastSlice, numberOfAddedSlices) == 2) {
					if (!touchesBottom.contains(column)) {
						touchesBottom.add(column);						
					}
				}
			}
			
			if (touchesTop.size() > 0 & touchesBottom.size() > 0) mCG.percolating.add(i);
			if (touchesTop.size() > 0 & touchesBottom.size() == 0) mCG.con2Top.add(i);
			if (touchesTop.size() == 0 & touchesBottom.size() > 0) mCG.con2Bot.add(i);
			
		}
		
		return mCG;
		
	}
	
	public ROIMorphoProps getDistances2AirFilledPores(InputOutput.MyFileCollection	mFC, ColumnRoi colRoi, MenuWaiter.ROISelectionOptions mRSO) {
		
		MorphologyAnalyzer morph = new MorphologyAnalyzer();
		
		//init output structure
		MorphologyAnalyzer.ROIMorphoProps myR = morph.new ROIMorphoProps();
		
		//colRoi.nowTiff.updateAndDraw();colRoi.nowTiff.show();
		
		//add black invert ROI
		ImageStack distStack = new ImageStack(colRoi.nowTiff.getWidth(), colRoi.nowTiff.getHeight());
		ImageProcessor whiteOne = colRoi.nowTiff.getProcessor().duplicate();
		whiteOne.fill();
		distStack.addSlice(whiteOne);
		//disp.displayIP(whiteOne, "test");
		for (int z = 1 ; z < colRoi.nowTiff.getNSlices() + 1; z++) {
			colRoi.nowTiff.setPosition(z);
			ImageProcessor nowIP = colRoi.nowTiff.getProcessor().duplicate();
			nowIP.invert();
			distStack.addSlice(nowIP);
		}
		
	    //monitor memory
	    Runtime rt = Runtime.getRuntime();
	    long currInUse = (rt.totalMemory() - rt.freeMemory()) / 1000000000;  //in GB		
		
		//distTiff.updateAndDraw();distTiff.show();
		
		//calculate distance map
		IJ.showStatus("Performing Euklidean distance transform ...");
		
		short[] shortWeights = ChamferWeights3D.BORGEFORS.getShortWeights();
		boolean normalize = true;
		DistanceTransform3D dt = new DistanceTransform3DShort(shortWeights, normalize);
		distStack = dt.distanceMap(distStack);
				
		//remove topmost slice and set to distTiff
		distStack.deleteSlice(1);
		ImagePlus distTiff = new ImagePlus("DistImage", distStack);

		//distTiff.updateAndDraw();distTiff.show();
		
		//calc average distance to aerated pore
		currInUse = (rt.totalMemory() - rt.freeMemory()) / 1000000000;  //in GB		
		DistanceProps myDP = calculatePhaseDistanceProperties(distTiff, colRoi);
		myR.averageDistance2PhaseBoundary = myDP.averageDistance2PhaseBoundaryInMM;
		myR.fractionLessThan3MMFromPhaseBoundary = myDP.fractionLessThan3MMFromPhaseBoundary;
		myR.fractionMoreThan3MMFromPhaseBoundary = myDP.fractionMoreThan3MMFromPhaseBoundary;
		
		//delete distTiff
		distTiff.unlock();distTiff.flush();
		System.gc();System.gc();
		currInUse = (rt.totalMemory() - rt.freeMemory()) / 1000000000;  //in GB	
		
//		//get histogram of distances
//		ImageStack binStack = new ImageStack(colRoi.nowTiff.getWidth(), colRoi.nowTiff.getHeight(),colRoi.nowTiff.getNSlices());			
//		for (int z = 1 ; z < colRoi.nowTiff.getNSlices() ; z++) {  //skip first slice since it was added to sumulate the soil surface..
//			distTiff.setPosition(z + 1);
//			ImageProcessor thIP = distTiff.getProcessor();
//			ImageProcessor cpIP = thIP.duplicate();			
//			ImageProcessor binIP = cpIP.convertToByte(false);
//			binStack.setProcessor(binIP, z);		
//		}
//		ImagePlus binTiff = new ImagePlus("distanceBin", binStack);
		
		//binTiff.updateAndDraw();binTiff.show();
		
		//get average distance to aerated pore
		//myR.distanceHistogram = hist.extractHistograms8(binTiff);
		
		return myR;		
		
	}
	
	public ROIMorphoProps getSomeSimpleMorphoProps(InputOutput.MyFileCollection	mFC, ColumnRoi colRoi, MenuWaiter.ROISelectionOptions mRSO) {
		
		MorphologyAnalyzer morph = new MorphologyAnalyzer();
		InputOutput jIO = new InputOutput();
		BoundingBox3D mBB = new BoundingBox3D();
		
		//init output structure
		MorphologyAnalyzer.ROIMorphoProps myR = morph.new ROIMorphoProps();
		
		//analyze air phase
		FractalProperties myFP = morph.calculateFractalProperties(colRoi, mRSO);
		myR.surfaceFractalDimension = myFP.surfaceFractalDim;
		IJ.showStatus("Calculating macroporosities ...");
		
		//macroporosity
		myR.phaseVolume = morph.macroPoreVolume(colRoi.nowTiff);
		double[] bulkSoilVolume = {0, 0, 0, 0};
		ImagePlus surfTiff = new ImagePlus();		
		if (mRSO.choiceOfRoi.equalsIgnoreCase("RealSample") & mRSO.includeSurfaceTopography) {
			String[] myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC);
			surfTiff = jIO.openTiff3D(mFC.mySurfaceFolder + mFC.pathSep + myGandS[1]);
			if (mFC.myCutSurfaceFolder != null) {
				surfTiff = jIO.openTiff3D(mFC.myCutSurfaceFolder + mFC.pathSep + myGandS[1]);
			}
			bulkSoilVolume = morph.calculateTheBulkSoilVolume(myGandS[0], surfTiff, mRSO, mFC);
			myR.roiBulkVolume = bulkSoilVolume[0];			
		}
		else { 
			bulkSoilVolume[0] = Math.round(mRSO.areaOfInterest) * (double)colRoi.nowTiff.getNSlices();
			myR.roiBulkVolume = bulkSoilVolume[0];
			surfTiff = null;
		}				
		myR.phaseVolumeFraction = myR.phaseVolume / myR.roiBulkVolume;	
		
		IJ.showStatus("Identifying connected pore-clusters ...");

		FloodFillComponentsLabeling3D myFCCL = new FloodFillComponentsLabeling3D(26, 32);
		ImageStack myLabelStack = myFCCL.computeLabels(colRoi.nowTiff.getStack());
										
		IJ.freeMemory();IJ.freeMemory();
		
		//calculate surface and Euler number.. using MorphoLibJ
		ImagePlus labelTiff = new ImagePlus();
		labelTiff.setStack(myLabelStack);				
		
		int[] labels = LabelImages.findAllLabels(myLabelStack);
		int numOfObjects = labels.length;
		Calibration calib = labelTiff.getCalibration();
		
		//also find bounding boxes 
		Box3D[] myBox = mBB.analyzeRegions(myLabelStack, labels, calib);
		
		//labelTiff.updateAndDraw();
		//labelTiff.show();
		
		//surface and curvature
		IntrinsicVolumesAnalyzer3D mIRA = new IntrinsicVolumesAnalyzer3D();
		mIRA.setConnectivity(26);
		Result[] myMorphos = mIRA.analyzeRegions(myLabelStack, labels, calib);
		
		double[] morphoVols = new double[myMorphos.length];
		double[] morphoSurf = new double[myMorphos.length];
		double[] morphoCurv = new double[myMorphos.length];
		double[] morphoEuler = new double[myMorphos.length];
		for (int i = 0 ; i < myMorphos.length ; i++) {
			morphoVols[i] = myMorphos[i].volume;
			morphoSurf[i] = myMorphos[i].surfaceArea;
			morphoCurv[i] = 2 * Math.PI * myMorphos[i].meanBreadth;
			morphoEuler[i] = myMorphos[i].eulerNumber;
		}
					
		//myP.phaseVolume = StatUtils.sum(morphoVols);
		myR.surfaceArea = StatUtils.sum(morphoSurf);
		myR.meanCurvature = StatUtils.sum(morphoCurv);			
		myR.eulerNumber = StatUtils.sum(morphoEuler);	
		
		//calculate percolating clusters ... needs the other parameters.. and is needed for following parameters	
		//calculate percolating clusters ... needs the other parameters.. and is needed for following parameters	
		ArrayList<Integer> conTop = new ArrayList<Integer>();
		ArrayList<Integer> conBot = new ArrayList<Integer>();
		ArrayList<Integer> isPercolating = new ArrayList<Integer>();
		if (numOfObjects > 0) {
			conTop = check4TouchingTheTop(labelTiff, colRoi.nowTiff, surfTiff);
			conBot = check4TouchingTheBottom(labelTiff, colRoi.nowTiff, surfTiff);
			isPercolating = checkForPercolatingClusters(numOfObjects, conTop, conBot);				
		}	
		
		//calculate percolating phase fraction
		double percolatingVolume = 0;
		for (int i : isPercolating) percolatingVolume += morphoVols[i-1];   
		myR.percolatingVolumeFraction = percolatingVolume / myR.roiBulkVolume;
		
		//check if volume percolates
		myR.phasePercolates = 0;
		if (myR.percolatingVolumeFraction > 0) myR.phasePercolates = 1;
		
		//calculate fraction of largest pore cluster
		double volOfLC = StatUtils.max(morphoVols);
		myR.largesPhaseClusterVolumeFraction = volOfLC / myR.roiBulkVolume;
		
		//calculate depth of phase penetration		
		double penetrationDepth = 0;
		for (int i = 0 ; i < myBox.length ; i++) if (myBox[i].getZMax() > penetrationDepth) penetrationDepth = myBox[i].getZMax();
		myR.depthOfPhasePenetration = penetrationDepth;
		
		//calculate global connection probability (gamma)
		double sumOfSquaredClusterSizes = 0;
		double rF = 1000000;
		for (int i = 0 ; i < myMorphos.length ; i++) {
			double squaredClusterSize = morphoVols[i] / rF * morphoVols[i] / rF;
			sumOfSquaredClusterSizes += squaredClusterSize;
		}
		myR.gamma = sumOfSquaredClusterSizes / (myR.phaseVolume / rF * myR.phaseVolume / rF);	
		
		return myR;
		
	}
}


