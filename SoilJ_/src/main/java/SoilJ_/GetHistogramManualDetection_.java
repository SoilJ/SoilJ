package SoilJ_;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

// import classes from SoilJ
import SoilJ.tools.InputOutput;
import SoilJ.tools.ObjectDetector;
import SoilJ.tools.RoiHandler;
import ij.IJ;
// import ImageJ and some libraries
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.EllipseFitter;
import ij.process.ImageProcessor;

// import other Java libraries


public class GetHistogramManualDetection_ extends ImagePlus implements PlugIn {

	/**
	 * Function to get the histogram of the column wall in order to calibrate the images (slice by slice)
	 * therefore, each slice has to be calibrated using a  gradient of greyvalues (air + wall material, even if it has holes) from top to bottom
	 * 
	 */
	@Override
	public void run(String arg) {
		// get the seaparator of the operating system ( / \ etc. )
		String pathSep = System.getProperty("file.separator");
			
		// Construct InputOutput >> Includes menu for choosing images
		InputOutput jIO = new InputOutput();
		// STruct for 3D Coordinates
		ObjectDetector jOD = new ObjectDetector();
		// Struct for InnerCircle
		RoiHandler roi = new RoiHandler();
		RoiHandler.ColumnRoi colRoi = roi.new ColumnRoi();

		// read file or files
		InputOutput.MyFileCollection mFC = jIO.fileSelector("Please choose the folder with your image data");
		
		// look for innerCircle Files
		String myGaugeFolder = jIO.findTheInnerCircle(mFC.myBaseFolder);
		String[] myGauges = jIO.listInnerCircleFiles(myGaugeFolder, ".tif");			
		mFC.myInnerCircleFolder = myGaugeFolder;
		mFC.myInnerCircleFiles = myGauges;

		// create output folder for the histograms
		String myOutFolderName = "Histograms";
		String myOutFolder = mFC.myBaseFolder + mFC.pathSep + myOutFolderName;
		new File(myOutFolder).mkdir();
		mFC.myOutFolder = myOutFolder;
		
		
		// loop through all the images / innerCircle Files
		for (int i = 0; i < mFC.myTiffs.length; i++) {
			
			// assign tiff file
			mFC.fileName = mFC.myTiffs[i];;
			mFC = jIO.addCurrentFileInfo(mFC);
			String filenameShort = mFC.fileName.replace(".tif", "");
			
			// load (current) image
			ImagePlus nowTiff = jIO.openTiff3D(mFC.myBaseFolder + mFC.pathSep + mFC.myTiffs[i]);
			
			// load corresponding InnerCircle File
			String[] GandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC);
			mFC.nowInnerCirclePath = mFC.myInnerCircleFolder + mFC.pathSep + GandS[0];
			
			//read InnerCircle file
			if (GandS[0].contains("Gauge")) {
				
				ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();				
				int versio = jIO.checkInnerCircleFileVersion(mFC.nowInnerCirclePath);			
				if (versio == 0) jCO = jIO.readInnerCircleVer0(mFC.nowInnerCirclePath);	
				else jCO = jIO.readInnerCircleVer1(mFC.nowInnerCirclePath);	
				colRoi.jCO = jCO;
			}			
			
			// first of al duplicate the image, so we can change it without changing the original
			ImagePlus inputTiff = new ImagePlus( nowTiff.getTitle(), nowTiff.getImageStack().duplicate() );
			System.out.println(inputTiff.getWidth());
			
			// ++++++++++++ DO THE MAGIC ++++++++++++++++++++
			// FUNCTION CREATING TWO ELLIPTIC POLYGONS FROM INNERCIRCLE AND SETTING EVERYTHIG ELSE THAN WALL TO 0		(colRoi, nowTiff)
			ImagePlus wallTiff = removeEverythingButWall(inputTiff, colRoi, myOutFolder, filenameShort);
			
			// FUNCTION EXTRACTING THE HISTOGRAM										(IMAGE WITH EVERYTHING ELSE THAN WALL 0)	
			
			// EVENTUALLY FUNCTION THAT ACTUALLY CALIBRATES THE IMAGES 					(HISTOGRAMS OR DERIVED VALUE, IMAGE
			
			
			
		} // loop through all images
	} // run
	
	
	private ImagePlus removeEverythingButWall (ImagePlus inputTiff, RoiHandler.ColumnRoi colRoi, String myOutFolder, String filenameShort) {
		
		ImageStack onlyWall =  new ImageStack(inputTiff.getWidth(), inputTiff.getHeight()); 
		ImagePlus wallTiff = new ImagePlus();
		
		// in each slice of the image: create the corresponding ellipse and set everything outside to 0 / everything inside to 0
		for (int z = 0; z < inputTiff.getNSlices(); z++) {
			
			inputTiff.setSlice(z+1); 										// no slice 0; start at 1 
			ImageProcessor ipNow = inputTiff.getProcessor();				// get the 2D image at slice z+1
			
			EllipseFitter outerEllipse = new EllipseFitter();				// OUTER ELLIPSE of InnerCircle
			outerEllipse.major = colRoi.jCO.outerMajorRadius[z]*2;	
			outerEllipse.minor = colRoi.jCO.outerMinorRadius[z]*2;
			outerEllipse.xCenter = colRoi.jCO.xmid[z];
			outerEllipse.yCenter = colRoi.jCO.ymid[z];
			outerEllipse.theta = colRoi.jCO.theta[z]+ Math.PI / 2.;	 		// (outerEllipse.angle would be degrees)
			outerEllipse.makeRoi(ipNow);
			
			EllipseFitter innerEllipse = new EllipseFitter();				// INNER ELLIPSE of InnerCircle
			innerEllipse.major = colRoi.jCO.innerMajorRadius[z]*2;
			innerEllipse.minor = colRoi.jCO.innerMinorRadius[z]*2;
			innerEllipse.xCenter = colRoi.jCO.ixmid[z];
			innerEllipse.yCenter = colRoi.jCO.iymid[z];
			innerEllipse.theta = colRoi.jCO.itheta[z]+ Math.PI / 2.;	 	// (innerEllipse.angle would be degrees)
			innerEllipse.makeRoi(ipNow);

			outerEllipse.makeRoi(ipNow);									// create the ROI in the RGB ImageProcessor	
			PolygonRoi outerP = new PolygonRoi(outerEllipse.xCoordinates, outerEllipse.yCoordinates, outerEllipse.nCoordinates, Roi.POLYGON);  // create a polygon Roi
			ImageProcessor mask = new ByteProcessor(ipNow.getWidth(), ipNow.getHeight()).convertToShort(true);
			mask.set(0);
			mask.setRoi(outerP);
			mask.setColor(1);
			mask.fill(outerP);
			
			PolygonRoi innerP = new PolygonRoi(innerEllipse.xCoordinates, innerEllipse.yCoordinates, innerEllipse.nCoordinates, Roi.POLYGON);  // create a polygon Roi
			mask.setRoi(innerP);
			mask.setColor(0);
			mask.fill(innerP);

			onlyWall.addSlice(mask);										// add the image to the stack

		} // loop through slices
		wallTiff.setStack(onlyWall);										// create 3d ImagePlus from image stack; either 1 for wall or 0 for non-wall
		
		ImageCalculator calc = new ImageCalculator();
        ImagePlus Wall = calc.run("Multiply create stack", wallTiff, inputTiff); // 16-bit image with all 0 except for wall with its original values (original image multiplied by either 1 or 0
        
//        IJ.saveAsTiff(Wall, myOutFolder + System.getProperty("file.separator") +  filenameShort + "tif");
        
        
        // For each slice in "Wall", create a histogram (slice, value, frequency)
        List<int[]> histlist = new ArrayList<>();
        
        for (int z = 0; z < Wall.getNSlices(); z ++) { // from top to bottom
        	Wall.setSlice(z+1); 											// no slice 0; start at 1 
        	
        	ImageProcessor ipWall = Wall.getProcessor();					// get the slices' 2D image
        	int[] sliceHisto = ipWall.getStats().histogram16;//ipWall.getStatistics().histogram16;
        	sliceHisto = Arrays.copyOfRange(sliceHisto, 1, sliceHisto.length+1);
        	histlist.add(sliceHisto);
        	
        }
        
        System.out.println(histlist.get(1).length);
		
//         later it would be better to write (and append to) file while creating the new slices, that would probably save a lot of time =)
//         but then it should also become an option to either save the histograms or not... because actually we don't really need them, just refer to the list in the calibration
        try (
        		FileWriter wF = new FileWriter(myOutFolder + System.getProperty("file.separator") + filenameShort + ".csv", false); 
        		BufferedWriter wFB = new BufferedWriter(wF);
        	) {
        	String nowLine = "";
			IJ.showStatus("Saving histogram of file " + filenameShort + " into histogram folder..");
			
			for (int i = 1; i < Math.pow(2,16)-1; i++) { // 
				for (int s = 0; s < Wall.getNSlices(); s++) {
					nowLine += Integer.toString(histlist.get(s)[i]) + "\t";
				}
				wFB.write(nowLine + "\n");
				nowLine = "";
			}
			wFB.flush();
			wFB.close();
			System.out.println("saved histogram =)");
        }
		catch(Exception e){
			IJ.error("Exception " + e.toString() + "\nSaving path: " + myOutFolder);
		}
        
        

		// GET THE SPLINE OF THE HISTOGRAM
        int maxF = 0;
        int maxG = 0;
        int minF = inputTiff.getWidth() * inputTiff.getHeight() * inputTiff.getNSlices();
        int minG = inputTiff.getWidth() * inputTiff.getHeight() * inputTiff.getNSlices();
        
        for (int[] slicevalues : histlist) {
        	int localMax = Arrays.stream(slicevalues).summaryStatistics().getMax();  // maximum in slice
        	if (localMax > maxF) {maxF = localMax;}
        	int localMin = Arrays.stream(slicevalues).summaryStatistics().getMin(); // minimum in slice
        	if (localMin < minF ) {minF = localMin;}
        	
        	int i = 0; int j = 1;
        	while (i == 0) {
        		i += slicevalues[j];
        		j += 1;
        	}
        	int q = slicevalues.length-1; j = slicevalues.length-1;
        	while(q == slicevalues.length-1) {
        		q += slicevalues[j];
        		j = j-1;
        	}
        	if (i > minG) {minF = i;}
        	if (q < maxG) {maxG = q;}

        }
        double range = maxF - minF;
        int min = minF;
        
        
        // also normalise grey values (max = largest with freq > 0, min = smallest with freq > 0 but not 0
        // not necessary, because anyawys almost 0 and 2^16...
        
        
        List<double[]> mostfrequent = new ArrayList<>();
        for (int i = 0; i < inputTiff.getNSlices(); i++) {
//        	int [] temp = histlist.get(i);
        	double [] temp = Arrays.stream(histlist.get(i)).asDoubleStream().toArray();  // must be double, as int everything would be floored to 0
        	double[] normal = Arrays.stream(temp).map(p -> (p - min )/range).toArray();  // normalize the array to values 0-1 
        	// loop through array of this slice to get all values > .1
        	for(int j = 6500; j < normal.length; j++) {
                if(normal[j] >= 0.1) {
                   double [] a = {i, j, normal[j]};   // ca. 6553.6 muss der Threshold sein...
                   mostfrequent.add(a);
                }
            }
        }
//        System.out.println(mostfrequent.size());
        
        double [] x_slice = new double[mostfrequent.size()];
        double [] y_grey = new double[mostfrequent.size()];
        double [] w_freq = new double[mostfrequent.size()];
        
        for (int i = 0; i < mostfrequent.size(); i++) {
        	x_slice[i] = (double) mostfrequent.get(i)[0];
        	y_grey[i] = (double) mostfrequent.get(i)[1];
        	w_freq[i] = (double) mostfrequent.get(i)[2];
        }
        
        
        ArrayList<GreyCoords> greylist = new ArrayList<GreyCoords>();
        for (int i = 0; i < mostfrequent.size(); i++) {
        	GreyCoords g = new GreyCoords();
        	g.slice = (int) mostfrequent.get(i)[0];
        	g.greyvalue = mostfrequent.get(i)[1];
        	greylist.add(g);
        }
        
        Map<Integer, Double> avgGrey = 
        	    greylist.stream()
        	          .collect(Collectors.groupingBy(e -> e.getSlice(),
        	                                         Collectors.averagingDouble(GreyCoords::getGreyvalue)));
        TreeMap<Integer, Double> sorted = new TreeMap<Integer, Double>(avgGrey);

        LoessInterpolator loess = new LoessInterpolator(.5, 20); // BANDWIDTH: Default 0.3, robustnessIters: default 2
        double [] xvals = new double[sorted.size()];
        double [] yvals = new double[sorted.size()];
        int c = 0;
        for(Map.Entry<Integer, Double> entry : sorted.entrySet()){
        	xvals[c] = entry.getKey().doubleValue();
            yvals[c] = entry.getValue();
            c++;
        }
        PolynomialSplineFunction b = loess.interpolate(xvals , yvals);
      
        System.out.println(sorted.lastKey()-sorted.firstKey());
        double [][] loesscurve = new double[sorted.lastKey()-sorted.firstKey()][2];
                
        for (int i = 0 ; i < sorted.lastKey()-sorted.firstKey(); i++) {
        	int n = i + sorted.firstKey();
        	loesscurve[i][0] = n;
        	loesscurve[i][1] = b.value(n);
        }
        
        
//        LoessInterpolator loess = new LoessInterpolator();
//        double [] xvals = avgGrey.keySet().toArray(new double[avgGrey.size()]);
//        double [] yvals = avgGrey.values().toArray(new double[avgGrey.size()]);
//        PolynomialSplineFunction b = loess.interpolate(  xvals  , yvals);
//        
        
//        CurveFitter curve = new CurveFitter(x_slice, y_grey);
//        curve.doFit(17); // 16: POLY6; 2: POLY3; 3: POLY4; 17: POLY5; 23: Chapman; 24: ERF; 8: GAMMA_VARIATE; 14: INV_RODBARD
//        curve.POLY3();
        
//       LoessInterpolator loess = new LoessInterpolator();
//       // PolynomialSplineFunction
//       double[] a = loess.smooth(x_slice, y_grey);
//       PolynomialSplineFunction b = loess.interpolate(x_slice, y_grey);  // muss wahrscheinlihc loess.interpolate(smooth) sein
       
       
        
//        System.out.println("\nSpline: ");
//        System.out.println(b.value(500.));
//        System.out.println(b.value(0.));
//        
//        double[][] interp = new double[inputTiff.getNSlices()][2];
//        for (int i = 0; i < inputTiff.getNSlices(); i++) {
//        	interp[i][0] = i;
//        	interp[i][1] = curve.f(i);
//        }
//        
        try (
        		FileWriter wF2 = new FileWriter(myOutFolder + System.getProperty("file.separator") + "LoessSpline_" + filenameShort + ".csv", false); 
        		BufferedWriter wFB2 = new BufferedWriter(wF2);
        	) {
        	String nowLine = "";
			for (int i = 1; i < loesscurve.length; i++) { // 
				nowLine += String.valueOf(loesscurve[i][0]) + "\t" + String.valueOf(loesscurve[i][1]);
				wFB2.write(nowLine + "\n");
				nowLine = "";
			}
			wFB2.flush();
			wFB2.close();
			System.out.println("saved profile =)");
        }
		catch(Exception e){
			IJ.error("Exception " + e.toString() + "\nSaving path: " + myOutFolder);
		}
//        System.out.println(mostfrequent.get(0)[0]);
//        System.out.println(mostfrequent.get(0)[1]);
//        System.out.println(mostfrequent.get(0)[2]);
        System.out.println("\nYay!");
        
        
        

        return Wall;
	}
	
	public class GreyCoords {
		public int slice;
		public double greyvalue;
		
		public int getSlice() {
			return slice;
		}
		public double getGreyvalue() {
			return greyvalue;
		}
	} // GreyCoords


	
} // class GetHistogramManualDetection_




