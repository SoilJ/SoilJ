package SoilJ_;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.Panel;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;

import SoilJ.tools.DisplayThings;
import SoilJ.tools.ImageManipulator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.InputOutput.MyFileCollection;
import SoilJ.tools.ObjectDetector;
import SoilJ.tools.TailoredMaths;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.ImageCanvas;
import ij.gui.Roi;
import ij.gui.StackWindow;
import ij.gui.Toolbar;
import ij.plugin.ContrastEnhancer;
import ij.plugin.PlugIn;


// ENTER COLUMN WALL THICKNESS (IF CONSTANT! OR ONLY DRAW OUTER LINE IN TOP AND BOTTOM SLICE AND GET THICKNESS FORM THERE)

/**
 * DefineColumnOutlinesManually is a SoilJ plugin for manual column selection.
 * The (inner & outer) outlines of the column are drawn using an ellipse. The inner
 * outline can also be determined based on the column thickness (only one ellipse has to be drawn)
 * 
 * Future: A side view of the column allows you to determine top and bottom of the column.
 * For now, the uppermost drawn ellipse corresponds to the top, the lowermost ellipse corresponds to the bottom of the column.
 * 
 * For now: Do this for each image separately
 * 
 * @author John Koestel & Maria Vorkauf
 *
 */

public class DefineColumnOutlinesManually2_ extends ImagePlus implements PlugIn {

	@Override
	public void run(String arg) {
		// get the seaparator of the operating system ( / \ etc. )
		String pathSep = System.getProperty("file.separator");
		
		// Construct InputOutput >> Includes menu for choosing images
		InputOutput jIO = new InputOutput();		

		// read file or files
		InputOutput.MyFileCollection mFC = jIO.fileSelectorOneFile("Please choose your image data");

		// create output folder
		String myOutFolderName = "WallCoordinateDefined";
		String myOutFolder = mFC.myBaseFolder + mFC.pathSep + myOutFolderName;
		new File(myOutFolder).mkdir();
		mFC.myOutFolder = myOutFolder;

		// create innerCircle Folder
		String myInnerCircleFolder = myOutFolder + mFC.pathSep + "InnerCircle";
		new File(myInnerCircleFolder).mkdir();
		mFC.myInnerCircleFolder = myInnerCircleFolder;

		// loop over 3D images
		for (int i = 0; i < mFC.myTiffs.length; i++) {
			
			// assign tiff file
			mFC.fileName = mFC.myTiffs[i];
			String filenameShort = mFC.fileName.replace(".tif", "");
			mFC = jIO.addCurrentFileInfo(mFC);

			// load image
			ImagePlus nowTiff = jIO.openVirtualStack3D(mFC.myBaseFolder + mFC.pathSep + mFC.myTiffs[i]);
			nowTiff.show();
			inputImage = nowTiff;			

			// Here I would like to open the image as virtual stack - 
			// only load the slices that are actually displayed. But doesn't work (yet)
			// ++ImagePlus nowTiff = jIO.openVirtualStack3D(mFC.myBaseFolder + mFC.pathSep + mFC.myTiffs[i]);
			// ++inputImage = nowTiff;
			
			inputStackCopy = inputImage.getImageStack().duplicate(); 	// duplicate the stack
			displayImage = new ImagePlus( inputImage.getTitle(), 		// generate the image to display (imagePlus)
								inputStackCopy );

			displayImage.setTitle(filenameShort);						// title to show
			displayImage.setColor(Color.CYAN);
			
			displayImage.setSlice( inputImage.getCurrentSlice() );		// set the slice of displayImage to the inputImage slice
			new ContrastEnhancer().stretchHistogram(displayImage, 0.5);	// adapt contrast (for the current slice, so sometimes when scrolling through there are "bumps")

			// set the 2D flag
			inputIs2D = inputImage.getImageStackSize() == 1;			// boolean whether 2D yes or no

			// correct Fiji error when the slices are read as frames
			
			if ( inputIs2D == false && 
					displayImage.isHyperStack() == false && 
					displayImage.getNSlices() == 1 )
			{
				// correct stack by setting number of frames as slices
				displayImage.setDimensions( 1, displayImage.getNFrames(), 1 );
			}
			
			
			//mFC hand-over to the action listener does only work as final! >>> re-initiate
			final InputOutput.MyFileCollection mFCfinal = mFC;
			// Build GUI
			
			// for each of the images: open the CustomWindow (below)
			SwingUtilities.invokeLater(
				new Runnable() {
					public void run() {	
						win = new CustomWindow( displayImage, myInnerCircleFolder, pathSep, jIO, 
								filenameShort, myOutFolder , mFCfinal);
						win.pack();
					} // run
				} // runnable
			); // SwingUtilities
		
		} // loop through all tiffs
	} // run
	
	
	// initiate variables, image(holders), windows, panels, buttons, etc. 
//	private CustomWindow win;					// Main GUI crate a CustomWindow (as defined below) which can then be filled with buttons, images, etc.

	boolean inputIs2D = false;					// set 2D to false, can be overwritten if actually 2D
	
	ImagePlus inputImage = null;				// initiate for Original input
	ImageStack inputStackCopy = null;			// initiate for Copy of Input (so to not overwrite)
	ImagePlus displayImage = null;				// initiate for Image displayed in GUI

	Panel mainPanel = new Panel();  			// Main Panel in GUI >> where we display the image
	JPanel paramsPanel = new JPanel();			// Parameters panel (segmentation + display options) ??????
	JPanel featureSelectionPanel = new JPanel(); // feature selection Panel (for rois?)
	JPanel radioPanel = new JPanel( new GridLayout(0, 1) ); //create a panel for storing the radio buttons
	JLabel inputImagePicture;					// 
	JPanel buttonsAndPicPanel = new JPanel();	// Panel for the image + the buttons together (will contain image panel and button panel in larger panel)
	JPanel roiCoordinatesPanel = new JPanel(); 	// create panel for ROI coordinates
	JPanel innerCirclePanel = new JPanel();		// create Panel for the InnerCircleButton
	
	ButtonGroup inputImageButtons;				// input options (depend on selected buttons)			
	JRadioButton sampleButton;					//radio-button for selecting soil + string for name
	static String sampleImageText = "Soil";
	JRadioButton wallButton;					// radio-button for selecting wall + string for name
	static String wallImageText = "Wall";
	JRadioButton calibrationButton;				// radio-button for selecting calibration (???) + string for name
	static String calibrationImageText = "Calibration";
	
	
		
	JButton addROIButton;						//Button to add ROI + text of button
	JCheckBox wallThickKnown;					//Checkbox to only draw one ellipse (but then enter wall thickness
	JTextField wallThicknessBox;					//Where to enter the wall thickness if known

	ArrayList<RoiCoords> rois = new ArrayList<RoiCoords>(); // initialize struct to store the outlines/ellipses (=ROIs)   
	
	// executor service to launch threads for the plugin methods and events
	final ExecutorService exec = Executors.newFixedThreadPool(1);
	
/** == > go create !! */
	
	private class CustomWindow extends StackWindow
	{
		/** serial version uid */  // ????????
		private static final long serialVersionUID = -6201855439028581892L;

		/** Listener for the GUI buttons */
		private ActionListener listener = new ActionListener() {

			@Override
			public void actionPerformed( final ActionEvent e ) {
				exec.submit(new Runnable() {					// listen to the buttons on separate threads not to block the event dispatch thread
					public void run(){
						// "add ROI"  button		
//						if( e.getSource() == addROIButton ){	
//						}						
					} // run
				}); // Runnable
			} // actionPerformed
		}; // listener

		ObjectDetector jOD = new ObjectDetector();			// includes struct for 3D coordinates
		
		
		/** Construct the plugin window
		 * @param imp input image [ImagePlus]
		 * @param ICpth Path to inner circle files [String]
		 * @param PathSep System path separator [String]
		 * @param jIO InputOutput Object [SoilJ.tools.InputOutput]
		 * @param filename current tiff file name (should acutally be taken out again, now with mFC) [String]
		 * @param mFC organizational infos [SoilJ.tools.InputOutput.MyFileCollection]
		 * @apiNote *Plugin window (GUI) to implement inner circle files manually*
		 * */
		CustomWindow( ImagePlus imp, String ICpth, String pathSep, 
				InputOutput jIO, String filename, 
				String myOutFolder, MyFileCollection mFC) {
			super(imp, new ImageCanvas(imp));
			
			final ImageCanvas canvas = (ImageCanvas) getCanvas(); 				//ImageCanvas canvas which will show the image
			Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();	// get screen dimensions for appropriate zoom 
			double screenWidth = screenSize.getWidth();
			double screenHeight = screenSize.getHeight();
			
			// Zoom in if image is too small
			while( ( ic.getWidth() < screenWidth/2 || ic.getHeight() < screenHeight/2 ) && ic.getMagnification() < 32.0 ){			
			    final int canvasWidth = ic.getWidth();
    			ic.zoomIn( 0, 0 );
    			// check if canvas size changed (otherwise stop zooming)
    			if( canvasWidth == ic.getWidth() ){
    				ic.zoomOut(0, 0);
    				break;
    			}
    		}
			// Zoom out if canvas is too large
			while( ( ic.getWidth() > 0.75 * screenWidth || ic.getHeight() > 0.75 * screenHeight ) && ic.getMagnification() > 1/72.0 ) {
				final int canvasWidth = ic.getWidth();
	    		ic.zoomOut( 0, 0 );
	    		// check if canvas size changed (otherwise stop zooming)
	    		if( canvasWidth == ic.getWidth() ){
	    			ic.zoomIn(0, 0);
	    			break;
	    		}
			}
			setTitle( filename );

			// select oval ROI tool to define outlines of features of interest (??? wie weiss es, dass es zu diesem Bild gehört?) Wo wird das Oval gezeichnet????
			// Mir noch nicht ganz klar, aber funktioniert soweit
			Toolbar.getInstance().setTool( Toolbar.OVAL );

	// FEATURE SELECTION PANEL ========================================================================

			// feature selection options (soil sample, inner wall perimeter or calibration rod)			
			// select inner circle of the wall
			sampleButton = new JRadioButton( sampleImageText );			
			sampleButton.setActionCommand( sampleImageText );
			sampleButton.addActionListener( listener );
			sampleButton.setToolTipText( "outline of soil sample" );
			sampleButton.setSelected(true);
			
			// select outer circle of the wall
			wallButton = new JRadioButton( wallImageText );
			wallButton.setActionCommand( wallImageText );
			wallButton.addActionListener( listener );
			wallButton.setToolTipText( "location of inner wall perimeter" );
			
			// select calibration thing location (???)
			calibrationButton = new JRadioButton( calibrationImageText );
			calibrationButton.setActionCommand( calibrationImageText );
			calibrationButton.addActionListener( listener );
			calibrationButton.setToolTipText( "location of calibration rod" );
			calibrationButton.setEnabled(false);		// vorerst mal deaktivieren

			inputImageButtons = new ButtonGroup();
			inputImageButtons.add( sampleButton );
			inputImageButtons.add( wallButton );
			inputImageButtons.add( calibrationButton );
			radioPanel.add( sampleButton );
			radioPanel.add( wallButton );
			radioPanel.add( calibrationButton );

			buttonsAndPicPanel.add(radioPanel);
			
			
			
			// add InnerCircle Button
			innerCirclePanel.setLayout(new BoxLayout(innerCirclePanel, BoxLayout.LINE_AXIS));
			GridBagConstraints innerCircleConstraints = new GridBagConstraints();
			innerCircleConstraints.anchor = GridBagConstraints.SOUTHWEST;
			innerCircleConstraints.fill = GridBagConstraints.NONE;
			JButton innerCircleButton = new JButton("Create InnerCircle File");
			innerCircleButton.setEnabled(false);								// Initially inactive, you should first draw and submit a ROI
		
		// What happens if you press the inner circle button:	
			innerCircleButton.addActionListener(new ActionListener() {			// >> set what happens if you press button >> basically create ROI file
				public void actionPerformed(ActionEvent e) {
					if (rois.size() > 0) {  // first check whether there is something to export!
						ObjectDetector.ColCoords3D roiCol = jOD.new ColCoords3D();			// initiate ColCoords3d roiCol
	
						int[] zValues = new int[rois.size()];	// just store the z Values to be able to go from top to bottom, get column height etc. 
						int[] wtValues =  new int[rois.size()];
						for (int i = 0; i < rois.size(); i++) {zValues[i] = (int) rois.get(i).z; wtValues[i] = (int) rois.get(i).wt;} // populate the zValues array with the z values (two (soil + wall) per z)
						ArrayList<Integer> zValueswt0 = new ArrayList<>();
						for (int i = 0; i < zValues.length; i ++) {
							if (wtValues[i] == 0) {zValueswt0.add(zValues[i]);}
						}
						int[] zValues_wt = zValueswt0.stream().mapToInt(i -> i).toArray();
						
						int [] uniqueZ = Arrays.stream(zValues_wt).distinct().toArray();		// only take the unique values.
						
						// create a frequency map for each z. Each slice (z) should have exactly 2 ROIs: 1 soil + 1 wall
						HashMap<Integer, Integer> freqMap = new HashMap<Integer, Integer>();
						for (int i=0; i<zValues_wt.length; i++) {
				            if (freqMap.containsKey(zValues_wt[i])) {freqMap.put(zValues_wt[i], freqMap.get(zValues_wt[i]) + 1);}
				            else {freqMap.put(zValues_wt[i], 1);}
				        }
						// check if exactly two per slice, otherwise do not create innerCircle
						for (Map.Entry entry : freqMap.entrySet()) {
				            int freq = (int) entry.getValue();
							if (freq != 2) {
								JOptionPane.showMessageDialog(null, "Exactly two ROIs (soil + wall) in each marked slice needed!");
								return;
							}
				        }
						
						//get top: 
						
						uniqueZ = Arrays.stream(zValues).distinct().toArray();
						int maxZ = Arrays.stream(uniqueZ).max().getAsInt();					// max Z corresponds to column top
						int minZ = Arrays.stream(uniqueZ).min().getAsInt();					// min Z corresponds to column bottom
						int heightZ = maxZ - minZ + 1;
						
						// create an empty list to fill
						// float[][] coordMatrix = new float[displayImage.getNSlices()][15];
						float[][] coordMatrix = new float[heightZ][15];						// create a matrix in which to store the values (for each Z between top & bottom one line), 15 columns for ellipse definition
						for (int i = 0; i < heightZ; i++) {									// for each line
							coordMatrix[i][0] = i+1;											// select the matrix line
							for (int j = 0; j < rois.size(); j++) {							// for each submitted ellipse
								RoiCoords roiNow = rois.get(j);								// 
								if ((roiNow.z-minZ) == i) {
									if (roiNow.type == "Soil") { 							// soil: inner circle
										coordMatrix[i][2] = (float) roiNow.z;				// z (=slice)
										coordMatrix[i][5] = (float) roiNow.x;				// 
										coordMatrix[i][6] = (float) roiNow.y;

										coordMatrix[i][9] = (float) roiNow.w/2;
										coordMatrix[i][10] = (float) roiNow.h/2;

										coordMatrix[i][12] = (float) roiNow.theta;
										coordMatrix[i][14] = (float) 1.0f;
										
										if (roiNow.wt > 0) {
											coordMatrix[i][1] = (float) roiNow.z;
											coordMatrix[i][3] = (float) roiNow.x;
											coordMatrix[i][4] = (float) roiNow.y;

											coordMatrix[i][7] = (float) roiNow.w/2 + roiNow.wt;
											coordMatrix[i][8] = (float) roiNow.h/2 + roiNow.wt;

											coordMatrix[i][11] = (float) roiNow.theta;
											coordMatrix[i][13] = (float) 1.0f;
										}
										

									} // if Soil >> Outer circle
									if (roiNow.type == "Wall") { 							// soil: inner circle
										coordMatrix[i][1] = (float) roiNow.z;
										coordMatrix[i][3] = (float) roiNow.x;
										coordMatrix[i][4] = (float) roiNow.y;

										coordMatrix[i][7] = (float) roiNow.w/2;
										coordMatrix[i][8] = (float) roiNow.h/2;

										coordMatrix[i][11] = (float) roiNow.theta;
										coordMatrix[i][13] = (float) 1.0f;
										
										if (roiNow.wt > 0) {
											coordMatrix[i][2] = (float) roiNow.z;				// z (=slice)
											coordMatrix[i][5] = (float) roiNow.x;				// 
											coordMatrix[i][6] = (float) roiNow.y;

											coordMatrix[i][9] = (float) roiNow.w/2 + roiNow.wt;
											coordMatrix[i][10] = (float) roiNow.h/2 + roiNow.wt;

											coordMatrix[i][12] = (float) roiNow.theta;
											coordMatrix[i][14] = (float) 1.0f;
										}
										
									} // if Soil >> Inner circle // if wall >> outer circle
								} // if there is a corresponding ROI entered by user
							} // for submitted ellipses j
						} // for each slice/line i 
						
						rois.sort(Comparator.comparing(RoiCoords::getZ)); 					// sort the submitted rois by Z (slice) top to bottom
						
						for (int i = 0; i < (uniqueZ.length - 1); i ++) {						// for slice with submitted ROI
							int a = uniqueZ[i] - minZ + 1;
							int b = uniqueZ[i+1] - minZ + 1;
							
							for (int j = a + 1; j < b; j++) {
								for (int k = 2; k < 13; k ++) {
									double start [] = {coordMatrix[a - 1][0],    coordMatrix[a - 1][k]};
									double end [] =   {coordMatrix[b - 1][0],      coordMatrix[b - 1][k]};
									coordMatrix[j-1][k] = (float) TailoredMaths.interpolate(start, end, j);
								} // fill matrix (interpolate x,y, diameter, inner, outer , etc.  BUT NOT THETA! THERE DIRECTION DEPENDENT!
								for (int k = 11; k < 13; k ++) {
									double t1 = coordMatrix[a - 1][k]; // theta start
									double t2 = coordMatrix[b - 1][k]; // theta end
									if (Math.abs((t1-t2)) > Math.abs(t1-(t2-Math.PI))){ t2 = t2 - Math.PI; }
									double start [] = {coordMatrix[a - 1][0],    t1};
									double end [] =   {coordMatrix[b - 1][0],    t2};
									coordMatrix[j-1][k] = (float) TailoredMaths.interpolate(start, end, j);
								}
								
							} // loop through slices without defined ROIs 
						} // loop through given z values (input by user)
						
						roiCol.xmid = TailoredMaths.getColumn(coordMatrix, 3);			//x out circle (cylinder wall)
						roiCol.ymid = TailoredMaths.getColumn(coordMatrix, 4);			//y midpoint
						roiCol.zmid = TailoredMaths.getColumn(coordMatrix, 2);
						
						roiCol.ixmid = TailoredMaths.getColumn(coordMatrix, 5);			//x midpoint (inner circle: soil)
						roiCol.iymid = TailoredMaths.getColumn(coordMatrix, 6);			//y midpoint (inner circle)
						
						roiCol.outerMajorRadius = TailoredMaths.getColumn(coordMatrix, 7);
						roiCol.outerMinorRadius = TailoredMaths.getColumn(coordMatrix, 8);
						roiCol.innerMajorRadius = TailoredMaths.getColumn(coordMatrix, 9);
						roiCol.innerMinorRadius = TailoredMaths.getColumn(coordMatrix, 10);
						roiCol.deltaMajorMinorRadius = meanWallThickness(roiCol.innerMajorRadius, roiCol.innerMinorRadius, roiCol.outerMajorRadius, roiCol.outerMinorRadius);
						roiCol.wallThickness = roiCol.deltaMajorMinorRadius;
						roiCol.theta = TailoredMaths.getColumn(coordMatrix, 11);  	//angle of major ellipse axis
						roiCol.itheta = TailoredMaths.getColumn(coordMatrix, 12);  	//angle of major ellipse axis (inner circle)
						
						roiCol.outerR2 = TailoredMaths.getColumn(coordMatrix, 13);
						roiCol.innerR2 = TailoredMaths.getColumn(coordMatrix, 14);
						
						roiCol.heightOfColumn = zValues[zValues.length-1] - zValues[0] + 1;
						roiCol.bottomOfColumn = zValues[zValues.length-1];
						roiCol.topOfColumn = minZ;
						
						// where the magic happens
						jIO.writeInnerCircleVer1(ICpth + File.separator + "Gauge" + filename + ".txt", roiCol);
						

	
						InputOutput jIO = new InputOutput();
						ImageManipulator jIM = new ImageManipulator();
						DisplayThings disp = new DisplayThings();
						
						// load image and cut out vertical ROI... // save a tiff with top and bottom cut off
						mFC.nowTiff = jIO.openTiff3D(mFC.myBaseFolder + System.getProperty("file.separator") + mFC.fileName);
						ImagePlus outTiff = jIM.cutOffTopAndBottomSteel(mFC.nowTiff, roiCol.topOfColumn, roiCol.bottomOfColumn);

						IJ.saveAsTiff(outTiff, mFC.myOutFolder + System.getProperty("file.separator") +  mFC.fileName);
						System.out.println("saved " + mFC.myOutFolder + System.getProperty("file.separator") +  mFC.fileName);
						
												
						// Review Column Outlines
						//disp.displayColumnOutlinesByZ(mFC.nowTiff, roiCol, mFC, null);
						//disp.displayInnerCircleInTiff(mFC.nowTiff, roiCol, mFC, null);
						
						
					} else { // first check whether there is something to export
						String noRoiMessage = "There is no roi to export!";
						JOptionPane.showMessageDialog(canvas, noRoiMessage);
					}
					win.dispose(); // when the inner Circle file is created, close the window!
				} // action event e for Inner Circle button
			});
			innerCirclePanel.add(innerCircleButton, innerCircleConstraints);
			
			// Before you can start "drawing" you should know what you are doing =)
			String infoText = 	"Manual: \n" +
								"For each slice (z) you draw TWO rois: one for the outer cylinder wall, one for the inner.\n\n" +
								"The uppermost two rois correspond to the cylinder top,\n" +
								"the lowermost two rois correspond to the cylinder bottom.\n\n" +
								"After pressing the space bar, the ellipse can be rotated\n" +
								"Careful: then it is not possible to switch the major and minor axis anymore.";
			JOptionPane.showMessageDialog(canvas, infoText);
			
			// Define Panel for displaying (later also deleting) ROI coordinates
			roiCoordinatesPanel.setBorder( BorderFactory.createTitledBorder( "List of ROIs" ) );
			GridBagLayout roiCoordinateLayout = new GridBagLayout();
			GridBagConstraints roiCoordinateConstraints = new GridBagConstraints();
			roiCoordinateConstraints.fill = GridBagConstraints.NONE;
			roiCoordinateConstraints.gridwidth = 1;
			roiCoordinateConstraints.gridheight = 1;
			roiCoordinateConstraints.gridx = 0;
			roiCoordinateConstraints.gridy = 0;
			roiCoordinateConstraints.insets = new Insets(0,0,0,0);
			roiCoordinatesPanel.setLayout( roiCoordinateLayout );	
			
			// field to enter the wall thickness
			JTextField wallThicknessBox = new JTextField(); 
			wallThicknessBox.setName("Wall Thickness (integer):");
			wallThicknessBox.setBounds(50,50,150,20); 
			wallThicknessBox.setInputVerifier(null);
			wallThicknessBox.setEditable(false);
			mainPanel.add(wallThicknessBox);
			
			// If the wall thickness is know, activate the possibility to only draw one ellipse
			JCheckBox wallThickKnown = new JCheckBox("Do you know the wall thickness and only want to draw one ellipse?");
			mainPanel.add(wallThickKnown);
			wallThickKnown.setSelected(false);									// only activate if the checkbox is activated!
			wallThickKnown.addItemListener(new ItemListener() {    
	             public void itemStateChanged(ItemEvent e) {
	            	 boolean wallThickNeeded = wallThicknessBox.isEditable();	// get current state (editable or not?)
	            	 wallThickNeeded = !wallThickNeeded; 						// change it to the opposite
	            	 wallThicknessBox.setEditable(wallThickNeeded);
	             }    
	          });    
			
			
			
			
			JButton addROIButton = new JButton("Addddd ROI");
			mainPanel.add(addROIButton);
			addROIButton.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					
		//**************************************************************************			
				// HERE I SHOULD FIRST ASK WHETHER THIS COMBINATION OF TYPE AND Z ALREADY EXISTS >> THEN DO NOT ACCEPT, FIRST DELETE	
					Roi currentRoi = displayImage.getRoi();
					RoiCoords r = new RoiCoords();
					
					Rectangle bb = currentRoi.getBounds();
					
					r.h = currentRoi.getFeretValues()[0];  //currentRoi.getFloatHeight();
					r.w = currentRoi.getFeretValues()[2];  //   .getFloatWidth();
					r.x = bb.x + bb.width/2.;  //currentRoi.getXBase() + r.w/2;
					r.y = bb.y + bb.height/2.; //currentRoi.getYBase() + r.h/2;
					r.z = displayImage.getCurrentSlice(); //>> mit displayImage funktioniert's logischerweise =)
					r.theta = Math.PI * 2 / 360 * currentRoi.getFeretValues()[1] ; // [0] maj. axis length, [1] angle (degree), [2] min. axis length, [3][4] x,y ouf outer point end of maj. axis
					r.s = rois.size();
					r.wt = 0;						// set wall thickness to 0
					if (wallThickKnown.isSelected()) {
						String enteredWT = wallThicknessBox.getText(); 
						try {
							r.wt = Integer.parseInt(enteredWT); 
						} catch (Exception e1){
							JOptionPane.showMessageDialog(null, "Only integers!");
							return;
						}
					}
					
					if (sampleButton.isSelected()) {r.type = "Soil";}
					else if (wallButton.isSelected()) {r.type = "Wall";}
					else if (calibrationButton.isSelected()) {r.type = "Cali";}
					
					// get the previous combinations of type and z >> should not occur more than once
					for (int q = 0; q < rois.size(); q ++) {
						RoiCoords qRoi = rois.get(q);
						if (r.z == qRoi.z && r.type == qRoi.type) {
							JOptionPane.showMessageDialog(null, "This combination of type (soil/wall) and slice already exists!");
							return;
						}
					}

					rois.add(r); 
					
					int roiNum = rois.size();
					roiCoordinatesPanel.removeAll(); // initialize the buttons again (not very good, better implement a listener and add new buttos on the fly!!!
					if (rois.size() > 0) {
						innerCircleButton.setEnabled(true);
						JButton[] selRoi = new JButton[roiNum];			
						
						for (int i = 0 ; i < roiNum ; i++) {
							if (i < rois.size()) {
								RoiCoords nowRoi = rois.get(i);
								String z = String.format("%1.0f",nowRoi.z);	
								String wtEntered = Integer.toString(nowRoi.wt);
								String coordinate = nowRoi.type + ", slice " + z + ", wall: " + wtEntered;
								selRoi[i] = new JButton(coordinate);								
								selRoi[i].setToolTipText( "click to remove from list" );

								// very unelegant solution. 
								// better not initiate selRoi newly each time a new roi is created and then add action listener to each button itself!
								selRoi[i].addActionListener(new ActionListener() {
									public void actionPerformed(ActionEvent e) {
										
										//loop through allbuttons to check if clicked >> if clicked disable it and when next drawn delete
									    for(int i = 0; i < selRoi.length; i++){
									       if(e.getSource() == selRoi[i]){
												rois.remove(i);
												selRoi[i].setEnabled(false);
									        }
									    }
									}
								}); // actionListener for Roi buttons (to remove them)

								
							} // i < rois.size()	
						} // for roiNum
						for (int i = 0 ; i < rois.size() ; i++) {
							roiCoordinateConstraints.gridy ++;
							roiCoordinatesPanel.add( selRoi[i], roiCoordinateConstraints );
							roiCoordinatesPanel.validate();
							paramsPanel.validate();
						} // for rois.size()
					} else {innerCircleButton.setEnabled(false);}  // if rois.size > 0
				} // actionPerformed (e)

			}); // actionListener addROIButton		
			
			// add components to input image panel   featureselectionpanel (heisst erst mal nur so)
			featureSelectionPanel.setBorder( BorderFactory.createTitledBorder( "Feature Selection" ) );
			GridBagLayout featureSelectionLayout = new GridBagLayout();
			GridBagConstraints featureSelectionConstraints = new GridBagConstraints();
			featureSelectionConstraints.anchor = GridBagConstraints.NORTH;
			featureSelectionConstraints.fill = GridBagConstraints.NONE;
			featureSelectionConstraints.gridwidth = 1;
			featureSelectionConstraints.gridheight = 1;
			featureSelectionConstraints.gridx = 0;
			featureSelectionConstraints.gridy = 0;
			featureSelectionConstraints.insets = new Insets(2, 2, 3, 3); // (insets (abstand) from top left bottom rigth)
			featureSelectionPanel.setLayout( featureSelectionLayout );						

			featureSelectionPanel.add( buttonsAndPicPanel, featureSelectionConstraints );
			featureSelectionConstraints.gridy++;		// gridy 0 specifies top, then goes gradually downwards through grid
			featureSelectionConstraints.anchor = GridBagConstraints.NORTHWEST; // Anchor at top left
			featureSelectionConstraints.fill = GridBagConstraints.HORIZONTAL; 
			featureSelectionPanel.add( wallThickKnown, featureSelectionConstraints );			
			featureSelectionConstraints.gridy++;
			featureSelectionPanel.add( wallThicknessBox, featureSelectionConstraints );			
			featureSelectionConstraints.gridy++;
			featureSelectionPanel.add( addROIButton, featureSelectionConstraints );			
			featureSelectionConstraints.gridy++;
			featureSelectionPanel.add( innerCircleButton, innerCircleConstraints );			
			featureSelectionConstraints.gridy++;


	
  // === final GUI tuning ===
			
			// Parameter panel (left side of the GUI, it includes the 
			// two main panels: feature selection and coordinate list
			GridBagLayout paramsLayout = new GridBagLayout();
			GridBagConstraints paramsConstraints = new GridBagConstraints();
			paramsConstraints.insets = new Insets( 5, 5, 6, 6 );
			paramsPanel.setLayout( paramsLayout );
			paramsConstraints.anchor = GridBagConstraints.NORTHWEST;
			paramsConstraints.fill = GridBagConstraints.HORIZONTAL;
			paramsConstraints.gridwidth = 1;
			paramsConstraints.gridheight = 1;
			paramsConstraints.gridx = 0;
			paramsConstraints.gridy = 0;
			paramsPanel.add( featureSelectionPanel, paramsConstraints);
			paramsConstraints.gridy++;
			paramsConstraints.weighty = 1.0;
			paramsPanel.add( roiCoordinatesPanel, paramsConstraints);
			paramsConstraints.gridy++;

			// main panel (including parameters panel and canvas)
			GridBagLayout layout = new GridBagLayout();
			GridBagConstraints allConstraints = new GridBagConstraints();
			mainPanel.setLayout(layout);

			// put parameter panel in place
			allConstraints.anchor = GridBagConstraints.NORTHWEST;
			allConstraints.fill = GridBagConstraints.BOTH;
			allConstraints.gridwidth = 1;
			allConstraints.gridheight = 1;
			allConstraints.gridx = 0;
			allConstraints.gridy = 0;
			allConstraints.weightx = 0;
			allConstraints.weighty = 0;

			// put canvas in place
			
			mainPanel.add( canvas, allConstraints );

			allConstraints.gridy++;
			allConstraints.weightx = 0;
			allConstraints.weighty = 0;

			// if the input image is 3d, put the
			// slice selectors in place
			if( null != super.sliceSelector ) {
				sliceSelector.setValue( inputImage.getCurrentSlice() );
				displayImage.setSlice( inputImage.getCurrentSlice() );
				new ContrastEnhancer().stretchHistogram(displayImage, 0.5);
				
				mainPanel.add( super.sliceSelector, allConstraints );

				if( null != super.zSelector )
					mainPanel.add( super.zSelector, allConstraints );
				if( null != super.tSelector )
					mainPanel.add( super.tSelector, allConstraints );
				if( null != super.cSelector )
					mainPanel.add( super.cSelector, allConstraints );

			}
			allConstraints.gridy--;

			GridBagLayout wingb = new GridBagLayout();
			GridBagConstraints winc = new GridBagConstraints();
			winc.anchor = GridBagConstraints.NORTHWEST;
			winc.fill = GridBagConstraints.BOTH;
			winc.weightx = 1;
			winc.weighty = 1;
			setLayout( wingb );
			add( mainPanel, winc );
			
			// Fix minimum size to the preferred size at this point
			pack();
			setMinimumSize( getPreferredSize() );

			// add special listener if the input image is a stack
			if(null != sliceSelector) {
				// add adjustment listener to the scroll bar
				sliceSelector.addAdjustmentListener(new AdjustmentListener() 
				{

					public void adjustmentValueChanged(final AdjustmentEvent e) {
						exec.submit(new Runnable() {
							public void run() {							
								if(e.getSource() == sliceSelector)
								{	
									new ContrastEnhancer().stretchHistogram(displayImage, 0.5);
									displayImage.updateAndDraw();	
									Roi nowRoi = displayImage.getRoi();
									if (nowRoi.isVisible()) {   //.isArea
										IJ.run("Fit Ellipse");
									}
								}

							}							
						});
					}
				});

				// mouse wheel listener to update the rois while scrolling
				addMouseWheelListener(new MouseWheelListener() {

					@Override
					public void mouseWheelMoved(final MouseWheelEvent e) {

						exec.submit(new Runnable() {
							public void run() 
							{
								new ContrastEnhancer().stretchHistogram(displayImage, 0.5);
								displayImage.updateAndDraw();
								Roi nowRoi = displayImage.getRoi();
								if (nowRoi.isVisible()) {   //.isArea
									IJ.run("Fit Ellipse");
								}
							}
						});

					}
				});
				
				
				// key listener to repaint the display image and the traces
				// when using the keys to scroll the stack
				KeyListener keyListener = new KeyListener() {

					@Override
					public void keyTyped(KeyEvent e) {}

					@Override
					public void keyReleased(final KeyEvent e) {
						exec.submit(new Runnable() {
							public void run() 
							{
								if(e.getKeyCode() == KeyEvent.VK_LEFT ||
										e.getKeyCode() == KeyEvent.VK_RIGHT ||
										e.getKeyCode() == KeyEvent.VK_LESS ||
										e.getKeyCode() == KeyEvent.VK_GREATER ||
										e.getKeyCode() == KeyEvent.VK_COMMA ||
										e.getKeyCode() == KeyEvent.VK_PERIOD ||
										e.getKeyCode() == KeyEvent.VK_SPACE)
								{
//									new ContrastEnhancer().stretchHistogram(displayImage, 0.5);
									
									// with the Space-Key you can make the ellipse editable
									if (e.getKeyCode() == KeyEvent.VK_SPACE) {
										Roi nowRoi = displayImage.getRoi();
										if (nowRoi.isVisible()) {   //.isArea
											IJ.run("Fit Ellipse");
										}
									}
									displayImage.updateAndDraw();
								}
							}
						});

					}

					@Override
					public void keyPressed(KeyEvent e) {}
				};
				// add key listener to the window and the canvas
				addKeyListener(keyListener);
				canvas.addKeyListener(keyListener);
			

			}
			
			allConstraints.gridx++;
			allConstraints.weightx = 1;
			allConstraints.weighty = 1;
			mainPanel.add( paramsPanel, allConstraints );	
			
//			mainPanel.add(innerCirclePanel);
			
		}// end CustomWindow constructor

		/**
		 * Overwrite windowClosing to display the input image after closing 
		 * the GUI and shut down the executor service
		 */
		@Override
		public void windowClosing( WindowEvent e ) 
		{							
			super.windowClosing( e );
			
			// remove listeners
			sampleButton.removeActionListener( listener );
			wallThickKnown.removeActionListener(listener);
			wallButton.removeActionListener( listener );			
			// addROIButton.removeActionListener( listener );			
			
			if( null != displayImage )
			{
				//displayImage.close();
				displayImage = null;
			}
			// shut down executor service
			exec.shutdownNow();
		}

		
		/**
		 * Update the display image with the gradient or the input image voxels.
		 */
	}
	
	/**
	 * Update the display image with the gradient or the input image voxels.
	 */
	public class RoiCoords {
		public int ID;
		public String type;				//soil / wall
		public double x;				// center 
		public double y;
		public float z;					//slice
		public double h;				// height>> axis
		public double w;				//width >> axis
		public int s;					
		public double theta;			// rot. angle
		public int wt;
		
		public String getType() {
			return type;
		}
		public float getZ() {
			return z;
		}
	} // RoiCoords

	
	
	/**
	 * calculate the mean wall thickness (mean of x and y direction per z-slice)
	 * @param innerMajorRadius 		double[z]
	 * @param innerMinorRadius		double[z]
	 * @param outerMajorRadius		double[z]
	 * @param outerMinorRadius		double[z]
	 * @return mean  Column wall thickness for each slice [double[] ]
	 */
	public static double[] meanWallThickness (double[] innerMajorRadius, double[] innerMinorRadius, double[] outerMajorRadius, double[] outerMinorRadius) {
		double mean [] = new double [innerMajorRadius.length];
		for (int i = 0; i < innerMajorRadius.length; i++) {
			mean[i] = (double) ((outerMajorRadius[i] - innerMajorRadius[i]) + (outerMinorRadius[i] - innerMinorRadius[i]) ) / 2;
		}
		return mean;
	}
	
}
