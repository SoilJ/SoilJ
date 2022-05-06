package SoilJ_;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.Panel;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.SwingUtilities;

import SoilJ.tools.DisplayThings;
import SoilJ.tools.ImageManipulator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.ObjectDetector;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.ImageCanvas;
import ij.gui.StackWindow;
import ij.gui.Toolbar;
import ij.plugin.PlugIn;

/**
 * DefineColumnOutlinesManually is a SoilJ plugin that is supposed to one fine
 * day contain a manual column outline detection.
 * 
 * @author John Koestel
 *
 */

public class DefineColumnOutlinesManually_ extends ImagePlus implements PlugIn {

	@Override
	public void run(String arg) {

					
		// probe operating system and adjust PathSep if necessary
		String pathSep;
		String myOS = System.getProperty("os.name");
		if (myOS.equalsIgnoreCase("Linux"))
			pathSep = "/";

		// construct biggish objects
		InputOutput jIO = new InputOutput();
		ObjectDetector jOD = new ObjectDetector();
		MenuWaiter menu = new MenuWaiter();
		DisplayThings disp = new DisplayThings();
		ImageManipulator jIM = new ImageManipulator();

		// init variables
		int i;

		// read file or files
		InputOutput.MyFileCollection mFC = jIO.fileSelector("Please choose a file or folder with your image data");

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
		int errorCounts = 0;
		for (i = 0; i < mFC.myTiffs.length; i++) {

			// assign tiff file
			mFC.fileName = mFC.myTiffs[i];
			mFC = jIO.addCurrentFileInfo(mFC);

			// load image
			ImagePlus nowTiff = jIO.openTiff3D(mFC.myBaseFolder + mFC.pathSep + mFC.myTiffs[i]);

			inputImage = nowTiff;			
					
			inputStackCopy = inputImage.getImageStack().duplicate();
			displayImage = new ImagePlus( inputImage.getTitle(), 
					inputStackCopy );
			displayImage.setTitle("Manual Feature Definition");
			displayImage.setSlice( inputImage.getCurrentSlice() );

			// hide input image (to avoid accidental closing)
			//inputImage.getWindow().setVisible( false );

			// set the 2D flag
			inputIs2D = inputImage.getImageStackSize() == 1;

			// correct Fiji error when the slices are read as frames
			if ( inputIs2D == false && 
					displayImage.isHyperStack() == false && 
					displayImage.getNSlices() == 1 )
			{
				// correct stack by setting number of frames as slices
				displayImage.setDimensions( 1, displayImage.getNFrames(), 1 );
			}
			
			// Build GUI
			SwingUtilities.invokeLater(
					new Runnable() {
						public void run() {
							win = new CustomWindow( displayImage );
							win.pack();
						}
					});
			
		}
	}
	
	/** main GUI window */
	private CustomWindow win;

	/** original input image */
	ImagePlus inputImage = null;
	
	ImageStack inputStackCopy = null;

	/** image to be displayed in the GUI */
	ImagePlus displayImage = null;

	/** gradient image stack */
	ImageStack gradientStack = null;

	/** image containing the final results of the watershed segmentation (basins with or without dams) */
	ImagePlus resultImage = null;		

	/** parameters panel (segmentation + display options) */
	JPanel paramsPanel = new JPanel();

	/** main panel */
	Panel all = new Panel();
	
	/** flag to indicate 2D input image */
	boolean inputIs2D = false;

	/** feature selection panel */
	JPanel featureSelectionPanel = new JPanel();
	/** input image options */
	ButtonGroup inputImageButtons;
	/** radio button to specify that soil should be selected */
	JRadioButton sampleButton;
	/** radio button to specify that wall should be selected */
	JRadioButton wallButton;
	/** radio button to specify that wall should be selected */
	JRadioButton calibrationButton;
	/** text for soil radio button */
	static String sampleImageText = "Soil";
	/** text for wall radio button */
	static String wallImageText = "Wall";
	/** text for calibration rod radio button */
	static String calibrationImageText = "Calibration";
	
	/** panel to store the radio buttons with the image type options */
	JPanel radioPanel = new JPanel( new GridLayout(0, 1) );
	JLabel inputImagePicture;
	JPanel buttonsAndPicPanel = new JPanel();
	
	/** Button to add ROIs */
	JButton addROIButton;
	/** text of addROI button */
	private String addROIText = "add ROI";
	/** tip text of addROI button */
	private String addROITip = "add ROI to InnerCircleFile";

	/** ROI coordinates panel */
	JPanel roiCoordinatesPanel = new JPanel();
	/** init struct to remember selected ROIs   */
	ArrayList<RoiCoords> rois = new ArrayList<RoiCoords>();
	
	/** executor service to launch threads for the plugin methods and events */
	final ExecutorService exec = Executors.newFixedThreadPool(1);
	/** opacity to display overlays */
	double opacity = 1.0/3.0;
	
	/** == > go create !! */
	
	private class CustomWindow extends StackWindow
	{
		/**
		 * serial version uid
		 */
		private static final long serialVersionUID = -6201855439028581892L;

		/**
		 * Listener for the GUI buttons
		 */
		private ActionListener listener = new ActionListener() {

			@Override
			public void actionPerformed( final ActionEvent e ) 
			{

				final String command = e.getActionCommand();

				// listen to the buttons on separate threads not to block
				// the event dispatch thread
				exec.submit(new Runnable() {

					public void run()
					{
						// "add ROI"  button		
						if( e.getSource() == addROIButton )
						{
							;						
						}						
//						else if( e.getSource()  == shuffleColorsButton )
//						{
//							shuffleColors();
//						}
					}

					
					
				});
			}

		};


		/**
		 * Construct the plugin window
		 * 
		 * @param imp input image
		 */
		CustomWindow( ImagePlus imp )
		{
			super(imp, new ImageCanvas(imp));

			final ImageCanvas canvas = (ImageCanvas) getCanvas();

			Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
			double screenWidth = screenSize.getWidth();
			double screenHeight = screenSize.getHeight();
			
			// Zoom in if image is too small
			while( ( ic.getWidth() < screenWidth/2 ||
				ic.getHeight() < screenHeight/2 ) &&
				ic.getMagnification() < 32.0 )			
			    {
    				final int canvasWidth = ic.getWidth();
    				ic.zoomIn( 0, 0 );
    				// check if canvas size changed (otherwise stop zooming)
    				if( canvasWidth == ic.getWidth() )
    				{
    					ic.zoomOut(0, 0);
    					break;
    				}
			    }
			// Zoom out if canvas is too large
			while( ( ic.getWidth() > 0.75 * screenWidth ||
				ic.getHeight() > 0.75 * screenHeight ) &&
				ic.getMagnification() > 1/72.0 )
				{
	    			final int canvasWidth = ic.getWidth();
	    			ic.zoomOut( 0, 0 );
	    			// check if canvas size changed (otherwise stop zooming)
	    			if( canvasWidth == ic.getWidth() )
	    			{
	    				ic.zoomIn(0, 0);
	    				break;
	    			}
				}

			setTitle( "Feature Location Definition" );

			// select oval ROI tool to define outlines of features of interest
			Toolbar.getInstance().setTool( Toolbar.OVAL );

			// === feature selection panel ===

			// feature selection options (soil sample, inner wall perimeter or calibration rod)			
			sampleButton = new JRadioButton( sampleImageText );			
			sampleButton.setActionCommand( sampleImageText );
			sampleButton.addActionListener( listener );
			sampleButton.setToolTipText( "outline of soil sample" );
						
			wallButton = new JRadioButton( wallImageText );
			wallButton.setActionCommand( wallImageText );
			wallButton.addActionListener( listener );
			wallButton.setToolTipText( "location of inner wall perimeter" );
			
			calibrationButton = new JRadioButton( calibrationImageText );
			calibrationButton.setActionCommand( calibrationImageText );
			calibrationButton.addActionListener( listener );
			calibrationButton.setToolTipText( "location of calibration rod" );

			inputImageButtons = new ButtonGroup();
			inputImageButtons.add( sampleButton );
			inputImageButtons.add( wallButton );
			inputImageButtons.add( calibrationButton );
			radioPanel.add( sampleButton );
			radioPanel.add( wallButton );
			radioPanel.add( calibrationButton );

			buttonsAndPicPanel.add( radioPanel);
			
			// add ROI button
			addROIButton = new JButton( addROIText );
			addROIButton.setToolTipText( addROITip );
			addROIButton.addActionListener( listener );
		
			// add components to input image panel
			featureSelectionPanel.setBorder( BorderFactory.createTitledBorder( "Feature Selection" ) );
			GridBagLayout featureSelectionLayout = new GridBagLayout();
			GridBagConstraints featureSelectionConstraints = new GridBagConstraints();
			featureSelectionConstraints.anchor = GridBagConstraints.NORTH;
			featureSelectionConstraints.fill = GridBagConstraints.NONE;
			featureSelectionConstraints.gridwidth = 1;
			featureSelectionConstraints.gridheight = 1;
			featureSelectionConstraints.gridx = 0;
			featureSelectionConstraints.gridy = 0;
			featureSelectionConstraints.insets = new Insets(2, 2, 3, 3);
			featureSelectionPanel.setLayout( featureSelectionLayout );						

			featureSelectionPanel.add( buttonsAndPicPanel, featureSelectionConstraints );
			featureSelectionConstraints.gridy++;
			featureSelectionConstraints.anchor = GridBagConstraints.NORTHWEST;
			featureSelectionConstraints.fill = GridBagConstraints.HORIZONTAL;
			featureSelectionPanel.add( addROIButton, featureSelectionConstraints );			
			featureSelectionConstraints.gridy++;

			// === selected coordinate list panel ===
			
			RoiCoords testRoi = new RoiCoords();
			testRoi.type = "gaga";
			testRoi.x=1;
			testRoi.y=1;
			testRoi.z=1;
		
			rois.add(testRoi);
			rois.add(testRoi);
			rois.add(testRoi);
			rois.add(testRoi);
			int roiNum = rois.size();
			
			JButton[] selRoi = new JButton[roiNum];			
			
			for (int i = 0 ; i < roiNum ; i++) {
				
				if (i < rois.size()) {
					
					RoiCoords nowRoi = rois.get(i);
					
					String x = String.format("%1.2f",nowRoi.x);
					String y = String.format("%1.2f",nowRoi.y);
					String z = String.format("%1.2f",nowRoi.z);
					String coordinate = nowRoi.type + ": x=" + x + " y=" + y + " z=" + z;					
					selRoi[i] = new JButton(coordinate);								
					selRoi[i].setToolTipText( "click to remove from list" );
					
				}				
			}
			
			roiCoordinatesPanel.setBorder( BorderFactory.createTitledBorder( "List of ROIs" ) );
			GridBagLayout roiCoordinateLayout = new GridBagLayout();
			GridBagConstraints roiCoordinateConstraints = new GridBagConstraints();
			roiCoordinateConstraints.anchor = GridBagConstraints.NORTHWEST;
			roiCoordinateConstraints.fill = GridBagConstraints.NONE;
			roiCoordinateConstraints.gridwidth = 1;
			roiCoordinateConstraints.gridheight = 1;
			roiCoordinateConstraints.gridx = 0;
			roiCoordinateConstraints.gridy = 0;
			roiCoordinateConstraints.insets = new Insets(2, 2, 3, 3);
			roiCoordinatesPanel.setLayout( roiCoordinateLayout );				
			
			for (int i = 0 ; i < roiNum ; i++) {
			
				roiCoordinatesPanel.add( selRoi[i], roiCoordinateConstraints );
				roiCoordinateConstraints.gridy++;
				
			}
			
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
			paramsPanel.add( roiCoordinatesPanel, paramsConstraints);
			paramsConstraints.gridy++;
//			paramsPanel.add( displayPanel, paramsConstraints);
//			paramsConstraints.gridy++;
//			paramsPanel.add( postProcessPanel, paramsConstraints);
//			paramsConstraints.gridy++;


			// main panel (including parameters panel and canvas)
			GridBagLayout layout = new GridBagLayout();
			GridBagConstraints allConstraints = new GridBagConstraints();
			all.setLayout(layout);

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
			
			all.add( canvas, allConstraints );

			allConstraints.gridy++;
			allConstraints.weightx = 0;
			allConstraints.weighty = 0;

			// if the input image is 3d, put the
			// slice selectors in place
			if( null != super.sliceSelector )
			{
				sliceSelector.setValue( inputImage.getCurrentSlice() );
				displayImage.setSlice( inputImage.getCurrentSlice() );
				
				all.add( super.sliceSelector, allConstraints );

				if( null != super.zSelector )
					all.add( super.zSelector, allConstraints );
				if( null != super.tSelector )
					all.add( super.tSelector, allConstraints );
				if( null != super.cSelector )
					all.add( super.cSelector, allConstraints );

			}
			allConstraints.gridy--;

			GridBagLayout wingb = new GridBagLayout();
			GridBagConstraints winc = new GridBagConstraints();
			winc.anchor = GridBagConstraints.NORTHWEST;
			winc.fill = GridBagConstraints.BOTH;
			winc.weightx = 1;
			winc.weighty = 1;
			setLayout( wingb );
			add( all, winc );
			
			// Fix minimum size to the preferred size at this point
			pack();
			setMinimumSize( getPreferredSize() );

			// add especial listener if the input image is a stack
			if(null != sliceSelector)
			{
				// add adjustment listener to the scroll bar
				sliceSelector.addAdjustmentListener(new AdjustmentListener() 
				{

					public void adjustmentValueChanged(final AdjustmentEvent e) {
						exec.submit(new Runnable() {
							public void run() {							
								if(e.getSource() == sliceSelector)
								{		
									displayImage.updateAndDraw();								
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
								displayImage.updateAndDraw();
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
										e.getKeyCode() == KeyEvent.VK_PERIOD)
								{
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
			all.add( paramsPanel, allConstraints );			
			
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
			wallButton.removeActionListener( listener );			
			addROIButton.removeActionListener( listener );			
			
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
		void updateDisplayImage()
		{
			int slice = displayImage.getCurrentSlice();
			displayImage.setStack( inputStackCopy );
			displayImage.setSlice( slice );
			displayImage.updateAndDraw();
		}		

	}
	
	
	public class RoiCoords {
		
		public String type;
		public float x;
		public float y;
		public float z;
	
	}
}
