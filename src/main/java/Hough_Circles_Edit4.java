/** houghCircles_.java:
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 @author Hemerson Pistori (pistori@ec.ucdb.br) and Eduardo Rocha Costa
 @created 18 de Mar�o de 2004
 
 The Hough Transform implementation was based on 
 Mark A. Schulze applet (http://www.markschulze.net/)
 
*/

//package sigus.templateMatching;
//import sigus.*;

import ij.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.*;
import java.awt.*;
import ij.gui.*;
import ij.plugin.*;

/**
*   This ImageJ plugin shows the Hough Transform Space and search for
*   circles in a binary image. The image must have been passed through
*   an edge detection module and have edges marked in white (background
*   must be in black).
*/
public class Hough_Circles_Edit4 implements PlugInFilter {

    public int radiusMin;  // Find circles with radius grater or equal radiusMin
    public int radiusMax;  // Find circles with radius less or equal radiusMax
    public int radiusInc;  // Increment used to go from radiusMin to radiusMax
    public int maxCircles; // Numbers of circles to be found
    public int threshold = -1; // An alternative to maxCircles. All circles with
    // a value in the hough space greater then threshold are marked. Higher thresholds
    // results in fewer circles.
    public int houghIndex; //Contains the current index (slice) position in the Hough radius series
    public int houghIndexCount; //Contains the total number of radii tested
    public double maxPixel; //Contains the brights pixel in the entire Hough array
    public boolean houghSeries; //Contains whether the user wants the Hough series stack as an output
    public boolean showCircles; //Contains whether the user wants the circles found as an output
    public boolean showCentroids; //Contains whether the user wants a map of centroids and radii outputed from search
    public boolean reverseSearch; //Search from largest sphere to smallest sphere
    public boolean recordHough; //Print the Hough Score below the centroid on the centroid map 
    public int radiusBegin; //Sets the starting radius for Hough search
    public int radiusStep; //Sets the increment used during the Hough Search
    byte imageValues[]; // Raw image (returned by ip.getPixels())
    double houghValues[][][]; // Hough Space Values
    public int width; // Hough Space width (depends on image width)
    public int height;  // Hough Space heigh (depends on image height)
    public int depth;  // Hough Space depth (depends on radius interval)
    public int offset; // Image Width
    public int offx;   // ROI x offset
    public int offy;   // ROI y offset
    Point centerPoint[]; // Center Points of the Circles Found.
    int centerRadii[]; //Corresponding radii of the cricles marked by the center points
    int houghScores[]; //Corresponding Hough scores for each centroid
    private final int vectorMaxSize = 500;
    boolean useThreshold = false;
    int lut[][][]; // LookUp Table for rsin e rcos values


    public int setup(String arg, ImagePlus imp) {
        if (arg.equals("about")) {
            showAbout();
            return DONE;
        }
        return DOES_8G+DOES_STACKS+SUPPORTS_MASKING;
    }

    public void run(ImageProcessor ip) {

        imageValues = (byte[])ip.getPixels();
        Rectangle r = ip.getRoi();


        offx = r.x;
        offy = r.y;
        width = r.width;
        height = r.height;
        offset = ip.getWidth();


        if( readParameters() ) { // Show a Dialog Window for user input of
            // radius and maxCircles.
            houghTransform();

            //Calculate the Hough transform image for each radius and input into the stack
            if(houghSeries)
               HoughSpaceSeries();
            
            ImageProcessor circlesip = new ByteProcessor(width, height);
            byte[] circlespixels = (byte[])circlesip.getPixels();

            // Mark the center of the found circles in a new image
            if(useThreshold)
                getCenterPointsByThreshold(threshold);
            else
                getCenterPoints(maxCircles);
            drawCircles(circlespixels);

            // Create image View for Marked Circles.
            if(showCircles){
                new ImagePlus(maxCircles+" Circles Found", circlesip).show();
            }
            
            //Create map of centroids where the intensity is the radius
            if (showCentroids){
                ImageProcessor centroidsip = new ByteProcessor(width, height);
                byte[] centroidpixels = (byte[])centroidsip.getPixels();
                drawCentroids(centroidpixels);
                new ImagePlus("Centroid Map contains " +maxCircles+ " circles", centroidsip).show();
            }
        }
    }
    
    //Clear out the memory after each iteration
        public void run(String arg) {
            for(int a = 1; a < 5; a++)
            System.gc();
        }
        
    void showAbout() {
        IJ.showMessage("About Circles_...",
                       "This plugin finds n circles\n" +
                       "using a basic HoughTransform operator\n." +
                       "For better results apply an Edge Detector\n" +
                       "filter and a binarizer before using this plugin\n"+
                       "\nAuthor: Hemerson Pistori (pistori@ec.ucdb.br)"
                      );
    }

    boolean readParameters() {
        //Get Hough search parameters
        GenericDialog gd = new GenericDialog("Hough Parameters", IJ.getInstance());
        gd.addNumericField("Minimum radius (in pixels) :", 10, 0);
        gd.addNumericField("Maximum radius (in pixels)", 20, 0);
        gd.addNumericField("Increment radius (in pixels) :", 2, 0);
        gd.addNumericField("Number of Circles (NC): (enter 0 if using threshold)", 10, 0);
        gd.addNumericField("Threshold: (not used if NC > 0)", 60, 0);
        
        //Allow user to choose desired output the plugin creates
        gd.addCheckbox("Hough Space Series", true);
        gd.addCheckbox("Circle centroids overlaid on original image", true);
        gd.addCheckbox("Map of centroids where intensity = radius (good for macros)", false);
        gd.addCheckbox("Reverse search direction (search large->small radius)", false);
        gd.addCheckbox("Record the Hough score on the map of centroids (1 pixel right of centroid)", false);

        gd.showDialog();

        if (gd.wasCanceled()) {
            return(false);
        }

        radiusMin = (int) gd.getNextNumber();
        radiusMax = (int) gd.getNextNumber();
        radiusInc = (int) gd.getNextNumber();
        depth = ((radiusMax-radiusMin)/radiusInc)+1;
        maxCircles = (int) gd.getNextNumber();
        threshold = (int) gd.getNextNumber();
        houghSeries = (boolean) gd.getNextBoolean();
        showCircles = (boolean) gd.getNextBoolean();
        showCentroids = (boolean) gd.getNextBoolean();
        reverseSearch = (boolean) gd.getNextBoolean();
        recordHough = (boolean) gd.getNextBoolean();

        
        if (maxCircles > 0) {
            useThreshold = false;
            threshold = -1;
        } else {
            useThreshold = true;
            if(threshold < 0) {
                IJ.showMessage("Threshold must be greater than 0");
                return(false);
            }
        }
        //If the reverse search is not selected then start at min, increment positively
        if(!reverseSearch){
            radiusBegin = radiusMin;
            radiusStep = radiusInc;
        }
        //If reverse search is selected then start at max, increment negatively
        else{
            radiusBegin = radiusMax;
            radiusStep = radiusInc * -1;
        }
        return(true);

    }

    /** The parametric equation for a circle centered at (a,b) with
        radius r is:

    a = x - r*cos(theta)
    b = y - r*sin(theta)

    In order to speed calculations, we first construct a lookup
    table (lut) containing the rcos(theta) and rsin(theta) values, for
    theta varying from 0 to 2*PI with increments equal to
    1/8*r. As of now, a fixed increment is being used for all
    different radius (1/8*radiusMin). This should be corrected in
    the future.

    Return value = Number of angles for each radius
       
    */
    private int buildLookUpTable() {

        int i = 0;
        int incDen = Math.round (8F * radiusMin);  // increment denominator

        lut = new int[2][incDen][depth];
        int radius = radiusBegin;
        while (radius >= radiusMin && radius <= radiusMax ) {
            i = 0;
            for(int incNun = 0; incNun < incDen; incNun++) {
                double angle = (2*Math.PI * (double)incNun) / (double)incDen;
                int indexR = (radius-radiusMin)/radiusInc;
                int rcos = (int)Math.round ((double)radius * Math.cos (angle));
                int rsin = (int)Math.round ((double)radius * Math.sin (angle));
                if((i == 0) | (rcos != lut[0][i][indexR]) & (rsin != lut[1][i][indexR])) {
                    lut[0][i][indexR] = rcos;
                    lut[1][i][indexR] = rsin;
                    i++;
                }
            }
            radius = radius + radiusStep;
        }

        return i;
    }

    private void houghTransform () {

        int lutSize = buildLookUpTable();

        houghValues = new double[width][height][depth];

        int k = width - 1;
        int l = height - 1;

        for(int y = 1; y < l; y++) {
            for(int x = 1; x < k; x++) {
                for(int radius = radiusMin;radius <= radiusMax;radius = radius+radiusInc) {
                    if( imageValues[(x+offx)+(y+offy)*offset] != 0 )  {// Edge pixel found
                        int indexR=(radius-radiusMin)/radiusInc;
                        for(int i = 0; i < lutSize; i++) {

                            int a = x + lut[1][i][indexR]; 
                            int b = y + lut[0][i][indexR]; 
                            if((b >= 0) & (b < height) & (a >= 0) & (a < width)) {
                                houghValues[a][b][indexR] += 1;
                            }
                        }

                    }
                }
            }

        }

    }

    //Find the largest Hough pixel in the 3D Hough transform array to scale the 8-bit conversion
    private double houghMaximum () {
        double d = -1D;
        for(int a=0; a<houghIndexCount; a++){
            for(int j = 0; j < height; j++) {
                for(int k = 0; k < width; k++){
                    if(houghValues[k][j][a] > d) {
                        d = houghValues[k][j][a];
                    }
                }
            }
        }
        return d;
    }
    
    //Create a stack of the full Hough tranform series
    private void HoughSpaceSeries(){
        // Create image View for Hough Transform.
        // Create a Hough Transform image for each radius increment
        //First, create a stack to deposit the transform
        houghIndexCount = (radiusMax-radiusMin)/radiusInc+1;
        ImageStack houghStack = new ImageStack(width, height, houghIndexCount);
            
        //Find the brightest pixel in the whole Hough transform array to normalize stack intensity
        maxPixel = houghMaximum();
        
        int radius = radiusBegin;
        while (radius >= radiusMin && radius <= radiusMax ) {
            //Create a new image for depositing the Hough pixels into
            ImageProcessor newip = new ByteProcessor(width, height);
            byte[] newpixels = (byte[])newip.getPixels();

            //Calculate the corresponding index
            houghIndex = (radius-radiusMin)/radiusInc;
                
            //Retrieve the pixel array for the current Hough radius image
            createHoughPixels(newpixels, houghIndex);
                
            //Deposit the array image into the corresponding slice in the stack
            houghStack.setPixels(newpixels, houghIndex+1);
                
            //Give the current slice the appropriate radius label
            houghStack.setSliceLabel("Hough Space [r="+radius+"]", houghIndex+1);
            
            radius = radius + radiusStep;
            }
        new ImagePlus("Hough Space Series", houghStack).show();
    }
    
    // Convert Values in Hough Space to an 8-Bit Image Space.
    private void createHoughPixels (byte houghPixels[], int index) {
	//Rescale all the Hough values to 8-bit to create the Hough image
        for(int l = 0; l < height; l++) {
            for(int i = 0; i < width; i++) {
                houghPixels[i + l * width] = (byte) Math.round ((houghValues[i][l][index] * 255D) / maxPixel);
            }

        }
    }

    // Draw the circles found in the original image.
    public void drawCircles(byte[] circlespixels) {
		
		// Copy original input pixels into output
		// circle location display image and
		// combine with saturation at 100
		int roiaddr=0;
		for( int y = offy; y < offy+height; y++) {
			for(int x = offx; x < offx+width; x++) {
				// Copy;
				circlespixels[roiaddr] = imageValues[x+offset*y];
				// Saturate
				if(circlespixels[roiaddr] != 0 )
					circlespixels[roiaddr] = 100;
				else
					circlespixels[roiaddr] = 0;
				roiaddr++;
			}
		}
		// Copy original image to the circlespixels image.
		// Changing pixels values to 100, so that the marked
		// circles appears more clear. Must be improved in
		// the future to show the resuls in a colored image.
		//for(int i = 0; i < width*height ;++i ) {
		//if(imageValues[i] != 0 )
		//if(circlespixels[i] != 0 )
		//circlespixels[i] = 100;
		//else
		//circlespixels[i] = 0;
		//}
		if(centerPoint == null) {
			if(useThreshold)
                            getCenterPointsByThreshold(threshold);
			else
                            getCenterPoints(maxCircles);
		}
		byte cor = -1;
		// Redefine these so refer to ROI coordinates exclusively
		int offset = width;
		int offx=0;
		int offy=0;
		
		for(int l = 0; l < maxCircles; l++) {
			int i = centerPoint[l].x;
			int j = centerPoint[l].y;
			// Draw a gray cross marking the center of each circle.
			for( int k = -10 ; k <= 10 ; ++k ) {
				int p = (j+k+offy)*offset + (i+offx);
				if(!outOfBounds(j+k+offy,i+offx))
					circlespixels[(j+k+offy)*offset + (i+offx)] = cor;
				if(!outOfBounds(j+offy,i+k+offx))
					circlespixels[(j+offy)*offset   + (i+k+offx)] = cor;
			}
			for( int k = -2 ; k <= 2 ; ++k ) {
				if(!outOfBounds(j-2+offy,i+k+offx))
					circlespixels[(j-2+offy)*offset + (i+k+offx)] = cor;
				if(!outOfBounds(j+2+offy,i+k+offx))
					circlespixels[(j+2+offy)*offset + (i+k+offx)] = cor;
				if(!outOfBounds(j+k+offy,i-2+offx))
					circlespixels[(j+k+offy)*offset + (i-2+offx)] = cor;
				if(!outOfBounds(j+k+offy,i+2+offx))
					circlespixels[(j+k+offy)*offset + (i+2+offx)] = cor;
			}
		}
	}
    
    // Draw the centroids found in the original image where intensity = radius.
    public void drawCentroids(byte[] centroidpixels) {
		
		for(int l = 0; l < maxCircles; l++) {
			int i = centerPoint[l].x;
			int j = centerPoint[l].y;
                        byte radius = (byte) centerRadii[l];
                        			
                        //Draw a point at the centroid and set the int to = radius
                        centroidpixels[j*offset + i] = radius;
                        
                        //Record the Hough score on the next pixel below the centroid if desired
                        if(recordHough){
                            centroidpixels[j*offset + i + 1] = (byte) houghScores[l];
                        }
                }
	}

    private boolean outOfBounds(int y,int x) {
        if(x >= width)
            return(true);
        if(x <= 0)
            return(true);
        if(y >= height)
            return(true);
        if(y <= 0)
            return(true);
        return(false);
    }

    public Point nthMaxCenter (int i) {
        return centerPoint[i];
    }


    /** Search for a fixed number of circles.

    @param maxCircles The number of circles that should be found.  
    */
    private void getCenterPoints (int maxCircles) {
        houghScores = new int [maxCircles];
        centerRadii = new int [maxCircles];
        centerPoint = new Point[maxCircles];
        int xMax = 0;
        int yMax = 0;
        int rMax = 0;




        for(int c = 0; c < maxCircles; c++) {
            double counterMax = -1;
            int radius = radiusBegin;
            while (radius >= radiusMin && radius <= radiusMax ) {


                int indexR = (radius-radiusMin)/radiusInc;
                for(int y = 0; y < height; y++) {
                    for(int x = 0; x < width; x++) {
                        if(houghValues[x][y][indexR] > counterMax) {
                            counterMax = houghValues[x][y][indexR];
                            xMax = x;
                            yMax = y;
                            rMax = radius;
                        }
                    }

                }
                radius = radius + radiusStep;
            }

            centerPoint[c] = new Point (xMax, yMax);
            centerRadii[c] = rMax;
            houghScores[c] = (int) counterMax;
            

            clearNeighbours(xMax,yMax,rMax);
        }
    }


    /** Search circles having values in the hough space higher than a threshold

    @param threshold The threshold used to select the higher point of Hough Space
    */
    private void getCenterPointsByThreshold (int threshold) {
        
        houghScores = new int [vectorMaxSize];
        centerRadii = new int [vectorMaxSize];
        centerPoint = new Point[vectorMaxSize];
        int xMax = 0;
        int yMax = 0;
        int countCircles = 0;

        int radius = radiusBegin;
        while (radius >= radiusMin && radius <= radiusMax ) {
            int indexR = (radius-radiusMin)/radiusInc;
            for(int y = 0; y < height; y++) {
                for(int x = 0; x < width; x++) {



                    if(houghValues[x][y][indexR] > threshold) {


                        if(countCircles < vectorMaxSize) {


                            centerPoint[countCircles] = new Point (x, y);
                            centerRadii[countCircles] = radius;
                            houghScores[countCircles] = (int) houghValues[x][y][indexR];

                            clearNeighbours(x,y,radius);

                            ++countCircles;
                        } else
                            break;
                    }
                }
            }
            radius = radius + radiusStep;
        }

        maxCircles = countCircles;
    }

    /** Clear, from the Hough Space, all the counter that are near (radius/2) a previously found circle C.
        
    @param x The x coordinate of the circle C found.
    @param x The y coordinate of the circle C found.
    @param x The radius of the circle C found.
    */
    private void clearNeighbours(int x,int y, int radius) {


        // The following code just clean the points around the center of the circle found.
	double radiusSquared = radius*radius;


        int y1 = (int)Math.floor ((double)y - radius);
        int y2 = (int)Math.ceil ((double)y + radius) + 1;
        int x1 = (int)Math.floor ((double)x - radius);
        int x2 = (int)Math.ceil ((double)x + radius) + 1;



        if(y1 < 0)
            y1 = 0;
        if(y2 > height)
            y2 = height;
        if(x1 < 0)
            x1 = 0;
        if(x2 > width)
            x2 = width;



        int r = radius;
        while (r >= radiusMin && r <= radiusMax )  {
            int indexR = (r-radiusMin)/radiusInc;
            for(int i = y1; i < y2; i++) {
                for(int j = x1; j < x2; j++) {	      	     
                    if(Math.pow (j - x, 2D) + Math.pow (i - y, 2D) < radiusSquared) {
                        houghValues[j][i][indexR] = 0.0D;
                    }
                }
            }
            r = r + radiusStep;
        }

    }

}

