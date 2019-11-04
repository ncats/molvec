package gov.nih.ncats.molvec.internal.image.binarization;

import java.awt.image.Raster;
import java.util.function.Consumer;

import gov.nih.ncats.molvec.internal.image.Bitmap;

/*
 * Adaptive Threshold based on local pixel average Loosely based on Bradley
 * & Roth Integral Image Adaptive Thresholding Added ability to use a
 * specific SIGMA minimum threshold. General Algorithm: 
 * 	1)Precompute the integral image (used to get local mean) 
 * 	2)Precompute the square integral image (used to get local 
 * 		standard dev) 
 * 	3)For each pixel: 
 * 		a)get the sum and sum of squares for those pixels in 
 * 		  an (wsize*2+1) sized box around the given pixel. 
 * 		  (using the precomputed integrals) 
 * 		b)calculate mean =  sum/ pixel count 
 * 		c)calculate stDEV = sqrt(sumSquares/(pixel count)-* mean^2) 
 * 		d)set threshold for given pixel at mean+stDEV*SIGMA 
 * 		  (where sigma is a provided constant)
 * UPDATE: Integral images calculated only for rows needed, on the fly
 */
public class AdaptiveThreshold implements Binarization {
    public static final double DEFAULT_SIGMA_THRESHOLD = 3;
    public static final int DEFAULT_ADAPTIVE_BOX_RADIUS = 200;
    public static final double DEFAULT_ADAPTIVE_MIN_THRESHOLD_RATIO = 0.2;
    public static final double DEFAULT_ADAPTIVE_MAX_THRESHOLD_RATIO = 0.8;
    public static final int DEFAULT_ADAPTIVE_MIN_STDDEV = 0;

    private int wsize;
    
    private double absMax;
    
    private double absMin, sigma, minSigma;

    public AdaptiveThreshold () {
        this (DEFAULT_ADAPTIVE_BOX_RADIUS);
    }

    public AdaptiveThreshold (int wsize) {
        this (wsize, DEFAULT_SIGMA_THRESHOLD, 
        		DEFAULT_ADAPTIVE_MIN_THRESHOLD_RATIO,
              DEFAULT_ADAPTIVE_MIN_STDDEV);
    }

    public AdaptiveThreshold (int wsize, double sigma, 
                              double absMin, double minSigma) {
        this.wsize = wsize;
        this.sigma = sigma;
        this.absMin = absMin;
        this.minSigma = minSigma;
        this.absMax = DEFAULT_ADAPTIVE_MAX_THRESHOLD_RATIO;
        
    }

   

    private static void addIntLines (Raster inRaster, int[] yMap,
                                     double[][] intLines,
                                     double[][] intSquareLines,
                                     int width, int y) {
        double sum = 0;
        double sumSquare = 0;
        double[] intLine = intLines[yMap[y]];
        double[] intSquareLine = intSquareLines[yMap[y]];
        
        double[] intPrevLine = null;
        double[] intSquarePrevLine = null;
        if (y > 0) {
            intPrevLine = intLines[yMap[y - 1]];
            intSquarePrevLine = intSquareLines[yMap[y - 1]];
        }

        for (int x = 0; x < width; ++x) {
            double sampleDouble = inRaster.getSampleDouble(x, y, 0);
            sum += sampleDouble;
            sumSquare += sampleDouble * sampleDouble;
            if (y == 0) {
                intLine[x] = sum;
                intSquareLine[x] = sumSquare;
            } else {
                intLine[x] = sum + intPrevLine[x];
                intSquareLine[x] = sumSquare + intSquarePrevLine[x];
            }
        }
    }

	@Override
	public Bitmap binarize(Raster inRaster, ImageStats stats, Consumer<ImageStats> cons) {

		if(stats==null){
		    stats = Binarization.computeImageStats(inRaster);
        }
		
		
        Bitmap bm = new Bitmap (inRaster.getWidth (), inRaster.getHeight ());
        wsize = (Math.min(wsize * 2 + 2,bm.height())-2)/2;
        wsize = (Math.min(wsize * 2 + 2,bm.width())-2)/2;
        
        int[] yMap = new int[bm.height()];
        double[][] intLines = new double[wsize * 2 + 2][bm.width()];
        double[][] intSquareLines = new double[wsize * 2 + 2][bm.width()];
        

        // make Integral-Image for quick averaging
        for (int y = 0; y < wsize * 2 + 2; ++y) {
            yMap[y] = y;
            addIntLines (inRaster, yMap, intLines, 
                         intSquareLines, bm.width(), y);
        }

        
        double range = stats.max-stats.min;
        
        
        double[] dl = new double[bm.width()];
        for (int y = 0; y < bm.height(); ++y) {
            int y1 = Math.max (y - wsize, 0);
            int y2 = Math.min (y + wsize, bm.height() - 1);
            int boxHeight = (y2 - y1 + 1);
            double[] topIntLine = null;
            double[] topSquareIntLine = null;
            double[] bottomIntLine;
            double[] bottomSquareIntLine;

            if (y1 >= 2) {
                if (boxHeight >= wsize * 2 + 1) {
                    yMap[y2] = yMap[y1 - 2];
                    addIntLines (inRaster, yMap, intLines, intSquareLines,
                                 bm.width(), y2);
                }
            }
            if (y1 > 0) {
                topIntLine = intLines[yMap[y1 - 1]];
                topSquareIntLine = intSquareLines[yMap[y1 - 1]];
            }
            bottomIntLine = intLines[yMap[y2]];
            bottomSquareIntLine = intSquareLines[yMap[y2]];
            inRaster.getSamples(0, y, bm.width(), 1, 0, dl);
            for (int x = 0; x < bm.width(); ++x) {
                // box to average over:
                int x1 = Math.max (x - wsize, 0);
                int x2 = Math.min (x + wsize, bm.width() - 1);
                int count = (x2 - x1 + 1) * boxHeight;
                double sum = bottomIntLine[x2];
                double sumSquare = bottomSquareIntLine[x2];
                if (y1 > 0) {
                    sum += -topIntLine[x2];
                    sumSquare += -topSquareIntLine[x2];
                }
                if (x1 > 0) {
                    sum += -bottomIntLine[x1 - 1];
                    sumSquare += -bottomSquareIntLine[x1 - 1];
                }
                if (y1 > 0 && x1 > 0) {
                    sum += topIntLine[x1 - 1];
                    sumSquare += topSquareIntLine[x1 - 1];
                }

                double mean = sum / count;
                double stdDEV = Math.sqrt 
                    (Math.abs (sumSquare / count- mean * mean));
                double threshold = Math.min (mean + stdDEV * sigma, stats.min + range*absMax);
                threshold = Math.max(threshold, stats.min+range*absMin);
                double pel = dl[x];
                bm.set (x, y, pel > threshold);
            }
        }

        return bm;
	}
}
