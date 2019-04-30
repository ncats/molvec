package gov.nih.ncats.molvec.image.binarization;

import java.awt.image.Raster;
import java.util.logging.Logger;

import gov.nih.ncats.molvec.Bitmap;
import gov.nih.ncats.molvec.image.Binarization;

/*
 * Simple threshold based on full image mean and standard deviation. 
 */
public class SigmaThreshold implements Binarization {

    static final Logger logger = Logger.getLogger(SigmaThreshold.class.getName());
    public static double DEFAULT_SIGMA_THRESHOLD = 1.2;
    private double sigma;
    
    public SigmaThreshold () {
        this (DEFAULT_SIGMA_THRESHOLD);
    }

    public SigmaThreshold (double sigma) {
        this.sigma = sigma;
    }

    public void setSigma (double sigma) { this.sigma = sigma; }
    public double getSigma () { return sigma; }

    public Bitmap binarize (Raster inRaster) {
        Bitmap bm = new Bitmap (inRaster.getWidth (), inRaster.getHeight ());
        double max = -1, min = Double.MAX_VALUE;
        double sum = 0;
        
        double sumSquare = 0;
        for (int y = 0; y < bm.height(); ++y) {
            for (int x = 0; x < bm.width(); ++x) {
                double pel = inRaster.getSampleDouble (x, y, 0);
                sum += pel;
                sumSquare += pel * pel;
                if (pel < min)
                    min = pel;
                if (pel > max)
                    max = pel;
            }
        }
        long tot = bm.height() * bm.width();
        double mean = sum / tot;
        double stdDEV = Math.sqrt (sumSquare / tot - mean * mean);
        
        double countAbove=0;
        double countBelow=0;
        
        double sumTop = 0;
        double sumBottom = 0;
        double sumSquareBottom = 0;
        double sumSquareTop = 0;
        
        
        for (int y = 0; y < bm.height(); ++y) {
            for (int x = 0; x < bm.width(); ++x) {
                double pel = inRaster.getSampleDouble (x, y, 0);
                if(pel>mean){
                	countAbove++;
                	sumTop+=pel;
                	sumSquareTop+=pel*pel;
                }else if(pel<mean){
                	countBelow++;
                	sumBottom+=pel;
                	sumSquareBottom+=pel*pel;
                }
                        
            }
        }
        double meanTop = sumTop / (double)countAbove;
        double stdDEVTop = Math.sqrt (sumSquareTop / countAbove - meanTop * meanTop);
        
        double meanBottom = sumBottom / (double)countBelow;
        double stdDEVBottom = Math.sqrt (sumSquareBottom / countBelow - meanBottom * meanBottom);
        
        double threshold=mean;
        double sign = 1;
        
        
        threshold = mean + stdDEV * sigma;
        
        if(threshold>max || threshold<min){
	        if(countAbove>countBelow){
	        	//System.out.println("More high than low");
	        	threshold = Math.min(meanBottom + stdDEV * sigma,max);
	        	
	        }else if(countAbove<countBelow){
	        	//System.out.println("More low than high");
	        	threshold = Math.max(meanTop - stdDEV * sigma, min);
	        }	
        }
        
        //double threshold =mean + stdDEV * sigma * sign;
        logger.info("## for threshold, mean ="+mean+", stDEV="+stdDEV + ", threshold=" + threshold + ", max=" + max + ", min=" + min);
        logger.info("## for threshold, countAbove ="+countAbove + ", countBelow="+countBelow);
        for (int y = 0; y < bm.height(); ++y) {
            for (int x = 0; x < bm.width(); ++x) {
                double pel = inRaster.getSampleDouble (x, y, 0);
                bm.set (x, y, pel >= threshold);
            }
        }

        return bm;
    }
}
