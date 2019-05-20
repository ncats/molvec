package gov.nih.ncats.molvec.image.binarization;

import java.awt.image.Raster;
import java.util.function.Consumer;
import java.util.logging.Logger;

import gov.nih.ncats.molvec.Bitmap;
import gov.nih.ncats.molvec.image.Binarization;

/*
 * Simple threshold based on full image mean and standard deviation. 
 */
public class SigmaThreshold implements Binarization {

    static final Logger logger = Logger.getLogger(SigmaThreshold.class.getName());
    public static double DEFAULT_SIGMA_THRESHOLD = 1.2;

    public static final double DEFAULT_MIN_THRESHOLD_RATIO = 0.1;
    public static final double DEFAULT_MAX_THRESHOLD_RATIO = 0.9;
    
    private double sigma;
    
    public SigmaThreshold () {
        this (DEFAULT_SIGMA_THRESHOLD);
    }

    public SigmaThreshold (double sigma) {
        this.sigma = sigma;
    }

    public void setSigma (double sigma) { this.sigma = sigma; }
    public double getSigma () { return sigma; }

    
    public static class ImageStats{
    	public double min;
    	public double max;
    	public double stdev;
    	public double mean;
    	public double threshold;
    	public int[] histogram;
    	public double count;
    	
    }
    public Bitmap binarize (Raster inRaster) {
    	return binarize(inRaster,(is)->{});
    }
    public Bitmap binarize (Raster inRaster, Consumer<ImageStats> cons) {
        Bitmap bm = new Bitmap (inRaster.getWidth (), inRaster.getHeight ());
        
        ImageStats stats = new ImageStats();
        
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
        
        
        stats.histogram = new int[101];
        for (int y = 0; y < bm.height(); ++y) {
            for (int x = 0; x < bm.width(); ++x) {
                double pel = inRaster.getSampleDouble (x, y, 0);
                double pct= (100.0*(pel-min))/(max-min);
                stats.histogram[(int)pct]++;
                if(pel>mean){
                	countAbove++;
                	sumTop+=pel;
                }else if(pel<mean){
                	countBelow++;
                	sumBottom+=pel;
                }       
            }
        }
        double meanTop = sumTop / (double)countAbove;
        double meanBottom = sumBottom / (double)countBelow;
        double threshold=mean;
        
        
        threshold = mean + stdDEV * sigma;
        
        if(threshold>max || threshold<min){
        	//determine whether inverted
	        if(countAbove>countBelow){
	        	threshold = Math.min(meanBottom + stdDEV * sigma,max);
	        }else if(countAbove<countBelow){
	        	threshold = Math.max(meanTop - stdDEV * sigma, min);
	        }	
        }
        
        threshold = Math.max(threshold,min+(max-min)*DEFAULT_MIN_THRESHOLD_RATIO);
        threshold = Math.min(threshold,min+(max-min)*DEFAULT_MAX_THRESHOLD_RATIO);
        
        for (int y = 0; y < bm.height(); ++y) {
            for (int x = 0; x < bm.width(); ++x) {
                double pel = inRaster.getSampleDouble (x, y, 0);
                bm.set (x, y, pel >= threshold);
            }
        }
        stats.min=min;
        stats.max=max;
        stats.mean=mean;
        stats.stdev=stdDEV;
        stats.threshold=threshold;
        stats.count=bm.width()*bm.height();
        cons.accept(stats);
        
        return bm;
    }
}
