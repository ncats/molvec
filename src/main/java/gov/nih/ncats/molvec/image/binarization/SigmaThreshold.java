package gov.nih.ncats.molvec.image.binarization;

import java.awt.image.Raster;
import java.util.function.Consumer;
import java.util.logging.Logger;

import gov.nih.ncats.molvec.image.Bitmap;

/*
 * Simple threshold based on full image mean and standard deviation. 
 */
public class SigmaThreshold implements Binarization {

    static final Logger logger = Logger.getLogger(SigmaThreshold.class.getName());
    public static double DEFAULT_SIGMA_THRESHOLD = 1.2;
    
    

    public static final double DEFAULT_MIN_THRESHOLD_RATIO = 0.1;
    public static final double DEFAULT_MAX_THRESHOLD_RATIO = 0.9;
    
    private double minThresholdRation = DEFAULT_MIN_THRESHOLD_RATIO;
    private double maxThresholdRation = DEFAULT_MAX_THRESHOLD_RATIO;
    
    private double sigma;
    
    public SigmaThreshold () {
        this (DEFAULT_SIGMA_THRESHOLD);
    }

    public SigmaThreshold (double sigma) {
        this.sigma = sigma;
    }
    public SigmaThreshold (double sigma, double minRat, double maxRat) {
        this.sigma = sigma;
        this.minThresholdRation=minRat;
        this.maxThresholdRation=maxRat;
    }
    

    public void setSigma (double sigma) { this.sigma = sigma; }
    public double getSigma () { return sigma; }

    
	@Override
	public Bitmap binarize(Raster inRaster, ImageStats stats, Consumer<ImageStats> cons) {
		 Bitmap bm = new Bitmap (inRaster.getWidth (), inRaster.getHeight ());
	        
	        if(stats==null)stats = Binarization.computeImageStats(inRaster);
	        
	        double countAbove=0;
	        double countBelow=0;
	        
	        double sumTop = 0;
	        double sumBottom = 0;
	        
	        for(int i=0;i<stats.histogramRaw.length;i++){
	        	if(i>stats.mean){
	        		countAbove+=stats.histogramRaw[i];
	        		sumTop+=i*stats.histogramRaw[i];
	        	}else if(i<stats.mean){
	        		countBelow+=stats.histogramRaw[i];
	        		sumBottom+=i*stats.histogramRaw[i];
	        	}
	        	
	        }
	        
	        
	        
	        
	        double meanTop = sumTop / (double)countAbove;
	        double meanBottom = sumBottom / (double)countBelow;
	        double threshold=stats.mean;
	        
	        
	        threshold = stats.mean + stats.stdev * sigma;
	        
	        if(threshold>stats.max || threshold<stats.min){
	        	//determine whether inverted
		        if(countAbove>countBelow){
		        	threshold = Math.min(meanBottom + stats.stdev * sigma,stats.max);
		        }else if(countAbove<countBelow){
		        	threshold = Math.max(meanTop - stats.stdev * sigma, stats.min);
		        }	
	        }
	        
	        threshold = Math.max(threshold,stats.min+(stats.max-stats.min)*minThresholdRation);
	        threshold = Math.min(threshold,stats.min+(stats.max-stats.min)*maxThresholdRation);
	        
	        Binarization.globalThreshold(inRaster,bm,threshold);
	        stats.threshold=threshold;
	        cons.accept(stats);
	        
	        return bm;
	}
}
