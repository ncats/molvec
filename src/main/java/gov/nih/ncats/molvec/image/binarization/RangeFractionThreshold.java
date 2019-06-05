package gov.nih.ncats.molvec.image.binarization;

import java.awt.image.Raster;
import java.util.function.Consumer;
import java.util.logging.Logger;

import gov.nih.ncats.molvec.Bitmap;
import gov.nih.ncats.molvec.image.Binarization;

/*
 * Simple threshold based on full image mean and standard deviation. 
 */
public class RangeFractionThreshold implements Binarization {
    public static double DEFAULT_FRACTION = 0.26;
    
    private double pct;
    
    public RangeFractionThreshold () {
        this (DEFAULT_FRACTION);
    }

    public RangeFractionThreshold (double pct) {
        this.pct = pct;
    }

    

	@Override
	public Bitmap binarize(Raster inRaster, ImageStats stats, Consumer<ImageStats> cons) {
		Bitmap bm = new Bitmap (inRaster.getWidth (), inRaster.getHeight ());
	        
		if(stats==null)stats = Binarization.computeImageStats(inRaster);
	     
		double threshold = stats.mean;
	        
        threshold = stats.min + (stats.max-stats.min)*pct;
        stats.threshold=threshold;
        
        
        Binarization.globalThreshold(inRaster, bm, threshold);
        
        cons.accept(stats);
        
        return bm;
	}
}
