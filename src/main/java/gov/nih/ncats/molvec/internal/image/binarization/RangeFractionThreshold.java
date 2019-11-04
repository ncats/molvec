package gov.nih.ncats.molvec.internal.image.binarization;

import java.awt.image.Raster;
import java.util.function.Consumer;

import gov.nih.ncats.molvec.internal.image.Bitmap;

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
	        
		if(stats==null){
		    stats = Binarization.computeImageStats(inRaster);
        }
	     

        stats.threshold=stats.min + (stats.max-stats.min)*pct;
        
        
        Binarization.globalThreshold(inRaster, bm, stats.threshold);
        
        cons.accept(stats);
        
        return bm;
	}
}
