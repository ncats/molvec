package gov.nih.ncats.molvec.image.binarization;

import java.awt.image.Raster;
import java.util.function.Consumer;

import gov.nih.ncats.molvec.Bitmap;
import gov.nih.ncats.molvec.image.Binarization;

/**
 * Simple threshold based on a constant cutoff
 */
public class ConstantThreshold implements Binarization {
    static final double DEFAULT_THRESHOLD = 128.;

    private double threshold;
    
    public ConstantThreshold () {
        this (DEFAULT_THRESHOLD);
    }

    public ConstantThreshold (double threshold) {
        this.threshold = threshold;
    }

    public void setThreshold (double threshold) { this.threshold = threshold; }
    public double getThreshold () { return threshold; }



	@Override
	public Bitmap binarize(Raster inRaster, ImageStats stats, Consumer<ImageStats> cons) {
		Bitmap bm = new Bitmap (inRaster.getWidth (), inRaster.getHeight ());
        
        Binarization.globalThreshold(inRaster, bm, threshold);
        
        return bm;
	}
}
