package tripod.molvec.image.binarization;

import java.awt.image.Raster;
import tripod.molvec.Bitmap;
import tripod.molvec.image.Binarization;

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

    public Bitmap binarize (Raster inRaster) {
        Bitmap bm = new Bitmap (inRaster.getWidth (), inRaster.getHeight ());
        for (int y = 0; y < bm.height(); ++y) {
            for (int x = 0; x < bm.width(); ++x) {
                double pel = inRaster.getSampleDouble (x, y, 0);
                bm.set (x, y, pel >= threshold);
            }
        }

        return bm;
    }
}
