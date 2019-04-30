package gov.nih.ncats.molvec.image.binarization;

import java.awt.image.Raster;
import gov.nih.ncats.molvec.Bitmap;
import gov.nih.ncats.molvec.image.Binarization;

/**
 * Simple threshold based on a constant cutoff
 */
public class WindowThreshold implements Binarization {
    private double low, high;
    
    public WindowThreshold () {
        this (0, 255);
    }

    public WindowThreshold (double low, double high) {
        this.low = low;
        this.high = high;
    }

    public void setLow (double low) { this.low = low; }
    public double getLow () { return low; }
    public void setHigh (double high) { this.high = high; }
    public double getHigh () { return high; }

    public Bitmap binarize (Raster inRaster) {
        Bitmap bm = new Bitmap (inRaster.getWidth (), inRaster.getHeight ());
        for (int y = 0; y < bm.height(); ++y) {
            for (int x = 0; x < bm.width(); ++x) {
                double pel = inRaster.getSampleDouble (x, y, 0);
                bm.set (x, y, pel >= low && pel <= high);
            }
        }

        return bm;
    }
}