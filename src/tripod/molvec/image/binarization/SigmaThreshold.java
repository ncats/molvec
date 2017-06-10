package tripod.molvec.image.binarization;

import java.awt.image.Raster;
import tripod.molvec.Bitmap;
import tripod.molvec.image.Binarization;

/*
 * Simple threshold based on full image mean and standard deviation. 
 */
public class SigmaThreshold implements Binarization {
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
        double threshold = mean + stdDEV * sigma;
        for (int y = 0; y < bm.height(); ++y) {
            for (int x = 0; x < bm.width(); ++x) {
                double pel = inRaster.getSampleDouble (x, y, 0);
                bm.set (x, y, pel >= threshold);
            }
        }

        return bm;
    }
}
