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
      
        
        double threshold = mean;
        
        stats.min=min;
        stats.max=max;
        stats.mean=mean;
        stats.stdev=stdDEV;
        stats.count=bm.width()*bm.height();
        
        threshold = stats.min + (stats.max-stats.min)*pct;
        stats.threshold=threshold;
        
        
        
        
        for (int y = 0; y < bm.height(); ++y) {
            for (int x = 0; x < bm.width(); ++x) {
                double pel = inRaster.getSampleDouble (x, y, 0);
                bm.set (x, y, pel >= threshold);
            }
        }
      
        cons.accept(stats);
        
        return bm;
    }
}
