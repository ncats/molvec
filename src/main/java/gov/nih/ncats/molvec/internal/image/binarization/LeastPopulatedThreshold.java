package gov.nih.ncats.molvec.internal.image.binarization;

import java.awt.image.Raster;
import java.util.function.Consumer;

import gov.nih.ncats.molvec.internal.image.Bitmap;

/*
 * Simple threshold based on full image mean and standard deviation. 
 */
public class LeastPopulatedThreshold implements Binarization {
    
    private int window;
    
    public LeastPopulatedThreshold () {
        this (10);
    }

    public LeastPopulatedThreshold (int window) {
        this.window = window;
    }

    
    public Bitmap binarize (Raster inRaster) {
    	return binarize(inRaster,(is)->{});
    }
    public Bitmap binarize (Raster inRaster, Consumer<ImageStats> cons) {

        
        return binarize(inRaster, null, cons);
    }

	@Override
	public Bitmap binarize(Raster inRaster, ImageStats stats, Consumer<ImageStats> cons) {
		Bitmap bm = new Bitmap (inRaster.getWidth (), inRaster.getHeight ());
        
		if(stats==null){
		    stats= Binarization.computeImageStats(inRaster);
        }
        

        
        int stot = 0;
        int mini = 55;
        int minVal = Integer.MAX_VALUE;
        for(int i=0;i<100;i++){
        	stot+=stats.histogram[i];
        	if(i>window){
        		stot-=stats.histogram[i-window];
        	}
        	if(i>=window){
        		if(stot<minVal){
        			minVal = stot;
        			mini=i;
        		}
        	}
        }

        stats.threshold=stats.min + (stats.max-stats.min)*(mini-window*0.5)/100D;
        
        
//        System.out.println("Threshold Num:" + mini);
//        System.out.println("Threshold:" + threshold);
        
        Binarization.globalThreshold(inRaster,bm,stats.threshold);
        
        cons.accept(stats);
        
        return bm;
	}
    
}
