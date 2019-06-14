package gov.nih.ncats.molvec.image.binarization;

import java.awt.image.Raster;
import java.util.Arrays;
import java.util.function.Consumer;

import gov.nih.ncats.molvec.image.Bitmap;

/**
 * Implementation of Sauvola threshold. 
 * @author tyler
 *
 */
public class SauvolaThreshold implements Binarization{

    public static final double DEFAULT_MIN_THRESHOLD_RATIO = 0.1;
    public static final double DEFAULT_MAX_THRESHOLD_RATIO = 0.9;
	int rad = 10;
	double k=-0.9;
	double r=128;
	
	public SauvolaThreshold(int rad, double k, double r){
		this.rad=rad;
		this.k=k;
		this.r=r;
		
	}
	
	public SauvolaThreshold(){
		
	}
	

	@Override
	public Bitmap binarize(Raster inRaster, ImageStats stats, Consumer<ImageStats> cons) {

		Bitmap bm = new Bitmap (inRaster.getWidth (), inRaster.getHeight ());
		if(stats==null)stats = Binarization.computeImageStats(inRaster);
        		
		double ek = k;
        if(stats.mean-stats.min > stats.max-stats.mean){
        	ek=ek*-1;
        }
        
        double t1=stats.min+(stats.max-stats.min)*DEFAULT_MIN_THRESHOLD_RATIO;
        double t2=stats.min+(stats.max-stats.min)*DEFAULT_MAX_THRESHOLD_RATIO;
        
		double[] vals = new double[rad*rad*4];
		
        for (int y = 0; y < bm.height(); ++y) {
            for (int x = 0; x < bm.width(); ++x) {
            	double pix = inRaster.getSampleDouble (x, y, 0);
            	boolean isOn = false;
            	if(pix>t2){
            		isOn=true;
            	}else if(pix<t1){
            		isOn=false;
            	}else{
	            	int minx = Math.max(x-rad, 0);
	            	int miny = Math.max(y-rad, 0);
	            	int maxx = Math.min(x+rad, bm.width()-1);
	            	int maxy = Math.min(y+rad, bm.height()-1);
	            	
	            	int tot = (maxx-minx)*(maxy-miny);
	            	
	            	inRaster.getSamples(minx, miny, maxx-minx, maxy-miny, 0, vals);
	            	
	            	
	            	double mean = Arrays.stream(vals).limit(tot).average().orElse(0);
	            	double var = Arrays.stream(vals).limit(tot).map(d->mean-d).map(d->d*d).average().orElse(0);
	            	
	            	isOn = pix > mean * (1 + ek * (Math.sqrt(var)/r-1.0));
            	}
            	
            	bm.set (x, y, isOn);
            }
        }
        cons.accept(stats);
        return bm;
	}

}
