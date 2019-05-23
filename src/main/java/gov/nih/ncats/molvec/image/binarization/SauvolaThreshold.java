package gov.nih.ncats.molvec.image.binarization;

import java.awt.image.Raster;
import java.util.Arrays;

import gov.nih.ncats.molvec.Bitmap;
import gov.nih.ncats.molvec.image.Binarization;

public class SauvolaThreshold implements Binarization{

    public static final double DEFAULT_MIN_THRESHOLD_RATIO = 0.1;
    public static final double DEFAULT_MAX_THRESHOLD_RATIO = 0.9;
	int rad = 10;
	double k=-0.9;
	double r=128;
	
	
	@Override
	public Bitmap binarize(Raster inRaster) {
		
		Bitmap bm = new Bitmap (inRaster.getWidth (), inRaster.getHeight ());
		double max = -1, min = Double.MAX_VALUE;
        for (int y = 0; y < bm.height(); ++y) {
            for (int x = 0; x < bm.width(); ++x) {
                double pel = inRaster.getSampleDouble (x, y, 0);
                if (pel < min)
                    min = pel;
                if (pel > max)
                    max = pel;
            }
        }
        
        double t1=min+(max-min)*DEFAULT_MIN_THRESHOLD_RATIO;
        double t2=min+(max-min)*DEFAULT_MAX_THRESHOLD_RATIO;
        
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
	            	
	            	isOn = pix > mean * (1 + k * (Math.sqrt(var)/r-1.0));
            	}
            	
            	bm.set (x, y, isOn);
            }
        }
        return bm;
	}

}
