package gov.nih.ncats.molvec.image;

import java.awt.image.Raster;
import java.util.function.Consumer;
import java.util.function.Predicate;

import gov.nih.ncats.molvec.image.binarization.ImageStats;

public interface Binarization {
    Bitmap binarize (Raster raster, ImageStats stats, Consumer<ImageStats> cons);
    
    default Binarization fallback(Binarization bb, Predicate<ImageStats> fallif){
    	Binarization _this=this;
    	return new Binarization(){

			@Override
			public Bitmap binarize(Raster raster,ImageStats stats, Consumer<ImageStats> cc) {
				ImageStats[] is = new ImageStats[]{null};
				Bitmap bm=_this.binarize(raster, stats, ss->is[0]=ss);
				
				if(is[0]!=null){
					if(fallif.test(is[0])){
						return bb.binarize(raster,is[0], cc);
					}else{
						cc.accept(is[0]);
					}
				}
				return bm;
			}
    		
    	};
    }
    
    
    public static ImageStats computeImageStats(Raster inRaster){
	    	int width = inRaster.getWidth();
	    	int height = inRaster.getHeight();
	    	ImageStats stats = new ImageStats();
	    	
	    	int[] fstats = new int[1000];
	    	
    	
    	
    	   double max = -1, min = Double.MAX_VALUE;
           double sum = 0;
           
           double sumSquare = 0;
           double[] nd = new double[width];
           
           for (int y = 0; y < height; ++y) {
        	   inRaster.getSamples(0, y, width, 1, 0, nd);
               for (int x = 0; x < width; ++x) {
                   double pel = nd[x];
                   sum += pel;
                   sumSquare += pel * pel;
                   if (pel < min)
                       min = pel;
                   if (pel > max)
                       max = pel;
                   fstats[(int)pel]++;
               }
           }
           long tot = width * height;
           double mean = sum / tot;
           double stdDEV = Math.sqrt (sumSquare / tot - mean * mean);
           
           stats.histogram = new int[101];
           stats.histogramRaw = new int[(int)max+1];
           
           
           for(int i=(int)min;i<=(int)max;i++){
        	   int ni = (int)((100*(i-min))/(max-min));
        	   stats.histogram[ni] +=fstats[i];
        	   stats.histogramRaw[i] = fstats[i];
           }
           
//           for(int i=0;i<stats.histogramRaw.length;i++){
//        	   System.out.println(i + "\t" + stats.histogramRaw[i]);
//           }
//         
           stats.min=min;
           stats.max=max;
           stats.mean=mean;
           stats.stdev=stdDEV;
           stats.count=width*height;
           stats.threshold=stats.mean;
           return stats;
    }
    

    public static void globalThreshold(Raster inRaster, Bitmap bm, double threshold){
    	 double[] nd = new double[bm.width()];
         
         for (int y = 0; y < bm.height(); ++y) {
      	   inRaster.getSamples(0, y, bm.width(), 1, 0, nd);
             for (int x = 0; x < bm.width(); ++x) {
                 bm.set (x, y, nd[x] >= threshold);
             }
         }
    }
}
