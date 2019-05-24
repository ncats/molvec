package gov.nih.ncats.molvec.image;

import java.awt.image.Raster;
import java.util.function.Consumer;
import java.util.function.Predicate;

import gov.nih.ncats.molvec.Bitmap;
import gov.nih.ncats.molvec.image.binarization.ImageStats;

public interface Binarization {
    Bitmap binarize (Raster raster);
    default Bitmap binarize(Raster raster, Consumer<ImageStats> cons){
    	return this.binarize(raster);
    }
    
    default Binarization fallback(Binarization bb, Predicate<ImageStats> fallif){
    	Binarization _this=this;
    	return new Binarization(){

			@Override
			public Bitmap binarize(Raster raster, Consumer<ImageStats> cc) {
				ImageStats[] is = new ImageStats[]{null};
				Bitmap bm=_this.binarize(raster, ss->is[0]=ss);
				
				if(is[0]!=null){
					if(fallif.test(is[0])){
						cc.accept(is[0]);
						return bb.binarize(raster);
					}else{
						cc.accept(is[0]);
					}
				}
				
				return bm;
			}

			@Override
			public Bitmap binarize(Raster raster) {
				return binarize(raster,(is)->{});
			}
    		
    	};
    }
}
