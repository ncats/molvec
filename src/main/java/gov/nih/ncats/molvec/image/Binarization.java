package gov.nih.ncats.molvec.image;

import java.awt.image.Raster;
import java.util.function.Consumer;

import gov.nih.ncats.molvec.Bitmap;
import gov.nih.ncats.molvec.image.binarization.ImageStats;

public interface Binarization {
    Bitmap binarize (Raster raster);
    default Bitmap binarize(Raster raster, Consumer<ImageStats> cons){
    	return this.binarize(raster);
    }
}
