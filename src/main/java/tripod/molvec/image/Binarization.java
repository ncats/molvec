package tripod.molvec.image;

import java.awt.image.Raster;
import tripod.molvec.Bitmap;

public interface Binarization {
    Bitmap binarize (Raster raster);
}
