package gov.nih.ncats.molvec.image;

import java.awt.image.Raster;
import gov.nih.ncats.molvec.Bitmap;

public interface Binarization {
    Bitmap binarize (Raster raster);
}
