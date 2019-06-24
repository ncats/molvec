package gov.nih.ncats.molvec.segmentation;

import java.util.Collection;
import java.awt.Shape;


/**
 * A generic segmentation interface
 */
public interface Segmentation {
    Collection<Shape> getZones ();
}
