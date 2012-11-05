package tripod.molvec.segmentation;

import java.awt.Shape;
import java.awt.geom.*;
import java.util.*;

import tripod.molvec.Zone;
import tripod.molvec.Segmentation;

/**
 * An implementation of the bottom-up docstrum page segmenation
 * algorithm
 */
public class Docstrum implements Segmentation {

    protected List<Zone> zones = new ArrayList<Zone>();

    public Docstrum () {
    }

    public Docstrum (Collection<Shape> ccomps) {
        analyze (ccomps);
    }

    /**
     * perform docstrum analysis from a collection of connected
     * components
     */
    public void analyze (Collection<Shape> ccomps) {
    }

    /**
     * Segmentation interface
     */
    public Collection<Zone> getZones () {
        return zones;
    }
}
