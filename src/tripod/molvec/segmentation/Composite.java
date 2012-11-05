package tripod.molvec.segmentation;

import java.io.Serializable;
import java.util.*;
import java.awt.Shape;
import java.awt.Polygon;
import java.awt.geom.*;

import java.util.logging.Logger;
import java.util.logging.Level;

import tripod.molvec.util.GeomUtil;

public class Composite implements Serializable {
    static final Logger logger = 
        Logger.getLogger(Composite.class.getName());

    protected Collection<Shape> bboxes; // connected component
    protected Collection<Shape> composites;
    protected double cutoff = -1.;

    public Composite (Collection<Shape> bboxes) {
        this.bboxes = bboxes;

        composites = bboxes;
        do {
            // this is extremely inefficient; not suitable for
            //  big images
            Collection<Shape> cc = createComposites (composites);
            if (cc.size() < composites.size()) {
                composites = cc;
            }
            else {
                break;
            }
        }
        while (true);
    }

    public Collection<Shape> getComposites () { return composites; }
    public Collection<Shape> getBoundingBoxes () { return bboxes; }

    /**
     * Given a collection of polygons, this routine generates composite
     * fragments based on containment and intersections of polygons.
     */
    public Collection<Shape> createComposites (Collection<Shape> polygons) {
        Shape[] poly = polygons.toArray(new Shape[0]);
        // first pass, collapse containments and intersections
        logger.info("## collapsing polygons...");
        do {
            int k = 0;
            for (int i = 0; i < poly.length; ++i) {
                if (poly[i] != null) {
                    k += composites (i, poly);
                }
            }

            if (k == 0) 
                break;
        }
        while (true);

        List<Shape> composites = new ArrayList<Shape>();
        for (int i = 0; i < poly.length; ++i) {
            if (poly[i] != null) {
                composites.add(poly[i]);
            }
        }

        logger.info("## calculating polygon nearest neighbor...");
        // calculate average nearest neighbor distance
        poly = composites.toArray(new Shape[0]);
        int[] nb = new int[poly.length]; // neighbors
        double[] dist = new double[poly.length];

        for (int i = 0; i < poly.length; ++i) {
            nb[i] = -1;
            dist[i] = Double.MAX_VALUE;
            for (int j = i+1; j < poly.length; ++j) {
                double d = GeomUtil.distance(poly[i], poly[j]);
                if (d < dist[i]) {
                    nb[i] = j;
                    dist[i] = d;
                }
            }
        }

        double[] median = new double[dist.length];
        System.arraycopy(dist, 0, median, 0, dist.length);
        Arrays.sort(median);
        
        if (cutoff < 0) {
            cutoff = median.length % 2 == 0 
                ? ((median[median.length/2] + median[median.length/2-1])/2.)
                : median[median.length/2];
        }
        
        logger.info("## merging nearest neighbor at cutoff "+cutoff+"..");
        composites.clear();
        for (int i = 0; i < poly.length; ++i) {
            if (poly[i] != null) {
                BitSet group = new BitSet (poly.length);
                dfs (cutoff, group, i, dist, nb);

                for (int j = group.nextSetBit(0); j >= 0; 
                     j = group.nextSetBit(j+1)) {
                    if (i != j && poly[j] != null) {
                        poly[i] = GeomUtil.add(poly[i], poly[j]);
                        poly[j] = null;
                        dist[j] = Double.MAX_VALUE;
                    }
                }

                composites.add(poly[i]);
            }
        }

        return composites;
    }

    static void dfs (double cutoff, BitSet group, 
                     int index, double[] dist, int[] nb) {
        group.set(index);
        int next = nb[index];
        if (next >= 0 && dist[index] < cutoff) {
            dfs (cutoff, group, next, dist, nb);
        }
    }

    static int composites (int index, Shape[] polygons) {
        Shape poly = polygons[index];
        Rectangle2D bounds = poly.getBounds2D();

        int k = 0;
        for (int i = 0; i < polygons.length; ++i) {
            Shape p = polygons[i];
            if (i != index && p != null) {
                if (p.contains(bounds)) {
                    polygons[index] = null;
                    ++k;
                    break;
                }
                else if (GeomUtil.intersects(poly, p)) {
                    polygons[index] = poly = GeomUtil.add(poly, p);
                    polygons[i] = null;
                    ++k;
                }
                else if (GeomUtil.contains(poly, p)) {
                    polygons[i] = null;
                    ++k;
                }
                else if (GeomUtil.contains(p, poly)) {
                    polygons[index] = null;
                    ++k;
                    break;
                }
            }
        }
        return k;
    }
}
