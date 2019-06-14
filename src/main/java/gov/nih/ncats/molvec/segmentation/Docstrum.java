package gov.nih.ncats.molvec.segmentation;

import java.awt.Shape;
import java.awt.geom.*;

import java.util.*;
import java.util.logging.Logger;

import java.io.*;

import gov.nih.ncats.molvec.algo.CentroidEuclideanMetric;
import gov.nih.ncats.molvec.algo.NearestNeighbors;
import gov.nih.ncats.molvec.algo.Peaks;
import gov.nih.ncats.molvec.algo.UnionFind;
import gov.nih.ncats.molvec.util.GeomUtil;

/**
 * An implementation of the bottom-up docstrum page segmenation
 * algorithm
 */
public class Docstrum implements Segmentation {
    private static final Logger logger = 
        Logger.getLogger(Docstrum.class.getName());

    static final boolean DEBUG;
    static {
        boolean debug = false;
        try {
            debug = Boolean.getBoolean("molvec.debug");
        }
        catch (Exception ex) {}
        DEBUG = debug;
    }

    interface Closure<T> {
        boolean connected (T v1, T v2);
    }

    static final int DEFAULT_K = 5;
    
    protected Collection<Shape> zones = new ArrayList<Shape>();
    protected int[] angle = new int[181]; // [0, 180]
    protected int[] dist;
    protected int[] docstrum;

    protected int[] withinLineDist;
    protected int[] betweenLineDist;
    protected int[] overallDist;

    protected int skew = 0; // estimated skew of image

    // +/- window around the peaks in angle histogram used 
    //  to estimate within line zones
    protected int withinRange = 10; 

    // +/- window around the peaks in angle histogram used
    //  to estimate between line zones
    protected int betweenRange = 10;

    public Docstrum () {
    }

    public Docstrum (Collection<Shape> ccomps) {
        analyze (ccomps);
    }

    public Docstrum (NearestNeighbors<Shape> knn) {
        estimateParameters (knn);
    }

    /**
     * perform docstrum analysis from a collection of connected
     * components
     */
    public void analyze (Collection<Shape> ccomps) {
        analyze (DEFAULT_K, ccomps);
    }

    public void analyze (int K, Collection<Shape> ccomps) {
        NearestNeighbors<Shape> knn = new NearestNeighbors<Shape>
            (K, new CentroidEuclideanMetric<Shape>());
        knn.addAll(ccomps);
        estimateParameters (knn);
    }

    /**
     * bounded angle within [-90, 90]
     */
    static int bound (int theta) {
        if (theta < 90) {
        }
        else if (theta < 270) {
            theta = 180 - theta;
        }
        else { // > 270
            theta = -(360 - theta);
        }
        return theta;
    }

    static int angle (double x0, double y0, double x1, double y1) {
        return (int)(Math.toDegrees(GeomUtil.angle(x0, y0, x1, y1))+.5);
    }

    static double toDegrees (Shape s1, Shape s2) {
        int x0 = (int)s1.getBounds2D().getCenterX();
        int y0 = (int)s1.getBounds2D().getCenterY();
        int x1 = (int)s2.getBounds2D().getCenterX();
        int y1 = (int)s2.getBounds2D().getCenterY();
        return Math.toDegrees(GeomUtil.angle(x0, y0, x1, y1));
    }

    public void estimateParameters (NearestNeighbors<Shape> knn) {
        for (int i = 0; i < angle.length; ++i) angle[i] = 0;
        /*
         * calculate neighbor statistics
         */
        int max = 0, k = 0;
        docstrum = new int[knn.size()*knn.getMaxNeighbors()];
        for (Shape s0 : knn.entries()) {
            double x0 = s0.getBounds2D().getCenterX();
            double y0 = s0.getBounds2D().getCenterY();
            for (NearestNeighbors.Neighbor<Shape> nb : knn.neighborList(s0)) {
                Shape s1 = nb.getNeighbor();
                double x1 = s1.getBounds2D().getCenterX();
                double y1 = s1.getBounds2D().getCenterY();
                int phi = angle (x0, y0, x1, y1); // [0,360]
                ++angle[bound (phi) + 90];
                int rho = (int)(nb.getValue() + .5); // distance
                if (rho > max) {
                    max = rho;
                }
                docstrum[k++] = (phi << 16) | rho;
            }
        }

        if (max > 0xffff) {
            max = 0xffff;
            logger.warning("Max nearest neighbor distance exceeds limit!");
        }

        dist = new int[max+1];
        for (int i = 0; i < docstrum.length; ++i) {
            ++dist[docstrum[i] & 0xffff];
        }

        Peaks p = new Peaks (5, 1.);
        int[] peaks = p.detect(dist);
        //System.out.print("** NNDist: "+peaks.length+" peaks detected;");
        for (int i = 0; i < peaks.length; ++i) {
            System.out.print(" "+peaks[i]);
        }
        System.out.println();

        /**
         * estimated skew
         */
        peaks = p.detect(angle);
        skew = peaks[0];
        if (skew < 0) {
            skew = 90 + skew;
        }

        //System.out.print("** NNAngle: "+peaks.length+" peaks detected;");
        for (int i = 0; i < peaks.length; ++i) {
            System.out.print(" "+(peaks[i] -90));
        }
        System.out.println();

        zones = transitiveClosure 
            (knn, new Closure<Shape> () {
                public boolean connected (Shape s1, Shape s2) {
                    if (s1 != null && s2 != null) {
                        int r = bound ((int)(toDegrees (s1, s2) + .5));
                        return Math.abs(r) <= 35;
                    }
                    return false;
                }
            });

        debug ();
    }

    /**
     * transitive closure based on nearest neighbor pairs within
     * designed lower and upper bounds on angle
     */
    Collection<Shape> transitiveClosure (NearestNeighbors<Shape> knn, 
                                         Closure<Shape> closure) {
        Map<Shape, Integer> g = new HashMap<Shape, Integer>();
        Map<Integer, Shape> ig = new HashMap<Integer, Shape>();

        int id = 0;
        UnionFind eqv = new UnionFind (knn.size());

        for (Shape s : knn.entries()) {
            Integer v0 = g.get(s);
            if (v0 == null) {
                g.put(s, v0 = id++);
                ig.put(v0, s);
            }

            Shape nb = knn.nearest(s);
            Integer v1 = g.get(nb);
            if (v1 == null) {
                g.put(nb, v1 = id++);
                ig.put(v1, nb);
            }

            if (closure.connected(s, nb)) {
                eqv.union(v0, v1);
            }
        }

        int[][] comps = eqv.getComponents();
        //logger.info(comps.length+" equivalence classes out of "+knn.size());

        Shape[] areas = new Shape[comps.length];
        for (int i = 0; i < comps.length; ++i) {
            Rectangle2D a = null;
            for (int j = 0; j < comps[i].length; ++j) {
                Rectangle2D r = ig.get(comps[i][j]).getBounds();
                if (a == null) {
                    a = r;
                }
                else {
                    a.add(r);
                }
                //System.out.println("  "+comps[i][j]+": "+a.getBounds());
            }
            areas[i] = a;
            //logger.info("## Component "+i+" "+a);
        }

        return Arrays.asList(areas);
    }

    void debug () {
        try { 
            PrintStream ps = new PrintStream
                (new FileOutputStream ("docstrum.txt"));
            ps.println("# "+docstrum.length+" data points");
            for (int i = 0; i < docstrum.length; ++i) {
                int phi = docstrum[i] >> 16;
                int rho = docstrum[i] & 0xffff;
                ps.println(phi+" "+rho);
            }
            ps.close();

            ps = new PrintStream (new FileOutputStream ("angles.txt"));
            for (int i = 0; i < angle.length; ++i) {
                if (angle[i] > 0) {
                    ps.println(String.format("%1$3d %2$d", i-90, angle[i]));
                }
            }
            ps.close();

            ps = new PrintStream (new FileOutputStream ("distances.txt"));
            for (int i = 0; i < dist.length; ++i) {
                if (dist[i] > 0) {
                    ps.println(i + " "+dist[i]);
                }
            }
            ps.close();
        }
        catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    /**
     * Segmentation interface
     */
    public Collection<Shape> getZones () {
        return zones;
    }
}
