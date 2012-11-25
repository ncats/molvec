package tripod.molvec.segmentation;

import java.awt.Shape;
import java.awt.geom.*;

import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;

import java.io.*;

import tripod.molvec.Zone;
import tripod.molvec.Segmentation;
import tripod.molvec.algo.*;
import tripod.molvec.util.GeomUtil;

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

    static final int DEFAULT_K = 5;
    
    protected List<Zone> zones = new ArrayList<Zone>();
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
        System.out.print("** NNDist: "+peaks.length+" peaks detected;");
        for (int i = 0; i < peaks.length; ++i) {
            System.out.print(" "+peaks[i]);
        }
        System.out.println();

        /**
         * estimate skew
         */
        peaks = p.detect(angle);
        System.out.print("** NNAngle: "+peaks.length+" peaks detected;");
        for (int i = 0; i < peaks.length; ++i) {
            System.out.print(" "+(peaks[i] -90));
        }
        System.out.println();

        debug ();
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
    public Collection<Zone> getZones () {
        return zones;
    }
}
