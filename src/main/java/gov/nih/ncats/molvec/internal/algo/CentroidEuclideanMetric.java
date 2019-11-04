package gov.nih.ncats.molvec.internal.algo;

import java.awt.Shape;

public class CentroidEuclideanMetric<T extends Shape> implements Metric<T> {
    public CentroidEuclideanMetric () {}

    public double evaluate (T s0, T s1) {
        double x0 = s0.getBounds2D().getCenterX();
        double y0 = s0.getBounds2D().getCenterY();
        double x1 = s1.getBounds2D().getCenterX();
        double y1 = s1.getBounds2D().getCenterY();
        double xx = x1 - x0, yy = y1 - y0;
        return Math.sqrt(xx*xx + yy*yy);
    }
}
