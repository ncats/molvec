package tripod.molvec.util;

import java.awt.Point;
import java.awt.Polygon;
import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.PathIterator;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import tripod.molvec.algo.LineUtil;
import tripod.molvec.algo.Tuple;


public class GeomUtil {
    private static final Logger logger = 
        Logger.getLogger (GeomUtil.class.getName ());

    private static final boolean DEBUG;
    public static final double EPS = 0.0001;

    static {
        boolean debug = false;
        try {
            debug = Boolean.getBoolean ("geomutil.debug");
        } catch (Exception ex) {
        }
        DEBUG = debug;
    }

    public static double ccw (Point2D p1, Point2D p2, Point2D p3) {
        return (p2.getX () - p1.getX ()) * (p3.getY () - p1.getY ())
            - (p2.getY () - p1.getY ()) * (p3.getX () - p1.getX ());
    }

    public static double angle (Point2D p0, Point2D p1) {
        return angle (p0.getX(), p0.getY(), p1.getX(), p1.getY());
    }

    /**
     * calculate angle between (x0,y0) and (x1,y1) in radian
     */
    public static double angle (double x0, double y0, double x1, double y1) {
        double dx = x1 - x0, dy = y1 - y0;
        if (dx > 0 && dy > 0) {
            return Math.atan (dy / dx);
        } else if (dx > 0 && dy == 0) {
            return 0.;
        } else if (dx < 0 && dy > 0) {
            return Math.PI - Math.atan (-1. * dy / dx);
        } else if (dx < 0 && dy == 0) {
            return Math.PI;
        } else if (dx == 0 && dy > 0) {
            return Math.PI / 2;
        } else if (dx == 0 && dy < 0) {
            return 3 * Math.PI / 2;
        } else if (dx < 0 && dy < 0) {
            return 3 * Math.PI / 2 - Math.atan (dy / dx);
        } else if (dx > 0 && dy < 0) {
            return 2 * Math.PI - Math.atan (-1. * dy / dx);
        }
        return 0.;
    }


    /**
     * Graham scan algorithm for convex hull
     */
    public static Polygon convexHull (Point2D... pts) {
        if (pts.length < 3) {
            Polygon poly = new Polygon ();
            for (Point2D pt : pts) {
                poly.addPoint ((int) (pt.getX () + .5), (int) (pt.getY () + .5));
            }
            return poly;
        }

        Point2D anchor = null;
        for (Point2D pt : pts) {
            if (anchor == null || pt.getY () < anchor.getY ()) {
                anchor = pt;
            } else if (pt.getY () == anchor.getY ()
                       && pt.getX () < anchor.getX ()) {
                anchor = pt;
            }
        }

        final Point2D p0 = anchor;
        Arrays.sort (pts, new Comparator<Point2D> () {
                         public int compare (Point2D a, Point2D b) {
                             double a0 = angle (p0, a), a1 = angle (p0, b);
                             /*
                               System.err.println("p0=("+p0.x+","+p0.y+") a=("+a.x+","+a.y+") "
                               +"b=("+b.x+","+b.y+") ccw="+ccw(p0,a,b)
                               +" theta(p0,a)="+a0+" theta(p0,b)="+a1);
                             */
                             if (a0 < a1) return -1;
                             if (a0 > a1) return 1;

                             double d0 = a.distance (p0), d1 = b.distance (p0);
                             if (d0 < d1) return -1;
                             if (d0 > d1) return 1;
                             return 0;
                         }
                     });

        if (DEBUG) {
            logger.info ("Starting point: " + p0);
            logger.info ("Points..." + pts.length);
            for (int i = 0; i < pts.length; ++i) {
                System.err.println (i + ": " + pts[i]);
            }
        }

        LinkedList<Point2D> stack = new LinkedList<Point2D> ();
        stack.push (pts[0]);
        stack.push (pts[1]);
        stack.push (pts[2]);
        for (int i = 3; i < pts.length; ++i) {
            Point2D pi = pts[i];
            while (true) {
                Point2D p2 = stack.pop ();
                Point2D p1 = stack.peek ();
                double dir = ccw (p1, p2, pi);
                if (DEBUG) {
                    System.err.println ("p1=(" + p1.getX () + "," + p1.getY () + ") p2=("
                                        + p2.getX () + "," + p2.getY () + ") p"
                                        + i + "=(" + pi.getX () + "," + pi.getY ()
                                        + ") ccw=" + dir);
                }
                if (dir >= 0.) { // push back
                    stack.push (p2);
                    break;
                }
                if (DEBUG) {
                    logger.info ("removing " + p2);
                }
            }

            stack.push (pi);
        }

        if (DEBUG) {
            logger.info ("Convex hull: " + stack.size ());
        }

        Polygon hull = new Polygon ();
        for (ListIterator<Point2D> it = stack.listIterator (stack.size ());
             it.hasPrevious (); ) {
            Point2D pt = it.previous ();
            // should prune collinear points
            if (DEBUG) {
                System.err.println (" " + pt);
            }
            hull.addPoint ((int) (pt.getX () + .5), (int) (pt.getY () + .5));
        }

        return hull;
    }

    public static boolean isNeighbor (Point p1, Point p2) {
        return (Math.abs (p1.x - p2.x) <= 1 && Math.abs (p1.y - p2.y) <= 1);
    }

    public static Polygon toPolygon (Shape shape) {
        Polygon poly = new Polygon ();
        for (Point2D p : vertices (shape)) {
            poly.addPoint ((int) (p.getX () + .5), (int) (p.getY () + .5));
        }
        return poly;
    }

    // create a new convex hull shape from two given shapes
    public static Polygon add (Shape s1, Shape s2) {
        ArrayList<Point2D> pts = new ArrayList<Point2D> ();
        for (Point2D p : vertices (s1)) {
            pts.add (p);
        }
        for (Point2D p : vertices (s2)) {
            pts.add (p);
        }
        return convexHull (pts.toArray (new Point2D[0]));
    }

    public static Point2D[] vertices (Shape shape) {
        return vertices (shape, null);
    }

    public static Point2D[] vertices (Shape shape, AffineTransform afx) {
        PathIterator p = shape.getPathIterator (afx);
        List<Point2D> vertices = new ArrayList<Point2D> ();
        double[] coord = new double[6];
        while (!p.isDone ()) {
            int type = p.currentSegment (coord);
            if (type != PathIterator.SEG_CLOSE) {
                vertices.add (new Point2D.Double (coord[0], coord[1]));
            }
            p.next ();
        }

        return vertices.toArray (new Point2D[0]);
    }

    public static Line2D[] lines (Shape shape) {
        return lines (shape, null);
    }

    public static Line2D[] lines (Shape shape, AffineTransform afx) {
        List<Line2D> lines = new ArrayList<Line2D> ();

        double[] coord = new double[6];
        Point2D start = null, prev = null;

        PathIterator p = shape.getPathIterator (afx);
        while (!p.isDone ()) {
            int type = p.currentSegment (coord);
            if (type == PathIterator.SEG_MOVETO) {
                start = prev = new Point2D.Double (coord[0], coord[1]);
            } else if (type == PathIterator.SEG_CLOSE) {
                if (prev == null || start == null) {
                    logger.warning ("Unexpected state while "
                                    + "iterating over shape!");
                } else {
                    lines.add (new Line2D.Double (prev, start));
                }
            } else {
                Point2D pt = new Point2D.Double (coord[0], coord[1]);
                lines.add (new Line2D.Double (prev, pt));
                prev = pt;
            }
            p.next ();
        }

        return lines.toArray (new Line2D[0]);
    }

    public static Point2D segmentIntersection(Line2D l1, Line2D l2){
    	 Point2D pt = intersection (l1, l2);
         if (pt != null && intersects (l1, pt) && intersects (l2, pt)) {
             return pt;
         }
         return null;
    }
    
    /**
     * return the intersection point (if any) between two lines
     */
    public static Point2D intersection (Line2D l1, Line2D l2) {
        Point2D p1 = l1.getP1 (), p2 = l1.getP2 ();
        Point2D p3 = l2.getP1 (), p4 = l2.getP2 ();

        double c = (p1.getX () - p2.getX ()) * (p3.getY () - p4.getY ())
            - (p1.getY () - p2.getY ()) * (p3.getX () - p4.getX ());

        if (Math.abs (c) < 0.0001)
            return null;


        double x = (p1.getX () * p2.getY () - p1.getY () * p2.getX ())
            * (p3.getX () - p4.getX ()) - (p1.getX () - p2.getX ())
            * (p3.getX () * p4.getY () - p3.getY () * p4.getX ());
        double y = (p1.getX () * p2.getY () - p1.getY () * p2.getX ())
            * (p3.getY () - p4.getY ()) - (p1.getY () - p2.getY ())
            * (p3.getX () * p4.getY () - p3.getY () * p4.getX ());

        return new Point2D.Double (x / c, y / c);

        /*
          double[] p1 = calcLineParams (l1);
          double[] p2 = calcLineParams (l2);
          if (p1[0] != p2[0]) {
          double x = (p2[1] - p1[1])/(p1[0] - p2[0]);
          double y = p1[0]*x + p1[1];
          return new Point2D.Double(x, y);
          }
          return null;
        */
    }
    public static List<List<Point2D>> groupPointsCloserThan(List<Point2D> points,double maxDistance){
    	int[] groups = IntStream.range(0, points.size())
    							.toArray();
    	
    	for(int i=0;i<points.size();i++){
    		Point2D p1=points.get(i);
    		for(int j=i+1;j<points.size();j++){
    			Point2D p2=points.get(j);
    			if(p1.distance(p2)<maxDistance){
    				int g1= groups[i];
    				int g2= groups[j];
    				groups[i]=Math.min(g1, g2);
    				groups[j]=Math.min(g1, g2);
    			}
    		}
    	}
    	
    	
    	return IntStream.range(0, points.size())
    	         .mapToObj(i->Tuple.of(groups[i],i))
    	         .map(Tuple.vmap(i->points.get(i)))
    	         .collect(Tuple.toGroupedMap())
    	         .values()
    	         .stream()
    	         .collect(Collectors.toList());
    	
    }
    
    public static class LineDistanceCalculator{
    	private Line2D line1;
    	private Line2D line2;
    	
    	private int line1CloserPoint=-1;
    	private int line2CloserPoint=-1;
    	
    	private Point2D[] pts1=new Point2D[2];
    	private Point2D[] pts2=new Point2D[2];
    	
    	private double[] lens = new double[4];
    	
    	public boolean isProcessed(){
    		return line1CloserPoint!=-1;
    	}
    	public LineDistanceCalculator process(){
    		pts1[0]=line1.getP1();
    		pts1[1]=line1.getP2();
    		pts2[0]=line2.getP1();
    		pts2[1]=line2.getP2();
    		
    		lens[0]=pts1[0].distance(pts2[0]);
    		lens[1]=pts1[0].distance(pts2[1]);
    		lens[2]=pts1[1].distance(pts2[0]);
    		lens[3]=pts1[1].distance(pts2[1]);
    		
    		double min=Double.POSITIVE_INFINITY;
    		int mini=0;
    		for(int i=0;i<lens.length;i++){
    			if(lens[i]<min){
    				min=lens[i];
    				mini=i;
    			}
    		}
    		
    		line1CloserPoint=(mini&2)>>1;
    		line2CloserPoint=(mini&1);
    		return this;
    	}
    	
    	public double getSmallestPointDistance(){
    		if(!isProcessed())process();
    		int i=(line1CloserPoint<<1) | line2CloserPoint;
    		
    		return lens[i];
    	}
    	public Point2D[] closestPoints(){
    		Point2D[] points= new Point2D[2];
    		if(!isProcessed())process();
    		points[0]=pts1[line1CloserPoint];
    		points[1]=pts2[line2CloserPoint];
    		return points;
    	}
    	
    	public Line2D getLineFromFarthestPoints(){
    		if(!isProcessed())process();
    		//System.out.println("L1 close:" + line1CloserPoint + " and " + line2CloserPoint);
    		Point2D p1 = (line1CloserPoint==0)?line1.getP2():line1.getP1();
    		Point2D p2 = (line2CloserPoint==0)?line2.getP2():line2.getP1();
    		Line2D nline= new Line2D.Double(p1,p2);
    		double ol1=LineUtil.length(line1);
    		double ol2=LineUtil.length(line2);
    		double nl=LineUtil.length(nline);
    		if(nl>ol1 && nl>ol2){
    			return nline;
    		}else{
    			if(ol1>ol2)return line1;
    			return line2;
    		}
    		
    	}
    	
    	public static LineDistanceCalculator from(Line2D line1, Line2D line2){
    		LineDistanceCalculator ldc = new LineDistanceCalculator();
    		ldc.line1=line1;
    		ldc.line2=line2;
    		return ldc;
    	}
    	
    }
    

    // calculate line parameters (slope & intercept) of a line
    public static double[] calcLineParams (Line2D line) {
        Point2D p1 = line.getP1 (), p2 = line.getP2 ();
        double intercept = 0, slope = 0;
        if (p1.getX () != p2.getX ()) {
            intercept = (p2.getX () * p1.getY () - p1.getX () * p2.getY ())
                / (p2.getX () - p1.getX ());

            slope = (p1.getY () - intercept) / p1.getX ();
        }

        return new double[]{slope, intercept};
    }

    public static boolean intersects (Line2D line, Point2D pt) {
        return line.ptSegDist (pt) < 0.0001;
    }

    public static boolean intersects (Shape s1, Shape s2) {
        Line2D[] lines1 = lines (s1);
        Line2D[] lines2 = lines (s2);
        for (Line2D l1 : lines1) {
            for (Line2D l2 : lines2) {
                Point2D pt = intersection (l1, l2);
                if (pt != null && intersects (l1, pt) && intersects (l2, pt)) {
                    return true;
                }
            }
        }
        return false;
    }
    
    public static Point2D getIntersection(Shape s1, Line2D l1){
    	
    	return Arrays.stream(lines(s1))
    	      .map(l->segmentIntersection(l,l1))
    	      .filter(p->p!=null)
    	      .findFirst()
    	      .orElse(null);
    }

    /**
     * Test whether s2 is contained within s1
     */
    public static boolean contains (Shape s1, Shape s2) {
        Point2D[] vertices = vertices (s2);
        int k = 0;
        for (Point2D p : vertices) {
            if (s1.contains (p))
                ++k;
        }
        return k == vertices.length;
    }

    /**
     * Return two closest vertices between two shapes
     */
    public static Point2D[] nearestNeighborVertices (Shape s1, Shape s2) {
        Point2D[] pts1 = vertices (s1);
        Point2D[] pts2 = vertices (s2);

        double minDist = Double.MAX_VALUE;
        Point2D p1 = null, p2 = null;
        for (int i = 0; i < pts1.length; ++i) {
            for (int j = 0; j < pts2.length; ++j) {
                double dist = length (pts1[i], pts2[j]);
                if (dist < minDist) {
                    p1 = pts1[i];
                    p2 = pts2[j];
                    minDist = dist;
                }
            }
        }

        return new Point2D[]{p1, p2};
    }
    
    public static double distanceTo(Shape s, Point2D pt){
    	if(s.contains(pt)){
    		return 0;
    	}
    	return Arrays.stream(lines(s))
    	      .mapToDouble(l->l.ptSegDist(pt))
    	      .min()
    	      .getAsDouble();
    }
    
    public static Tuple<Shape,Double> distanceTo(List<Shape> shapes, Point2D pt){
    	return shapes.stream()
    	      .map(s->Tuple.of(s,distanceTo(s,pt)).withVComparator())
    	      .sorted()
    	      .findFirst()
    	      .get();
    }
    
  
  

    /**
     * Euclidean distance between two points
     */
    public static double length (Point2D pt1, Point2D pt2) {
        double dx = pt1.getX () - pt2.getX (), dy = pt1.getY () - pt2.getY ();
        return Math.sqrt (dx * dx + dy * dy);
    }

    /**
     * Euclidean between two polygons based on nearest vertices distance
     */
    public static double distance (Shape s1, Shape s2) {
        Point2D[] vertex = nearestNeighborVertices (s1, s2);
        return length (vertex[0], vertex[1]);
    }
    
    public static double[] centerVector(double[] v){
    	double a1=(IntStream.range(0,v.length)
   			 .mapToDouble(i->v[i])
   			 .average()
   			 .getAsDouble());
    	
    	
    	return IntStream.range(0,v.length)
			 .mapToDouble(i->v[i]-a1)
			 .toArray();
    }
    
    public static double l2Norm(double[] v){
    	return Math.sqrt(IntStream.range(0,v.length)
   			 .mapToDouble(i->v[i])
   			 .map(d->d*d)
   			 .sum());
    }
    public static double dot(double[] a, double[] b){
    	return IntStream.range(0,a.length)
      			 .mapToDouble(i->a[i]*b[i])
      			 .sum();
    }
    
    public static double pearsonCorrel(double[] l1, double[] l2){
    	l1=centerVector(l1);
    	l2=centerVector(l2);
    	return dot(l1,l2)/(l2Norm(l1)*l2Norm(l2));
    }
    public static double pearsonCorrel(int[] l1, int[] l2){
    	return pearsonCorrel(Arrays.stream(l1).mapToDouble(i->i).toArray(),
    						 Arrays.stream(l2).mapToDouble(i->i).toArray());
    }
    
    public static double ordinalCorrel(int[] l1){
    	return pearsonCorrel(IntStream.range(0, l1.length).toArray(),
    						 l1);
    }
    
    public static double rankedCorrel(int[] l1){
    	int[] rank=IntStream.range(0, l1.length)
    	         .mapToObj(i->Tuple.of(i,l1[i]).withVComparator())
    	         .sorted()
    	         .mapToInt(t->t.k())
    	         .toArray();
    	return ordinalCorrel(rank);
    }
    
    
}
