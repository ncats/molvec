package tripod.molvec.util;

import java.awt.Point;
import java.awt.Polygon;
import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Optional;
import java.util.function.Predicate;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import tripod.molvec.algo.Tuple;


public class GeomUtil {
    private static final double ZERO_DISTANCE_TOLERANCE = 0.0001;

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

    public static Optional<Point2D> segmentIntersection(Line2D l1, Line2D l2){
    	 Point2D pt = intersection (l1, l2);
         if (pt != null && intersects (l1, pt) && intersects (l2, pt)) {
             return Optional.of(pt);
         }
         return Optional.empty();
    }
    
    /**
     * return the intersection point (if any) between two lines
     */
    public static Point2D intersection (Line2D l1, Line2D l2) {
        Point2D p1 = l1.getP1 (), p2 = l1.getP2 ();
        Point2D p3 = l2.getP1 (), p4 = l2.getP2 ();

        double c = (p1.getX () - p2.getX ()) * (p3.getY () - p4.getY ())
            - (p1.getY () - p2.getY ()) * (p3.getX () - p4.getX ());

        if (Math.abs (c) < ZERO_DISTANCE_TOLERANCE)
            return null;


        double x = (p1.getX () * p2.getY () - p1.getY () * p2.getX ())
            * (p3.getX () - p4.getX ()) - (p1.getX () - p2.getX ())
            * (p3.getX () * p4.getY () - p3.getY () * p4.getX ());
        double y = (p1.getX () * p2.getY () - p1.getY () * p2.getX ())
            * (p3.getY () - p4.getY ()) - (p1.getY () - p2.getY ())
            * (p3.getX () * p4.getY () - p3.getY () * p4.getX ());

        return new Point2D.Double (x / c, y / c);
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
    
    
    
    public static List<List<Shape>> groupShapesIfClosestPointsMatchCriteria(List<Shape> shapes, Predicate<Point2D[]> merge){
    	return groupShapes(shapes,(t)->merge.test(closestPoints(t.k(),t.v())));
    }
    
    public static List<List<Shape>> groupShapes(List<Shape> points, Predicate<Tuple<Shape,Shape>> merge){
    	int[] groups = IntStream.range(0, points.size())
    							.toArray();
    	
    	for(int i=0;i<points.size();i++){
    		Shape p1=points.get(i);
    		for(int j=i+1;j<points.size();j++){
    			Shape p2=points.get(j);
    			if(merge.test(Tuple.of(p1,p2))){
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
    		double ol1=GeomUtil.length(line1);
    		double ol2=GeomUtil.length(line2);
    		double nl=GeomUtil.length(nline);
    		if(nl>ol1 && nl>ol2){
    			return nline;
    		}else{
    			if(ol1>ol2)return line1;
    			return line2;
    		}
    		
    	}
    	
    	public Point2D[] getAbsoluteClosestPointsOnEachLine(){
    		return GeomUtil.closestPointsOnLines(line1, line2);
    	}
    	
    	public double getAbsoluteClosestDistance(){
    		Point2D[] pts=getAbsoluteClosestPointsOnEachLine();
    		return pts[0].distance(pts[1]);
    	}
    	
    	public static LineDistanceCalculator from(Line2D line1, Line2D line2){
    		LineDistanceCalculator ldc = new LineDistanceCalculator();
    		ldc.line1=line1;
    		ldc.line2=line2;
    		return ldc;
    	}
    	
    	
    	
    }
    
    public static boolean intersects (Line2D line, Point2D pt) {
        return line.ptSegDist (pt) < ZERO_DISTANCE_TOLERANCE;
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
    
    public static Optional<Point2D> getIntersection(Shape s1, Line2D l1){
    	
    	
    	return Arrays.stream(lines(s1))
    	      .map(l->segmentIntersection(l,l1))
    	      .filter(p->p.isPresent())
    	      .findFirst()
    	      .map(p->p.get());
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
    
    /**
     * Return two closest points on shape pair
     */
    public static Point2D[] closestPoints (Shape s1, Shape s2) {
        return closestPoints(lines(s1),lines(s2));
    }
    
    /**
     * Return two closest points on shape and line
     */
    public static Point2D[] closestPoints (Shape s1, Line2D line) {
        return closestPoints(lines(s1),new Line2D[]{line});
    }
    
    public static Point2D[] closestPoints (Line2D[] l1s, Line2D[] l2s) {
        
        Point2D[] closest = null;
        double dist = Double.MAX_VALUE;
      
        for(Line2D l1: l1s){
        	for(Line2D l2: l2s){
        		Point2D[] p=closestPointsOnLines(l1,l2);
        		double d=p[0].distance(p[1]);
        		if(d<dist){
        			dist=d;
        			closest=p;
        		}
        	}
        }
        return closest;
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
    
    public static Point2D findCenterOfShape(Shape s){
    	Rectangle2D rect=s.getBounds2D();
    	return new Point2D.Double(rect.getCenterX(),rect.getCenterY());
    }
    
    public static Point2D findCenterMostPoint(List<Point2D> points){
    	double[] avg=new double[]{0,0};
    	points.forEach(p->{
    		avg[0]+=p.getX();
    		avg[1]+=p.getY();
    	});
    	avg[0]=avg[0]/points.size();
    	avg[1]=avg[1]/points.size();
    	Point2D avgPt=new Point2D.Double(avg[0],avg[1]);
    	
    	
    	return findClosestPoint(points,avgPt);
    }
    
    public static Point2D findClosestPoint(List<Point2D> points, Point2D avgPt){
    	return points.stream()
	      .map(p->Tuple.of(p,p.distance(avgPt)).withVComparator())
	      .sorted()
	      .map(t->t.k())
	      .findFirst()
	      .orElse(null);
    }
    
    
    
    public static Tuple<Shape,Double> findClosestShapeTo(List<Shape> shapes, Point2D pt){
    	return shapes.stream()
    	      .map(s->Tuple.of(s,distanceTo(s,pt)).withVComparator())
    	      .sorted()
    	      .findFirst()
    	      .orElse(null);
    }
    
    public static Point2D closestPointOnLine(Line2D l, Point2D p){
    	Point2D pp=GeomUtil.projectPointOntoLine(l,p);
    	double len = GeomUtil.length(l);
    	if(pp.distance(l.getP1())>len){
    		return l.getP2();
    	}else if(pp.distance(l.getP2())>len){
    		return l.getP1();
    	}
    	return pp;
    }
    public static Point2D closestPointOnLines(Line2D[] ls, Point2D p){
    	return Arrays.stream(ls)
    	      .map(l->closestPointOnLine(l,p))
    	      .map(p1->Tuple.of(p1,p.distance(p1)).withVComparator())
    	      .sorted()
    	      .map(t->t.k())
    	      .findFirst()
    	      .orElse(null);
    }
    public static Point2D closestPointOnShape(Shape s, Point2D p){
    	return closestPointOnLines(lines(s),p);
    }
    
    public static Point2D[] closestPointsOnLines(Line2D l1, Line2D l2){
    	Point2D inter=GeomUtil.intersection(l1, l2);
    	if(inter!=null && l1.ptSegDist(inter) < ZERO_DISTANCE_TOLERANCE && l2.ptSegDist(inter) < ZERO_DISTANCE_TOLERANCE){
    		return new Point2D[]{inter,inter};
    	}
    	Point2D[] closestProj1 = new Point2D[2];
    	Point2D[] closestProj2 = new Point2D[2];
    	
    	Point2D l21Closest=closestPointOnLine(l1,l2.getP1());
    	Point2D l22Closest=closestPointOnLine(l1,l2.getP2());
    	Point2D l11Closest=closestPointOnLine(l2,l1.getP1());
    	Point2D l12Closest=closestPointOnLine(l2,l1.getP2());
    	
    	
    	
    	
    	if(l21Closest.distance(l2.getP1()) < l22Closest.distance(l2.getP2())){
    		closestProj1[0]=l21Closest;
    		closestProj1[1]=l2.getP1();
    	}else{
    		closestProj1[0]=l22Closest;
    		closestProj1[1]=l2.getP2();
    	}
    	
    	if(l11Closest.distance(l1.getP1()) < l12Closest.distance(l1.getP2())){
    		closestProj2[0]=l1.getP1();
    		closestProj2[1]=l11Closest;
    	}else{
    		closestProj2[0]=l1.getP2();
    		closestProj2[1]=l12Closest;
    	}
    	
    	if(closestProj1[0].distance(closestProj1[1]) <closestProj2[0].distance(closestProj2[1]) ){
    		return closestProj1;
    	}else{
    		return closestProj2;
    	}
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

	public static List<Line2D> asLines(Collection<Path2D> segments){
			List<Line2D> lines= new ArrayList<Line2D>();
	        for(Path2D p2:segments){
	            PathIterator pi=p2.getPathIterator(null);
	            double[] prevPt=null;
	            while(!pi.isDone()){
	        		
	                double[] coord= new double[2];
	                pi.currentSegment(coord);
	                if(prevPt!=null){
	                    Line2D line = new Line2D.Double
	                        (coord[0], coord[1], prevPt[0], prevPt[1]);
	//                    double lineLength=Math.sqrt
	//                        ((coord[0]-prevPt[0])*(coord[0]-prevPt[0])
	//                         +(coord[1]-prevPt[1])*(coord[1]-prevPt[1]));
	                    lines.add(line);
	                    //System.out.println(lineLength);
	                }
	                prevPt=coord;
	                pi.next();
	            }
	        }
	        return lines;
		}

	public static Path2D fromLine(Line2D line){
		GeneralPath gp = new GeneralPath();
		gp.moveTo(line.getX1(),line.getY1());
		gp.lineTo(line.getX2(),line.getY2());
		return gp;
	}

	public static List<Path2D> fromLines(List<Line2D> lines){
		return lines.stream().map(l->fromLine(l)).collect(Collectors.toList());
	}

	public static double length(Line2D l){
		return l.getP1().distance(l.getP2());
	}

	/**
		 * <p>Currently, this method takes in a set of lines and attempts to group the lines based on being</p>
		 * <ol>
		 * <li>Sufficiently parallel</li>
		 * <li>Sufficiently close together</li>
		 * <li>Sufficiently "overlapping" in projection of little line to big line</li>
		 * <li>Sufficiently "overlapping" in projection of big line to little line</li>  
		 * </ol>
		 * 
		 * 
		 * @param lines
		 * @param maxDeltaTheta
		 * @param maxDeltaOffset
		 * @param minProjectionRatio
		 * @param minLargerProjectionRatio
		 * @return
		 */
		public static List<List<Line2D>> groupMultipleBonds(List<Line2D> lines, double maxDeltaTheta, double maxDeltaOffset, double minProjectionRatio, double minLargerProjectionRatio){
			double[] recipLengths=new double[lines.size()];
			Map<Integer,Integer> groups = new HashMap<>();
			for(int i=0;i<lines.size();i++){
				Line2D line=lines.get(i);
				recipLengths[i]=1.0/line.getP1().distance(line.getP2());
				groups.put(i, i);
			}
			double cosThetaCutoff=Math.cos(maxDeltaTheta);
			
			
			//So, this will be an n^2 operation
			for(int i=0;i<lines.size();i++){
				Line2D line1=lines.get(i);
				double[] vec1=asVector(line1);
				for(int j=i+1;j<lines.size();j++){
					Line2D line2=lines.get(j);
					double[] vec2=asVector(line2);
					double dot=vec1[0]*vec2[0]+vec1[1]*vec2[1];
					double cosTheta=Math.abs(dot*recipLengths[i]*recipLengths[j]);
					//Sufficiently parallel
					if(cosTheta>cosThetaCutoff){
						double[] vec2Off = new double[2];
						vec2Off[0] = (line2.getX1()+line2.getX2())/2 -line1.getX1();
						vec2Off[1] = (line2.getY1()+line2.getY2())/2 -line1.getY1();
						double antiDot=(vec2Off[0]*vec1[1]-vec2Off[1]*vec1[0]);
						double rej=Math.abs(antiDot*recipLengths[i]);
						//Sufficiently close to line
						if(rej<maxDeltaOffset){
							
							double len1=1/recipLengths[i];
							double proj1=0;
							double proj2=0;
							vec2Off[0] = line2.getX1() -line1.getX1();
							vec2Off[1] = line2.getY1() -line1.getY1();
							proj1=(vec2Off[0]*vec1[0]+vec2Off[1]*vec1[1])*recipLengths[i];
							vec2Off[0] = line2.getX2() -line1.getX1();
							vec2Off[1] = line2.getY2() -line1.getY1();
							proj2=(vec2Off[0]*vec1[0]+vec2Off[1]*vec1[1])*recipLengths[i];
							
							if(proj1<0)proj1=0;
							if(proj1>len1){
								proj1=len1;
							}
							if(proj2<0)proj2=0;
							if(proj2>len1){
								proj2=len1;
							}
							double intersectLength=Math.abs(proj1-proj2);
							double cutoffRat1=Math.max(recipLengths[i], recipLengths[j]);
							double cutoffRat2=Math.min(recipLengths[i], recipLengths[j]);
							if(intersectLength*cutoffRat1>minProjectionRatio 
									&& intersectLength*cutoffRat2>minLargerProjectionRatio
									){
	//							System.out.println("WORKED");
	//							System.out.println((intersectLength*cutoffRat1) + " >? " + minProjectionRatio);
	//							System.out.println((intersectLength*cutoffRat2) + " >? " + minLargerProjectionRatio);
								int g1=groups.get(i);
								int g2=groups.get(j);
								groups.put(i, Math.min(g1, g2));
								groups.put(j, Math.min(g1, g2));
							}else{
	//							System.out.println("Didn't work");
	//							System.out.println((intersectLength*cutoffRat1) + " >? " + minProjectionRatio);
	//							System.out.println((intersectLength*cutoffRat2) + " >? " + minLargerProjectionRatio);
							}
						}
					}else{
	//					System.out.println("Not par");
	//					System.out.println((cosTheta) + " >? " + cosThetaCutoff);
					}
				}
			}
			
			return groups.entrySet().stream()
			      .map(Tuple::of)
			      .map(Tuple.kmap(i->lines.get(i)))
			      .map(t->t.swap())
			      .collect(Tuple.toGroupedMap())
			      .entrySet()
			      .stream()
			      .map(Tuple::of)
			      .map(t->t.v())
			      .collect(Collectors.toList());
			
			
		}

	/**
	 * This method, like {@link #groupMultipleBonds(List, double, double, double, double)}, attempts to group bonds, but also 
	 * goes one step further to select the longest line from each group while also returning the count of lines
	 * found in the group. This is an estimate of the bond order.
	 * 
	 * @param lines
	 * @param maxDeltaTheta
	 * @param maxDeltaOffset
	 * @param minProjectionRatio
	 * @param minLargerProjectionRatio
	 * @return
	 */
	public static List<Tuple<Line2D, Integer>> reduceMultiBonds(List<Line2D> lines, double maxDeltaTheta, double maxDeltaOffset, double minProjectionRatio, double minLargerProjectionRatio){
		
		return groupMultipleBonds(lines,maxDeltaTheta,maxDeltaOffset,minProjectionRatio,minLargerProjectionRatio)
				.stream()
				.map(l->Tuple.of(l,l.size()))
				.map(Tuple.kmap(l->l.stream()
								   .map(s->Tuple.of(s,-length(s)).withVComparator()) //always choose longer line
								   .sorted()
								   .findFirst()
						           .get()
						           .k()
						           ))
				.collect(Collectors.toList());
	}

	public static double[] asVector(Line2D line){
		double dx=line.getX2()-line.getX1();
		double dy=line.getY2()-line.getY1();
		return new double[]{dx,dy};
		
	}

	public static Shape getClosestShapeTo(List<Shape> shapes, Point2D pnt){
		return shapes.stream()
				      .map(s->Tuple.of(s,distanceTo(s, pnt)).withVComparator())
				      .sorted()
				      .findFirst()
				      .map(t->t.k())
				      .orElse(null);
	}

	public static double cosTheta(Line2D l1, Line2D l2){
		double[] vec1=asVector(l1);
		double[] vec2=asVector(l2);
		double dot=vec1[0]*vec2[0]+vec1[1]*vec2[1];
		double cosTheta=Math.abs(dot/(length(l1)*length(l2)));
		return cosTheta;
	}

	public static Line2D longestLineFromOneVertexToPoint(Line2D line, Point2D pnt){
		double d1= pnt.distance(line.getP1());
		double d2= pnt.distance(line.getP2());
		if(d1>d2){
			return new Line2D.Double(line.getP1(), pnt);
		}else{
			return new Line2D.Double(line.getP2(), pnt);
		}
	}

	public static ConnectionTable getConnectionTable(List<Tuple<Line2D,Integer>> linest, List<Shape> likelyNodes, double maxDistanceRatioNonLikely,double maxDistanceRatioLikely, double maxDistanceRatioPerLine, Predicate<Line2D> acceptNewLine){
				List<Tuple<Line2D,Integer>> lines = new ArrayList<>(linest);
				
				for(int i=0;i<lines.size();i++){
					Tuple<Line2D, Integer> lineOrder1=lines.get(i);
					double distance1 = GeomUtil.length(lineOrder1.k());
					for(int j=i+1;j<lines.size();j++){
						
						
						Tuple<Line2D, Integer> lineOrder2=lines.get(j);
						double distance2 = GeomUtil.length(lineOrder2.k());
						double totalDistance = distance1+distance2;
						
						Point2D intersect = intersection(lineOrder1.k(),lineOrder2.k());
						if(intersect==null)continue;
						double ndistance1 = Math.max(intersect.distance(lineOrder1.k().getP1()), intersect.distance(lineOrder1.k().getP2()));
						double ndistance2 = Math.max(intersect.distance(lineOrder2.k().getP1()), intersect.distance(lineOrder2.k().getP2()));
						double totalDistanceAfter = ndistance1+ndistance2;
						
						double ratioTotal = Math.max(totalDistanceAfter,totalDistance)/Math.min(totalDistanceAfter, totalDistance);
						
						double ratioLine1 = Math.max(ndistance1,distance1)/Math.min(ndistance1, distance1);
						double ratioLine2 = Math.max(ndistance2,distance2)/Math.min(ndistance2, distance2);
						
						boolean merge = false;
						
					
								
						if(ratioTotal<maxDistanceRatioLikely){
							if(ratioTotal<maxDistanceRatioNonLikely){
								if(ratioLine1<maxDistanceRatioPerLine && ratioLine2<maxDistanceRatioPerLine){
								merge=true;
								}
							}else{
								boolean inLikelyNode=likelyNodes.stream().filter(s->s.contains(intersect)).findAny().isPresent();
								if(inLikelyNode){
									merge=true;
								}
							}
						}
						
						if(merge){
							Line2D newLine1 = longestLineFromOneVertexToPoint(lineOrder1.k(),intersect);
							Line2D newLine2 = longestLineFromOneVertexToPoint(lineOrder2.k(),intersect);
							if(acceptNewLine.test(newLine1) && acceptNewLine.test(newLine2)){
								Tuple<Line2D,Integer> norder1 = Tuple.of(newLine1,lineOrder1.v());
								Tuple<Line2D,Integer> norder2 = Tuple.of(newLine2,lineOrder2.v());
								lines.set(i, norder1);
								lines.set(j, norder2);
								lineOrder1 = norder1;
							}
						}				
					}
				}
					
				ConnectionTable ct=new ConnectionTable();
				
				for(int i=0;i<lines.size();i++){
					Tuple<Line2D, Integer> lineOrder1=lines.get(i);
					ct.addNode(lineOrder1.k().getP1());
					ct.addNode(lineOrder1.k().getP2());
					ct.addEdge(i*2, i*2+1, lineOrder1.v());
				}
				return ct;
				
				
	}

	public static Point2D projectPointOntoLine(Line2D l, Point2D p){
		double[] vec=asVector(l);
		double[] pvec = new double[]{p.getX()-l.getX1(),p.getY()-l.getY1()};
		double dot=dot(vec, pvec);
		double proj=dot/Math.pow(l2Norm(vec),2);
		double[] projVec= new double[]{proj*vec[0],proj*vec[1]};
		return new Point2D.Double(projVec[0]+l.getX1(), projVec[1]+l.getY1());
	}
    
    
}
