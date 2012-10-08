
import java.util.*;
import java.awt.Shape;
import java.awt.Polygon;
import java.awt.Point;
import java.awt.geom.*;

import java.util.logging.Logger;
import java.util.logging.Level;


public class GeomUtil {
    private static final Logger logger = Logger.getLogger(GeomUtil.class.getName());

    private static final boolean DEBUG;
    static {
	boolean debug = false;
	try {
	     debug = Boolean.getBoolean("geomutil.debug");
	}
	catch (Exception ex) {
	}
	DEBUG = debug;
    }

    public static double ccw (Point2D p1, Point2D p2, Point2D p3) {
	return (p2.getX() - p1.getX())*(p3.getY()-p1.getY()) 
            - (p2.getY()-p1.getY())*(p3.getX()-p1.getX());
    }

    public static double angle (Point2D p0, Point2D p1) {
	double dx = p1.getX() - p0.getX(), dy = p1.getY() - p0.getY();
	if (dx > 0 && dy > 0) {
	    return Math.atan(dy/dx);
	}
	else if (dx > 0 && dy == 0) {
	    return 0.;
	}
	else if (dx < 0 && dy > 0) {
	    return Math.PI - Math.atan(-1.*dy/dx);
	}
	else if (dx < 0 && dy == 0) {
	    return Math.PI;
	}
	else if (dx == 0 && dy > 0) {
	    return Math.PI/2;
	}
	else if (dx == 0 && dy < 0) {
	    return 3*Math.PI/2;
	}
	else if (dx < 0 && dy < 0) {
	    return 3*Math.PI/2 - Math.atan(dy/dx);
	}
	else if (dx > 0 && dy < 0) {
	    return 2*Math.PI - Math.atan(-1.*dy/dx);
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
                poly.addPoint((int)(pt.getX()+.5), (int)(pt.getY()+.5));
            }
            return poly;
	}

	Point2D anchor = null;
	for (Point2D pt : pts) {
	    if (anchor == null || pt.getY() < anchor.getY()) {
		anchor = pt;
	    }
	    else if (pt.getY() == anchor.getY() 
                     && pt.getX() < anchor.getX()) {
		anchor = pt;
	    }
	}

	final Point2D p0 = anchor;
	Arrays.sort(pts, new Comparator<Point2D>() {
			public int compare (Point2D a, Point2D b) {
			    double a0 = angle (p0, a), a1 = angle (p0, b);
			    /*
			    System.err.println("p0=("+p0.x+","+p0.y+") a=("+a.x+","+a.y+") "
					       +"b=("+b.x+","+b.y+") ccw="+ccw(p0,a,b)
					       +" theta(p0,a)="+a0+" theta(p0,b)="+a1);
			    */
			    if (a0 < a1) return -1;
			    if (a0 > a1) return 1;

			    double d0 = a.distance(p0), d1 = b.distance(p0);
			    if (d0 < d1) return -1;
			    if (d0 > d1) return 1;
			    return 0;
			}
		    });

	if (DEBUG) {
	    logger.info("Starting point: "+p0);
	    logger.info("Points..."+pts.length);
	    for (int i = 0; i < pts.length; ++i) {
		System.err.println(i+": "+pts[i]);
	    }
	}

	LinkedList<Point2D> stack = new LinkedList<Point2D>();
	stack.push(pts[0]);
	stack.push(pts[1]);
	stack.push(pts[2]);
	for (int i = 3; i < pts.length; ++i) {
	    Point2D pi = pts[i];
	    while (true) {
		Point2D p2 = stack.pop();
		Point2D p1 = stack.peek();
		double dir = ccw (p1, p2, pi);
		if (DEBUG) {
		    System.err.println("p1=("+p1.getX()+","+p1.getY()+") p2=("
				       +p2.getX()+","+p2.getY()+") p"
                                       +i+"=("+pi.getX()+","+pi.getY()
                                       +") ccw="+dir);
		}
		if (dir >= 0.) { // push back
		    stack.push(p2);
		    break;
		}
		if (DEBUG) {
		    logger.info("removing "+p2);
		}
	    }

	    stack.push(pi);
	}

	if (DEBUG) {
	    logger.info("Convex hull: "+stack.size());
	}

	Polygon hull = new Polygon ();
	for (ListIterator<Point2D> it = stack.listIterator(stack.size()); 
             it.hasPrevious(); ) {
	    Point2D pt = it.previous();
	    // should prune collinear points
	    if (DEBUG) {
		System.err.println(" "+pt);
	    }
	    hull.addPoint((int)(pt.getX()+.5), (int)(pt.getY()+.5));
	}

	return hull;
    }    

    public static boolean isNeighbor (Point p1, Point p2) {
	return (Math.abs(p1.x-p2.x) <= 1 && Math.abs(p1.y-p2.y) <= 1);
    }

    public static Point2D[] vertices (Shape shape) {
        return vertices (shape, null);
    }

    public static Polygon toPolygon (Shape shape) {
        Polygon poly = new Polygon ();
        for (Point2D p : vertices (shape)) {
            poly.addPoint((int)(p.getX()+.5), (int)(p.getY()+.5));
        }
        return poly;
    }

    // create a new convex hull shape from two given shapes
    public static Polygon add (Shape s1, Shape s2) {
        ArrayList<Point2D> pts = new ArrayList<Point2D>();
        for (Point2D p : vertices (s1)) {
            pts.add(p);
        }
        for (Point2D p : vertices (s2)) {
            pts.add(p);
        }
        return convexHull (pts.toArray(new Point2D[0]));
    }

    public static Point2D[] vertices (Shape shape, AffineTransform afx) {
        PathIterator p = shape.getPathIterator(afx);
        List<Point2D> vertices = new ArrayList<Point2D>();
        double[] coord = new double[6];
        while (!p.isDone()) {
            int type = p.currentSegment(coord);
            if (type != PathIterator.SEG_CLOSE) {
                vertices.add(new Point2D.Double (coord[0], coord[1]));
            }
            p.next();
        }

        return vertices.toArray(new Point2D[0]);
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
                double dist = distance (pts1[i], pts2[j]);
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
     * Euclidean distance between two points
     */
    public static double distance (Point2D pt1, Point2D pt2) {
        double dx = pt1.getX() - pt2.getX(), dy = pt1.getY() - pt2.getY();
        return Math.sqrt(dx*dx + dy*dy);
    }

    /**
     * Euclidean between two polygons
     */
    public static double distance (Shape s1, Shape s2) {
        Point2D[] vertex = nearestNeighborVertices (s1, s2);
        return distance (vertex[0], vertex[1]);
    }

    /**
     * Given a collection of polygons, this routine performs transitive 
     * closure based on containment and neighborhood distance threshold.
     * Here distance is measured as the distance between the closest
     * vertices between two polygons.
     */
    public static List<Shape> transitiveClosure 
        (Collection<Shape> polygons, double distance) {

        Shape[] poly = polygons.toArray(new Shape[0]);
        List<Shape> closure = new ArrayList<Shape>();
        
        // first pass remove containments and merge intersections
        for (int i = 0; i < poly.length; ++i) {
            for (int j = 0; j < poly.length; ++j) {
                if (i != j && poly[i] != null && poly[j] != null) {
                    Rectangle2D bounds = poly[j].getBounds2D();
                    if (poly[i].contains(bounds)) {
                        poly[j] = null;
                    }
                    else if (poly[i].intersects(bounds)) {
                        poly[i] = add (poly[i], poly[j]);
                        poly[j] = null;
                    }
                }
            }
        }

        for (int i = 0; i < poly.length; ++i) {
            if (poly[i] != null) {
                closure.add(poly[i]);
            }
        }

        return closure;
    }


    static int transitiveClosure (int index, Polygon[] polygons) {
        Polygon poly = polygons[index];
        if (poly == null) {
            return 0;
        }

        Rectangle2D bounds = poly.getBounds2D();
        double cx = bounds.getCenterX();
        double cy = bounds.getCenterY();

        double[] coord = new double[6];
        int k = 0;
        for (int i = 0; i < polygons.length; ++i) {
            if (i == index) continue;
            Polygon p = polygons[i];
            PathIterator path = polygons[i].getPathIterator(null);
            while (!path.isDone()) {
                int type = path.currentSegment(coord);
                if (type != PathIterator.SEG_CLOSE) {
                    double x = coord[0], y = coord[1];
                    if (Math.abs(cx-x) <= bounds.getWidth()
                        && Math.abs(cy-y) <= bounds.getHeight()) {
                        poly.addPoint((int)x, (int)y);
                        ++k;
                    }
                }
                path.next();
            }
        }
        return 0;
    }
}
