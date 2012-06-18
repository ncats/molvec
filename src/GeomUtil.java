
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

    public static int ccw (Point p1, Point p2, Point p3) {
	return (p2.x - p1.x)*(p3.y-p1.y) - (p2.y-p1.y)*(p3.x-p1.x);
    }

    public static double angle (Point p0, Point p1) {
	int dx = p1.x - p0.x, dy = p1.y - p0.y;
	if (dx > 0 && dy > 0) {
	    return Math.atan((double)dy/dx);
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
	    return 3*Math.PI/2 - Math.atan((double)dy/dx);
	}
	else if (dx > 0 && dy < 0) {
	    return 2*Math.PI - Math.atan(-1.*dy/dx);
	}
	return 0.;
    }

    /**
     * Graham scan algorithm for convex hull
     */
    public static Polygon convexHull (Point... pts) {
	if (pts.length < 3) {
	    return new Polygon ();
	}

	Point anchor = null;
	for (Point pt : pts) {
	    if (anchor == null || pt.y < anchor.y) {
		anchor = pt;
	    }
	    else if (pt.y == anchor.y && pt.x < anchor.x) {
		anchor = pt;
	    }
	}

	final Point p0 = anchor;
	Arrays.sort(pts, new Comparator<Point>() {
			public int compare (Point a, Point b) {
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

	LinkedList<Point> stack = new LinkedList<Point>();
	stack.push(pts[0]);
	stack.push(pts[1]);
	stack.push(pts[2]);
	for (int i = 3; i < pts.length; ++i) {
	    Point pi = pts[i];
	    while (true) {
		Point p2 = stack.pop();
		Point p1 = stack.peek();
		int dir = ccw (p1, p2, pi);
		if (DEBUG) {
		    System.err.println("p1=("+p1.x+","+p1.y+") p2=("
				       +p2.x+","+p2.y+") p"+i+"=("+pi.x+","+pi.y+") ccw="+dir);
		}
		if (dir >= 0) { // push back
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
	for (ListIterator<Point> it = stack.listIterator(stack.size()); it.hasPrevious(); ) {
	    Point pt = it.previous();
	    // should prune collinear points
	    if (DEBUG) {
		System.err.println(" "+pt);
	    }
	    hull.addPoint(pt.x, pt.y);
	}

	return hull;
    }    

    public static boolean isNeighbor (Point p1, Point p2) {
	return (Math.abs(p1.x-p2.x) <= 1 && Math.abs(p1.y-p2.y) <= 1);
    }
}
