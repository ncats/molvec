package gov.nih.ncats.molvec.internal.util;

import java.awt.Point;
import java.awt.Polygon;
import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Area;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.NoninvertibleTransformException;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.Stack;
import java.util.concurrent.ThreadLocalRandom;
import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.Supplier;
import java.util.logging.Logger;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import gov.nih.ncats.molvec.internal.algo.Tuple;
import gov.nih.ncats.molvec.internal.util.GeomUtil.ShapeWrapper;


public class GeomUtil {
    private static final double ZERO_DISTANCE_TOLERANCE = 0.0001;
    private static final double ZERO_DISTANCE_OPTIMISTIC_TOLERANCE = 0.000001;

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
     * calculate angle between (x0,y0) and (x1,y1) in radians
     */
    public static double angle (double x0, double y0, double x1, double y1) {
        double dx = x1 - x0, dy = y1 - y0;
        return Math.atan2(dy, dx);
    }

    public static Shape booleanAddShape(Shape s1, Shape s2){
    	Area a1 = new Area(s1);
    	Area a2 = new Area(s2);
    	a1.add(a2);
    	return a1;
    }
    
    public static Shape growShapeHex(Shape s, double rad){
    	return growShapeNPoly(s,rad,6);
    }
    
    public static Shape growShapeNPoly(Shape s, double rad, int n){
    	return Arrays.stream(vertices(s))
    				 .flatMap(v->Arrays.stream(makeNPolyCenteredAt(v,n,rad)))
    				 .collect(convexHull());
    }
    
    public static Point2D[] makeNPolyCenteredAt(Point2D p, int n, double rad){
    	double delta = Math.PI*2.0/n;
    	
    	return IntStream.range(0, n)
    	         .mapToDouble(i->i*delta)
    	         .mapToObj(t->new Point2D.Double(Math.cos(t)*rad, Math.sin(t)*rad))
    	         .map(p1->new Point2D.Double(p.getX()+p1.getX(),p1.getY()+p.getY()))
    	         .toArray(i->new Point2D[i]);
    }
    
    public static ShapeWrapper makeNPolygonOriginCenter(int n, double rad){
    	
    	return ShapeWrapper.of(GeomUtil.convexHull2(makeNPolyCenteredAt(new Point2D.Double(0,0),n,rad)));
    }

    /**
     * Graham scan algorithm for convex hull
     */
    public static Polygon convexHullOldIntPrecision (Point2D... pts) {
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

                             double d0 = a.distanceSq(p0), d1 = b.distanceSq(p0);
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
                //this is really a check to make sure dir >=0
                //but because of floating point math we have to have a delta
                if (dir >= -0.02) { // push back
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
        int px=Integer.MIN_VALUE;
        int py=Integer.MIN_VALUE;
        
        boolean hasDelta=false;
        boolean hasPrevious=false;
        
        int pdx=Integer.MIN_VALUE;
        int pdy=Integer.MIN_VALUE;
        int rej=0;
        boolean error=false;
        
        Stack<int[]> pointsToAdd = new Stack<int[]>();
        
        for (ListIterator<Point2D> it = stack.listIterator (stack.size ());
             it.hasPrevious (); ) {
            Point2D pt = it.previous ();
            // should prune collinear points
            if (DEBUG) {
                System.err.println (" " + pt);
            }
            int nx= (int) (pt.getX () + .5);
            int ny= (int) (pt.getY () + .5);
            
            
            if(hasPrevious && (nx==px && ny==py)){
            	
            }else{

            	int ndx=nx-px;
        		int ndy=ny-py;
        		
            	
            	if(hasDelta){
            		int nrej = (int)Math.signum(ndx*pdy - ndy*pdx);
            		
            		if(nrej==0){
            			pointsToAdd.pop();
            		}else{
	            		if(rej!=0){
	            			if(rej!=nrej){
	            				error=true;
	            			}
	            		}
	            		rej=nrej;
            		}
            	}
           		pointsToAdd.push(new int[]{nx,ny});
           		
       			if(hasPrevious){
            		pdx=ndx;
            		pdy=ndy;
            		hasDelta=true;
            	}
                px=nx;
                py=ny;
                hasPrevious=true;
            }
        }
        for(int[] xy : pointsToAdd){
        	hull.addPoint(xy[0], xy[1]);
        }
    	
        if(error){
        	//TODO: There must be a bug in this code, because this should never get called, but it is called sometimes.
//        	System.err.println("Not hull:" + Arrays.toString(pts));
//        	System.err.println("Not hull:" + Arrays.toString(vertices(hull)));

        	
        	hull= convexHullOldIntPrecision(vertices(hull));
        }
        
        

        return hull;
    }
    
    /**
     * Graham scan algorithm for convex hull
     */
    public static Shape convexHull2 (Point2D... pts) {
    	//if(true)return convexHullOldIntPrecision(pts);
    	
    	AffineTransform at = new AffineTransform();
    	at.scale(5, 5);
    	
    	Point2D[] tpts = Arrays.stream(pts)
    			               .map(p->at.transform(p, null))
    			               .toArray(i->new Point2D[i]);
    	
    	Shape oShape=GeomUtil.convexHullOldIntPrecision(tpts);
    	AffineTransform atInv=null;
		try {
			atInv = at.createInverse();
		} catch (NoninvertibleTransformException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//TODO: right now this actually just uses the integer-precision
		//algorithm, and sclaes it up/down for slightly better precision.
		//The code below this statement does work, in general,
		//but certain aspects of the image prediction expect the slightly
		//worse precision
    	if(true)return atInv.createTransformedShape(oShape);
    	
    	if(pts.length==0){
    		//GeneralPath poly = new GeneralPath ();
    		return convexHullOldIntPrecision(pts);
    	}
    	if (pts.length < 3) {
            GeneralPath poly = new GeneralPath ();
            int i=0;
            
            for (Point2D pt : pts) {
            	if(i==0)poly.moveTo(pt.getX(), pt.getY());
            	else poly.lineTo(pt.getX(), pt.getY());
            	i++;
            }
            poly.closePath();
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

        GeneralPath hull = new GeneralPath ();
        double px=Integer.MIN_VALUE;
        double py=Integer.MIN_VALUE;
        
        boolean hasDelta=false;
        boolean hasPrevious=false;
        
        double pdx=Integer.MIN_VALUE;
        double pdy=Integer.MIN_VALUE;
        int rej=0;
        boolean error=false;
        
        Stack<double[]> pointsToAdd = new Stack<double[]>();
        
        for (ListIterator<Point2D> it = stack.listIterator (stack.size ());
             it.hasPrevious (); ) {
            Point2D pt = it.previous ();
            // should prune collinear points
            if (DEBUG) {
                System.err.println (" " + pt);
            }
            double nx= pt.getX ();
            double ny= pt.getY ();
            
            
            if(hasPrevious && (Math.abs(nx-px)<ZERO_DISTANCE_TOLERANCE && Math.abs(ny-py)<ZERO_DISTANCE_TOLERANCE)){
            	
            }else{

            	double ndx=nx-px;
        		double ndy=ny-py;
        		
            	
            	if(hasDelta){
            		int nrej = (int)Math.signum(ndx*pdy - ndy*pdx);
            		
            		if(nrej==0){
            			pointsToAdd.pop();
            		}else{
	            		if(rej!=0){
	            			if(rej!=nrej){
	            				error=true;
	            			}
	            		}
	            		rej=nrej;
            		}
            	}
           		pointsToAdd.push(new double[]{nx,ny});
           		
       			if(hasPrevious){
            		pdx=ndx;
            		pdy=ndy;
            		hasDelta=true;
            	}
                px=nx;
                py=ny;
                hasPrevious=true;
            }
        }
        int i=0;
        
        for(double[] xy : pointsToAdd){
        	if(i==0){
        		hull.moveTo(xy[0], xy[1]);	
        	}else{
        		hull.lineTo(xy[0], xy[1]);
        	}
        	i++;
        }
        hull.closePath();
    	
        if(error){
        	//TODO: There must be a bug in this code, because this should never get called, but it is called sometimes.
        	
        	System.err.println("Not hull:" + Arrays.toString(pts));
        	System.err.println("Not hull:" + Arrays.toString(vertices(hull)));

        	
        	hull= (GeneralPath)convexHull2(vertices(hull));
        }
        
        

        return hull;
    }
    
    /**
     * Find the center of mass (centroid) of a shape. This is done by using the triangle area method:
     * <pre>
     * X = SUM[(Xi + Xi+1) * (Xi * Yi+1 - Xi+1 * Yi)] / (6 * A)
	 * Y = SUM[(Yi + Yi+1) * (Xi * Yi+1 - Xi+1 * Yi)] / (6 * A)
	 * </pre>
	 * 
	 * Where A is the signed area of the polygon. For shapes with near zero area, the center of the bounding box is returned instead.
	 * 
     * @param s1
     * @return
     */
    public static Point2D centerOfMass(Shape s1){
    	return ShapeWrapper.of(s1).centerOfMass();
    }
    
    
    
    
    /**
     * Finds the center point of a set of lines, weighted by the line lengths
     * @param lines
     * @return
     */
    public static Point2D centerOfMass(List<Line2D> lines){
    	double sumLength = lines.stream().mapToDouble(l->length(l)).sum()*2;
    	
    	Point2D centerPoint = lines.stream()
    							  .flatMap(l->Stream.of(Tuple.of(l.getP1(), length(l)),Tuple.of(l.getP2(), length(l))))
    							  .map(t->Tuple.of(t.k().getX()*t.v(), t.k().getY()*t.v()))
    							  .reduce((t1,t2)->Tuple.of(t1.k() + t2.k(), t1.v()+t2.v()))
    							  .map(Tuple.vmap(d->d/sumLength))
    							  .map(Tuple.kmap(d->d/sumLength))
    							  .map(t->new Point2D.Double(t.k(),t.v()))
    							  .orElse(new Point2D.Double(0, 0));
    	return centerPoint;
    }
    
    
    /**
     * <p>Given a line and a set of points, this will return a new set of points that are transformed
     * by the same transformation that would take the line to be from (0,0)->(1,0). That is,
     * it will scale the line to be of length 1, center it at the origin, and rotate it so that the line
     * points to (1,0). All points transformed by the same scale,rotation and translation to bring them
     * either below or above the line (y-axis), and between 0 and 1 on the x-axis if their projection would
     * be ON the line segment.
     * </p>
     * 
     * <p>
     * This transformation is useful when trying to detect some kind of signal along a line. For example,
     * seeing if the points projected near a line are oscillating, or have some other periodic nature.
     * </p>
     * @param pts
     * @param line
     * @return
     */
    public static List<Point2D> getOffsetsOfPointsOntoLine(List<Point2D> pts, Line2D line){
    	
    	AffineTransform at = GeomUtil.getTransformFromLineToLine(line, new Line2D.Double(0,0,1,0),false);
    	
    	return pts.stream()
    	   .map(p->at.transform(p, null))
    	   .map(p->Tuple.of(p.getX(),p).withKComparator())
    	   .sorted()
    	   .map(t->t.v())
    	   .collect(Collectors.toList());
    }
    
    //TODO: implement
//    public static Line2D getLongestPartitioningLine(Shape s){
//    	//To do this, you actually want to find the line which, when passing through the shape, would
//    	//produce the smallest sum of square residual lengths for each point (vertex or otherwise) onto that line
//    	//This is effectively the same as a PCA.
//    
//    	//The tricky part here is that you need to have this work for all points along a side of a shape,
//        //not just the vertices
//    	
//    	
//    	
//    	
//    	
//    }
    
    
    public static Line2D findMaxSeparationLine(List<Point2D> pts){
    	double a = pts.stream().mapToDouble(p->p.getX()*p.getX()-p.getY()*p.getY()).sum();
    	double b = pts.stream().mapToDouble(p->p.getX()*p.getY()).sum();
    	Function<Double,Double> deriv=(theta)->{
    		double cos=Math.cos(theta);
    		double sin=Math.sin(theta);
    		return cos*sin*a+(sin*sin-cos*cos)*b;
    	};
    	
    	Tuple<Double,Double> dub = findZeroGrid(deriv,Math.PI/2,0,5);
    	
    	double btheta= findZeroInterp(deriv,dub.v(),dub.k(),0.001,0.001);
    	
    	Line2D lp1 = new Line2D.Double(0, 0, Math.cos(btheta), Math.sin(btheta));
    	Line2D lp2 = new Line2D.Double(0, 0, -Math.sin(btheta),Math.cos(btheta));
    	
    	double[] v1= asVector(lp1);
    	double[] v2= asVector(lp2);
    	
    	Tuple<Double,Double> resids=pts.stream()
    	   .map(p->asVector(p))
    	   .map(v->Tuple.of(orthoDot(v1,v), orthoDot(v2,v)))
    	   .map(t->Tuple.of(t.k()*t.k(), t.v()*t.v()))
    	   .reduce((t1,t2)->Tuple.of(t1.k()+t2.k(),t1.v()+t2.v())).orElse(null);
    	
    	if(resids==null)return new Line2D.Double(0,0,1,0);
    	
    	if(resids.k()>resids.v()){
    		Line2D lt=lp2;
    		lp2=lp1;
    		lp1=lt;
    	}
    	
    	return lp1;
    }
    
    public static Line2D findMaxSeparationLinePCA(List<Point2D> pts){
    	
    	double[] xs = new double[pts.size()];
    	double[] ys = new double[pts.size()];
    	for(int i=0;i<pts.size();i++){
    		Point2D p=pts.get(i);
    		xs[i]=p.getX();
    		ys[i]=p.getY();
    	}
    	double[] vec=getPCALikeUnitVector(xs, ys);
    	
    	return new Line2D.Double(0,0,vec[0],vec[1]);
    }
    
    public static LineWrapper findLongestSplittingLine(Shape s){
    	if(s instanceof Line2D)return LineWrapper.of((Line2D)s);
    	Point2D center = GeomUtil.centerOfMass(s);
    	int ptNum = 1000;
    	List<Point2D> pts = 
    			//Arrays.stream(vertices(s))
    			getNEquallySpacedPointsAroundShape(s,ptNum).stream()
    			                .map(p->new Point2D.Double(p.getX()-center.getX(),p.getY()-center.getY()))
    			                .collect(Collectors.toList());
    	
    	Line2D mline=findMaxSeparationLine(pts);
    	
    	Line2D mlineReal = new Line2D.Double(mline.getX1() + center.getX(), mline.getY1()+center.getY(), mline.getX2() + center.getX(), mline.getY2()+center.getY());
    	
    	
    	List<Point2D> ptsIntersect=Stream.of(lines(s))
    	      .map(l->Tuple.of(l,GeomUtil.intersection(mlineReal, l)))
    	      .filter(t->t.v()!=null)
    	      .filter(t->t.k().ptSegDist(t.v())<0.001)
    	      .map(t->t.v())
    	      .collect(Collectors.toList());
    	
    	Point2D[] ptsDist=GeomUtil.getPairOfFarthestPoints(ptsIntersect);
    	
    	return LineWrapper.of(new Line2D.Double(ptsDist[0], ptsDist[1]));
    	
    	
    }
    
    public static List<Point2D> getNEquallySpacedPointsAroundShape(Shape s, int n){
    	double p=getPerimeter(s);
    	
    	//assuming this always gives back in CW or CCW order
    	List<LineWrapper> lines = Arrays.stream(lines(s)).map(l->LineWrapper.of(l)).collect(Collectors.toList());
    	
    	if(lines.size()==0){
    		Rectangle2D rect=s.getBounds2D();
    		Point2D p1=new Point2D.Double(rect.getMinX(),rect.getMinY());
    		Point2D p2=new Point2D.Double(rect.getMaxX(),rect.getMaxY());
    		return Stream.of(p1,p2).collect(Collectors.toList());
    	}
    	
    	List<Point2D> pts = new ArrayList<Point2D>();
    	
    	double mult = p / ((double)n);
    	
    	for(int i=0;i<n;i++){
    		double pc = i*mult;
    		
    		double lenSoFar=0;
    		
    		LineWrapper gline=null;
    		for(int j=0;j<lines.size();j++){
    			double nl = lines.get(j).length();
    			if(nl+lenSoFar>pc){
    				//found the right one!
    				gline=lines.get(j);
    				break;
    			}
    			lenSoFar+=nl;
    			gline=lines.get(j);
    		}
    		
    		
    		
    		double res=pc-lenSoFar;
    		
    		double pres = res*gline.recipLength();
    		
    		double[] v1=gline.vector();
    		double[] start=gline.offset();
    		double[] end =addVectors(negate(start),new double[]{pres*v1[0],pres*v1[1]});
    		pts.add(new Point2D.Double(end[0],end[1]));
    	}
    	
    	return pts;
    	
    }
    
    public static double getPerimeter(Shape s){
    	return Arrays.stream(lines(s))
    	      .mapToDouble(l->length(l))
    	      .sum();
    }
    
    
	 /**
	  * 
	  * Returns the y-axis, x-axis, and product second moment of area for a set
	  * of CCW vertices defining a shape.
	  * 
	  * This is effectively the same as returning the variance in the x-direction,
	  * y-direction, and along y=x.
	  * 
	  * 
	  * For more information:
	  * https://apps.dtic.mil/dtic/tr/fulltext/u2/a183444.pdf
	  * 
	  * @param xs
	  * @param ys
	  * @return
	  */
    public static double[] getSecondMomentXYandCross(double[] xs, double[] ys){
    	
    	double ix = 0;
    	double iy = 0;
    	double ixy = 0;
    	double c1 = 1/12.0;
    	double c2 = 1/24.0;
    	
    	
    	for(int i=0;i<xs.length;i++){
    		int ni = i+1;
    		if(ni>=xs.length)ni=ni%xs.length;
    		double cp = xs[i]*ys[ni] - xs[ni]*ys[i];
    		double xp = xs[i]*xs[i] + xs[i]*xs[ni] + xs[ni]*xs[ni];
    		double yp = ys[i]*ys[i] + ys[i]*ys[ni] + ys[ni]*ys[ni];
    		double xyp = xs[i]*ys[ni] + 2*xs[i]*ys[i]+ 2*xs[ni]*ys[ni] + xs[ni]*ys[i];
    		ix+=cp*xp;
    		iy+=cp*yp;
    		ixy+=cp*xyp;
    	}
    	
    	return new double[]{ix*c1,iy*c1,ixy*c2};
    }
    
    public static double[] getUnitVectorFromVariance(double vx, double vy, double vxy){
    	double cos;
    	
//    	double cos = Math.sqrt(vx)/hyp;
    	if(Math.abs(vxy)<ZERO_DISTANCE_OPTIMISTIC_TOLERANCE){
    		vxy=ZERO_DISTANCE_OPTIMISTIC_TOLERANCE;
    	}
    	
    	double k = (vx-vy)/vxy;
    	double ki = k / (Math.sqrt(k*k+4)); // bounded between -1 and 1
    	
    	
		if (ki <= 1) {
			cos = Math.sqrt((1 - ki) / 2.0);
		}else {
			cos = 0; // if everything else fails, default to +x-direction
		}
		
    	//we now know that the best cos(theta) is either cos1, cos2 or -cos1,-cos2
    	//If we assume we always report sin(theta) as positive    	
    	double sin=Math.sqrt(1-cos*cos);
    	double keepCos = cos;
    	double keepSin = sin;
    	
    	if(vx>vy){
    			//more variance in the x direction, so we want a larger |cos(theta)|
    			keepCos = Math.max(Math.abs(cos),Math.abs(sin));
    			keepSin = Math.min(Math.abs(cos),Math.abs(sin));
    	}else{
    			keepCos = Math.min(Math.abs(cos),Math.abs(sin));
    			keepSin = Math.max(Math.abs(cos),Math.abs(sin));
    	}
    	//aligned in the 1st and 3rd quad, so cos(theta) is positive
    	if(vxy<0){
    		keepCos=-keepCos;
    	}
    	
    	return new double[]{keepCos, keepSin};
    }
    
    public static double[] getPCALikeUnitVector(double[] xs, double[] ys){
    	double vx=0;
    	double vy=0;
    	double vxy=0;  //orientation. >0 means 1st and 3rd quad, <0 means 2nd and 4th quad
    	
    	for(int i=0;i<xs.length;i++){
    		vx+=xs[i]*xs[i];
    		vy+=ys[i]*ys[i];
    		vxy+=xs[i]*ys[i];
    	}
    	
    	return getUnitVectorFromVariance(vx,vy, vxy);
    	    	
    }
    
    
    /**
     * Grid-based search for a zero in a supplied function, using a set number of steps, and returning the first bounds
     * where there is a sign change in the supplied function.  
     * @param f
     * @param high
     * @param low
     * @param steps
     * @return
     */
    private static Tuple<Double,Double> findZeroGrid(Function<Double,Double> f, double high, double low, int steps){
    	double delta = (high-low)/steps;
    	
    	boolean hasLast=false;
    	
    	int psign=1;
    	double px=0;
    	for(int i=0;i<steps;i++){
    		double x=low+delta*i;
    		double y=f.apply(x);
    		int snum=(int)Math.signum(y);
    		if(hasLast){
    			if(snum!=psign){
    				return Tuple.of(px,x);
    			}
    		}
    		
    		psign = snum;
    		px=x;
    		hasLast=true;
    	}
    	return Tuple.of(low,high);
    }
    
    /**
     * Finds the zero of a function by linear interpolation. That is, a high and low x value is specified, and the midpoint
     * x value is also sampled. I
     * @param f
     * @param high
     * @param low
     * @param tolError
     * @param minStep
     * @return
     */
    private static double findZeroInterp(Function<Double,Double> f, double high, double low, double tolError, double minStep){
    	
    	double hv=f.apply(high);
    	double lv=f.apply(low);
    	
    	if(Math.abs(hv)< tolError){
    		return high;
    	}else if(Math.abs(lv)< tolError){
    		return low;
    	}
    	
    	double mid = (high+low)/2;
    	
    	double mv=f.apply(mid);
    	
    	boolean hpos=(hv>0);
    	boolean lpos=(lv>0);
    	
    	if(hpos && !lpos){
    		if(mv>0){
    			high=mid;
    		}else{
    			low=mid;
    		}
    	}else if(!hpos && lpos){
    		if(mv>0){
    			low=mid;
    		}else{
    			high=mid;
    		}
    	}
    	
    	
    	if(Math.abs(high-low)<minStep){
    		return (high+low)/2;
    	}
    	return findZeroInterp(f,high,low,tolError,minStep);
    	
    }
    
    
    public static Tuple<Point2D,double[]> getCircumscribedAndInscribedCircles(ShapeWrapper s1){
    	Point2D centerOfMass=s1.centerOfMass();
    	

    	Line2D[] lines=s1.getLines();
    	
    	double min=Arrays.stream(lines)
    			.map(l->GeomUtil.projectPointOntoLine(l, centerOfMass))
			    	      .mapToDouble(lp->lp.distanceSq(centerOfMass))
			    	      .min()
			    	      .orElse(0);
    	double max=Arrays.stream(s1.getVerts())
	    	      .mapToDouble(v->v.distanceSq(centerOfMass))
	    	      .max()
	    	      .orElse(0);
    	
    	return Tuple.of(centerOfMass,new double[]{Math.sqrt(min),Math.sqrt(max)});
    }
    
    public static double getCircleLikeScore(Shape s){
    	Tuple<Point2D,double[]> circles=getCircumscribedAndInscribedCircles(ShapeWrapper.of(s));
    	double[] rads=circles.v();
    	return rads[0]*rads[0]/(rads[1]*rads[1]);
    }
    

    public static boolean isNeighbor (Point p1, Point p2) {
        return (Math.abs (p1.x - p2.x) <= 1 && Math.abs (p1.y - p2.y) <= 1);
    }

   

    // create a new convex hull shape from two given shapes
    public static Shape add (Shape s1, Shape s2) {
        ArrayList<Point2D> pts = new ArrayList<Point2D> ();
        for (Point2D p : vertices (s1)) {
            pts.add (p);
        }
        for (Point2D p : vertices (s2)) {
            pts.add (p);
        }
        return convexHull2(pts.toArray (new Point2D[0]));
    }
    
    public static ShapeWrapper add(ShapeWrapper s1, ShapeWrapper s2){
    	ArrayList<Point2D> pts = new ArrayList<Point2D> ();
        for (Point2D p : s1.getVerts()) {
            pts.add (p);
        }
        for (Point2D p : s2.getVerts()) {
            pts.add (p);
        }
        return ShapeWrapper.of(convexHull2(pts.toArray (new Point2D[0])));
    }

    public static Point2D[] vertices (Shape shape) {
        return vertices (shape, null);
    }
    
    public static Point2D[] vertices(List<Line2D> lines){
    	return lines.stream().flatMap(l->Stream.of(l.getP1(),l.getP2())).toArray(i->new Point2D[i]);
    }
    
    

    public static Point2D[] vertices (Shape shape, AffineTransform afx) {
    	if(shape ==null){
    		return new Point2D[0];
    	}
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
    
    public static class ShapeWrapper{
    	private Shape s;
    	private Rectangle2D r;
    	private Line2D[] lines;
    	private Point2D[] verts;
    	private Double signedArea = null;
    	
    	private Double inscribedRadius = null;
    	private Double circumscribedRadius = null;
    	
    	private Point2D[] extremes = null;
    	
    	
    	private LineWrapper longestSplittingLine = null;
    	
    	public ShapeWrapper(Shape s){
    		this.s=s;
    	}
    	
    	public Shape getShape(){
    		return this.s;
    	}
    	
    	
    	public LineWrapper findLongestSplittingLine(){
    		if(longestSplittingLine==null){
	        	if(s instanceof Line2D)return LineWrapper.of((Line2D)s);
	        	
	        	Point2D[] verts = this.getVerts();
	        	if(verts.length==2){
	        		longestSplittingLine =  LineWrapper.of(new Line2D.Double(verts[0], verts[1]));
	        		return longestSplittingLine;
	        	}
	        	
	        	Point2D center = centerOfMass();
	        	
	        	double[] xs = new double[verts.length];
	        	double[] ys = new double[verts.length];
	        	for(int i=0;i<verts.length;i++){
	        		xs[i]=verts[i].getX()-center.getX();
	        		ys[i]=verts[i].getY()-center.getY();
	        	}
	        	
	        	double[] stats=getSecondMomentXYandCross(xs,ys);
	        	
	        	double[] vec =getUnitVectorFromVariance(stats[0],stats[1],stats[2]);
	        	
	        	Line2D mline = new Line2D.Double(0,0,vec[0],vec[1]);
	        	
	        	Line2D mlineReal = new Line2D.Double(mline.getX1() + center.getX(), mline.getY1()+center.getY(), mline.getX2() + center.getX(), mline.getY2()+center.getY());
	        	
	        	
	        	List<Point2D> ptsIntersect=Stream.of(getLines())
	        	      .map(l->Tuple.of(l,GeomUtil.intersection(mlineReal, l)))
	        	      .filter(t->t.v()!=null)
	        	      .filter(t->t.k().ptSegDist(t.v())<0.001)
	        	      .map(t->t.v())
	        	      .collect(Collectors.toList());
	        	
	        	if(ptsIntersect.size()<2){
	        		Point2D[] pfpts = this.getPairOfFarthestPoints();
	        		longestSplittingLine= LineWrapper.of(new Line2D.Double(pfpts[0], pfpts[1]));
	        	}else{
		        	Point2D[] ptsDist=GeomUtil.getPairOfFarthestPoints(ptsIntersect);
		        	longestSplittingLine= LineWrapper.of(new Line2D.Double(ptsDist[0], ptsDist[1]));
		        	if(longestSplittingLine.length()<2){
		        		Point2D[] pfpts = this.getPairOfFarthestPoints();
		        		longestSplittingLine= LineWrapper.of(new Line2D.Double(pfpts[0], pfpts[1]));
		        	}
	        	}
	        	
    		}
    		return longestSplittingLine;
        	
        	
        }
    	public double distanceTo(Point2D pt){
    		if(s.contains(pt)){
    	    	return 0;
    	    }
    		
    	    return Math.sqrt(Arrays.stream(getLines())
    	    	      .mapToDouble(l->l.ptSegDistSq(pt))
    	    	      .min()
    	    	      .getAsDouble());    	    
    	}
    	
    	public Point2D[] closestPointsTo(ShapeWrapper sw){
    		return GeomUtil.closestPoints(getLines(), sw.getLines());
    	}
    	
    	public Point2D[] getPairOfFarthestPoints(){
    		if(extremes==null)extremes=GeomUtil.getPairOfFarthestPoints(getVerts());
    		return extremes;
    	}
    	
    	private void calcRadii(){
    		Tuple<Point2D, double[]> circles = getCircumscribedAndInscribedCircles(this);
    		this.circumscribedRadius = circles.v()[1];
    		this.inscribedRadius = circles.v()[0];
    	}
    	
    	public double getInscribedRadius(){
    		if(inscribedRadius==null)calcRadii();
    		return inscribedRadius;
    	}
    	public double getCircumscribedRadius(){
    		if(circumscribedRadius==null)calcRadii();
    		return circumscribedRadius;
    	}
    	
    	public double getCircleLikeScore(){
    		double cr=getCircumscribedRadius();
    		double ir = getInscribedRadius();
    		return ir*ir/(cr*cr);
    	}
    	
    	private void cacheAll(){
    		this.r=s.getBounds2D();
    		this.lines = GeomUtil.lines(s);
    		this.verts=vertices(s);
    		this.signedArea = areaVerticesCW(verts);    		
    	}
    	
    	
    	public List<Line2D> getLinesNotInside(Line2D l){

    		boolean p1Inside = contains(l.getP1());
    		boolean p2Inside = contains(l.getP2());
    		
    		//empty
    		if(p1Inside && p2Inside)return new ArrayList<Line2D>();
    		
    		
    		
    		List<Point2D> ips=getAllIntersections(l);
    		if(ips.isEmpty())return Arrays.asList(l);
    		if(ips.size()>2){
//    			System.out.println("Line:" + l);
//    			System.out.println("Shape:" + Arrays.toString(vertices(s)));
//    			System.out.println("Shape:" + Arrays.toString(vertices(convexHull2(vertices(s)))));
    			throw new IllegalStateException("Line should not intersect with convex hull more than twice, maybe the shape isn't a convux hull?");
    		}
    		
    		
    		Point2D ip1 = ips.get(0);
    		Point2D ip2 = (ips.size()==2)?ips.get(1):null;
    		
    		if(p1Inside && ip2==null){			
    			Line2D ln = new Line2D.Double(ip1,l.getP2());
    			return Arrays.asList(ln);
    		}else if(p2Inside && ip2==null){
    			Line2D ln = new Line2D.Double(l.getP1(),ip1);
    			return Arrays.asList(ln);
    		}else{
    			if(ip2==null){
    				//Could be tangent to one point actually, in this case just return whole line
    				return Arrays.asList(l);
    				
    				//throw new IllegalStateException("Line should not have both points inside convux hull and have no intersections. Something is wrong.");
    			}else{
    				double p1Distance1 = ip1.distanceSq(l.getP1());
    				double p1Distance2 = ip2.distanceSq(l.getP1());
    				
    				Point2D kp1 = l.getP1();
    				Point2D kp2 = l.getP2();
    				
    				if(p1Distance1>p1Distance2){
    					kp1 = l.getP2();
    					kp2 = l.getP1();
    				}
    				return Arrays.asList(new Line2D.Double(kp1,ip1),new Line2D.Double(kp2,ip2));
    			}
    		}
    	}
    	
    	public Line2D[] getLines(){
    		if(lines==null)lines = GeomUtil.lines(s);
    		return lines;
    	}
    	
    	public List<Point2D> getAllIntersections(Line2D l1){
    		   	List<Point2D> plist= Arrays.stream(getLines())
    								    	      .map(l->segmentIntersection(l,l1))
    								    	      .filter(p->p.isPresent())
    								    	      .map(p->p.get())
    								    	      .collect(Collectors.toList());
    		    	
		    	return GeomUtil.groupPointsCloserThan(plist, ZERO_DISTANCE_TOLERANCE)
		    	        .stream()
		    	        .map(pl->pl.get(0))
		    	        .collect(Collectors.toList());
    	}
    	
    	public Point2D[] getVerts(){
    		if(verts==null)verts = GeomUtil.vertices(s);
    		return verts;
    	}
    	
    	public double getWidth(){
    		return getBounds().getWidth();
    	}
    	public double getHeight(){
    		return getBounds().getHeight();
    	}
    	
    	public Rectangle2D getBounds(){
    		if(r==null)r=s.getBounds2D();
    		return r;
    	}
    	
    	public double getSignedArea(){
    		if(signedArea==null)signedArea=areaVerticesCW(getVerts());
    		return signedArea;
    	}
    	
    	public double getArea(){
    		return Math.abs(getSignedArea());
    	}
    	
    	public static ShapeWrapper of(Shape s){
    		return new ShapeWrapper(s);
    	}

		public boolean contains(Point2D p) {
			return s.contains(p);
		}
		
	    /**
	     * Find the center of mass (centroid) of a shape. This is done by using the triangle area method:
	     * <pre>
	     * X = SUM[(Xi + Xi+1) * (Xi * Yi+1 - Xi+1 * Yi)] / (6 * A)
		 * Y = SUM[(Yi + Yi+1) * (Xi * Yi+1 - Xi+1 * Yi)] / (6 * A)
		 * </pre>
		 * 
		 * Where A is the signed area of the polygon. For shapes with near zero area, the center of the bounding box is returned instead.
		 * 
	     * @param s1
	     * @return
	     */
		public Point2D centerOfMass(){

	    	Point2D[] vts = this.getVerts();
	    	
	    	double sarea = this.getSignedArea();
	    	
	    	if(Math.abs(sarea)<ZERO_DISTANCE_TOLERANCE){
	    		return centerOfBounds();
	    	}
	    	
	    	double sumX=0;
	    	double sumY=0;
	    	double rArea = 1/(6*sarea);
	    	
	    	for(int i=0;i<vts.length;i++){
	    		Point2D p1=vts[i];
				Point2D p2=vts[(i+1+vts.length)%vts.length];
				
	    		double w=(p1.getX()*p2.getY()-p2.getX()*p1.getY());
	    		sumX+= (p1.getX()+p2.getX())*w;
	    		sumY+= (p1.getY()+p2.getY())*w;    		
	    	}
	    	
	    	return new Point2D.Double(sumX*rArea, sumY*rArea);
		}
		
		
		public Point2D centerOfBounds(){
			return new Point2D.Double(getBounds().getCenterX(),getBounds().getCenterY());
		}
		
		/**
		 * This method scales the shape based on its bounding box. The value given will scale the shape
		 * such that the new bounding box width is 2*dr larger than it was before. Note that this
		 * may scale a shape quite differently depending on whether they are "tall" or "wide"
		 * @param s
		 * @param dr
		 * @return
		 */
		public ShapeWrapper growShapeBounds(double dr){
			AffineTransform at = new AffineTransform();
			Rectangle2D rect = getBounds();
			
			double scale = (rect.getWidth()+2*dr)/rect.getWidth();
			
			Point2D pt=centerOfBounds();
			at.translate(pt.getX(), pt.getY());
			at.scale(scale, scale);
			at.translate(-pt.getX(), -pt.getY());
			return ShapeWrapper.of(at.createTransformedShape(s));			
		}
		
		public boolean intersects(ShapeWrapper other){
			return GeomUtil.intersects(this, other);
		}
		
		public ShapeWrapper and(ShapeWrapper s2){
			return add(this, s2);
		}
		
		public double distanceSq(ShapeWrapper sother){
			Point2D[] pts = GeomUtil.nearestNeighborVertices(this, sother);
			return pts[0].distanceSq(pts[1]);
		}
		
		public ShapeWrapper growShapeNPoly(double r,int n){
			return ShapeWrapper.of(GeomUtil.growShapeNPoly(getShape(), r, n));
		}
		
		public Optional<Line2D> getLineInside(Line2D l){
			boolean p1Inside = s.contains(l.getP1());
			boolean p2Inside = s.contains(l.getP2());
			
			if(p1Inside && p2Inside)return Optional.of(l);
			
			
			List<Point2D> ips=this.getAllIntersections(l);
			
			if(ips.isEmpty())return Optional.empty();
			
			if(ips.size()>2){
//				System.out.println("Line:" + l);
//				System.out.println("Shape:" + Arrays.toString(vertices(s)));
//				System.out.println("Shape:" + Arrays.toString(vertices(convexHull2(vertices(s)))));
				throw new IllegalStateException("Line should not intersect with convex hull more than twice, maybe the shape isn't a convux hull?");
			}

			Point2D ip1 = ips.get(0);
			Point2D ip2 = (ips.size()==2)?ips.get(1):null;
			
			if(p1Inside && ip2==null){			
				Line2D ln = new Line2D.Double(l.getP1(),ip1);
				return Optional.of(ln);
			}else if(p2Inside && ip2==null){
				Line2D ln = new Line2D.Double(ip1,l.getP2());
				return Optional.of(ln);
			}else{
				if(ip2==null){
					//Could be tangent to one point actually, in this case just return nothing
					return Optional.empty();
				}else{
					return Optional.of(new Line2D.Double(ip1,ip2));
				}
			}
		}
		
		public boolean contains(ShapeWrapper sother){
			return GeomUtil.contains(this.s, sother.s);
		}

		public ShapeWrapper getTransformed(AffineTransform at) {
			Shape ts=at.createTransformedShape(this.getShape());
			return ShapeWrapper.of(ts);
		}

		/**
		 * Return a {@link ShapeWrapper} scaled, translated and rotated so that
		 * the center-of-mass is at the origin and the point farthest 
		 * from the center of mass is at point (1.0,0.0).
		 * @return
		 */
		public ShapeWrapper normalize() {
			Point2D cent=this.centerOfMass();
			Point2D far = Arrays.stream(this.getVerts())
					.map(p->Tuple.of(p.distanceSq(cent), p).withKComparator())
					.max(Comparator.naturalOrder())
					.map(t->t.v())
					.get();
			Line2D l = new Line2D.Double(cent, far);
			Line2D unitLine = new Line2D.Double(0,0,1,0);
			
			AffineTransform aft=GeomUtil.getTransformFromLineToLine(l, unitLine, this.getSignedArea()<0);
			
			return this.getTransformed(aft);
		}
		
		public Optional<ShapeWrapper> getIntersection(ShapeWrapper s2){
			return getIntersectionShape(this,s2);
		}
		
		public double similarity(ShapeWrapper s2){
			
			double iarea=this.getIntersection(s2)
			    .map(is->is.getArea())
			    .orElse(0.0);
			
			double tarea=this.and(s2).getArea();
			
			return iarea/tarea;
			
		}
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

        return lines.stream().filter(l->!l.getP1().equals(l.getP2())).toArray(i->new Line2D[i]);
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
    
    
    
    public static List<List<Shape>> groupShapesIfClosestPointsMatchCriteria(Collection<Shape> shapes, Predicate<Tuple<Shape[],Point2D[]>> merge){
    	return groupShapes(shapes,(t)->{
    		Point2D[] far=closestPoints(t.k(),t.v());
    		Shape[] s= new Shape[]{t.k(),t.v()};
    		return merge.test(Tuple.of(s,far));	
    	});
    }
    
    public static List<List<ShapeWrapper>> groupShapesIfClosestPointsMatchCriteriaSW(Collection<ShapeWrapper> shapes, Predicate<Tuple<ShapeWrapper[],Point2D[]>> merge){
    	return groupThings(shapes,(t)->{
    		Point2D[] far=t.k().closestPointsTo(t.v());
    		ShapeWrapper[] s= new ShapeWrapper[]{t.k(),t.v()};
    		return merge.test(Tuple.of(s,far));	
    	});
    }
    
    public static List<List<Shape>> groupShapes(Collection<Shape> points, Predicate<Tuple<Shape,Shape>> merge){
    	return groupThings(points,merge);
    }
    
    public static <T> List<List<T>> groupThings(Collection<T> things, Predicate<Tuple<T,T>> merge){
    	int[] groups = IntStream.range(0, things.size())
    							.toArray();
    	BitSet bs = new BitSet(things.size());
    	
    	List<T> asList = (things instanceof List)?(List<T>)things:things.stream().collect(Collectors.toList());
    	
    	for(int i=0;i<things.size();i++){
    		T p1=asList.get(i);
    		for(int j=i+1;j<things.size();j++){
    			T p2=asList.get(j);
    			if(merge.test(Tuple.of(p1,p2))){
    				int g1= groups[i];
    				int g2= groups[j];
    				int ng=Math.min(g1, g2);
    				int rg=Math.max(g1, g2);
    				groups[i]=ng;
    				groups[j]=ng;
    				bs.set(ng);
    				if(bs.get(rg)){
    					for(int k=0;k<groups.length;k++){
    						if(groups[k]==rg){
    							groups[k]=ng;
    						}
    					}
    				}
    			}
    		}
    	}
    	
    	
    	return IntStream.range(0, things.size())
    	         .mapToObj(i->Tuple.of(groups[i],i))
    	         .map(Tuple.vmap(i->asList.get(i)))
    	         .collect(Tuple.toLinkedGroupedMap())
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
        return intersects(ShapeWrapper.of(s1), ShapeWrapper.of(s2));
    }
    
    public static boolean intersects (ShapeWrapper s1, ShapeWrapper s2) {
        if(!s1.getBounds().intersects(s2.getBounds()))return false;
    	Line2D[] lines1 = s1.getLines();
        Line2D[] lines2 = s2.getLines();
        for (Line2D l1 : lines1) {
            for (Line2D l2 : lines2) {
                Point2D pt = intersection (l1, l2);
                if (pt != null && intersects (l1, pt) && intersects (l2, pt)) {
                    return true;
                }
            }
        }
        Point2D[] p1s=s1.getVerts();
        Point2D[] p2s=s2.getVerts();
        
        for(Point2D p:p1s){
        	if(s2.contains(p))return true;
        	break;
        }
        for(Point2D p:p2s){
        	if(s1.contains(p))return true;
        	break;
        }
        return false;
    }
    
    public static Optional<Point2D> getIntersection(Shape s1, Line2D l1){
    	
    	return getAllIntersections(s1,l1).stream().findFirst();
    }
    
    public static List<Point2D> getAllIntersections(Shape s1, Line2D l1){
    	return ShapeWrapper.of(s1).getAllIntersections(l1);
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
       
        return nearestNeighborVertices(ShapeWrapper.of(s1),ShapeWrapper.of(s2));
    }
    
    
    public static Point2D[] nearestNeighborVertices (ShapeWrapper s1, ShapeWrapper s2) {
        Point2D[] pts1 = s1.getVerts();
        Point2D[] pts2 = s2.getVerts();

        double minDist = Double.MAX_VALUE;
        Point2D p1 = null, p2 = null;
        for (int i = 0; i < pts1.length; ++i) {
            for (int j = 0; j < pts2.length; ++j) {
                double dist =pts1[i].distanceSq(pts2[j]);
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
        		double d=p[0].distanceSq(p[1]);
        		if(d<dist){
        			dist=d;
        			closest=p;
        		}
        	}
        }
        return closest;
    }
    
    public static List<Point2D> intersectingPoints(List<Line2D> l1s, List<Line2D> l2s){
    	return l1s.stream()
    			  .flatMap(l->l2s.stream().map(l2->segmentIntersection(l2,l)).filter(ip->ip.isPresent()))
    			  .map(op->op.get())
    			  .collect(Collectors.toList());
    }
    public static List<Point2D> intersectingPoints(Shape l1s, Shape l2s){
    	return intersectingPoints(Arrays.asList(lines(l1s)),Arrays.asList(lines(l2s)));
    }
    
    public static List<Point2D> intersectingPoints(ShapeWrapper l1s, ShapeWrapper l2s){
    	return intersectingPoints(Arrays.asList(l1s.getLines()),Arrays.asList(l2s.getLines()));
    }
    
    
    public static double distanceTo(Shape s, Point2D pt){
    	if(s.contains(pt)){
    		return 0;
    	}
    	return Math.sqrt(Arrays.stream(lines(s))
    	      .mapToDouble(l->l.ptSegDistSq(pt))
    	      .min()
    	      .getAsDouble());
    }
    
    public static Point2D findCenterOfShape(Shape s){
    	Rectangle2D rect=s.getBounds2D();
    	return new Point2D.Double(rect.getCenterX(),rect.getCenterY());
    }
    
    public static Point2D findCenterOfVertices(List<Point2D> points){
    	double[] avg=new double[]{0,0};
    	points.forEach(p->{
    		avg[0]+=p.getX();
    		avg[1]+=p.getY();
    	});
    	avg[0]=avg[0]/points.size();
    	avg[1]=avg[1]/points.size();
    	
    	return new Point2D.Double(avg[0],avg[1]);
    }
    
    
    
   
    
    public static Point2D findCenterMostPoint(List<Point2D> points){
    	Point2D avgPt=findCenterOfVertices(points);
    	return findClosestPoint(points,avgPt);
    }
    
    public static Point2D findClosestPoint(List<Point2D> points, Point2D avgPt){
    	return points.stream()
	      .map(p->Tuple.of(p,p.distance(avgPt)).withVComparator())
	      .min(CompareUtil.naturalOrder())
	      .map(t->t.k())
	      .orElse(null);
    }
    
    
    
    public static Optional<Tuple<Shape,Double>> findClosestShapeTo(Collection<Shape> shapes, Point2D pt){
    	return shapes.stream()
    	      .map(s->Tuple.of(s,distanceTo(s,pt)).withVComparator())
    	      .min(CompareUtil.naturalOrder());
    }
    
    public static Optional<Tuple<ShapeWrapper,Double>> findClosestShapeWTo(Collection<ShapeWrapper> shapes, Point2D pt){
    	return shapes.stream()
    	      .map(s->Tuple.of(s,s.distanceTo(pt)).withVComparator())
    	      .min(CompareUtil.naturalOrder());
    }
    
    public static Shape growLine(Line2D l, double width){
    	double[] v= normalize(asVector(l));
    	double wc=width/2;
    	Point2D p1 = l.getP1();
    	Point2D p2 = l.getP2();
    	
    	Point2D p11 = new Point2D.Double(p1.getX()-v[1]*wc,p1.getY()+v[0]*wc);
    	Point2D p12 = new Point2D.Double(p1.getX()+v[1]*wc,p1.getY()-v[0]*wc);
    	Point2D p21 = new Point2D.Double(p2.getX()-v[1]*wc,p2.getY()+v[0]*wc);
    	Point2D p22 = new Point2D.Double(p2.getX()+v[1]*wc,p2.getY()-v[0]*wc);
    	
    	return convexHull2(p11,p12,p21,p22);
    }
    
    public static double[] normalize(double[] v){
    	double lr=1/l2Norm(v);
    	return new double[]{v[0]*lr,v[1]*lr};
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
    	      .min(CompareUtil.naturalOrder())
    	      .map(t->t.k())
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
        return distance(ShapeWrapper.of(s1), ShapeWrapper.of(s2));
    }
    
    /**
     * Euclidean between two polygons based on nearest vertices distance
     */
    public static double distance (ShapeWrapper s1, ShapeWrapper s2) {
        Point2D[] vertex = nearestNeighborVertices (s1, s2);
        return length (vertex[0], vertex[1]);
    }
    
    
    
    /**
     * Euclidean Squared distance between two polygons based on nearest vertices distance
     */
    public static double distanceSq (Shape s1, Shape s2) {
        Point2D[] vertex = nearestNeighborVertices (s1, s2);
        return vertex[0].distanceSq(vertex[1]);
    }
    
    public static double[] centerVector(double[] v){
    	double a1=(IntStream.range(0,v.length)
   			 .mapToDouble(i->v[i])
   			 .average()
   			 .orElse(0));
    	
    	
    	return IntStream.range(0,v.length)
			 .mapToDouble(i->v[i]-a1)
			 .toArray();
    }
    
    public static double l2Norm(double[] v){
    	if(v.length==2)return Math.sqrt(v[0]*v[0]+v[1]*v[1]);
    	return Math.sqrt(IntStream.range(0,v.length)
   			 .mapToDouble(i->v[i])
   			 .map(d->d*d)
   			 .sum());
    }
    public static double dot(double[] a, double[] b){
    	if(a.length==2)return a[0]*b[0]+a[1]*b[1];
    	return IntStream.range(0,a.length)
      			 .mapToDouble(i->a[i]*b[i])
      			 .sum();
    }
    
    public static double orthoDot(double[] a, double[] b){
    	if(a.length!=2 || b.length!=2)throw new IllegalStateException("Cannot do an orthoganal dot product on vectors longer than 2");
    	return dot(a, new double[]{b[1],-b[0]});
    }
    
    /**
     * Get the scalar projection of B onto A
     * @param a
     * @param b
     * @return
     */
    public static double projection(double[] a, double[] b){
    	return dot(a,b)/(l2Norm(a));
    }
    
    /**
     * Get the scalar rejection of B onto A
     * @param a
     * @param b
     * @return
     */
    public static double rejection(double[] a, double[] b){
    	return orthoDot(a,b)/(l2Norm(a));
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
    
    public static double variance(double[] l1){
    	double[] cl=centerVector(l1);
    	return Arrays.stream(cl)
    			     .map(d->d*d)
    			     .average()
    			     .orElse(0);
    }
    public static double variance(int[] l1){
    	return variance(Arrays.stream(l1).mapToDouble(i->i).toArray());
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
	public static double lengthSquared(Line2D l){
		return l.getP1().distanceSq(l.getP2());
	}
	
	public static Predicate<Line2D> longerThan(double d){
		double sq=d*d;
		return (l)->lengthSquared(l)>sq;
	}


	public static class LineWrapper implements Comparable<LineWrapper>{
		private Line2D line;
		private double[] vec = null;
		private double len = 0;
		private double rlen = 0;
		private double[] cvec= null;
		private double[] offset = null;
		private Point2D center=null;
		
		private boolean process=false;
		
		public Line2D getLine(){
			return this.line;
		}
		
		public Point2D centerPoint(){
			process();
			return center;
		}
		
		public List<Point2D> splitIntoNPieces(int n){

			double[] v=this.vector();
			double[] off=this.offset();
			double rlen=1.0/(n);
			return IntStream.range(0,n+1)
					 .mapToObj(i->new double[]{i*v[0]*rlen-off[0], i*v[1]*rlen-off[1]})
					 .map(d->new Point2D.Double(d[0],d[1]))
					 .collect(Collectors.toList());
		}
		
		public Stream<Point2D> streamPoints(){
			return Stream.of(line.getP1(),line.getP2());
		}
		
		public Shape growLine(double d){
			return GeomUtil.growLine(this.line,d);
		}
		
		public LineWrapper process(){
			if(process){
			    return this;
            }
			vec=asVector(line);
			len=l2Norm(vec);
			center=GeomUtil.findCenterOfShape(line);
			rlen=1/len;
			cvec = asVector(findCenterOfShape(line));
			offset= negate(asVector(line.getP1()));
			process=true;
			return this;
		}
		
		public double cosTheta(LineWrapper lw2){
			double d=dot(this.vector(), lw2.vector());
			return d*this.rlen*lw2.rlen;
		}
		
		public double absCosTheta(LineWrapper lw2){
			return Math.abs(cosTheta(lw2));
		}
		
		
		public double centerRejection(LineWrapper lw2){
			double[] l1center=this.centerVector();
			double[] offset = this.offset();
			double[] l2center=lw2.centerVector();
			double[] sVec=addVectors(l1center,offset);
			double[] rVec=addVectors(l2center,offset);
			return rejection(sVec,rVec);
		}
		
		public double projectionOffset(Point2D p){
			double[] pvec = asVector(p);
			double[] nvec = addVectors(pvec,this.offset());
			return dot(vector(),nvec)*this.recipLength();
		}
		
		public double rejectionOffset(Point2D p){
			double[] pvec = asVector(p);
			double[] nvec = addVectors(pvec,this.offset());
			return orthoDot(vector(),nvec)*this.recipLength();
		}
		
		public double[] vector(){
			process();
			return this.vec;
		}
		public double length(){
			process();
			return this.len;
		}
		public double lengthSq(){
			return GeomUtil.lengthSquared(this.line);
		}
		
		public double recipLength(){
			process();
			return this.rlen;
		}
		public double[] offset(){
			process();
			return this.offset;
		}
		public double[] centerVector(){
			process();
			return this.cvec;
		}
		public LineWrapper(Line2D l){
			this.line=l;
		}
		public static LineWrapper of(Line2D l){
			return new LineWrapper(l);
		}

		public Point2D projectPointOntoLine(Point2D p) {
			//return GeomUtil.projectPointOntoLine(this.getLine(), p);
			double[] vec=this.vector();
			double[] off=this.offset();
			
			double[] pvec = new double[]{p.getX()+off[0],p.getY()+off[1]};
			double dot=dot(vec, pvec);
			double proj=dot*this.recipLength()*this.recipLength();
			double[] projVec= new double[]{proj*vec[0],proj*vec[1]};
			return new Point2D.Double(projVec[0]-off[0], projVec[1]-off[1]);
		}

		@Override
		public int compareTo(LineWrapper arg0) {
			return Double.compare(length(),arg0.length());
		}

		public boolean intersectsLine(LineWrapper line2) {
			return this.getLine().intersectsLine(line2.getLine());
		}
	}
	
	
	public static class BoundingBox{
		private Rectangle2D rect;
		private List<Tuple<Line2D,Line2D>> lines;

		private CachedSupplier<Double> totLen = CachedSupplier.of(()->{
			return lines.stream().mapToDouble(l->length(l.k())).sum();
		});

		private CachedSupplier<Shape> convexHull = CachedSupplier.of(()->{
			return lines.stream().map(t->t.k()).flatMap(l->Stream.of(l.getP1(),l.getP2())).collect(GeomUtil.convexHull());
		});
		
		
		public double getLineDensity(){
			return Math.pow(totLen.get(),2)/(rect.getWidth()*2+rect.getHeight()*2);
		}
		
		
		public Rectangle2D getRect(){ //son!
			return rect;
		}
		
		public Rectangle2D getCenterOfMassRect(){
			Point2D center=getCenterOfMass();
			return new Rectangle2D.Double(center.getX()-rect.getWidth()/2, center.getY()-rect.getHeight()/2, rect.getWidth(), rect.getHeight());
		}
		public Point2D getCenterOfMass(){
			return centerOfMass(lines.stream().map(l->l.k()).collect(Collectors.toList()));
		}
		
		public double getTotalLineLength(){
			return totLen.get();
		}
		public static BoundingBox of(Rectangle2D rect, List<Line2D> lines){
			return ofCombo(rect, lines.stream().map(l->Tuple.of(l,l)).collect(Collectors.toList()));
		}
		public static BoundingBox ofCombo(Rectangle2D rect, List<Tuple<Line2D,Line2D>> lines){
			BoundingBox bb= new BoundingBox();
			bb.rect=rect;
			bb.lines=lines;
			return bb;
		}
		
		public BoundingBox resize(double width, double height){
			double cx= rect.getCenterX();
			double cy= rect.getCenterY();
			double x1= cx-width/2;
			double y1= cy-height/2;
			
			return BoundingBox.ofCombo(new Rectangle2D.Double(x1,y1,width,height), lines);
		}
		
		
		public int numberOfSplitLines(){
			
			return (int)lines.stream()
							 .flatMap(l->Stream.of(l.v().getP1(), l.v().getP2()))
							 .filter(p->!rect.contains(p))
							 .count();
		}
		
		public double getPercentageOfSplitLines(){
			return numberOfSplitLines()/((double)lines.size());
		}
		
		public Shape getConvuxHull(){
			return convexHull.get();
		}
		
		public Rectangle2D getCenteredConvexRect(){
			Shape s=getConvuxHull();
			
			Rectangle2D bb2=s.getBounds2D();
			return new Rectangle2D.Double(bb2.getCenterX()-rect.getWidth()/2,bb2.getCenterY()-rect.getHeight()/2, rect.getWidth(),rect.getHeight());
		}
		
		public Shape getConvuxHullOnlyInside(){
			return lines.stream()
						.filter(l->rect.contains(l.v().getP1()) && rect.contains(l.v().getP2()))
					    .map(t->t.k())
					    .flatMap(l->Stream.of(l.getP1(),l.getP2()))
					    .collect(GeomUtil.convexHull());
		}
	
		
	}
	
	public static List<BoundingBox> getBoundingBoxesContaining(double width, double height, List<Line2D> lines, double minDistance){
		double dx=1;
		double dy=1;

		List<Point2D> pts = lines.stream().flatMap(l->Stream.of(l.getP1(),l.getP2())).collect(Collectors.toList());
		double[] xs = pts.stream()
				           .flatMap(p->Stream.of(p.getX()-dx, p.getX()-width+dx))
				           .mapToDouble(d->d)
				           .sorted()
				           .distinct()
				           .toArray();
		double[] ys = pts.stream()
				           .flatMap(p->Stream.of(p.getY()-dy, p.getY()-height+dy))
				           .mapToDouble(d->d)
				           .sorted()
				           .distinct()
				           .toArray();
		
		
		
		//for n lines, that's 2n points. That's 4n x's and 4n y's. That's 16n^2 combinations. 
		//For n=100, that's 160,000 shapes (probably too many), still, keeping it naive for now
		
		
		//there are 3 criteria we care about for triage:
		//1. Should have lots of points in it
		//2. Should have long length of line inside
		//3. Internal bounding box area of lines should be pretty close to the bounding box area of the actual box.

		//We're going to focus on 1 and 2 for now.
		//For 1, we're going to make sure it has at least 4 points inside
		
		
		List<Point2D> alreadyGot = new ArrayList<Point2D>();
		
		
		
		return Arrays.stream(xs)
		      .mapToObj(x->x)
		      .flatMap(x->Arrays.stream(ys).mapToObj(y->new double[]{x,y}))
		      .map(xy->new Rectangle2D.Double(xy[0], xy[1], width, height))
		      .map(r->Tuple.of(r,pts.stream().filter(p->r.contains(p)).collect(Collectors.toList())))
		      .filter(t->t.v().size()>=4)
		      .map(Tuple.vmap(pl->Tuple.of(pl,-pl.size()).withVComparator()))
		      .map(t->t.withVComparator())
		      .sorted()
		      .map(Tuple.vmap(v->v.k()))
		      .map(t->Tuple.of(t.k(), lines.stream().map(l->Tuple.of(l,getLineInside(l,t.k()))).filter(o->o.v().isPresent()).map(Tuple.vmap(o->o.get())).map(Tuple::swap).collect(Collectors.toList())))
		      .map(t->BoundingBox.ofCombo(t.k(), t.v()))
		      .limit(1000) //shouldn't do this, but trying to stop an explosion
		      .filter(b->area(b.getConvuxHull()) > area(b.getConvuxHull().getBounds2D())*0.6)
		      .filter(b->area(b.getConvuxHull()) > area(b.getRect())*0.5)
		      .filter(b->area(b.getConvuxHullOnlyInside()) > area(b.getConvuxHull())*0.4)
		      .map(b->Tuple.of(b,b.numberOfSplitLines()))
		      .map(t->t.swap())
		      .collect(Tuple.toGroupedMap())
		      .entrySet()
		      .stream()
		      .map(Tuple::of)
		      .map(t->t.withKComparator())
		      .sorted()
		      .flatMap(t->t.v().stream().map(v->Tuple.of(v,-v.getTotalLineLength()).withVComparator()).sorted().map(t1->t1.k()))
		      .filter(r->{
		    	 Point2D c=r.getCenterOfMass();
		    	 if(alreadyGot.stream().filter(p->p.distance(c)<minDistance).findFirst().isPresent()){
		    		 return false;
		    	 }
		    	 alreadyGot.add(c);
		    	 return true;
		      })
		      .map(b->Tuple.of(b,-b.getLineDensity()).withVComparator())
		      .sorted()
		      .map(t->t.k())
		      .limit(5)
		      .collect(Collectors.toList());
		
	}
	
	public static List<List<Line2D>> groupMultipleBondsPreGrouped(List<List<LineWrapper>> lineLists, double maxDeltaTheta, double maxDeltaOffset, double minProjectionRatio, double minLargerProjectionRatio){
		return lineLists.stream()
					    .flatMap(ll->groupMultipleBonds(ll, maxDeltaTheta,  maxDeltaOffset,  minProjectionRatio,  minLargerProjectionRatio).stream())
					    .collect(Collectors.toList());
		
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
		
		public static List<List<Line2D>> groupMultipleBonds(List<LineWrapper> lines, double maxDeltaTheta, double maxDeltaOffset, double minProjectionRatio, double minLargerProjectionRatio){
			

			double cosThetaCutoff=Math.cos(maxDeltaTheta);
			
			List<LineWrapper> lws=lines;
			
			
			
			return GeomUtil.groupThings(lws, (lt)->{
				LineWrapper line1=lt.k();
				LineWrapper line2=lt.v();
								
				
				if(Math.abs(line1.cosTheta(line2))>cosThetaCutoff){
					double rej = Math.abs(line1.centerRejection(line2));
					if(rej<maxDeltaOffset){
						if(minProjectionRatio<=0 && minLargerProjectionRatio<=0)return true;
						double proj1=line1.projectionOffset(line2.getLine().getP1());
						double proj2=line1.projectionOffset(line2.getLine().getP2());
						
						if(proj1<0)proj1=0;
						if(proj1>line1.length()){
							proj1=line1.length();
						}
						if(proj2<0)proj2=0;
						if(proj2>line1.length()){
							proj2=line1.length();
						}
						double intersectLength=Math.abs(proj1-proj2);
						double cutoffRat1=Math.max(line1.recipLength(), line2.recipLength());
						double cutoffRat2=Math.min(line1.recipLength(), line2.recipLength());
						if(        intersectLength*cutoffRat1>minProjectionRatio 
								&& intersectLength*cutoffRat2>minLargerProjectionRatio
								){
							return true;
						}
					}
				}
				return false;
			})
			.stream()
			.map(l->l.stream().map(lw->lw.getLine()).collect(Collectors.toList()))
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
	public static List<Tuple<Line2D, Integer>> reduceMultiBonds(List<List<LineWrapper>> lines, double maxDeltaTheta, double maxDeltaOffset, double minProjectionRatio, double minLargerProjectionRatio, double maxStitchLineDistanceDelta, Consumer<Line2D> rejected){
		
		return groupMultipleBondsPreGrouped(lines,maxDeltaTheta,maxDeltaOffset,minProjectionRatio,minLargerProjectionRatio)
				.stream()
				//.map(l->GeomUtil.stitchSufficientlyStitchableLines(l, maxStitchLineDistanceDelta))
				.map(l->GeomUtil.getLineOffsetsToLongestLine(l))
				.map(l->{
					double doff=l.stream()
								 .filter(t->Math.abs(t.v())>GeomUtil.ZERO_DISTANCE_TOLERANCE)
								 .mapToDouble(t->t.v())
								 .min()
								 .orElse(1.0);
					List<Line2D> nlines= l.stream()
						 .map(Tuple.vmap(d->(int)Math.round(d/doff)))
						 
						 .map(t->t.swap())
						 .collect(Tuple.toGroupedMap())
						 .entrySet()
						 .stream()
						 .map(Tuple::of)
						 .map(Tuple.vmap(v1->GeomUtil.getPairOfFarthestPoints(vertices(v1))))
						 .map(Tuple.vmap(v1->new Line2D.Double(v1[0], v1[1])))
						 .map(t->t.v())
						 .collect(Collectors.toList());
					
					return nlines;
				})
				//.map(l->l.stream().map(t->t.k()).collect(Collectors.toList()))
				.map(l->Tuple.of(l,l.size()))
				.map(Tuple.kmap(l->{
					
					List<Line2D> sort= l.stream()
					   .map(s->Tuple.of(s,length(s)).withVComparator()) //always choose longer line
					   .sorted(Comparator.reverseOrder())
			           .map(l1->l1.k())
			           .collect(Collectors.toList());
					
					Line2D longest=sort.remove(0);
					sort.forEach(l1->rejected.accept(l1));
					return longest;
				}))
				.collect(Collectors.toList());
	}
	
	
	
	public static List<Line2D> stitchEdgesInMultiBonds(List<LineWrapper> lines, double maxDeltaTheta, double maxDeltaOffset, double minProjectionRatio, double minLargerProjectionRatio, double maxStitchLineDistanceDelta){
		
		return groupMultipleBonds(lines,maxDeltaTheta,maxDeltaOffset,minProjectionRatio,minLargerProjectionRatio)
				.stream()
				.flatMap(l->GeomUtil.stitchSufficientlyStitchableLines(l, maxStitchLineDistanceDelta).stream())
				.collect(Collectors.toList());
	}

	public static double[] asVector(Line2D line){
		double dx=line.getX2()-line.getX1();
		double dy=line.getY2()-line.getY1();
		return new double[]{dx,dy};
		
	}

	public static Shape getClosestShapeTo(Collection<Shape> shapes, Point2D pnt){
		return shapes.stream()
				      .map(s->Tuple.of(s,distanceTo(s, pnt)).withVComparator())
				      .min(CompareUtil.naturalOrder())
				      .map(t->t.k())
				      .orElse(null);
	}
	
	public static ShapeWrapper getClosestShapeToSW(Collection<ShapeWrapper> shapes, Point2D pnt){
		return shapes.stream()
				      .map(s->Tuple.of(s,s.distanceTo(pnt)).withVComparator())
				      .min(CompareUtil.naturalOrder())
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
			return new Line2D.Double(pnt,line.getP2());
		}
	}

	public static ConnectionTable getConnectionTable(List<Tuple<Line2D, Integer>> linest, List<ShapeWrapper> likelyNodes,
			double maxDistanceRatioNonLikely, 
			double maxDistanceRatioLikely, 
			double maxDistanceRatioPerLine,
			double minPerLineDistanceRatioForIntersection,
			double maxCandidateRatioForIntersectionWithNeighbor,
			Predicate<Line2D> acceptNewLine) {
			
			return ConnectionTable.fromLinesAndOrders(linest)
			                      .getNewConnectionTable(likelyNodes, 
			                    		  maxDistanceRatioNonLikely, 
			                    		  maxDistanceRatioLikely, 
			                    		  maxDistanceRatioPerLine, 
			                    		  minPerLineDistanceRatioForIntersection,
			                    		  
			                    		  maxCandidateRatioForIntersectionWithNeighbor,
			                    		  acceptNewLine);
			
				
	}

	public static Point2D projectPointOntoLine(Line2D l, Point2D p){
		double[] vec=asVector(l);
		double[] pvec = new double[]{p.getX()-l.getX1(),p.getY()-l.getY1()};
		double dot=dot(vec, pvec);
		double proj=dot/Math.pow(l2Norm(vec),2);
		double[] projVec= new double[]{proj*vec[0],proj*vec[1]};
		return new Point2D.Double(projVec[0]+l.getX1(), projVec[1]+l.getY1());
	}
	public static Line2D resizeLine(Line2D l1, double size){
		double[] v=normalize(asVector(l1));
		v[0]=v[0]*size;
		v[1]=v[1]*size;
		double[] off=new double[]{l1.getX1(),l1.getY1()};
		
		return new Line2D.Double(l1.getP1(), new Point2D.Double(v[0]+off[0], v[1]+off[1]));
	}
	
	public static Tuple<Point2D, Double> projectPointOntoLineWithRejection(Line2D l, Point2D p){
		double[] vec=asVector(l);
		double[] pvec = new double[]{p.getX()-l.getX1(),p.getY()-l.getY1()};
		double dot=dot(vec, pvec);
		double proj=dot/Math.pow(l2Norm(vec),2);
		double[] projVec= new double[]{proj*vec[0],proj*vec[1]};
		Point2D pr=new Point2D.Double(projVec[0]+l.getX1(), projVec[1]+l.getY1());;
		double r=GeomUtil.rejection(vec, pvec);
		
		return Tuple.of(pr,r);
	}
	
	public static Point2D[] getPairOfFarthestPoints(Shape s){
		return getPairOfFarthestPoints(vertices(s));
	}
	public static Point2D[] getPairOfFarthestPoints(Point2D[] pts){
		return getPairOfFarthestPoints(Arrays.asList(pts));
	}
	public static Point2D[] getPairOfFarthestPoints(List<Point2D> pts){
		//not sure if this is really what should be done
		if(pts.size()==0)return new Point2D[]{new Point2D.Double(0,0),new Point2D.Double(0,0)};
		return eachCombination(pts)
		         .map(t->Tuple.of(t,t.k().distanceSq(t.v())).withVComparator())
		         .max(CompareUtil.naturalOrder())
		         .map(t->t.k())
		         .map(t->new Point2D[]{t.k(),t.v()})
		         .orElse(new Point2D[]{pts.get(0),pts.get(0)});
	}
	
	public static List<Shape> mergeOverlappingShapes(List<Shape> shapes, double minOverlapRatio){
		return GeomUtil.groupShapes(shapes, t->{
							if(intersects(t.k(),t.v())){
								double areaCombined=area(add(t.k(),t.v()));
								double areaIntersect=getIntersectionShape(t.k(), t.v()).map(s->area(s)).orElse(0.0);
								if(areaIntersect/areaCombined>minOverlapRatio){
									return true;
								}
							}
							return false;
						})
					.stream()
					.map(ll->ll.stream().reduce((s1,s2)->add(s1,s2)).orElse(null))
					.filter(s->s!=null)
				    .collect(Collectors.toList());
	}
	
	public static List<ShapeWrapper> mergeOverlappingShapesSW(List<ShapeWrapper> shapes, double minOverlapRatio){
		return GeomUtil.groupThings(shapes, t->{
							double areaIntersect=getIntersectionShape(t.k(), t.v()).map(s->s.getArea()).orElse(0.0);
							if(areaIntersect>ZERO_DISTANCE_TOLERANCE){
								double areaCombined=add(t.k(),t.v()).getArea();
								
								if(areaIntersect/areaCombined>minOverlapRatio){
									return true;
								}
							}
							return false;
						})
					.stream()
					.map(ll->ll.stream().reduce((s1,s2)->add(s1,s2)).orElse(null))
					.filter(s->s!=null)
				    .collect(Collectors.toList());
	}
	
	public static <T> Stream<Tuple<T,T>> eachCombination(List<T> list){
		return IntStream.range(0, list.size())
		         .mapToObj(i->IntStream.range(i+1, list.size())
		        		               .mapToObj(j->Tuple.of(i,j))
		        		  )
		         .flatMap(s->s)
		         .map(Tuple.vmap(i->list.get(i)))
		         .map(Tuple.kmap(i->list.get(i)));
		         
		         
	}
	
	public static Shape makeShapeAround(Point2D p, double rad){
		return new Rectangle2D.Double(p.getX()-rad, p.getY()-rad, rad*2, rad*2);
	}
	
	/**
	 * This method scales the shape based on its bounding box. The value given will scale the shape
	 * such that the new bounding box width is 2*dr larger than it was before. Note that this
	 * may scale a shape quite differently depending on whether they are "tall" or "wide"
	 * @param s
	 * @param dr
	 * @return
	 */
	public static Shape growShape(Shape s, double dr){
		return ShapeWrapper.of(s).growShapeBounds(dr).getShape();
	}
	public static List<Line2D> splitLineIn2(Line2D l){
		double cx=l.getBounds2D().getCenterX();
		double cy=l.getBounds2D().getCenterY();
		Point2D cp = new Point2D.Double(cx,cy);
		Line2D l1 = new Line2D.Double(l.getP1(),cp);
		Line2D l2 = new Line2D.Double(l.getP2(),cp);
		return Stream.of(l1,l2).collect(Collectors.toList());
		
	}
	
	public static Line2D stitchLines(Line2D l1,Line2D l2){
		LineDistanceCalculator ldc=LineDistanceCalculator.from(l1,l2);
		return ldc.getLineFromFarthestPoints();
	}
	
	public static Line2D longestLineFrom(List<Point2D> points){
		Point2D[] pts=getPairOfFarthestPoints(points);
		return new Line2D.Double(pts[0], pts[1]);
	}
	
	public static List<Line2D> stitchSufficientlyStitchableLines(List<Line2D> lines, double maxDistance){
		
		//Look through each line,
		List<Tuple<Line2D,Double>> lineLengthSet = lines.stream()
				 .map(l->Tuple.of(l,length(l)).withVComparator())
				 .sorted()
				 .collect(Collectors.toList());
		
		Map<Line2D,Integer> group = new HashMap<>();
		
		for(int i=0;i<lines.size();i++){
			group.put(lines.get(i),i);
		}
		
		eachCombination(lineLengthSet).forEach(t->{
			double dl1=t.k().v();
			double dl2=t.v().v();
			double sum=dl1+dl2;
			Line2D l1=t.k().k();
			Line2D l2=t.v().k();
			Line2D nl=stitchLines(l1,l2);
			double dln=length(nl);
			if(Math.abs(sum-dln)<maxDistance){
				int g1=group.get(l1);
				int g2=group.get(l2);
				int ng=Math.min(g1, g2);
				group.put(l1, ng);
				group.put(l2, ng);
			}
		});
		
		return group.entrySet()
		     .stream()
		     .map(Tuple::of)
		     .map(t->t.swap())
		     .collect(Tuple.toGroupedMap())
		     .values()
		     .stream()
		     .map(ll->vertices(ll))
		     .map(vt->longestLineFrom(Arrays.asList(vt)))
		     .collect(Collectors.toList());
		
				 
		
	}
	
	public static double[] asVector(Point2D p){
		return new double[]{p.getX(),p.getY()};
	}
	public static double[] addVectors(double[] a, double[] b){
		return IntStream.range(0, a.length)
			         .mapToDouble(i->a[i]+b[i])
			         .toArray();
	}
	public static double[] negate(double[] a){
		return IntStream.range(0, a.length)
		         .mapToDouble(i->-a[i])
		         .toArray();
	}
	
	
	public static List<Tuple<Line2D,Double>> getLineOffsetsToLongestLine(List<Line2D> lines){
		Line2D longest= lines.stream()
				             .map(l->Tuple.of(l,length(l)).withVComparator())
				             .max(CompareUtil.naturalOrder())
				             .map(t->t.k())
				             .orElse(null);
		if(longest==null)return null;
		double[] lvec= GeomUtil.asVector(longest);
		double[] offset = negate(new double[]{longest.getX1(),longest.getY1()});
		
		
		return lines.stream()
			         .map(l->Tuple.of(l,findCenterOfShape(l)))
			         .map(Tuple.vmap(p->asVector(p)))
			         .map(Tuple.vmap(b->addVectors(b,offset)))
			         .map(Tuple.vmap(b->rejection(lvec,b)))
			         .collect(Collectors.toList());
		
	}
	
	public static Optional<Line2D> getLineInside(Line2D l, Shape s){

		boolean p1Inside = s.contains(l.getP1());
		boolean p2Inside = s.contains(l.getP2());
		
		if(p1Inside && p2Inside)return Optional.of(l);
		
		
		List<Point2D> ips=GeomUtil.getAllIntersections(s, l);
		
		if(ips.isEmpty())return Optional.empty();
		
		if(ips.size()>2){
//			System.out.println("Line:" + l);
//			System.out.println("Shape:" + Arrays.toString(vertices(s)));
//			System.out.println("Shape:" + Arrays.toString(vertices(convexHull2(vertices(s)))));
			throw new IllegalStateException("Line should not intersect with convex hull more than twice, maybe the shape isn't a convux hull?");
		}

		Point2D ip1 = ips.get(0);
		Point2D ip2 = (ips.size()==2)?ips.get(1):null;
		
		if(p1Inside && ip2==null){			
			Line2D ln = new Line2D.Double(l.getP1(),ip1);
			return Optional.of(ln);
		}else if(p2Inside && ip2==null){
			Line2D ln = new Line2D.Double(ip1,l.getP2());
			return Optional.of(ln);
		}else{
			if(ip2==null){
				//Could be tangent to one point actually, in this case just return nothing
				return Optional.empty();
			}else{
				return Optional.of(new Line2D.Double(ip1,ip2));
			}
		}
		
	}
	
	public static List<Line2D> getLinesNotInside(Line2D l, Shape s){
		return ShapeWrapper.of(s).getLinesNotInside(l);
	}
	
	public static List<Line2D> getLinesNotInside(Line2D l, Collection<Shape> shapes){
		List<Line2D> start = Arrays.asList(l);
		
		for(Shape s: shapes){
			start=start.stream()
					   .flatMap(ls->getLinesNotInside(ls,s).stream())
					   .collect(Collectors.toList());
		}
		return start;
	}
	
	public static List<Line2D> getLinesNotInsideSW(Line2D l, Collection<ShapeWrapper> shapes){
		List<Line2D> start = Arrays.asList(l);
		
		for(ShapeWrapper s: shapes){
			start=start.stream()
					   .flatMap(ls->s.getLinesNotInside(ls).stream())
					   .collect(Collectors.toList());
			if(start.isEmpty())break;
		}
		return start;
	}
	
	public static Optional<Line2D> getLongestLineNotInside(Line2D l, Collection<Shape> shapes){
		//if(true)return Optional.of(l);
		
		return getLinesNotInside(l,shapes)
		       .stream()
		       .map(l1->Tuple.of(l1,length(l1)).withVComparator())
		       .max(Comparator.naturalOrder())
		       .map(t->t.k());
		       
	}
	public static Optional<Shape> getIntersectionShape(Shape s1, Shape s2){
		return getIntersectionShape(ShapeWrapper.of(s1), ShapeWrapper.of(s2))
				.map(s->s.getShape());
	}
	
	public static Optional<ShapeWrapper> getIntersectionShape(ShapeWrapper s1, ShapeWrapper s2){
		if(!s1.getBounds().intersects(s2.getBounds()))return Optional.empty();
		   
		List<Point2D> pointsInside1 = Arrays.stream(s1.getVerts())
				                            .filter(p->s2.contains(p))
				                            .collect(Collectors.toList());
		List<Point2D> pointsInside2 = Arrays.stream(s2.getVerts())
                .filter(p->s1.contains(p))
                .collect(Collectors.toList());
		
		List<Point2D> pp=intersectingPoints(s1,s2);
		
		List<Point2D> allPoints = new ArrayList<Point2D>();
		allPoints.addAll(pointsInside1);
		allPoints.addAll(pointsInside2);
		allPoints.addAll(pp);
		
		if(allPoints.isEmpty() || allPoints.size()<3)return Optional.empty();
		return Optional.ofNullable(convexHull2(allPoints.toArray(new Point2D[0])))
				.map(s->ShapeWrapper.of(s));
	}
	
	/**
	 * Returns the signed area of a triangle, assuming that a CW orientation is positive,
	 * and a CCW orientation is negative. This is accomplished by centering the triangle
	 * so that the first point (P1) is on the origin, and then taking the signed "rejection" of the the
	 * vector from P1->P3 onto the vector P1->P2. That value, multiplied by the length of P1->P2 gives
	 * the area of the parallelogram. Divide the answer by 2, and it's the area of the triangle.
	 * @param verts
	 * @return
	 */
	public static double areaTriangle(Point2D p1, Point2D p2, Point2D p3){
		//base x height
		//base x rejection of non-base onto base
		double[] vp1=asVector(p1);
		double[] vp2=asVector(p2);
		double[] vp3=asVector(p3);
		vp2=addVectors(vp2,negate(vp1));
		vp3=addVectors(vp3,negate(vp1));
		
		//orthoDot is just the non-normalized scalar rejection of vector 2 onto vector 1
		//which is the same as the area of the parallelogram formed by V1 and V2 
		return orthoDot(vp2,vp3)/2;
		
	}
	/**
	 * Returns the signed area of a triangle, assuming that a CW orientation is positive,
	 * and a CCW orientation is negative. This just delegates to the {@link #areaTriangle(Point2D, Point2D, Point2D)}
	 * method.
	 * @param s
	 * @return
	 */
	public static double areaTriangle(Shape s){
		Point2D[] pts=vertices(s);
		if(pts.length!=3)throw new IllegalStateException("Triangles must have 3 points");
		return areaTriangle(pts[0],pts[1],pts[2]);
	}
	
	/**
	 * Returns the signed area of a set of vertices, assuming that a CW orientation is positive,
	 * and a CCW orientation is negative. This is accomplished by breaking the vertices into triangles
	 * by fetching every-other triangle formed by the vertices. Once all areas that way are calculated,
	 * it then recursively calculates the area for the vertices "left over" and adds that area. All of this
	 * eventually delegates to {@link #areaTriangle(Point2D, Point2D, Point2D)}. Note that this method
	 * works well for convex hulls and for shapes that do not have "crossing" edge. More specifically, the 
	 * results of this method will not be well-defined if the provided vertices form an edge that 'crosses"
	 * another edge, and therefore has an inconsistent definition of what it means to be "inside" of the shape. 
	 * @param verts
	 * @return
	 */
	public static double areaVerticesCW(Point2D[] verts){
		if(verts.length<3)return 0;
		if(verts.length==3)return areaTriangle(verts[0],verts[1],verts[2]);
		List<Point2D> oddPoints=new ArrayList<Point2D>();
		double areaSoFar=0;
		oddPoints.add(verts[0]);
		for(int i=1;i<verts.length;i+=2){
			Point2D p1=verts[i-1];
			
			Point2D p2=verts[i];
			
			Point2D p3=verts[(i+1)%verts.length];
			oddPoints.add(p3);
			areaSoFar+=areaTriangle(p1,p2,p3);
		}
		return areaSoFar + areaVerticesCW(oddPoints.toArray(new Point2D[]{}));
		
	}
	
	public static Shape shapeFromVertices(Point2D[] pnts){
		Path2D path = new Path2D.Double();
		
		for(int i=0;i<pnts.length;i++){
			Point2D p=pnts[i];
			if(i>0){
				path.lineTo(p.getX(), p.getY());
			}else{
				path.moveTo(p.getX(), p.getY());
			}
		}
		path.closePath();
		return path;
	}
	
	/**
	 * Returns the absolute (positive) area of a shape.
	 * @param s
	 * @return
	 */
	public static double area(Shape s){
		if(s instanceof Rectangle2D){
			return ((Rectangle2D)s).getWidth()*((Rectangle2D)s).getHeight();
		}
		return Math.abs(areaVerticesCW(vertices(s)));
	}
	
	public static double signedArea(Shape s){
		return areaVerticesCW(vertices(s));
	}
	
	
	@Deprecated
	public static double poorMansArea(Shape s){
		Rectangle2D r = s.getBounds2D();
		int hits=0;
		int tot=10000;
		double wid=r.getWidth();
		double hi=r.getHeight();
		double tx;
		double ty;
		for(int i=0;i<tot;i++){
			tx=r.getMinX()+ ThreadLocalRandom.current().nextDouble(0, r.getWidth());
			ty=r.getMinY()+ ThreadLocalRandom.current().nextDouble(0, r.getHeight());
			if(s.contains(tx, ty))hits++;
		}
		return wid*hi*(hits)/(double)tot;
	}
	public static AffineTransform getTransformFromLineToLine(Line2D source,Line2D destination, boolean flip){
		AffineTransform at = new AffineTransform();
		
		double[] v1 = asVector(source);
		double[] v2 = asVector(destination);
		double l1=l2Norm(v1);
		double l2=l2Norm(v2);
		
		double scale  =l2/l1;
		
		double thetap1 = angleFromVec1ToVec2(v1,new double[]{1,0});
		double thetap2 = angleFromVec1ToVec2(new double[]{1,0},v2);
		
		double scaleY = (flip)?-scale:scale;
		
		double sy=-source.getP1().getY();
		double dy=destination.getP1().getY();
		
		at.translate(destination.getP1().getX(), dy);

		at.rotate(thetap2);
		at.scale(scale, scaleY);
		at.rotate(thetap1);
		
		at.translate(-source.getP1().getX(), sy);
		
		//at.translate(-destination.getP1().getX(), -destination.getP1().getY());
		
		return at;
	}
	
	public static double angleFromVec1ToVec2(double[] v1, double[] v2){
		double l1=l2Norm(v1);
		double l2=l2Norm(v2);
		
		double dot=dot(v1,v2);
		double cosTheta=dot/(l1*l2);
		double thetap = Math.acos(cosTheta);
		double rej = rejection(v1,v2);
		if(rej<0){
			thetap=thetap*-1;
		}
		return thetap;
	}
	
	
	public static List<Point2D> splitIntoNPieces(Line2D l, int n){
		
		double[] v=asVector(l);
		double[] off=new double[]{l.getP1().getX(),l.getP1().getY()};
		double rlen=1.0/(n);
		
		return IntStream.range(0,n+1)
				 .mapToObj(i->new double[]{i*v[0]*rlen+off[0], i*v[1]*rlen+off[1]})
				 
				 .map(d->new Point2D.Double(d[0],d[1]))
				 .peek(d->System.out.println(d))
				 .collect(Collectors.toList());
	}
	public static Collector<Point2D,List<Point2D>,Point2D> averagePoint(){
		
		return new Collector<Point2D,List<Point2D>,Point2D>(){

			@Override
			public BiConsumer<List<Point2D>, Point2D> accumulator() {

				return (l,p)->{
					l.add(p);
				};
			}

			@Override
			public Set<java.util.stream.Collector.Characteristics> characteristics() {
				return new HashSet<java.util.stream.Collector.Characteristics>();
			}

			@Override
			public BinaryOperator<List<Point2D>> combiner() {
				return (l1,l2)->{
					l1.addAll(l2);
					return l1;
				};
			}

			@Override
			public Function<List<Point2D>, Point2D> finisher() {
				return (pts)->GeomUtil.findCenterOfVertices(pts);
			}

			@Override
			public Supplier<List<Point2D>> supplier() {
				return ()->new ArrayList<>();
			}
			
		};
	}
	
	public static Collector<Point2D,List<Point2D>,Shape> convexHull(){
		
		return new Collector<Point2D,List<Point2D>,Shape>(){

			@Override
			public BiConsumer<List<Point2D>, Point2D> accumulator() {

				return (l,p)->{
					l.add(p);
				};
			}

			@Override
			public Set<java.util.stream.Collector.Characteristics> characteristics() {
				return new HashSet<java.util.stream.Collector.Characteristics>();
			}

			@Override
			public BinaryOperator<List<Point2D>> combiner() {
				return (l1,l2)->{
					l1.addAll(l2);
					return l1;
				};
			}

			@Override
			public Function<List<Point2D>, Shape> finisher() {
				return (pts)->GeomUtil.convexHull2(pts.stream().toArray(i->new Point2D[i]));
			}

			@Override
			public Supplier<List<Point2D>> supplier() {
				return ()->new ArrayList<>();
			}
			
		};
	}
	

	public static Collector<Shape,List<Shape>,Shape> joined(){
		
		return new Collector<Shape,List<Shape>,Shape>(){

			@Override
			public BiConsumer<List<Shape>, Shape> accumulator() {

				return (l,p)->{
					l.add(p);
				};
			}

			@Override
			public Set<java.util.stream.Collector.Characteristics> characteristics() {
				return new HashSet<java.util.stream.Collector.Characteristics>();
			}

			@Override
			public BinaryOperator<List<Shape>> combiner() {
				return (l1,l2)->{
					l1.addAll(l2);
					return l1;
				};
			}

			@Override
			public Function<List<Shape>, Shape> finisher() {
				return (pts)->pts.stream().reduce((s1,s2)->add(s1,s2)).orElse(null);
			}

			@Override
			public Supplier<List<Shape>> supplier() {
				return ()->new ArrayList<>();
			}
			
		};
	}
	
public static Collector<ShapeWrapper,List<ShapeWrapper>,ShapeWrapper> joinedSW(){
		
		return new Collector<ShapeWrapper,List<ShapeWrapper>,ShapeWrapper>(){

			@Override
			public BiConsumer<List<ShapeWrapper>, ShapeWrapper> accumulator() {

				return (l,p)->{
					l.add(p);
				};
			}

			@Override
			public Set<java.util.stream.Collector.Characteristics> characteristics() {
				return new HashSet<java.util.stream.Collector.Characteristics>();
			}

			@Override
			public BinaryOperator<List<ShapeWrapper>> combiner() {
				return (l1,l2)->{
					l1.addAll(l2);
					return l1;
				};
			}

			@Override
			public Function<List<ShapeWrapper>, ShapeWrapper> finisher() {
				return (pts)->pts.stream().reduce((s1,s2)->add(s1,s2)).orElse(null);
			}

			@Override
			public Supplier<List<ShapeWrapper>> supplier() {
				return ()->new ArrayList<>();
			}
			
		};
	}
	
	public static Collector<Shape,List<Shape>,Optional<Shape>> union(){
		
		return new Collector<Shape,List<Shape>,Optional<Shape>>(){

			@Override
			public BiConsumer<List<Shape>, Shape> accumulator() {

				return (l,p)->{
					l.add(p);
				};
			}

			@Override
			public Set<java.util.stream.Collector.Characteristics> characteristics() {
				return new HashSet<java.util.stream.Collector.Characteristics>();
			}

			@Override
			public BinaryOperator<List<Shape>> combiner() {
				return (l1,l2)->{
					l1.addAll(l2);
					return l1;
				};
			}

			@Override
			public Function<List<Shape>, Optional<Shape>> finisher() {
				return (pts)->pts.stream().map(s->new Area(s)).reduce((s1,s2)->{s1.add(s2);return s1;}).map(s->(Shape)s);
			}

			@Override
			public Supplier<List<Shape>> supplier() {
				return ()->new ArrayList<>();
			}
			
		};
	}
	
	public static Collector<Double,List<Double>,Double> median(){
		
		return new Collector<Double,List<Double>,Double>(){

			@Override
			public BiConsumer<List<Double>, Double> accumulator() {

				return (l,p)->{
					l.add(p);
				};
			}

			@Override
			public Set<java.util.stream.Collector.Characteristics> characteristics() {
				return new HashSet<java.util.stream.Collector.Characteristics>();
			}

			@Override
			public BinaryOperator<List<Double>> combiner() {
				return (l1,l2)->{
					l1.addAll(l2);
					return l1;
				};
			}

			@Override
			public Function<List<Double>, Double> finisher() {
				return (pts)->{
					double[] arr=pts.stream().mapToDouble(d->d).sorted().toArray();
					if(arr.length==1)return arr[0];
					if(arr.length==0)return 0.0;
					if(arr.length%2==1){
						return arr[arr.length/2];
					}else{
						return (arr[arr.length/2-1] +arr[arr.length/2])/2;
					}
				};
			}

			@Override
			public Supplier<List<Double>> supplier() {
				return ()->new ArrayList<>();
			}
			
		};
	}
	

	public static <T> Collector<T,List<T>,List<List<T>>> groupThings(Predicate<Tuple<T,T>> accept){
		
		
		
		return new Collector<T,List<T>,List<List<T>>>(){

			@Override
			public BiConsumer<List<T>, T> accumulator() {
				return (l,t)->l.add(t);
			}

			@Override
			public Set<java.util.stream.Collector.Characteristics> characteristics() {
				return new HashSet<java.util.stream.Collector.Characteristics>();
			}

			@Override
			public BinaryOperator<List<T>> combiner() {

				return (l1,l2)->{
					l1.addAll(l2);
					return l1;
				};
			}

			@Override
			public Function<List<T>, List<List<T>>> finisher() {
				return l->groupThings(l,accept);
			}

			@Override
			public Supplier<List<T>> supplier() {
				return ()->new ArrayList<T>();
			}

		};
	}



    
}
