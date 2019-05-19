package gov.nih.ncats.molvec;

import java.awt.Point;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.DataBuffer;
import java.awt.image.MultiPixelPackedSampleModel;
import java.awt.image.Raster;
import java.awt.image.RenderedImage;
import java.awt.image.SampleModel;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.BiFunction;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import javax.imageio.ImageIO;



import gov.nih.ncats.molvec.algo.Tuple;
import gov.nih.ncats.molvec.image.ImageUtil;
import gov.nih.ncats.molvec.image.TiffTags;
import gov.nih.ncats.molvec.image.binarization.AdaptiveThreshold;
import gov.nih.ncats.molvec.image.binarization.SigmaThreshold;
import gov.nih.ncats.molvec.image.binarization.SigmaThreshold.ImageStats;
import gov.nih.ncats.molvec.util.GeomUtil;
import gov.nih.ncats.molvec.util.GeomUtil.LineDistanceCalculator;
import gov.nih.ncats.molvec.util.GeomUtil.LineWrapper;

/**
 * A bitmap image
 */
public class Bitmap implements Serializable, TiffTags {
    private static final long serialVersionUID = 0x5f1f54d8fed49ab3l;
    private static final Logger logger =
        Logger.getLogger (Bitmap.class.getName ());

    private static final double EPS = 0.000001;
    public static final double DEFAULT_AEV_THRESHOLD = 1.5;

    private static final boolean DEBUG;
    static {
        boolean debug = false;
        try {
            debug = Boolean.getBoolean ("bitmap.debug");
        } catch (Exception ex) { }
        DEBUG = debug;
    }
    
    public static class BitmapScaled{
    	public int twidth;
    	public int theight;
    	public int[][] ccount;
    	public List<int[]> xys;
    	public int tcount;
    	
    	
    	
    	
    	public static BitmapScaled of(Bitmap r, int DEF_WIDTH, int DEF_HEIGHT){
	    	int twidth = r.width();
			int theight = r.height();
			
			int[][] ccount = new int[DEF_WIDTH][DEF_HEIGHT];
			
			List<int[]> xys=r.getXYOnPoints()
							    .map(xy->{
							    	int cx = (xy[0] * DEF_WIDTH) / twidth;
							    	int cy = (xy[1] * DEF_HEIGHT) / theight;
							    	return new int[]{cx,cy,1};
							    })
							    .collect(Collectors.groupingBy(i->i[0]+","+i[1]))
							    .values()
							    .stream()
							    .map(il->{
							    	int[] r1=il.get(0);
							    	r1[2]=il.size();
							    	return r1;
							    })
							    .collect(Collectors.toList());
			
			int tcount= xys.stream().mapToInt(r1->r1[2]).sum();
			
			for(int i=0;i<twidth;i++){
				for(int j=0;j<theight;j++){
				   	int cx = (i * DEF_WIDTH) / twidth;
			    	int cy = (j * DEF_HEIGHT) / theight;
			    	ccount[cx][cy]++;
				}
			}
			BitmapScaled bms=new BitmapScaled();
			bms.ccount=ccount;
			bms.twidth=twidth;
			bms.theight=theight;
			bms.xys=xys;
			bms.tcount=tcount;
			return bms;
    	}
    }
    
    private Map<String,BitmapScaled> _scaleMap = new HashMap<>();
    
    public BitmapScaled getScaled(int nwid, int nhit){
    	return _scaleMap.computeIfAbsent(nwid+"_"+nhit, k->{
    		return BitmapScaled.of(this, nwid, nhit);
    	});
    }
    
    
    public static class CropBackedBitmap extends Bitmap{

    	private Rectangle r = null;
    	private Shape cropShape = null;
    	private Bitmap real = null;
    	private int x0;
    	private int y0;
    	private int x1;
    	private int y1;
    	private int w;
    	private int h;

    	private CachedSupplier<List<int[]>> onXYs = CachedSupplier.of(()->{
    		   		
    		if(w*h<400){
    			return real.crop(cropShape).getXYOnPoints().collect(Collectors.toList());
    		}
    		
    		return real.getXYOnPoints()
    				   .filter(xy->(xy[0]>=x0 && xy[0]<=x1) && (xy[1]>=y0 && xy[1]<=y1))
    				   .filter(xy->xy[0]==x1||xy[1]==y1||cropShape.contains(xy[0], xy[1]))
    				   .map(xy->new int[]{(int)(xy[0]-x0),(int)(xy[1]-y0)})
    				   .collect(Collectors.toList());
    	});
     	
    	
    	
		public CropBackedBitmap(Bitmap copy, Shape cropShape) {
			real = copy;
			this.cropShape=cropShape;
			
			r = cropShape.getBounds ();
		     
			w=r.width+1;
			h=r.height+1;
		    x1 = Math.min (real.width, r.x + r.width);
		    y1 = Math.min (real.height, r.y + r.height);
		    x0 = r.x;
		    y0 = r.y;
		}

		public int height() {
			return h;
		}
		
		public int width() {
			return w;
		}
		
		public Stream<int[]> getXYOnPoints() {
			return onXYs.get().stream();
		}
    	public double fractionPixelsOn(){
    		return onXYs.get().size()/((double)(width()*height()));
    	}
    }
    
    
    public CropBackedBitmap getLazyCrop(Shape c){
    	
    	Rectangle r = c.getBounds();
    	if(r.width==0||r.height==0)return null;
    	return new CropBackedBitmap(this,c);
    }
    
    public static class BitmapBuilder{
    	private Bitmap source;
    	
    	private int vblurRad=0;
    	private int hblurRad=0;
    	private int thresh=1;
    	
    	
    	private List<Shape> toRemove=new ArrayList<>();
    	
    	private int nwidth;
    	private int nheight;
    	private int blurRep=1;
    	
    	private double scale=1;
    	
    	public BitmapBuilder(Bitmap bm){
    		this.source=bm;
    		this.nwidth=source.width;
    		this.nheight=source.height;
    	}
    	
    	public BitmapBuilder scale(double s){
    		this.scale=s;
    		this.nwidth=(int)(source.width*scale);
    		this.nheight=(int)(source.height*scale);
    		return this;
    	}
    	
    	
    	public BitmapBuilder boxBlur(int rad){
    		return this.vblur(rad).hblur(rad);
    	}
    	public BitmapBuilder gaussBlur(int rad){
    		return this.vblur(2).hblur(2).blurRepeats(rad);
    	}
    	
    	public BitmapBuilder blurRepeats(int n){
    		blurRep=n;
    		return this;
    	}
    	
    	public BitmapBuilder vblur(int rad){
    		vblurRad =rad;
    		
    		return this;
    	}
    	public BitmapBuilder hblur(int rad){
    		hblurRad =rad;
    		return this;
    	}
    	public BitmapBuilder threshold(int t){
    		this.thresh=t;
    		return this;
    	}
    	
    	public Bitmap build(){
    		int[][] raw = new int[nwidth][nheight];
    		double iscale=1/scale;
    		for(int i=0;i<nwidth;i++){
    			for(int j=0;j<nheight;j++){
    				int x=(int)Math.round(iscale*i);
    				int y=(int)Math.round(iscale*j);
    				boolean ignore=toRemove.stream()
    				        .filter(r->r.contains(x, y))
    				        .findAny()
    				        .isPresent();
    				if(!ignore){
    					raw[i][j]=source.getAsInt((int)Math.round(iscale*i), (int)Math.round(iscale*j));
    				}
    			}
    		}
    		for(int i=0;i<this.blurRep;i++){
    			if(vblurRad>0)vblur(raw,vblurRad);
    			if(hblurRad>0)hblur(raw,hblurRad);
    		}
    		
    		Bitmap bm2 = new Bitmap(nwidth,nheight);
    		
    		for(int i=0;i<nwidth;i++){
    			for(int j=0;j<nheight;j++){
    				if(raw[i][j]>=thresh)bm2.set(i, j, true);
    			}
    		}
    		return bm2;
    	}
    	
    	
    	/**
    	 * Vertical "motion" Blur: This is a one-dimensional blurring of an image,
    	 * without renormalizing.
    	 * 
    	 * @param bmap
    	 *            bitmap as 2d int array
    	 * @param rad
    	 *            radius for vertical blur
    	 */
    	protected static void vblur(int[][] bmap, int rad) {
    		int sofar = 0;
    		List<Integer> added = new ArrayList<Integer>();
    		for (int i = 0; i < bmap.length; i++) {
    			for (int j = 0; j < rad * 2 + bmap[0].length; j++) {
    				// System.out.println(sofar + ":" + added.size());
    				if (j < bmap[0].length) {
    					added.add(bmap[i][j]);
    					sofar += bmap[i][j];
    				} else {
    					if (j < bmap[0].length + rad) {
    						added.add(bmap[i][bmap[0].length - 1]);
    						sofar += bmap[i][bmap[0].length - 1];
    					}
    				}
    				if (j >= rad) {
    					if (j - rad < bmap[0].length) {
    						bmap[i][j - rad] = sofar;
    					}
    					sofar -= added.get(0);
    					added.remove(0);
    				}
    			}
    		}
    	}

    	/**
    	 * Horizontal "motion" Blur: This is a one-dimensional blurring of an image,
    	 * without renormalizing.
    	 * 
    	 * @param bmap
    	 *            bitmap as 2d int array
    	 * @param rad
    	 *            radius for horizontal blur
    	 */
    	protected static void hblur(int[][] bmap, int rad) {
    		int sofar = 0;
    		List<Integer> added = new ArrayList<Integer>();
    		for (int j = 0; j < bmap[0].length; j++) {
    			for (int i = 0; i < rad * 2 + bmap.length; i++) {
    				// System.out.println(sofar + ":" + added.size());
    				if (i < bmap.length) {
    					added.add(bmap[i][j]);
    					sofar += bmap[i][j];
    				} else {
    					if (i < bmap.length + rad) {
    						added.add(bmap[bmap.length - 1][j]);
    						sofar += bmap[bmap.length - 1][j];
    					}
    				}
    				if (i >= rad) {
    					if (i - rad < bmap.length) {
    						bmap[i - rad][j] = sofar;
    					}
    					sofar -= added.get(0);
    					added.remove(0);
    				}
    			}
    		}
    	}

		public BitmapBuilder remove(List<Shape> toRemove) {	
			this.toRemove=toRemove;
			return this;
		}
    }

    /**
     * bounding box shape
     */
    public static enum Bbox {
        Rectangular,
            Polygon,
            DoublePolygon
            }

    static final int[] MASK = new int[]{
        0x80,
        0x40,
        0x20,
        0x10,
        0x08,
        0x04,
        0x02,
        0x01
    };


    
    public List<int[]> findHollowPoints(){
    	return IntStream.range(0, width)
    			.mapToObj(i->IntStream.range(0, height).mapToObj(j->new int[]{i,j}))
    			.flatMap(t->t)
    	    .filter(xy->this.get(xy[0]-1, xy[1]) &&
    	    		    this.get(xy[0], xy[1]-1) &&
    	    		    this.get(xy[0]+1, xy[1]) &&
    	    		    this.get(xy[0], xy[1]+1)
    	    		)
    	    .collect(Collectors.toList());
    }
    
    
    
    
    static class SparseArray<T> {
        Map<Integer, Map<Integer, T>> data =
            new HashMap<Integer, Map<Integer, T>> ();

        T set (int i, int j, T v) {
            Map<Integer, T> m = data.get (i);
            if (m == null) {
                data.put (i, m = new HashMap<Integer, T> ());
            }
            T x = m.put (j, v);
            return x;
        }

        T get (int i, int j) {
            Map<Integer, T> m = data.get (i);
            if (m != null) {
                return m.get (j);
            }
            return null;
        }
    }

    public enum ChainCode {
        E (1, 0, '0'), // 0
            NE (1, -1, '1'), // 1
            N (0, -1, '2'), // 2
            NW (-1, -1, '3'), // 3
            W (-1, 0, '4'), // 4
            SW (-1, 1, '5'), // 5
            S (0, 1, '6'), // 6
            SE (1, 1, '7'); // 7

        final int dx, dy;
        final char ch;

        ChainCode (int dx, int dy, char ch) {
            this.dx = dx;
            this.dy = dy;
            this.ch = ch;
        }

        public int dx () {
            return dx;
        }

        public int dy () {
            return dy;
        }

        // angle (in radians) measured ccw
        public double angle () {
            if (dy == 0 && dx == 0) return 0.;
            if (dy == 1 && dx == 1) return Math.PI / 4;
            if (dy == 1 && dx == 0) return Math.PI / 2;
            if (dy == 1 && dx == -1) return 3 * Math.PI / 4;
            if (dy == 0 && dx == -1) return Math.PI;
            if (dy == -1 && dx == -1) return 3 * Math.PI / 2;
            if (dy == -1 && dx == 0) return 5 * Math.PI / 4;
            if (dy == -1 && dx == 1) return 7 * Math.PI / 4;
            return -1.;
        }

        public char ch () {
            return ch;
        }
    }

    /**
     * this only generate for the first connected component
     * found in the bitmap.  The chain code is based on
     * the following 8-neighbor definition:
     * 3 2 1
     * 4 * 0
     * 5 6 7
     */
    public static class ChainCodeSequence {
        Point2D start; // starting x & y
        List<ChainCode> codes = new ArrayList<ChainCode> ();
        LinkedList<Point2D> coords = new LinkedList<Point2D> ();

        public ChainCodeSequence (int x, int y) {
            start = new Point (x, y);
            coords.add (start);
        }

        // return the new coordinate correspond to this
        public Point2D add (ChainCode code) {
            Point2D pt = coords.getLast ();
            //logger.info(pt.toString()+" "+code);
            Point newPt = new Point ((int) (pt.getX () + code.dx () + .5),
                                     (int) (pt.getY () + code.dy () + .5));
            if (!contains (newPt)) {
                coords.add (newPt);
                codes.add (code);
            } else {
                newPt = null;
            }
            return newPt;
        }

        public boolean contains (double x, double y) {
            for (Point2D pt : coords) {
                if (Math.abs (pt.getX () - x) < EPS
                    && Math.abs (pt.getY () - y) < EPS) {
                    return true;
                }
            }
            return false;
        }

        public boolean contains (Point2D pt) {
            return contains (pt.getX (), pt.getY ());
        }

        public int getStartX () {
            return (int) start.getX ();
        }

        public int getStartY () {
            return (int) start.getY ();
        }

        public Point2D getStartPt () {
            return start;
        }

        public int length () {
            return codes.size ();
        }

        public Point2D[] getCoords () {
            return coords.toArray (new Point2D[0]);
        }

        public ChainCode getCode (Point2D pt) {
            return getCode (pt.getX (), pt.getY ());
        }

        public ChainCode getCode () { // last code
            if (codes.isEmpty ())
                return null;
            return codes.get (codes.size () - 1);
        }

        public ChainCode getCode (double x, double y) {
            return getCode ((int) x, (int) y);
        }

        public ChainCode getCode (int x, int y) {
            int xi = (int) start.getX (), yi = (int) start.getY ();
            for (Iterator<ChainCode> it = codes.iterator (); it.hasNext (); ) {
                ChainCode c = it.next ();

                if (xi == x && yi == y) {
                    return c;
                }

                xi += c.dx ();
                yi += c.dy ();
            }
            return null;
        }

        public ChainCode[] getSequence () {
            return codes.toArray (new ChainCode[0]);
        }

        /**
         * identify dominant points from a chain code sequence
         */
        static class AEV implements Comparable<AEV> {
            int k;
            double dist;

            AEV () {
            }

            AEV (int k, double dist) {
                this.k = k;
                this.dist = dist;
            }

            AEV (int k) {
                this.k = k; // chain code index
            }

            public int compareTo (AEV x) {
                if (dist < x.dist) return -1;
                if (dist > x.dist) return 1;
                return 0;
            }
        }

        ;

        public Point2D[] dominantPoints (double threshold) {
            // break points are candiates for dominant points
            List<AEV> breaks = new ArrayList<AEV> ();
            Point2D[] cc = getCoords ();

            int i = 1;
            for (int j = 0; i < codes.size (); ++i, ++j) {
                if (codes.get (i) != codes.get (j)) {
                    breaks.add (new AEV (i));
                }
            }

            boolean open = cc[i].getX () != cc[0].getX ()
                || cc[i].getY () != cc[0].getY ();
            if (open) {
                // if not closed curve, add the starting and terminating
                //   points
                breaks.add (0, new AEV (0));
                breaks.add (new AEV (i));
            }

            calcAEV (breaks, cc, threshold, !open);

            if (DEBUG) {
                System.out.println
                    ("## " + breaks.size () + " dominant points!");
                for (i = 0; i < breaks.size (); ++i) {
                    AEV aev = breaks.get (i);
                    Point2D pt = cc[aev.k];
                    System.out.println
                        ("** dominant point at " + pt + "; aev = " + aev.dist);
                    if (i + 1 < breaks.size ()) {
                        for (int k = aev.k; k <= breaks.get (i + 1).k; ++k) {
                            System.out.println ("  ++ " + cc[k]);
                        }
                    }
                }
            }

            Point2D[] pts = new Point2D[breaks.size ()];
            for (i = 0; i < breaks.size (); ++i) {
                AEV aev = breaks.get (i);
                pts[i] = cc[aev.k];
            }

            return pts;
        }

        // calculate associated error value for each break point
        void calcAEV (List<AEV> breaks, Point2D[] cc,
                      double threshold, boolean closed) {
            if (breaks.size () < 3) {
                return;
            }

            if (DEBUG) {
                System.out.println ("## " + breaks.size () + " break points!");
            }

            AEV min = new AEV (-1, Double.MAX_VALUE);
            for (int i = 0; i < breaks.size (); ++i) {
                AEV aev = calcAEV (i, breaks, cc, closed);
                if (aev.dist < min.dist) {
                    min.dist = aev.dist;
                    min.k = i;
                }

                if (DEBUG) {
                    Point2D pt = cc[aev.k];
                    System.out.println
                        ("** break point at " + pt + "; aev = " + aev.dist);
                    if (i + 1 < breaks.size ()) {
                        for (int k = aev.k; k <= breaks.get (i + 1).k; ++k) {
                            System.out.println ("  ++ " + cc[k]);
                        }
                    }
                }
            }

            while (min.k >= 0 && min.dist <= threshold) {
                //logger.info("removing break "+min.k+" "+min.dist);
                breaks.remove (min.k);

                min.dist = Double.MAX_VALUE;
                min.k = -1;
                for (int i = 0; i < breaks.size (); ++i) {
                    AEV b = calcAEV (i, breaks, cc, closed);
                    if (min.k < 0 || b.dist < min.dist) {
                        min.dist = b.dist;
                        min.k = i;
                    }
                }
                //logger.info("next min "+min.k+" "+min.dist);
            }
        }


        AEV calcAEV (int i, List<AEV> breaks, Point2D[] cc, boolean closed) {
            int size = breaks.size ();
            AEV aev = breaks.get (i);
            Point2D pt = cc[aev.k];
            aev.dist = Double.MAX_VALUE;
            if (closed && (i == 0 || (i + 1) == size)) {
                if (i == 0) {
                    Point2D pi = cc[breaks.get (size - 1).k];
                    Point2D pj = cc[breaks.get (i + 1).k];
                    aev.dist = sqDist (pt, pi, pj);
                } else { // i+1 == breaks.size()
                    Point2D pi = cc[breaks.get (i - 1).k];
                    Point2D pj = cc[breaks.get (0).k];
                    aev.dist = sqDist (pt, pi, pj);
                }
            } else if (i > 0 && (i + 1) < size) {
                Point2D pi = cc[breaks.get (i - 1).k];
                Point2D pj = cc[breaks.get (i + 1).k];
                aev.dist = sqDist (pt, pi, pj);
            }
            return aev;
        }

        // calculate squared distance from pk to line pj-pi
        static double sqDist (Point2D pk, Point2D pi, Point2D pj) {
            double a = (pk.getX () - pi.getX ()) * (pj.getY () - pi.getY ())
                - (pk.getY () - pi.getY ()) * (pj.getX () - pi.getX ());
            double b = (pi.getX () - pj.getX ());
            double c = (pi.getY () - pj.getY ());
            return a * a / (b * b + c * c);
        }

        public String toString () {
            StringBuilder sb = new StringBuilder
                (getClass () + "{x=" + start.getX () + ",y=" + start.getY ());
            if (!codes.isEmpty ()) {
                sb.append (",length=" + codes.size () + ",");
                for (ChainCode c : codes) {
                    sb.append (c.ch ());
                }
                sb.deleteCharAt (sb.length () - 1);
            }
            sb.append ("}");
            return sb.toString ();
        }
    }

    private byte[] data; // pixel values

    
    
    
    private int width, height;
    private int scanline;
    private SampleModel sampleModel;
    
    private Map<Integer, Integer> scanlineCache = new ConcurrentHashMap<>();
    
    private int getScanlineFor(int y){
    	return scanlineCache.computeIfAbsent(y, k->scanline*k);
    }
    
    private static byte[] sqrtCache = CachedSupplier.of(()->{
    	int maxReal = 2*Byte.MAX_VALUE*Byte.MAX_VALUE; // this is like 127*127*2 = 15-bits ~ 32kB of 4x sqrts for each input
    	byte[] cache= new byte[maxReal];
    	
    	for(int i=0;i<maxReal;i++){
    		cache[i] = (byte)Math.min(Byte.MAX_VALUE, Math.round(Math.sqrt(i)*4));
    	}
    	return cache;
    }).get();
    
    private CachedSupplier<byte[]> distanceData = CachedSupplier.of(()->{
    	byte[] distanceX=new byte[data.length*8];
    	byte[] distanceY=new byte[data.length*8];
    	
    	
    	//TODO may want to fiddle with this number/algorithm
    	int growUntil = 2;
    	
    	
    	int nscan=scanline*8;

    	
    	for (int y = 0; y < this.height; ++y) {
             for (int x = 0; x < this.width; ++x) {
            	 int loc = y * nscan + x;
            	 if(get(x,y)){
            		 distanceX[loc] =0;
            		 distanceY[loc] =0;
            	 }else{
            		 
            		 distanceX[loc] =Byte.MAX_VALUE;
            		 distanceY[loc] =Byte.MAX_VALUE;
            	 }
             }
        }
    	
    	BiFunction<Integer,Integer,Integer> translate = (x,y)->y*nscan + x;
    	
    	int[] dx = new int[]{-1, 0, 1,-1,1,-1,0,1};
    	int[] dy = new int[]{-1,-1,-1, 0,0, 1,1,1};
    	
    	boolean changedSomething=true;
    	for (int r = 0; changedSomething && r<growUntil;r++){
    		changedSomething=false;
	    	for (int y = 0; y < this.height; ++y) {
	    		int yoff=y*nscan;
	    		
	            for (int x = 0; x < this.width; ++x) {
	          		 int loc = yoff + x;
	          		 
	          		 if(distanceX[loc]==Byte.MAX_VALUE && distanceY[loc] == Byte.MAX_VALUE){
	          			 //It's unset
	          			 int mini=-1;
	          			 byte ndx=0;
	          			 byte ndy=0;
	     				 double minsqdist=Byte.MAX_VALUE*Byte.MAX_VALUE+Byte.MAX_VALUE*Byte.MAX_VALUE;
	          			 for(int i=0;i<dx.length;i++){
	          				 int nx = dx[i]+x;
	          				 int ny = dy[i]+y;
	          				 if(nx<this.width && nx>=0 && ny<this.height && ny>=0){
	          					int nloc=translate.apply(nx, ny);
	          					double bdx=distanceX[nloc]*1.0 + Math.abs(dx[i]);
	          					double bdy=distanceY[nloc]*1.0 + Math.abs(dy[i]);
	          					
	          					double d = bdx*bdx+bdy*bdy;
	          					if(d<minsqdist){
	          						minsqdist=d;
	          						mini=i;
	          						ndx=(byte) Math.min(Byte.MAX_VALUE, bdx);
	          						ndy=(byte) Math.min(Byte.MAX_VALUE, bdy);
	          					}
	          				 }	 
	          			 }
	          			 if(mini>=0){
	          				changedSomething=true;
	          				distanceX[loc] =ndx;
	      					distanceY[loc] =ndy;
	          			 }
	          		 }
	            }
	       }
    	}
    	
    	
    	
    	for (int y = 0; y < this.height; ++y) {
    		int yoff=y*nscan;
    		
            for (int x = 0; x < this.width; ++x) {
          		 int loc = yoff + x;
          		 byte ddx=distanceX[loc];
         		 byte ddy=distanceY[loc];
          		 if(ddx>=Byte.MAX_VALUE || ddy>=Byte.MAX_VALUE){
          			distanceX[loc] = (byte)Byte.MAX_VALUE;
          			continue;
          		 }          		 
          		 int dd=ddx*ddx+ddy*ddy;
          		 distanceX[loc]=sqrtCache[dd];
            }
       }
    	
    	
    	return distanceX;
    }); //distance to nearest pixel*4
    
    
    public static Bitmap createBitmap (Raster raster, double sigma) {
        SampleModel model = raster.getSampleModel();
        int band = model.getNumBands ();
        if (band > 1) {
            throw new IllegalArgumentException
                ("Can't handle sample with multiple channels");
        }


        ImageStats[] is = new ImageStats[]{null};


        Bitmap bm= new SigmaThreshold (sigma).binarize(raster, stat->{
        	is[0]=stat;
        });




        double tPct=100*(is[0].threshold-is[0].min)/(is[0].max-is[0].min);
        double tPctMinSig=100*((is[0].threshold-is[0].stdev)-is[0].min)/(is[0].max-is[0].min);
        double tPctMaxSig=100*((is[0].threshold+is[0].stdev)-is[0].min)/(is[0].max-is[0].min);

        double count = 0;
        double countOn = 0;
        for(int i=(int)Math.max(1, tPctMinSig);i<=Math.min(tPctMaxSig, 100);i++){
        	count+=is[0].histogram[i];
        }
        for(int i=(int)Math.max(1, tPct);i<=100;i++){
        	countOn+=is[0].histogram[i];
        }

//print histogram
//      for(int i=0;i<is[0].histogram.length;i++){
//  		System.out.println(i + "\t" + is[0].histogram[i]);
//  	}
//  	System.out.println("Threshold:" + tPct);
//  	System.out.println("Min Threshold:" + tPctMinSig);
//  	System.out.println("Max Threshold:" + tPctMaxSig);
//
        //If there's a little uncertainty about where to draw the threshold line, try
        //the adaptive threshold, possibly
        if(count> countOn*0.1 || count> (is[0].count-countOn)*0.1){

//        	System.out.println("testing adaptive");
            List<Shape> polys2= bm.connectedComponents(Bitmap.Bbox.DoublePolygon);

            if(polys2.size()<4000){

    	        Bitmap bm1= new AdaptiveThreshold().binarize(raster);
    	        List<Shape> polys1= bm1.connectedComponents(Bitmap.Bbox.DoublePolygon);

    	        if(polys1.size()<4000){
    	            long sum1=polys1.stream()
    		              .mapToLong(s->polys1.stream().filter(s2->GeomUtil.contains(s, s2)).count())
    		              .sum();
    		        long sum2=polys2.stream()
    		                .mapToLong(s->polys2.stream().filter(s2->GeomUtil.contains(s, s2)).count())
    		                .sum();
//    		        //if there are at least 3 more shapes inside other shapes, it's
//    		        //probably a thresholding issue that should use the one with more shapes
//    		        //The logic here is that aromatic double bonds are quite common, and if
//    		        //the thresholding washes them out, then you'd expect to see about 3 or more shapes
//    		        //inside other shapes with good thresholding than with bad thresholding
    		        if(sum1>=sum2+3){
//    		        	System.out.println("Using adaptive");
    		        	return bm1;
    		        }
    	        }
    	        //System.out.println("time:" + (System.currentTimeMillis()-tstart));
            }
        }



        
        return bm;
        
    }




    
    public Bitmap clean(){
    	double pon=fractionPixelsOn();
    	
    	if(pon>.5){
    		return this.invert();
    	}else{
    		return this;
    	}
    }
    
    private CachedSupplier<Double> fractionOn = CachedSupplier.of(()->{
    	return this._fractionPixelsOn();
    });
    
    public double fractionPixelsOn(){
    	return this.fractionOn.get();
    }
    private double _fractionPixelsOn(){
       int on=0;
       for (int i = 0; i < data.length; ++i) {
  		 on+= Integer.bitCount((data[i] & 0xff));
      }
       return ((double)on) / (data.length*8);
    }
    
    
    
    
    public Bitmap invert(){
    	Bitmap clone = new Bitmap(this);
    	for (int i = 0; i < clone.data.length; ++i) {
    		 clone.data[i] = (byte) (~clone.data[i] & 0xff);
        }
    	return clone;
    }

    public static Bitmap read (byte[] file, double sigma) throws IOException {

               return createBitmap (ImageUtil.grayscale(file).getData(),sigma);

    }
    public static Bitmap read (byte[] file) throws IOException {
    	return read(file,1.2);
    }
    public static Bitmap read (File file, double sigma) throws IOException {
            return createBitmap (ImageUtil.grayscale(file).getData(), sigma);

    }
    
    public static Bitmap read (BufferedImage bi, double sigma) {
    	 return createBitmap (ImageUtil.grayscale(bi).getData(),sigma);
    }
    
    




    // create an empty image
    public Bitmap (Bitmap copy) {
        this (copy.width, copy.height);
        System.arraycopy (copy.data, 0, this.data, 0, this.data.length);
    }

    public Bitmap (int width, int height) {
        this.width = width;
        this.height = height;

        scanline = (width + 7) >> 3;
        data = new byte[scanline * height];
        sampleModel = new MultiPixelPackedSampleModel
            (DataBuffer.TYPE_BYTE, width, height, 1, scanline, 0);
    }
    
    public Bitmap(){
    	
    }

    public Object clone () {
        return new Bitmap (this);
    }

    public int width () {
        return width;
    }

    public int height () {
        return height;
    }

    public int getAsInt (int x, int y) {
    	return get(x,y) ? 1: 0;
    }
    public boolean get (int x, int y) {
        if (x >= 0 && x < width && y >= 0 && y < height) {
            return isOn(x,y);
        }
        return false;
    }

    // same as get() but without the bound checking
    public boolean isOn (int x, int y) {
        return (data[getScanlineFor(y) + x / 8] & MASK[x % 8]) != 0;
    }

    public void set (int x, int y, boolean on) {
        int loc = getScanlineFor(y) + x / 8;
        if (on) {
            data[loc] |= MASK[x % 8];
        } else {
            data[loc] &= ~MASK[x % 8];
        }
    }


    /*
     * 8-neighbor of p
     *   p(7)  p(0)  p(1)
     *   p(6)   p    p(2)
     *   p(5)  p(4)  p(3)
     */
    public int p0 (int x, int y) {
        return y > 0 && isOn (x, y-1) ? 1 : 0;
    }

    public int p1 (int x, int y) {
        return y > 0 && (x+1 < width) && isOn (x+1, y-1) ? 1 : 0;
    }

    public int p2 (int x, int y) {
        return x+1 < width && isOn (x+1, y) ? 1 : 0;
    }

    public int p3 (int x, int y) {
        return (y+1 < height) && (x+1 < width) && isOn (x+1, y+1) ? 1 : 0;
    }

    public int p4 (int x, int y) {
        return (y+1 < height) && isOn (x, y+1) ? 1 : 0;
    }

    public int p5 (int x, int y) {
        return x > 0 && (y+1 < height) && isOn (x-1, y+1) ? 1 : 0;
    }

    public int p6 (int x, int y) {
        return x > 0 && isOn (x-1, y) ? 1 : 0;
    }

    public int p7 (int x, int y) {
        return x > 0 && y > 0 && isOn (x-1, y-1) ? 1 : 0;
    }

    /*
     * count number of 8-neighbor pixels
     */
    public int neighbor8 (int x, int y) {
        int nb = 0;
        if (p0 (x, y) == 1) ++nb;
        if (p1 (x, y) == 1) ++nb;
        if (p2 (x, y) == 1) ++nb;
        if (p3 (x, y) == 1) ++nb;
        if (p4 (x, y) == 1) ++nb;
        if (p5 (x, y) == 1) ++nb;
        if (p6 (x, y) == 1) ++nb;
        if (p7 (x, y) == 1) ++nb;
        return nb;
    }

    /*
     * number of 8-neighbor pixels that transition from off to on
     */
    public int transition8 (int x, int y) {
        int nb = 0;
        if (p0 (x, y) - p7 (x, y) == 1) ++nb;
        if (p1 (x, y) - p0 (x, y) == 1) ++nb;
        if (p2 (x, y) - p1 (x, y) == 1) ++nb;
        if (p3 (x, y) - p2 (x, y) == 1) ++nb;
        if (p4 (x, y) - p3 (x, y) == 1) ++nb;
        if (p5 (x, y) - p4 (x, y) == 1) ++nb;
        if (p6 (x, y) - p5 (x, y) == 1) ++nb;
        if (p7 (x, y) - p6 (x, y) == 1) ++nb;
        return nb;
    }

    public void dump (OutputStream os) {
        PrintStream ps = new PrintStream (os, true);
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                ps.print (get (x, y) ? '*' : '.');
            }
            if (y % 10 == 0) {
                ps.print (" " + y);
            }
            ps.println ();
        }
    }

    public SampleModel getSampleModel () {
        return sampleModel;
    }

    public WritableRaster createRaster () {
        WritableRaster raster =
            Raster.createWritableRaster (sampleModel, null);
        for (int y = 0; y < height; ++y) {
            int band = getScanlineFor(y);
            for (int x = 0; x < width; ++x) {
                // the default IndexColorModel is 0 for black and 1 white
                raster.setSample
                    (x, y, 0, (data[band + x / 8] & MASK[x % 8]) == 0 ? 1 : 0);
            }
        }
        return raster;
    }

    public BufferedImage createBufferedImage () {
        BufferedImage image = new BufferedImage
            (width, height, BufferedImage.TYPE_BYTE_BINARY);
        image.setData (createRaster ());
        return image;
    }

    public boolean write (String format, OutputStream os) throws IOException {
        return ImageIO.write (createBufferedImage (), format, os);
    }

    public boolean write (String format, File output) throws IOException {
        return ImageIO.write (createBufferedImage (), format, output);
    }

    public boolean write (File output) throws IOException {
        return write ("png", output);
    }

    public boolean write (OutputStream os) throws IOException {
        return write ("png", os);
    }
    
    
    private CachedSupplier<List<int[]>> onInts = CachedSupplier.of(()->{
    	List<int[]> on = new ArrayList<>();
    	for(int i=0;i<this.width();i++){
    		for(int j=0;j<this.height();j++){
    			if(this.get(i, j)){
    				on.add(new int[]{i,j});
    			}
    		}
    	}
    	return on;
    });
    
    public Stream<int[]> getXYOnPoints(){
    	return onInts.get().stream();
    }

    public Bitmap crop (Shape s) {
        Rectangle r = s.getBounds ();
        if (r.width == 0 || r.height == 0) {
            return null;
        }

        Bitmap dst = new Bitmap (r.width + 1, r.height + 1);
        int x1 = Math.min (width, r.x + r.width);
        int y1 = Math.min (height, r.y + r.height);
        int x0 = r.x, y0 = r.y;
        
        int area=r.width*r.height;
        
        //if the area is large, it's usually easier to deal with the sparse array
        //instead of a bit-for-but copy
        if(area > 500){
	        getXYOnPoints().filter(xy->(xy[0]>=x0 && xy[0]<=x1) && (xy[1]>=y0 && xy[1]<=y1))
	        			   .filter(xy->xy[0]==x1 || xy[1] == y1 || s.contains(xy[0],xy[1]))
	        			   
	        			   .map(xy->new int[]{xy[0]-x0,xy[1]-y0})
	        			   .forEach(xy->{
	        				   
	        				   dst.set(xy[0], xy[1], true);
	        			   });
        }else{
	        int i, j = 0;
	        for (int y = y0; y <= y1; ++y, ++j) {
	            i = 0;
	            for (int x = x0; x <= x1; ++x, ++i) {
	                if (x == x1 || y == y1 || s.contains (x, y))
	                    dst.set (i, j, get (x, y));
	            }
	        }
        }
        return dst;
    }
    
    

    public Bitmap crop (int x, int y, int w, int h) {
        Bitmap dst = new Bitmap (w, h);
        int x1 = Math.min (width, x + w);
        int y1 = Math.min (height, y + h);
        // this is pretty slow... we should be using a lut
        //   (look-up-table) here
        int i, j = 0;
        for (int y0 = y; y0 < y1; ++y0, ++j) {
            i = 0;
            for (int x0 = x; x0 < x1; ++x0, ++i)
                dst.set (i, j, get (x0, y0));
        }
        return dst;
    }
    
    public long onPixelsInShape(Shape s){
    	Rectangle r = s.getBounds ();
        if (r.width == 0 || r.height == 0) {
            return 0;
        }

        
        int x1 = Math.min (width, r.x + r.width);
        int y1 = Math.min (height, r.y + r.height);
        int x0 = r.x, y0 = r.y;
        long total =0;
        int i, j = 0;
        for (int y = y0; y <= y1; ++y, ++j) {
            i = 0;
            for (int x = x0; x <= x1; ++x, ++i) {
                if (x == x1 || y == y1 || s.contains (x, y)){
                	if(get(x,y))total++;
                }
                	
                    
            }
        }
        return total;
    }

    /* Thinning algorithm based on Nagendraprasad, Wang, and Gupta.  The 
       following description is based on 
       Gonzalez and Woods, Digital Image Processing, Addison Wesley, 1992.
       
       The algorithm is as follows:  First, all contour pixels are identified: 
       a pixel is a contour pixel if it's a foreground pixel and one of its 
       8-connected neighbors is a background pixel.  Then for each contour 
       pixel, the following two steps are iteratively applied to the image 
       until no change occurs, in which case we're done.
       
       Step 1:  All contour pixels satisfying the following conditions are
       marked for deletion.
       
       (a) 2 <= N(p1) <= 6
       (b) S(p1) = 1
       (c) p2 * p4 * p6 = 0
       (d) p4 * p6 * p8 = 0
       
       where N(p1) is the number of 8-connected foreground pixels around p1.  
       S(p1) is the number of 0-1 transitions from p2, p3, p4, ..., p8, p9, p2
       according to the following labeling scheme.
       p9  p2  p3
       p8  p1  p4
       p7  p6  p5
       After all contour pixels have been processed, those that were marked for
       deletion are deleted from the image.  Next, step 2 is applied to the 
       image.
       
       Step 2:  This step is almost identical to step 1.  The only difference 
       here is that conditions (c) and (d) are changed to
       
       (c') p2 * p4 * p8 = 0
       (d') p2 * p6 * p8 = 0
       
       respectively.  Both steps 1 and 2 are repeated until there is no change
       in the image. */
    public Bitmap skeleton () {
        int N, S;
        boolean changed;

        Bitmap thin = new Bitmap (this);
        byte[] copy = new byte[this.data.length];
        System.arraycopy (this.data, 0, copy, 0, this.data.length);

        int[] p = new int[10];
        boolean flag, step = false;
        while (true) {
            changed = false;

            for (int y = 0; y < height; ++y)
                for (int x = 0; x < width; ++x) {
                    /* only process contour pixels */
                    if (thin.get (x, y)
                        && ((x == 0 || y == 0      /* boundary */
                             || (y + 1) == height || (x + 1) == width)
                            /* checking 8-neighbor of white pixel */
                            || !thin.get(x - 1, y) || !thin.get(x + 1, y)
                            || !thin.get(x, y - 1) || !thin.get(x, y + 1)
                            || !thin.get(x - 1, y - 1) 
                            || !thin.get(x + 1, y - 1)
                            || !thin.get(x - 1, y + 1) 
                            || !thin.get(x + 1, y + 1))) {
                        /* count the number of pixels in the mask */
                        p[2] = p[3] = p[4] = p[5] =
                            p[6] = p[7] = p[8] = p[9] = 0;
                        N = 0;
                        if (x > 0) {
                            N += p[8] = thin.get(x - 1, y) ? 1 : 0;
                            N += p[9] = thin.get(x - 1, y - 1) ? 1 : 0;
                            N += p[7] = thin.get(x - 1, y + 1) ? 1 : 0;
                        }
                        if (x + 1 < width) {
                            N += p[4] = thin.get(x + 1, y) ? 1 : 0;
                            N += p[3] = thin.get(x + 1, y - 1) ? 1 : 0;
                            N += p[5] = thin.get(x + 1, y + 1) ? 1 : 0;
                        }
                        N += p[2] = thin.get(x, y - 1) ? 1 : 0;
                        N += p[6] = thin.get(x, y + 1) ? 1 : 0;

                        /* count the number of 0-1 transition */
                        S = p[3] - p[2] == 1 ? 1 : 0;
                        S += p[4] - p[3] == 1 ? 1 : 0;
                        S += p[5] - p[4] == 1 ? 1 : 0;
                        S += p[6] - p[5] == 1 ? 1 : 0;
                        S += p[7] - p[6] == 1 ? 1 : 0;
                        S += p[8] - p[7] == 1 ? 1 : 0;
                        S += p[9] - p[8] == 1 ? 1 : 0;
                        S += p[2] - p[9] == 1 ? 1 : 0;

                        /* step 1 or step 2 and does the proper connectivity
                           checking */
                        flag = step
                            ? ((p[2] & p[4] & p[8]) == 0) 
                            && ((p[2] & p[6] & p[8]) == 0)
                            : ((p[2] & p[4] & p[6]) == 0) 
                            && ((p[4] & p[6] & p[8]) == 0);

                        if ((N >= 2 && N <= 6) && S == 1 && flag) {
                            /* flag this pixel to be deleted */
                            copy[getScanlineFor(y) + x / 8] &= ~MASK[x % 8];
                            changed = true;
                        }
                    } /* endif boundary pixel */
                } /* endfor x */

            /* check if there is any change; if there isn't we're done */
            if (!changed)
                break;
            /* copy the image back */
            System.arraycopy (copy, 0, thin.data, 0, copy.length);
            step = !step; // toggle the step
        } /* endwhile (1) */

        return thin;
    }

    /**
     * This version is a slight improvement to the NWG algorithm above. 
     * It's based on the following paper:
     * R. Carrsco, M. Forcada, A note on the Nagendraprasad-Wang-Gupta
     * thinning algorithm, Pattern Recognition Letters, 16, 539-541, 1995.
     */
    public Bitmap thin () {
        Bitmap thin = new Bitmap (this);
        byte[] copy = new byte[this.data.length];
        System.arraycopy (thin.data, 0, copy, 0, copy.length);

        int parity = 1;
        boolean changed;
        
        do {
            changed = false;
            parity = 1 - parity;

            for (int y = 0; y < height; ++y)
                for (int x = 0; x < width; ++x) {
                    if (thin.isOn(x, y)) {
                        int nb = thin.neighbor8(x, y);
                        int ap = thin.transition8(x, y);
                        int cp = ((thin.p0(x, y) == 0 
                                   && thin.p1(x, y) == 0
                                   && thin.p2(x, y) == 0
                                   && thin.p5(x, y) == 0
                                   && thin.p4(x, y) == 1
                                   && thin.p6(x, y) == 1)
                                  || (thin.p2(x, y) == 0
                                      && thin.p3(x, y) == 0 
                                      && thin.p4(x, y) == 0
                                      && thin.p7(x, y) == 0
                                      && thin.p6(x, y) == 1
                                      && thin.p0(x, y) == 1)) ? 1 : 0; 
                        int dp = ((thin.p1(x, y) == 0 
                                   && thin.p4(x, y) == 0
                                   && thin.p5(x, y) == 0
                                   && thin.p6(x, y) == 0
                                   && thin.p0(x, y) == 1
                                   && thin.p2(x, y) == 1)
                                  || (thin.p0(x, y) == 0 
                                      && thin.p3(x, y) == 0
                                      && thin.p6(x, y) == 0
                                      && thin.p7(x, y) == 0
                                      && thin.p2(x, y) == 1 
                                      && thin.p4(x, y) == 1)) ? 1 : 0;
                        if (nb > 1 && nb < 7 
                            && (ap == 1 || ((1-parity)*cp + parity*dp) == 1)) {
                            int ep = (thin.p2(x, y) + thin.p4(x, y))
                                *thin.p0(x, y) * thin.p6(x, y);
                            int fp = (thin.p6(x, y) + thin.p0(x, y))
                                *thin.p4(x, y) * thin.p2(x, y);
                            if ((parity == 0 && ep == 0) 
                                || (parity == 1 && fp == 0)) {
                                // delete this pixel
                                copy[getScanlineFor(y) + x / 8] &= ~MASK[x % 8];
                                changed = true;
                            }
                        }
                    } // if pixel is on
                } // endfor each pixel

            // update the image
            if (changed) {
                System.arraycopy (copy, 0, thin.data, 0, copy.length);
            }
        }
        while (changed);
        copy = null;

        return thin;
    }

    void union (short[] eqvtab, short cls1, short cls2) {
        short i = cls1, j = cls2, k;
        //logger.info("union "+cls1+" "+cls2);

        while (eqvtab[i] > 0) i = eqvtab[i];
        while (eqvtab[j] > 0) j = eqvtab[j];

        while (eqvtab[cls1] > 0) {
            k = cls1;
            cls1 = eqvtab[cls1];
            eqvtab[k] = i;
        }

        while (eqvtab[cls2] > 0) {
            k = cls2;
            cls2 = eqvtab[cls2];
            eqvtab[k] = j;
        }

        if (i != j) {
            if (eqvtab[j] < eqvtab[i]) {
                eqvtab[j] += eqvtab[i] - 1;
                eqvtab[i] = j;
            } else {
                eqvtab[i] += eqvtab[j] - 1;
                eqvtab[j] = i;
            }
        }
    }

    /*
     * return connected components as rectangular bounding boxes
     */
    public List<Shape> rectConnectedComponents () {
        return connectedComponents (Bbox.Rectangular);
    }

    /*
     * return connected components as convex hull polygons
     */
    public List<Shape> polyConnectedComponents () {
        return connectedComponents (Bbox.Polygon);
    }

    public List<Shape> connectedComponents () {
        return connectedComponents (Bbox.Rectangular);
    }

    private Map<Bbox,List<Shape>> _cacheShapes = new HashMap<>();



    public List<Shape> connectedComponents (Bbox shape) {
    	return _cacheShapes.computeIfAbsent(shape, (ss)->{
    		short label = 0; // current label
            short[][] labels = new short[height][width + 1];

            // equivalence class
            short[] eqvtab = new short[500]; // some initial default
            short[] L = new short[4];

            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    /* check to see if we have a black pixel */
                    if (!get (x, y))
                        /* do nothing */ ;
                    /* boundary conditions */
                    else if (y == 0 && x == 0) {
                        labels[y][x] = ++label;
                    } else if (y == 0) {
                        short label1 = labels[y][x - 1];
                        if (label1 == 0) {
                            label1 = ++label;
                        }
                        labels[y][x] = label1;
                    } else if (x == 0) {
                        int label1 = labels[y - 1][x];
                        int label2 = labels[y - 1][x + 1];
                        if (label1 != 0 && label2 != 0)
                            label1 = Math.min (label1, label2);
                        else if (label1 == 0 && label2 == 0) {
                            label1 = ++label;
                        } else
                            label1 = Math.max (label1, label2);
                        labels[y][x] = (short) label1;
                    }
                    /* assign new label */
                    else if (labels[y][x - 1] == 0
                             && labels[y - 1][x] == 0
                             && labels[y - 1][x - 1] == 0
                             && labels[y - 1][x + 1] == 0) {
                        labels[y][x] = ++label;
                    } else {
                        L[0] = labels[y - 1][x - 1];
                        L[1] = labels[y - 1][x];
                        L[2] = labels[y - 1][x + 1];
                        L[3] = labels[y][x - 1];

                        Arrays.sort (L);

                        /* skip all non-labeled pixels */
                        int n;
                        for (n = 0; n < 4 && L[n] == 0; ++n)
                            ;
                        /* n should not be 4 */
                        if (n == 4) {
                            throw new IllegalStateException ("n == 4");
                        }

                        labels[y][x] = L[n];
                        /* now enumerate from n to 4 - 1 */
                        for (int i = n; i < 4; ++i)
                            for (int j = i + 1; j < 4; ++j)
                                /* update equivalence table */
                                union (eqvtab, L[i], L[j]);
                    }

                    if (label == Short.MAX_VALUE) {
                        logger.log (Level.SEVERE, "Max number of labels reached: "
                                    + label + "; truncating search!");
                        break;
                    }
                    // ensure there's enough space in the eqvtab
                    else if (label >= eqvtab.length) {
                        short[] newtab = new short[label + 100];
                        System.arraycopy (eqvtab, 0, newtab, 0, eqvtab.length);
                        eqvtab = newtab;
                    }
                }
            }

            if (DEBUG) {
                System.err.print ("eqvtab:");
                for (int i = 1; i <= label; ++i) {
                    System.err.print (" " + i + ":" + eqvtab[i]);
                }
                System.err.println ();
                System.err.println ("label: " + label);

                System.err.println ("eqv class labels...");
                for (int y = 0; y < height; ++y) {
                    for (int x = 0; x < width; ++x) {
                        System.err.print
                            (get (x, y) ? String.valueOf (labels[y][x]) : '.');
                    }
                    System.err.println ();
                }
            }


            List<Shape> comps;
            switch (shape) {
            case Polygon:
                comps = connectedComponentPolygonShapes (eqvtab, labels);
                break;
            case DoublePolygon:
                comps = connectedComponentDoublePrecisionPolygonShapes (eqvtab, labels);
                break;

            case Rectangular:
            default:
                comps = connectedComponentRectangularShapes (eqvtab, labels);
            }

            if (DEBUG) {
                System.err.println ("merged labels...");
                for (int y = 0; y < height; ++y) {
                    for (int x = 0; x < width; ++x) {
                        System.err.print
                            (get (x, y) ? String.valueOf (labels[y][x]) : '.');
                    }
                    System.err.println ();
                }
            }
            eqvtab = null;
            labels = null;

            return comps;
    	});

    }


    List<Shape> connectedComponentPolygonShapes
        (short[] eqvtab, short[][] labels) {

        Map<Short, List<Point>> coords = new LinkedHashMap<Short, List<Point>> ();
        for (int y = 0; y < height; ++y)
            for (int x = 0; x < width; ++x) {
                short label = labels[y][x];
                if (label != 0) {
                    short l = label;
                    /* find equivalence class */
                    while (eqvtab[l] > 0)
                        l = eqvtab[l];

                    labels[y][x] = l;

                    List<Point> pts = coords.get (l);
                    if (pts == null) {
                        coords.put (l, pts = new ArrayList<Point> ());
                    }
                    pts.add (new Point (x, y));
                }
            }

        List<Shape> comps = new ArrayList<Shape> ();
        for (List<Point> pts : coords.values ()) {
            Polygon hull = GeomUtil.convexHullOldIntPrecision (pts.toArray (new Point[0]));
            comps.add (hull);
        }

        return comps;
    }
    List<Shape> connectedComponentDoublePrecisionPolygonShapes
	    (short[] eqvtab, short[][] labels) {
	
	    Map<Short, List<Point>> coords = new LinkedHashMap<Short, List<Point>> ();
	    for (int y = 0; y < height; ++y)
	        for (int x = 0; x < width; ++x) {
	            short label = labels[y][x];
	            if (label != 0) {
	                short l = label;
	                /* find equivalence class */
	                while (eqvtab[l] > 0)
	                    l = eqvtab[l];
	
	                labels[y][x] = l;
	
	                List<Point> pts = coords.get (l);
	                if (pts == null) {
	                    coords.put (l, pts = new ArrayList<Point> ());
	                }
	                pts.add (new Point (x, y));
	            }
	        }
	
	    List<Shape> comps = new ArrayList<Shape> ();
	    for (List<Point> pts : coords.values ()) {
	    	Point2D[] ptsadjusted=pts.stream()
	    	    .flatMap(pt->{
	    	    	return Stream.of(new Point2D.Double(pt.getX(),pt.getY())
	    	    					 //,new Point2D.Double(pt.getX()+1,pt.getY()),
	    	    					 //,new Point2D.Double(pt.getX()+1,pt.getY()+1),
	    	    					 //,new Point2D.Double(pt.getX(),pt.getY()+1
	    	    			);
	    	    })
	    	    .toArray(i->new Point2D[i]);
	    	Shape hull = GeomUtil.convexHull2 (ptsadjusted);
	    	
	        comps.add (hull);
	    }
	
	    return comps;
	}

    List<Shape> connectedComponentRectangularShapes
        (short[] eqvtab, short[][] labels) {

        Map<Short, Rectangle> ltab = new LinkedHashMap<Short, Rectangle> ();
        List<Shape> comps = new ArrayList<Shape> ();

        for (int y = 0; y < height; ++y)
            for (int x = 0; x < width; ++x) {
                short label = labels[y][x];
                if (label != 0) {
                    short l = label;
                    /* find equivalence class */
                    while (eqvtab[l] > 0)
                        l = eqvtab[l];

                    labels[y][x] = l;
                    /* create bounding box for each class and make
                       sure that it does not go outside of the image
                       boundary */
                    Rectangle r = ltab.get (l);
                    if (r == null) {
                        ltab.put (l, r = new Rectangle (x, y, 1, 1));
                        comps.add (r);
                    }
                    int x0 = Math.min (r.x, x);
                    int y0 = Math.min (r.y, y);
                    int x1 = Math.min (width, Math.max (r.x + r.width, x + 1));
                    int y1 = Math.min (height, Math.max (r.y + r.height, y + 1));
                    r.setBounds (x0, y0, x1 - x0, y1 - y0);
                }
            }

        return comps;
    }

    static EnumSet<ChainCode> getNeighbors (Bitmap b, int x, int y) {
        EnumSet<ChainCode> Nb = EnumSet.noneOf (ChainCode.class);

        // boundary cases
        if (x == 0 && y == 0) { // Nb in {0, 6, 7}
            if (b.get (x + 1, y)) Nb.add (ChainCode.E);
            if (b.get (x, y + 1)) Nb.add (ChainCode.S);
            if (b.get (x + 1, y + 1)) Nb.add (ChainCode.SE);
        } else if (x == 0) { // Nb in {2,1,0,7,6}
            if (b.get (x + 1, y)) Nb.add (ChainCode.E);
            if (b.get (x + 1, y - 1)) Nb.add (ChainCode.NE);
            if (b.get (x, y - 1)) Nb.add (ChainCode.N);
            if (b.get (x, y + 1)) Nb.add (ChainCode.S);
            if (b.get (x + 1, y + 1)) Nb.add (ChainCode.SE);
        } else if (y == 0) { // Nb in {4,5,6,7,0}
            if (b.get (x - 1, y)) Nb.add (ChainCode.W);
            if (b.get (x - 1, y - 1)) Nb.add (ChainCode.SW);
            if (b.get (x, y + 1)) Nb.add (ChainCode.S);
            if (b.get (x + 1, y + 1)) Nb.add (ChainCode.SE);
            if (b.get (x + 1, y)) Nb.add (ChainCode.E);
        } else if (x == b.width - 1 && y == b.height - 1) {
            // Nb in {2,3,4}
            if (b.get (x, y - 1)) Nb.add (ChainCode.N);
            if (b.get (x - 1, y - 1)) Nb.add (ChainCode.NW);
            if (b.get (x - 1, y)) Nb.add (ChainCode.W);
        } else if (x == b.width - 1) {
            // Nb in {2,3,4,5,6}
            if (b.get (x, y - 1)) Nb.add (ChainCode.N);
            if (b.get (x - 1, y - 1)) Nb.add (ChainCode.NW);
            if (b.get (x - 1, y)) Nb.add (ChainCode.W);
            if (b.get (x - 1, y + 1)) Nb.add (ChainCode.SW);
            if (b.get (x, y + 1)) Nb.add (ChainCode.S);
        } else if (y == b.height - 1) {
            // Nb in {0,1,2,3,4}
            if (b.get (x + 1, y)) Nb.add (ChainCode.E);
            if (b.get (x + 1, y - 1)) Nb.add (ChainCode.NE);
            if (b.get (x, y - 1)) Nb.add (ChainCode.N);
            if (b.get (x - 1, y - 1)) Nb.add (ChainCode.NW);
            if (b.get (x - 1, y)) Nb.add (ChainCode.W);
        } else {
            if (b.get (x + 1, y)) Nb.add (ChainCode.E);
            if (b.get (x + 1, y - 1)) Nb.add (ChainCode.NE);
            if (b.get (x, y - 1)) Nb.add (ChainCode.N);
            if (b.get (x - 1, y - 1)) Nb.add (ChainCode.NW);
            if (b.get (x - 1, y)) Nb.add (ChainCode.W);
            if (b.get (x - 1, y + 1)) Nb.add (ChainCode.SW);
            if (b.get (x, y + 1)) Nb.add (ChainCode.S);
            if (b.get (x + 1, y + 1)) Nb.add (ChainCode.SE);
        }
        return Nb;
    }

    public List<ChainCodeSequence> chainCodes () {
        return chainCodes (5);
    }

    public List<ChainCodeSequence> chainCodes (int minsize) {
        Bitmap clone = new Bitmap (this);

        List<ChainCodeSequence> seqs = new ArrayList<ChainCodeSequence> ();
        for (ChainCodeSequence seq; (seq = chainCode (clone)) != null; ) {
            if (DEBUG) {
                System.out.println ("-- " + seq);
                for (int y = 0; y < clone.height (); ++y) {
                    for (int x = 0; x < clone.width (); ++x) {
                        String s = clone.get (x, y) ? "*" : ".";
                        ChainCode c = seq.getCode (x, y);
                        if (c != null) {
                            s = "" + c.ch ();
                        }
                        System.out.print (s);
                    }
                    System.out.println ();
                }
                seq.dominantPoints (DEFAULT_AEV_THRESHOLD);
            }
            if (seq.length () >= minsize) {
                seqs.add (seq);
            }
        }

        return seqs;
    }

    public List<Point2D> dominantPoints () {
        return dominantPoints (5, DEFAULT_AEV_THRESHOLD);
    }

    public List<Point2D> dominantPoints (int minsize, double threshold) {
        List<Point2D> dps = new ArrayList<Point2D> ();
        for (ChainCodeSequence ccs : chainCodes (minsize)) {
            for (Point2D pt : ccs.dominantPoints (threshold)) {
                dps.add (pt);
            }
        }
        return dps;
    }

    
    
    
    public List<Path2D> segments () {
        return segments (1, DEFAULT_AEV_THRESHOLD);
    }
    
    //0 means not wedge like
    //+1 means perfect wedge-like from A to B, with B being wider
    //-1 means perfect wedge-like from B to A, with A being wider
    public double getWedgeLikeScore(Line2D line){
    	
    	double sx=line.getX1();
		double sy=line.getY1();
		double dx=line.getX2()-line.getX1();
		double dy=line.getY2()-line.getY1();
		
    	double len=GeomUtil.length(line);
    	if(len<1)return 0;
    	double mult=1/len;
    	
    	int widthDistance=(int)(Math.round(len/4));
    	
    	double stepX=dx*mult;
    	double stepY=dy*mult;
    	
    	int[] c = new int[(int)len];
		for(int d=0;d<c.length;d++){
			double ddx = stepX*d+sx;
			double ddy = stepY*d+sy;
			
			for(int i=-widthDistance;i<widthDistance;i++){
				double iddx = i*stepY+ddx;
				double iddy = -i*stepX+ddy;
				if(this.get((int)Math.round(iddx), (int)Math.round(iddy))){
					c[d]++;
				}
			}
		}
		int[] rc = Arrays.stream(c).filter(cr->cr>0).toArray();
		
		return GeomUtil.ordinalCorrel(rc);
    }
    
    
    public static class WedgeInfo{
    	private Shape hull;
    	int onPixels;
		double area;
    	double correl;
    	Line2D line;
    	Line2D longSplit;
    	double padding;
    	
    	public Line2D getLine(){
    		return line;
    	}
    	public int getOnPixels() {
			return onPixels;
		}

		public void setOnPixels(int onPixels) {
			this.onPixels = onPixels;
		}

		public double getArea() {
			return area;
		}

		public void setArea(double area) {
			this.area = area;
		}

		public double getCorrel() {
			return correl;
		}

		public void setCorrel(double correl) {
			this.correl = correl;
		}

    	
    	public WedgeInfo(Shape s, int pix, double area, double c){
    		this.setHull(s);
    		this.onPixels=pix;
    		this.area=area;
    		this.correl=c;
    	}

		public Shape getHull() {
			return hull;
		}

		public void setHull(Shape hull) {
			this.hull = hull;
			this.longSplit=GeomUtil.findLongestSplittingLine(hull).getLine();
		}
    	
		
		public double getAverageThickness(){
			return this.onPixels/(GeomUtil.length(this.line)-padding*2);
		}
		
		
		public double pctOfHull(){
			return GeomUtil.length(this.longSplit)/GeomUtil.length(this.line);
		}
    	
    }
    
    public Optional<WedgeInfo> getconfexHullAlongLine(Line2D line){
    	
    	double sx=line.getX1();
		double sy=line.getY1();
		double dx=line.getX2()-line.getX1();
		double dy=line.getY2()-line.getY1();
		
    	double len=GeomUtil.length(line);
    	if(len<1)return Optional.empty();
    	double mult=1/len;
    	
    	int widthDistance=(int)(Math.round(len/4));
    	
    	List<Point2D> pts= new ArrayList<>();
    	
    	double stepX=dx*mult;
    	double stepY=dy*mult;
    	
    	int pad=(int)(len/6);
    	
    	int c=0;
    	int[] cl = new int[(int)Math.ceil(len-2*pad)];
    	
    	int k=0;
    	
    	for(int d=pad;d<len-pad;d++){
			double ddx = stepX*d+sx;
			double ddy = stepY*d+sy;
			int f=0;
			
			int fx=(int)Math.round(ddx);
	    	int fy=(int)Math.round(ddy);
	    	int lx=fx;
	    	int ly=fy;
	    	
	    	boolean found1=false;
	    	
	    	if(this.get(fx, fy)){
	    		pts.add(new Point2D.Double(fx, fy));
	    		c++;
	    		f++;
	    	}
	    	
			for(int i=1;i<widthDistance;i++){
				double iddx = i*stepY+ddx;
				double iddy = -i*stepX+ddy;
				//TODO
				if(this.get((int)Math.round(iddx), (int)Math.round(iddy))){
					c++;
					f++;
					fx=(int)Math.round(iddx);
					fy=(int)Math.round(iddy);
					if(!found1){
						lx=fx;
						ly=fy;
					}
					found1=true;
				}else{
					break;
				}
			}
			for(int i=-1;i>-widthDistance;i--){
				double iddx = i*stepY+ddx;
				double iddy = -i*stepX+ddy;
				//TODO
				if(this.get((int)Math.round(iddx), (int)Math.round(iddy))){
					c++;
					f++;
					lx=(int)Math.round(iddx);
					ly=(int)Math.round(iddy);
					if(!found1){
						fx=lx;
						fy=ly;
					}
					found1=true;
				}else{
					break;
				}
			}
			if(found1){
				pts.add(new Point2D.Double(fx, fy));
				pts.add(new Point2D.Double(fx+1, fy));
				pts.add(new Point2D.Double(fx, fy+1));
				pts.add(new Point2D.Double(fx+1, fy+1));
				pts.add(new Point2D.Double(lx, ly));
				pts.add(new Point2D.Double(lx+1, ly));
				pts.add(new Point2D.Double(lx, ly+1));
				pts.add(new Point2D.Double(lx+1, ly+1));
				cl[k]+=f;
			}
			k++;
		}
		Shape shull=GeomUtil.convexHull2(pts.stream().toArray(i->new Point2D[i]));
		if(shull==null)return Optional.empty();
		double area=GeomUtil.area(shull);
		if(area<1)return Optional.empty();
		double correl = GeomUtil.ordinalCorrel(cl);
		
		WedgeInfo wi = new WedgeInfo(shull,c,area,correl);
		wi.line=line;
		wi.padding=pad;
		
		return Optional.of(wi);
    }
    
    //0 means not wedge like
    //+1 means perfect wedge-like from A to B, with B being wider
    //-1 means perfect wedge-like from B to A, with A being wider
    public double getDashLikeScore(Line2D line){
    	
    	double sx=line.getX1();
		double sy=line.getY1();
		double dx=line.getX2()-line.getX1();
		double dy=line.getY2()-line.getY1();
		
    	double len=GeomUtil.length(line);
    	if(len<1)return 0;
    	double mult=1/len;
    	
    	int widthDistance=(int)(Math.round(len/4));
    	
    	double stepX=dx*mult;
    	double stepY=dy*mult;
    	
    	int[] c = new int[(int)len];
		for(int d=0;d<c.length;d++){
			double ddx = stepX*d+sx;
			double ddy = stepY*d+sy;
			
			for(int i=-widthDistance;i<widthDistance;i++){
				double iddx = i*stepY+ddx;
				double iddy = -i*stepX+ddy;
				if(this.get((int)Math.round(iddx), (int)Math.round(iddy))){
					c[d]++;
				}
			}
		}
		
		return Math.sqrt(GeomUtil.variance(c))/len;
    }
    
    public double getLineLikeScore(Line2D line){
    	byte[] distMet=distanceData.get();
    	BiFunction<Double,Double,Double> sample = (x,y)->{
    		//TODO: do interp eventually
    		int ix=(int)Math.round(x);
    		int iy=(int)Math.round(y);
    		int ni=getScanlineFor(iy)*8+ix;
    		if(ni>distMet.length || ni<0)return (double)Byte.MAX_VALUE/4;
    		return (double) (distMet[ni]/4);    		
    	};
    	double sx=line.getX1();
		double sy=line.getY1();
		double dx=line.getX2()-line.getX1();
		double dy=line.getY2()-line.getY1();
		
    	double len=GeomUtil.length(line);
    	double mult=1/len;
    	double sumDist = 0;
		for(int d=0;d<len;d++){
			double ddx = mult*d*dx+sx;
			double ddy = mult*d*dy+sy;
			double dist=sample.apply(ddx, ddy);
			sumDist+=dist;
		}
		return sumDist/len;
    }
    
    //TODO: worry about this
    public static int MAX_REPS=1000;
    
    
    public List<LineWrapper> combineLines(List<LineWrapper> ilines, double maxMinDistance, double maxAvgDeviation, double maxDistanceToConsiderSamePoint, double maxAngle, double minLengthForAngleCompare){
    	List<LineWrapper> lines = ilines;
    	int[] reps1=new int[]{0};
    	
    	boolean gotOne=true;
    	while(gotOne){
    		gotOne=false;
    		List<LineWrapper> nlines=combineLines2(lines, maxMinDistance, maxAvgDeviation, maxDistanceToConsiderSamePoint,maxAngle,minLengthForAngleCompare,reps1);
    		if(nlines.size()!=lines.size()){
    			gotOne=true;
    			lines=nlines;
    		}
    		if(reps1[0]>MAX_REPS)break;
    	}
    	return lines;
    }
    
    private List<LineWrapper> combineLines2(List<LineWrapper> ilines, double maxMinDistance, double maxAvgDeviation, double maxDistanceToConsiderSamePoint,double maxAngle, double minLengthForAngleCompare, int[] reps){
    	byte[] distMet=distanceData.get();
    	
    	List<LineWrapper> lines=ilines.stream()
    	     .sorted()
    	     .collect(Collectors.toList());
    	
    	List<Point2D> allPoints=lines.stream()
       	     .flatMap(l->l.streamPoints())
       	     .collect(Collectors.toList());
    	
    	List<Point2D> dontmerge=GeomUtil.groupPointsCloserThan(allPoints,maxDistanceToConsiderSamePoint)
    	        .stream()
    	        .filter(pL->pL.size()>2)
    	        .flatMap(l->l.stream())
    	        .collect(Collectors.toList());
    	
    	
    	
    	BiFunction<Double,Double,Double> sample = (x,y)->{
    		//do interp eventually
    		int ix=(int)Math.round(x);
    		int iy=(int)Math.round(y);
    		return (double) (distMet[getScanlineFor(iy)*8+ix]/4);    		
    	};
    	
    	double maxCosAng = Math.abs(Math.cos(maxAngle));

    	
    	
    	int[] ii = new int[]{0};
    	
		for (int i = 0; i < lines.size(); i++) {
			if (reps[0] >= MAX_REPS)
				break;    
			LineWrapper line1 = lines.get(i);
			ii[0]=i;			
			
			boolean[] foundOne= new boolean[]{false};
			
			IntStream.range(i+1,lines.size())
			         .mapToObj(k->Tuple.of(k,k))
			         .map(Tuple.kmap(k->lines.get(k)))
			         .map(Tuple.kmap(l->LineDistanceCalculator.from(line1.getLine(), l.getLine())))
			         .map(t->t.withKComparatorMap(lu->lu.getAbsoluteClosestDistance()))
			         .sorted()
			         .filter(t->{
			        	 LineWrapper line2=lines.get(t.v());
			        	 double dist=t.k().getAbsoluteClosestDistance();
			        	 //if(alreadyMerged.clines.get(t.v())
			        	
			        	 return dist<maxMinDistance; 
			         })
			         .filter(t->{
			        	 
			        	 LineWrapper line2=lines.get(t.v());
			        	 if(line1.length()>minLengthForAngleCompare && 
			        	    line2.length()>minLengthForAngleCompare){
				        	 if(line1.absCosTheta(line2)<maxCosAng){
				        		 return false;
				        	 }
			        	 }
			        	 return true;
			         })
			         .filter(t->{
			        	 	Point2D[] closest=t.k().closestPoints();
							
							boolean partOfTriple=dontmerge.stream()
							         .flatMap(p->Stream.of(p.distance(closest[0]),p.distance(closest[1])))
							         .filter(d->(d<maxDistanceToConsiderSamePoint))
							         .findAny()
							         .isPresent();
							if(partOfTriple)return false;
							return true;
			         })
			         .filter(t->{
			        	 	int j=t.v();
			        	 	LineDistanceCalculator ldc = t.k();
			        	 	LineWrapper line2 = lines.get(j);
							Line2D combined = ldc.getLineFromFarthestPoints();
							double sx = combined.getX1();
							double sy = combined.getY1();
							double dx = combined.getX2() - combined.getX1();
							double dy = combined.getY2() - combined.getY1();
							double len = GeomUtil.length(combined);

							double mult = 1 / len;

							double sumSqDist = 0;
							for (int d = 0; d < len; d++) {
								double ddx = mult * d * dx + sx;
								double ddy = mult * d * dy + sy;
								double dist = sample.apply(ddx, ddy);
								sumSqDist += dist * dist;
							}
							double sqrtSTDErr = Math.sqrt(sumSqDist / len);
							if (sqrtSTDErr <= maxAvgDeviation) {
								if (len > line1.length() && len > line2.length()) {
									return true;
								}else if(line1.intersectsLine(line2)){
										return true;
								}
							} else {
//								System.out.println("No Good:" + sqrtSTDErr);
							}
							return false;
			         })
			         .findFirst()
			         .ifPresent(t->{
			        	 	int j=t.v();
			        	 	LineDistanceCalculator ldc = t.k();
			        	 	LineWrapper line2 = lines.get(j);
			        	 	Line2D combined = ldc.getLineFromFarthestPoints();
			        	 	
//			        	 	System.out.println("Old line 1:" + LineUtil.length(line1));
//							System.out.println("Old line 2:" + LineUtil.length(line2));
//							System.out.println("New line:" + LineUtil.length(combined));
//							System.out.println("Dist:" + ldc.getSmallestPointDistance());
//							
							lines.set(ii[0], LineWrapper.of(combined));
							lines.remove(j);
							reps[0]++;
							foundOne[0]=true;
//							System.out.println("reps:" + reps[0] + " of " + MAX_REPS);
			         });
			if(foundOne[0])i--;
			
			    
			
		}
	    	
    	
    	return lines;
    }
    

    public List<Path2D> segments (int minsize, double threshold) {
        List<Path2D> segs = new ArrayList<Path2D> ();
        for (ChainCodeSequence ccs : chainCodes (minsize)) {
            GeneralPath gp = null;
            for (Point2D pt : ccs.dominantPoints (threshold)) {
                if (gp == null) {
                    gp = new GeneralPath ();
                    gp.moveTo (pt.getX (), pt.getY ());
                } else {
                    gp.lineTo (pt.getX (), pt.getY ());
                }
            }
            segs.add (gp);
        }
        return segs;
    }

    public static ChainCodeSequence chainCode (Bitmap bitmap) {
        int x = -1, y = -1; // locate the first point

        boolean done = false;
        for (int j = 0; j < bitmap.height; ++j) {
            int band = j * bitmap.scanline;
            for (int i = 0; i < bitmap.width; ++i) {
                if ((bitmap.data[band + i / 8] & MASK[i % 8]) != 0) {
                    x = i;
                    y = j;
                    done = true;
                    break;
                }
            }
            if (done) break;
        }

        if (x < 0 || y < 0) {
            return null;
        }

        ChainCodeSequence seq = new ChainCodeSequence (x, y);

        do {
            EnumSet<ChainCode> Nb = getNeighbors (bitmap, x, y);

            Point2D pt = null;
            if (Nb.isEmpty ()) {
            } else if (Nb.size () == 1) { //
                pt = seq.add (Nb.iterator ().next ());
            } else {
                // multiple choice; pick best one based on the following
                //  rule: select the one for which
                ChainCode best = null;
                EnumSet<ChainCode> bestNq = null;

                for (ChainCode c : Nb) {
                    int xp = x + c.dx (), yp = y + c.dy ();
                    if (!seq.contains (xp, yp)) {
                        EnumSet<ChainCode> Nq = getNeighbors (bitmap, xp, yp);
                        if (bestNq == null
                            || (!Nq.isEmpty () && Nq.size () < bestNq.size ())
                            // choose the least change in direction
                            || (Nq.size () == bestNq.size ()
                                && c.angle () < best.angle ())) {
                            best = c;
                            bestNq = Nq;
                        }
                    }
                }

                if (best != null) {
                    pt = seq.add (best);
                }
            }

            if (pt != null) {
                // continue
                x = (int) pt.getX ();
                y = (int) pt.getY ();
            } else {
                break; // we're done
            }
        }
        while (true);

        /*
         * all pixels connected to any chain code should be deleted
         */
        // remove all the pixels that make up the chain code
        x = seq.getStartX ();
        y = seq.getStartY ();
        for (ChainCode c : seq.getSequence ()) {
            bitmap.set (x, y, false);
            x += c.dx ();
            y += c.dy ();
        }
        bitmap.set (x, y, false);

        // do one more pass to remove all singleton pixels
        //  that are left behind from this chain code
        x = seq.getStartX ();
        y = seq.getStartY ();
        for (ChainCode c : seq.getSequence ()) {
            EnumSet<ChainCode> Nb = getNeighbors (bitmap, x, y);
            for (ChainCode n : Nb) {
                int xp = x + n.dx (), yp = y + n.dy ();
                EnumSet<ChainCode> Nq = getNeighbors (bitmap, xp, yp);
                if (Nq.isEmpty ()) {
                    bitmap.set (xp, yp, false);
                }
            }
            x += c.dx ();
            y += c.dy ();
        }

        return seq;
    }

    public static ChainCodeSequence chainCode2 (Bitmap bitmap) {
        int x = -1, y = -1; // locate the first point

        boolean done = false;
        for (int j = 0; j < bitmap.height; ++j) {
            int band = j * bitmap.scanline;
            for (int i = 0; i < bitmap.width; ++i) {
                if ((bitmap.data[band + i / 8] & MASK[i % 8]) != 0) {
                    x = i;
                    y = j;
                    done = true;
                    break;
                }
            }
            if (done) break;
        }
        Bitmap visited = new Bitmap (bitmap.width, bitmap.height);

        ChainCodeSequence seq = new ChainCodeSequence (x, y);
        chainCode2 (seq, visited, bitmap, x, y);

        return seq;
    }

    public static void chainCode2
        (ChainCodeSequence seq, Bitmap visited, Bitmap bitmap, int x, int y) {
        if (x < 0 || y < 0 || x >= bitmap.width || y >= bitmap.height) {
            return;
        }

        bitmap.set (x, y, false);
        EnumSet<ChainCode> Nb = getNeighbors (bitmap, x, y);

        logger.info ("+ x:" + x + " y:" + y + " N:" + Nb.size ());
        if (Nb.isEmpty ()) {
        } else if (Nb.size () == 1) { //
            Point2D pt = seq.add (Nb.iterator ().next ());
            if (pt != null) {
                chainCode2 (seq, visited, bitmap,
                            (int) pt.getX (), (int) pt.getY ());
            }
        } else {
            // multiple choice; pick best one based on the following
            //  rule: select the one for which
            ChainCode best = null;
            EnumSet<ChainCode> bestNq = null;

            for (ChainCode c : Nb) {
                int xp = x + c.dx (), yp = y + c.dy ();
                if (!seq.contains (xp, yp)) {
                    EnumSet<ChainCode> Nq = getNeighbors (bitmap, xp, yp);
                    if (bestNq == null
                        || (!Nq.isEmpty () && Nq.size () < bestNq.size ())
                        // choose the least change in direction
                        || (Nq.size () == bestNq.size ()
                            && c.angle () < best.angle ())) {
                        best = c;
                        bestNq = Nq;
                    }
                }
            }

            if (best != null) {
                Point2D pt = seq.add (best);
                if (pt != null) {
                    chainCode2 (seq, visited, bitmap,
                                (int) pt.getX (), (int) pt.getY ());
                    bitmap.set ((int) pt.getX (), (int) pt.getY (), true);
                }
            }

            for (ChainCode c : Nb) {
                if (c != best) {
                    Point2D pt = seq.add (c);
                    if (pt != null) {
                        chainCode2 (seq, visited, bitmap,
                                    (int) pt.getX (), (int) pt.getY ());
                    }
                }
            }
        }
        logger.info ("- x:" + x + " y:" + y);
    }


    /**
     * detect line segments in the bitmap using the Hough transform
     * thetaDelta - angle partition in degree
     * rhoDelta - radius partition in degree
     */
    public List<Line2D> segments (int thetaDelta, int rhoDelta) {
        if (thetaDelta <= 0) {
            throw new IllegalArgumentException
                ("Invalid delta value: " + thetaDelta);
        }

        // rho = x cos(theta) + y sin(theta)
        int nsteps = 180 / thetaDelta;
        int rmax = (int) (0.5 + Math.sqrt (width * width + height * height));
        int nrho = 2 * rmax / rhoDelta;
        List[][] H = new List[nsteps + 1][nrho];

        List<List> lines = new ArrayList<List> ();
        for (int y = 0; y < height; ++y) {
            int band = y * scanline;
            for (int x = 0; x < width; ++x) {
                if ((data[band + x / 8] & MASK[x % 8]) != 0) {
                    //System.err.println("x="+x+" y="+y);
                    for (int n = 0; n < nsteps; ++n) {
                        int t = n * thetaDelta;
                        double theta = Math.toRadians (t);
                        int rho = (int) (0.5 + x * Math.cos (theta)
                                         + y * Math.sin (theta));
                        int r = (rho + rmax) / rhoDelta;
                        List l = H[n][r];
                        if (l == null) {
                            H[n][r] = l = new ArrayList ();
                            lines.add (l);
                        }
                        l.add (new Point (x, y));
                    }
                }
            }
        }
        H = null;

        Collections.sort (lines, new Comparator<List> () {
                              public int compare (List l1, List l2) {
                                  return l2.size () - l1.size ();
                              }
                          });

        List<Line2D> segments = new ArrayList<Line2D> ();
        for (List l : lines) {
            if (l.size () < 2) {
                break;
            }

            Point start = null, p = null;

            //System.err.println("** path "+l.size());
            next: for (int i = 1; i < l.size (); ++i) {
                Point pt = (Point) l.get (i);

                //System.out.println("  "+pt);
                if (start == null) {
                    start = pt;
                } else if (GeomUtil.isNeighbor (pt, p)) {

                } else {
                    if (start.distance (p) > 2.) {
                        Line2D ln = new Line2D.Float (start, p);
                        for (Line2D s : segments) {
                            if (ln.intersectsLine (s)) {
                                ln = null;
                                break;
                            }
                        }
                        if (ln != null) {
                            segments.add (ln);
                        }
                    }
                    start = pt;
                }
                p = pt;
            }

            if (start.distance (p) > 2.) {
                Line2D ln = new Line2D.Float (start, p);
                for (Line2D s : segments) {
                    if (ln.intersectsLine (s)) {
                        ln = null;
                        break;
                    }
                }
                if (ln != null) {
                    segments.add (ln);
                }
            }
        }

        logger.info (segments.size () + " segments!");
        for (Line2D l : segments) {
            System.err.println (l.getP1 () + " - " + l.getP2 ()
                                + " length=" + l.getP1 ().distance (l.getP2 ()));
        }

        return segments;
    }



}
