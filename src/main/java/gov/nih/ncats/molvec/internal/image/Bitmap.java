package gov.nih.ncats.molvec.internal.image;

import java.awt.Point;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.Transparency;
import java.awt.color.ColorSpace;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.ComponentColorModel;
import java.awt.image.DataBuffer;
import java.awt.image.DataBufferByte;
import java.awt.image.MultiPixelPackedSampleModel;
import java.awt.image.Raster;
import java.awt.image.RenderedImage;
import java.awt.image.SampleModel;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;
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
import java.util.function.BiFunction;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import javax.imageio.ImageIO;

import gov.nih.ncats.molvec.internal.image.binarization.Binarization;
import gov.nih.ncats.molvec.internal.util.CachedSupplier;
import gov.nih.ncats.molvec.internal.algo.StructureImageExtractor;
import gov.nih.ncats.molvec.internal.algo.Tuple;
import gov.nih.ncats.molvec.internal.image.binarization.AdaptiveThreshold;
import gov.nih.ncats.molvec.internal.image.binarization.ImageStats;
import gov.nih.ncats.molvec.internal.util.GeomUtil;
import gov.nih.ncats.molvec.internal.util.GeomUtil.LineDistanceCalculator;
import gov.nih.ncats.molvec.internal.util.GeomUtil.LineWrapper;
import gov.nih.ncats.molvec.internal.util.GeomUtil.ShapeWrapper;

/**
 * A bitmap image
 */
public class Bitmap implements Serializable {
    private static final long serialVersionUID = 0x5f1f54d8fed49ab3l;
    private static final Logger logger =
        Logger.getLogger (Bitmap.class.getName ());

    private static final double EPS = 0.000001;
    public static final double DEFAULT_AEV_THRESHOLD = 1.5;

    
    
    
    public static class Grid{
    	private static final int MIN_WIDTH=16;
    	int x;
    	int y;
    	int wid;
    	int swid;
    	int count=0;
    	
//    	int minx=-1;
//    	int miny=-1;
//    	int maxx=-1;
//    	int maxy=-1;
    	
    	boolean giveChildren=true;
    	
    	// 0 1
    	// 2 3
    
    	Grid[] children = new Grid[4];
    	
    	public Grid(int x, int y, int wid){
    		this.x=x;
    		this.y=y;
    		this.wid=wid;
    		this.swid=wid/2;
    		
    		if(wid<=MIN_WIDTH){
    			giveChildren=false;
    		}
    	}
    	
    	public Grid add(int x1, int y1, int c){
    		count+=c;
    		if(giveChildren){
    			int idx = getIndex(x1,y1);
    			Grid gg=children[idx];
    			if(gg==null){
    				gg= new Grid(x+(idx/2)*swid,y+(idx%2)*swid,swid);
    				children[idx]=gg;
    			}
    			gg.add(x1, y1,c);    			
    		}
    		return this;
    	}
    	public Grid add(int x1, int y1){
    		return add(x1,y1,1);
    	}
    	
    	public Grid remove(int x1, int y1){
    		return add(x1,y1,-1);
    	}
    	
    	public int getIndex(int x1, int y1){
    		int idx=0;
    		if(x1>=x+swid)idx+=2;
    		if(y1>=y+swid)idx+=1;
    		return idx;
    	}
    	
    	public List<Grid> getLeafGrids(int c){
    		List<Grid> all = new ArrayList<Grid>();
    		addGridsAbove(c,all);
    		return all;
    		
    	}
    	
    	public List<Rectangle> getLeafBoundsInside(int c, int maxx, int maxy){
    		 return this.getLeafGrids(c)
    				 .stream()
    				 .map(tgrid->{
    					int widy = Math.min(tgrid.y+tgrid.wid, maxy);
     	               	int widx = Math.min(tgrid.x+tgrid.wid, maxx);
     	               	return new Rectangle(tgrid.x,tgrid.y,widx,widy);
    				 })
    				 .collect(Collectors.toList());
    	}
    	
    	public boolean addGridsAbove(int c, List<Grid> gridList){
    		if(count>c){
    			boolean added=false;
    			for(int i=0;i<4;i++){
    				Grid child = children[i];
    				if(child!=null){
    					added=child.addGridsAbove(c,gridList) || added;
    				}
    			}
    			if(!added){
    				gridList.add(this);
    			}
    			return true;
    		}    		
    		return false;
    	}
    	
    	public Rectangle2D getBounds(){
    		return new Rectangle2D.Double(x,y, wid, wid);
    	}
    	
    	
    	public Grid invert(){
    		int tot = wid*wid;
    		count = tot -count;
    		for(int i=0;i<4;i++){
				Grid child = children[i];
				if(child!=null){
					child.invert();
				}else if(this.giveChildren){
					children[i]=new Grid(x+(i/2)*swid,y+(i%2)*swid,swid);
					children[i].invert();
				}
			}
    		return this;
    	}
    	
    	public Grid clone(){
    		Grid gclone = new Grid(x,y,wid);
    		gclone.count=this.count;
    		for(int i=0;i<4;i++){
				Grid child = children[i];
				if(child!=null){
					gclone.children[i]=child.clone();
				}
			}
    		return gclone;
    	}
    }
    
    private Grid onGrid=null;
    
    
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
    public enum Bbox {
        Rectangular{
            @Override
            List<Shape> computeConnectedComponentShapes(short[] eqvtab, short[][] labels) {
                Map<Short, Rectangle> ltab = new LinkedHashMap<Short, Rectangle> ();
                List<Shape> comps = new ArrayList<Shape> ();
//labels = new short[height][width + 1];
                int height = labels.length;
                int width = labels[0].length -1;
                for (int y = 0; y < height; ++y) {
                    //TODO is this an off by 1 error since label array goes +1 ?
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
                            Rectangle r = ltab.get(l);
                            if (r == null) {
                                ltab.put(l, r = new Rectangle(x, y, 1, 1));
                                comps.add(r);
                            }
                            int x0 = Math.min(r.x, x);
                            int y0 = Math.min(r.y, y);
                            int x1 = Math.min(width, Math.max(r.x + r.width, x + 1));
                            int y1 = Math.min(height, Math.max(r.y + r.height, y + 1));
                            r.setBounds(x0, y0, x1 - x0, y1 - y0);
                        }
                    }
                }
                return comps;
            }
        },
            Polygon{
                @Override
                List<Shape> computeConnectedComponentShapes(short[] eqvtab, short[][] labels) {
                    Map<Short, List<Point>> coords = new LinkedHashMap<Short, List<Point>> ();
                    //labels = new short[height][width + 1];
                    int height = labels.length;
                    int width = labels[0].length -1;
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
            },
        /**
         * Double Precision Polygon ?
         */
            DoublePolygon{
                @Override
                List<Shape> computeConnectedComponentShapes(short[] eqvtab, short[][] labels) {
                    Map<Short, List<Point>> coords = new LinkedHashMap<Short, List<Point>> ();
                    //labels = new short[height][width + 1];
                    int height = labels.length;
                    int width = labels[0].length -1;
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
                                    return Stream.of(new Point2D.Double(pt.getX(),pt.getY()));

//	    	    	return Stream.of(new Point2D.Double(pt.getX()-0.5,pt.getY()-0.5)
//	    	    					 ,new Point2D.Double(pt.getX()-0.5,pt.getY()+0.5)
//	    	    					 ,new Point2D.Double(pt.getX()+0.5,pt.getY()-0.5)
//	    	    					 ,new Point2D.Double(pt.getX()+0.5,pt.getY()+0.5)
//	    	    			);
                                })
                                .toArray(i->new Point2D[i]);
                        Shape hull = GeomUtil.convexHull2 (ptsadjusted);

                        comps.add (hull);
                    }

                    return comps;
                }
            }
        ;

        abstract List<Shape> computeConnectedComponentShapes(short[] eqvtab, short[][] labels);


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
//    	if(true)return new ArrayList<>();
    	
    	return IntStream.range(1, width-1)
    			.mapToObj(i->IntStream.range(1, height-1).mapToObj(j->new int[]{i,j}))
    			.flatMap(t->t)
    	    .filter(xy->this.isOn(xy[0]-1, xy[1]) &&
    	    		    this.isOn(xy[0], xy[1]-1) &&
    	    		    this.isOn(xy[0]+1, xy[1]) &&
    	    		    this.isOn(xy[0], xy[1]+1) &&
    	    		    !this.isOn(xy[0], xy[1]) 
    	    		)
    	    .collect(Collectors.toList());
    }
    
    
    
    //y+ is down
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
        final double len;
        final double rlen;
        

        ChainCode (int dx, int dy, char ch) {
            this.dx = dx;
            this.dy = dy;
            this.ch = ch;
            len = Math.sqrt(this.dx*this.dx + this.dy*this.dy);
            rlen = 1/len;
        }

        public int dx () {
            return dx;
        }

        public int dy () {
            return dy;
        }

        // angle (in radians) measured ccw
        public double angle () {
        	
        	//cardinal directions
        	if (dy == 0 && dx == 0) return 0.; //not a line, shouldn't happen
        	
        	if (dy == 0 && dx == 1) return 0.; //0 degrees
        	if (dy == 1 && dx == 0) return Math.PI / 2; // 90 degrees
        	if (dy == 0 && dx == -1) return Math.PI; //180 degrees
        	if (dy == -1 && dx == 0) return 3 * Math.PI / 2; //270 degrees
            
        	//diagonals
        	if (dy == 1 && dx == 1) return Math.PI / 4; // 45 degrees
            if (dy == 1 && dx == -1) return 3 * Math.PI / 4; //135 degrees
            if (dy == -1 && dx == -1) return 5 * Math.PI / 4;  //225 degrees
            if (dy == -1 && dx == 1) return 7 * Math.PI / 4;  //315 degrees
            
            return -1.;
        }
        

        public char ch () {
            return ch;
        }
        
        public ChainCode inverse(){
        	return ChainCode.values()[(this.ordinal()+4)%8];
        }
        
        public double cosine(ChainCode cc){
        	return (this.dx*cc.dx + this.dy*cc.dy)*(this.rlen*cc.rlen);
        }
        
        //this is the || priority
        private static final int[] priority = new int[] {0,1,2,3,4,3,2,1};
        
        public int priorityChange(ChainCode cc){
        	//if(true)return priority[cc.ordinal()];
        	int delta= this.ordinal();
        	int ni=(cc.ordinal()-delta+8)%8;
        	return priority[ni];
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
        
        public ChainCode peek(){
        	if(codes.isEmpty())return ChainCode.E;
        	return codes.get(codes.size()-1);
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
                return Integer.compare(this.k,x.k);
            }
        }

        ;

        public Point2D[] dominantPoints (double threshold) {
            // break points are candidates for dominant points
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

            
            //This part can be optimized a little, as the 
            //as it does a few unnecessary recalculations
            while (min.k >= 0 && min.dist <= threshold) {
                //logger.info("removing break "+min.k+" "+min.dist);
                breaks.remove (min.k);
                
                int size = breaks.size();
            	int nii=(min.k)%size;
            	int pii=(min.k-1+size)%size;
                
            	//recalculate neighbor AEVs
            	calcAEV (nii, breaks, cc, closed);
            	calcAEV (pii, breaks, cc, closed);
                
                min.dist = Double.MAX_VALUE;
                min.k = -1;

                for (int i = 0; i < breaks.size (); ++i) {
                	AEV b = breaks.get(i);
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
            
            int nii=(i+1)%size;
            int pii=(i-1+size)%size;
            
            Point2D pt = cc[aev.k];
            aev.dist = Double.MAX_VALUE;
            if ((i > 0 && (i + 1) < size) || closed) {
                Point2D pi = cc[breaks.get (pii).k];
                Point2D pj = cc[breaks.get (nii).k];
                aev.dist = sqDist (pt, pi, pj);
            }
            return aev;
        }

        // calculate squared distance from pk to line pj-pi
        // this is the same as a rejection
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
    

    private CachedSupplier<List<int[]>> onInts = CachedSupplier.of(()->{
    	List<int[]> on = new ArrayList<>(onGrid.count);
        
    	
    	
    	List<Rectangle> lookGrids = onGrid.getLeafBoundsInside(0,width,height);
    	
        for(int i=0;i<lookGrids.size();i++){
               	Rectangle tgrid = lookGrids.get(i);
   	            for (int y = tgrid.y; y < tgrid.height; ++y){
   	                for (int x = tgrid.x; x < tgrid.width; ++x) {
   	                		if(this.isOn(x, y)){
   	                			on.add(new int[]{x,y});
   	                		}
   	                }
   	            }
        }
        on.sort(Comparator.comparing(xy->xy[1]*width+xy[0]));
    	return on;
    });
    
    private int getScanlineFor(int y){
    	return scanline*y;
    }
    
    private static byte[] sqrtCache = CachedSupplier.of(()->{
    	int maxReal = 2*Byte.MAX_VALUE*Byte.MAX_VALUE; // this is like 127*127*2 = 15-bits ~ 32kB of 4x sqrts for each input
    	byte[] cache= new byte[maxReal];
    	
    	for(int i=0;i<maxReal;i++){
    		cache[i] = (byte)Math.min(Byte.MAX_VALUE, Math.round(Math.sqrt(i)*4));
    	}

    	return cache;
    }).get();
    
    
    
    //Prime for optimization
    //This is just a map of the l2 distances from each pixel location to the nearest
    //feature pixel. This is sometimes called the "distance transform", and is useful
    //for a few things like thickening, etc. Here, it's used mostly for giving some 
    //tolerance when walking line segments through the bitmap.
    
    //This algorithm below is a fairly intensive one, and only works up to a fixed
    //radius from feature pixels.
    //A better implementation would likely be:
    //http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.507.6791
    //
    private CachedSupplier<byte[]> distanceData = CachedSupplier.of(()->{
    	//TODO may want to fiddle with this number/algorithm
    	long start=System.currentTimeMillis();
    	int growUntil = 5;

//    	
//    	byte[] distanceX=new byte[width*height];
//    	byte[] distanceY=new byte[width*height];
//    	
//    	
//    	
//    	
//    	int nscan=width;
//
//    	Arrays.fill(distanceX, Byte.MAX_VALUE);
//    	Arrays.fill(distanceY, Byte.MAX_VALUE);
//    	int[] dx = new int[]{-1, 0, 1,-1,1,-1,0,1};
//    	int[] dy = new int[]{-1,-1,-1, 0,0, 1,1,1};
//    	
//    	HashSet<Integer> currentUpdates = new HashSet<>();
//    	for(int[] xy:onInts.get()){
//    		int loc = xy[1] * nscan + xy[0];
//           	distanceX[loc] =0;
//       		distanceY[loc] =0;
//       		for(int i=0;i<dx.length;i++){
//       			int nx = xy[0]+dx[i];
//       			int ny = xy[1]+dy[i];
//       			if(nx>=0 && nx<width && ny>=0 && ny<height){
//	       			int locn = ny * nscan + nx;
//	       			currentUpdates.add(locn);
//       			}
//       		}
//    	}
//    	
//    	
//    	HashSet<Integer> nextUpdates = new HashSet<>();
//       
//    	Stack<int[]> toChange = new Stack<>();
//    	
//    	
//    	double MIN_SQ=Byte.MAX_VALUE*Byte.MAX_VALUE+Byte.MAX_VALUE*Byte.MAX_VALUE;
//    	for (int r = 0; r<growUntil;r++){
//    		int p=r%2;
//    		
//    		HashSet<Integer> cur;
//    		HashSet<Integer> nex;
//    		
//    		
//    		if(p==0){
//    			cur=currentUpdates;
//    			nex=nextUpdates;
//    		}else{
//    			nex=currentUpdates;
//    			cur=nextUpdates;
//    		}
//    		
//    		for(int loc:cur){
//    			if(distanceX[loc]==Byte.MAX_VALUE && distanceY[loc] == Byte.MAX_VALUE){
//   				 int x = loc %nscan;
//   				 int y = loc /nscan;
//   				//It's unset
//         			 int mini=-1;
//         			 byte ndx=0;
//         			 byte ndy=0;
//    				 double minsqdist=MIN_SQ;
//         			 for(int i=0;i<dx.length;i++){
//         				 int nx = dx[i]+x;
//         				 int ny = dy[i]+y;
//         				 if(nx<this.width && nx>=0 && ny<this.height && ny>=0){
//         					int nloc=ny*nscan + nx;
//         					double bdx=distanceX[nloc] + Math.abs(dx[i]);
//         					double bdy=distanceY[nloc] + Math.abs(dy[i]);
//         					
//         					double d = bdx*bdx+bdy*bdy;
//         					if(d<MIN_SQ){
//	          					if(d<minsqdist){
//	          						minsqdist=d;
//	          						mini=i;
//	          						ndx=(byte) Math.min(Byte.MAX_VALUE, bdx);
//	          						ndy=(byte) Math.min(Byte.MAX_VALUE, bdy);
//	          					}
//         					}else{
//         						nex.add(nloc);
//         					}
//         				 }	 
//         			 }
//         			 if(mini>=0){
//     					toChange.add(new int[]{loc,ndx,ndy});
//     					
//         			 }
//   			 }
//    		}
//    		
//    		while(!toChange.empty()){
//    			int[] v=toChange.pop();
//    			distanceX[v[0]] =(byte)v[1];
//				distanceY[v[0]] =(byte)v[2];
//				nex.remove(v[0]);
//				
//    		}
//    		cur.clear();
//    	}
//    	
////    	System.out.println("Total Time:" + (System.currentTimeMillis() - start));
//    	
//    	for (int y = 0; y < this.height; ++y) {
//    		int yoff=y*nscan;
//    		
//            for (int x = 0; x < this.width; ++x) {
//          		 int loc = yoff + x;
//          		 byte ddx=distanceX[loc];
//         		 byte ddy=distanceY[loc];
//          		 if(ddx>=Byte.MAX_VALUE || ddy>=Byte.MAX_VALUE){
//          			distanceX[loc] = (byte)Byte.MAX_VALUE;
//          			continue;
//          		 }          		 
//          		 int dd=ddx*ddx+ddy*ddy;
//          		 distanceX[loc]=sqrtCache[dd];
//            }
//       }
//	BufferedImage gg=getGrayscale(width, distanceX);
//    	
//    	try {
//			ImageIO.write(gg, "PNG", new File("loWARBrute.png"));
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//    	System.out.println("Total Time:" + (System.currentTimeMillis() - start));
//    	return distanceX;
    	
    	
    	//This method should be faster, but isn't yet ... probably because it
    	//has to compute the real thing, not just the first few pixels
    	
    	return getNNPixelMapX(growUntil);
    	
    	
    }); //distance to nearest pixel*4
    

    
    //Based roughly on this publication:
    //http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.507.6791
    public byte[] getNNPixelMapX(int limit){
    	int BAD= Integer.MAX_VALUE;
    	int pow2widt=Integer.highestOneBit(width);
		if(pow2widt<width){
			pow2widt=pow2widt*2;
		}
		int pow2wid = pow2widt;
		
		int nshift = (int)Math.round(Math.log(pow2wid)/Math.log(2.0));
		
    	int[] bestIndex = new int[pow2wid*height];
    	byte[] distanceY=new byte[width*height];
    	
    	IntStream.range(0, height)
//    	   .parallel()
    	   .forEach(y->{
	       		int lOnx = -1;
	       		//each column
	       		for(int x=0;x<width;x++){	
	       			boolean on=isOn(x,y);
	       			int locn = y * pow2wid + x;
	       			if(on){
	       				
	       				bestIndex[locn] = locn;
	       				int upto=x;
	       				int setx=lOnx;
	       				if(lOnx>-1){
	       					upto=((lOnx+x)/2);
	       				}else{
	       					setx=x;
	       				}
	       				if(x>lOnx+1){
	   	    				//update nn for last point
	   	    				for(int i=lOnx+1;i<=upto;i++){
	   	    					int locn2 = (y) * pow2wid + i;
	   	    					bestIndex[locn2] = y*pow2wid +setx;
	   	    				}
	   	    				for(int i=upto+1;i<x;i++){
	   	    					int locn2 = (y) * pow2wid + i;
	   	    					bestIndex[locn2] = locn;
	   	    				}
	       				}
	       				lOnx=x;
	       			}else{
	       				bestIndex[locn] = BAD;
	       			}
	       		}
	       		if(lOnx>-1){
	       			for(int i=lOnx+1;i<width;i++){
	   					int locn2 = (y) * pow2wid + i;
	   					bestIndex[locn2] = y*pow2wid+lOnx;
	   				}
	       		}
    	   });
    	
//
////		//each column
    	IntStream.range(0, width)
//    	.parallel()
 	   .forEach(x->{
 		  List<Double> cutoffs = new ArrayList<>();
 		  List<Integer> locationsX = new ArrayList<>();
 		  List<Integer> locationsY = new ArrayList<>();
 		  
 		  
 		  
 		  
 		  
 		  int startY=0;
 		  
 		  //find first real candidate
 		  for(int y=0;y<height;y++){
 			int locn = y*pow2wid + x;
				
			
 		  	if(bestIndex[locn] != BAD){
 		  		int bx = bestIndex[locn]& (pow2wid - 1); 
				int by = bestIndex[locn]>>nshift; // unnecessary probably
 		  		locationsX.add(bx);
 		  		locationsY.add(by);
 		  		
 		  		cutoffs.add(0.0);
 		  		startY=y+1;
 		  		break;	
 		  		
 		  	}
 		  }
 		  if(!locationsX.isEmpty()){
 			 //construct voronoi
 			  for(int y=startY;y<height;y++){
 	 			int locn = y*pow2wid + x;
 	 			
 	 		  	if(bestIndex[locn] != BAD){
 	 		  		
 	 		  		//need to iterate here until it's okay
 	 		  		int bx = bestIndex[locn]& (pow2wid - 1); 
 	 		  		int by = bestIndex[locn]>>nshift; // unnecessary probably
 			 		int dx2= (x-bx);
 			 		
 			 		
 			 		
 			 		double iy=Double.NEGATIVE_INFINITY;
 			 		
 			 		while(locationsX.size()>0){
	 	 		  		int plocX = locationsX.get(locationsX.size()-1);
	 	 		  		int plocY = locationsY.get(locationsY.size()-1);
	 	 		  		double pcut = cutoffs.get(cutoffs.size()-1);
	 	 		  		double dy = by -plocY;
	 	 		  		int dx1= (x-plocX);
	 	 		  		
	 	 		  		//intersection point
	 	 		  		iy = (dx2*dx2-dx1*dx1+dy*dy) / ((double)2*dy) + plocY;
	 	 		  		
	 	 		  		if(iy<pcut){
	 	 		  			//need to remove and loop
	 	 		  			locationsX.remove(locationsX.size()-1);
	 	 		  			locationsY.remove(locationsY.size()-1);
	 	 		  			cutoffs.remove(cutoffs.size()-1);
	 	 		  		}else{
		 	 		  		break;
	 	 		  		}
 			 		}
 			 		locationsX.add(bx);
 	 		  		locationsY.add(by);
 	 		  		cutoffs.add(iy);
 	 		  	}
 	 		 }
 			 cutoffs.add(Double.POSITIVE_INFINITY);
 			 locationsX.add(Integer.MAX_VALUE);
		  	 locationsY.add(Integer.MAX_VALUE);
 			  			 
 			 int py=0;
 			 double pco=cutoffs.get(0);
 			 for(int ci=1;ci<cutoffs.size();ci++){
 				 double co=Math.min(cutoffs.get(ci),height);
 				 int vx=locationsX.get(ci-1);
	 		  	 int vy=locationsY.get(ci-1);
	 		  	 
 				 int from=(int)Math.round(Math.max(pco, py));
 				 for(int y=from;y<co;y++){
 					 int locn = y*pow2wid + x;
 		 		  	 bestIndex[locn]=vy*pow2wid + vx;
 					 py=y+1;
 				 } 				 
 				 pco=co;
 			 }
 			 
 		  }
 		 for(int y=0;y<height;y++){
 			int locn = y*pow2wid + x;
 			int bx = bestIndex[locn]& (pow2wid - 1); 
		  	int by = bestIndex[locn]>>nshift; 
 			int ddx=Math.abs(bx-x);
   			int ddy=Math.abs(by-y);
   			
   		 	int bloc =width*y+x;
   		 
    		 if(ddx>=Byte.MAX_VALUE || ddy>=Byte.MAX_VALUE){
    			distanceY[bloc] = (byte)Byte.MAX_VALUE;
    			continue;
    		 }          		 
    		 int dd=ddx*ddx+ddy*ddy;
    		 
    		 if(dd < limit*limit*2){
    			distanceY[bloc]=(byte) Math.min(sqrtCache[dd], Byte.MAX_VALUE);	 
    		 }else{
    			distanceY[bloc] = (byte)Byte.MAX_VALUE;
    		 }
 		 }
 		   
 	   });
//
//    	BufferedImage gg=getGrayscale(width, distanceY);
//    	
//    	try {
//			ImageIO.write(gg, "PNG", new File("loWARvoronFull.png"));
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}

    	
    	return distanceY;
    	
    }
    
    
    
    
    /**
    *
    * @param width The image width (height derived from buffer length)
    * @param buffer The buffer containing raw grayscale pixel data
    *
    * @return THe grayscale image
    */
   private static BufferedImage getGrayscale(int width, byte[] buffer) {
       int height = buffer.length / width;
       ColorSpace cs = ColorSpace.getInstance(ColorSpace.CS_GRAY);
       int[] nBits = { 8 };
       ColorModel cm = new ComponentColorModel(cs, nBits, false, true,
               Transparency.OPAQUE, DataBuffer.TYPE_BYTE);
       SampleModel sm = cm.createCompatibleSampleModel(width, height);
       DataBufferByte db = new DataBufferByte(buffer, width * height);
       WritableRaster raster = Raster.createWritableRaster(sm, db, null);
       BufferedImage result = new BufferedImage(cm, raster, false, null);

       return result;
   }
    
    public static Bitmap createBitmap (Raster raster, Binarization bb) {
        SampleModel model = raster.getSampleModel();
        int band = model.getNumBands ();
        if (band > 1) {
            throw new IllegalArgumentException
                ("Can't handle sample with multiple channels");
        }


        ImageStats[] is = new ImageStats[]{null};
        
        
        
        
        
    	Bitmap bm= bb.binarize(raster, null, stat->{
        	is[0]=stat;
        });


    	if(is[0]!=null){
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
	
	        //If there's a little uncertainty about where to draw the threshold line, try
	        //the adaptive threshold, possibly
	        if(count> countOn*0.1 || count> (is[0].count-countOn)*0.1){
	
	            List<Shape> polys2= bm.connectedComponents(Bitmap.Bbox.DoublePolygon);
	
	            if(polys2.size()<4000){
	
	    	        Bitmap bm1= new AdaptiveThreshold().binarize(raster, is[0],(ist)->{});
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
	    		        	return bm1;
	    		        }
	    	        }
	    	        //System.out.println("time:" + (System.currentTimeMillis()-tstart));
	            }
	        }
    	}


        
        return bm;
        
    }


	public static RenderedImage readToImage(File file) throws IOException{

			return ImageUtil.grayscale(file);

	}
	public static RenderedImage readToImage(byte[] file) throws IOException{

			return ImageUtil.grayscale(file);

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
       if(onGrid!=null){
    	   return onGrid.count / (double)(width*height);
       }else{
    	   int on=0;
	       
	       for (int i = 0; i < data.length; ++i) {
	  		 on+= Integer.bitCount((data[i] & 0xff));
	       }
	       return ((double)on) / (data.length*8);
       }
    }
    
    
    
    
    public Bitmap invert(){
    	Bitmap clone = new Bitmap(this);
    	for (int i = 0; i < clone.data.length; ++i) {
    		 clone.data[i] = (byte) (~clone.data[i] & 0xff);
        }
    	clone.onGrid.invert();
    	
    	return clone;
    }


    public static Bitmap read (byte[] file, Binarization bin) throws IOException {

               return createBitmap (ImageUtil.grayscale(file).getData(),bin);

    }
    public static Bitmap read (byte[] file) throws IOException {
    	return read(file,StructureImageExtractor.DEF_BINARIZATION);
    }
    public static Bitmap read (File file) throws IOException {
    	return read(file,StructureImageExtractor.DEF_BINARIZATION);
    }
    public static Bitmap read (File file, Binarization bin) throws IOException {
    	
            return createBitmap (ImageUtil.grayscale(file).getData(), bin);

    }
    
    public static Bitmap read (BufferedImage bi, Binarization bb) {
    	 return createBitmap (ImageUtil.grayscale(bi).getData(),bb);
    }
    
    




    // create an empty image
    public Bitmap (Bitmap copy) {
        this (copy.width, copy.height);
        System.arraycopy (copy.data, 0, this.data, 0, this.data.length);
        this.onGrid=copy.onGrid.clone();
    	
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
    	if(onGrid==null){
    		int pow2wid=Integer.highestOneBit(width);
    		int pow2hit=Integer.highestOneBit(height);
    		int twid = Math.max(pow2wid, pow2hit);
    		if(twid<width || twid<height){
    			twid=twid*2;
    		}
    		
    		onGrid = new Grid(0,0,twid);
    	}
    	
    	
    	
        int loc = getScanlineFor(y) + x / 8;
        boolean wasOn = ((data[loc] & MASK[x % 8]) !=0);
        if(on!=wasOn){
	        if(on){
	    		onGrid.add(x, y,1);
	    	}else{
	    		onGrid.add(x, y,-1);
	    	}
        }
        
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
    
    public int neighbor8Index(int x, int y){
    	int nb = 0;
        if (p0 (x, y) == 1) nb|=1;
        if (p1 (x, y) == 1) nb|=2;
        if (p2 (x, y) == 1) nb|=4;
        if (p3 (x, y) == 1) nb|=8;
        if (p4 (x, y) == 1) nb|=16;
        if (p5 (x, y) == 1) nb|=32;
        if (p6 (x, y) == 1) nb|=64;
        if (p7 (x, y) == 1) nb|=128;
        return nb;
    }

    public boolean shouldThin(int nindex, int parity){
    	if(parity==0){
    		return thinCache1[nindex];
    	}else{
    		return thinCache2[nindex];
    	}
    }
    
    
    
    
    private static final boolean[] thinCache1 = new boolean[256];
    private static final boolean[] thinCache2 = new boolean[256];
    
    static{
    	 for(int i=0;i<256;i++){
    		 thinCache1[i] =rule1(i);
    		 thinCache2[i] =rule2(i);
    	 }
    }
    
    
   	
	private static boolean rule1(int b){
		int tot = Integer.bitCount(b);
		if(tot>=2 && tot<=6){
			boolean[] on = new boolean[8];
			for(int i=0;i<8;i++){
				on[i] = ((b & (1<<i))!=0)?true:false;
			}
			boolean ex1 = !on[0] && !on[1] && !on[2] && !on[5] && on[4] && on[6];
			boolean ex2 = !on[2] && !on[3] && !on[4] && !on[7] && on[6] && on[0];
			
			
			
			if(countTransition(b)==1 || ex1 || ex2){				
				
				return ((!on[2] || !on[0] || !on[6])) && ((!on[4] || !on[0] || !on[6]));
				
			}
		}
		return false;
		
	}
	
	private static boolean rule2(int b){
		int tot = Integer.bitCount(b);
		if(tot>=2 && tot<=6){
			boolean[] on = new boolean[8];
			for(int i=0;i<8;i++){
				on[i] = ((b & (1<<i))!=0)?true:false;
			}
		     
//          if (nb > 1 && nb < 7 
//              && (ap == 1 || ((1-parity)*cp + parity*dp) == 1)) {
//              int ep = (thin.p2(x, y) + thin.p4(x, y))
//                  *thin.p0(x, y) * thin.p6(x, y);
//              int fp = (thin.p6(x, y) + thin.p0(x, y))
//                  *thin.p4(x, y) * thin.p2(x, y);
			boolean ex1 = !on[1] && !on[4] && !on[5] && !on[6] && on[0] && on[2];
			boolean ex2 = !on[0] && !on[3] && !on[6] && !on[7] && on[2] && on[4];
			
			
			
			if(countTransition(b)==1 || ex1 || ex2){				
				return ((!on[2] || !on[4] || !on[6])) && ((!on[4] || !on[0] || !on[2]));
			}
		}
		return false;
		
	}
	
	public static int countTransition(int b){
		int tc = 0;
		int[] on = new int[8];
		for(int i=0;i<8;i++){
			on[i] = ((b & (1<<i))!=0)?1:0;
		}
		
		for(int i=0;i<8;i++){
			if(on[(i+1)%8] >on[i]){
				tc++;
			}
		}
		return tc;
		
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

    /**
     * This version is a slight improvement to the NWG algorithm. 
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
        
        Grid gg = thin.onGrid;
        List<Rectangle> lookGrids = gg.getLeafBoundsInside(0,width,height);
        do {
        	changed = false;
            parity = 1 - parity;
            
            
            for(int i=0;i<lookGrids.size();i++){
            	Rectangle tgrid = lookGrids.get(i);
            	for (int y = tgrid.y; y < tgrid.height; ++y){
	                for (int x = tgrid.x; x < tgrid.width; ++x) {
	                    if (thin.isOn(x, y)) {
	                    	int ni = thin.neighbor8Index(x, y);
	                    	
	                    	boolean should = thin.shouldThin(ni, parity);
	                        if(should){
	                        	copy[getScanlineFor(y) + x / 8] &= ~MASK[x % 8];
                                gg.remove(x, y);
                                changed = true;
	                        }
	                        
	                        
	                        
	                    } // if pixel is on
	                } // endfor each pixel
	            }
        	}
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
    		short[] label = new short[]{0}; // current label
            final short[][] labels = new short[height][width + 1];

            // equivalence class
            short[][] eqvtab = new short[1][]; 
            short[] L = new short[4];
            
            eqvtab[0] = new short[500]; //some initial default

            List<int[]> xys = onInts.get();
            
            for(int p=0;p<xys.size();p++){
            	int[] xy = xys.get(p);
            	int x = xy[0];
            	int y = xy[1];
            	if (y == 0 && x == 0) {
                    labels[y][x] = ++label[0];
                } else if (y == 0) {
                    short label1 = labels[y][x - 1];
                    if (label1 == 0) {
                        label1 = ++label[0];
                    }
                    labels[y][x] = label1;
                } else if (x == 0) {
                    int label1 = labels[y - 1][x];
                    int label2 = labels[y - 1][x + 1];
                    if (label1 != 0 && label2 != 0)
                        label1 = Math.min (label1, label2);
                    else if (label1 == 0 && label2 == 0) {
                        label1 = ++label[0];
                    } else
                        label1 = Math.max (label1, label2);
                    labels[y][x] = (short) label1;
                }
                /* assign new label */
                else if (labels[y][x - 1] == 0
                         && labels[y - 1][x] == 0
                         && labels[y - 1][x - 1] == 0
                         && labels[y - 1][x + 1] == 0) {
                    labels[y][x] = ++label[0];
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
                            union (eqvtab[0], L[i], L[j]);
                }

                if (label[0] == Short.MAX_VALUE) {
                    logger.log (Level.SEVERE, "Max number of labels reached: "
                                + label + "; truncating search!");
                    break;
                }
                // ensure there's enough space in the eqvtab
                else if (label[0] >= eqvtab[0].length) {
                    short[] newtab = new short[label[0] + 100];
                    System.arraycopy (eqvtab[0], 0, newtab, 0, eqvtab[0].length);
                    eqvtab[0] = newtab;
                }
            }
            

            if (DEBUG) {
                System.err.print ("eqvtab:");
                for (int i = 1; i <= label[0]; ++i) {
                    System.err.print (" " + i + ":" + eqvtab[0][i]);
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


            List<Shape> comps = ss.computeConnectedComponentShapes(eqvtab[0], labels);

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
//            eqvtab = null;
//            labels = null;

            return comps;
    	});

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
			this.longSplit=ShapeWrapper.of(hull).findLongestSplittingLine().getLine();
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
			int ix=(int)Math.round(ddx);
    		int iy=(int)Math.round(ddy);
    		int ni=width*iy+ix;
    		if(ni>distMet.length || ni<0)return (double)Byte.MAX_VALUE*0.25;
    		double dist= (distMet[ni]*0.25);    		
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
    	
    	//The commented section looks for intersections to help in heuristics. This doesn't
    	//seem to be necessary anymore, so it and the filter have been commented out.
    	
//    	List<Point2D> allPoints=lines.stream()
//       	     .flatMap(l->l.streamPoints())
//       	     .collect(Collectors.toList());
    	
//    	List<Point2D> dontmerge=GeomUtil.groupPointsCloserThan(allPoints,maxDistanceToConsiderSamePoint)
//    	        .stream()
//    	        .filter(pL->pL.size()>2)
//    	        .flatMap(l->l.stream())
//    	        .collect(Collectors.toList());
    	
    	
    	
    	BiFunction<Double,Double,Double> sample = (x,y)->{
    		//TODO: do interp eventually
    		int ix=(int)Math.round(x);
    		int iy=(int)Math.round(y);
    		int ni=width*iy+ix;
    		if(ni>distMet.length || ni<0)return (double)Byte.MAX_VALUE*0.25;
    		return (double) (distMet[ni]*0.25);    		
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
			        	 double dist=t.k().getAbsoluteClosestDistance();
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
//			         .filter(t->{
//			        	 	Point2D[] closest=t.k().closestPoints();
//							
//							boolean partOfTriple=dontmerge.stream()
//							         .flatMap(p->Stream.of(p.distance(closest[0]),p.distance(closest[1])))
//							         .filter(d->(d<maxDistanceToConsiderSamePoint))
//							         .findAny()
//							         .isPresent();
//							if(partOfTriple)return false;
//							return true;
//			         })
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
        	ChainCode pcode = seq.peek();
            EnumSet<ChainCode> Nb = getNeighbors (bitmap, x, y);

            Point2D pt = null;
            if (Nb.isEmpty ()) {
            } else if (Nb.size () == 1) { //
                pt = seq.add (Nb.iterator ().next ());
            } else {
                // multiple choice; pick best one based on the following
                //  rule: select the yet-unseen one for which
            	//  the number of neighbor neighbors is minimized. 
            	//  If there's a tie, chose the neighbor in the following order:
            	//  E,SE,S,SW,W,NW,N,NE
                ChainCode best = null;
                EnumSet<ChainCode> bestNq = null;
               
                for (ChainCode c : Nb) {
                    int xp = x + c.dx (), yp = y + c.dy ();
                    if (!seq.contains (xp, yp)) {
                        EnumSet<ChainCode> Nq = getNeighbors (bitmap, xp, yp);
                        if (bestNq == null
                            || (Nq.size () < bestNq.size ())
                            // it's actually just choosing directions in the priority of:
                            // E,SE,S,SW,W,NW,N,NE
                            // So it prefers to go to the right and up rather than
                            // to the left and down. Many heuristics were tweaked 
                            // based on this implicit rule, so it won't be changed now
                            // but could probably be better                            
                            
//                            || (Nq.size () == bestNq.size ()
//                                && c.angle() < best.angle())
                        	
                        	
//                    	    || (Nq.size () == bestNq.size ()
//                                && c.cosine(pcode) > best.cosine(pcode))
                    	    
                    	    || (Nq.size () == bestNq.size ()
                            && pcode.priorityChange(c) < pcode.priorityChange(best))
                        	
                        ){
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

	public Bitmap paste(Bitmap bm2, Shape ss) {

		Rectangle2D bounds = ss.getBounds2D();


		int minX = (int) Math.max(bounds.getMinX(), 0);
		int minY = (int) Math.max(bounds.getMinY(), 0);

		int maxX = (int) Math.min(bounds.getMaxX(), this.width());
		int maxY = (int) Math.min(bounds.getMaxY(), this.height());

		for(int x=minX;x<maxX;x++){
			for(int y=minY;y<maxY;y++){
				if(ss.contains(x, y)){
					this.set(x, y, bm2.get(x, y));
				}
			}
		}


		return this;

	}


	private CachedSupplier<Bitmap> _halfer=CachedSupplier.of(()->_half());
	
	
	private Bitmap _half() {

		Bitmap bm2 = new Bitmap(width/2,height/2);
		
		for(int x=0;x<bm2.width;x++){
			for(int y=0;y<bm2.height;y++){
				bm2.set(x, y, this.get(x*2, y*2));
			}
		}
		
		
		return bm2;
	}
	
	public Bitmap half(){
		return _halfer.get();
	}

}
