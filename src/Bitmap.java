import java.io.*;

import java.awt.geom.*;
import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;

import java.awt.Shape;
import java.awt.Rectangle;
import java.awt.Point;
import java.awt.Polygon;
import java.awt.image.*;

import javax.imageio.*;

import com.sun.media.jai.codec.*;

/**
 * A bitmap image
 */
public class Bitmap implements Serializable {
    private static final long serialVersionUID = 0x5f1f54d8fed49ab3l;

    private static final Logger logger =
        Logger.getLogger (Bitmap.class.getName ());

    private static final double EPS = 0.000001;

    private static final boolean DEBUG;

    static {
        boolean debug = false;
        try {
            debug = Boolean.getBoolean ("bitmap.debug");
        } catch (Exception ex) {
        }
        DEBUG = debug;
    }

    static final double DEFAULT_AEV_THRESHOLD = 1.5;
    static final double DEFAULT_SIGMA_THRESHOLD = 1.2;
    static final int DEFAULT_ADAPTIVE_BOX_RADIUS = 20;
    static final int DEFAULT_ADAPTIVE_MIN_THRESHOLD = 20;
    static final int DEFAULT_ADAPTIVE_MIN_STDDEV = 20;


    /**
     * bounding box shape
     */
    public static enum Bbox {
        Rectangular,
            Polygon
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


    /*
     * commonly used TIFF tags
     */
    static final int TAG_IMAGELENGTH = 257;
    static final int TAG_IMAGEWIDTH = 256;
    static final int TAG_RESOLUTIONUNIT = 296;
    static final int RESOLUTIONUNIT_NONE = 1; // no absolute unit
    static final int RESOLUTIONUNIT_INCH = 2; // inch
    static final int RESOLUTIONUNIT_CENT = 3; // centimeter
    static final int TAG_XRESOLUTION = 282;
    static final int TAG_YRESOLUTION = 283;
    static final int TAG_BITSPERSAMPLE = 258;
    static final int TAG_SAMPLESPERPIXEL = 277;
    static final int TAG_PHOTOMETRIC = 262;
    static final int PHOTOMETRIC_WHITEISZERO = 0;
    static final int PHOTOMETRIC_BLACKISZERO = 1;
    static final int PHOTOMETRIC_RGB = 2;
    static final int PHOTOMETRIC_PALETTE = 3;
    static final int PHOTOMETRIC_MASK = 4;
    static final int TAG_STRIPOFFSETS = 273;
    static final int TAG_ROWSPERSTRIP = 278;
    static final int TAG_STRIPBYTECOUNTS = 279;
    static final int TAG_COMPRESSION = 259;
    static final int TAG_IMAGEDESCRIPTION = 270;
    static final int TAG_MAXSAMPLEVALUE = 281;
    static final int TAG_MINSAMPLEVALUE = 280;
    static final int TAG_PLANARCONFIGURATION = 284;
    static final int PLANARCONFIGURATION_CHUNKY = 1;
    static final int PLANARCONFIGURATION_PLANAR = 2;

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

    public static Bitmap createBitmap (Raster raster) {
        SampleModel model = raster.getSampleModel ();
        int band = model.getNumBands ();
        if (band > 1) {
            throw new IllegalArgumentException
                ("Can handle sample with multiple channels");
        }

        return adaptiveThreshold
            (raster, Math.max((DEFAULT_ADAPTIVE_BOX_RADIUS*raster.getWidth())/350,DEFAULT_ADAPTIVE_BOX_RADIUS),
             DEFAULT_SIGMA_THRESHOLD,
             DEFAULT_ADAPTIVE_MIN_THRESHOLD,
             DEFAULT_ADAPTIVE_MIN_STDDEV);
    }

    /*
     * Simple threshold based on full image mean and standard deviation. 
     */
    public static Bitmap sigmaThreshold (Raster inRaster, double sigma) {
        Bitmap bm = new Bitmap (inRaster.getWidth (), inRaster.getHeight ());
        double max = -1, min = Double.MAX_VALUE;
        double sum = 0;
        double sumSquare = 0;
        for (int y = 0; y < bm.height; ++y) {
            for (int x = 0; x < bm.width; ++x) {
                double pel = inRaster.getSampleDouble (x, y, 0);
                sum += pel;
                sumSquare += pel * pel;
                if (pel < min)
                    min = pel;
                if (pel > max)
                    max = pel;
            }
        }
        long tot = bm.height * bm.width;
        double mean = sum / tot;
        double stdDEV = Math.sqrt (sumSquare / tot - mean * mean);
        double threshold = mean + stdDEV * sigma;
        for (int y = 0; y < bm.height; ++y) {
            for (int x = 0; x < bm.width; ++x) {
                double pel = inRaster.getSampleDouble (x, y, 0);
                bm.set (x, y, pel >= threshold);
            }
        }
        return bm;
    }

    /*
     * Adaptive Threshold based on local pixel average Loosely based on Bradley
     * & Roth Integral Image Adaptive Thresholding Added ability to use a
     * specific SIGMA minimum threshold. General Algorithm: 
     * 	1)Precompute the integral image (used to get local mean) 
     * 	2)Precompute the square integral image (used to get local 
     * 		standard dev) 
     * 	3)For each pixel: 
     * 		a)get the sum and sum of squares for those pixels in 
     * 		  an (nsize*2+1) sized box around the given pixel. 
     * 		  (using the precomputed integrals) 
     * 		b)calculate mean =  sum/ pixel count 
     * 		c)calculate stDEV = sqrt(sumSquares/(pixel count)-* mean^2) 
     * 		d)set threshold for given pixel at mean+stDEV*SIGMA 
     * 		  (where sigma is a provided constant)
     * UPDATE: Integral images calculated only for rows needed, on the fly
     */
    public static Bitmap adaptiveThreshold (Raster inRaster, int nsize,
                                            double sigma, int absMin,
                                            double minSigma) {
        Bitmap bm = new Bitmap (inRaster.getWidth (), inRaster.getHeight ());
        int[] yMap = new int[bm.height];
        double[][] intLines = new double[nsize * 2 + 2][bm.width];
        double[][] intSquareLines = new double[nsize * 2 + 2][bm.width];
        

        // make Integral-Image for quick averaging
        for (int y = 0; y < nsize * 2 + 2; ++y) {
            yMap[y] = y;
            addIntLines (inRaster, yMap, intLines, intSquareLines, bm.width, y);
        }

        for (int y = 0; y < bm.height; ++y) {
            int y1 = Math.max (y - nsize, 0);
            int y2 = Math.min (y + nsize, bm.height - 1);
            int boxHeight = (y2 - y1 + 1);
            double[] topIntLine = null;
            double[] topSquareIntLine = null;
            double[] bottomIntLine = null;
            double[] bottomSquareIntLine = null;

            if (y1 >= 2) {
                if (boxHeight >= nsize * 2 + 1) {
                    yMap[y2] = yMap[y1 - 2];
                    addIntLines (inRaster, yMap, intLines, intSquareLines,
                                 bm.width, y2);
                }
            }
            if (y1 > 0) {
                topIntLine = intLines[yMap[y1 - 1]];
                topSquareIntLine = intSquareLines[yMap[y1 - 1]];
            }
            bottomIntLine = intLines[yMap[y2]];
            bottomSquareIntLine = intSquareLines[yMap[y2]];

            for (int x = 0; x < bm.width; ++x) {
                // box to average over:
                int x1 = Math.max (x - nsize, 0);
                int x2 = Math.min (x + nsize, bm.width - 1);
                int count = (x2 - x1 + 1) * boxHeight;
                double sum = bottomIntLine[x2];
                double sumSquare = bottomSquareIntLine[x2];
                if (y1 > 0) {
                    sum += -topIntLine[x2];
                    sumSquare += -topSquareIntLine[x2];
                }
                if (x1 > 0) {
                    sum += -bottomIntLine[x1 - 1];
                    sumSquare += -bottomSquareIntLine[x1 - 1];
                }
                if (y1 > 0 && x1 > 0) {
                    sum += topIntLine[x1 - 1];
                    sumSquare += topSquareIntLine[x1 - 1];
                }

                double mean = sum / ((double) count);
                double stdDEV = Math.sqrt (Math.abs (sumSquare / ((double) count)
                                                     - mean * mean));
                double threshold = Math.max (mean + stdDEV * sigma, absMin);
                double pel = inRaster.getSampleDouble (x, y, 0);
                bm.set (x, y, pel > threshold && stdDEV > minSigma);
            }
        }
        return bm;
    }

    private static void addIntLines (Raster inRaster, int[] yMap,
                                     double[][] intLines,
                                     double[][] intSquareLines,
                                     int width, int y) {
        double sum = 0;
        double sumSquare = 0;
        double[] intLine = intLines[yMap[y]];
        double[] intSquareLine = intSquareLines[yMap[y]];

        double[] intPrevLine = null;
        double[] intSquarePrevLine = null;
        if (y > 0) {
            intPrevLine = intLines[yMap[y - 1]];
            intSquarePrevLine = intSquareLines[yMap[y - 1]];
        }

        for (int x = 0; x < width; ++x) {
            sum += inRaster.getSampleDouble (x, y, 0);
            sumSquare += inRaster.getSampleDouble (x, y, 0)
                * inRaster.getSampleDouble (x, y, 0);
            if (y == 0) {
                intLine[x] = sum;
                intSquareLine[x] = sumSquare;
            } else {
                intLine[x] = sum + intPrevLine[x];
                intSquareLine[x] = sumSquare + intSquarePrevLine[x];
            }
        }
    }

    public static Bitmap readtif (String file) throws IOException {
        return readtif (new FileInputStream (file));
    }

    public static Bitmap readtif (File file) throws IOException {
        return readtif (new FileInputStream (file));
    }

    public static Bitmap readtif (InputStream is) throws IOException {
        ImageDecoder decoder = ImageCodec.createImageDecoder
            ("TIFF", is, new TIFFDecodeParam ());

        TIFFDirectory tif = new TIFFDirectory
            (decoder.getInputStream (), 0);
        TIFFField[] fields = tif.getFields ();

        String unit = "";
        double xres = 0., yres = 0.;
        int rows = -1, photometric = -1, bpp = -1;
        for (int j = 0; j < fields.length; ++j) {
            TIFFField f = fields[j];
            int tag = f.getTag ();
            switch (tag) {
            case TAG_RESOLUTIONUNIT: {
                int u = f.getAsInt (0);
                if (u == RESOLUTIONUNIT_NONE) {
                } else if (u == RESOLUTIONUNIT_INCH) {
                    unit = "in";
                } else if (u == RESOLUTIONUNIT_CENT) {
                    unit = "cm";
                }
            }
                break;

            case TAG_XRESOLUTION:
                xres = f.getAsFloat (0);
                break;

            case TAG_YRESOLUTION:
                yres = f.getAsFloat (0);
                break;

            case TAG_ROWSPERSTRIP:
                //rows = f.getAsInt(0);
                break;

            case TAG_PHOTOMETRIC:
                photometric = f.getAsInt (0);
                break;

            case TAG_BITSPERSAMPLE:
                bpp = f.getAsInt (0);
                break;
                /*
                  case TAG_IMAGEWIDTH:
                  width = f.getAsFloat(0);
                  break;

                  case TAG_IMAGELENGTH:
                  height = f.getAsFloat(0);
                  break;
                */
            }
        }

        if (bpp != 1) {
            throw new IllegalArgumentException ("BitsPerSample != 1");
        }

        RenderedImage decodedImage = decoder.decodeAsRenderedImage ();
        Raster raster = decodedImage.getData ();
        if (bpp == 1) {
            Bitmap bm = new Bitmap (raster.getWidth (), raster.getHeight ());
            for (int y = 0; y < bm.height; ++y) {
                int band = bm.scanline * y;
                for (int x = 0; x < bm.width; ++x) {
                    int pel = raster.getSample (x, y, 0);
                    if (pel == 1) {
                        bm.data[band + x / 8] |= MASK[x % 8];
                    }
                }
            }

            if (photometric == PHOTOMETRIC_BLACKISZERO) {
                // flip
                for (int i = 0; i < bm.data.length; ++i) {
                    bm.data[i] = (byte) (~bm.data[i] & 0xff);
                }
            }

            return bm;
        }

        logger.info ("TIFF image: bpp=" + bpp
                     + " channels=" + raster.getNumBands ());
        return adaptiveThreshold (raster.getNumBands () == 1
                                  ? raster : createGrayscaleRaster (raster),
                                  DEFAULT_ADAPTIVE_BOX_RADIUS,
                                  DEFAULT_SIGMA_THRESHOLD,
                                  DEFAULT_ADAPTIVE_MIN_THRESHOLD,
                                  DEFAULT_ADAPTIVE_MIN_STDDEV);
    }

    /*
     * convert rgb raster into grayscale
     */
    static Raster createGrayscaleRaster (Raster inRaster) {
        int height = inRaster.getHeight ();
        int width = inRaster.getWidth ();
        WritableRaster outRaster = Raster.createWritableRaster
            (new BandedSampleModel
             (DataBuffer.TYPE_BYTE, width, height, 1), null);

        double[] sample = new double[3];
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int s = rgb2grayscale (inRaster.getPixel (x, y, sample));
                outRaster.setSample (x, y, 0, s & 0xff);
            }
        }


        try {
            BufferedImage img = new BufferedImage
                (width, height, BufferedImage.TYPE_BYTE_GRAY);
            img.setData (outRaster);
            ImageIO.write (img, "png", new File ("gray.png"));
        } catch (Exception ex) {
            ex.printStackTrace ();
        }

        return outRaster;
    }

    static int rgb2grayscale (double[] rgb) {
        return (int) (0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] + .5);
    }

    public void writetif (String file) throws IOException {
        writetif (new FileOutputStream (file));
    }

    public void writetif (File file) throws IOException {
        writetif (new FileOutputStream (file));
    }

    public void writetif (OutputStream os) throws IOException {
        TIFFEncodeParam param = new TIFFEncodeParam ();
        param.setCompression (TIFFEncodeParam.COMPRESSION_GROUP4);
        ImageEncoder encoder = ImageCodec.createImageEncoder
            ("TIFF", os, param);
        encoder.encode (createBufferedImage ());
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

    public Object clone () {
        return new Bitmap (this);
    }

    public int width () {
        return width;
    }

    public int height () {
        return height;
    }

    public boolean get (int x, int y) {
        if (x >= 0 && x < width && y >= 0 && y < height) {
            return ((data[y * scanline + x / 8] & 0xff) & MASK[x % 8]) != 0;
        }
        return false;
    }

    // same as get() but without the bound checking
    public boolean isOn (int x, int y) {
        return ((data[y * scanline + x / 8] & 0xff) & MASK[x % 8]) != 0;
    }

    public void set (int x, int y, boolean on) {
        int loc = y * scanline + x / 8;
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
            int band = y * scanline;
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

    public Bitmap crop (Shape s) {
        Rectangle r = s.getBounds ();
        if (r.width == 0 || r.height == 0) {
            return null;
        }

        Bitmap dst = new Bitmap (r.width + 1, r.height + 1);
        int x1 = Math.min (width, r.x + r.width);
        int y1 = Math.min (height, r.y + r.height);
        int x0 = r.x, y0 = r.y;

        int i, j = 0;
        for (int y = y0; y <= y1; ++y, ++j) {
            i = 0;
            for (int x = x0; x <= x1; ++x, ++i) {
                if (x == x1 || y == y1 || s.contains (x, y))
                    dst.set (i, j, get (x, y));
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
                            copy[scanline * y + x / 8] &= ~MASK[x % 8];
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
     * This version is similar to the NWG algorithm above. It's based on
     * R. Carrsco, M. Forcada, A note on the Nagendraprasad-Wang-Gupta
     * thinning algorithm, Pattern Recognition Letters, 16, 539-541, 1995.
     */
    public Bitmap thin () {
        Bitmap thin = new Bitmap (this);
        byte[] copy = new byte[this.data.length];
        System.arraycopy (thin.data, 0, copy, 0, copy.length);

        int[] p = new int[10]; // window
        int parity = 1;
        boolean changed;
        
        do {
            changed = false;
            parity = 1 - parity;

            for (int y = 0; y < height; ++y)
                for (int x = 0; x < width; ++x) {
                    if (!thin.isOn(x, y))
                        continue;
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
                        int ep = (thin.p2(x, y) + thin.p4(x, y))*thin.p0(x, y)
                            * thin.p6(x, y);
                        int fp = (thin.p6(x, y) + thin.p0(x, y))*thin.p4(x, y)
                            * thin.p2(x, y);
                        if ((parity == 0 && ep == 0) 
                            || (parity == 1 && fp == 0)) {
                            // delete this pixel
                            copy[scanline * y + x / 8] &= ~MASK[x % 8];
                            changed = true;
                        }
                    }
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

    public List<Shape> connectedComponents (Bbox shape) {
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
    }


    List<Shape> connectedComponentPolygonShapes
        (short[] eqvtab, short[][] labels) {

        Map<Short, List<Point>> coords = new HashMap<Short, List<Point>> ();
        for (int y = 0; y < height; ++y)
            for (int x = 0; x < width; ++x) {
                short label = labels[y][x];
                if (label != 0) {
                    short l = label;
                    /* find equivalence class */
                    while (eqvtab[l] > 0)
                        l = eqvtab[l];

                    labels[y][x] = l;
                    Point pt = new Point (x, y);

                    List<Point> pts = coords.get (l);
                    if (pts == null) {
                        coords.put (l, pts = new ArrayList<Point> ());
                    }
                    pts.add (pt);
                }
            }

        List<Shape> comps = new ArrayList<Shape> ();
        for (List<Point> pts : coords.values ()) {
            Polygon hull = GeomUtil.convexHull (pts.toArray (new Point[0]));
            comps.add (hull);
        }

        return comps;
    }

    List<Shape> connectedComponentRectangularShapes
        (short[] eqvtab, short[][] labels) {

        Map<Short, Rectangle> ltab = new HashMap<Short, Rectangle> ();
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
        return segments (2, DEFAULT_AEV_THRESHOLD);
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

    public static void main (String[] argv) throws Exception {
        Bitmap bm = new Bitmap (16, 16);
        java.util.Random rand = new java.util.Random ();
        int n = 0;
        for (int y = 0; y < bm.height (); ++y) {
            for (int x = 0; x < bm.width (); ++x) {
                boolean on = rand.nextDouble () > .5;
                if (on) {
                    //System.out.println("pixel on at "+x+" "+y);
                    ++n;
                }
                bm.set (x, y, on);
            }
        }

        /*
          for (int x = 4; x < 6; ++x) {
          //bm.set(x, 4, true);
          bm.set(x, 5, true);
          }

          for (int y = 6; y < 10; ++y) {
          for (int x = 2; x < 8; ++x) {
          boolean on = rand.nextDouble() > .5;
          if (on) {
          //System.out.println("pixel on at "+x+" "+y);
          ++n;
          }
          bm.set(x, y, true);
          }
          }

          bm.dump(System.out);
          System.out.println(n + " pixels on!");
          ChainCodeSequence ccs = bm.chainCode();
          System.out.println("chain code "+ccs.length());

          Bitmap thin = bm.skeleton();
          thin.dump(System.out);

          List<Rectangle> bboxes = bm.connectedComponents();
          System.out.println(bboxes.size()+ " connected components!");
          for (Rectangle r : bboxes) {
          System.out.println(r);
          Bitmap c = bm.crop(r);
          c.dump(System.out);
          }
          bm.write("png", new File ("bitmap.png"));

          List<Point> corners = bm.corners(3);
          System.out.println(corners.size()+ " corners!");
          for (Point pt : corners) {
          System.out.println(pt);
          }
        */

        if (argv.length == 0) {
            System.err.println ("Usage: Bitmap FILE.tif");
            System.exit (1);
        }

        Bitmap tif = readtif (new File (argv[0]));
        Bitmap ske = tif.skeleton ();
        //tif.write("png", new File ("bitmap.png"));
        tif.dump (System.out);
        /*
          List<Shape> cc = tif.connectedComponents();
          for (Shape s : cc) {
          System.out.println("crop: "+s.getBounds());
          Bitmap c = tif.crop(s);
          c.dump(System.out);
          }
        */
        System.out.println ("Image " + tif.width () + "x" + tif.height ());

        char[][] ascii = new char[tif.height ()][tif.width ()];
        for (int y = 0; y < tif.height (); ++y) {
            for (int x = 0; x < tif.width (); ++x) {
                ascii[y][x] = tif.get (x, y) ? '*' : '.';
            }
        }

        List<ChainCodeSequence> seqs = tif.chainCodes ();
        System.out.println ("Num chain codes: " + seqs.size ());
        for (ChainCodeSequence s : seqs) {
            int x = s.getStartX (), y = s.getStartY ();
            for (ChainCode c : s.getSequence ()) {
                ascii[y][x] = c.ch ();
                x += c.dx ();
                y += c.dy ();
            }
            ascii[y][x] = '#'; // mark the end of this chain code
        }

        for (int y = 0; y < tif.height (); ++y) {
            for (int x = 0; x < tif.width (); ++x)
                System.out.print (ascii[y][x]);
            System.out.println ();
        }

        List<Point2D> dps = tif.dominantPoints ();
        System.out.println ("** " + dps.size () + " dominant points");
        for (int i = 0; i < dps.size (); ++i) {
            Point2D pt = dps.get (i);
            System.out.println (" ++ " + (i + 1) + ": " + pt);
            ascii[(int) pt.getY ()][(int) pt.getX ()] = '@';
        }


        for (int y = 0; y < tif.height (); ++y) {
            for (int x = 0; x < tif.width (); ++x)
                System.out.print (ascii[y][x]);
            System.out.println ();
        }
        //tif.segments(10, 2);
    }

    public static class Skeleton {
        public static void main (String[] argv) throws Exception {
            if (argv.length < 2) {
                System.err.println ("Usage: Bitmap$Skeleton INTIF OUTTIF");
                System.exit (1);
            }

            Bitmap bm = Bitmap.readtif (new File (argv[0]));
            bm.skeleton ().writetif (new File (argv[1]));
        }
    }
}
