package gov.nih.ncats.molvec.internal.image;

import java.awt.image.BandedSampleModel;
import java.awt.image.BufferedImage;
import java.awt.image.DataBuffer;
import java.awt.image.Raster;
import java.awt.image.WritableRaster;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.Optional;
import java.util.logging.Logger;

import javax.imageio.ImageIO;

/**
 * Converts an image into a 256 gray scale image.
 *
 *
 */
public class Grayscale {
    static final Logger logger = Logger.getLogger(Grayscale.class.getName());

    
    private Raster grayscale;


    /**
     * Listener interface that accepts pixel values for each image.
     * The order of method calls is:
     * newImage()
     * multipleCalls to accept()
     * finishImage()
     */
    interface GrayscaleListener{
        /**
         * Begin listening to a new image.
         * @param height the height of the new image.
         * @param width the weight of the new image.
         * @param nband the number of bands of the new image.
         */
        void newImage(int height, int width, int nband);

        /**
         * accept the given pixel values in the current image.
         * @param s the current pixel value.
         */
        void accept(int s);

        /**
         * The current image is finished.
         */
        void finishImage();
    }

    private enum NoOpGrayscaleListner implements GrayscaleListener{
        INSTANCE;

        @Override
        public void newImage(int height, int width, int nband) {

        }

        @Override
        public void accept(int s) {

        }

        @Override
        public void finishImage() {

        }

    }

    /**
     * Listener that computes min, max, mean, stddev and a histogram
     * of all pixel values.
     */
    public static class GrayscaleStats implements GrayscaleListener{

        private int min, max;
        private double mean, stddev;
        private int[] histogram;
        @Override
        public void newImage(int height, int width, int nband) {
            max = 0;
            min = Integer.MAX_VALUE;
            histogram = new int[256];
        }

        @Override
        public void accept(int s) {
            histogram[s]++;
        }

        @Override
        public void finishImage() {
            mean = 0;
            stddev = 0;
            int cnt = 0;
            for (int i = 0; i < histogram.length; ++i) {
                int p = histogram[i];
                if (p > 0) {
                    mean += p;
                    ++cnt;
                    max=i;
                    if(i<min)min=i;
                }
            }

            if (cnt > 0) {
                mean /= cnt;
                for (int i = 0; i < histogram.length; ++i) {
                    int p = histogram[i];
                    if (p > 0) {
                        double x = p - mean;
                        stddev += x*x;
                    }
                }
                stddev = Math.sqrt(stddev/cnt);
            }
        }

        public int getMin() {
            return min;
        }

        public int getMax() {
            return max;
        }

        public double getMean() {
            return mean;
        }

        public double getStddev() {
            return stddev;
        }

        public int[] getHistogram() {
            return histogram;
        }
    }
    public Grayscale (Raster raster){
        this(raster, NoOpGrayscaleListner.INSTANCE);
    }
    public Grayscale (Raster raster, GrayscaleListener listener) {
        if (raster == null) {
            throw new IllegalArgumentException ("Input raster is null");
        }
        grayscale = createRaster (raster, listener);
       
    }

    public BufferedImage asNewBufferedImage(){
        BufferedImage image = new BufferedImage(grayscale.getWidth(), grayscale.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
        //this does a copy operation so it's safe
        image.setData(grayscale);
        return image;
    }
    public Raster getRaster () {
        //TODO we never call method with isRescale == true ?
//        return getRaster (false);
        return grayscale;
    }


    /**
     * convert rgb raster into grayscale
     */
    private Raster createRaster (Raster raster, GrayscaleListener listener) {
        int height = raster.getHeight();
        int width = raster.getWidth();
        int nband = raster.getNumBands();

        Grayscaler grayscaler = Grayscaler.getFor(raster);

        listener.newImage(height, width, nband);
//        max = 0;
//        min = Integer.MAX_VALUE;

    
        Optional<AlphaInfo> alphaInfo = grayscaler.computeAlphaInfo(raster);
//        for (int i = 0; i < histogram.length; ++i) {
//            histogram[i] = 0;
//        }

        WritableRaster outRaster = Raster.createWritableRaster
                (new BandedSampleModel
                 (DataBuffer.TYPE_BYTE, width, height, 1), null);
        
        double[] sampleRow = new double[width*raster.getNumBands()];
        double[] resultRow = new double[width];
        for (int y = 0; y < height; ++y) {
        	raster.getPixels(0, y, width, 1, sampleRow);
        	//System.out.println(ff.length + "  vs " + width);
        	for (int x = 0; x < width; x++) {
        		
        		double[] pp=Arrays.copyOfRange(sampleRow,x*nband,x*nband+nband);
        		grayscaler.adjustAlpha(pp, alphaInfo);

        		int s = grayscaler.grayscaleValue (pp) & 0xff;
                resultRow[x]=s;
                listener.accept(s);

            }
        	outRaster.setSamples(0, y, width, 1, 0, resultRow);
        }
        raster = outRaster;

        listener.finishImage();
        
        

        return raster;
    }

    public BufferedImage getImage () {
        if (grayscale == null) {
            throw new IllegalStateException ("No buffer available");
        }

        BufferedImage img = new BufferedImage
            (grayscale.getWidth(), grayscale.getHeight(), 
             BufferedImage.TYPE_BYTE_GRAY);
        img.setData(grayscale);
        return img;
    }

    public void write (OutputStream out) throws IOException {
        ImageIO.write(getImage (), "png", out);
    }

    public static int grayscale (double[] rgb) {
    	if(rgb.length==4){
    		return (int) ((0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] + .5) * rgb[3]/255D);
    	}else if(rgb.length==1){
            return (int)rgb[0];
    	}else if(rgb.length==2){

            return (int)(rgb[0] * (rgb[1]/255D));
        }else{
    		return (int) (0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] + .5);
    	}
    }
    
    public static double[] hsv (double[] rgb) {
    	
    	if(rgb.length>=3){
    		double rp = rgb[0]/255;
        	double gp = rgb[1]/255;
        	double bp = rgb[2]/255;
        	double cmax= Math.max(Math.max(rp, gp),bp);
        	double cmin= Math.max(Math.min(rp, gp),bp);
        	double delta=cmax-cmin;
        	
        	double s = 0;
        	if(cmax>0){
        		s=delta/cmax;
        	}
        	double v = cmax;
        	
        	double h=0;
        	if(rp>gp && rp>bp){
        		h=60*(gp-bp)/delta;
        		if(h<0){
        			h=360+h;
        		}
        	}else if(gp>rp && gp>bp){
        		h=60*((bp-rp)/delta+2);
        	}else if(bp>gp && bp>rp){
        		h=60*((rp-gp)/delta+4);
        	}
        		
        	return new double[]{h,s,v};
    	}
    	return new double[]{0,0,rgb[0]};
    }
    
    public static int sOnly255 (double[] rgb) {
    	
    	if(rgb.length>=3){
    		double rp = rgb[0];
        	double gp = rgb[1];
        	double bp = rgb[2];
        	double cmax= Math.max(Math.max(rp, gp),bp);
        	return (int)cmax;
    	}
    	return (int)rgb[0];
    }

    private enum Grayscaler{
    	//USES hard-coded ratios for going to grey
        RGBA_G(3){
            @Override
            protected int grayscaleValue(double[] rgb) {
            	int start=(int) ((0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] + 0.5) );
            	return (int)(start * rgb[3]/255D);
            }

        },
        RGB_G{
            @Override
            protected int grayscaleValue(double[] rgb) {
            	int start=(int) (0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]+0.5);
            	return start;
            }
        },
        
        //Uses HSV model with "V" component
        RGBA_V(3){
            @Override
            protected int grayscaleValue(double[] rgb) {
//            	double[] hsv1=hsv(rgb);
//            	int vstart=(int) ((hsv1[2])*255);
//            	return (int)(vstart* rgb[3]/255D);
            	return (int)(sOnly255(rgb)* rgb[3]/255D);
            }

        },
        RGB_V{
            @Override
            protected int grayscaleValue(double[] rgb) {
//            	double[] hsv1=hsv(rgb);
//            	int vstart=(int) ((hsv1[2])*255);
//            	return vstart;
            	return sOnly255(rgb);
            }
        },
        GRAY{
            @Override
            protected int grayscaleValue(double[] rgb) {
                return (int) rgb[0];
            }

        },
        GRAY_ALPHA(1){
            @Override
            protected int grayscaleValue(double[] rgb) {
                return (int)(rgb[0] * (rgb[1]/255D));
            }

        }        
        ;

        private final int alphaBand;
        Grayscaler(){
            this(-1);
        }
        Grayscaler(int alphaBand){
            this.alphaBand = alphaBand;
        }

        protected abstract int grayscaleValue(double[] rgb);

        protected int getAlphaBand(){
            return alphaBand;
        }

        protected Optional<AlphaInfo> computeAlphaInfo(Raster raster){
            if(alphaBand == -1){
                return Optional.empty();
            }
            int width = raster.getWidth();
            int height = raster.getHeight();
            int[] row = new int[width];
            AlphaInfo info = new AlphaInfo();
            for (int j = 0; j < height; ++j) {
                raster.getSamples(0, j, width, 1, getAlphaBand(), row);
                for (int i = 0; i < width; ++i) {
                    info.accept(row[i]);
                }
            }
            return Optional.of(info);
        }

        public static Grayscaler getFor(Raster raster){
            int nbands = raster.getNumBands();
            switch(nbands){
                case 1: return Grayscaler.GRAY;
                case 2 :  return Grayscaler.GRAY_ALPHA;
                case 3 : return Grayscaler.RGB_V;
                default : return Grayscaler.RGBA_V;
            }
        }


        public void adjustAlpha(double[] pp, Optional<AlphaInfo> alphaInfo){
            if(!alphaInfo.isPresent()){
                return;
            }
            AlphaInfo info = alphaInfo.get();

            if(info.maxAlpha<=info.minAlpha) {
                pp[alphaBand]=255;
            }else {
                pp[alphaBand]=(pp[alphaBand]-info.minAlpha)/(info.maxAlpha-info.minAlpha);
                pp[alphaBand]*=255;
                //assume that there's supposed to be more
                //transparent things than opaque things.
                //If that's not the case, invert the alpha
                //channel
                if(info.moreAlphaCount>info.lessAlphaCount) {
                    pp[alphaBand]=255-pp[alphaBand];
                }
            }
        }
    }
    private static class AlphaInfo{
        private int maxAlpha=0;
        private int minAlpha=255;
        private int moreAlphaCount=0;
        private int lessAlphaCount=0;

        public void accept(int pixel){
            if (pixel > maxAlpha) {
                maxAlpha = pixel;
            }
            if (pixel < minAlpha) {
                minAlpha = pixel;
            }
            if(pixel > 128) {
                moreAlphaCount++;
            }
            if(pixel <= 128) {
                lessAlphaCount++;
            }
        }
    }
}
