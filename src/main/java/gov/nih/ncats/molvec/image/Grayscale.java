package gov.nih.ncats.molvec.image;

import java.io.IOException;
import java.io.OutputStream;
import java.util.logging.Logger;

import java.awt.image.BufferedImage;
import java.awt.image.DataBuffer;
import java.awt.image.Raster;
import java.awt.image.WritableRaster;
import java.awt.image.BandedSampleModel;
import java.awt.image.RescaleOp;

import javax.imageio.ImageIO;
import gov.nih.ncats.molvec.algo.Peaks;

/**
 * 256 gray scale 
 */
public class Grayscale {
    static final Logger logger = Logger.getLogger(Grayscale.class.getName());

    private Raster grayscale;
    private int[] histogram = new int[256];
    private double mean, stddev;
    private int max, min;

    public Grayscale () {
    }

    public Grayscale (Raster raster) {
        setRaster (raster);
    }

    public void setRaster (Raster raster) {
        if (raster == null) {
            throw new IllegalArgumentException ("Input raster is null");
        }
        grayscale = createRaster (raster);
    }

    public Raster getRaster () { 
        return getRaster (false);
    }

    public Raster getRaster (boolean isRescale) { 
        if (isRescale) {
            // rescale to 8-bit
            double scale = Math.max(256./(max-min+1),1);
            RescaleOp rescale = new RescaleOp 
                ((float)scale, -(float)scale*min, null);
            Raster raster = rescale.filter(grayscale, null);
            rescale = new RescaleOp (-1, 255, null);
            return rescale.filter(raster, null);
        }
        return grayscale; 
    }

    /*
     * convert rgb raster into grayscale
     */
    protected Raster createRaster (Raster raster) {
        int height = raster.getHeight();
        int width = raster.getWidth();

        max = 0;
        min = Integer.MAX_VALUE;
    	for (int i = 0; i < width; ++i) {
    	    for (int j = 0; j < height; ++j) {
                int pixel = raster.getSample(i, j, 0);
                if (pixel > max) max = pixel;
                if (pixel < min) min = pixel;
    	    }
    	}
        int range = max - min;
        logger.info("## range="+range+", min="+min+", max="+max);

        for (int i = 0; i < histogram.length; ++i)
            histogram[i] = 0;

        WritableRaster outRaster = Raster.createWritableRaster
                (new BandedSampleModel
                 (DataBuffer.TYPE_BYTE, width, height, 1), null);
        double[] sample = new double[Math.max(raster.getNumBands(),3)];
        
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
            	
                int s = grayscale (raster.getPixel(x, y, sample)) & 0xff;
                //System.out.println(Arrays.toString(sample));
                outRaster.setSample (x, y, 0, s);
                ++histogram[s];
            }
        }
        raster = outRaster;

        mean = 0;
        stddev = 0;
        int cnt = 0;
        for (int i = 0; i < histogram.length; ++i) {
            int p = histogram[i];
            if (p > 0) {
                mean += p;
                ++cnt;
            }
        }

        /*
          try {
          BufferedImage img = new BufferedImage
          (raster.getWidth(), raster.getHeight(), 
          BufferedImage.TYPE_BYTE_GRAY);
          img.setData(raster);
          ImageIO.write(img, "png", new java.io.FileOutputStream ("gray.png"));
          }
        catch (IOException ex) {
            ex.printStackTrace();
        }
        */

        { Peaks peaks = new Peaks ();
            int[] p = peaks.detect(histogram);
            for (int i = 0; i < p.length; ++i) {
                System.out.print(" "+p[i]);
            }
            System.out.println();
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

        logger.info("mean: "+mean+", std: "+stddev);

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

    public int[] histogram () { return histogram; }
    public double mean () { return mean; }
    public double stddev () { return stddev; }

    public void write (OutputStream out) throws IOException {
        ImageIO.write(getImage (), "png", out);
    }

    public static int grayscale (double[] rgb) {
    	if(rgb.length==4){
    		
    		return (int) ((0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] + .5) * (1-(1.0/255.0)*rgb[3]));
    	}else{
    		return (int) (0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] + .5);
    	}
    }
}
