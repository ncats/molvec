package gov.nih.ncats.molvec.image;

import java.awt.image.BandedSampleModel;
import java.awt.image.BufferedImage;
import java.awt.image.DataBuffer;
import java.awt.image.Raster;
import java.awt.image.RescaleOp;
import java.awt.image.WritableRaster;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.logging.Logger;

import javax.imageio.ImageIO;

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
        int nband = raster.getNumBands();
        

        max = 0;
        min = Integer.MAX_VALUE;
        
        int maxAlpha=0;
        int minAlpha=255;
        int moreAlphaCount=0;
        int lessAlphaCount=0;
        
        //Assumes RGBA
        if(nband>=4) {
        	int[] row = new int[width];
            for (int j = 0; j < height; ++j) {
            	raster.getSamples(0, j, width, 1, 3, row);
            	for (int i = 0; i < width; ++i) {
        		    int pixel = row[i];
                    if (pixel > maxAlpha) maxAlpha = pixel;
                    if (pixel < minAlpha) minAlpha = pixel;
                    if(pixel > 128) moreAlphaCount++;
                    if(pixel <= 128) lessAlphaCount++;
                    
        	    }
        	}
        }
        
    
        
        for (int i = 0; i < histogram.length; ++i)
            histogram[i] = 0;

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
        		
        		if(pp.length==4 ) {
        			if(maxAlpha<=minAlpha) {
        				pp[3]=255;
        			}else {
        				pp[3]=(pp[3]-minAlpha)/(maxAlpha-minAlpha);
        				pp[3]*=255;
        				//assume that there's supposed to be more
        				//transparent things than opaque things.
        				//If that's not the case, invert the alpha
        				//channel
        				if(moreAlphaCount>lessAlphaCount) {
        					pp[3]=255-pp[3];
        				}
        			}
        		}
        		
        		int s = grayscale (pp) & 0xff;
                resultRow[x]=s;
                ++histogram[s];
            }
        	outRaster.setSamples(0, y, width, 1, 0, resultRow);
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
    		return (int) ((0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] + .5) * ((1.0/255.0)*rgb[3]));
    	}else if(rgb.length==1){
    		return (int)rgb[0];
    	}else{
    		return (int) (0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] + .5);
    	}
    }
}
