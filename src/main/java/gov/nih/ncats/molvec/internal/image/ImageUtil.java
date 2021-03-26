package gov.nih.ncats.molvec.internal.image;

import com.twelvemonkeys.imageio.stream.ByteArrayImageInputStream;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.IndexColorModel;
import java.awt.image.Raster;
import java.awt.image.RenderedImage;
import java.awt.image.RescaleOp;
import java.io.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.logging.Logger;

import javax.imageio.ImageIO;
import javax.imageio.ImageReader;
import javax.imageio.stream.ImageInputStream;


public class ImageUtil {
    private static final Logger logger = Logger.getLogger
	(ImageUtil.class.getName());


    private static BufferedImage toRGBColorModel(BufferedImage bi2){
		
		int nwidth=bi2.getWidth();
		int nheight=bi2.getHeight();
		
        // creates output image
        BufferedImage outputImage = new BufferedImage(nwidth, nheight,ColorModel.BITMASK);
 
        // scales the input image to the output image
        Graphics2D g2d = outputImage.createGraphics();
        
        g2d.drawImage(bi2, 0, 0, nwidth, nheight, null);
        g2d.dispose();
        Map<Integer, Integer> rgbConverterMap = new HashMap<>();
        for (int x = 0; x < outputImage.getWidth(); x++) {
            for (int y = 0; y < outputImage.getHeight(); y++) {
                int rgba = outputImage.getRGB(x, y);
                int newValue = rgbConverterMap.computeIfAbsent(rgba, k->{
                    Color col = new Color(rgba, true);
                    col = new Color(255 - col.getRed(),
                            255 - col.getGreen(),
                            255 - col.getBlue());
                    return col.getRGB();
                });

                outputImage.setRGB(x, y, newValue);
            }
        }
        

		
		return outputImage;
	}

    public static BufferedImage decode (BufferedImage bi) {
    	
    	if(bi.getColorModel() instanceof IndexColorModel){
    		bi=toRGBColorModel(bi);
    	}
        Raster raster = bi.getData();

        int bands=raster.getNumBands();
    	int max = 0;
    	int min = Integer.MAX_VALUE;
    	int[] pix = new int[raster.getWidth()];
    	for (int j = 0; j < raster.getHeight(); ++j) {
//    		for (int i = 0; i < raster.getWidth(); ++i) {
    		for(int k=0;k<bands;k++){
    			raster.getSamples(0, j, raster.getWidth(), 1, k, pix);
    	    	for(int i=0;i<pix.length;i++){
	                int pixel = pix[i];
	                //System.out.print(pixel%10);
	                if (pixel > max) max = pixel;
	                if (pixel < min) min = pixel;
    	    	}
    	    	
    	    }
    	    //System.out.println();
    	}

//        logger.info("## dynamic range: "+(max-min)+" color model: "+bi.getColorModel());
        
        // rescale to 8-bit
        double scale = Math.max(256./(max-min+1),1);
        RescaleOp rescale = new RescaleOp 
            ((float)scale, -(float)scale*min, null);
        Raster scaled = rescale.filter(raster, null);
        rescale = new RescaleOp (-1, 255, null);
        raster = rescale.filter(scaled, null);

        Grayscale grayscale = new Grayscale (raster);
    	return grayscale.getImage();
    }
    private static boolean isTiff(byte[] f) throws IOException{

        //0x4949 or 0x4d4d
        return ((f[0] == 0x49 && f[1] == 0x49) || (f[0] == 0x4d && f[1] == 0x4d));
    }
    private static boolean isPng(byte[] f) throws IOException{

        //0x89PNG
        return f[0] == 0x89 && f[1] == 0x50 && f[1] == 0x4E && f[1] == 0x47;
    }
    public static BufferedImage grayscale (BufferedImage bi) {
        Grayscale grayscale = new Grayscale (bi.getData());
    	return grayscale.asNewBufferedImage();
    }


    public static BufferedImage grayscale (byte[] file) throws IOException {
        // Create input stream

        //this whole thing is to avoid using ImageIO.read(byte[])
        //since that method writes the array to a temp file
        //and saves it in some kind of cache which is slow(er)
        //this avoids that and eliminates the need for temp files
        //for a little it of a speed improvement and less hassle managing files
        //or worrying if the caller wants to use the ImageIO cache and has set it properly.
        try(ImageInputStream input = new ByteArrayImageInputStream(file)) {
            // Get the reader
            Iterator<ImageReader> readers;
            if(isPng(file)){
                readers = ImageIO.getImageReadersByFormatName("png");
            }else if(isTiff(file)){
                readers = ImageIO.getImageReadersByFormatName("tiff");
            }else{
                readers = ImageIO.getImageReaders(input);
            }

            if (!readers.hasNext()) {
                throw new IOException("No reader found for format provided in byte array");
            }

            ImageReader reader = readers.next();

            try {
                reader.setInput(input);

                // Finally read the image, using settings from param
                return  decode(reader.read(0));

                // ...
            }
            finally {
                // Dispose reader in finally block to avoid memory leaks
                reader.dispose();
            }
        }
//        return decode(ImageIO.read(new ByteArrayImageInputStream(file)));
    }
    public static BufferedImage grayscale (File file) throws IOException {

        return decode(ImageIO.read(file));
    }

    public static BufferedImage decode (File file) throws IOException {
//    	System.out.println(file.getName());
        return ImageIO.read(file);
    }

    public static BufferedImage fuse (BufferedImage img0, BufferedImage img1) {
	if (img0.getWidth() != img1.getWidth() 
	    || img0.getHeight() != img1.getHeight()) {
	    throw new IllegalArgumentException 
		("Images are not of the same dimension");
	}

	int width = img0.getWidth(), height = img0.getHeight();

	BufferedImage fused = new BufferedImage 
	    (width, height, BufferedImage.TYPE_INT_RGB);

	/*
	BufferedImage green = new BufferedImage 
	    (width, height, BufferedImage.TYPE_INT_RGB);
	BufferedImage blue = new BufferedImage 
	    (width, height, BufferedImage.TYPE_INT_RGB);
	for (int x = 0; x < width; ++x) {
	    for (int y = 0; y < height; ++y) {
		green.setRGB(x, y, (img0.getRGB(x,y) & 0xff) << 8);
		blue.setRGB(x, y, grayscale(img1.getRGB(x,y)) & 0xff);
	    }
	}
	
	Graphics2D g = fused.createGraphics();
	g.drawImage(green, 0, 0, null);
	AlphaComposite c = AlphaComposite.getInstance
	    (AlphaComposite.SRC_OVER, .35f);
	g.setComposite(c);
	g.drawImage(blue, 0, 0, null);
	g.dispose();
	*/	
	int p, q, g, b;
	for (int x = 0; x < width; ++x) {
	    for (int y = 0; y < height; ++y) {
		p = img0.getRGB(x, y) & 0xff00;
		q = img1.getRGB(x, y) & 0xff;
		fused.setRGB(x, y, p | q);
	    }
	}

	return fused;
    }

    protected static int grayscale (int rgb) {
	int r = (rgb & 0x00ffffff) >> 16;
	int g = (rgb & 0x0000ffff) >> 8;
	int b = (rgb & 0x000000ff);
	/*
	return (int)(0.3*r + 0.59*g + 0.11*b + 0.5);
	*/
	return Math.max(r, Math.max(g, b));
    }

}
