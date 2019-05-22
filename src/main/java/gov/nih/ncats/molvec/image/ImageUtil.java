package gov.nih.ncats.molvec.image;

import java.io.*;
import java.awt.image.*;

import javax.imageio.*;
import java.util.Map;
import java.util.logging.Logger;



public class ImageUtil implements TiffTags {
    private static final Logger logger = Logger.getLogger
	(ImageUtil.class.getName());


    public static BufferedImage decode (BufferedImage bi) {
        Raster raster = bi.getData();

        int bands=raster.getNumBands();
    	int max = 0;
    	int min = Integer.MAX_VALUE;
    	for (int i = 0; i < raster.getWidth(); ++i) {
    	    for (int j = 0; j < raster.getHeight(); ++j) {
    	    	for(int k=0;k<bands;k++){
	                int pixel = raster.getSample(i, j, k);
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

    private static boolean isTiff(PushbackInputStream pushbackInputStream) throws IOException{
       byte[] magicNumber = new byte[2];
        int len = pushbackInputStream.read(magicNumber,0,2);
//        System.out.println("MAGIC NUMBER = " + Arrays.toString(magicNumber));
        pushbackInputStream.unread(magicNumber,0,len);
        if(len !=2){
            return false;
        }

        //0x4949 or 0x4d4d
        return ((magicNumber[0] == 0x49 && magicNumber[1] == 0x49) || (magicNumber[0] == 0x4d && magicNumber[1] == 0x4d));
    }
    private static boolean isTiff(byte[] f) throws IOException{

        //0x4949 or 0x4d4d
        return ((f[0] == 0x49 && f[1] == 0x49) || (f[0] == 0x4d && f[1] == 0x4d));
    }
    public static BufferedImage grayscale (BufferedImage bi) {
        Grayscale grayscale = new Grayscale (bi.getData());
    	Raster raster = grayscale.getRaster();
    	BufferedImage image = new BufferedImage 
    	    (raster.getWidth(), raster.getHeight(), 
    	     BufferedImage.TYPE_BYTE_GRAY);
    	image.setData(raster);
        return image;
    }

    public static BufferedImage grayscale (byte[] file) throws IOException {

        return decode(ImageIO.read(new ByteArrayInputStream(file)));
    }
    public static BufferedImage grayscale (File file) throws IOException {

        return decode(ImageIO.read(file));
    }

    public static BufferedImage decode (File file) throws IOException {
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

    public static void main (String[] argv) throws Exception {
	if (argv.length == 0) {
	    System.out.println("Usage: ImageUtil IMAGE.tif");
	    System.exit(1);
	}

	BufferedImage img = decode (new File (argv[0]));
    }
}
