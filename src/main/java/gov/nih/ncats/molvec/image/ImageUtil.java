package gov.nih.ncats.molvec.image;

import java.io.*;
import java.awt.image.*;

import javax.imageio.*;

import java.util.logging.Logger;

import com.sun.media.jai.codec.*;


public class ImageUtil implements TiffTags {
    private static final Logger logger = Logger.getLogger
	(ImageUtil.class.getName());
    private static BufferedImage decodeTIFF (byte[] file) throws IOException {
        return decodeTiff( ImageCodec.createImageDecoder("TIFF", new ByteArraySeekableStream(file), new TIFFDecodeParam ()));
    }

    private static BufferedImage decodeTIFF (InputStream file) throws IOException {
        return decodeTiff(ImageCodec.createImageDecoder("TIFF", file, new TIFFDecodeParam ()));
    }

    private static BufferedImage decodeTiff(ImageDecoder decoder) throws IOException{
        int ndirs = decoder.getNumPages();
        TIFFDirectory tif = new TIFFDirectory
                (decoder.getInputStream(), 0);
        TIFFField[] fields = tif.getFields();

        double width = 0, height = 0;
        String unit = "";
        double xres = 0., yres = 0.;
        double rows = -1;
        int photometric = -1, bpp = -1;
        for (int j = 0; j < fields.length; ++j) {
            TIFFField f = fields[j];
            int tag = f.getTag();
            try {
                switch (tag) {
                    case TAG_RESOLUTIONUNIT:
                    {
                        int u = f.getAsInt(0);
                        if (u == RESOLUTIONUNIT_NONE) {
                        }
                        else if (u == RESOLUTIONUNIT_INCH) {
                            unit = "in";
                        }
                        else if (u == RESOLUTIONUNIT_CENT) {
                            unit = "cm";
                        }
                    }
                    break;

                    case TAG_XRESOLUTION:
                        xres = f.getAsFloat(0);
                        break;

                    case TAG_YRESOLUTION:
                        yres = f.getAsFloat(0);
                        break;

                    case TAG_ROWSPERSTRIP:
                        rows = f.getAsFloat(0);
                        break;

                    case TAG_PHOTOMETRIC:
                        photometric = f.getAsInt(0);
                        break;

                    case TAG_BITSPERSAMPLE:
                        bpp = f.getAsInt(0);
                        break;

                    case TAG_IMAGEWIDTH:
                        width = f.getAsFloat(0);
                        break;

                    case TAG_IMAGELENGTH:
                        height = f.getAsFloat(0);
                        break;
                }
            }
            catch (Exception ex) {
                logger.warning("## TIFF decoder tag="
                        +tag+"; "+ex.getMessage());
            }
        }

	/*
	if (xres > 0) {
	    width /= xres;
	}
	if (yres > 0) {
	    height /= yres;
	}
	*/

        logger.info( tif.toString() + " has " + ndirs + " image; width="+width
                +" height="+height +" xres="+xres+unit
                +" yres="+yres+unit+" bpp="+bpp
                +" photometric="+photometric+" rows="+rows);

        RenderedImage decodedImage = decoder.decodeAsRenderedImage();
        //ImageIO.write(decodedImage, "png", new File("tmp.png"));
        Raster raster = decodedImage.getData();
	/*
	if (raster.getNumBands() > 1) {
	    throw new IllegalArgumentException
		("Sorry, can't support multiband image at the moment!");
	}
	*/

        logger.info("sample model: nbands="+raster.getNumBands()
                +" "+raster.getSampleModel().getClass());
        /*
	MultiPixelPackedSampleModel packed =
	    (MultiPixelPackedSampleModel)raster.getSampleModel();
	logger.info("scanline: "+packed.getScanlineStride());
	logger.info("bit stride: "+packed.getPixelBitStride());
        */

        Grayscale grayscale = new Grayscale (raster);
        BufferedImage big=grayscale.getImage();
        //ImageIO.write(decodedImage, "png", new File("tmp-gray.png"));
        return big;
    }

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

        logger.info("## dynamic range: "+(max-min)+" color model: "+bi.getColorModel());
        
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
        if(isTiff(file)) {
            try {
                return decodeTIFF(file);
            } catch (Exception ex) {
                logger.info("## bytearray not a TIFF image ("
                        + ex.getMessage() + "); trying generic decoding... ");

            }
        }
        return decode(ImageIO.read(new ByteArrayInputStream(file)));
    }
    public static BufferedImage grayscale (File file) throws IOException {
        BufferedInputStream in = new BufferedInputStream(new FileInputStream(file));
        try(PushbackInputStream pushbackInputStream = new PushbackInputStream(in,2);
        ) {
           if (isTiff(pushbackInputStream)) {
                try {
                    return decodeTIFF(pushbackInputStream);
                } catch (Exception ex) {
                    logger.info("## " + file.getName() + " not a TIFF image ("
                            + ex.getMessage() + "); trying generic decoding... ");
                }
           }
        }
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
