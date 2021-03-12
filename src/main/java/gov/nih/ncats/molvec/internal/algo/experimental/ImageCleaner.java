package gov.nih.ncats.molvec.internal.algo.experimental;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.awt.image.BufferedImageOp;
import java.awt.image.ConvolveOp;
import java.awt.image.Kernel;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

public class ImageCleaner {
	public static BufferedImage rotateCw(BufferedImage img) {
		int width = img.getWidth();
		int height = img.getHeight();
		BufferedImage newImage = new BufferedImage(height, width, img.getType());
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				newImage.setRGB(height - 1 - j, i, img.getRGB(i, j));
			}
		}
		return newImage;
	}
	
	public static BufferedImage preCleanImageResize(BufferedImage biIn, double scale, boolean qblur, boolean rotateWidest) throws IOException{

		if(biIn.getWidth()>biIn.getHeight()){
			rotateWidest=!rotateWidest;
		}
		if(rotateWidest){
			biIn=rotateCw(biIn);
		}
		 
		if(qblur){	
		    float[] matrix = {
		            0.25f, 0.25f, 
		            0.25f, 0.25f, 
		        };
	
	        BufferedImageOp op = new ConvolveOp( new Kernel(2, 2, matrix) );
	        BufferedImage biIn2= new BufferedImage(biIn.getWidth(),
	        		biIn.getHeight(),BufferedImage.TYPE_3BYTE_BGR);
	        
	        biIn2 = op.filter(biIn, biIn2);
	        biIn=biIn2;

	        for (int x = 0; x < biIn.getWidth(); x++) {
	        	for (int y = 0; y < biIn.getHeight(); y++) {
	        		int rgba = biIn.getRGB(x, y);
	        		Color col = new Color(rgba, false);
	        		col = new Color(col.getRed(),
	        				col.getRed(),
	        				col.getRed());
	        		if(col.getRed()>253 || x<3||y<3||x>biIn.getWidth()-3||y>biIn.getHeight()-3){
	        			biIn.setRGB(x, y, Color.WHITE.getRGB());	
	        		}else{
	        			biIn.setRGB(x, y, Color.BLACK.getRGB());
	        		}
	        	}
	        }

		}
        
		int nwidth=(int) (biIn.getWidth() *scale);
		int nheight=(int) (biIn.getHeight() *scale);
		BufferedImage outputImage=null;
		
		
        // creates output image
        outputImage = new BufferedImage(nwidth,
                nheight,BufferedImage.TYPE_3BYTE_BGR);
 
        // scales the input image to the output image
        Graphics2D g2d = outputImage.createGraphics();
        
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
        

    	g2d.setRenderingHint(RenderingHints.KEY_INTERPOLATION,
        	       RenderingHints.VALUE_INTERPOLATION_NEAREST_NEIGHBOR);

        g2d.scale(scale, scale);
        g2d.drawImage(biIn, 0, 0, null);
        g2d.dispose();

        for (int x = 0; x < outputImage.getWidth(); x++) {
            for (int y = 0; y < outputImage.getHeight(); y++) {
                int rgba = outputImage.getRGB(x, y);
                Color col = new Color(rgba, false);
                col = new Color(col.getRed(),
                		col.getRed(),
                		col.getRed());
                outputImage.setRGB(x, y, col.getRGB());
                
            }
        }

//        ImageIO.write(outputImage, "png", new File("tmp.2.png"));
		
		return outputImage;
	}

}
