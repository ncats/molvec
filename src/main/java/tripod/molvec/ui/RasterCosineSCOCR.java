package tripod.molvec.ui;

import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.font.GlyphVector;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.DataBuffer;
import java.awt.image.MultiPixelPackedSampleModel;
import java.awt.image.Raster;
import java.awt.image.SampleModel;
import java.awt.image.WritableRaster;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

import javax.imageio.ImageIO;

public class RasterCosineSCOCR extends SCOCR {
	Set<Character> _alphabet;
	Map<Character, Collection<int[][]>> charVal = new HashMap<Character, Collection<int[][]>>();
	Map<Character, Double> WH_RATIO = new HashMap<Character, Double>();

	static int DEF_WIDTH = 20;
	static int DEF_HEIGHT = 20;
	private static BufferedImage bi = new BufferedImage(DEF_WIDTH, DEF_HEIGHT,
			BufferedImage.TYPE_BYTE_GRAY);
	private List<Font> fontList = new ArrayList<Font>();
	
	public static List<Font> DEFAULT_FONTS(){
		return SANS_SERIF_FONTS();
	}
	
	public static List<Font> SERIF_FONTS(){
		return Arrays.asList(new Font[] {
				new Font(Font.SERIF, Font.PLAIN, 8),
				new Font(Font.SERIF, Font.BOLD, 8) 
				});
	}
	public static List<Font> SANS_SERIF_FONTS(){
		return Arrays.asList(new Font[] {
				new Font(Font.SANS_SERIF, Font.BOLD, 8),
				new Font(Font.SANS_SERIF, Font.PLAIN, 8),
				new Font("DejaVu Sans", Font.BOLD, 8),
				new Font("Khmer OS", Font.PLAIN, 8),
				new Font("DejaVu Sans", Font.PLAIN, 8),
				new Font("Khmer OS", Font.BOLD, 8) 
				});
	}
	
	public RasterCosineSCOCR(List<Font> flist) {
		fontList.addAll(flist);
	}
	
	public RasterCosineSCOCR() {
		this(Arrays.asList(new Font[] {
				new Font(Font.SANS_SERIF, Font.BOLD, 8),
				new Font(Font.SANS_SERIF, Font.PLAIN, 8),
				new Font("DejaVu Sans", Font.BOLD, 8),
				new Font("Khmer OS", Font.PLAIN, 8),
				new Font("DejaVu Sans", Font.PLAIN, 8),
				new Font("Khmer OS", Font.BOLD, 8), 
				new Font(Font.SERIF, Font.PLAIN, 8),
				new Font(Font.SERIF, Font.BOLD, 8) 
				}));
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
	private static void vblur(int[][] bmap, int rad) {
		int sofar = 0;
		List<Integer> added = new ArrayList<Integer>();
		for (int i = 0; i < bmap.length; i++) {
			for (int j = 0; j < rad * 2 + bmap[0].length; j++) {
				// System.out.println(sofar + ":" + added.size());
				if (j < bmap.length) {
					added.add(bmap[i][j]);
					sofar += bmap[i][j];
				} else {
					if (j < bmap.length + rad) {
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
	private static void hblur(int[][] bmap, int rad) {
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

	private static void blurbmap(int[][] bmap) {
		vblur(bmap, 2);
		hblur(bmap, 2);

	}

	/**
	 * 
	 * @param c
	 *            Character to generate bitmap for
	 * @param f
	 *            Font to use in generation
	 * @param toRet
	 *            Integer array to put result into
	 * @return Returns the width/height of the character
	 */
	private static double getBitmap(Character c, Font f, int[][] toRet) {
		if (toRet == null) {
			toRet = new int[DEF_WIDTH][DEF_HEIGHT];
		}

		Graphics2D g = (Graphics2D) bi.getGraphics();
		GlyphVector gv = f.createGlyphVector(g.getFontRenderContext(),
				new char[] { c });
		Shape outline=gv.getOutline();
		Rectangle2D r = gv.getVisualBounds();

		double w = r.getWidth();
		double h = r.getHeight();

		Point2D testPoint = new Point2D.Double();
		
		for (int j = 0; j < DEF_HEIGHT; j++) {
			for (int i = 0; i < DEF_WIDTH; i++) {
				testPoint.setLocation(r.getMinX() + (w * i) / DEF_WIDTH, 
									  r.getMinY() + (h * j) / DEF_HEIGHT);
				toRet[i][j] = outline.contains(testPoint) ? 1
						: 0;
			}
		}
		
		
		        
		blurbmap(toRet);
		return w / h;
	}
	
	
	private void saveModel(File dir){
		if(!dir.exists())dir.mkdirs();
		if(!dir.isDirectory()) throw new IllegalStateException("Must specify a directory to save model");
		
		for(Character ch : charVal.keySet()){
			AtomicInteger ai = new AtomicInteger(0);
			double rat=WH_RATIO.get(ch);
			
			charVal.get(ch).stream()
	                        .forEach(data->{
	                        	int i=ai.getAndIncrement();
	                        	try{
	                        		Raster rs=makeRaster(data);
	                        		 BufferedImage image = new BufferedImage
	                        		            (DEF_WIDTH, DEF_HEIGHT, BufferedImage.TYPE_BYTE_GRAY);
	                        		        image.setData (rs);
	                        		        
	                        		        String fullpath = dir.getAbsolutePath() + "/" + "glyph." + ((int)ch.charValue()) +"." +  i + "." + rat + ".png";
	                        		        File file = new File(fullpath);
	                        		        ImageIO.write (image, "png", file);
	                        		}catch(Exception e){
	                        			e.printStackTrace();
	                        		}
	                        });
		}
		
	}
	
	private static Raster makeRaster(int[][] data){
		int width=data.length;
		int height=data[0].length;
		int scanline=1*width;
		SampleModel sampleModel = new MultiPixelPackedSampleModel(DataBuffer.TYPE_BYTE, width, height, 8, scanline, 0);
		WritableRaster raster =
		            Raster.createWritableRaster (sampleModel, null);
		        for (int y = 0; y < height; ++y) {
		            for (int x = 0; x < width; ++x) {
		                // the default IndexColorModel is 0 for black and 1 white
		                raster.setSample
		                    (x, y, 0, data[x][y]*8);
		            }
		        }
		        return raster;
	}
	

	private void makeAlphabetMaps() {
		for (Font f : fontList) {
			for (char c : _alphabet) {
				int[][] charMat = new int[DEF_WIDTH][DEF_HEIGHT];
				double wh = getBitmap(c, f, charMat);
				WH_RATIO.put(c, wh);
				Collection<int[][]> bmap = charVal.get(c);
				if (bmap == null) {
					bmap = new ArrayList<int[][]>();
					charVal.put(c, bmap);
				}
				bmap.add(charMat);
			}
		}
		//saveModel(new File("tmp1"));
	}
	
	
	

	@Override
	public void setAlphabet(Set<Character> charSet) {
		_alphabet = charSet;
		makeAlphabetMaps();
	}

	@Override
	public Set<Character> getAlphabet() {
		return _alphabet;
	}

	@Override
	public Map<Character, Number> getRanking(Raster r) {
		Map<Character, Number> myMap = new HashMap<Character, Number>();
		int[][] bmap = new int[r.getHeight()][r.getWidth()];
		for (int i = 0; i < r.getHeight(); i++) {
			r.getPixels(0, i, r.getWidth(), 1, bmap[i]);
		}
		int[][] bmap2 = new int[r.getWidth()][r.getHeight()];
		int twidth = bmap2.length;
		int theight = bmap2[0].length;
		for (int j = 0; j < theight; j++) {
			for (int i = 0; i < twidth; i++) {
				bmap2[i][j] = (bmap[j][i] > 0) ? 0 : 1;
			}
		}
		// debugPrintBmap(bmap2);
		for (char c2 : _alphabet) {
			myMap.put(c2, correlation(bmap2, c2));
		}
		return myMap;
	}

	public static void debugPrintBmap(int[][] test) {
		int twidth = test.length;
		int theight = test[0].length;

		for (int j = 0; j < theight; j++) {
			for (int i = 0; i < twidth; i++) {
				System.out.print(test[i][j]);
			}
			System.out.println();
		}

	}

	public double correlation(int[][] test, Character c) {
		double maxCor = Double.MIN_VALUE;
		for (int[][] cM : charVal.get(c)) {
			int twidth = test.length;
			int theight = test[0].length;
			int total = 0;
			int totalC = 0;
			double cor = 0;
			for (int i = 0; i < twidth; i++) {
				for (int j = 0; j < theight; j++) {
					int cx = Math.round((i * DEF_WIDTH) / twidth);
					int cy = Math.round((j * DEF_HEIGHT) / theight);
					int val = cM[cx][cy];
					cor += val * test[i][j];
					total += test[i][j] * test[i][j];
					totalC += val * val;
				}
			}
			double whrat = (double) twidth / (double) theight;
			whrat = WH_RATIO.get(c) / whrat;
			if (whrat > 1)
				whrat = 1 / whrat;
			whrat = 1 - Math.pow(1 - whrat, 2);

			double tcor = whrat * cor / (Math.sqrt(total) * Math.sqrt(totalC));
			maxCor = Math.max(tcor, maxCor);
		}
		return maxCor;
	}

}
