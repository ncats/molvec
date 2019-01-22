package tripod.molvec.ui;

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
import java.util.Base64;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import javax.imageio.ImageIO;

import tripod.molvec.algo.Tuple;

public abstract class RasterBasedCosineSCOCR implements SCOCR{
	Set<Character> _alphabet;
	Map<Character, List<RasterChar>> charVal = new HashMap<Character, List<RasterChar>>();
	//Map<Character, Rectangle2D> WH_RATIO = new HashMap<Character, Rectangle2D>();

	static int DEF_WIDTH = 20;
	static int DEF_HEIGHT = 20;
	protected static BufferedImage bi = new BufferedImage(DEF_WIDTH, DEF_HEIGHT,
			BufferedImage.TYPE_BYTE_GRAY);
	
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

	protected static void blurbmap(int[][] bmap) {
		vblur(bmap, 2);
		hblur(bmap, 2);
	}
	
	
	public void saveModel(File dir){
		if(!dir.exists())dir.mkdirs();
		if(!dir.isDirectory()) throw new IllegalStateException("Must specify a directory to save model");
		
		for(Character ch : charVal.keySet()){
			AtomicInteger ai = new AtomicInteger(0);
			//RasterChar rc=charVal.get(ch);
			
			

			
			charVal.get(ch).stream()
	                        .forEach(rc->{
	                        		Rectangle2D rat=rc.rect;
	                    			String recSTR=rat.getWidth() + "x" + rat.getHeight();
	                        		int i=ai.getAndIncrement();
	                        		try{
	                        				Raster rs=makeRaster(rc.data);
	                        		 		BufferedImage image = new BufferedImage(DEF_WIDTH, DEF_HEIGHT, BufferedImage.TYPE_BYTE_GRAY);
	                        		        image.setData (rs);
	                        		        String fullpath = dir.getAbsolutePath() + "/" + "glyph." + ((int)ch.charValue()) +"." +  i + "." + recSTR + ".png";
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
	
	public static class RasterChar{
		public int[][] data;
		public Rectangle2D rect;
		
		public RasterChar(int[][] dat, Rectangle2D rect){
			this.data=dat;
			this.rect=rect;
		}
		
		//hacky way to do this for now
		public String rawDataAsString(){
			return data.length + "x" + data[0].length + "\n" +
				   rect.getWidth() + "x" + rect.getHeight() + "\n" 
		           + Arrays.stream(data)
					     .map(line->Arrays.stream(line)
					    		           .mapToObj(i->Integer.toString(i, 16))
					    		           .collect(Collectors.joining(",")))
					     .collect(Collectors.joining("\n"));      
		}
		
		public RasterChar readDataFromString(String raw){
			String[] rawlines = raw.split("\n");
			int[] dim =Arrays.stream(rawlines[0].split("x")).mapToInt(s->Integer.parseInt(s)).toArray();
			double[] dimB =Arrays.stream(rawlines[1].split("x")).mapToDouble(s->Double.parseDouble(s)).toArray();
			this.data=Arrays.stream(rawlines)
							.skip(2)
							.map(l->Arrays.stream(l.split(","))
			    		        .mapToInt(s->Integer.parseInt(s,16))
			    		        .toArray())
							.toArray(i->new int[dim[0]][dim[1]]);
			this.rect=new Rectangle2D.Double(0,0,dimB[0],dimB[1]);
			return this;
		}		
	}
	
	
	public void debug(){
		charVal.entrySet()
		.stream()
		.map(Tuple::of)
		.map(Tuple.vmap(rcs->rcs.stream().map(rc->Base64.getEncoder().encodeToString(rc.rawDataAsString().getBytes()))))
		.flatMap(t->t.v().map(rc->t.k() + "\t" + rc))
		.forEach(l->{
			System.out.println(l);
		})
		;
	
		//Base64.getEncoder()
	}
	
	
	
	
	public abstract void getBitmapsForChar(Character c, Consumer<RasterChar> rconsumer);
	

	private void makeAlphabetMaps() {
		for (char c : _alphabet) {
			getBitmapsForChar(c,rc->{
				charVal.computeIfAbsent(c, (k)->new ArrayList<>())
				       .add(rc);
			});
		}
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
		for (RasterChar rc: charVal.get(c)) {
			int[][] cM = rc.data;
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
			whrat = rc.rect.getWidth()/rc.rect.getHeight();
			if (whrat > 1)
				whrat = 1 / whrat;
			whrat = 1 - Math.pow(1 - whrat, 2);

			double tcor = whrat * cor / (Math.sqrt(total) * Math.sqrt(totalC));
			maxCor = Math.max(tcor, maxCor);
		}
		return maxCor;
	}

}
