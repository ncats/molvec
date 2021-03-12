package gov.nih.ncats.molvec.ui;

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
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.imageio.ImageIO;

import gov.nih.ncats.molvec.internal.image.Bitmap;
import gov.nih.ncats.molvec.internal.image.Bitmap.BitmapScaled;
import gov.nih.ncats.molvec.internal.util.CachedSupplier;
import gov.nih.ncats.molvec.internal.algo.Tuple;

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
			
			

			
			charVal.get(ch)
//                    .stream()
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

		private static final Pattern LINE_SEP_PATTERN = Pattern.compile("\n");
		private static final Pattern COORD_SEP_PATTERN = Pattern.compile("x");
		private static final Pattern VALUE_SEP_PATTERN = Pattern.compile(",");

		/**
		 * Cache for Ints to their hex String representation this
		 * is used to speed up performance of deserializing the rasters
		 * in {@link #readDataFromString(String)}
		 */
		private static Map<String, Integer> INT_16_CACHE = new ConcurrentHashMap<>();
		private static Map<String, Integer> INT_CACHE = new ConcurrentHashMap<>();
		private static Map<String, Double> DOUBLE_CACHE = new ConcurrentHashMap<>();

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
		static RasterChar  parseFromString(String raw){
			RasterChar rs = new RasterChar(null, null);
			return rs.readDataFromString(raw);
		}
		private RasterChar readDataFromString(String raw){
			String[] rawlines = LINE_SEP_PATTERN.split(raw);

			String[] dimLine = COORD_SEP_PATTERN.split(rawlines[0]);
			int x = INT_CACHE.computeIfAbsent(dimLine[0], v-> Integer.parseInt(v));
			int y = INT_CACHE.computeIfAbsent(dimLine[1], v-> Integer.parseInt(v));


			String[] rectLine = COORD_SEP_PATTERN.split(rawlines[1]);
			double rectX = DOUBLE_CACHE.computeIfAbsent(rectLine[0], v->Double.parseDouble(v));
			double rectY = DOUBLE_CACHE.computeIfAbsent(rectLine[1], v->Double.parseDouble(v));

			this.data = new int[x][y];
			for(int i=2; i< rawlines.length; i++){
				int[] row = new int[y];
				int j=0;
				for(String v : VALUE_SEP_PATTERN.split(rawlines[i])){
					row[j] = INT_16_CACHE.computeIfAbsent(v, value ->Integer.parseInt(value, 16));
					j++;
				}
				data[i-2] = row;

			}

			this.rect=new Rectangle2D.Double(0,0,rectX,rectY);
			return this;
		}		

		
		public static RasterChar from(Bitmap bm, int wid, int hit){
			int[][] bmap = new int[wid][hit];
			
			double scalex = ((double)wid)/bm.width();
			double scaley = ((double)hit)/bm.height();
			
			bm.getXYOnPoints()
			  .map(xy->new int[]{(int)(xy[0]*scalex),(int)(xy[1]*scaley)})
			  .forEach(xy->{
				  bmap[xy[0]][xy[1]]=1;
			  });  
			Rectangle2D rect = new Rectangle2D.Double(0, 0, bm.width(),bm.height());
			return new RasterChar(bmap,rect);
		}
		
		public RasterChar blur(int r){
			vblur(this.data, r);
			hblur(this.data, r);
			return this;
		}
		
		public static RasterChar fromDefault(Bitmap bm){
			return from(bm,DEF_WIDTH,DEF_HEIGHT);
		}
		
		
		private CachedSupplier<Integer> _totalC = CachedSupplier.of(()->{
			int totalc = 0;
			for(int i=0;i<data.length;i++){
				for(int j=0;j<data[i].length;j++){
					totalc+=data[i][j]*data[i][j];
				}
			}
			return totalc;
		});
		
		private CachedSupplier<int[]> _totalCH = CachedSupplier.of(()->{
			int[] vec = new int[data.length];
			for(int i=0;i<data.length;i++){
				int c=0;
				for(int j=0;j<data[i].length;j++){
					c+=data[i][j]*data[i][j];
				}
				vec[i]=c;
			}
			return vec;
		});
		
		private CachedSupplier<int[]> _totalCV = CachedSupplier.of(()->{
			int[] vec = new int[data[0].length];
			for(int i=0;i<data.length;i++){
				for(int j=0;j<data[i].length;j++){
					vec[j]+=data[i][j]*data[i][j];
				}
			}
			return vec;
		});
		
		public int getSqLength(){
			return Arrays.stream(_totalCH.get()).sum();
		}
		public int[] getSqLengthHorizontal(){
			return _totalCH.get();
		}
		
		public int[] getSqLengthVertical(){
			return _totalCV.get();
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
		_alphabet = charSet.stream().collect(Collectors.toSet());
		
		//hacky control chars
		_alphabet.add('~');
		_alphabet.add('$');
		_alphabet.add('!');
		_alphabet.add('%');
		_alphabet.add('`');
		_alphabet.add('(');
		_alphabet.add(')');
//		_alphabet.add('?');
		makeAlphabetMaps();
	}

	@Override
	public Set<Character> getAlphabet() {
		return _alphabet;
	}

	@Override
	public Map<Character, Number> getRanking(Bitmap r) {
		BitmapScaled bms=r.getScaled(DEF_WIDTH, DEF_HEIGHT);
		
		return _alphabet
				.parallelStream()
				.collect(Collectors.toMap(Function.identity(), c2->correlation(bms, c2)));
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
	
	
	
	private double correlation(BitmapScaled bms, Character c){

		List<RasterChar> rcl = charVal.get(c);
		if(rcl==null){
			return 0;
		}
		List<int[]> xys = bms.xys;
		int[][] ccount = bms.ccount;
		int total = bms.tcount;
		int twidth = bms.twidth;
		int theight = bms.theight;

		double maxCor = Double.MIN_VALUE;
		for (RasterChar rc: rcl) {
			int[][] cM = rc.data;
			double cor = 0;
			int sum=0;
			for(int[] xy : xys){
				cor +=cM[xy[0]][xy[1]]*xy[2];
			}

			
			for(int i=0;i<DEF_WIDTH;i++){
				for(int j=0;j<DEF_HEIGHT;j++){
					int val=cM[i][j];
					sum+=val*val*ccount[i][j];
				}
			}
			int totalC=sum;
			
			double whrat = (double) twidth / (double) theight;
			whrat = (rc.rect.getWidth()/rc.rect.getHeight()) / whrat;
			if (whrat > 1) {
				whrat = 1 / whrat;
			}
			whrat = 1 - Math.pow(1 - whrat, 2);

			double tcor = whrat * cor / (Math.sqrt(total) * Math.sqrt(totalC));
			maxCor = Math.max(tcor, maxCor);
		}
		return maxCor;
	}

//	public double correlation(Bitmap test, Character c) {
//		int twidth = test.width();
//		int theight = test.height();
//
//		int[][] ccount = new int[DEF_WIDTH][DEF_HEIGHT];
//
//		List<int[]> xys=test.getXYOnPoints()
//						    .map(xy->{
//						    	int cx = (xy[0] * DEF_WIDTH) / twidth;
//						    	int cy = (xy[1] * DEF_HEIGHT) / theight;
//						    	return new int[]{cx,cy};
//						    })
//						    .collect(Collectors.toList());
//
//		for(int i=0;i<twidth;i++){
//			for(int j=0;j<theight;j++){
//			   	int cx = (i * DEF_WIDTH) / twidth;
//		    	int cy = (j * DEF_HEIGHT) / theight;
//		    	ccount[cx][cy]++;
//			}
//		}
//		return correlation(xys,ccount,xys.size(),twidth,theight,c);
//	}
	
	

}
