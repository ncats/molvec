package gov.nih.ncats.molvec.algo;

import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.ColorConvertOp;
import java.awt.image.Raster;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import gov.nih.ncats.molvec.Bitmap;
import gov.nih.ncats.molvec.Bitmap.WedgeInfo;
import gov.nih.ncats.molvec.CachedSupplier;
import gov.nih.ncats.molvec.ui.FontBasedRasterCosineSCOCR;
import gov.nih.ncats.molvec.ui.SCOCR;
import gov.nih.ncats.molvec.ui.StupidestPossibleSCOCRSansSerif;
import gov.nih.ncats.molvec.ui.StupidestPossibleSCOCRSerif;
import gov.nih.ncats.molvec.util.CompareUtil;
import gov.nih.ncats.molvec.util.ConnectionTable;
import gov.nih.ncats.molvec.util.ConnectionTable.Edge;
import gov.nih.ncats.molvec.util.ConnectionTable.Node;
import gov.nih.ncats.molvec.util.ConnectionTable.Ring;
import gov.nih.ncats.molvec.util.GeomUtil;
import gov.nih.ncats.molvec.util.GeomUtil.BoundingBox;
import gov.nih.ncats.molvec.util.GeomUtil.LineWrapper;
import gov.nih.ncats.molvec.util.RunningAverage;

public class StructureImageExtractor {

	static final SCOCR OCR_DEFAULT=new StupidestPossibleSCOCRSansSerif();
	//static final SCOCR OCR_DEFAULT=new FontBasedRasterCosineSCOCR(FontBasedRasterCosineSCOCR.SANS_SERIF_FONTS());
	//static final SCOCR OCR_BACKUP=new FontBasedRasterCosineSCOCR(FontBasedRasterCosineSCOCR.SERIF_FONTS())
	static final SCOCR OCR_BACKUP=new StupidestPossibleSCOCRSerif()
			.adjustWeights(t->{
				double ov=t.v().doubleValue();
				Tuple<Character,Number> ret=t;
				if("N".equals(t.k()+"")){
					ret= Tuple.of(t.k(),(Number)(Math.max(1-(1-ov)*1.1,0)));
				}
				return ret;
			});


	static final SCOCR OCR_ALL=new FontBasedRasterCosineSCOCR();

	static{
		Set<Character> alpha=SCOCR.SET_COMMON_CHEM_ALL();
		alpha.add(Character.valueOf('/'));
		alpha.add(Character.valueOf('\\'));
		OCR_DEFAULT.setAlphabet(alpha);
		OCR_BACKUP.setAlphabet(alpha);

		Set<Character> alphaAll=SCOCR.SET_COMMON_PUCTUATION();
		alphaAll.addAll(SCOCR.SET_ALPHANUMERIC());
		OCR_ALL.setAlphabet(alphaAll);

		//((FontBasedRasterCosineSCOCR)OCR_BACKUP).debug();
	}

	private Bitmap bitmap; // original bitmap
	private Bitmap thin; // thinned bitmap


	private List<Shape> polygons;
	private List<Shape> rescueOCRShapes;
	public List<Shape> getRescueOCRShapes() {
		return rescueOCRShapes;
	}

	private List<LineWrapper> lines;
	private List<LineWrapper> linesJoined;
	private List<Tuple<Line2D,Integer>> linesOrder;    
	private Map<Shape,List<Tuple<Character,Number>>> ocrAttempt = new ConcurrentHashMap<>();
	private Map<Shape,String> bestGuessOCR = new LinkedHashMap<>();
	

	private ConnectionTable ctab;
	private List<ConnectionTable> ctabRaw = new ArrayList<ConnectionTable>();

	
	private final int MAX_OCR_FULL_REPEATS=6;
	private final int MAX_REPS = 2;
	
	
	private final double INITIAL_MAX_BOND_LENGTH=Double.MAX_VALUE;
	private final double MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE = 1/3.6;
	private final double MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE_AFTER_SPLIT = 1/6.0;
	private final double MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE_INITIAL = 1/2.9;

	private final double MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL = 1.3;
	private final double MIN_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL = 0.4;

	
	private final boolean REMOVE_NONSENSE_OCR_LINES = false;

	private final double MAX_TOLERANCE_FOR_DASH_BONDS = 3.0;
	private final double MAX_TOLERANCE_FOR_SINGLE_BONDS = 0.4;

	private final double OCRcutoffCosine=0.65;
	private final double OCRcutoffCosineRescueInitial=0.55;
	private final double OCRcutoffCosineRescue=0.50;
	

	private final double WEDGE_LIKE_PEARSON_SCORE_CUTOFF=.60;
	private final double WEDGE_LIKE_PEARSON_SCORE_CUTOFF_DOUBLE=.80;

	private final double MAX_BOND_TO_AVG_BOND_RATIO_TO_KEEP = 1.8;
	private final double MAX_DISTANCE_BEFORE_MERGING_NODES = 4.0;
	private final double maxRatioForIntersection = 1.10;
	private final double maxPerLineDistanceRatioForIntersection = 1.6;
	private final double minPerLineDistanceRatioForIntersection = 0.7;
	private final double OCR_TO_BOND_MAX_DISTANCE=3.0;
	private final double maxCandidateRatioForIntersection = 1.5;
	private final double maxCandidateRatioForIntersectionWithNeighbor = 1.3;       
	private final double MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS_THIN = 1;
	private final double MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS_FULL = 0.5;
	private final double MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS = 6;
	private final double MIN_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS = 2;
	private final double MAX_BOND_TO_AVG_BOND_RATIO_FOR_INTERSECTION= 0.8;

	private final double MAX_ANGLE_FOR_JOINING_SEGMENTS=25 * Math.PI/180.0;
	private final double MIN_SIZE_FOR_ANGLE_COMPARE_JOINING_SEGMENTS=8.0;

	private final double MAX_DISTANCE_TO_MERGE_PARALLEL_LINES=2;
	private final double MAX_POINT_DISTANCE_TO_BE_PART_OF_MULTI_NODE= 1;

	private final double MAX_BOND_RATIO_FOR_OCR_CHAR_SPACING=0.3;
	private final double MAX_THETA_FOR_OCR_SEPERATION=80 * Math.PI/180.0;


	//This number is likely one of the most important to adjust.
	//It may have to have some changes done to the algorithm using it too
	private final double MAX_BOND_RATIO_FOR_MERGING_TO_OCR=0.28;
	
	private final double MAX_BOND_RATIO_FOR_LINES_CONSIDERED_FOR_POSITIONING_OCR=MAX_BOND_RATIO_FOR_MERGING_TO_OCR;


	private final double MIN_LONGEST_WIDTH_RATIO_FOR_OCR_TO_AVERAGE=0.5;
	private final double MIN_AREA_RATIO_FOR_OCR_TO_AVERAGE=0.6;
	private final double MAX_AREA_RATIO_FOR_OCR_TO_AVERAGE=2.5;
	private final double MIN_AREA_RATIO_FOR_HULL_TO_BBOX_OCR=0.5;



	//For finding high order bonds
	private final double MIN_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS=.5;
	private final double MIN_BIGGER_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS=.25;
	private final double MAX_ANGLE_FOR_PARALLEL=10.0 * Math.PI/180.0;

	private final double MIN_ST_DEV_FOR_KEEPING_DASHED_LINES=0.07;



	//Parallel lines

	private final double MAX_DELTA_LENGTH_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC=7.0;
	private final double MAX_ANGLE_FOR_PARALLEL_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC=15.0 * Math.PI/180.0;


	
	
	//This is a newish feature, and it slows things down. Turn off if speed is needed.
	private final boolean PRE_RESCUE_OCR = true;

	/**
	 * Create a new extractor from the given bufferedImage.
	 * @param bufferedImage
	 * @return
	 * @throws IOException
	 */
	public static StructureImageExtractor createFromImage(BufferedImage bufferedImage)throws Exception{
		BufferedImage img = bufferedImage;
		if(BufferedImage.TYPE_BYTE_GRAY != bufferedImage.getType()){
			img = toGrayScale(bufferedImage);
		}
		return new StructureImageExtractor(img.getRaster());

	}

	/**
	 * Convert the given BufferdImage into a new GrayScaled image.
	 * the input image is not modified.
	 * @param image the image to convert;
	 * @return a new BufferedImage of the same size as the input
	 * but in grayscale.
	 */
	private static BufferedImage toGrayScale(BufferedImage image){
		BufferedImage grayScaled = new BufferedImage(image.getWidth(),image.getHeight(),	BufferedImage.TYPE_BYTE_GRAY);

		ColorConvertOp op = new ColorConvertOp(
				image.getColorModel().getColorSpace(),
				grayScaled.getColorModel().getColorSpace(),null);
		op.filter(image,grayScaled);
		return grayScaled;

	}
	public StructureImageExtractor(Raster raster)throws Exception{
		this(raster, false);
	}
	public StructureImageExtractor(Raster raster, boolean debug )throws Exception{
		this.DEBUG = debug;
		load(Bitmap.createBitmap(raster));
	}
	public StructureImageExtractor(byte[] file, boolean debug) throws Exception{
		this.DEBUG=debug;
		load(file);
	}
	public StructureImageExtractor(File file, boolean debug) throws Exception{
		this.DEBUG=debug;
		load(file);
	}
	
	public StructureImageExtractor(byte[] file) throws Exception{
		this(file,false);
	}
	public StructureImageExtractor(File file) throws Exception{
		this(file,false);
	}
	
	
	
	
	
	

	private static Tuple<Character,Number> adjustConfidence(Tuple<Character,Number> tup){
		String ch=tup.k()+"";
		double invScore=1-tup.v().doubleValue();
		if(ch.equals("K") || ch.equals("k") || ch.equals("f")){
			invScore=invScore*3.5; // penalize "K"
		}
		if(ch.equals("R")||
				//ch.equalsIgnoreCase("A")||
				//ch.equalsIgnoreCase("Z")||
				ch.equalsIgnoreCase("-")||
				ch.equals("m")||
				ch.equalsIgnoreCase("W")||

				ch.equals("n")){
			invScore=invScore*3; // penalize
		}
		if(ch.equalsIgnoreCase("X") || ch.equalsIgnoreCase("+") || ch.equalsIgnoreCase("D") || ch.equalsIgnoreCase("Z")){
			invScore=invScore*1.5; // penalize
		}else if(ch.equals("N") 
				|| ch.equalsIgnoreCase("C") 
				|| ch.equalsIgnoreCase("O")				
				){
			invScore=invScore*(0.87); // promote
		}
		
		if(ch.equals("H")){
			invScore=invScore*(0.92); // promote
		}
//		
		if(ch.equals("h")){
			invScore=invScore*(1.04); // penalize
		}
		

		return Tuple.of(tup.k(),Math.max(0,1-invScore));

	}

	public static enum CharType{
		ChemLikely,
		NumericLikely,
		BondLikely,
		VerticalBondLikely,
		
	}
	
	private static CharType OCRIsLikely(Tuple<Character,Number> tup){
		String t=tup.k()+"";
		if(     "I".equalsIgnoreCase(t) || 
				"L".equalsIgnoreCase(t) ||
				"1".equalsIgnoreCase(t) ||
				"t".equalsIgnoreCase(t) ||
				"(".equalsIgnoreCase(t) ||
				")".equalsIgnoreCase(t) ||
				"f".equals(t)){
			return CharType.VerticalBondLikely;
		}else if(
				"-".equalsIgnoreCase(t) ||
				"/".equalsIgnoreCase(t) ||
				"K".equalsIgnoreCase(t) ||
				"Y".equalsIgnoreCase(t) ||
				"W".equalsIgnoreCase(t) ||
				"\\".equalsIgnoreCase(t)
				){
			return CharType.BondLikely;
			
		}else if(
				"2".equalsIgnoreCase(t) ||
				"3".equalsIgnoreCase(t) ||
				"4".equalsIgnoreCase(t) ||
				"5".equalsIgnoreCase(t) ||
				"6".equalsIgnoreCase(t) ||
				"7".equalsIgnoreCase(t) ||
				"9".equalsIgnoreCase(t)){
			return CharType.NumericLikely;
		}
		return CharType.ChemLikely;
	}



	private void processOCR(SCOCR socr, List<Shape> polygons,Bitmap bitmap, Bitmap thin, BiConsumer<Shape,List<Tuple<Character,Number>>> onFind) throws Exception{
		
		boolean[] interupt=new boolean[]{false};
		/*
		 * Looks at each polygon, and gets the likely OCR chars.
		 */   
		List<Shape> toAddShapes = Collections.synchronizedList(new ArrayList<Shape>());
		List<Shape> toRemoveShapes = Collections.synchronizedList(new ArrayList<Shape>());

		polygons
				.parallelStream()
//				.stream()
//		.collect(Collectors.toList())
//		.stream()
		
//		.parallel()
		.forEach(s->{
			if (Thread.interrupted())  // Clears interrupted status!
				interupt[0]=true;
			if(interupt[0]){
				return;
			}
			
			Rectangle2D bounds2d = s.getBounds2D();
			if(bounds2d.getWidth()>0 && bounds2d.getHeight()>0){
				List<Tuple<Character,Number>> ll = new ArrayList<>();
				
				processOCRShape(socr,s,bitmap,thin,(sf,lf)->{
					ll.addAll(lf);
					//onFind.accept(sf, lf);
				});
				
				double bestMatch1 = ll.stream().findFirst().map(t->t.v().doubleValue()).orElse(0.0);
				double[] bestMatch = new double[]{bestMatch1,bestMatch1};
				
				
				// if the width is too wide, it might be two chars pushed together
				// but that's only really likely if the shape is pretty close to being a box
				if(bounds2d.getWidth() >  bounds2d.getHeight()){
					double sarea=GeomUtil.area(s);
					double bbarea=GeomUtil.area(bounds2d);
					
					
					
					
					if(sarea/bbarea > 0.8 && bounds2d.getHeight()>4){
						
						double[] ratios=new double[]{0.5};
						
						
						if(bounds2d.getWidth()>2.3*bounds2d.getHeight() &&bounds2d.getWidth()<3.4*bounds2d.getHeight()){
							ratios=new double[]{0.58,0.44};
						}
						
						boolean better=false;
						Shape gshape1=null;
						Shape gshape2=null;
						
						List<Tuple<Shape,List<Tuple<Character,Number>>>> bsplitMatches=null;
						
						for(double ratio:ratios){
							List<Tuple<Shape,List<Tuple<Character,Number>>>> splitMatches = new ArrayList<>();
							
							Shape box1= new Rectangle2D.Double(bounds2d.getMinX(), bounds2d.getMinY(), bounds2d.getWidth()*ratio, bounds2d.getHeight());
							Shape box2= new Rectangle2D.Double(bounds2d.getMinX() + bounds2d.getWidth()*ratio, bounds2d.getMinY(), bounds2d.getWidth()*(1-ratio), bounds2d.getHeight());
							
							Shape cropShape1=GeomUtil.getIntersectionShape(box1, s).get();
							Shape cropShape2=GeomUtil.getIntersectionShape(box2, s).get();
							
							processOCRShape(socr,cropShape1,bitmap,thin,(sf,lf)->{
								if(lf.get(0).v().doubleValue()>=bestMatch[0]-0.0){
									splitMatches.add(Tuple.of(sf,lf));
								}
								//onFind.accept(sf, lf);
							});
							processOCRShape(socr,cropShape2,bitmap,thin,(sf,lf)->{
								if(lf.get(0).v().doubleValue()>=bestMatch[1]-0.0){
									splitMatches.add(Tuple.of(sf,lf));
								}
							});
							if(splitMatches.size()>1){
								bsplitMatches=splitMatches;
								gshape1=cropShape1;
								gshape2=cropShape2;
								
								bestMatch[0]=splitMatches.get(0).v().get(0).v().doubleValue();
								bestMatch[1]=splitMatches.get(1).v().get(0).v().doubleValue();
								
								better=true;
							}
						}
						
						if(better){
							bsplitMatches.stream().forEach(t->{
								onFind.accept(t.k(), t.v());
							});
							toAddShapes.add(gshape1);
							toAddShapes.add(gshape2);
							toRemoveShapes.add(s);
						}else{
							onFind.accept(s, ll);
						}
						
					}else{
						onFind.accept(s, ll);
					}
					// if the height is too tall, it might be two chars pushed together
					// but that's only really likely if the shape is pretty close to being a box
				}else if(bounds2d.getHeight() >  bounds2d.getWidth()*2.4 && bounds2d.getHeight() <  bounds2d.getWidth()*3){
						double sarea=GeomUtil.area(s);
						double bbarea=GeomUtil.area(bounds2d);
						
						if(sarea/bbarea > 0.8 && bounds2d.getHeight()>4 && bounds2d.getWidth()>8){
							
							double[] ratios=new double[]{0.5};
							
							boolean better=false;
							Shape gshape1=null;
							Shape gshape2=null;
							
							List<Tuple<Shape,List<Tuple<Character,Number>>>> bsplitMatches=null;
							
							for(double ratio:ratios){
								List<Tuple<Shape,List<Tuple<Character,Number>>>> splitMatches = new ArrayList<>();
								
								Shape box1= new Rectangle2D.Double(bounds2d.getMinX(), bounds2d.getMinY(), bounds2d.getWidth(), bounds2d.getHeight()*ratio);
								Shape box2= new Rectangle2D.Double(bounds2d.getMinX(), bounds2d.getMinY() + bounds2d.getHeight()*ratio, bounds2d.getWidth(), bounds2d.getHeight()*(1-ratio));
								
								Shape cropShape1=GeomUtil.getIntersectionShape(box1, s).get();
								Shape cropShape2=GeomUtil.getIntersectionShape(box2, s).get();
								
								processOCRShape(socr,cropShape1,bitmap,thin,(sf,lf)->{
									if(lf.get(0).v().doubleValue()>=bestMatch[0]-0.0){
										splitMatches.add(Tuple.of(sf,lf));
									}
									//onFind.accept(sf, lf);
								});
								processOCRShape(socr,cropShape2,bitmap,thin,(sf,lf)->{
									if(lf.get(0).v().doubleValue()>=bestMatch[1]-0.0){
										splitMatches.add(Tuple.of(sf,lf));
									}
								});
								if(splitMatches.size()>1){
									bsplitMatches=splitMatches;
									gshape1=cropShape1;
									gshape2=cropShape2;
									
									bestMatch[0]=splitMatches.get(0).v().get(0).v().doubleValue();
									bestMatch[1]=splitMatches.get(1).v().get(0).v().doubleValue();
									
									better=true;
								}
							}
							
							if(better){
								bsplitMatches.forEach(t->{
									onFind.accept(t.k(), t.v());
								});
								toAddShapes.add(gshape1);
								toAddShapes.add(gshape2);
								toRemoveShapes.add(s);
							}else{
								onFind.accept(s, ll);
							}
							
						}else{
							onFind.accept(s, ll);
						}
					}else{
					onFind.accept(s, ll);
				}				
			}
		});
		
		if(interupt[0])throw new InterruptedException();
		
		polygons.removeAll(toRemoveShapes);
		polygons.addAll(toAddShapes);
	}

	private void processOCRShape(SCOCR socr, Shape s, Bitmap bitmap, Bitmap thin,BiConsumer<Shape,List<Tuple<Character,Number>>> onFind){
		//this is a critical section that is called thousands of times
		//so some work as been put in to optimize it
		//shaving off a few ms really adds up!

		if(s.getBounds2D().getWidth()>0 && s.getBounds2D().getHeight()>0){
			
			
			double areareal=GeomUtil.area(s);
			//we compute this a couple of times in the if statements below so cache it.
			CachedSupplier<Double> areaRealDivByAreaBox = CachedSupplier.of(() ->areareal/ GeomUtil.area(s.getBounds2D()));

			if(areareal<=5)return;
			//we only really care about the "best" character and never update that value
			//even after filtering...
			//the rest of the time char lookups just look for contains without worrying about order
			Character[] best =new Character[1]; //this is done to set it in a lambda
			char[] asciiCache = new char[128]; // we only check against ASCII values

			Bitmap cropped = bitmap.getLazyCrop(s);
			
			if(cropped.fractionPixelsOn()>0.9){
				if(s.getBounds2D().getWidth()>s.getBounds2D().getHeight()){
					onFind.accept(s, Stream.of('-','-','-','-').map(c->Tuple.of(c, (Number)0.8)).collect(Collectors.toList()));
				}else{
					onFind.accept(s, Stream.of('I','l','t','i').map(c->Tuple.of(c, (Number)0.8)).collect(Collectors.toList()));
				}
				return;
				
			}
			
			List<Tuple<Character,Number>> potential = socr.getNBestMatches(4,
					cropped,
					thin.getLazyCrop(s))
					.stream()
					.map(Tuple::of)
					.map(t->adjustConfidence(t))
					.map(t->t.withVComparator())
					.sorted(Comparator.reverseOrder())
					.peek(t->{
						if(best[0] ==null){
							best[0] = t.k();
						}
						char c = t.k();
						if(c < 128){
							asciiCache[c]=1;
						}
					})
					.collect(Collectors.toList());

			if(asciiCache['N']==1 || asciiCache['S']==1|| asciiCache['s']==1){


				//this usually means it's not a real "N" or S
				boolean alreadyFiltered=false;
				if(areaRealDivByAreaBox.get() <0.5){
					if(asciiCache['\\']==1 || asciiCache['X']==1 || asciiCache['K']==1 ||
							asciiCache['k']==1 || asciiCache['-']==1){
						potential = potential.stream()
								.filter(t->!t.k().equals('N') && !t.k().equals('S') && !t.k().equals('s'))
								.collect(Collectors.toList());
						alreadyFiltered=true;
					}
				}
				if(!alreadyFiltered){
					Rectangle2D rbox = s.getBounds2D();
					//probably not an N or S
					if(rbox.getWidth()>rbox.getHeight()*1.3){
						potential = potential.stream()
								.filter(t->!t.k().equals('N') && !t.k().equals('S') && !t.k().equals('s'))
								.collect(Collectors.toList());
					}
				}
				
				//This is the least justified tweak, just a strange thing about N+ that I've noticed
			}else if(asciiCache['M']==1 && asciiCache['m']==1 && (asciiCache['P']==1 || (asciiCache['K']==1 && asciiCache['-']==1))){

				boolean goodEnough=potential.stream().filter(t->t.v().doubleValue()>0.4).findAny().isPresent();
				boolean badEnough=potential.stream().filter(t->t.v().doubleValue()>0.3).findAny().isPresent();
				
				if(!goodEnough && badEnough){
					if(areaRealDivByAreaBox.get() >0.5){
						potential.add(0,Tuple.of('?',0.7));
					}
				}
			}else if((asciiCache['9']==1 && asciiCache['p']==1 && asciiCache['b']==1 && asciiCache['3']==1 ) ||
					 (asciiCache['9']==1 && asciiCache['I']==1 && asciiCache['i']==1 && asciiCache['3']==1 )
					 ){

				boolean goodEnough=potential.stream().filter(t->t.v().doubleValue()>0.4).findAny().isPresent();
				boolean badEnough=potential.stream().filter(t->t.v().doubleValue()>0.3).findAny().isPresent();
				
				if(!goodEnough && badEnough){
					if(areaRealDivByAreaBox.get() >0.5){
						potential.add(0,Tuple.of(')',0.7));
					}
				}
			}else if(asciiCache['C']==1 && asciiCache['c']==1 && asciiCache['6']==1 && (asciiCache['0']==1 )){

				boolean goodEnough=potential.stream().filter(t->t.v().doubleValue()>0.45).findAny().isPresent();
				boolean badEnough=potential.stream().filter(t->t.v().doubleValue()>0.3).findAny().isPresent();
				
				if(!goodEnough && badEnough){
					if(areaRealDivByAreaBox.get() >0.5){
						potential.add(0,Tuple.of('(',0.7));
					}
				}
			}

			if(Character.valueOf('L').equals(best[0])){
				Rectangle2D rbox = s.getBounds2D();
				if(rbox.getWidth()>rbox.getHeight()*0.6){
					//Too wide for an L, but since it's the highest confidence, it's
					//probably just part of a bond system. 
					potential = potential.stream()
							.map(Tuple.vmap(n->(Number)Double.valueOf(0D)))
							.collect(Collectors.toList());
				}
			}
			if(asciiCache['K']==1  && asciiCache['X']==1){
				if(areaRealDivByAreaBox.get() <0.5){

					potential = potential.stream()
							.map(Tuple.vmap(n->(Number)0.0))
							.collect(Collectors.toList());

				}
			}

			if(asciiCache['S']==1  && asciiCache['s']==1 && asciiCache['8']==1){
				//It's probably an S in this case, slightly adjust numbers for those
				potential = potential.stream()
						.map(t->{
							if(t.k().equals('S') || t.k().equals('s')){
								return Tuple.of(t.k(),(Number)Math.max(0, (1-(1-t.v().doubleValue())*0.8)));
							}
							return t;
						})
						.map(t->t.withVComparator())
						.sorted(Comparator.reverseOrder())
						.collect(Collectors.toList());
			}

			if(asciiCache['D']==1 && (asciiCache['U']==1 ||asciiCache['u']==1)){

				if(best[0] != null && (best.equals('D') || best.equals('U') || best.equals('u'))){
					//It's probably an O, just got flagged wrong
					potential = potential.stream()
							.map(Tuple.kmap(c->'O'))
							.collect(Collectors.toList());
				}
			}

			onFind.accept(s, potential);
		}
	}

	private void load(byte[] file) throws Exception{
		load(bitmap = Bitmap.read(file).clean());

	}
	private void load(File file) throws Exception{
		load(bitmap = Bitmap.read(file).clean());

	}
	
	private boolean DEBUG=false;
	
	
	private void rescueOCR(List<LineWrapper> lines, List<Shape> polygons, Set<Shape> likelyOCRAll, SCOCR scocr, BiConsumer<Shape,List<Tuple<Character,Number>>> cons){
		
			double averageWidthOCRFinal=likelyOCRAll.stream()
					.map(Shape::getBounds2D)
					.filter(Objects::nonNull)
					.mapToDouble(Rectangle2D::getWidth)
					.filter(Objects::nonNull) //sometimes get NPEs here not sure why or how
					.average()
					.orElse(0);
			double averageHeightOCRFinal=likelyOCRAll.stream()
					.map(Shape::getBounds2D)
					.filter(Objects::nonNull)
					.mapToDouble(Rectangle2D::getHeight)
					.filter(Objects::nonNull) //sometimes get NPEs here not sure why or how
					.average()
					.orElse(0);
			
			
		
			//this is a spotty attempt at final rescue, not particularly good
			//only do this sometimes
			List<Tuple<Point2D,Line2D>> centerPoints = lines.stream()
															.map(l->Tuple.of(l.centerPoint(),l.getLine()))
															.collect(Collectors.toList());
			
			double ldensityCutoff=likelyOCRAll.stream()
			        .filter(p->p.getBounds2D().getWidth()>averageWidthOCRFinal*0.8)
			        .filter(p->p.getBounds2D().getHeight()>averageHeightOCRFinal*0.8)
			        .map(p->Tuple.of(p,GeomUtil.centerOfMass(p)))
			        .map(t->t.k())
			        .map(s->Tuple.of(Tuple.of(s,GeomUtil.growShape(s, 2)), centerPoints))
			        .map(t->Tuple.of(t.k().k(),t.v().stream()
			        		     .filter(p->t.k().v().contains(p.k()))
			        		     .map(p->p.v())
			        		     .collect(Collectors.toList())
			        		))
			        .map(t->{
			        	BoundingBox ocrBounds=BoundingBox.of(t.k().getBounds2D(), t.v());
			        	return ocrBounds;
			        })
			        .mapToDouble(o->o.getLineDensity())
			        .average()
			        .orElse(0);
			
			List<BoundingBox> rescueShapes1 = new ArrayList<>();
			
			polygons.stream()
			        .filter(p->p.getBounds2D().getWidth()>averageWidthOCRFinal*0.8)
			        .filter(p->p.getBounds2D().getHeight()>averageHeightOCRFinal*0.8)
			        .filter(p->GeomUtil.area(p.getBounds2D())>averageHeightOCRFinal*averageWidthOCRFinal)
			        .map(p->Tuple.of(p,GeomUtil.centerOfMass(p)))
			        .filter(t-> !likelyOCRAll.stream().filter(oc->oc.contains(t.v())).findFirst().isPresent())
			        .map(t->t.k())
			        .map(s->Tuple.of(GeomUtil.growShape(s, 2), centerPoints))
			        .map(t->Tuple.of(t.k(),t.v().stream()
			        		     .filter(p->t.k().contains(p.k()))
			        		     .map(p->p.v())
			        		     .collect(Collectors.toList())
			        		))
			        .map(t->t.v())
			        
			        .flatMap(t->{
			        	List<Tuple<Line2D,Point2D>> mlist= t.stream().map(t1->Tuple.of(t1, GeomUtil.findCenterOfShape(t1))).collect(Collectors.toList());
			        	
			        	
			        	
			        	return GeomUtil.groupThings(mlist, t1->{
					        		return t1.k().v().distance(t1.v().v())<averageHeightOCRFinal*1.2;
					        	})
					        	.stream()
					        	.filter(lpts->GeomUtil.area(lpts.stream()
					        			                        .flatMap(lt->Stream.of(lt.k().getP1(),lt.k().getP2()))
					        			                        .collect(GeomUtil.convexHull()))>averageHeightOCRFinal*averageWidthOCRFinal*0.4
					        			);
			        	
					        	
			        })
			        .map(t->t.stream().map(t1->t1.k()).collect(Collectors.toList()))
			        .filter(t->t.size()>3)
			        .forEach(t->{
			        	List<BoundingBox> bblist=GeomUtil.getBoundingBoxesContaining(averageWidthOCRFinal*1.1, averageHeightOCRFinal*1.1, t, averageWidthOCRFinal);
			        	bblist.stream()
			        		.map(bb->bb.resize(averageWidthOCRFinal,averageHeightOCRFinal))
			        		.filter(bb->bb.getLineDensity()>ldensityCutoff*0.9)
				        	.forEach(b->{
				        		Point2D center = b.getCenterOfMass();
				        		if(!likelyOCRAll.stream().filter(oc->oc.contains(center)).findFirst().isPresent()){
				        			
				        			rescueShapes1.add(b);
				        		}
				        	});
			        });
			
			rescueShapes1.forEach(bb->{

				boolean[] got = new boolean[]{false};

				if(!got[0]){
					//I'm not sure 1 is the right number here
					Shape nshape = GeomUtil.growShape(bb.getCenteredConvexRect(),1);
					
					
					
					
					processOCRShape(scocr,nshape,bitmap,thin,(s,potential)->{
						String st=potential.get(0).k().toString();
						double cos=potential.get(0).v().doubleValue();
						if(cos>OCRcutoffCosineRescueInitial){
							BranchNode n=BranchNode.interpretOCRStringAsAtom2(st);
							if(n!=null && n.isRealNode()){
								if(cos>OCRcutoffCosine || bb.getLineDensity()>ldensityCutoff*1.05){
									cons.accept(nshape, potential);
									got[0]=true;
								}
							}
						}	
					});
				}
				if(!got[0]){
					Shape nshape =  GeomUtil.growShape(bb.getRect(),1);
					processOCRShape(scocr,nshape,bitmap,thin,(s,potential)->{
						
						String st=potential.get(0).k().toString();
						double cos=potential.get(0).v().doubleValue();
						
						if(cos>OCRcutoffCosineRescueInitial){
							BranchNode n=BranchNode.interpretOCRStringAsAtom2(st);
							if(n!=null && n.isRealNode()){
								if(cos>OCRcutoffCosine || bb.getLineDensity()>ldensityCutoff*1.05){
									cons.accept(nshape, potential);
								}
							}
						}	
					});
				}
			});
	}

	private void load(Bitmap aBitMap) throws Exception{

		
		
		
		ctabRaw.clear();
		ocrAttempt.clear();
		bitmap = aBitMap;

		SCOCR[] socr=new SCOCR[]{OCR_DEFAULT.orElse(OCR_BACKUP, OCRcutoffCosine)};

		double[] maxBondLength=new double[]{INITIAL_MAX_BOND_LENGTH};    




		
		

		thin = bitmap.thin();
		
		{
			List<int[]> hollow =thin.findHollowPoints();
			
			if(hollow.size()> 0.002*thin.fractionPixelsOn()*thin.width()*thin.height()){
				bitmap=new Bitmap.BitmapBuilder(bitmap).boxBlur(1).threshold(2).build();
				thin=bitmap.thin();
			}
			
			
		}
		
		polygons = bitmap.connectedComponents(Bitmap.Bbox.DoublePolygon);

		boolean isLarge = false;
		if (!polygons.isEmpty()) {
			isLarge = polygons.size() > 4000;
			
			//it could be that there's just a lot of really tiny pieces
//			
//			if(isLarge){
//				System.out.println("Looks large");
//				List<Shape> toRemove=polygons.stream()
//				        					 .map(p->Tuple.of(p,GeomUtil.area(p)))
//				        					 .filter(t->t.v()<5)
//				        					 .map(t->t.k())
//				        					 .map(s->GeomUtil.growShapeNPoly(s, 1, 16))
//				        					 .collect(Collectors.toList());
//				if(!toRemove.isEmpty() && (polygons.size()-toRemove.size())<4000){
//					System.out.println("removing noise");
//					Bitmap bm=new Bitmap.BitmapBuilder(bitmap)
//					          .remove(toRemove)
//					          .build();
//					load(bm);
//					return;
//					
//				}
//			}
			
		}
		


		// segments are generated for thinned bitmap only, since
		//  it can quite noisy on normal bitmap!
		if (isLarge) {
			throw new IllegalStateException("Cannot support images with over 4000 polygons at this time");
		}

		Set<Shape> likelyOCR= Collections.synchronizedSet(new LinkedHashSet<Shape>());
		Set<Shape> likelyOCRNumbers= Collections.synchronizedSet(new LinkedHashSet<Shape>());
		Set<Shape> likelyOCRNonBond= Collections.synchronizedSet(new LinkedHashSet<Shape>());
		Set<Shape> likelyOCRAll=Collections.synchronizedSet(new LinkedHashSet<Shape>());
		Set<Shape> likelyOCRIgnore = Collections.synchronizedSet(new LinkedHashSet<Shape>());
		
		List<Shape> ocrRescues = new ArrayList<Shape>();
		
		/*
		 * Looks at each polygon, and gets the likely OCR chars.
		 */   

		Set<Shape> verticalShapes = Collections.synchronizedSet(new LinkedHashSet<Shape>());

		Set<Shape> NPlusShapes = Collections.synchronizedSet(new LinkedHashSet<Shape>());
		
		processOCR(socr[0],polygons,bitmap,thin,(s,potential)->{
			ocrAttempt.put(s, potential);
			
			if(potential.stream().filter(e->e.v().doubleValue()>OCRcutoffCosine).findAny().isPresent()){
				CharType ct=OCRIsLikely(potential.get(0));
				
				if(potential.get(0).k().toString().equals("?")){
					NPlusShapes.add(s);
				}else{
					if(ct.equals(CharType.ChemLikely)){
						likelyOCR.add(s);
						likelyOCRNonBond.add(s);
						
					}else if(ct.equals(CharType.NumericLikely)){
						likelyOCRNonBond.add(s);
						likelyOCRNumbers.add(s);
					}else if(ct.equals(CharType.VerticalBondLikely)){
						verticalShapes.add(s);
					}
					likelyOCRAll.add(s);
				}
			}
			
		});
		
		
		double averageHeightOCR1=likelyOCR.stream()
				.map(Shape::getBounds2D)
				.filter(Objects::nonNull)
				.mapToDouble(Rectangle2D::getHeight)
				.filter(Objects::nonNull) //sometimes get NPEs here not sure why or how
				.average()
				.orElse(10000);
		
		double averageWidthOCR1=likelyOCR.stream()
				.map(Shape::getBounds2D)
				.filter(Objects::nonNull)
				.mapToDouble(Rectangle2D::getWidth)
				.filter(Objects::nonNull) //sometimes get NPEs here not sure why or how
				.average()
				.orElse(0);
		

		NPlusShapes
			.stream()
			.filter(s->s.getBounds2D().getWidth()<averageWidthOCR1*2)
			.forEach(s->{
				
				
				
				CharType ct=OCRIsLikely(ocrAttempt.get(s).get(0));
				
				if(ct.equals(CharType.ChemLikely)){
					likelyOCR.add(s);
					likelyOCRNonBond.add(s);
					
				}else if(ct.equals(CharType.NumericLikely)){
					likelyOCRNonBond.add(s);
					likelyOCRNumbers.add(s);
				}else if(ct.equals(CharType.VerticalBondLikely)){
					verticalShapes.add(s);
				}
				likelyOCRAll.add(s);				
				
			});
		
		
		
		verticalShapes.stream()
		              .filter(s->s.getBounds2D().getHeight()>averageHeightOCR1*1.5)
		              .forEach(s->{
		            	  likelyOCR.remove(s);
		            	  likelyOCRNumbers.remove(s);
		            	  likelyOCRAll.remove(s);
		              });
		
		

		List<Shape> circles =  polygons.stream()
		        .filter(p->!likelyOCRAll.contains(p))
		        .map(s->Tuple.of(s,GeomUtil.getCircleLikeScore(s)))
		        .filter(t->t.v()>0.9)
		        .map(t->t.k())
		        .map(s->GeomUtil.growShape(s, 2))
		        .collect(Collectors.toList());
		
		
		lines= GeomUtil.asLines(thin.segments())
					   .stream()
					   .filter(l->!circles.stream()
								           .filter(s->s.contains(l.getP1()) || s.contains(l.getP1()))
								           .findFirst()
								           .isPresent())
					   .map(l->GeomUtil.LineWrapper.of(l))
					   .collect(Collectors.toList());
		
		
		
		
		
		
		
		
		ctabRaw.clear();

		boolean[] foundNewOCR=new boolean[]{true};
	
		int repeats=0;
		List<Shape> realRescueOCRCandidates = Collections.synchronizedList(new ArrayList<>());
		
		double[] averageHeightOCRFinal = new double[]{0};
		double[] averageWidthOCRFinal = new double[]{0};
		
		List<Point2D> intersectionNodes = new ArrayList<>();
		
		
		if(PRE_RESCUE_OCR){
			rescueOCR(lines,polygons,likelyOCR,socr[0],(s,potential)->{
				String ss = potential.get(0).k().toString();
				if(ss.equalsIgnoreCase("C") || ss.equalsIgnoreCase("P") || ss.equalsIgnoreCase("S")){
					return;
				}
////				//never get B, probably an H
//				if(ss.equals("H")){
//					System.out.println("It's a b");
////					potential.set(0, Tuple.of('H',potential.get(0).v()));
//					RasterChar rc=RasterChar.fromDefault(bitmap.crop(s)).blur(2);
//	            	System.out.println(Base64.getEncoder().encodeToString(rc.rawDataAsString().getBytes()));
//				}				
				
				Point2D cent=GeomUtil.findCenterOfShape(s);
				boolean cont=likelyOCRAll.stream()
				            .filter(s1->s1.contains(cent))
				            .findAny()
				            .isPresent();
				if(cont)return;
				
				ocrAttempt.put(s, potential);
				
				CharType ct=OCRIsLikely(potential.get(0));
				if(ct.equals(CharType.ChemLikely)){
					likelyOCR.add(s);
					likelyOCRNonBond.add(s);
					
				}else if(ct.equals(CharType.NumericLikely)){
					likelyOCRNonBond.add(s);
					likelyOCRNumbers.add(s);
				}
				likelyOCRAll.add(s);
				
				//realRescueOCRCandidates.add(s);
			});
		}

		double[] ignoreTooSmall=new double[]{0.0};
		
		while(foundNewOCR[0] && repeats<MAX_OCR_FULL_REPEATS){
			
			if (Thread.currentThread().isInterrupted()){
			      throw new InterruptedException();
			}
			
			repeats++;
			foundNewOCR[0]=false;
			intersectionNodes.clear();

			double averageLargestOCR=likelyOCR.stream()
					.map(s->GeomUtil.getPairOfFarthestPoints(s))
					.filter(p -> p !=null && p.length ==2)
					.mapToDouble(p->p[0].distance(p[1]))
					.average()
					.orElse(0);
			double averageAreaOCR=likelyOCR.stream()
					.mapToDouble(s->GeomUtil.area(s))
					.average()
					.orElse(0);

			double averageWidthOCR=likelyOCR.stream()
					.map(Shape::getBounds2D)
					.filter(Objects::nonNull)
					.mapToDouble(Rectangle2D::getWidth)
					.filter(Objects::nonNull) //sometimes get NPEs here not sure why or how
					.average()
					.orElse(0);
			double averageHeightOCR=likelyOCR.stream()
					.map(Shape::getBounds2D)
					.filter(Objects::nonNull)
					.mapToDouble(Rectangle2D::getHeight)
					.filter(Objects::nonNull) //sometimes get NPEs here not sure why or how
					.average()
					.orElse(0);
			

			
			
			
			
			
			
			
			double averageWidthNumberOCR=likelyOCRNumbers.stream()
					.map(Shape::getBounds2D)
					.filter(Objects::nonNull)
					.mapToDouble(Rectangle2D::getWidth)
					.filter(Objects::nonNull) //sometimes get NPEs here not sure why or how
					.average()
					.orElse(0);
			
			
			double averageHeightNumberOCR=likelyOCRNumbers.stream()
					.map(Shape::getBounds2D)
					.filter(Objects::nonNull)
					.mapToDouble(Rectangle2D::getHeight)
					.filter(Objects::nonNull) //sometimes get NPEs here not sure why or how
					.average()
					.orElse(0);
			//likelyOCRNumbers
			averageHeightOCRFinal[0] = averageHeightOCR;
			averageWidthOCRFinal[0] = averageWidthOCR;


			likelyOCRAll.retainAll(likelyOCRAll.stream()
					.filter(Objects::nonNull)
					.map(s->Tuple.of(s,GeomUtil.getPairOfFarthestPoints(s)))
					.filter(t->t.v()[0].distance(t.v()[1]) > averageLargestOCR*MIN_LONGEST_WIDTH_RATIO_FOR_OCR_TO_AVERAGE)
					.map(t->t.k())
					.collect(Collectors.toList()));

			likelyOCR.retainAll(likelyOCR.stream()
					.map(s->Tuple.of(s,GeomUtil.getPairOfFarthestPoints(s)))
					.filter(t->t.v()[0].distance(t.v()[1]) > averageLargestOCR*MIN_LONGEST_WIDTH_RATIO_FOR_OCR_TO_AVERAGE)
					.map(t->t.k())
					.collect(Collectors.toList()));



			Predicate<Line2D> isInOCRShape = (l)->{
				if(likelyOCR.isEmpty())return false;
				Optional<Tuple<Shape,Double>> shape1=GeomUtil.findClosestShapeTo(likelyOCRNonBond, l.getP1());
				
				if(!shape1.isPresent() || shape1.get().v()>OCR_TO_BOND_MAX_DISTANCE){
					return false;
				}
				Optional<Tuple<Shape,Double>> shape2=GeomUtil.findClosestShapeTo(likelyOCRNonBond, l.getP2());
				if(!shape2.isPresent() || shape2.get().v()>OCR_TO_BOND_MAX_DISTANCE){
					return false;
				}
				if(shape1.get().k().equals(shape2.get().k())){
					return true;
				}
				
				boolean anyOutside=GeomUtil.getLinesNotInside(l, Arrays.asList(shape1.get().k(),shape2.get().k()))
				        .stream()
				        .filter(Objects::nonNull)
				        .filter(GeomUtil.longerThan(1))
				        .findAny()
				        .isPresent();
				if(!anyOutside)return true;
				
				return false;
			};

			Predicate<Line2D> tryToMerge = isInOCRShape.negate().and((l)->{
				return true;

				//return LineUtil.length(l)<largestBond;
			});

			List<Line2D> useLines = lines.stream()
										 .filter(t->!likelyOCRIgnore.stream()
									    		    .filter(s->s.contains(t.centerPoint()))
									    		    .findAny().isPresent())
										 .map(l->l.getLine())
									     .collect(Collectors.toList());
			


			List<LineWrapper> smallLines=useLines.stream()
					.filter(tryToMerge)
					.map(l->LineWrapper.of(l))
					.collect(Collectors.toList());

			List<Line2D> bigLines=useLines.stream()
					.filter(tryToMerge.negate())
					.collect(Collectors.toList());


			// smallLines= thin.combineLines(smallLines, MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS, MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS_THIN, MAX_POINT_DISTANCE_TO_BE_PART_OF_MULTI_NODE,MAX_ANGLE_FOR_JOINING_SEGMENTS,MIN_SIZE_FOR_ANGLE_COMPARE_JOINING_SEGMENTS);

			smallLines= bitmap.combineLines(smallLines, MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS, MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS_FULL, MAX_POINT_DISTANCE_TO_BE_PART_OF_MULTI_NODE,MAX_ANGLE_FOR_JOINING_SEGMENTS,MIN_SIZE_FOR_ANGLE_COMPARE_JOINING_SEGMENTS);
			smallLines= thin.combineLines(smallLines, MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS, MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS_THIN, MAX_POINT_DISTANCE_TO_BE_PART_OF_MULTI_NODE,MAX_ANGLE_FOR_JOINING_SEGMENTS,MIN_SIZE_FOR_ANGLE_COMPARE_JOINING_SEGMENTS);

			List<Line2D> removedTinyLines =smallLines.stream()
					.map(l->l.getLine())
					.filter(GeomUtil.longerThan(MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS).negate())
					.collect(Collectors.toList());

			List<Point2D> removedTinyVertices = removedTinyLines.stream()
					.flatMap(l->Stream.of(l.getP1(),l.getP2()))
					.collect(Collectors.toList());



			smallLines=smallLines.stream()
					.filter(l->l.length()>MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS)
					.collect(Collectors.toList());



			List<Point2D> verts = smallLines.stream()
					.flatMap(l->l.streamPoints())
					.collect(Collectors.toList());
			//find average spacing from OCR shapes to closest vertex         
			double[] lDistOCRToLine=likelyOCR.stream()
					.map(s->Tuple.of(s,GeomUtil.findCenterOfShape(s)))
					.map(Tuple.vmap(p->Tuple.of(p,GeomUtil.findClosestPoint(verts, p))))
					.map(Tuple.vmap(t->t.k().distance(t.v())))
					.mapToDouble(t->t.v())
					.sorted()
					.toArray();
			
			List<Shape> extendableOCR = likelyOCR.stream()
					                             .map(s->Tuple.of(s,ocrAttempt.get(s)))
					                             .filter(t->t.v()!=null)
					                             .filter(t->t.v().size()>0)
					                             .map(Tuple.vmap(l->l.get(0).k().toString()))
					                             .filter(t->!t.v().equals("H"))
					                             .map(t->t.k())
					                             .collect(Collectors.toList());

			OptionalDouble avgDistOCRToLine = Optional.of(0)
					.filter(d->lDistOCRToLine.length>0)
					.map(d->lDistOCRToLine[lDistOCRToLine.length/2])
					.map(d->OptionalDouble.of(d))
					.orElse(OptionalDouble.empty());



			linesJoined=Stream.concat(
					bigLines.stream().map(l->LineWrapper.of(l)),
					smallLines.stream())
					.collect(Collectors.toList());
			
			


			double largestBond=smallLines.stream()
					.mapToDouble(l->l.length())
					.max()
					.orElse(0);

			double averageLine=smallLines.stream()
					.mapToDouble(l->l.length())
					.filter(d->d>ignoreTooSmall[0])
					.average()
					.orElse(0);

			if(largestBond>2.0*averageLine){
				largestBond=1.4*averageLine;
			}

			List<List<LineWrapper>> preprocess= GeomUtil.reduceMultiBonds(Arrays.asList(linesJoined), MAX_ANGLE_FOR_PARALLEL, MAX_DISTANCE_TO_MERGE_PARALLEL_LINES, MIN_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS,0,MAX_DELTA_LENGTH_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC, (l)->{})
					.stream()
					.map(t->t.k())
					.map(t->LineWrapper.of(t))
					.map(t->Tuple.of(t,likelyOCR.stream()
							                    .filter(s->s.contains(t.centerPoint()))
							                    .findAny()
							                    .isPresent()))
					.collect(Collectors.groupingBy(t->t.v()))
					.values()
					.stream()
					.map(tl->tl.stream().map(t->t.k()).collect(Collectors.toList()))
					.collect(Collectors.toList());


			List<LineWrapper> rejBondOrderLines=new ArrayList<>();
			
			
			
			
			linesOrder=GeomUtil.reduceMultiBonds(preprocess, MAX_ANGLE_FOR_PARALLEL, largestBond/3, MIN_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS, MIN_BIGGER_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS,MAX_DELTA_LENGTH_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC,
					(rejLine->rejBondOrderLines.add(LineWrapper.of(rejLine)))
					);

			
			List<Shape> growLines = linesOrder.stream()
											  .map(t->t.k())
											  .map(l->GeomUtil.growLine(l, 5))
											  .collect(Collectors.toList());
			

			List<Shape> rescueOCRCandidates = new ArrayList<>();


			List<Shape> connectedComponents = polygons.stream()
					.map(s->GeomUtil.growShape(s, 2))
					.collect(Collectors.toList());

			int reps=0;
			boolean tooLongBond=true;
			
			double maxRatioInitial=0.5;
			double maxTotalRatioInitial=1.4;
			
			double maxRatio=0.5;
			double maxTotalRatio=1.4;
			
			List<LineWrapper> dottedLines = new ArrayList<>();
			
			
			while(tooLongBond){
				if (Thread.currentThread().isInterrupted()){
				      throw new InterruptedException();
				}
				
				
				rescueOCRCandidates.clear();
				List<Tuple<Line2D,Integer>> linesOrderRestricted =linesOrder.stream()
						.filter(t->{
							Line2D l=t.k();
							return isInOCRShape.negate().test(l);
						})
						.collect(Collectors.toList());

				ctab = GeomUtil.getConnectionTable(linesOrderRestricted, extendableOCR, 
						maxRatioForIntersection, 
						maxCandidateRatioForIntersection,
						maxPerLineDistanceRatioForIntersection,
						minPerLineDistanceRatioForIntersection,
						maxCandidateRatioForIntersectionWithNeighbor,
						GeomUtil.longerThan(maxBondLength[0]).negate())
						.mergeNodesCloserThan(MAX_DISTANCE_BEFORE_MERGING_NODES);

				if(DEBUG)ctabRaw.add(ctab.cloneTab());
				
				ctab.getEdgesWhichMightBeWiggleLines()
					.forEach(t->{
						
						realRescueOCRCandidates.add(t.k().getLine());
						
						Point2D p1=t.k().getLine().getP1();
						Point2D p2=t.k().getLine().getP2();
						List<Node> rnodes= new ArrayList<>();
				    	List<Node> lnodes= new ArrayList<>();
				    	t.v().stream()
				    		 .flatMap(e->e.streamNodes())
				    	     .distinct()
				    	     .forEach(n->{
				    	    	if(n.getPoint().distance(p1) <n.getPoint().distance(p2)){
				    	    		n.setPoint(p1);
				    	    		rnodes.add(n);
				    	    	}else{
				    	    		n.setPoint(p2);
				    	    		lnodes.add(n);
				    	    	}
				    	     });
				    	
				    	if(rnodes.size()>0 && lnodes.size()>0){
				    		Edge e=ctab.addEdge(rnodes.get(0).getIndex(), lnodes.get(0).getIndex(), 1);
				    		e.setDashed(true);
				    		foundNewOCR[0]=true;
				    	}
					});
				ctab.mergeNodesCloserThan(MAX_DISTANCE_BEFORE_MERGING_NODES);
				ctab.standardCleanEdges();
				
				System.out.println("Looking for dashes:" + ctabRaw.size());
				if(DEBUG)ctabRaw.add(ctab.cloneTab());

				RunningAverage allDashLengths = new RunningAverage(2);
				
				ctab.getEdgesWhichMightBeDottedLines()
				    .stream()
				    .map(e->Tuple.of(e,GeomUtil.vertices(e.stream().map(e1->e1.getLine()).collect(Collectors.toList()))))
				    .map(Tuple.vmap(vts->GeomUtil.getPairOfFarthestPoints(vts)))
				    .forEach(s->{
				    	List<Node> rnodes= new ArrayList<>();
				    	List<Node> lnodes= new ArrayList<>();
				    	s.k().stream()
				    		 .peek(e->allDashLengths.add(e.getEdgeLength()))
				    	     .flatMap(e->e.streamNodes())
				    	     .distinct()
				    	     .forEach(n->{
				    	    	if(n.getPoint().distance(s.v()[0]) <n.getPoint().distance(s.v()[1])){
				    	    		n.setPoint(s.v()[0]);
				    	    		rnodes.add(n);
				    	    	}else{
				    	    		n.setPoint(s.v()[1]);
				    	    		lnodes.add(n);
				    	    	}
				    	     });
				    	
				    	dottedLines.add(LineWrapper.of(new Line2D.Double(s.v()[0],s.v()[1])));
				    	if(rnodes.size()>0 && lnodes.size()>0){
				    		Edge e=ctab.addEdge(rnodes.get(0).getIndex(), lnodes.get(0).getIndex(), 1);
				    		e.setDashed(true);
				    		
				    		foundNewOCR[0]=true;
				    	}
				    	
				    });
				ctab.standardCleanEdges();
				ctab.mergeNodesCloserThan(MAX_DISTANCE_BEFORE_MERGING_NODES);
				
				double avgDot=allDashLengths.computeAvg();
				
				if(foundNewOCR[0] && avgDot>ignoreTooSmall[0]){
					ignoreTooSmall[0]=avgDot*1.1;
					break;
				}else{
					foundNewOCR[0]=false;
				}

				for(Shape s: likelyOCR){
					ctab.mergeAllNodesInsideCenter(s, OCR_TO_BOND_MAX_DISTANCE);
				}

				if(DEBUG)ctabRaw.add(ctab.cloneTab());
				
			

				List<List<Node>> newNodesForMerge = new ArrayList<>();

				Function<List<Node>,Point2D> bestIntersectionPoint = (nl)->{
					Point2D center= GeomUtil.findCenterOfVertices(nl.stream().map(n->n.getPoint()).collect(Collectors.toList()));
					double radSq=Math.pow(ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE,2);
					

					List<Edge> el=nl.stream()
							.flatMap(n->n.getEdges().stream()) //Note that this double-counts edges sometimes. That's actually good in this case
							.collect(Collectors.toList());


					List<Point2D> intersections=GeomUtil.eachCombination(el)
							.flatMap(t->{
								if(t.k()==t.v())return Stream.of(t.k().getPoint1(),t.k().getPoint2());
								return Stream.of(GeomUtil.intersection(t.k().getLine(),t.v().getLine()));
							})
							.filter(p->p!=null)
							.filter(p->p.distanceSq(center)<radSq)
							.collect(Collectors.toList());

					if(!intersections.isEmpty()){
						//return GeomUtil.findCenterOfVertices(intersections);
						return GeomUtil.findCenterMostPoint(intersections);
					}else{
						return center;
					}
				};

				List<Point2D> mergedPoints = new ArrayList<Point2D>();

				ctab.mergeNodesCloserThan(ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE_INITIAL, (nl)->{
					Point2D[] far=GeomUtil.getPairOfFarthestPoints(nl.stream().map(n->n.getPoint()).collect(Collectors.toList()));

					//This is pretty hacky (as if most of this code isn't)

					if(nl.size()==2 && far[0].distance(far[1])>ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE){
						return null;
					}

					if(far[0].distance(far[1])>0.9*ctab.getAverageBondLength()){
						//something is wrong here, there must be some bad nodes
						//flag for merge
						
						//Maybe still add as candidate? IDK.

						List<Node> group1=new ArrayList<>();
						List<Node> group2=new ArrayList<>();
						nl.forEach(n->{
							if(n.getPoint().distance(far[0])<n.getPoint().distance(far[1])){
								group1.add(n);
							}else{
								group2.add(n);
							}
						});

						newNodesForMerge.add(group1);
						newNodesForMerge.add(group2);
						//cancel the merge
						return null;
					}

					//find any vertices that are close to this one left over

					Point2D cpt=GeomUtil.findCenterOfVertices(nl.stream().map(n->n.getPoint()).collect(Collectors.toList()));
					List<Point2D> missingPoints=removedTinyVertices.stream()
							.filter(pt->pt.distance(cpt)<ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE)
							.collect(Collectors.toList());
					missingPoints.addAll(nl.stream().map(n->n.getPoint()).collect(Collectors.toList()));



					if(missingPoints.size()>3){
						Shape candidate=GeomUtil.convexHull2(missingPoints.stream().toArray(i->new Point2D[i]));
						//	        		double area=GeomUtil.area(candidate);
						
						if(GeomUtil.area(candidate)>0.5*averageAreaOCR){
							Point2D center=GeomUtil.findCenterOfVertices(missingPoints);;

							candidate=GeomUtil.growShape(candidate,4);
							rescueOCRCandidates.add(candidate);
							
							return center;
						}
					}

					mergedPoints.addAll(missingPoints);
					
					return bestIntersectionPoint.apply(nl);

				});

				newNodesForMerge.forEach(ln->{
					Point2D center= bestIntersectionPoint.apply(ln);
					ctab.mergeNodes(ln.stream().map(n->n.getIndex()).collect(Collectors.toList()), (ll)->{
						return center;
					});
				});


				System.out.println("Merging:" + ctabRaw.size());
				if(DEBUG)ctabRaw.add(ctab.cloneTab());


				//ctab.mergeAllNodesOnParLines();
				ctab.removeOrphanNodes();
				ctab.standardCleanEdges();

				ctab.mergeNodesCloserThan(ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE, (nl)->{
					Point2D cpt=GeomUtil.findCenterOfVertices(nl.stream().map(n->n.getPoint()).collect(Collectors.toList()));
					List<Point2D> missingPoints=removedTinyVertices.stream()
							.filter(pt->pt.distance(cpt)<ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE)
							.collect(Collectors.toList());
					missingPoints.addAll(nl.stream().map(n->n.getPoint()).collect(Collectors.toList()));



					if(missingPoints.size()>3){
						Shape candidate=GeomUtil.convexHull2(missingPoints.stream().toArray(i->new Point2D[i]));
						double area=GeomUtil.area(candidate);
						if(area>0.5*averageAreaOCR){
							candidate=GeomUtil.growShape(candidate,4);
							rescueOCRCandidates.add(candidate);
							//polygons.add(candidate);
							
							return GeomUtil.findCenterOfVertices(missingPoints);
						}
						//polygons.add(candidate);
						//return center;
					}

					mergedPoints.addAll(missingPoints);

					return bestIntersectionPoint.apply(nl);

				});


				GeomUtil.groupThings(mergedPoints, (tp)->{
		
							Point2D p1=tp.k();
							Point2D p2=tp.v();
							if(p1.distance(p2)<ctab.getAverageBondLength()*0.6){
								return connectedComponents.stream()
										.filter(s->s.contains(p1) && s.contains(p2))
										.findAny()
										.isPresent();
							}
							return false;
		
						})
						.forEach(ll->{
							Point2D[] pts=ll.toArray(new Point2D[0]);
							if(pts.length>3){
								Shape candidate=GeomUtil.convexHull2(pts);
								if(GeomUtil.area(candidate)>0.5*averageAreaOCR){
									candidate=GeomUtil.growShape(candidate,4);
									//polygons.add(candidate);
									rescueOCRCandidates.add(candidate);
								}
							}
						});

				ctab.removeOrphanNodes();
				ctab.standardCleanEdges();


				if(DEBUG)ctabRaw.add(ctab.cloneTab());
				
				
				Set<Edge> splitEdges = new HashSet<Edge>();

				ctab.createNodesOnIntersectingLines(2, elist->{
					splitEdges.addAll(elist);
					return true;
				}, (nn)->{
					intersectionNodes.add(nn.getPoint());
				});
				List<Shape> expectedLineZones =Stream.concat(growLines.stream(), likelyOCR.stream()).collect(Collectors.toList());
				
				ctab.getEdges()
				    .stream()
				    .filter(e->splitEdges.contains(e))
				    .collect(Collectors.toList())
				    .forEach(e->{
				    	Line2D ll=e.getLine();
				    	
				    	double totLen=GeomUtil.getLinesNotInside(ll,expectedLineZones).stream().mapToDouble(l->GeomUtil.length(l)).sum();
				    	if(totLen>GeomUtil.length(ll)*0.8){
				    		ctab.removeEdge(e);
				    	}
				    });

				if(DEBUG)ctabRaw.add(ctab.cloneTab());
				
				

				ctab.mergeFilteredNodesCloserThan(ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE_AFTER_SPLIT, n->true);
				ctab.standardCleanEdges();
				

				if(DEBUG)ctabRaw.add(ctab.cloneTab());

				

				if(avgDistOCRToLine.isPresent()){
					double nmaxRatio=(avgDistOCRToLine.getAsDouble())/ctab.getAverageBondLength();

					if(nmaxRatio>maxRatioInitial){
						//TODO: sometimes there are only pairs of bonds, in such a case the heuristics don't
						//work well
						maxRatio=nmaxRatio;

						double maxlen=ctab.getEdges()
								.stream()
								.mapToDouble(e->e.getEdgeLength())
								.max()
								.orElse(1);

						maxlen=Math.max(maxlen, averageWidthOCR);

						maxTotalRatio = Math.max(maxTotalRatioInitial, maxlen/ctab.getAverageBondLength());
					}
				}

				ctab.mergeNodesExtendingTo(likelyOCR,maxRatio,maxTotalRatio);
				
				if(DEBUG)ctabRaw.add(ctab.cloneTab());

				ctab.removeOrphanNodes();

				ctab.mergeFilteredNodesCloserThan(ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE_AFTER_SPLIT, n->true);
				ctab.standardCleanEdges();
				ctab.mergeFilteredNodesCloserThan(ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE, n->{
					if(intersectionNodes.stream().filter(p->p.distance(n.getPoint())<0.01*ctab.getAverageBondLength()).findAny().isPresent()){
						if(n.getEdgeCount()==4)return false;
						return true;
					}
					return true;
				});
				if(DEBUG)ctabRaw.add(ctab.cloneTab());
				
				for(Shape s: likelyOCR){
					ctab.mergeAllNodesInsideCenter(s, OCR_TO_BOND_MAX_DISTANCE);
				}

				
				System.out.println("Finding new bonds:" + ctabRaw.size());
				if(DEBUG)ctabRaw.add(ctab.cloneTab());

				ctab.makeMissingNodesForShapes(likelyOCR,MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL,MIN_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL);

				Set<Node> toRemove = new LinkedHashSet<Node>();
				Set<Edge> toRemoveEdges = new LinkedHashSet<Edge>();
				Set<Edge> toRemoveEdgesImmediately = new LinkedHashSet<Edge>();
				


				ctab.makeMissingBondsToNeighbors(bitmap,MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL,MAX_TOLERANCE_FOR_DASH_BONDS,likelyOCR,OCR_TO_BOND_MAX_DISTANCE, (t)->{
					
					//It could be that there is already a bond between the two nodes through a bad intermediate
					Edge e=t.v();
					Node n1=e.getRealNode1();
					Node n2=e.getRealNode2();
					if(intersectionNodes.stream().filter(p->p.distance(n1.getPoint())<3).findAny().isPresent()){
						if(n1.getEdgeCount()>=4){
							toRemoveEdgesImmediately.add(e);
							return;
						}
					}
					if(intersectionNodes.stream().filter(p->p.distance(n2.getPoint())<3).findAny().isPresent()){
						if(n2.getEdgeCount()>=4){
							toRemoveEdgesImmediately.add(e);
							return;
						}
					}
					List<Edge> existingEdges1=n1.getEdges();
					List<Edge> existingEdges2=n2.getEdges();
					Set<Node> n1Neigh=existingEdges1.stream()
							.flatMap(ne->Stream.of(ne.getRealNode1(),ne.getRealNode2()))
							.filter(n->!n.equals(n1))
							.collect(Collectors.toSet());

					Set<Node> n2Neigh=existingEdges2.stream()
							.flatMap(ne->Stream.of(ne.getRealNode1(),ne.getRealNode2()))
							.filter(n->!n.equals(n2))
							.collect(Collectors.toSet());
					List<Node> commonNeigh = n1Neigh.stream().filter(nn->n2Neigh.contains(nn)).collect(Collectors.toList());

					boolean alreadyExists = false;

					if(!commonNeigh.isEmpty()){
						for(Node cn:commonNeigh){
							Point2D cp=cn.getPoint();
							double distance1=n1.getPoint().distance(cp);
							double distance2=n2.getPoint().distance(cp);
							double sumd=distance1+distance2;
							double ddelta=Math.abs(sumd-e.getEdgeLength());
							List<Edge> edges=cn.getEdges();
							if(edges.size()==2){
								if(ddelta<MAX_DELTA_LENGTH_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC){
									if(!toRemove.contains(n1) && !toRemove.contains(n2)){
										toRemove.add(cn);
	
									}
									double o2=edges.stream().map(et->Tuple.of(et,et.getEdgeLength()))
											.mapToDouble(e1->(e1.k().getOrder() * e1.v()))
											.sum();
									int o=(int)Math.round(((o2/sumd)+0.05));
									t.v().setOrder(o);

								}
								if(!edges.stream().anyMatch(e2->e2.getDashed())){
									alreadyExists=true;	
								}
							}else if(edges.size()==4){ //might be intersection for cage
								boolean isIntersection = intersectionNodes.stream()
								                 .filter(in->in.distance(cp)<2)
								                 .findAny()
								                 .isPresent();
								if(isIntersection && Math.abs(sumd-t.v().getEdgeLength())<ctab.getAverageBondLength()*0.05){
									toRemoveEdges.add(t.v());
								}
								        
								
							}

						}
					}


					
					if(t.k()>MAX_TOLERANCE_FOR_SINGLE_BONDS){
						if(!alreadyExists){
							t.v().setDashed(true);
						}
					}
				});

				toRemoveEdgesImmediately.forEach(e->ctab.removeEdge(e));
				
				
				System.out.println("Looking for rings:" + ctabRaw.size());
				if(DEBUG)ctabRaw.add(ctab.cloneTab());
				
				ctab.getRings()
			    .stream()
			    .filter(r->r.size()>4 && r.size()<7)
			    .forEach(r->{
			    	Point2D center=GeomUtil.centerOfMass(r.getConvexHull());
			    	Shape p=GeomUtil.convexHull2(GeomUtil.makeNPolyCenteredAt(new Point2D.Double(0,0), r.size(), 100));
			    	//Shape p2=GeomUtil.convexHull2(GeomUtil.makeNPolyCenteredAt(center, 6, ctab.getAverageBondLength()));
			    	Point2D anchor = r.getNodes().stream()
			    			          .map(n->Tuple.of(n, n.getPoint().distance(center)).withVComparator())
			    			          .max(Comparator.naturalOrder())
			    			          .map(t->t.k().getPoint())
			    			          .orElse(null);
			    	Line2D nline = new Line2D.Double(center,anchor);
			    	AffineTransform at=GeomUtil.getTransformFromLineToLine(new Line2D.Double(new Point2D.Double(0,0),new Point2D.Double(100,0)),nline,false);
			    	Shape ns=at.createTransformedShape(p);
			    	
			    	double ll=GeomUtil.length(nline)*0.02;
			    	
			    	Point2D[] verts2 = GeomUtil.vertices(ns);
			    	
			    	boolean looksOkay = r.getNodes()
			    	 .stream()
			    	 .map(n->{
			    		double dd=Arrays.stream(verts2)
			    		      .map(v->Tuple.of(v, v.distance(n.getPoint())).withVComparator())
			    		      .min(Comparator.naturalOrder())
			    		      .map(t->t.v())
			    		      .orElse(0.0);
			    		return dd;
			    	 })
			    	 .filter(d->d>ll)
			    	 .findAny()
			    	 .isPresent();
			    	
			    	if(looksOkay){
			    		//this means the rings are likely to be real.
			    		//with that in mind, let's take a look and see if there are things that should be merged
			    		//specifically, we want edges that are dashes
			    		Set<Edge> eset = r.getEdges().stream().collect(Collectors.toSet());
			    		toRemoveEdges.removeAll(eset);
			    		toRemove.removeAll(r.getNodes());
			    		
			    		eset
			    		 .stream()
			    		 //.filter(ee->ee.getDashed())
			    		 .forEach(ed->{
			    			 
			    			 boolean[] changed=new boolean[]{false};
			    			 
			    			 ed.getRealNode1()
			    			   .getNeighborNodes()
			    			   .stream()
			    			   .filter(ne->!eset.contains(ne.v()))
			    			   .filter(ne->ne.k().distanceTo(ed.getRealNode2()) < ed.getEdgeLength())
			    			   .filter(ne->{
			    				   return GeomUtil.cosTheta(ne.v().getLine(), ed.getLine())>Math.cos(20*Math.PI/180.0);
			    			   })
			    			   .filter(ne->{
			    				   return ns.contains(ne.k().getPoint());
			    			   })
			    			   .forEach(ne->{
			    				   ed.setDashed(ne.v().getDashed());
			    				   if(ed.getOrder()==1 && ne.v().getOrder()!=1){
			    					   ed.setOrder(ne.v().getOrder());   
			    				   }
			    				   toRemoveEdges.add(ne.v());
			    				   toRemove.add(ne.k());
			    				   changed[0]=true;
			    			   });
			    			 ed.getRealNode2()
			    			   .getNeighborNodes()
			    			   .stream()
			    			   .filter(ne->!eset.contains(ne.v()))
			    			   .filter(ne->ne.k().distanceTo(ed.getRealNode1()) < ed.getEdgeLength())
			    			   .filter(ne->{
			    				   return GeomUtil.cosTheta(ne.v().getLine(), ed.getLine())>Math.cos(20*Math.PI/180.0);
			    			   })
			    			   .filter(ne->{
			    				   return ns.contains(ne.k().getPoint());
			    			   })
			    			   .forEach(ne->{
			    				   ed.setDashed(ne.v().getDashed());
			    				   if(ed.getOrder()==1 && ne.v().getOrder()!=1){
			    					   ed.setOrder(ne.v().getOrder());   
			    				   }
			    				   toRemoveEdges.add(ne.v());
			    				   toRemove.add(ne.k());
			    				   changed[0]=true;
			    			   });
			    			 
			    			 if(changed[0] && ed.getOrder()!=2){
			    				 LineWrapper lw1= LineWrapper.of(ed.getLine());
			    				 Point2D centerpt=lw1.centerPoint();
			    				 linesOrder.stream()
			    				           .filter(t->t.v()==2)
			    				           .map(t->t.k())
			    				           .map(t->LineWrapper.of(t))
			    				           .filter(lw->lw.length()>ctab.getAverageBondLength()*0.4)
			    				           .filter(lw->lw.absCosTheta(lw1)>Math.cos(Math.PI*10.0/180.0))
			    				           .map(lw->lw.growLine(ctab.getAverageBondLength()*0.4))
			    				           .filter(s->s.contains(centerpt))
			    				           .findFirst()
			    				           .ifPresent(s->{
			    				        	  ed.setOrder(2);
			    				        	  ed.setDashed(false);
			    				        	  ed.setWedge(false);
			    				           });
			    			 }
			    			 
			    			  
			    			 
			    			 
			    		 });
			    		
			    		
			    		
			    	}			    	
			    });

				
				
				toRemoveEdges.forEach(e->ctab.removeEdge(e));
				
				//fuzzy adding missing stuff
				System.out.println("Removing bad edges:" + ctabRaw.size());
				if(DEBUG)ctabRaw.add(ctab.cloneTab());
				
				
				toRemove.forEach(n->ctab.removeNodeAndEdges(n));
				//ctab.removeOrphanNodes();

			

				System.out.println("Checking if bond length looks okay:" + ctabRaw.size());
				if(DEBUG)ctabRaw.add(ctab.cloneTab());

				double avgBondLength=ctab.getAverageBondLength();
				maxBondLength[0]=avgBondLength*MAX_BOND_TO_AVG_BOND_RATIO_TO_KEEP;



				
				Predicate<Line2D> longerThanMax = GeomUtil.longerThan(maxBondLength[0]);
				tooLongBond = ctab.getEdges()
						.stream()
						.filter(e->longerThanMax.test(e.getLine()))
						.findAny()
						.isPresent();
				if(tooLongBond){
					reps++;
				}
				if(reps>MAX_REPS)break;
			}
			//realRescueOCRCandidates.addAll(rescueOCRCandidates);


			AtomicBoolean anyOtherIntersections = new AtomicBoolean(false);

			ctab.createNodesOnIntersectingLines(3, elist->{

				long longEnoughBonds=elist.stream()
						
						.filter(e->e.getEdgeLength()>MAX_BOND_TO_AVG_BOND_RATIO_FOR_INTERSECTION*ctab.getAverageBondLength())
						.count();
				if(longEnoughBonds<3)return false;
				anyOtherIntersections.set(true);
				return true;
			}, (nn)->{
				intersectionNodes.add(nn.getPoint());
			});

			if(anyOtherIntersections.get()){
				
				if(DEBUG)ctabRaw.add(ctab.cloneTab());
				ctab.mergeNodesCloserThan(ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE);
				ctab.standardCleanEdges();
			}




			double shortestRealBondRatio = .3;
			ctab.fixBondOrders(likelyOCR,shortestRealBondRatio, e->{
				e.setOrder(1);
			});





			//This is probably where we try to add some missed OCR based on the nodes

			//1. Find all nodes that aren't in likely OCR shapes
			//2. Create a bounding circle around the node that's ~half a bond width
			//3. See how many initial segment points were found inside. We want there to have been more than the number of edges to that node
			//4. Take all those segment points found, and find the center.
			//5. repeat step 2 until no changes are made to the point set
			//6. make a convex hull from the shapes, and "grow" it a little bit (configurable)
			//7. Crop the original bitmap (and thin, I guess) and run OCR on that chunk
			//8. Add everything that is OKAY to the OCR sets

			double avgbond=ctab.getAverageBondLength();
			double SEED_BOND_RATIO_FOR_OCR_WIDTH=0.0;

			double SEED_BOND_RATIO_FOR_OCR_WIDTH_FOR_CENTROID=0.5;

			List<Node> unmatchedNodes=ctab.getNodesNotInShapes(likelyOCR, OCR_TO_BOND_MAX_DISTANCE + avgbond*SEED_BOND_RATIO_FOR_OCR_WIDTH);

			List<LineWrapper> verticesJl=linesJoined.stream()
					.collect(Collectors.toList());
			List<Point2D> verticesJ=
					verticesJl.stream()
					.flatMap(l->l.streamPoints())
					.collect(Collectors.toList());



			List<Shape> toAddAllOCR=new ArrayList<Shape>();

			Map<Shape,List<Tuple<Character,Number>>> gotCache=new HashMap<>();


			unmatchedNodes.forEach(n->{

				Point2D cpt=n.getPoint();

				//Do another step for likely candidates


				Point2D centerRescue=rescueOCRCandidates.stream()
						.filter(sc->sc.contains(n.getPoint()))
						.map(sc->Tuple.of(sc,GeomUtil.area(sc)).withVComparator())
						.sorted(Comparator.reverseOrder())
						.map(t->t.k())
						.findAny()
						.map(sc->GeomUtil.findCenterOfShape(sc))
						.orElse(null);



				int numEdges=n.getEdges().size();
				Shape[] area=new Shape[]{null};
				Shape nshape=null;
				double radius=Math.max(avgbond*SEED_BOND_RATIO_FOR_OCR_WIDTH_FOR_CENTROID,averageLargestOCR/2);

				List<Point2D> allVertices = verticesJ;

				if(centerRescue!=null){
					cpt=centerRescue;
				}

				boolean keep=true;
				for(int i=0;i<3;i++){
					keep=true;
					area[0]=GeomUtil.makeShapeAround(cpt,radius);

					//
					List<Point2D> insideVertices=allVertices.stream()
							.filter(v->area[0].contains(v))
							.collect(Collectors.toList());
					if(insideVertices.size()<=numEdges+1){
						keep=false;
					}
					List<Point2D> insideVertices2=allVertices.stream()
							.filter(v->area[0].contains(v))
							.collect(Collectors.toList());

					Point2D center=GeomUtil.findCenterOfVertices(insideVertices2);;
					//remove outliers

					double distanceMean= insideVertices2.stream()
							.mapToDouble(pt->center.distance(pt))
							.average()
							.orElse(0);
					double distanceVar= insideVertices2.stream()
							.mapToDouble(pt->Math.pow(distanceMean - center.distance(pt),2))
							.sum();
					double distanceSTDEV=Math.sqrt(distanceVar/(insideVertices2.size()-1));
					//distanceStDev = Math.sqrt(distanceMean*distanceMean-distanceStDev);

					List<Point2D> realMissing=insideVertices2.stream()
							.filter(pt->center.distance(pt)<distanceMean+distanceSTDEV*2.5)
							.collect(Collectors.toList());

					nshape = GeomUtil.convexHull2(realMissing.toArray(new Point2D[0]));


					Point2D[] far=GeomUtil.getPairOfFarthestPoints(nshape);

					double arean = GeomUtil.area(nshape);

					double r=0;
					if(far!=null){
						r=far[0].distance(far[1]);	
					}
					if(r < averageLargestOCR*MIN_LONGEST_WIDTH_RATIO_FOR_OCR_TO_AVERAGE){
						keep=false;
					}

					if(arean < GeomUtil.area(nshape.getBounds2D())*MIN_AREA_RATIO_FOR_HULL_TO_BBOX_OCR){
						keep=false;
					}
					if(GeomUtil.area(nshape.getBounds2D())>avgbond*avgbond*0.5){
						keep=false;
					}
					if(GeomUtil.area(nshape.getBounds2D()) < averageAreaOCR*MIN_AREA_RATIO_FOR_OCR_TO_AVERAGE){
						keep=false;
					}
					
					if(GeomUtil.area(nshape.getBounds2D()) > averageAreaOCR*MAX_AREA_RATIO_FOR_OCR_TO_AVERAGE){
						keep=false;
					}

					radius=Math.max(averageLargestOCR/2,r/2);
					//cpt=GeomUtil.findCenterOfVertices(Arrays.asList(GeomUtil.vertices(nshape)));
					cpt=GeomUtil.findCenterOfShape(nshape);
				}



				if(keep){



					Bitmap nmap=bitmap.getLazyCrop(nshape);
					Bitmap nthinmap=thin.getLazyCrop(nshape);
					if(nmap!=null && nthinmap!=null){
						
						nshape=GeomUtil.growShape(nshape, 2);
						
						nmap=bitmap.crop(nshape);

						List<Shape> slist=nmap.connectedComponents(Bitmap.Bbox.DoublePolygon);


						Shape bshape=slist.stream()
								.map(s->Tuple.of(s,s.getBounds2D().getWidth()*s.getBounds2D().getHeight()).withVComparator())
								.max(CompareUtil.naturalOrder())
								.map(t->t.k())
								.orElse(nshape);
						Rectangle2D rect1 = nshape.getBounds2D();
						AffineTransform at = new AffineTransform();
						at.translate(rect1.getX(),rect1.getY());
						nshape=at.createTransformedShape(bshape).getBounds2D();
						//                    

						nmap=bitmap.getLazyCrop(nshape);
						nthinmap=thin.getLazyCrop(nshape);
						
						if(GeomUtil.area(nshape)<averageAreaOCR*MIN_AREA_RATIO_FOR_OCR_TO_AVERAGE){
							return;
						}

						if(nmap!=null && nthinmap!=null){
							processOCRShape(socr[0],nshape,bitmap,thin,(s,potential)->{

								if(potential.get(0).v().doubleValue()>OCRcutoffCosineRescue){

									String st=potential.get(0).k().toString();
									if(BranchNode.interpretOCRStringAsAtom2(st)!=null){
										toAddAllOCR.add(s);	
										gotCache.put(s,potential);
									}
								}	
							});
						}
					}
				}

				//}

			});

			
			if(DEBUG)ctabRaw.add(ctab.cloneTab());


			GeomUtil.mergeOverlappingShapes(toAddAllOCR,0.75)
			.forEach(nshape->{
				//if(true)return;
				boolean sigOverlap = 
						likelyOCRAll.stream()
						.map(s->Tuple.of(s,GeomUtil.getIntersectionShape(nshape, s)))
						.filter(os->os.v().isPresent())
						.map(Tuple.vmap(os->os.get()))
						.map(Tuple.vmap(s->GeomUtil.area(s)))
						.map(Tuple.kmap(s->GeomUtil.area(s)))
						.mapToDouble(t->t.v()/t.k())
						.filter(areaFraction->areaFraction>0.0)
						.findAny()
						.isPresent();

				if(sigOverlap){
					return;
				}
				//if(ctab.getNodesInsideShape(nshape, 0).isEmpty())return;
				List<Tuple<Character,Number>> matches=gotCache.getOrDefault(nshape, new ArrayList<>());
				if(matches.isEmpty()){
					Bitmap nmap=bitmap.getLazyCrop(nshape);
					Bitmap nthinmap=thin.getLazyCrop(nshape);
					if(nmap!=null && nthinmap!=null){
						processOCRShape(socr[0],nshape,bitmap,thin,(s,potential)->{
							matches.addAll(potential);
						});
					}
				}
				ocrAttempt.put(nshape, matches);
				
				//polygons.add(nshape);
				if (matches.get(0).v().doubleValue() > OCRcutoffCosineRescue) {
					
					
					CharType ct=OCRIsLikely(matches.get(0));
					if(ct.equals(CharType.ChemLikely)){
						likelyOCR.add(nshape);
						likelyOCRNonBond.add(nshape);
						
					}else if(ct.equals(CharType.NumericLikely)){
						likelyOCRNonBond.add(nshape);
						likelyOCRNumbers.add(nshape);
					}
					likelyOCRAll.add(nshape);
					ocrRescues.add(nshape);
					
					foundNewOCR[0]=true;
				}
			});


			
			ctab.mergeNodesExtendingTo(likelyOCR,maxRatio,maxTotalRatio);

//			System.out.println("Extended to OCR:" + ctabRaw.size());
			if(DEBUG)ctabRaw.add(ctab.cloneTab());

			double cosThetaOCRShape =Math.cos(MAX_THETA_FOR_OCR_SEPERATION);
			
		

			List<Shape> growLikelyOCRNonBond=likelyOCRNonBond.stream().map(s->GeomUtil.growShape(s, 2)).collect(Collectors.toList());
			List<Shape> growLikelyOCR=likelyOCR.stream().map(s->GeomUtil.growShape(s, 2)).collect(Collectors.toList());

			List<Shape> maybeDash = polygons.stream()
					.filter(s->GeomUtil.area(s)<averageAreaOCR)
			        .filter(s->!likelyOCR.contains(s))
			        .map(s->Tuple.of(s,GeomUtil.findCenterOfShape(s)))
			        .filter(st->!growLikelyOCR.stream().filter(g->g.contains(st.v())).findFirst().isPresent())
			        .map(t->t.k())
			        .map(t->Tuple.of(t,GeomUtil.findLongestSplittingLine(t)))
			        
			        .filter(t->{
			        	if(t.v()!=null){
			        		return t.v().length()<ctab.getAverageBondLength()*0.5;
			        	}
			        	return true;
			        })
			        .map(t->t.k())
			        .collect(Collectors.toList());

			        
			List<Shape> maybeDashCollection = GeomUtil.groupThings(maybeDash, t->{
									Shape s1=t.k();
									Shape s2=t.v();
									Point2D p1=GeomUtil.findCenterOfShape(s1);
									Point2D p2=GeomUtil.findCenterOfShape(s2);
									return p1.distance(p2)<ctab.getAverageBondLength()/3;
								})
								.stream()
								.filter(sl->sl.size()>=3)
								.map(l->{
									LineWrapper splitting=GeomUtil.findLongestSplittingLine(l.stream().collect(GeomUtil.joined()));
									return l.stream()
											.map(l1->Tuple.of(l1, GeomUtil.findLongestSplittingLine(l1)))
									        .filter(l2->splitting.cosTheta(l2.v()) < Math.cos(45*Math.PI/180))
									        .map(t->t.k())
									        .collect(Collectors.toList());
									
								})
								.map(l->l.stream().collect(GeomUtil.joined()))
								.filter(b->b!=null)
								.filter(bshape->{
									Point2D[] pts=GeomUtil.getPairOfFarthestPoints(bshape);
									double dist=pts[0].distance(pts[1]);
									if(dist < ctab.getAverageBondLength()*1.3 && dist>ctab.getAverageBondLength()*0.6){
										return true;
									}else{
										return false;
									}
								})
								.collect(Collectors.toList());
			
			likelyOCRAll.stream()
			            .filter(p->maybeDashCollection.stream().filter(ds->ds.contains(GeomUtil.findCenterOfShape(p))).findAny().isPresent())
			            .collect(Collectors.toList())
			            .forEach(ocs->{
			            	likelyOCRAll.remove(ocs);
			            	likelyOCR.remove(ocs);
			            	likelyOCRNumbers.remove(ocs);
			            	likelyOCRNonBond.remove(ocs);
			            });
			
			double wid=averageWidthOCR;

			List<List<Shape>> ocrGroupList=GeomUtil.groupShapesIfClosestPointsMatchCriteria(likelyOCRAll, t->{
				Point2D[] pts=t.v();
				Shape[] shapes =t.k();

				//Line2D l2 = new Line2D.Double(pts[0],pts[1]);
				double dist=pts[0].distanceSq(pts[1]);
				double cutoff=Math.max(ctab.getAverageBondLength()*MAX_BOND_RATIO_FOR_OCR_CHAR_SPACING,wid);
				
				if(dist>cutoff*cutoff){
					return false;
				}
				
				List<Tuple<Character, Number>> attempt0 = ocrAttempt.get(shapes[0]);
				List<Tuple<Character, Number>> attempt1 = ocrAttempt.get(shapes[1]);
				String v1= (attempt0 ==null || attempt0.isEmpty())? "" : attempt0.get(0).k().toString();
				String v2= (attempt1 ==null || attempt1.isEmpty())? "" : attempt1.get(0).k().toString();
				if(v1.equals("\\") || v1.equals("/") || 
						v2.equals("\\") || v2.equals("/")){
					return false;
				}
				
				
				

				
				double s1y1=shapes[0].getBounds2D().getMinY();
				double s1y2=shapes[0].getBounds2D().getMaxY();
				double s2y1=shapes[1].getBounds2D().getMinY();
				double s2y2=shapes[1].getBounds2D().getMaxY();
				
				if(!((s1y1>=s2y1 && s1y1<=s2y2) || 
				     (s1y2>=s2y1 && s1y2<=s2y2) ||
				     ((s1y2+s1y1)/2>=s2y1 && (s1y2+s1y1)/2<=s2y2))
						){
					return false;
				}
				
				double pvoverlapT = Math.max(s1y2, s2y2) - Math.min(s1y1, s2y1);
				double pvoverlapI = Math.min(s1y2, s2y2) - Math.max(s1y1, s2y1);
				
				if(pvoverlapI<pvoverlapT*0.4){
					if(likelyOCR.contains(shapes[0]) && likelyOCR.contains(shapes[1])){
						if(!v1.equalsIgnoreCase("S") && !v2.equalsIgnoreCase("S")){
							return false;	
						}
						
					}
				}
				
				

				Point2D cs1=GeomUtil.findCenterOfShape(shapes[0]);
				Point2D cs2=GeomUtil.findCenterOfShape(shapes[1]);

				LineWrapper cenLine = LineWrapper.of(new Line2D.Double(cs1, cs2));
				double[] vec=cenLine.vector();
				double cosTheta=Math.abs(vec[0]*cenLine.recipLength());

				if(cosTheta>cosThetaOCRShape){
					return true;
				}
				return false;
			});

			bestGuessOCR.clear();
			
						
			boolean areLikelyNumbers=(likelyOCRNumbers.size()>=4);

			ocrGroupList.stream()
			//.filter(l->l.size()>1)
			.forEach(g->{
				List<Shape> sorted=g.stream()
						.map(s->Tuple.of(s,s))
						.map(Tuple.vmap(s->s.getBounds2D().getMinX()))
						.map(t->t.withVComparator())
						.sorted()
						.map(t->t.k())
						.collect(Collectors.toList());



				String soFar="";
				Shape making=null;
				
				
				List<Tuple<Shape,Tuple<List<Shape>,String>>> toAdd = new ArrayList<>();
				List<Shape> lsofar = new ArrayList<>();

				for(Shape s: sorted){
					lsofar.add(s);
					List<Tuple<Character, Number>> list = ocrAttempt.get(s);
					String v= (list ==null || list.isEmpty())? "":list.get(0).k().toString();
					
					if(s.getBounds2D().getHeight()<averageHeightOCR*0.8 || (s.getBounds2D().getHeight() <=averageHeightNumberOCR*1.1 && areLikelyNumbers)){
						if(v.equalsIgnoreCase("S")){
							
							//probably a 3 or 8
							Tuple<Character,Number> tc=list.stream()
									                       .filter(c->c.k().toString().equals("3")||c.k().toString().equals("8"))
									                       .findFirst().orElse(null);
							if(tc!=null){
								if(tc.v().doubleValue()>OCRcutoffCosineRescue){
									v=tc.k().toString();
								}
							}
						}
					}
					if(v.equals("?")){
						v="N+";
					}

					if(v.equals("-")){
						if(!soFar.equals("t")){
							if(making!=null){
								toAdd.add(Tuple.of(making,Tuple.of(lsofar,soFar)));
								lsofar.remove(lsofar.size()-1);
								lsofar= new ArrayList<>();
							}
							soFar="";
							making=null;
							continue;
						}
					}
					if(making==null){
						making=s;
					}else{
						making=GeomUtil.add(making,s);
					}
					soFar+=v;
				}

				if(making!=null){
					toAdd.add(Tuple.of(making,Tuple.of(lsofar,soFar)));
				}
				
				Map<String,List<String>> dontMerge = new HashMap<>();
				
				dontMerge.put("OO", Arrays.asList("O","O"));
				dontMerge.put("Oo", Arrays.asList("O","o"));
				dontMerge.put("oO", Arrays.asList("o","O"));
				dontMerge.put("oo", Arrays.asList("o","o"));
				
				dontMerge.put("SO", Arrays.asList("S","O"));
				dontMerge.put("OS", Arrays.asList("O","S"));
				dontMerge.put("so", Arrays.asList("s","o"));
				dontMerge.put("os", Arrays.asList("o","s"));
				
				
				
				dontMerge.put("OF", Arrays.asList("O","F"));
				dontMerge.put("oF", Arrays.asList("o","F"));
				dontMerge.put("Fo", Arrays.asList("F","o"));
				dontMerge.put("FO", Arrays.asList("F","O"));
				dontMerge.put("FF", Arrays.asList("F","F"));
				dontMerge.put("CH3CH3", Arrays.asList("CH3","CH3"));
				dontMerge.put("cH3cH3", Arrays.asList("cH3","cH3"));
				dontMerge.put("HOOH", Arrays.asList("HO","OH"));
				dontMerge.put("OHOH", Arrays.asList("OH","OH"));
				dontMerge.put("OHHO", Arrays.asList("OH","HO"));
				dontMerge.put("BrBr", Arrays.asList("Br","Br"));

				for(Tuple<Shape,Tuple<List<Shape>,String>> tt: toAdd){
					boolean removeBad=false;
					String val = tt.v().v();
					List<Shape> contains = tt.v().k();
					Shape parent=tt.k();
					
					
					if(val.contains("~")){
						if(val.equals("~")){
							val="Cl";
						}else{
							val=val.replace("~", "O");
						}
					}
					
					if(val.contains("$")){
						val=val.replace("$", "O2");
					}
					if(val.contains("!")){
						val=val.replace("!", "H3");
					}
					
					//This is pretty hacky
					if(val.contains("`")){
						if(val.contains("``")){
							val=val.replace("``","`");
						}
						
						val=val.replace("`", "HO");
						
					}
					if(val.contains("%")){
						val=val.replace("%", "OC");
					}
					
					
					
					
					if(val.matches("[cC][h][1ilt][r][a][1ilt]")){
						removeBad=true;
					}
					
					if(val.equals("IN") || val.equals("tN") || val.equals("lN")){
						bestGuessOCR.put(contains.get(1), "N");
						continue;
					}
					
					
					if(dontMerge.containsKey(val)){
						List<String> keepAs=dontMerge.get(val);
						
						int findex=0;
						
						for(int i=0;i<keepAs.size();i++){
							String keep=keepAs.get(i);
							String g1=val.substring(findex,findex+keep.length());
							Shape parts= contains.stream().skip(findex).limit(keep.length()).collect(GeomUtil.joined());
							findex=findex+keep.length();
							if(keep.equals(g1)){
								bestGuessOCR.put(parts, keep);
							}
						}
						
						
//						bestGuessOCR.put(contains.get(0), val.substring(0,1));
//						bestGuessOCR.put(contains.get(1), val.substring(0,1));
						continue;
					}
					
					if(val.contains("F8")){
						val=val.replace("F8", "F3");
					}
					if(val.matches("[cC]H[Ss]")){
						val=val.replaceAll("[cC]H[Ss]", "CH3");
					}
					if(val.matches(".*[cC][tlI][Ss].*")){
						val=val.replaceAll("[cC][tlI][Ss]", "Cl3");
					}
					if(val.matches(".*[cC][tlI][8].*")){
						val=val.replaceAll("[cC][tlI][8]", "Cl3");
					}
					
					
					BranchNode bn = BranchNode.interpretOCRStringAsAtom2(val);
					
					//sometimes the letters it chooses are kinda weird ... and it should know better based on context
					//for example OO
					
					
					if(val.length()>5){
						if(bn==null){
							removeBad=true;
						}
					}
					
					if(parent.getBounds2D().getHeight()<=averageHeightNumberOCR*1.1 && areLikelyNumbers){
						//might be number
						if(val.matches("[0-9t][0-9t]*") ){
							val="#" + val.replace("t", "1");
							removeBad=true;
						}
					}
					
					if(bn==null || !bn.isRealNode()){
						if(REMOVE_NONSENSE_OCR_LINES){
							removeBad=true;
						}
					}
					
					if(removeBad){
						contains.stream()
					      //.filter(s2->likelyOCR.contains(s2))
					      .forEach(s2->{
					    	  likelyOCR.remove(s2);
					    	  likelyOCRNumbers.remove(s2);
					    	  likelyOCRNonBond.remove(s2);
					    	  likelyOCRAll.remove(s2);
					    	  
					      });
						foundNewOCR[0]=true;
						likelyOCRIgnore.add(GeomUtil.growShape(parent,2));
						
						//we need to clear up things that might have been found by accident due to bad
						//OCR
						likelyOCR.removeAll(ocrRescues);
						likelyOCRNumbers.removeAll(ocrRescues);
						likelyOCRNonBond.removeAll(ocrRescues);
						likelyOCRAll.removeAll(ocrRescues);
					}
				
					bestGuessOCR.put(parent, val);
				}
			});
			

			
			
			//if(foundNewOCR[0] && repeats<MAX_OCR_FULL_REPEATS)continue;
			
			bestGuessOCR.entrySet()
						.stream()
						.map(Tuple::of)
						.filter(t->t.v().equals("H"))
						.collect(Collectors.toList())
						.forEach(t->{
							double cutoff=Math.max(ctab.getAverageBondLength()*MAX_BOND_RATIO_FOR_OCR_CHAR_SPACING,averageWidthOCR);
							
							Tuple<Shape,String> toConnect=bestGuessOCR.entrySet()
							            .stream()
							            .map(Tuple::of)
							            .filter(t1->t1.k()!=t.k())
							            .filter(t1->t1.v().equals("N") || t1.v().equals("Nt")  || 
							            		t1.v().equals("N+") || t1.v().equals("NI") || 
							            		t1.v().equals("Nl") || t1.v().equals("O") || 
							            		t1.v().equals("S"))
							            .filter(t1->GeomUtil.distance(t.k(), t1.k())<cutoff)
							            .filter(t1->(Math.abs(t1.k().getBounds2D().getMinX()-t.k().getBounds2D().getMinX())< cutoff/3.0))
							            .filter(t1->(Math.abs(t1.k().getBounds2D().getMaxX()-t.k().getBounds2D().getMaxX())< cutoff/1.5))
							            .findFirst()
							            .orElse(null);
							            ;
							            
							 if(toConnect!=null){
								 Shape nshape=GeomUtil.add(t.k(),toConnect.k());
								 String old=toConnect.v();
								 if(old.contains("N") && !old.equals("N+"))old="N";
								 String nstring=old + t.v();
								 bestGuessOCR.put(nshape, nstring);
								 bestGuessOCR.remove(t.k());
								 bestGuessOCR.remove(toConnect.k());
							 }
						});
			
			bestGuessOCR.entrySet()
			.stream()
			.map(Tuple::of)
			.filter(t->t.v().equals("O2") || t.v().equalsIgnoreCase("Boc"))
			.collect(Collectors.toList())
			.stream()
			.forEach(t->{
				double cutoff=Math.max(ctab.getAverageBondLength()*MAX_BOND_RATIO_FOR_OCR_CHAR_SPACING,averageWidthOCR);
				
				Tuple<Shape,String> toConnect=bestGuessOCR.entrySet()
				            .stream()
				            .map(Tuple::of)
				            .filter(t1->t1.k()!=t.k())
				            .filter(t1->t1.v().equals("S") || t1.v().equals("N"))
				            .filter(t1->GeomUtil.distance(t.k(), t1.k())<cutoff)
				            .filter(t1->(Math.abs(t1.k().getBounds2D().getMinX()-t.k().getBounds2D().getMinX())< cutoff/3.0))
				            .findFirst()
				            .orElse(null);
				            ;
				            
				 if(toConnect!=null){
					 Shape nshape=GeomUtil.add(t.k(),toConnect.k());
					 String old=toConnect.v();
					 String nstring=old + t.v();
					 bestGuessOCR.put(nshape, nstring);
					 bestGuessOCR.remove(t.k());
					 bestGuessOCR.remove(toConnect.k());
				 }
			});
			
			
			
			

			ctab.standardCleanEdges();
			
			


			List<Shape> ocrMeaningful=bestGuessOCR.keySet()
					.stream()
					.filter(s->BranchNode.interpretOCRStringAsAtom2(bestGuessOCR.get(s))!=null)
					.collect(Collectors.toList());


			System.out.println("Assigning OCR:" + ctabRaw.size());
			if(DEBUG)ctabRaw.add(ctab.cloneTab());
			//ctab.removeOrphanNodes();

			Set<Node> alreadyFixedNodes = new HashSet<Node>();

			Function<Shape, Double> priorityNumber = (s)->{
				Rectangle2D rect = s.getBounds2D();
				return rect.getWidth()*rect.getHeight() + (1+rect.getX())/1000.0 + (1+rect.getY())/1000000.0; 
			};
			
			bestGuessOCR.entrySet()
						.stream()
						.map(Tuple::of)
						.map(t->Tuple.of(t.k(),
								         Tuple.of(t.v(),(t.v().equals("H"))?1.0:0.0 + 1.0/(priorityNumber.apply(t.k())))
								              .withVComparator()))
						.map(t->t.withVComparator())
						.sorted()
						.filter(t->!t.v().k().startsWith("#")) //??? if it isn't a number, this should never be the case now, but might be in the future
						.map(Tuple.vmap(t->t.k()))
						.collect(Collectors.toList())
						.forEach(shapeString->{
							Shape s= shapeString.k();
							String sym = shapeString.v();
							
							Point2D cen = GeomUtil.findCenterOfShape(s);
							Point2D centert = cen;
			
							BranchNode actual1;
							if(sym.equals("I")){
								Shape gs=GeomUtil.growShapeNPoly(s, 2, 10);
								//realRescueOCRCandidates.add(gs);
								List<Node> ln=ctab.getNodesInsideShape(gs, 0.2);
								
								if(ln.isEmpty())return;
								boolean isLinker=ln.stream().allMatch(n->n.getEdgeCount()>1);
								if(isLinker){
									
									return;
								}
								List<Node> tooClose=ctab.getNodes()
								    .stream()
								    .filter(n->!ln.contains(n))
								    .filter(n->n.getPoint().distance(cen)<ctab.getAverageBondLength()*0.8)
								    .collect(Collectors.toList());
								if(!tooClose.isEmpty()){
									if(tooClose.size()>1)return;
									Node tooCloseN = tooClose.get(0);
									if(tooCloseN.getEdgeCount()!=2)return;
									double smallestX = gs.getBounds2D().getMinX();
									double largestX = gs.getBounds2D().getMaxX();
									
									if(tooCloseN.getPoint().getX()>=smallestX && tooCloseN.getPoint().getX()<=largestX){
										boolean works=tooCloseN.getNeighborNodes().stream()
																 .map(t->t.k().getPoint().getX())
														         .allMatch(xx->xx>=smallestX && xx<=largestX);
										if(!works)return;
									}else{
										return;
									}
									
								}
								boolean doublePoss=rejBondOrderLines.stream()
												  .filter(lw->lw.growLine(ctab.getAverageBondLength()*0.1).contains(cen))
												  .findAny()
												  .isPresent();
								
								if(doublePoss)return;
								
								
								
								if(s.getBounds2D().getHeight()>Math.max(averageHeightOCR,ctab.getAverageBondLength())*1.1){
									//realRescueOCRCandidates.add(gs);
									return;
								}else{
									Line2D ll = new Line2D.Double(cen.getX(),cen.getY()-ctab.getAverageBondLength(),
																  cen.getX(),cen.getY()+ctab.getAverageBondLength());
									Shape llw = LineWrapper.of(ll).growLine(1);
									
									//realRescueOCRCandidates.add(ll);
									boolean inter=likelyOCRNonBond.stream()
									                .filter(ss->GeomUtil.intersects(ss, llw))
									                .findAny()
									                .isPresent();
									if(inter)return;
									
									
								}
								actual1=new BranchNode("I");
								//bestGuessOCR.remove(s);
								bestGuessOCR.put(s, "t"); //t is used for 'I' for bad legacy reasons I'm too lazy to fix
							}else{
								actual1 = BranchNode.interpretOCRStringAsAtom2(sym);	
							}
							BranchNode actual = actual1;
							
							
							if(actual!=null && actual.isRealNode()){
								
								if(sym.length()>1){
									List<Line2D> externalLines=ctab.getAllEdgesEntering(s, MAX_BOND_RATIO_FOR_LINES_CONSIDERED_FOR_POSITIONING_OCR*ctab.getAverageBondLength())
											.stream()
											.filter(t->!bestGuessOCR.keySet()
													               .stream()
													               .filter(ss->ss!=s)
													               .filter(ss->ss.contains(t.v().k().getPoint()))
													               .findAny()
													               .isPresent())
//											.filter(t->!bestGuessOCR.keySet()
//										               .stream()
//										               .filter(ss->ss!=s)
//										               .filter(ss->ss.contains(t.v().v().getPoint()))
//										               .findAny()
//										               .isPresent())
											.map(t->t.k().getLine())
											.collect(Collectors.toList());
									if(externalLines.size()==1){
										Line2D exl=externalLines.get(0);
										Point2D tc=centert;
										Point2D cnew=likelyOCR.stream()
												.map(s1->GeomUtil.findCenterOfShape(s1))
												.filter(spt->s.contains(spt))
												.map(cpt->GeomUtil.projectPointOntoLineWithRejection(exl, cpt))
												.map(Tuple.vmap(d->Math.abs(d)))
												.map(t->t.withVComparator())
												.min(Comparator.naturalOrder())
												.map(t->t.k())
												.orElseGet(()->{
													return GeomUtil.projectPointOntoLine(externalLines.get(0), tc);
												});
										if(s.contains(cnew)){
											centert=cnew;
			
										}
									}else{
										List<Point2D> intersections =GeomUtil.eachCombination(externalLines)
												.map(t->GeomUtil.intersection(t.k(),t.v()))
												.filter(p->p!=null)
												.filter(p->s.contains(p))
												.collect(Collectors.toList());
										if(intersections.size()==1){
											centert=intersections.get(0);
										}else if(intersections.size()>1){
											centert=GeomUtil.findCenterOfVertices(intersections);
										}
									}
									if(!s.contains(centert)){
										centert=GeomUtil.findCenterOfShape(s);
									}
								}
								Point2D center = centert;
								
								
								
								//This is likely the source of lots of problems
								ctab.mergeAllNodesInside(s, MAX_BOND_RATIO_FOR_MERGING_TO_OCR*ctab.getAverageBondLength(),(n)->{
									if(sym.equals("H")){
										Optional<Tuple<Shape, Double>> findClosestShapeTo = GeomUtil.findClosestShapeTo(ocrMeaningful, n.getPoint());
										if(!findClosestShapeTo.isPresent() || findClosestShapeTo.get().k() !=s){
											return false;
										}
									}
									boolean contains=s.contains(n.getPoint());
									
									if(!contains){
										boolean bb=bestGuessOCR.keySet()
												               .stream()
										                       .filter(ss->ss.contains(n.getPoint()))
										                       .findAny()
										                       .isPresent();
										if(bb)return false;
										
										if(actual.isTerminal() && n.getEdgeCount()>1){
											//System.out.println("Term?");
											long cc=n.getNeighborNodes()
											 .stream()
											 .map(t->t.k())
											 .filter(nn->GeomUtil.distanceTo(s, nn.getPoint())<2)
											 .count();
											
											if(cc==0)return false;
										}
										if(n.getEdgeCount()>1){
											boolean nhas=n.getNeighborNodes().stream()
											                    .map(n1->n1.k().getPoint())
											                    .filter(p->s.contains(p))
											                    .findAny()
											                    .isPresent();
											if(!nhas){
												boolean edgeHas=n.getEdges()
																 .stream()
																 .map(e->e.getLine())
																 .map(l->GeomUtil.getLineInside(l, s))
																 .filter(l->l.isPresent())
																 .findAny()
																 .isPresent();
												if(!edgeHas){
													
													double avgNDist=GeomUtil.eachCombination(n.getNeighborNodes()
															                  .stream()
																			  .map(n1->n1.k().getPoint())
																			  .collect(Collectors.toList()))
													        .map(t->new Line2D.Double(t.k(), t.v()))
													        .mapToDouble(l->GeomUtil.length(l))
													        .average()
													        .orElseGet(()->ctab.getAverageBondLength());
													if(avgNDist>ctab.getAverageBondLength()){
														if(!sym.equals("H")){
															Tuple<Shape,Double> bs=GeomUtil.findClosestShapeTo(likelyOCR, n.getPoint()).orElse(null);
															if(bs!=null){
																char tt=Optional.ofNullable(ocrAttempt.get(bs.k()))
																        .map(l->l.get(0).k())
																        .orElse(null);
																if(tt=='H'){
																	return false;
																}
																		
															}
														}
													}
													
												}
											}
										}
									}
									
									
									if(alreadyFixedNodes.contains(n))return false;
									if(n.getEdgeCount()==0)return false;
									return true;
								},(l)->{
			
									boolean matchesOthers=l.stream()
											.map(pt->GeomUtil.findClosestShapeTo(ocrMeaningful, pt))
											.filter(Optional::isPresent)
											.map(o-> o.get().k())
											
											.filter(sb->(sb!=s))
											.findAny()
											.isPresent();
									if(!matchesOthers){
										return center;
									}else{
										//return center;
										return GeomUtil.findCenterMostPoint(l);
									}
			
			
			
								});
								List<Node> mergedNodes=ctab.getAllNodesInsideShape(s,MAX_BOND_RATIO_FOR_MERGING_TO_OCR*ctab.getAverageBondLength());
								alreadyFixedNodes.addAll(mergedNodes);
							}
						});

			ctab.standardCleanEdges();


			if(DEBUG)ctabRaw.add(ctab.cloneTab());

			


			List<Tuple<Line2D,Point2D>> lj =Stream.concat(removedTinyLines.stream(), linesJoined.stream().map(l->l.getLine()))
								  // 
								   .flatMap(l->GeomUtil.getLinesNotInside(l, growLikelyOCRNonBond).stream())
								   //.filter(l->GeomUtil.length(l)>2) // needed?
								   .map(l->Tuple.of(l,GeomUtil.findCenterOfShape(l)))
								   
								   .collect(Collectors.toList());

			//lines=lj;
			


			Set<Line2D> taken = new HashSet<Line2D>();


			List<Tuple<List<Line2D>,Tuple<Node,Node>>> edgesToMake = new ArrayList<>();

			GeomUtil.groupShapesIfClosestPointsMatchCriteria(likelyOCR, (t)->{
				Point2D[] pts=t.v();
				if(pts[0].distance(pts[1])<ctab.getAverageBondLength()*.9){
					return true;
				}
				return false;
			})
			.stream()
			.filter(ls->ls.size()>=2)
			.flatMap(ls->{


				return GeomUtil.eachCombination(ls)
						.filter(t->GeomUtil.distance(t.k(), t.v())<ctab.getAverageBondLength()*1.2)
						.map(t->{
							return Tuple.of(t,GeomUtil.add(GeomUtil.growShape(t.k(),1), GeomUtil.growShape(t.v(),1)));
						})
						.map(Tuple.vmap(s->Tuple.of(s,GeomUtil.area(s)).withVComparator()))
						.map(t->t.withVComparator())
						.sorted()
						.map(Tuple.vmap(st->st.k()))
						.map(t->{
							//It should also have at least 1 line segment between the two
							Shape cshape = t.v();

							Line2D makeLine = new Line2D.Double(GeomUtil.findCenterOfShape(t.k().k()), 
									GeomUtil.findCenterOfShape(t.k().v()));
							List<Line2D> opl=lj.stream()
									.filter(l1->cshape.contains(l1.v()))
									.map(t1->t1.k())
									.filter(l->GeomUtil.cosTheta(l, makeLine) > Math.cos(25*Math.PI/180))
									//.map(l1->Tuple.of(l1,taken.add(l1)))
									//.map(l1->Tuple.of(l1,true))
									.collect(Collectors.toList());

							return Tuple.of(t.k(),opl);
						});

			})
			.forEach(lst->{

				List<Node> nodes = Stream.of(lst.k().k(),lst.k().v())
						.map(s->ctab.getNodesInsideShape(s, 2))
						.flatMap(nds->nds.stream())
						.distinct()
						.collect(Collectors.toList());
				if(nodes.size()==2){
					Node n1=nodes.get(0);
					Node n2=nodes.get(1);
					Edge alreadyEdge=ctab.getEdgeBetweenNodes(n1, n2).orElse(null);
					boolean haspossibleLine = !lst.v().isEmpty();


					boolean already=(alreadyEdge!=null);
					


					if(!already && haspossibleLine){
						//sometimes adds wrong bonds
						edgesToMake.add(Tuple.of(lst.v(),Tuple.of(nodes.get(0),nodes.get(1))));
					}else if(!already && !haspossibleLine){
						//do nothing
					}else if(already && !haspossibleLine){
						ctab.removeEdge(alreadyEdge);
					}else{
						
						taken.addAll(lst.v());

						int order=GeomUtil.groupThings(lst.v(), tlines->{
							Line2D l1=tlines.k(); 
							Line2D l2=tlines.v();
							if(GeomUtil.cosTheta(l1, l2) > Math.cos(10*Math.PI/180)){
								return true;
							}
							return false;
						})
								.stream()
								.map(l->GeomUtil.getLineOffsetsToLongestLine(l))
								.map(l->{
									OptionalDouble opdoff=l.stream()
												 .filter(t->Math.abs(t.v())>ctab.getAverageBondLength()*0.1)
												 .mapToDouble(t->t.v())
												 .min();
									if(!opdoff.isPresent()){
										return l.stream().map(l1->l1.k()).limit(1).collect(Collectors.toList());
									}
									double doff=opdoff.getAsDouble();
									
									List<Line2D> nlines= l.stream()
										 .map(Tuple.vmap(d->(int)Math.round(d/doff)))
										 .map(t->t.swap())
										 .collect(Tuple.toGroupedMap())
										 .entrySet()
										 .stream()
										 .map(Tuple::of)
										 .map(Tuple.vmap(v1->GeomUtil.getPairOfFarthestPoints(GeomUtil.vertices(v1))))
										 .map(Tuple.vmap(v1->new Line2D.Double(v1[0], v1[1])))
										 .map(t->t.v())
										 .collect(Collectors.toList());
									
									return nlines;
								})
								.mapToInt(ll->ll.size())
								.max().getAsInt();

						if(order>alreadyEdge.getOrder()){
							alreadyEdge.setOrder(order);
						}
					}
				}
			});


			edgesToMake.stream()
						.filter(t->{
							return taken.addAll(t.k());   
						})
						.forEach(t->{
			
							List<Edge> crossingEdges=ctab.getBondsThatCross(t.v().k(),t.v().v());
							if(!crossingEdges.isEmpty())return;
			
							Line2D nline = new Line2D.Double(t.v().k().getPoint(),t.v().v().getPoint());
							
							List<Line2D> k = t.k().stream().filter(l->GeomUtil.cosTheta(nline, l)> Math.cos(45*Math.PI/180)).collect(Collectors.toList());
							
							if(k.isEmpty())return;
							//don't add cross bonds
			
							int order=GeomUtil.groupThings(k, tlines->{
								Line2D l1=tlines.k(); 
								Line2D l2=tlines.v();
								if(GeomUtil.cosTheta(l1, l2) > Math.cos(10*Math.PI/180)){
									return true;
								}
								return false;
							})
									.stream()
									.map(l->GeomUtil.getLineOffsetsToLongestLine(l))
									.map(l->{
										OptionalDouble opdoff=l.stream()
													 .filter(tb->Math.abs(tb.v())>ctab.getAverageBondLength()*0.1)
													 .mapToDouble(tb->tb.v())
													 .min();
										if(!opdoff.isPresent()){
											return l.stream().map(l1->l1.k()).limit(1).collect(Collectors.toList());
										}
										double doff=opdoff.getAsDouble();
										
										List<Line2D> nlines= l.stream()
											 .map(Tuple.vmap(d->(int)Math.round(d/doff)))
											 .map(tb->tb.swap())
											 .collect(Tuple.toGroupedMap())
											 .entrySet()
											 .stream()
											 .map(Tuple::of)
											 .map(Tuple.vmap(v1->GeomUtil.getPairOfFarthestPoints(GeomUtil.vertices(v1))))
											 .map(Tuple.vmap(v1->new Line2D.Double(v1[0], v1[1])))
											 .map(tb->tb.v())
											 .collect(Collectors.toList());
										
										return nlines;
									})
									.mapToInt(ll->ll.size())
									.max().getAsInt();
							ctab.addEdge(t.v().k().getIndex(), t.v().v().getIndex(), order);
						});

			ctab.standardCleanEdges();

			

			if(DEBUG)ctabRaw.add(ctab.cloneTab());

			ctab.makeMissingBondsToNeighbors(bitmap,MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL,MAX_TOLERANCE_FOR_SINGLE_BONDS,likelyOCR,OCR_TO_BOND_MAX_DISTANCE, (t)->{

				if(t.k()>MAX_TOLERANCE_FOR_SINGLE_BONDS){
					t.v().setDashed(true);
				}
			});

			
			if(DEBUG)ctabRaw.add(ctab.cloneTab());
			//remove triangles that are obviously wrong
			{
				double avgL = ctab.getAverageBondLength();
				Set<Edge> skip= new HashSet<Edge>();
				ctab.getEdges()
					.stream()
					.map(l->Tuple.of(l,LineWrapper.of(l.getLine())).withVComparator())
					.filter(e->e.v().length()>avgL)
					.sorted(Comparator.reverseOrder())
					.map(t->t.k())
					.filter(t->!skip.contains(t))
					.forEach(e->{
						Node n1= e.getRealNode1();
						Node n2= e.getRealNode2();
						List<Tuple.KEqualityTuple<Node,Edge>> neigh1=n1.getNeighborNodes();
						List<Tuple.KEqualityTuple<Node,Edge>> neigh2=n2.getNeighborNodes();
						List<Tuple.KEqualityTuple<Node,Edge>> things1=neigh1.stream()
								.filter(ne->neigh2.contains(ne))
								.collect(Collectors.toList());
						List<Tuple.KEqualityTuple<Node,Edge>> things2=neigh2.stream()
								.filter(ne->neigh1.contains(ne))
								.collect(Collectors.toList());
						List<Tuple.KEqualityTuple<Node,Edge>> things = Stream.concat(things1.stream(), things2.stream())
								.collect(Collectors.toList());
	
						if(things.size()>0){
							
							Point2D p1=n1.getPoint();
							Point2D p2=n2.getPoint();
							Node n3=things.get(0).k();
							
							Point2D p3=n3.getPoint();
							
							if(n3.getEdgeCount()==4 && intersectionNodes.stream().filter(in->in.distance(p3)<2).findAny().isPresent()){
								things.stream()
									.map(t->t.v())
									.forEach(e2->{
										skip.add(e2);
									});
								return;
							}
							
							if(n1.getEdgeCount()==4 && intersectionNodes.stream().filter(in->in.distance(p1)<2).findAny().isPresent()){
								things.stream()
									.map(t->t.v())
									.forEach(e2->{
										skip.add(e2);
									});
								return;
							}
							if(n2.getEdgeCount()==4 && intersectionNodes.stream().filter(in->in.distance(p2)<2).findAny().isPresent()){
								things.stream()
									.map(t->t.v())
									.forEach(e2->{
										skip.add(e2);
									});
								return;
							}
							
							
	
							Edge oedge1=things.get(0).v();
							Edge oedge2=things.get(1).v();
	
							double tarea=Math.abs(GeomUtil.areaTriangle(p1,p2,p3));
	
							double expected = Math.sqrt(3)/4*Math.pow(e.getEdgeLength(),2);
							if(tarea<expected*0.5){
								
								boolean removeLong=true;
	
								if((e.getEdgeLength()<avgL*1.8) &&
										(!e.getDashed() && (oedge1.getDashed() && oedge2.getDashed()) ||
										(oedge1.getDashed() || oedge2.getDashed() && e.getOrder()>1))
										){
									removeLong=false;
								}else{
									if(oedge1.getEdgeLength()<avgL*0.7 && oedge2.getEdgeLength()<avgL*0.7){
										removeLong=false;
									}else{
										removeLong=true;
									}
								}
								
								if(removeLong){
									ctab.removeEdge(e);
								}else{
									ctab.removeEdge(oedge1);
									ctab.removeEdge(oedge2);
								}
							}else{
								things.stream()
								.map(t->t.v())
								.forEach(e2->{
									skip.add(e2);
								});
							}
						}
	
	
					});
				//for each edge, 

			}
			ctab.removeOrphanNodes();

			if(DEBUG)ctabRaw.add(ctab.cloneTab());

			//final cleanup
			{
				ctab.getDashLikeScoreForAllEdges(bitmap,likelyOCR)
				.forEach(t->{
					if(t.v()<MIN_ST_DEV_FOR_KEEPING_DASHED_LINES && t.k().getDashed()){
						t.k().setDashed(false);
						//Maybe it shouldn't even be here?
						//Try to remove it
						double tol=ctab.getToleranceForEdge(t.k(),bitmap,likelyOCR);
						if(tol>MAX_TOLERANCE_FOR_DASH_BONDS){
							ctab.removeEdge(t.k());
						}
					}
				});

			}
			
			System.out.println("cleaning triples:" + ctabRaw.size());
			if(DEBUG)ctabRaw.add(ctab.cloneTab());
			double fbondlength=ctab.getAverageBondLength();

			

			
			//clean bad triple bonds
			ctab.getEdges().stream()
					.filter(e->e.getOrder()==3)
					.collect(Collectors.toList())
					.forEach(e->{
						LineWrapper lb = LineWrapper.of(e.getLine());
						
						double len = e.getEdgeLength();
						
						int c=ctab.getEdges().size();
						
						double otherBondAverage = (c*ctab.getAverageBondLength()-len) 
																/
															  (c-1);
						
						int n = (int)Math.round(len/otherBondAverage);
						
						if(n>1){
							Node n1=e.getRealNode1();
							Node n2=e.getRealNode2();
							
							Shape bigLineShape = lb.growLine(len/3);
							
							Point2D apt=rejBondOrderLines.stream()
											 .filter(l->l.length()>otherBondAverage*0.4)
							                 .filter(p->bigLineShape.contains(p.centerPoint()))
							                 .map(p->lb.projectPointOntoLine(p.centerPoint()))
							                 .collect(GeomUtil.averagePoint());
							
							
							List<Point2D> pts=lb.splitIntoNPieces(n);
							
							Node pnode=n1;
							Edge closestEdge = null;
							double closestD = Double.POSITIVE_INFINITY;
							
							for(int i=1;i<pts.size();i++){
								
								Node nn=null;
								if(i<pts.size()-1){
									Point2D np=pts.get(i);
									nn=ctab.addNode(np);	
								}else{
									nn=n2;
								}
								
								Edge ne=ctab.addEdge(pnode.getIndex(),nn.getIndex(),1);
								Point2D cpt = Stream.of(ne.getPoint1(),ne.getPoint2()).collect(GeomUtil.averagePoint());
								double dpt=cpt.distance(apt);
								if(dpt<closestD){
									closestD=dpt;
									closestEdge=ne;
								}
								pnode=nn;
							}
							if(closestEdge!=null){
								
								closestEdge.setOrder(3);
								ctab.removeEdge(e);
							}
							
						}
					});
			
			if(DEBUG)ctabRaw.add(ctab.cloneTab());
			
			List<Shape> appliedOCR = new ArrayList<Shape>();

			for(Shape s: bestGuessOCR.keySet()){
				String sym=bestGuessOCR.get(s);
				
				BranchNode actual=BranchNode.interpretOCRStringAsAtom2(sym);
				if(actual!=null && actual.isRealNode()){
					
					//This means it's a free floating group, don't connect it. Probably
					//don't even keep it for now
					if(!actual.isLinkable()){
						for(Node n:ctab.getAllNodesInsideShape(s, 0.1)){
							ctab.removeNodeAndEdges(n);	
						}
						
						//eventually add it back
						//ctab.addNode(p)
						continue;
					}
					
					appliedOCR.add(s);
					
					List<Node> nlist=ctab.getNodesInsideShape(s, 0.1).stream().filter(n->!n.isInvented()).collect(Collectors.toList());
					if(actual.getSymbol().equals("I") && !nlist.isEmpty() && nlist.get(0).getEdgeCount()>1){
					
						continue;
					}
					
					
					//Didn't merge, maybe merge now?
					if(nlist.size()>1){
						Point2D np=nlist.get(0).getPoint();
						ctab.mergeAllNodesInside(s, 0.1, (n)->!n.isInvented(), (l)->np);
						ctab.standardCleanEdges();
						nlist=ctab.getNodesInsideShape(s, 0.1).stream().filter(n->!n.isInvented()).collect(Collectors.toList());
					}
					nlist.forEach(n->{
						n.setSymbol(actual.getSymbol());
					});
					
					//ctab.setNodeToSymbol(s, actual.getSymbol());

					if(nlist.size()==1){
						Node pnode=nlist.get(0);
						pnode.setCharge(actual.getCharge());
						Point2D ppoint=pnode.getPoint();
						
						Node lneigh=null;
						Node rneigh=null;
						
						
						if(pnode.getEdgeCount()>=2){
							//non terminal
							lneigh=pnode.getNeighborNodes().stream().map(n->n.k())
									//.filter(n->n.getPoint().getX()<pnode.getPoint().getX())
									.map(n->Tuple.of(n,n.getPoint().getX()))
									.map(t->t.withVComparator())
									.min(Comparator.naturalOrder())
									.map(t->t.k())
									.orElse(null);
							rneigh=pnode.getNeighborNodes().stream().map(n->n.k())
									//.filter(n->n.getPoint().getX()>pnode.getPoint().getX())
									.map(n->Tuple.of(n,n.getPoint().getX()))
									.map(t->t.withVComparator())
									.max(Comparator.naturalOrder())
									.map(t->t.k())
									.orElse(null);
						}
						
						if(actual.hasChildren()){
							actual.generateCoordinates();
							AffineTransform at = new AffineTransform();
							at.translate(ppoint.getX(), ppoint.getY());
							at.scale(fbondlength, fbondlength);
							if(pnode.getEdges().size()>0){
								Point2D otherPoint =pnode.getEdges().stream().map(e->e.getOtherNode(pnode)).map(n->n.getPoint()).collect(GeomUtil.averagePoint());
								double ang=GeomUtil.angle(ppoint, otherPoint);
								at.rotate(ang+Math.PI);
							}
							actual.applyTransform(at);

							Map<BranchNode,Node> parentNodes = new HashMap<BranchNode,Node>();

							parentNodes.put(actual, pnode);

							actual.forEachBranchNode((parN,curN)->{
								if(parN==null) return;
								Node mpnode=pnode;
								if(parN!=null){
									mpnode=parentNodes.get(parN);
								}

								Node n= ctab.addNode(curN.getSuggestedPoint())
										    .setSymbol(curN.getSymbol())
										    .setCharge(curN.getCharge())
										    .setInvented(true);
								
								Edge e=ctab.addEdge(mpnode.getIndex(), n.getIndex(), curN.getOrderToParent());
								if(curN.getWedgeType()==1){
									 e.setWedge(true);
								}else if(curN.getWedgeType()==-1){
									 e.setDashed(true);
								}
								
								curN.getRingBond().ifPresent(t->{
									Node rn=parentNodes.get(t.k());
									ctab.addEdge(rn.getIndex(), n.getIndex(), t.v());
								});
								
								parentNodes.put(curN, n);
							});
							
							//might be a chain
							if(actual.canBeChain()){
								if(lneigh!=null&&rneigh!=null){
									Node l=lneigh;
									Node r=rneigh;
									Node nl=parentNodes.get(actual.getLeftBranchNode());
									Node nr=parentNodes.get(actual.getRightBranchNode());
									pnode.getNeighborNodes().stream()
															.filter(t->t.k()==l || t.k()==r)
															.collect(Collectors.toList())
									                        .forEach(t->{
									                        	Edge e=t.v();
									                        	Node on = t.k();
									                        	Node nn = (on==l)?nr:nl;
									                        	ctab.removeEdge(e);
									                        	ctab.addEdge(on.getIndex(), nn.getIndex(), e.getOrder());
									                        });		
									Line2D oldLine = new Line2D.Double(nl.getPoint(), nr.getPoint());
									Point2D[] far = GeomUtil.getPairOfFarthestPoints(s);
								
									double minx = Math.min(far[0].getX(),far[1].getX());
									double maxx = Math.max(far[0].getX(),far[1].getX());
									
									Line2D newLine = new Line2D.Double(maxx-averageWidthOCR/2,pnode.getPoint().getY(),minx+averageWidthOCR/2,pnode.getPoint().getY());
									
									
									Point2D navg=Stream.of(l.getPoint(),r.getPoint()).collect(GeomUtil.averagePoint());
									boolean flip=false;
									if(navg.getY()>newLine.getY1()){
										flip=true;
									}
									
									AffineTransform att=GeomUtil.getTransformFromLineToLine(oldLine, newLine, flip);
									
									parentNodes.values().stream().forEach(pn->{
										pn.setPoint(att.transform(pn.getPoint(),null));
									});
									
								}
							}
							
						}
					}
				}

			}
			
			if(DEBUG)ctabRaw.add(ctab.cloneTab());
			
			// see what a bond would look like if it were in a ring
			// To do this, first make a 5-membered ring, and then transform it to be on the bond
			// in 2 ways, then look for atoms near where they're supposed to be. Only do this for edges which
			// are not ring edges
			{
				Shape nshape=GeomUtil.convexHull2(GeomUtil.makeNPolyCenteredAt(new Point2D.Double(0,0), 5, 100));
				Line2D sl = GeomUtil.lines(nshape)[0];
				
				List<Tuple<Tuple<Node,Node>,Integer>> toAdd = new ArrayList<>();
				
				
				ctab.getEdges()
					.stream()
				    .filter(e->!e.isRingEdge())
				    .filter(e->!e.isInventedBond())
				    .forEach(e->{
				    	Line2D l = e.getLine();
				    	
				    	
				    	for(int i=0;i<2;i++){
					    	AffineTransform at1=GeomUtil.getTransformFromLineToLine(sl,l,i==0);
					    	Shape nnshape = at1.createTransformedShape(nshape);
					    	
					    	Point2D[] verts2 = GeomUtil.vertices(nnshape);
					    	
					    	
					    	List<Node> nnodes=Arrays.stream(verts2)
					    	      .map(v->{
					    	    	  return ctab.getNodes()
						    	    	  .stream()
						    	    	  .filter(n->!n.isInvented())
					    	    		  .map(n1->Tuple.of(n1,n1))
					    	    		  .map(Tuple.vmap(n->n.getPoint().distance(v)))
					    	    		  .map(t->t.withVComparator())
					    	    		  .min(Comparator.naturalOrder())
					    	    		  .orElse(null);
					    	      })
					    	      .filter(t->t!=null)
					    	      .filter(t->t.v()<ctab.getAverageBondLength()*0.05)
					    	      .map(t->t.k())
					    	      .collect(Collectors.toList());
					    	if(nnodes.size()==5){
					    		realRescueOCRCandidates.add(nnshape);
					    		
					    		Node[] nodes = nnodes.stream().toArray(s->new Node[s]);
					    		
					    		for(int j=0;j<nodes.length;j++){
					    			Node nextn = nodes[(j+1)%nodes.length];
					    			Node thisn = nodes[(j)%nodes.length];
					    			
					    			if(!thisn.connectsTo(nextn)){
					    				if(!thisn.getSymbol().equals("C") || !nextn.getSymbol().equals("C")){
					    					int o=0;
					    					LineWrapper nline = LineWrapper.of(new Line2D.Double(nextn.getPoint(), thisn.getPoint()));
					    					Point2D avg=Stream.of(nextn.getPoint(),thisn.getPoint()).collect(GeomUtil.averagePoint());
					    					//probably should make an edge
					    					o+=linesOrder.stream()
					    									 .map(Tuple.kmap(ll->LineWrapper.of(ll)))
					    					                 .filter(ll->ll.k().centerPoint().distance(avg)<nline.length()*0.3)
					    					                 .filter(ll->ll.k().absCosTheta(nline)>Math.cos(10*Math.PI/180.0))
					    					                 .mapToInt(t->t.v())
					    					                 .sum();
					    					if(o==0)o=1;
					    					toAdd.add(Tuple.of(Tuple.of(thisn,nextn),o));
					    				}
					    			}
					    		}
					    	}
				    	}
				    	
				    });
				
				toAdd.forEach(t->{
					ctab.addEdge(t.k().k().getIndex(), t.k().v().getIndex(), t.v());
				});
				
				ctab.standardCleanEdges();
			}
			
			


			if(true){
				ConnectionTable biggestSection = ctab.getDisconnectedComponents()
					    .stream()
					    .map(ct->Tuple.of(ct,GeomUtil.area(ct.getConvexHull())).withVComparator())
					    .max(Comparator.naturalOrder())
					    .map(t->t.k())
					    .orElse(null);
				if(biggestSection!=null){
					if(biggestSection.getAverageBondLength()<20){
							bitmap=new Bitmap.BitmapBuilder(bitmap).scale(2).gaussBlur(1).threshold(6).build();
							load(bitmap);
							return;
					}
				}
			}
			

			if(ctab.getDisconnectedComponents().size()>1){
				
				Tuple<ConnectionTable,Shape> ctshape = ctab.getDisconnectedComponents()
				    .stream()
				    .map(ct->Tuple.of(ct,ct.getConvexHull()))
				    .map(Tuple.vmap(h->Tuple.of(h,-GeomUtil.area(h)).withVComparator()))
				    .map(t->t.withVComparator())
				    .sorted()
				    .limit(1)
				    .map(Tuple.vmap(t->t.k()))
				    .findFirst()
				    .orElse(null);
				if(ctshape!=null){
					
					double bestbond=ctshape.k().getAverageBondLength();
					List<ConnectionTable> ctabs=ctab.getDisconnectedComponents();
					Shape crop=ctabs
					    .stream()
					    .filter(ct->ct.getAverageBondLength()<0.7*bestbond || ct.getAverageBondLength()>1.33*bestbond)
					    .map(ct->ct.getAreaAround(bestbond))
					    .filter(a->a!=null)
					    .collect(GeomUtil.union())
					    .orElse(null);
					
					Shape keep=ctabs
						    .stream()
						    .filter(ct->ct.getAverageBondLength()>=0.7*bestbond && ct.getAverageBondLength()<=1.33*bestbond)
						    .map(ct->ct.getAreaAround(bestbond))
						    .collect(GeomUtil.union())
						    .orElse(null);
					
					//rescueOCRShapes.add(crop);
					if(crop!=null && keep!=null){
						polygons.stream()
								.filter(s->!likelyOCR.contains(s))
								.map(s->Tuple.of(s,GeomUtil.findCenterOfShape(s)))
								.filter(t->!keep.contains(t.v()))
						        .filter(t->crop.contains(t.v()))
						        .map(t->t.k())
						        .forEach(p->{
						        	likelyOCR.remove(p);
									likelyOCRNumbers.remove(p);
									likelyOCRNonBond.remove(p);
									likelyOCRAll.remove(p);
									likelyOCRIgnore.add(GeomUtil.growShapeNPoly(p, 2, 16));
									//realRescueOCRCandidates.add(GeomUtil.growShapeNPoly(p, 2, 16));
									foundNewOCR[0]=true;
						        });
					}
				}
				//realRescueOCRCandidates.add(crop);
			}

			if(foundNewOCR[0] && repeats<MAX_OCR_FULL_REPEATS){
				System.out.println("Starting over...");
				continue;
			}
//			

			//Now get all the nodes with 2 edges which have shorter than average bond length,
			//and which are not in OCR, and where the sum of the distances of the bonds is within 95%
			//of the distance of the two neighbors. Those are nodes to remove
			Set<Node> toRemove = new HashSet<Node>();
			do{
				toRemove.clear();

				if(Thread.currentThread().isInterrupted()){
					throw new InterruptedException();
				}


				ctab.getNodesNotInShapes(appliedOCR, 0)
				.stream()
				.map(n->Tuple.of(n,n.getNeighborNodes()))
				.filter(t->t.v().size()==2)
				.filter(t->t.v().get(0).v().getEdgeLength()<ctab.getAverageBondLength())
				.filter(t->t.v().get(1).v().getEdgeLength()<ctab.getAverageBondLength())
				.filter(t->t.v().get(0).v().getOrder()==1)
				.filter(t->t.v().get(1).v().getOrder()==1)
				.collect(Collectors.toList())
				.forEach(t->{
					Node n1=t.k();
					Node n2=t.v().get(0).k();
					Node n3=t.v().get(1).k();
					if(toRemove.contains(n2) || toRemove.contains(n3))return;




					double d1=n1.distanceTo(n2) + n1.distanceTo(n3);
					double d2=n2.distanceTo(n3);
					if(d2/d1>.95){
						Line2D longLine=new Line2D.Double(n2.getPoint(),n3.getPoint());
						Point2D np=GeomUtil.projectPointOntoLine(longLine, n1.getPoint());
						if(np.distance(n1.getPoint())<3){
							
							
							ctab.addEdge(n2.getIndex(), n3.getIndex(), 1);
							toRemove.add(n1);
						}
					}
				});        
				toRemove.forEach(n->ctab.removeNodeAndEdges(n));
				ctab.standardCleanEdges();
				
				
				
				
				
				
			}while(!toRemove.isEmpty());


			if(DEBUG)ctabRaw.add(ctab.cloneTab());


			//Cleanup "duplicate" lines that are probably problems. 

			ctab.getNodes()
				.stream()
				.filter(n->n.getEdgeCount()>=2)
				.map(n->{
					return Tuple.of(n,GeomUtil.eachCombination(n.getEdges())
							.map(t->{
								if(t.k().getEdgeLength()>t.v().getEdgeLength()){
									t=t.swap();
								}
								return t;
							})
							.filter(t->t.v().getEdgeLength()>=ctab.getAverageBondLength())
							.collect(Collectors.toList()));
				})
				.filter(ed->!ed.v().isEmpty())
				.collect(Collectors.toList())
				.forEach(te->{
					Node n=te.k();
	
					te.v().forEach(t->{
						Node tnode=t.k().getOtherNode(n);
						Node otherNode=t.v().getOtherNode(n);
						Point2D ppnt=GeomUtil.projectPointOntoLine(t.v().getLine(), tnode.getPoint());
						if(ppnt.distance(tnode.getPoint())<0.1*ctab.getAverageBondLength()){
							double sd1=ppnt.distance(otherNode.getPoint());
							if(sd1<t.v().getEdgeLength()){
								//remove long bond
								//change point
								//add edge to long bond other node
	
								
								tnode.setPoint(ppnt);
								ctab.addEdge(tnode.getIndex(), t.v().getOtherNode(n).getIndex(), t.v().getOrder());
								ctab.removeEdge(t.v());
							}
						}
					});
				});

			ctab.standardCleanEdges();


			if(DEBUG)ctabRaw.add(ctab.cloneTab());


			ctab.getEdges()
				.stream()
				.filter(e->e.getEdgeLength()<ctab.getAverageBondLength()*0.55)
				.collect(Collectors.toList())
				.forEach(e->{
					//look at small bonds
					//Maybe merge the atoms?
					Node n1=e.getRealNode1();
					Node n2=e.getRealNode2();
					
					
					boolean wasIntersection =intersectionNodes.stream()
					                .filter(p->p.distance(n1.getPoint())<ctab.getAverageBondLength()*0.2 || p.distance(n2.getPoint())<ctab.getAverageBondLength()*0.2 )
					                .findAny()
					                .isPresent();
					
					if(wasIntersection){
						//maybe it's a cage structure. If it is, try removing the node and connecting the other edges
						return;
					}
					                
					
					
					List<Node> neigh1=n1.getEdges().stream()
							     .filter(e1->e1!=e)
					             .map(e1->e1.getOtherNode(n1))
					             .collect(Collectors.toList());
					List<Node> neigh2=n2.getEdges().stream()
						     .filter(e1->e1!=e)
				             .map(e1->e1.getOtherNode(n2))
				             .collect(Collectors.toList());
					
					//Don't bother trying to merge terminal groups
					if(neigh1.isEmpty()||neigh2.isEmpty())return;
					
					//what would be the bond lengths if they were merged to the first node?
					
					double avgBondLengthIfN1Merge = neigh2.stream().mapToDouble(nn->nn.distanceTo(n1)).average().orElse(0);
					double avgBondLengthIfN2Merge = neigh1.stream().mapToDouble(nn->nn.distanceTo(n2)).average().orElse(0);
					
					double maxBondLengthIfN1Merge = neigh2.stream().mapToDouble(nn->nn.distanceTo(n1)).max().orElse(0);
					double maxBondLengthIfN2Merge = neigh1.stream().mapToDouble(nn->nn.distanceTo(n2)).max().orElse(0);
					
					double minBondLengthIfN1Merge = neigh2.stream().mapToDouble(nn->nn.distanceTo(n1)).min().orElse(0);
					double minBondLengthIfN2Merge = neigh1.stream().mapToDouble(nn->nn.distanceTo(n2)).min().orElse(0);
					
					Predicate<Double> isReasonable = d->{
						return d< ctab.getAverageBondLength()*1.3 && d >ctab.getAverageBondLength()*0.7; 
					};
					
					boolean n1Merge=false;
					boolean n2Merge=false;
					
					if(isReasonable.test(avgBondLengthIfN1Merge) && isReasonable.test(maxBondLengthIfN1Merge) && isReasonable.test(minBondLengthIfN1Merge)){
						n1Merge=true;
					}
					if(isReasonable.test(avgBondLengthIfN2Merge) && isReasonable.test(maxBondLengthIfN2Merge) && isReasonable.test(minBondLengthIfN2Merge)){
						n2Merge=true;
					}
					
					if(!n1Merge && !n2Merge)return;
					
					if(n1Merge && !n2Merge){
						Point2D mp=n1.getPoint();
						ctab.mergeNodes(Stream.of(n1.getIndex(),n2.getIndex()).collect(Collectors.toList()), (l)->mp);
					}else if(n1Merge && !n2Merge){
						Point2D mp=n2.getPoint();
						ctab.mergeNodes(Stream.of(n1.getIndex(),n2.getIndex()).collect(Collectors.toList()), (l)->mp);
					}else{
						ctab.mergeNodesAverage(n1.getIndex(), n2.getIndex());
					}
				});
			
			ctab.standardCleanEdges();

			System.out.println("Cleaning cages:" + ctabRaw.size());
			if(DEBUG)ctabRaw.add(ctab.cloneTab());
			
			List<Node> toRemoveNodesCage = new ArrayList<>();
			
			ctab.getNodes()
			    .stream()
			    .filter(n->n.getEdgeCount()==4)
			    .filter(n->n.getSymbol().equals("C")) 
			    .filter(n->!n.isInvented())
			    .filter(n->intersectionNodes.stream().filter(p->p.distance(n.getPoint())<ctab.getAverageBondLength()*0.1).findAny().isPresent())
			    .forEach(n->{

    				
			    	//might be cage node / cross bond. To resolve this, you should see if
			    	//1. If node is in a ring
			    	//2. The variability in bond length will decrease if you split them
			    	//3. Removing the node wouldn't stop ALL neighbors from being in a ring
			    	//4. There is a pretty big gap somewhere (sufficient) 
			    	//5. The node is in a 4-membered ring, a 3-membered ring, and a 5 membered ring
			    	if(n.isInRing(7)){
			    		
			    		//probably cage then
			    		List<Node> neigh =n.getNeighborNodes().stream().map(t->t.k()).collect(Collectors.toList());
			    		
//			    		neigh.stream().map(n1->n1.getPoint()).map(p->GeomUtil.makeShapeAround(p, 10))
//			    		.forEach(s->{
//			    			realRescueOCRCandidates.add(s);
//			    		});
			    		
			    		//realRescueOCRCandidates
			    		
			    		if(neigh.stream().filter(nn->toRemoveNodesCage.contains(nn)).findAny().isPresent()){
			    			return;
			    		}
			    		List<List<Node>> pairs=GeomUtil.groupThings(neigh, (t1)->{
			    			Line2D nline = new Line2D.Double(t1.k().getPoint(),t1.v().getPoint());
			    			Point2D pp=GeomUtil.projectPointOntoLine(nline,n.getPoint());
			    			if(pp.distance(n.getPoint())< ctab.getAverageBondLength()*0.02){
			    				return true;
			    			}
			    			return false;
			    		});
			    		
			    		if(pairs.size()==2){
			    			List<Node> npair1=pairs.get(0);
			    			List<Node> npair2=pairs.get(1);
			    			if(npair1.size()==2 && npair2.size()==2){
			    				Line2D l1p = new Line2D.Double(npair1.get(0).getPoint(),npair1.get(1).getPoint());
			    				Line2D l2p = new Line2D.Double(npair2.get(0).getPoint(),npair2.get(1).getPoint());
			    				Shape l1s=GeomUtil.growLine(l1p, ctab.getAverageBondLength()*0.2);
			    				Shape l2s=GeomUtil.growLine(l2p, ctab.getAverageBondLength()*0.2);
			    				
			    				boolean l1Has=lj.stream()
			    				  .filter(t->l1s.contains(t.v()))
			    				  .filter(t->GeomUtil.cosTheta(t.k(), l1p)>0.6)
			    				  .map(l->l.k())
			    				  .map(l->GeomUtil.growLine(l, ctab.getAverageBondLength()*0.2))
			    				  .filter(ls->ls.contains(n.getPoint()))
			    				  .findFirst()
			    				  .isPresent();
			    				
			    				boolean l2Has=lj.stream()
					    				  .filter(t->l2s.contains(t.v()))
					    				  .filter(t->GeomUtil.cosTheta(t.k(), l2p)>0.6)
					    				  .map(l->l.k())
					    				  .map(l->GeomUtil.growLine(l, ctab.getAverageBondLength()*0.2))
					    				  .filter(ls->ls.contains(n.getPoint()))
					    				  .findFirst()
					    				  .isPresent();
			    				
			    				
			    				boolean doit=false;
			    				
			    				if((l1Has && !l2Has) || (l2Has && !l1Has)){
			    					doit=true;
			    				}else{		    				
				    				//Really do it, but first make sure that the average bond length makes sense
			    					//OR that it meets special ring criteria
			    					
			    					int[] rings = new int[7];
			    					
			    					n.getAllRings()
			    					 .forEach(r->{
			    						 rings[r.size()]++;
			    					 });
			    					
			    					if(rings[3]>0 && rings[4]>0 && rings[5]>0){
			    						doit=true;
			    					}else{
				    					
					    				
					    				double averageOtherBondLength=ctab.getEdges().stream()
					    				               .filter(e2->!e2.hasNode(n))
					    				               .mapToDouble(e2->e2.getEdgeLength())
					    				               .average()
					    				               .orElse(ctab.getAverageBondLength());
					    				
					    				double averageCurrentBondLength=n.getEdges().stream().mapToDouble(e->e.getEdgeLength()).average().orElse(0);
					    				
					    				double averageNewBondLength=0.5*(npair1.get(0).distanceTo(npair1.get(1)) +
					    											npair2.get(0).distanceTo(npair2.get(1)));
					    				
					    				if(Math.abs(averageNewBondLength-averageOtherBondLength)<Math.abs(averageCurrentBondLength-averageOtherBondLength)){
					    						doit=true;
					    				}
			    					}
			    				}
			    				
			    				if(doit){
			    					List<Edge> redges=n.getEdges();
			    					List<Node> nn = n.getNeighborNodes().stream().map(t->t.k()).collect(Collectors.toList());
			    					Edge nedge1=ctab.addEdge(npair1.get(0).getIndex(), npair1.get(1).getIndex(),1);
			    					Edge nedge2=ctab.addEdge(npair2.get(0).getIndex(), npair2.get(1).getIndex(),1);
				    				
				    				redges.forEach(e->ctab.removeEdge(e));
				    				//ctab.removeNode(n.getIndex());
				    				
				    				boolean cancel=false;
				    				if(!nn.stream().filter(nn1->nn1.isInRing(7)).findAny().isPresent()){
				    					cancel=true;
				    				}
				    				
				    				if(cancel){
				    					ctab.removeEdge(nedge1);
				    					ctab.removeEdge(nedge2);
				    				}else{
				    					toRemoveNodesCage.add(n);	
				    				}
				    				redges.forEach(e->ctab.addEdge(e.getRealNode1().getIndex(),e.getRealNode2().getIndex(),e.getOrder())
				    						.setDashed(e.getDashed())
				    						.setWedge(e.getWedge())
				    						);
				    				
				    				
				    				
			    				}
			    			}
			    		}
			    	}
			    });
			    
			for(Node r:toRemoveNodesCage){
				ctab.removeNodeAndEdges(r);
			}
			
			
			toRemoveNodesCage.clear();
			
			ctab.getNodes()
			    .stream()
			    .filter(n->n.getEdgeCount()==2)
			    .filter(n->n.getSymbol().equals("C"))
			    .filter(n->!n.isInvented())
			    .filter(n->n.getEdges().stream().filter(e->e.getOrder()==1).count()==2)
			    .filter(n->GeomUtil.findClosestShapeTo(likelyOCR, n.getPoint()).map(t->t.v()).orElse(100.0)>ctab.getAverageBondLength()*0.2)
			    //.filter(n->n.isInRing(8))
			    .forEach(n->{
			    	List<Tuple<Node,Node>> tn=GeomUtil.eachCombination(n.getNeighborNodes())
			    	        .filter(t->{
			    	        	Node n1=t.k().k();
			    	        	Node n2=t.v().k();
			    	        	
			    	        	Line2D l = new Line2D.Double(n1.getPoint(), n2.getPoint());
			    	        	Point2D pp=GeomUtil.projectPointOntoLine(l, n.getPoint());
			    	        	if(pp.distance(n.getPoint())<ctab.getAverageBondLength()*0.1){
			    	        		return true;
			    	        	}
			    	        	return false;
			    	        })
			    	        .map(t->Tuple.of(t.k().k(),t.v().k()))
			    	        .collect(Collectors.toList());
			    	if(tn.size()==1){
			    		toRemoveNodesCage.add(n);
			    		ctab.addEdge(tn.get(0).k().getIndex(), tn.get(0).v().getIndex(), 1);
			    	}
			    	//System.out.println("Cage like:" + tn.size());
			    	        
			    });
			if(DEBUG)ctabRaw.add(ctab.cloneTab());
			toRemoveNodesCage.forEach(n->{
				ctab.removeNodeAndEdges(n);	
			});
			ctab.standardCleanEdges();
			
			
			//System.out.println("MAde new nodes on triples:" + ctabRaw.size());
			if(DEBUG)ctabRaw.add(ctab.cloneTab());
			ctab.removeOrphanNodes();
			
			
			List<Shape> singleBondInfluenceAreas = ctab.getEdges()
					                                   .stream()
					                                   .filter(e->e.getOrder()==1)
					                                   .filter(e->!e.isInventedBond())
					                                   .map(e->Tuple.of(e,GeomUtil.growLine(e.getLine(), e.getEdgeLength()/5)))
					                                   .map(t->t.v())
					                                   .collect(Collectors.toList());
			Set<Edge> wasDouble = new HashSet<Edge>();
			
			//clean bad double bonds
			ctab.getEdges().stream()
				.filter(e->e.getOrder()==2)
				.filter(e->!e.isInventedBond())
				.filter(e->e.getRealNode1().getSymbol().equals("C") && e.getRealNode2().getSymbol().equals("C"))
				.peek(e->wasDouble.add(e))
				.collect(Collectors.toList())
				.forEach(e->{
						Line2D lb = e.getLine();
						Point2D apnt=Stream.of(lb.getP1(),lb.getP2()).collect(GeomUtil.averagePoint());
						
						List<Point2D> possibleOtherDoubleBonds = rejBondOrderLines.stream()
								 .filter(l->l.length()>ctab.getAverageBondLength()*0.4)
				                 .map(l->l.centerPoint())
				                 .filter(p->p.distance(apnt)<ctab.getAverageBondLength()*0.8)
				                 .collect(Collectors.toList());
						
						if(possibleOtherDoubleBonds.isEmpty()){
							//This makes me suspicious ... could be that it shouldn't have been flagged in the first place
							//If it's colinear, probably a mistake
							e.getNeighborEdges()
							 .stream()
							 .filter(e2->GeomUtil.cosTheta(e.getLine(),e2.getLine())>Math.cos(2*Math.PI/180))
							 .filter(e2->wasDouble.contains(e2))
							 .findAny()
							 .ifPresent(eo->{
								
								e.setOrder(1); 
							 });
							return;
						}
						
						boolean couldBeAnother=possibleOtherDoubleBonds.stream()
						                 .flatMap(p->singleBondInfluenceAreas.stream().filter(s->s.contains(p)))
						                 .findAny()
						                 .isPresent();
						                 //.collect(GeomUtil.averagePoint());
						if(couldBeAnother){
							e.setOrder(1);
						}
						//ctab.removeEdge(e);
				});
			
			//clean bad triple bonds
			ctab.getEdges().stream()
				.filter(e->e.getOrder()==3)
				.filter(e->!e.isInventedBond())
				.filter(e->e.getRealNode1().getSymbol().equals("C") && e.getRealNode2().getSymbol().equals("C"))
				.collect(Collectors.toList())
				.forEach(e->{
						
						LineWrapper lb = LineWrapper.of(e.getLine());
						Point2D apnt=lb.centerPoint();
						
						List<Tuple<Line2D,Point2D>> possibleOtherDoubleBonds = rejBondOrderLines.stream()
								 .filter(l->l.length()>ctab.getAverageBondLength()*0.4)
								 .filter(l->l.absCosTheta(lb)>0.8)
				                 .map(l->Tuple.of(l,l.centerPoint()))
				                 .filter(p->p.v().distance(apnt)<ctab.getAverageBondLength()*0.8)
				                 .map(Tuple.kmap(l->l.getLine()))
				                 .collect(Collectors.toList());
						
						if(possibleOtherDoubleBonds.size() ==1){
							//probably should have been a double bond
							e.setOrder(2);
							return;
						}
						List<Tuple<Edge,Shape>> edgeShapes = ctab.getEdges()
								                                 .stream()
								                                 .filter(e1->e1!=e)
								                                 .map(e1->Tuple.of(e1,GeomUtil.growLine(e1.getLine(), ctab.getAverageBondLength()*0.5)))
								                                 .collect(Collectors.toList());
						
						List<Tuple<Line2D,Edge>> bestEdge = possibleOtherDoubleBonds.stream()
						                        .map(t->Tuple.of(t.k(),edgeShapes.stream()
						                        		                     .filter(es->es.v().contains(t.v()))
						                        		                     .filter(es->GeomUtil.cosTheta(t.k(),es.k().getLine())>0.8)
						                        		                     .map(ee->ee.k())
						                        		                     .findAny()
						                        		                     .orElse(null)
						                        		))
						                        .collect(Collectors.toList());
						
						
						long c=bestEdge.stream()
						        .filter(te->{
						        	if(te.v()!=null){
						        		if(te.v().getOrder()==1){
						        			te.v().setOrder(2);
						        		}
						        		return true;
						        	}
						        	return false;
						        })
						        .count();
						if(c==1){
							e.setOrder(2);
						}else if(c==2){
							e.setOrder(1);
						}
						
						//ctab.removeEdge(e);
				});
			
			
			
			//look for pentavalent Carbons
			ctab.getNodes()
			    .stream()
			    .filter(n->n.getSymbol().equals("C"))
			    .filter(n->n.getValanceTotal()>=5)
			    .forEach(n->{
			    	//only deal with bad double bonds
			    	if(n.getEdgeCount()<=4){
			    		int toomuch = n.getValanceTotal()-4;
			    		if(toomuch==1 && n.getEdgeCount()==4){
			    			Optional<Edge> ope=n.getEdges()
			    			 .stream()
			    			 .filter(e->e.getOrder()==1)
			    			 .filter(e->e.getDashed())
			    			 .findFirst();
			    			if(ope.isPresent()){
			    				ope.ifPresent(e->{
			    					ctab.removeEdge(e);
			    				});
			    			}else{
			    				n.getEdges()
				    			 .stream()
				    			 .filter(e->e.getOrder()>1)
				    			 .forEach(e->e.setOrder(1));
			    			}
			    		}else if(toomuch==1 && n.getEdgeCount()==3){
			    			Optional<Edge> opEdge=n.getEdges()
								    			   .stream()
								    			   .filter(e->e.getOrder()>2)
								    			   .findFirst();
			    			if(opEdge.isPresent()){
			    				opEdge.ifPresent(e->{
			    					e.setOrder(2);
			    				});
			    			}else{
			    				n.getEdges()
				    			   .stream()
				    			   .filter(e->e.getOrder()==2)
				    			   .map(e->Tuple.of(e, GeomUtil.growLine(e.getLine(), ctab.getAverageBondLength()/3.0)))
				    			   .map(t->Tuple.of(t.k(),rejBondOrderLines.stream()
				    					                        .map(l->Tuple.of(l,l.centerPoint()))
				    					                        .filter(p->t.v().contains(p.v()))
				    					                        .mapToDouble(t1->t1.k().length())
				    					                        .findAny()))
				    			   
				    			   //.filter(t->t.v().isPresent())
				    			   .map(Tuple.vmap(v->v.orElse(0.0)))
				    			   .map(t->t.withVComparator())
				    			   .min(Comparator.naturalOrder())
				    			   .map(e->e.k())
				    			   .ifPresent(e->{
				    				   e.setOrder(1);
				    			   });
				    			   ;
			    				
			    				
			    				
			    			}
			    		}
			    	}
			    });
			
			//Here, we should remove some bad cage bonds
			
			
			//In this case, it's all very short bonds in 2 4-membered ring 
			//and all other edges in that ring are more than 30% longer than 
			//the edge in question
			ctab.getEdges()
			    .stream()
			    .filter(e->!e.isInventedBond())
			    .filter(e->e.getOrder()==1)
			    .filter(e->e.getEdgeLength()<ctab.getAverageBondLength()*0.5)
			    .collect(Collectors.toList())
			    .forEach(e->{
			    	List<Ring> rings = e.getAllRings().stream().filter(r->r.size() ==4).collect(Collectors.toList());
			    	
			    	if(rings.size()==2){
			    		boolean smallestEdge=rings.stream()
			    		     .allMatch(r->r.getEdges().stream().filter(e2->e2!=e).allMatch(e2->e2.getEdgeLength()>1.3*e.getEdgeLength()));
			    		if(smallestEdge){
			    			ctab.removeEdge(e);
			    		}
			    	}
			    });
			
			
			
			
			//Find floating methyls
			
			
			List<Shape> dashShapes = new ArrayList<Shape>();
			
			List<List<Node>> toMergeNodes = new ArrayList<>();
			
			
			
			maybeDashCollection
					.forEach(bshape->{
								//realRescueOCRCandidates.add(GeomUtil.growShapeNPoly(bshape,2,12));
							
							
								//Looks very promising
								dashShapes.add(GeomUtil.growShapeNPoly(bshape,2,12));
								Point2D[] pts=GeomUtil.getPairOfFarthestPoints(bshape);
								List<Node> forN1=ctab.getNodes()
								    .stream()
								    .filter(n->n.getPoint().distance(pts[0])< ctab.getAverageBondLength()*0.6)
								    .collect(Collectors.toList());
								List<Node> forN2=ctab.getNodes()
									    .stream()
									    .filter(n->n.getPoint().distance(pts[1])< ctab.getAverageBondLength()*0.6)
									    .filter(n->!forN1.contains(n))
									    .collect(Collectors.toList());
								
								
								
								if(forN1.size()+forN2.size()==0)return;
								LineWrapper splitLine=GeomUtil.findLongestSplittingLine(bshape);
								
								
								
								
								if(forN1.size()+forN2.size()>1){
									if(forN1.size()>1){
										//probably merge
										List<Node> nadd = forN1.stream().filter(nn->!forN2.contains(nn))
												.collect(Collectors.toList());
										toMergeNodes.add(nadd);
									}else if(forN2.size()>1){
										//probably merge
										List<Node> nadd = forN2.stream().filter(nn->!forN1.contains(nn))
												.collect(Collectors.toList());
										toMergeNodes.add(nadd);
									}else{
										if(forN1.size()==1 && forN2.size()==1){
											//probably reset the nodes now
											//realRescueOCRCandidates.add(GeomUtil.growShapeNPoly(bshape,2,12));
											
											realRescueOCRCandidates.add(splitLine.getLine());
											Node n1=forN1.get(0);	
											Node n2=forN2.get(0);
											Point2D np1=splitLine.projectPointOntoLine(n1.getPoint());
											Point2D np2=splitLine.projectPointOntoLine(n2.getPoint());
											n1.setPoint(np1);
											n2.setPoint(np2);
										}
										
									}
									return;
								}
								
								Node pnode=null;
								Point2D newPoint=pts[0];
								if(!forN1.isEmpty()){
									newPoint=pts[1];
									pnode=forN1.get(0);
								}else{
									pnode=forN2.get(0);
								}
								
								double ndist=newPoint.distance(pnode.getPoint());
								
								if(ndist<ctab.getAverageBondLength()*1.3 && ndist>ctab.getAverageBondLength()*0.6){
									//looks good
									Point2D cShape=GeomUtil.centerOfMass(bshape);
									Line2D nline = new Line2D.Double(pnode.getPoint(),cShape);
									double len = ctab.getAverageBondLength();
									Point2D op = GeomUtil.resizeLine(nline,len).getP2();
									
									Node otherNode=ctab.getNodes()
									    .stream()
									    .map(n->Tuple.of(n,n.getPoint().distance(op)).withVComparator())
									    .filter(t->t.v()<len*0.3)
									    .max(Comparator.reverseOrder())
									    .map(t->t.k())
									    .orElse(null);
									
									if(otherNode==null){
										Node realNode=ctab.addNode(op);
										BranchNode bn=bestGuessOCR.entrySet().stream()
																 .map(Tuple::of)
																 .filter(t->GeomUtil.growShapeNPoly(t.k(), len*0.3, 12).contains(op))
																 .map(Tuple.vmap(s1->BranchNode.interpretOCRStringAsAtom2(s1)))
																 .findFirst()
																 .map(t->t.v())
																 .orElse(null);
										if(bn!=null && bn.isRealNode()){
											realNode.setSymbol(bn.getSymbol());
										}
										
									}else{
										
										otherNode.setPoint(op);
									}
								}
						
					});
			
			
			
			if(!toMergeNodes.isEmpty()){
				toMergeNodes.forEach(nm->{
					//nervous about this
					ctab.mergeNodes(nm.stream().map(n->n.getIndex()).collect(Collectors.toList()), pl->pl.stream().collect(GeomUtil.averagePoint()));		
				});
				ctab.standardCleanEdges();
			}
			if(DEBUG)ctabRaw.add(ctab.cloneTab());
			        
			//Not sure about this, sometimes want to do a final merge
			ctab.getNodes()
			    .stream()
			    .filter(n->!n.isInvented())
			    .collect(GeomUtil.groupThings(t->{
			    	Node n1=t.k();
			    	Node n2=t.v();
			    	if(n1.distanceTo(n2)<ctab.getAverageBondLength()*0.2){
			    		if(n1.connectsTo(n2)){
			    			return true;
			    		}
			    	}
			    	return false;
			    }))
			    .stream()
			    .filter(nl->nl.size()==2)
			    .forEach(nlist->{
			    	Point2D p1 = nlist.stream().map(n->n.getPoint()).collect(GeomUtil.averagePoint());
			    	ctab.mergeNodes(nlist.stream().map(n->n.getIndex()).collect(Collectors.toList()), (pl)->p1);
			    	ctab.standardCleanEdges();
			    });
			

			@SuppressWarnings("unchecked")
			List<Tuple<Edge,WedgeInfo>> winfo=(List<Tuple<Edge, WedgeInfo>>) ctab.getEdges()
					.stream()
					.filter(e->!e.isInventedBond())
					.map(e->{
						LineWrapper useLine=GeomUtil.getLinesNotInside(e.getLine(), growLikelyOCR)
								.stream()
								.map(l->LineWrapper.of(l))
								.max(Comparator.naturalOrder())
								.orElse(null);
						if(useLine!=null){
							
							Point2D c=useLine.centerPoint();
							
							Optional<Shape> isDash=dashShapes.stream()
									  .filter(d->d.contains(c))
									  .findFirst();
							
							boolean isDotted=dottedLines.stream().map(dl->dl.growLine(ctab.getAverageBondLength()*.3))
							                   .filter(sl->sl.contains(c))
							                   .findAny()
							                   .isPresent();
							if(isDotted){
								e.setDashed(true);
							}else{
							
								if(isDash.isPresent()){
									
									e.setDashed(true);
									if(e.getOrder()!=1){
										if(e.getRealNode1().getSymbol().equals("C") && e.getRealNode2().getSymbol().equals("C")){
												//do nothing
										}else{
											e.setOrder(1);
										}							
									}
									
									Point2D cmass=GeomUtil.centerOfMass(isDash.get());
									if(e.getRealNode1().getPoint().distance(cmass)< e.getRealNode2().getPoint().distance(cmass)){
										e.switchNodes();
									}
								}else{
									//IDK
									if(e.getRealNode1().getSymbol().equals("C") && e.getRealNode2().getSymbol().equals("C")){
										//well, sometimes it's still okay to be a dash, if we _really_ think that
										//it's a "regular" dashed line "- - - - "
										//instead of the wedge-like dash line "| | | | |"
										//Not sure how to check for this
										if(e.getDashed() && e.getOrder()==1){
											//Keep it _only_ if there are 
											//at least 3 lines that are found in the area around the bond
											double grow=ctab.getAverageBondLength()*0.1;
											Shape lshape=useLine.growLine(grow);
											List<Shape> lineShapes=lines.stream()
											     .filter(l->!bestGuessOCR.keySet().stream().filter(ss->ss.contains(l.centerPoint())).findAny().isPresent())
											     .filter(l->l.absCosTheta(useLine)>0.8)
											     .filter(l->lshape.contains(l.centerPoint()))
											     .map(l->l.growLine(grow))
											     .collect(Collectors.toList());
											
											if(lineShapes.size()<3){
												e.setDashed(false);
											}else{
												double tarea=lineShapes.stream()
												                       .mapToDouble(s->GeomUtil.area(s))
												                       .sum();
												if(tarea>0.7*GeomUtil.area(lshape)){
													e.setDashed(false);
												}
											}
										}
									}
									return bitmap.getconfexHullAlongLine(useLine.getLine())
											.map(w->Tuple.of(e,w));
								}
							}
						}
						return Optional.empty();
					})
					.filter(t->t.isPresent())
					.map(o->o.get())
					.collect(Collectors.toList());
			
			Predicate<Node> couldBeStereoCenter = (n1)->n1.getEdgeCount()>=3 && n1.getSymbol().equals("C") && !n1.getEdges().stream().filter(e1->e1.getOrder()>1).findAny().isPresent();
			
			double averageThickness = winfo.stream()
										   .filter(t->Math.abs(t.v().getCorrel()) < WEDGE_LIKE_PEARSON_SCORE_CUTOFF)
										   .filter(t->t.v().pctOfHull()>0.5)
					                       .mapToDouble(t->t.v().getAverageThickness())
					                       .average()
					                       .orElse(2);
			
			Set<Edge> thickEdges = new HashSet<Edge>();
			Set<Edge> wedgeEdges = new HashSet<Edge>();
			
			
			winfo.forEach(t->{
				Edge e=t.k();
				
				WedgeInfo s = t.v();
				
				//realRescueOCRCandidates.add(s.getHull());
				
				
				double wl=s.getCorrel();
				//TODO: we need something for thick bonds which are _not_ strictly wedges, but are meant to be
				//wedges.
				if(s.getAverageThickness()>averageThickness*1.9){
					if(s.getOnPixels()>s.getArea()*0.5 && s.pctOfHull()>0.5){
						
						double cutoff=WEDGE_LIKE_PEARSON_SCORE_CUTOFF;
						if(e.getRealNode1().getEdges().stream().filter(ed->ed.getOrder()>1).findAny().isPresent()){
							cutoff=WEDGE_LIKE_PEARSON_SCORE_CUTOFF_DOUBLE;
						}
						if(e.getRealNode2().getEdges().stream().filter(ed->ed.getOrder()>1).findAny().isPresent()){
							cutoff=WEDGE_LIKE_PEARSON_SCORE_CUTOFF_DOUBLE;
						}
						
						
						if(wl>cutoff){
							e.setWedge(true);
							wedgeEdges.add(e);
						}else if(wl<-cutoff){
							e.setWedge(true);
							wedgeEdges.add(e);
							e.switchNodes();
						}else if(s.getAverageThickness()>averageThickness*2.3){
							//very thick line
							thickEdges.add(e);
							e.setWedge(true);
						}
						if(e.getWedge()){
							if(e.getRealNode1().getEdgeCount()<3 && e.getRealNode2().getEdgeCount()>=3){
								e.switchNodes();
							}
							//not the best yet
							if(couldBeStereoCenter.test(e.getRealNode2()) && ! couldBeStereoCenter.test(e.getRealNode1())){
								e.switchNodes();
							}
						}
					}
				}
			});
			
			
			if(thickEdges.size()>0 && wedgeEdges.size()>0){
				//There are 2 kinds of wedges here, disable all thick ones,
				//at least that share a node with a wedge/
				thickEdges.stream()
						  .filter(e->{
							  if(e.getRealNode1().getEdges().stream().filter(e2->e2!=e).filter(e2->e2.getWedge()||e2.getDashed()).findAny().isPresent() ||
								 e.getRealNode2().getEdges().stream().filter(e2->e2!=e).filter(e2->e2.getWedge()||e2.getDashed()).findAny().isPresent()){
								  return true;
							  }  
							  return false;
						  })
				          .forEach(w->w.setWedge(false));
			}
			
			

			GeomUtil.eachCombination(ctab.getNodes().stream().filter(n->!n.isInvented()).collect(Collectors.toList()))
					.filter(t->t.k().distanceTo(t.v())<1.5*ctab.getAverageBondLength())
					.filter(t->!t.k().getBondTo(t.v()).isPresent())
					.forEach(t1->{
						Line2D l2 = new Line2D.Double(t1.k().getPoint(),t1.v().getPoint());
						LineWrapper useLine=GeomUtil.getLinesNotInside(l2, growLikelyOCR)
								.stream()
								.map(l->LineWrapper.of(l))
								.max(Comparator.naturalOrder())
								.orElse(null);
						
						if(useLine==null)return;
						List<LineWrapper> lt=polygons.stream()
								.filter(s->GeomUtil.getIntersection(s, useLine.getLine()).isPresent())
								.map(s->GeomUtil.findLongestSplittingLine(s))
								.filter(l->l.length()<ctab.getAverageBondLength())
								.collect(Collectors.toList());
						boolean found=false;
						if(lt.size()>2){
							found=true;
							
						}else if(lt.size()>0){
							LineWrapper lwo=LineWrapper.of(l2);
							found=lt.stream()
							  .filter(lw->lw.absCosTheta(lwo)>0.9)
							  .filter(lw->lw.length()<useLine.length())
							  .filter(lw->lw.centerPoint().distance(useLine.centerPoint())<ctab.getAverageBondLength()/4.0)
							  .findAny()
							  .isPresent();
						}
						if(found){
							ctab.addEdge(t1.k().getIndex(), t1.v().getIndex(), 1);
							Edge e=ctab.getEdges().get(ctab.getEdges().size()-1);
							e.setDashed(true);
						}
					});
			
			ctab.getEdges()
			    .stream()
			    .filter(e->e.getDashed())
			    .map(e->Tuple.of(e,GeomUtil.findCenterOfShape(e.getLine())))
			    .filter(e->!dashShapes.stream().filter(d->d.contains(e.v())).findAny().isPresent())
			    .map(t->t.k())
			    .collect(Collectors.toList())
			    .forEach(t->{
			    	//might be bad edges
			    	Shape sl=GeomUtil.growLine(t.getLine(),ctab.getAverageBondLength()/3);
			    	boolean findLines=lines.stream()
			    						   .map(l->l.centerPoint())
			    	                       .filter(l->sl.contains(l))
			    	                       .findAny()
			    	                       .isPresent();
			    	//realRescueOCRCandidates.add(sl);
			    	if(!findLines){
			    		ctab.removeEdge(t);
			    	}
			    });

			ctab.getEdges()
				.stream()
				.filter(e->e.getDashed())
				.forEach(e->{
					if(couldBeStereoCenter.test(e.getRealNode2()) && ! couldBeStereoCenter.test(e.getRealNode1())){
						e.switchNodes();
					}
				});


			
	
//			
			

			
			ctab.getRings()
			    .stream()
			    .filter(r->r.size()==6)
			    .forEach(r->{
			    	List<Edge> doubles=r.getEdges()
			    	 .stream()
			    	 .filter(e->e.getOrder()==2)
			    	 .collect(Collectors.toList());
			    	
			    	if(doubles.size()==2){
			    		Set<Node> sp2s = doubles.stream().flatMap(e->Stream.of(e.getRealNode1(),e.getRealNode2())).collect(Collectors.toSet());
			    		
			    		Edge maybeDouble=r.getEdges()
			    		 .stream()
			    		 .filter(e->!sp2s.contains(e.getRealNode1()))
			    		 .filter(e->!sp2s.contains(e.getRealNode2()))
			    		 .findFirst()
			    		 .orElse(null);
			    		if(maybeDouble!=null){
			    			
			    			if(maybeDouble.getRealNode1().getValanceTotal()==4 || maybeDouble.getRealNode2().getValanceTotal()==4){
			    				return;
			    			}
			    			
			    			Shape ls = GeomUtil.growLine(maybeDouble.getLine(), ctab.getAverageBondLength()*0.3);
			    			
			    			long c=lines.stream()
			    			     .map(l->Tuple.of(l, l.centerPoint()))
			    			     .filter(t->ls.contains(t.v()))
			    			     .filter(t->t.k().absCosTheta(LineWrapper.of(maybeDouble.getLine()))>0.8)
			    			     .count();
			    			
			    			if(c>1){
			    				maybeDouble.setOrder(2);
			    			}
			    			     
			    		}
			    		
			    	}
			    });
			
			ctab.getRings()
			    .stream()
			    .filter(r->r.size()==6)
			    .filter(r->r.isConjugated())
			    .forEach(r->{
			    	Point2D center=GeomUtil.centerOfMass(r.getConvexHull());
			    	Shape p=GeomUtil.convexHull2(GeomUtil.makeNPolyCenteredAt(new Point2D.Double(0,0), 6, 100));
			    	//Shape p2=GeomUtil.convexHull2(GeomUtil.makeNPolyCenteredAt(center, 6, ctab.getAverageBondLength()));
			    	Point2D anchor = r.getNodes().stream()
			    			          .map(n->Tuple.of(n, n.getPoint().distance(center)).withVComparator())
			    			          .max(Comparator.naturalOrder())
			    			          .map(t->t.k().getPoint())
			    			          .orElse(null);
			    	Line2D nline = new Line2D.Double(center,anchor);
			    	AffineTransform at=GeomUtil.getTransformFromLineToLine(new Line2D.Double(new Point2D.Double(0,0),new Point2D.Double(100,0)),nline,false);
			    	Shape ns=at.createTransformedShape(p);
			    	
			    	Point2D[] verts2 = GeomUtil.vertices(ns);
			    	
			    	
			    	r.getNodes()
			    	 .forEach(n->{
			    		Point2D np=Arrays.stream(verts2)
			    		      .map(v->Tuple.of(v, v.distance(n.getPoint())).withVComparator())
			    		      .min(Comparator.naturalOrder())
			    		      .map(t->t.k())
			    		      .orElse(null);
			    		n.setPoint(np);
			    	 });
			    	
			    });
//			ctab.getRings()
//			    .stream()
//			    .filter(r->r.size()==4)
//			    .peek(r->System.out.println("Found ring:" + r.size()))
//			    .forEach(r->{
//			    	Point2D center=GeomUtil.centerOfMass(r.getConvexHull());
//			    	Shape p=GeomUtil.convexHull2(GeomUtil.makeNPolyCenteredAt(new Point2D.Double(0,0), 4, 100));
//			    	//Shape p2=GeomUtil.convexHull2(GeomUtil.makeNPolyCenteredAt(center, 6, ctab.getAverageBondLength()));
//			    	Point2D anchor = r.getNodes().stream()
//			    			          .map(n->Tuple.of(n, n.getPoint().distance(center)).withVComparator())
//			    			          .max(Comparator.naturalOrder())
//			    			          .map(t->t.k().getPoint())
//			    			          .orElse(null);
//			    	Line2D nline = new Line2D.Double(center,anchor);
//			    	AffineTransform at=GeomUtil.getTransformFromLineToLine(new Line2D.Double(new Point2D.Double(0,0),new Point2D.Double(100,0)),nline,false);
//			    	Shape ns=at.createTransformedShape(p);
//			    	
//			    	Point2D[] verts2 = GeomUtil.vertices(ns);
//			    	
//			    	
//			    	r.getNodes()
//			    	 .forEach(n->{
//			    		Point2D np=Arrays.stream(verts2)
//			    		      .map(v->Tuple.of(v, v.distance(n.getPoint())).withVComparator())
//			    		      .min(Comparator.naturalOrder())
//			    		      .map(t->t.k())
//			    		      .orElse(null);
//			    		n.setPoint(np);
//			    	 });
//			    	
//			    	realRescueOCRCandidates.add(ns);
//			    	
//			    	//realRescueOCRCandidates.add(p2);
//			    	
//			    });


			if(DEBUG)ctabRaw.add(ctab.cloneTab());
			

		}
		if(Thread.currentThread().isInterrupted()){
			throw new InterruptedException();
		}
		rescueOCRShapes=realRescueOCRCandidates;
	
		if(DEBUG)ctabRaw.add(ctab.cloneTab());
		
		
		
		
		//fix bad Sulfurs
		ctab.getNodes().stream()
		    .filter(n->n.getSymbol().equals("S"))
		    .filter(n->n.getValanceTotal()==7)
		    .forEach(n->{
		    		n.getNeighborNodes().stream()
		    						    .filter(t->t.v().getOrder()==2)
		    		                    .filter(nn->nn.k().getSymbol().equals("N"))
		    		                    .findFirst()
		    		                    .ifPresent(nn->{
		    		                    	nn.v().setOrder(1);
		    		                    });
		    });
		
		
		//fix bad double-bond to stupid OCR
				ctab.getNodes().stream()
				    .filter(n->n.getSymbol().equals("H"))
				    .filter(n->n.getValanceTotal()>2)
				    .forEach(n->{
				    		//probably a C
				    		n.setSymbol("C");
				    });
				

		
		
		//charge bad nitrogens		
		ctab.getNodes().stream()
		    .filter(n->n.getSymbol().equals("N"))
		    .filter(n->n.getCharge()==0)
		    .forEach(n->{
		    		int so=n.getEdges().stream().mapToInt(e->e.getOrder()).sum();
		    		if(so>3){
		    			n.setCharge(so-3);
		    			
		    			n.getNeighborNodes().stream()
		    								.filter(t->t.v().getOrder()==1)
		    								.map(t->t.k())
		    								.filter(nn->nn.getSymbol().equals("O"))
		    								.filter(nn->nn.getCharge()==0)
		    								.filter(nn->nn.getEdgeCount()==1)
		    								.findFirst()
		    								.ifPresent(nn->{
		    									nn.setCharge(-1);
		    								});
		    			
		    		}
		    });
		
		//charge bad sulfurs
		ctab.getNodes().stream()
		    .filter(n->n.getSymbol().equals("S"))
		    .filter(n->n.getCharge()==0)
		    .forEach(n->{
		    		int so=n.getEdges().stream().mapToInt(e->e.getOrder()).sum();
		    		if(so==3){
		    			n.getNeighborNodes().stream()
		    								.map(t->t.k())
		    								.filter(nn->nn.getSymbol().equals("O"))
		    								.filter(nn->nn.getCharge()==0)
		    								.filter(nn->nn.getEdgeCount()==1)
		    								.findFirst()
		    								.ifPresent(nn->{
		    									nn.setCharge(-1);
		    									n.setCharge(1);
		    								});
		    			
		    		}
		    });
		
		if(DEBUG)ctabRaw.add(ctab.cloneTab());
		
		List<Shape> mightBeNegative=polygons.stream()
										    .filter(s->s.getBounds2D().getHeight()<ctab.getAverageBondLength()/10)
										    .filter(s->s.getBounds2D().getWidth()>1)
										    .filter(s->s.getBounds2D().getWidth()<ctab.getAverageBondLength()/2)
										    .filter(s->!likelyOCRAll.contains(s) || Optional.ofNullable(ocrAttempt.get(s)).map(l->l.get(0)).filter(t->t.k().toString().equals("-")).isPresent())
										    .collect(Collectors.toList());
		
		if(!mightBeNegative.isEmpty()){
			Set<Shape> already = new HashSet<Shape>();
			BiConsumer<Node, Integer> maybeCharge=(n,o)->{
		  		int so=n.getEdges().stream().mapToInt(e->e.getOrder()).sum();
	    		if(so==o){
		    		Point2D p=n.getPoint();
		    		Shape neg=mightBeNegative.stream().filter(s->!already.contains(s)).filter(s->GeomUtil.distanceTo(s, p)<ctab.getAverageBondLength()*0.7).findAny().orElse(null);
		    		
		    		if(neg!=null){
		    			Point2D ncenter=GeomUtil.findCenterOfShape(neg);	
		    			boolean inLine=n.getEdges().stream()
						    			.map(e->e.getLine())
						    			.map(l->GeomUtil.growLine(l, averageWidthOCRFinal[0]))
						    			.filter(s->s.contains(ncenter))
						    			.findAny()
						    			.isPresent();
		    			if(!inLine){
		    				already.add(neg);
		    				n.setCharge(-1);
		    			}else{
		    				if(ctab.getSumCharge()>0){
		    					already.add(neg);
		    					n.setCharge(-1);
		    				}
		    			}
		    		}
	    		}
			};
			
			Map<String,Integer> neededValances = new HashMap<String,Integer>();
			neededValances.put("O",1);
			neededValances.put("N",2);
			neededValances.put("Cl",1);
			
			ctab.getNodes().stream()
				.filter(n->neededValances.get(n.getSymbol())!=null)
				.filter(n->n.getCharge()==0)
			    .filter(n->!n.isInvented())
			    .forEach(n->maybeCharge.accept(n, neededValances.get(n.getSymbol())));
			
			int sc=ctab.getSumCharge();
			
			if(sc>0){
				ctab.getNodes().stream()
							    .filter(n->n.getSymbol().equals("C"))
							    .filter(n->n.getCharge()==0)
							    .filter(n->n.getEdgeCount()==2)
							    .filter(n->!n.isInvented())
							    .filter(n->!n.getEdges().stream().filter(e->e.getDashed()).findAny().isPresent())
							    .forEach(n->maybeCharge.accept(n, 2));
			}
			if(ctab.getSumCharge()>0){
				mightBeNegative.stream()
							   .filter(s->!already.contains(s))
							   .map(s->Tuple.of(s,s))
							   .map(Tuple.vmap(s->bestGuessOCR.keySet()
									                          .stream()
									                          .map(os->Tuple.of(os,GeomUtil.distance(os, s)).withVComparator())
									                          .min(Comparator.naturalOrder())))
							   .filter(t->t.v().isPresent())
							   .map(Tuple.vmap(v1->v1.get()))
							   .filter(t->t.v().v()<ctab.getAverageBondLength()*0.3)
							   .forEach(t->{
								  Shape neg=t.k();
								  Shape ocr=t.v().k();
								  //for right now, only charge OH's
								  ctab.getAllNodesInsideShape(ocr, 3)
								      .stream()
								      .map(n->Tuple.of(Stream.of(n),n.getNeighborNodes().stream().map(t1->t1.k())))
								      .flatMap(nt->Stream.concat(nt.k(), nt.v()))
								      .distinct()
								      .filter(n->n.getSymbol().equals("O"))
								      .filter(n->n.getValanceTotal()==1)
								      .filter(n->n.getCharge()==0)
								      .findFirst()
								      .ifPresent(ncharge->{
								    	 ncharge.setCharge(-1); 
								      });
							   });
							   
			}
		}
		
		ctab.getNodes()
			.stream()
			.filter(ca->ca.getSymbol().equals("S"))
			.filter(ca->ca.getValanceTotal()==5)
			.forEach(ca->{
				ca.getNeighborNodes().stream()
				  .filter(cn->cn.k().getSymbol().equals("O"))
				  .filter(cn->cn.k().getCharge()==0)
				  .filter(cn->cn.v().getOrder()==1)
				  .limit(1)
				  .forEach(co->{
					  co.v().setOrder(2);
				  });
			});
		
		
		ctab.getNodes()
			.stream()
			.filter(ca->ca.getSymbol().equals("N"))
			.filter(ca->ca.getValanceTotal()==4)
			.filter(ca->ca.getEdgeCount()==2)
			.forEach(ca->{
				ca.getNeighborNodes().stream()
				  .filter(cn->!cn.k().getSymbol().equals("C"))
				  .filter(cn->cn.k().getCharge()==0)
				  .filter(cn->cn.v().getOrder()==2)
				  .limit(1)
				  .forEach(co->{
					  co.v().setOrder(1);
				  });
			});
		
		 
		if(DEBUG)ctabRaw.add(ctab.cloneTab());
		
		ctab.simpleClean();
		//Make aromatic bonds
		circles.stream()
		       .map(c->GeomUtil.findCenterOfShape(c))
		       .forEach(cp->{
		    	   ctab.getEdgesWithCenterWithin(cp,ctab.getAverageBondLength())
				       .stream()
				       .filter(e->e.isRingEdge())
				       .forEach(Edge::setToAromatic);
		       });
		
		if(DEBUG)ctabRaw.add(ctab.cloneTab());
		
		
		    
		

	}

	public String toMol(){
		return ctab.toMol();
	}

	public Bitmap getBitmap() {
		return bitmap;
	}

	public Bitmap getThin() {
		return thin;
	}

	public List<Shape> getPolygons() {
		return polygons;
	}
	public Map<Shape,String> getBestGuessOCR(){
		return bestGuessOCR;
	}

	public List<Line2D> getLines() {
		return lines.stream().map(lw->lw.getLine()).collect(Collectors.toList());
	}

	public List<Line2D> getLinesJoined() {
		return linesJoined.stream().map(lw->lw.getLine()).collect(Collectors.toList());
	}

	public List<Tuple<Line2D, Integer>> getLinesOrder() {
		return linesOrder;
	}

	public Map<Shape, List<Tuple<Character, Number>>> getOcrAttmept() {
		return ocrAttempt;
	}

	public ConnectionTable getCtab() {
		return ctab;
	}



	public List<ConnectionTable> getCtabRaw() {
		return this.ctabRaw;
	}





}
