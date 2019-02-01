package tripod.molvec.algo;

import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
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

import gov.nih.ncats.chemkit.api.Chemical;
import tripod.molvec.Bitmap;
import tripod.molvec.Bitmap.WedgeInfo;
import tripod.molvec.CachedSupplier;
import tripod.molvec.algo.Tuple.KEqualityTuple;
import tripod.molvec.ui.FontBasedRasterCosineSCOCR;
import tripod.molvec.ui.SCOCR;
import tripod.molvec.ui.StupidestPossibleSCOCRSansSerif;
import tripod.molvec.ui.StupidestPossibleSCOCRSerif;
import tripod.molvec.util.CompareUtil;
import tripod.molvec.util.ConnectionTable;
import tripod.molvec.util.ConnectionTable.Edge;
import tripod.molvec.util.ConnectionTable.Node;
import tripod.molvec.util.GeomUtil;

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

	private List<Line2D> lines;
	private List<Line2D> linesJoined;
	private List<Tuple<Line2D,Integer>> linesOrder;    
	private Map<Shape,List<Tuple<Character,Number>>> ocrAttempt = new ConcurrentHashMap<>();
	private Map<Shape,String> bestGuessOCR = new HashMap<>();
	

	private ConnectionTable ctab;
	private List<ConnectionTable> ctabRaw = new ArrayList<ConnectionTable>();


	private final double MAX_REPS = 10;
	private final double INITIAL_MAX_BOND_LENGTH=Double.MAX_VALUE;
	private final double MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE = 1/3.5;
	private final double MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE_INITIAL = 1/2.9;

	private final double MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL = 1.3;
	private final double MIN_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL = 0.5;

	
	private final boolean REMOVE_NONSENSE_OCR_LINES = false;

	private final double MAX_TOLERANCE_FOR_DASH_BONDS = 3.0;
	private final double MAX_TOLERANCE_FOR_SINGLE_BONDS = 0.4;

	private final double OCRcutoffCosine=0.65;
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
	private final double MAX_THETA_FOR_OCR_SEPERATION=45 * Math.PI/180.0;


	//This number is likely one of the most important to adjust.
	//It may have to have some changes done to the algorithm using it too
	private final double MAX_BOND_RATIO_FOR_MERGING_TO_OCR=0.31;


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






	public StructureImageExtractor(byte[] file) throws IOException{
		load(file);
	}
	public StructureImageExtractor(File file) throws IOException{
		load(file);
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
		BondLikely
	}
	
	private static CharType OCRIsLikely(Tuple<Character,Number> tup){
		String t=tup.k()+"";
		if(     "I".equalsIgnoreCase(t) || 
				"L".equalsIgnoreCase(t) ||
				"1".equalsIgnoreCase(t) ||
				"-".equalsIgnoreCase(t) ||
				"/".equalsIgnoreCase(t) ||
				"K".equalsIgnoreCase(t) ||
				"t".equalsIgnoreCase(t) ||
				"(".equalsIgnoreCase(t) ||
				")".equalsIgnoreCase(t) ||
				"Y".equalsIgnoreCase(t) ||
				"W".equalsIgnoreCase(t) ||
				"f".equals(t) ||
				"\\".equalsIgnoreCase(t)){
			return CharType.BondLikely;
		}else if(
				"2".equalsIgnoreCase(t) ||
				"3".equalsIgnoreCase(t) ||
				"4".equalsIgnoreCase(t) ||
//				!"5".equalsIgnoreCase(t) &&
				"6".equalsIgnoreCase(t) ||
				"7".equalsIgnoreCase(t) ||
				"9".equalsIgnoreCase(t)){
			return CharType.NumericLikely;
		}
		return CharType.ChemLikely;
	}



	private void processOCR(SCOCR socr, List<Shape> polygons,Bitmap bitmap, Bitmap thin, BiConsumer<Shape,List<Tuple<Character,Number>>> onFind){
		/*
		 * Looks at each polygon, and gets the likely OCR chars.
		 */   
		List<Shape> toAddShapes = Collections.synchronizedList(new ArrayList<Shape>());
		List<Shape> toRemoveShapes = Collections.synchronizedList(new ArrayList<Shape>());

		polygons.stream()
//		.collect(Collectors.toList())
//		.stream()
		
		.parallel()
		.forEach(s->{
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
						
						
						if(bounds2d.getWidth()>2.3*bounds2d.getHeight()){
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
				}else{
					onFind.accept(s, ll);
				}
				
				
			}
		});
		
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

			List<Tuple<Character,Number>> potential = socr.getNBestMatches(4,
					bitmap.crop(s),
					thin.crop(s))
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

	private void load(byte[] file) throws IOException{
		load(bitmap = Bitmap.read(file).clean());

	}
	private void load(File file) throws IOException{
		load(bitmap = Bitmap.read(file).clean());

	}

	private void load(Bitmap aBitMap) throws IOException{


		ctabRaw.clear();
		ocrAttempt.clear();

		SCOCR[] socr=new SCOCR[]{OCR_DEFAULT.orElse(OCR_BACKUP, OCRcutoffCosine)};

		double[] maxBondLength=new double[]{INITIAL_MAX_BOND_LENGTH};    




		polygons = bitmap.connectedComponents(Bitmap.Bbox.DoublePolygon);

		boolean isLarge = false;
		if (!polygons.isEmpty()) {
			isLarge = polygons.size() > 4000;
		}
		

		

		thin = bitmap.thin();

		// segments are generated for thinned bitmap only, since
		//  it can quite noisy on normal bitmap!
		if (isLarge) {
			throw new IllegalStateException("Cannot support images with over 4000 line segments at this time");
		}

		List<Shape> likelyOCR= Collections.synchronizedList(new ArrayList<Shape>());
		List<Shape> likelyOCRNumbers= Collections.synchronizedList(new ArrayList<Shape>());
		List<Shape> likelyOCRNonBond= Collections.synchronizedList(new ArrayList<Shape>());
		List<Shape> likelyOCRAll=Collections.synchronizedList(new ArrayList<Shape>());
		List<Shape> likelyOCRIgnore = Collections.synchronizedList(new ArrayList<Shape>());
		
		List<Shape> ocrRescues = new ArrayList<Shape>();
		
		/*
		 * Looks at each polygon, and gets the likely OCR chars.
		 */   

		List<Shape> initialDebug = Collections.synchronizedList(new ArrayList<>());
		
		
		
		

		processOCR(socr[0],polygons,bitmap,thin,(s,potential)->{
			ocrAttempt.put(s, potential);
			initialDebug.add(s);
			
			if(potential.stream().filter(e->e.v().doubleValue()>OCRcutoffCosine).findAny().isPresent()){
				CharType ct=OCRIsLikely(potential.get(0));
				if(ct.equals(CharType.ChemLikely)){
					likelyOCR.add(s);
					likelyOCRNonBond.add(s);
					
				}else if(ct.equals(CharType.NumericLikely)){
					likelyOCRNonBond.add(s);
					likelyOCRNumbers.add(s);
				}
				likelyOCRAll.add(s);

			}
		});

		//if(true)return;
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
				.collect(Collectors.toList());
		
		
		
		
		ctabRaw.clear();

		AtomicBoolean foundNewOCR=new AtomicBoolean(true);
		int maxFullRepeats=5;
		int repeats=0;
		List<Shape> realRescueOCRCandidates = new ArrayList<>();
		
		double[] averageHeightOCRFinal = new double[]{0};
		double[] averageWidthOCRFinal = new double[]{0};
		
		while(foundNewOCR.get() && repeats<maxFullRepeats){
			repeats++;
			foundNewOCR.set(false);

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
				Tuple<Shape,Double> shape1=GeomUtil.findClosestShapeTo(likelyOCRNonBond, l.getP1());
				if(shape1.v()>OCR_TO_BOND_MAX_DISTANCE){
					return false;
				}
				Tuple<Shape,Double> shape2=GeomUtil.findClosestShapeTo(likelyOCRNonBond, l.getP2());
				if(shape2.v()>OCR_TO_BOND_MAX_DISTANCE){
					return false;
				}
				if(shape1.k()==shape2.k()){
					return true;
				}
				
				boolean anyOutside=GeomUtil.getLinesNotInside(l, Arrays.asList(shape1.k(),shape2.k()))
				        .stream()
				        .filter(l1->l1!=null)
				        .filter(l1->GeomUtil.length(l1)>1)
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
										 .map(l->Tuple.of(l,GeomUtil.findCenterOfShape(l)))
									     .filter(t->!likelyOCRIgnore.stream()
									    		    .filter(s->s.contains(t.v()))
									    		    .findAny().isPresent())
									     .map(t->t.k())
									     .collect(Collectors.toList());
			


			List<Line2D> smallLines=useLines.stream()
					.filter(tryToMerge)
					.collect(Collectors.toList());

			List<Line2D> bigLines=useLines.stream()
					.filter(tryToMerge.negate())
					.collect(Collectors.toList());

			smallLines=smallLines.stream()
					.flatMap(l->{
						//			                	  double len=GeomUtil.length(l);
						//			                	  if(len<MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS && len>MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS){
						//			                		  return GeomUtil.splitLineIn2(l).stream();  
						//			                	  }
						return Stream.of(l);
					})
					.collect(Collectors.toList());

			// smallLines= thin.combineLines(smallLines, MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS, MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS_THIN, MAX_POINT_DISTANCE_TO_BE_PART_OF_MULTI_NODE,MAX_ANGLE_FOR_JOINING_SEGMENTS,MIN_SIZE_FOR_ANGLE_COMPARE_JOINING_SEGMENTS);

			smallLines= bitmap.combineLines(smallLines, MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS, MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS_FULL, MAX_POINT_DISTANCE_TO_BE_PART_OF_MULTI_NODE,MAX_ANGLE_FOR_JOINING_SEGMENTS,MIN_SIZE_FOR_ANGLE_COMPARE_JOINING_SEGMENTS);
			smallLines= thin.combineLines(smallLines, MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS, MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS_THIN, MAX_POINT_DISTANCE_TO_BE_PART_OF_MULTI_NODE,MAX_ANGLE_FOR_JOINING_SEGMENTS,MIN_SIZE_FOR_ANGLE_COMPARE_JOINING_SEGMENTS);

			List<Line2D> removedTinyLines =smallLines.stream()
					.filter(l->GeomUtil.length(l)<=MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS)
					.collect(Collectors.toList());

			List<Point2D> removedTinyVertices = removedTinyLines.stream()
					.flatMap(l->Stream.of(l.getP1(),l.getP2()))
					.collect(Collectors.toList());



			smallLines=smallLines.stream()
					.filter(l->GeomUtil.length(l)>MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS)
					.collect(Collectors.toList());



			List<Point2D> verts = smallLines.stream()
					.flatMap(l->Stream.of(l.getP1(),l.getP2()))
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



			linesJoined=Stream.concat(bigLines.stream(),
					smallLines.stream())
					.collect(Collectors.toList());


			double largestBond=smallLines.stream()
					.mapToDouble(l->GeomUtil.length(l))
					.max()
					.orElse(0);

			double averageLine=smallLines.stream()
					.mapToDouble(l->GeomUtil.length(l))
					.average()
					.orElse(0);

			if(largestBond>2*averageLine){
				largestBond=1.4*averageLine;
			}

			List<Line2D> preprocess= GeomUtil.reduceMultiBonds(linesJoined, MAX_ANGLE_FOR_PARALLEL, MAX_DISTANCE_TO_MERGE_PARALLEL_LINES, MIN_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS,0,MAX_DELTA_LENGTH_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC, (l)->{})
					.stream()
					.map(t->t.k())
					.collect(Collectors.toList());


			List<Line2D> rejBondOrderLines=new ArrayList<Line2D>();
			
			linesOrder=GeomUtil.reduceMultiBonds(preprocess, MAX_ANGLE_FOR_PARALLEL, largestBond/3, MIN_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS, MIN_BIGGER_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS,MAX_DELTA_LENGTH_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC,
					(rejLine->rejBondOrderLines.add(rejLine))
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
			
			while(tooLongBond){
				
				
				
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
						l-> (GeomUtil.length(l) < maxBondLength[0]))
						.mergeNodesCloserThan(MAX_DISTANCE_BEFORE_MERGING_NODES);
				ctabRaw.add(ctab.cloneTab());
				



				for(Shape s: likelyOCR){
					ctab.mergeAllNodesInsideCenter(s, OCR_TO_BOND_MAX_DISTANCE);
				}

				ctabRaw.add(ctab.cloneTab());


				List<List<Node>> newNodesForMerge = new ArrayList<>();

				Function<List<Node>,Point2D> bestIntersectionPoint = (nl)->{
					Point2D center= GeomUtil.findCenterOfVertices(nl.stream().map(n->n.getPoint()).collect(Collectors.toList()));
					double rad=ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE;

					List<Edge> el=nl.stream()
							.flatMap(n->n.getEdges().stream())
							.collect(Collectors.toList());


					List<Point2D> intersections=GeomUtil.eachCombination(el)
							.flatMap(t->{
								if(t.k()==t.v())return Stream.of(t.k().getPoint1(),t.k().getPoint2());
								return Stream.of(GeomUtil.intersection(t.k().getLine(),t.v().getLine()));
							})
							.filter(p->p!=null)
							.filter(p->p.distance(center)<rad)
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
							//polygons.add(candidate);
							
							return center;
						}
						//polygons.add(candidate);
						//return center;
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


				ctabRaw.add(ctab.cloneTab());


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


				ctabRaw.add(ctab.cloneTab());
				
				List<Edge> splitEdges = new ArrayList<Edge>();

				ctab.createNodesOnIntersectingLines(2, elist->{
					splitEdges.addAll(elist);
					return true;
				});

				ctab.getEdges()
			    .stream()
			    //.filter(e->e.getOrder()==1)
			    .filter(e->splitEdges.contains(e))
			    .collect(Collectors.toList())
			    .forEach(e->{
			    	Line2D ll=e.getLine();
			    	double totLen=GeomUtil.getLinesNotInside(ll,growLines).stream().mapToDouble(l->GeomUtil.length(l)).sum();
			    	if(totLen>GeomUtil.length(ll)*0.8){
			    		ctab.removeEdge(e);
			    	}
			    });

				ctabRaw.add(ctab.cloneTab());
				
				
				//This part needs work
//				ctab.getEdges()
//			    .stream()
//			    .filter(e->e.getOrder()==1)
//			    .collect(Collectors.toList())
//			    .forEach(e->{
//			    	Line2D ll=GeomUtil.getLinesNotInside(e.getLine(),likelyOCR).stream()
//			    						.map(l1->Tuple.of(l1,GeomUtil.length(l1)).withVComparator())
//			    						.max(Comparator.naturalOrder())
//			    						.map(t->t.k())
//			    						.orElse(null);
//			    	if(ll==null)return;
//			    	double totLen=GeomUtil.getLinesNotInside(ll,growLines).stream().mapToDouble(l->GeomUtil.length(l)).sum();
//			    	if(totLen>GeomUtil.length(ll)*0.9){
//			    		ctab.removeEdge(e);
//			    	}
//			    });


				ctab.mergeNodesCloserThan(ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE);
				ctab.standardCleanEdges();

				ctabRaw.add(ctab.cloneTab());

				

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
				
				ctabRaw.add(ctab.cloneTab());

				ctab.removeOrphanNodes();
				ctab.mergeNodesCloserThan(ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE);
				ctabRaw.add(ctab.cloneTab());
				
				for(Shape s: likelyOCR){
					ctab.mergeAllNodesInsideCenter(s, OCR_TO_BOND_MAX_DISTANCE);
				}

				
				ctabRaw.add(ctab.cloneTab());

				ctab.makeMissingNodesForShapes(likelyOCR,MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL,MIN_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL);

				Set<Node> toRemove = new HashSet<Node>();
				


				ctab.makeMissingBondsToNeighbors(bitmap,MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL,MAX_TOLERANCE_FOR_DASH_BONDS,likelyOCR,OCR_TO_BOND_MAX_DISTANCE, (t)->{
					
					//It could be that there is already a bond between the two nodes through a bad intermediate
					Edge e=t.v();
					Node n1=e.getRealNode1();
					Node n2=e.getRealNode2();
					List<Edge> existingEdges1=n1.getEdges();
					List<Edge> existingEdges2=n2.getEdges();
					List<Node> n1Neigh=existingEdges1.stream()
							.flatMap(ne->Stream.of(ne.getRealNode1(),ne.getRealNode2()))
							.filter(n->!n.equals(n1))
							.collect(Collectors.toList());

					List<Node> n2Neigh=existingEdges2.stream()
							.flatMap(ne->Stream.of(ne.getRealNode1(),ne.getRealNode2()))
							.filter(n->!n.equals(n2))
							.collect(Collectors.toList());
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


							}

						}
					}


					
					if(t.k()>MAX_TOLERANCE_FOR_SINGLE_BONDS){
						if(!alreadyExists){
							t.v().setDashed(true);
						}
					}
				});



				
				//fuzzy adding missing stuff
				ctabRaw.add(ctab.cloneTab());
				toRemove.forEach(n->ctab.removeNodeAndEdges(n));
				//ctab.removeOrphanNodes();

			



//				System.out.println("Removed bad edges:" + ctabRaw.size());
				ctabRaw.add(ctab.cloneTab());

				double avgBondLength=ctab.getAverageBondLength();
				maxBondLength[0]=avgBondLength*MAX_BOND_TO_AVG_BOND_RATIO_TO_KEEP;



				

				tooLongBond = ctab.getEdges()
						.stream()
						
						.filter(e->e.getEdgeLength()>maxBondLength[0])
						.findAny()
						.isPresent();
				if(tooLongBond){
					
					reps++;
				}
				if(reps>MAX_REPS)break;
			}
			realRescueOCRCandidates.addAll(rescueOCRCandidates);


			AtomicBoolean anyOtherIntersections = new AtomicBoolean(false);

			ctab.createNodesOnIntersectingLines(3, elist->{

				long longEnoughBonds=elist.stream()
						
						.filter(e->e.getEdgeLength()>MAX_BOND_TO_AVG_BOND_RATIO_FOR_INTERSECTION*ctab.getAverageBondLength())
						.count();
				if(longEnoughBonds<3)return false;
				anyOtherIntersections.set(true);
				return true;
			});

			if(anyOtherIntersections.get()){
				
				ctabRaw.add(ctab.cloneTab());
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

			List<Line2D> verticesJl=linesJoined.stream()
					.collect(Collectors.toList());
			List<Point2D> verticesJ=
					//Stream.concat(
					//			  removedTinyVertices.stream(), 
					verticesJl.stream()
					.flatMap(l->Stream.of(l.getP1(),l.getP2()))
					//.filter(pt->!likelyOCRAll.stream().filter(s->s.contains(pt)).findAny().isPresent())
					//		                 									     )
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

				boolean isresc=false;


				if(centerRescue!=null){
					cpt=centerRescue;
					isresc=true;
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



					Bitmap nmap=bitmap.crop(nshape);
					Bitmap nthinmap=thin.crop(nshape);
					if(nmap!=null && nthinmap!=null){
						
						nshape=GeomUtil.growShape(nshape, 2);
						nmap=bitmap.crop(nshape);
						nthinmap=thin.crop(nshape);

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

						nmap=bitmap.crop(nshape);
						nthinmap=thin.crop(nshape);
						
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

			ctabRaw.add(ctab.cloneTab());


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
					Bitmap nmap=bitmap.crop(nshape);
					Bitmap nthinmap=thin.crop(nshape);
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
					
					foundNewOCR.set(true);
				}
			});


			
			ctab.mergeNodesExtendingTo(likelyOCR,maxRatio,maxTotalRatio);

//			System.out.println("Extended to OCR:" + ctabRaw.size());
			ctabRaw.add(ctab.cloneTab());

			double cosThetaOCRShape =Math.cos(MAX_THETA_FOR_OCR_SEPERATION);
			
		
			double wid=averageWidthOCR;

			List<List<Shape>> ocrGroupList=GeomUtil.groupShapesIfClosestPointsMatchCriteria(likelyOCRAll, t->{
				Point2D[] pts=t.v();
				Shape[] shapes =t.k();

				List<Tuple<Character, Number>> attempt0 = ocrAttempt.get(shapes[0]);
				List<Tuple<Character, Number>> attempt1 = ocrAttempt.get(shapes[1]);
				String v1= (attempt0 ==null || attempt0.isEmpty())? "" : attempt0.get(0).k().toString();
				String v2= (attempt1 ==null || attempt1.isEmpty())? "" : attempt1.get(0).k().toString();
				if(v1.equals("\\") || v1.equals("/") || 
						v2.equals("\\") || v2.equals("/")){
					return false;
				}
				
				

				Line2D l2 = new Line2D.Double(pts[0],pts[1]);
				double dist=GeomUtil.length(l2);
				double cutoff=Math.max(ctab.getAverageBondLength()*MAX_BOND_RATIO_FOR_OCR_CHAR_SPACING,wid);
				
				if(dist>cutoff){
					return false;
				}

				Point2D cs1=GeomUtil.findCenterOfShape(shapes[0]);
				Point2D cs2=GeomUtil.findCenterOfShape(shapes[1]);

				Line2D cenLine = new Line2D.Double(cs1, cs2);
				double[] vec=GeomUtil.asVector(cenLine);
				double cosTheta=Math.abs(vec[0]/GeomUtil.length(cenLine));

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
				dontMerge.put("FF", Arrays.asList("F","F"));
				dontMerge.put("HOOH", Arrays.asList("HO","OH"));
				dontMerge.put("OHOH", Arrays.asList("OH","OH"));
				dontMerge.put("OHHO", Arrays.asList("OH","HO"));
				dontMerge.put("BrBr", Arrays.asList("Br","Br"));

				for(Tuple<Shape,Tuple<List<Shape>,String>> tt: toAdd){
					boolean removeBad=false;
					String val = tt.v().v();
					List<Shape> contains = tt.v().k();
					Shape parent=tt.k();
					
					
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
							String g1=val.substring(findex,keep.length());
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
					
					BranchNode bn = BranchNode.interpretOCRStringAsAtom2(val);
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
						foundNewOCR.set(true);
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
			
			bestGuessOCR.entrySet()
						.stream()
						.map(Tuple::of)
						.filter(t->t.v().equals("H"))
						.collect(Collectors.toList())
						.stream()
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
					.peek(t->System.out.println(bestGuessOCR.get(t)))
					.filter(s->BranchNode.interpretOCRStringAsAtom2(bestGuessOCR.get(s))!=null)
					.collect(Collectors.toList());

			
			ctabRaw.add(ctab.cloneTab());
			//ctab.removeOrphanNodes();

			List<Node> alreadyFixedNodes = new ArrayList<Node>();


			bestGuessOCR.entrySet()
			.stream()
			.map(Tuple::of)
			.map(Tuple.vmap(s->Tuple.of(s,(s.equals("H"))?1:0).withVComparator()))
			.map(t->t.withVComparator())
			.sorted()
			.filter(t->!t.v().k().startsWith("#")) //???
			.map(Tuple.vmap(t->t.k()))
			.forEach(shapeString->{
				Shape s= shapeString.k();
				String sym = shapeString.v();
				BranchNode actual = BranchNode.interpretOCRStringAsAtom2(sym);
				Point2D centert = GeomUtil.findCenterOfShape(s);

				if(actual!=null && actual.isRealNode()){
					if(sym.length()>1){
						List<Line2D> externalLines=ctab.getAllEdgesEntering(s, MAX_BOND_RATIO_FOR_MERGING_TO_OCR*ctab.getAverageBondLength())
								.stream()
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
							if(GeomUtil.findClosestShapeTo(ocrMeaningful, n.getPoint()).k() !=s){
								return false;
							}
						}
						if(!s.contains(n.getPoint()) && actual.isTerminal() && n.getEdgeCount()>1){
							//System.out.println("Term?");
							long cc=n.getNeighborNodes()
							 .stream()
							 .map(t->t.k())
							 .filter(nn->GeomUtil.distanceTo(s, nn.getPoint())<2)
							 .count();
							
							if(cc==0)return false;
						}
						if(alreadyFixedNodes.contains(n))return false;
						if(n.getEdgeCount()==0)return false;
						return true;
					},(l)->{

						boolean matchesOthers=l.stream()
								.map(pt->GeomUtil.findClosestShapeTo(ocrMeaningful, pt).k())
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


			
			ctabRaw.add(ctab.cloneTab());

			List<Shape> growLikelyOCRNonBond=likelyOCRNonBond.stream().map(s->GeomUtil.growShape(s, 2)).collect(Collectors.toList());
			List<Shape> growLikelyOCR=likelyOCR.stream().map(s->GeomUtil.growShape(s, 2)).collect(Collectors.toList());



			List<Line2D> lj =Stream.concat(removedTinyLines.stream(), linesJoined.stream())
					.flatMap(l->GeomUtil.getLinesNotInside(l, growLikelyOCRNonBond).stream())
					.filter(l->GeomUtil.length(l)>2) // needed?
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
						.map(t->{
							return Tuple.of(t,GeomUtil.add(t.k(), t.v()));
						})
						.map(Tuple.vmap(s->Tuple.of(s,GeomUtil.area(s)).withVComparator()))
						.map(t->t.withVComparator())
						.sorted()
						.map(Tuple.vmap(st->st.k()))
						.map(t->{
							//It should also have at least 1 line segment between the two
							Shape cshape = t.v();

							List<Line2D> opl=lj.stream()
									.filter(l1->cshape.contains(GeomUtil.findCenterOfShape(l1)))
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

				//don't add cross bonds


				int order=GeomUtil.groupThings(t.k(), tlines->{
					Line2D l1=tlines.k(); 
					Line2D l2=tlines.v();
					if(GeomUtil.cosTheta(l1, l2) > Math.cos(10*Math.PI/180)){
						return true;
					}
					return false;
				})
						.stream()
						.mapToInt(ll->ll.size())
						.max().getAsInt();
				ctab.addEdge(t.v().k().getIndex(), t.v().v().getIndex(), order);
			});

			ctab.standardCleanEdges();

			

			ctabRaw.add(ctab.cloneTab());

			ctab.makeMissingBondsToNeighbors(bitmap,MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL,MAX_TOLERANCE_FOR_SINGLE_BONDS,likelyOCR,OCR_TO_BOND_MAX_DISTANCE, (t)->{

				if(t.k()>MAX_TOLERANCE_FOR_SINGLE_BONDS){
					t.v().setDashed(true);
				}
			});

			
			ctabRaw.add(ctab.cloneTab());
			//remove triangles that are obviously wrong
			{
				double avgL = ctab.getAverageBondLength();
				Set<Edge> skip= new HashSet<Edge>();
				ctab.getEdges()
				.stream()
				.filter(e->e.getEdgeLength()>avgL)
				.map(l->Tuple.of(l,GeomUtil.length(l.getLine())).withVComparator())
				.sorted(Comparator.reverseOrder())
				.map(t->t.k())
				.filter(t->!skip.contains(t))
				.forEach(e->{
					Node n1= e.getRealNode1();
					Node n2= e.getRealNode2();
					List<KEqualityTuple<Node,Edge>> neigh1=n1.getNeighborNodes();
					List<KEqualityTuple<Node,Edge>> neigh2=n2.getNeighborNodes();
					List<KEqualityTuple<Node,Edge>> things1=neigh1.stream()
							.filter(ne->neigh2.contains(ne))
							.collect(Collectors.toList());
					List<KEqualityTuple<Node,Edge>> things2=neigh2.stream()
							.filter(ne->neigh1.contains(ne))
							.collect(Collectors.toList());
					List<KEqualityTuple<Node,Edge>> things = Stream.concat(things1.stream(), things2.stream())
							.collect(Collectors.toList());

					if(things.size()>0){
						
						Point2D p1=n1.getPoint();
						Point2D p2=n2.getPoint();
						Point2D p3=things.get(0).k().getPoint();

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

			ctabRaw.add(ctab.cloneTab());

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
			ctabRaw.add(ctab.cloneTab());
			double fbondlength=ctab.getAverageBondLength();

			

			
			//clean bad triple bonds
			ctab.getEdges().stream()
			.filter(e->e.getOrder()==3)
			.collect(Collectors.toList())
			.forEach(e->{
				Line2D lb = e.getLine();
				
				double len = e.getEdgeLength();
				
				double otherBondAverage = ctab.getEdges().stream().filter(e1->e1!=e).mapToDouble(e1->e1.getEdgeLength()).average().orElse(1);
				
				int n = (int)Math.round(len/otherBondAverage);
				
				if(n>1){
					Node n1=e.getRealNode1();
					Node n2=e.getRealNode2();
					
					Shape bigLineShape = GeomUtil.growLine(lb,len/3);
					
					Point2D apt=rejBondOrderLines.stream()
									 .filter(l->GeomUtil.length(l)>ctab.getAverageBondLength()*0.5)
					                 .map(l->GeomUtil.findCenterOfShape(l))
					                 .filter(p->bigLineShape.contains(p))
					                 .map(p->GeomUtil.projectPointOntoLine(lb,p))
					                 .collect(GeomUtil.averagePoint());
					
					
					List<Point2D> pts=GeomUtil.splitIntoNPieces(lb,n);
					Node pnode=n1;
					Edge closestEdge = null;
					double closestD = 999999;
					
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
					closestEdge.setOrder(3);
					//System.out.println("Setting order to 3");
					ctab.removeEdge(e);
					
				}
				
				
				//rejBondOrderLines
				
				
			});
			
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
							lneigh=pnode.getNeighborNodes().stream().map(n->n.k()).filter(n->n.getPoint().getX()<pnode.getPoint().getX())
									.map(n->Tuple.of(n,n.getPoint().getX()))
									.map(t->t.withVComparator())
									.max(Comparator.naturalOrder())
									.map(t->t.k())
									.orElse(null);
							rneigh=pnode.getNeighborNodes().stream().map(n->n.k()).filter(n->n.getPoint().getX()>pnode.getPoint().getX())
									.map(n->Tuple.of(n,n.getPoint().getX()))
									.map(t->t.withVComparator())
									.min(Comparator.naturalOrder())
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

								Node n= ctab.addNode(curN.suggestedPoint)
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
			
			ctabRaw.add(ctab.cloneTab());



			//Now get all the nodes with 2 edges which have shorter than average bond length,
			//and which are not in OCR, and where the sum of the distances of the bonds is within 97%
			//of the distance of the two neighbors. Those are nodes to remove
			List<Node> toRemove = new ArrayList<Node>();
			do{
				toRemove.clear();



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


			ctabRaw.add(ctab.cloneTab());


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

			});;

			ctab.standardCleanEdges();


			ctabRaw.add(ctab.cloneTab());


			ctab.getEdges()
			.stream()
			.filter(e->e.getEdgeLength()<ctab.getAverageBondLength()*0.55)
			.collect(Collectors.toList())
			.forEach(e->{
				//look at small bonds
				//Maybe merge the atoms?
				Node n1=e.getRealNode1();
				Node n2=e.getRealNode2();
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

			ctabRaw.add(ctab.cloneTab());
			
			List<Node> toRemoveNodesCage = new ArrayList<>();
			
			ctab.getNodes()
		    .stream()
		    .filter(n->n.getEdgeCount()==2)
		    .filter(n->n.getSymbol().equals("C"))
		    .filter(n->!n.isInvented())
		    .filter(n->n.getEdges().stream().filter(e->e.getOrder()==1).count()==2)
		    .filter(n->Optional.ofNullable(GeomUtil.findClosestShapeTo(likelyOCR, n.getPoint())).map(t->t.v()).orElse(100.0)>2)
		    .forEach(n->{
		    	List<Tuple<Node,Node>> tn=GeomUtil.eachCombination(n.getNeighborNodes())
		    	        .filter(t->{
		    	        	Node n1=t.k().k();
		    	        	Node n2=t.v().k();
		    	        	
		    	        	Line2D l = new Line2D.Double(n1.getPoint(), n2.getPoint());
		    	        	Point2D pp=GeomUtil.projectPointOntoLine(l, n.getPoint());
		    	        	if(pp.distance(n.getPoint())<ctab.getAverageBondLength()*0.01){
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
			ctabRaw.add(ctab.cloneTab());
			toRemoveNodesCage.forEach(n->{
				ctab.removeNodeAndEdges(n);	
			});
			ctab.standardCleanEdges();
			
			
			//System.out.println("MAde new nodes on triples:" + ctabRaw.size());
			ctabRaw.add(ctab.cloneTab());
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
							 .filter(l->GeomUtil.length(l)>ctab.getAverageBondLength()*0.4)
			                 .map(l->GeomUtil.findCenterOfShape(l))
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
			
			
			
			
			
			//Find floating methyls
			
			List<Shape> maybeDash = polygons.stream()
					.filter(s->GeomUtil.area(s)<averageAreaOCR)
			        .filter(s->!likelyOCR.contains(s))
			        .map(s->Tuple.of(s,GeomUtil.findCenterOfShape(s)))
			        .filter(st->!growLikelyOCR.stream().filter(g->g.contains(st.v())).findFirst().isPresent())
			        .map(t->t.k())
			        .collect(Collectors.toList());
			
			List<Shape> dashShapes = new ArrayList<Shape>();
			
			GeomUtil.groupThings(maybeDash, t->{
				Shape s1=t.k();
				Shape s2=t.v();
				Point2D p1=GeomUtil.findCenterOfShape(s1);
				Point2D p2=GeomUtil.findCenterOfShape(s2);
				
				return p1.distance(p2)<ctab.getAverageBondLength()/3;
			})
			.stream()
			.filter(sl->sl.size()>=3)
			.forEach(sl->{
				
				Shape bshape=sl.stream()
							   .reduce((s1,s2)->GeomUtil.add(s1, s2)).orElse(null);
				if(bshape!=null){
					Point2D[] pts=GeomUtil.getPairOfFarthestPoints(bshape);
					double dist=pts[0].distance(pts[1]);
					if(dist < ctab.getAverageBondLength()*1.3 && dist>ctab.getAverageBondLength()*0.6){
						//Looks very promising
						System.out.println("Found a possible dash:" + sl.size());
						dashShapes.add(bshape);
						List<Node> forN1=ctab.getNodes()
						    .stream()
						    .filter(n->n.getPoint().distance(pts[0])< ctab.getAverageBondLength()*0.3)
						    .collect(Collectors.toList());
						List<Node> forN2=ctab.getNodes()
							    .stream()
							    .filter(n->n.getPoint().distance(pts[1])< ctab.getAverageBondLength()*0.3)
							    .collect(Collectors.toList());
						
						
						if(forN1.size()+forN2.size()>1)return;
						if(forN1.size()+forN2.size()==0)return;
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
								 .filter(t->t.k().contains(op))
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
						
					}
				}
				
				
			});
			
			        
			
			

			List<Tuple<Edge,WedgeInfo>> winfo=(List<Tuple<Edge, WedgeInfo>>) ctab.getEdges()
					.stream()
			.map(e->{
				if(e.getRealNode1().isInvented() || e.getRealNode2().isInvented())return Optional.empty();
				Line2D useLine=GeomUtil.getLinesNotInside(e.getLine(), growLikelyOCR)
						.stream()
						.map(l->Tuple.of(l, GeomUtil.length(l)).withVComparator())
						.max(Comparator.naturalOrder())
						.map(t->t.k())
						.orElse(null);
				if(useLine!=null){
					
					Point2D c=GeomUtil.findCenterOfShape(useLine);
					
					Optional<Shape> isDash=dashShapes.stream()
							  .filter(d->d.contains(c))
							  .findFirst();
					
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
							e.setDashed(false);
						}
	
						int mult=1;
						if(e.getPoint1().distance(useLine.getP1())< e.getPoint2().distance(useLine.getP1())){
							mult=-1;
						}
						return bitmap.getconfexHullAlongLine(useLine)
								.map(w->Tuple.of(e,w));
					}
				}
				return Optional.empty();
			})
			.filter(t->t.isPresent())
			.map(o->o.get())
			.collect(Collectors.toList());
			
			Predicate<Node> couldBeStereoCenter = (n1)->n1.getEdgeCount()>=3 && n1.getSymbol().equals("C") && !n1.getEdges().stream().filter(e1->e1.getOrder()>1).findAny().isPresent();
			
			double averageThickness = winfo.stream().mapToDouble(t->t.v().getAverageThickness()).average().orElse(2);
			
			winfo.forEach(t->{
				Edge e=t.k();
				
				WedgeInfo s = t.v();
				
				//realRescueOCRCandidates.add(s.getHull());
				
				
				double wl=s.getCorrel();
				//TODO: we need something for thick bonds which are _not_ strictly wedges, but are meant to be
				//wedges.
				if(s.getAverageThickness()>averageThickness*1.9){
					if(s.getOnPixels()>s.getArea()*0.5){
						double cutoff=WEDGE_LIKE_PEARSON_SCORE_CUTOFF;
						if(e.getRealNode1().getEdges().stream().filter(ed->ed.getOrder()>1).findAny().isPresent()){
							cutoff=WEDGE_LIKE_PEARSON_SCORE_CUTOFF_DOUBLE;
						}
						if(e.getRealNode2().getEdges().stream().filter(ed->ed.getOrder()>1).findAny().isPresent()){
							cutoff=WEDGE_LIKE_PEARSON_SCORE_CUTOFF_DOUBLE;
						}
						
						
						if(wl>cutoff){
							e.setWedge(true);
						}else if(wl<-cutoff){
							e.setWedge(true);
							e.switchNodes();
						}else if(s.getAverageThickness()>averageThickness*2.5){
							//very thick line
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
				//System.out.println("Pixels on:" + s.getOnPixels() + " and area=" + s.getArea() + " and correl=" + s.getCorrel() + " and thickness:" + s.getAverageThickness());
			});
			
			

			GeomUtil.eachCombination(ctab.getNodes())
			.filter(t->t.k().distanceTo(t.v())<1.2*ctab.getAverageBondLength())
			.filter(t->!t.k().getBondTo(t.v()).isPresent())
			.forEach(t1->{
				Line2D l2 = new Line2D.Double(t1.k().getPoint(),t1.v().getPoint());
				Line2D useLine=GeomUtil.getLinesNotInside(l2, growLikelyOCR)
						.stream()
						.map(l->Tuple.of(l, GeomUtil.length(l)).withVComparator())
						.max(Comparator.naturalOrder())
						.map(t->t.k())
						.orElse(null);
				if(useLine==null)return;
				long c=polygons.stream()
						.filter(s->GeomUtil.getIntersection(s, useLine).isPresent())
						.count();
				if(c>2){
					ctab.addEdge(t1.k().getIndex(), t1.v().getIndex(), 1);
					Edge e=ctab.getEdges().get(ctab.getEdges().size()-1);
					e.setDashed(true);
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





			ctabRaw.add(ctab.cloneTab());
			
			
			
			
			
			/*
			 * long c=polygons.stream()
		        		        .filter(s->GeomUtil.getIntersection(s, useLine).isPresent())
		        		        .count();
			 */



			//ctab=ctabRaw;

		}
		
		rescueOCRShapes=realRescueOCRCandidates;
	
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
		
		ctab.getNodes().stream()
	    .filter(n->n.getSymbol().equals("S"))
	    .filter(n->n.getCharge()==0)
	    .forEach(n->{
	    		int so=n.getEdges().stream().mapToInt(e->e.getOrder()).sum();
	    		if(so==3){
	    			
	    			
	    			
	    			//n.setCharge(1);
	    			
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
		
		List<Shape> mightBeNegative=polygons.stream()
										    .filter(s->s.getBounds2D().getHeight()<ctab.getAverageBondLength()/10)
										    .filter(s->s.getBounds2D().getWidth()>1)
										    .filter(s->!likelyOCRAll.contains(s) || Optional.ofNullable(ocrAttempt.get(s)).map(l->l.get(0)).filter(t->t.k().toString().equals("-")).isPresent())
										    .collect(Collectors.toList());
		
		if(!mightBeNegative.isEmpty()){
			ctab.getNodes().stream()
		    .filter(n->n.getSymbol().equals("O"))
		    .filter(n->n.getCharge()==0)
		    .filter(n->n.getEdgeCount()==1)
		    .filter(n->!n.isInvented())
		    .forEach(n->{
		    		int so=n.getEdges().stream().mapToInt(e->e.getOrder()).sum();
		    		if(so==1){
			    		Point2D p=n.getPoint();
			    		Shape neg=mightBeNegative.stream().filter(s->GeomUtil.distanceTo(s, p)<ctab.getAverageBondLength()/2).findAny().orElse(null);
			    		
			    		if(neg!=null){
			    			Point2D ncenter=GeomUtil.findCenterOfShape(neg);	
			    			boolean inLine=n.getEdges().stream()
							    			.map(e->e.getLine())
							    			.map(l->GeomUtil.growLine(l, averageWidthOCRFinal[0]))
							    			.filter(s->s.contains(ncenter))
							    			.findAny()
							    			.isPresent();
			    			if(!inLine){
			    				n.setCharge(-1);
			    			}
			    		}
		    		}
		    });
			
			ctab.getNodes().stream()
		    .filter(n->n.getSymbol().equals("N"))
		    .filter(n->n.getCharge()==0)
		    .filter(n->!n.isInvented())
		    .forEach(n->{
		    		int so=n.getEdges().stream().mapToInt(e->e.getOrder()).sum();
		    		if(so==2){
			    		Point2D p=n.getPoint();
			    		Shape neg=mightBeNegative.stream().filter(s->GeomUtil.distanceTo(s, p)<ctab.getAverageBondLength()/2).findAny().orElse(null);
			    		
			    		if(neg!=null){
			    			Point2D ncenter=GeomUtil.findCenterOfShape(neg);
			    			
			    			boolean inLine=n.getEdges().stream()
							    			.map(e->e.getLine())
							    			.map(l->GeomUtil.growLine(l, averageWidthOCRFinal[0]))
							    			.filter(s->s.contains(ncenter))
							    			.findAny()
							    			.isPresent();
			    			if(!inLine){
			    				n.setCharge(-1);
			    			}
			    		}
		    		}
		    });
			
			ctab.getNodes().stream()
		    .filter(n->n.getSymbol().equals("Cl"))
		    .filter(n->n.getCharge()==0)
		    .filter(n->!n.isInvented())
		    .forEach(n->{
		    		int so=n.getValanceTotal();
		    		if(so<=1){
			    		Point2D p=n.getPoint();
			    		Shape neg=mightBeNegative.stream().filter(s->GeomUtil.distanceTo(s, p)<ctab.getAverageBondLength()/2).findAny().orElse(null);
			    		
			    		if(neg!=null){
			    			Point2D ncenter=GeomUtil.findCenterOfShape(neg);
			    			
			    			boolean inLine=n.getEdges().stream()
							    			.map(e->e.getLine())
							    			.map(l->GeomUtil.growLine(l, averageWidthOCRFinal[0]))
							    			.filter(s->s.contains(ncenter))
							    			.findAny()
							    			.isPresent();
			    			if(!inLine){
			    				n.setCharge(-1);
			    			}
			    		}
		    		}
		    });
			
			int sc=ctab.getSumCharge();
			
			if(sc>0){
				ctab.getNodes().stream()
			    .filter(n->n.getSymbol().equals("C"))
			    .filter(n->n.getCharge()==0)
			    .filter(n->n.getEdgeCount()==2)
			    .filter(n->!n.isInvented())
			    .filter(n->!n.getEdges().stream().filter(e->e.getDashed()).findAny().isPresent())
			    .forEach(n->{
			    		int so=n.getEdges().stream().mapToInt(e->e.getOrder()).sum();
			    		if(so==2){
				    		Point2D p=n.getPoint();
				    		Shape neg=mightBeNegative.stream().filter(s->GeomUtil.distanceTo(s, p)<ctab.getAverageBondLength()/2).findAny().orElse(null);
				    		
				    		if(neg!=null){
				    			Point2D ncenter=GeomUtil.findCenterOfShape(neg);
				    			boolean inLine=n.getEdges().stream()
								    			.map(e->e.getLine())
								    			.map(l->GeomUtil.growLine(l, averageWidthOCRFinal[0]))
								    			.filter(s->s.contains(ncenter))
								    			.findAny()
								    			.isPresent();
				    			if(!inLine){
				    				n.setCharge(-1);
				    			}
				    		}
			    		}
			    });
			}
		}
		 
		
		ctab.simpleClean();
		//Make aromatic bonds
		circles.stream()
		       .map(c->GeomUtil.findCenterOfShape(c))
		       .forEach(cp->{
		    	   ctab.getEdgesWithCenterWithin(cp,ctab.getAverageBondLength())
				       .stream()
				       .forEach(Edge::setToAromatic);
		       });
		

	}


	public Chemical getChemical(){
		return ctab.toChemical();
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
		return lines;
	}

	public List<Line2D> getLinesJoined() {
		return linesJoined;
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
