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
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.Predicate;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import gov.nih.ncats.chemkit.api.Chemical;
import tripod.molvec.Bitmap;
import tripod.molvec.CachedSupplier;
import tripod.molvec.algo.Tuple.KEqualityTuple;
import tripod.molvec.ui.RasterCosineSCOCR;
import tripod.molvec.ui.SCOCR;
import tripod.molvec.util.CompareUtil;
import tripod.molvec.util.ConnectionTable;
import tripod.molvec.util.ConnectionTable.Edge;
import tripod.molvec.util.ConnectionTable.Node;
import tripod.molvec.util.GeomUtil;

public class StructureImageExtractor {
	static final SCOCR OCR_DEFAULT=new RasterCosineSCOCR(RasterCosineSCOCR.SANS_SERIF_FONTS());
	static final SCOCR OCR_BACKUP=new RasterCosineSCOCR(RasterCosineSCOCR.SERIF_FONTS());
	
	
	static final SCOCR OCR_ALL=new RasterCosineSCOCR();
   
	static{
		Set<Character> alpha=SCOCR.SET_COMMON_CHEM_ALL();
		alpha.add(Character.valueOf('/'));
		alpha.add(Character.valueOf('\\'));
		OCR_DEFAULT.setAlphabet(alpha);
		OCR_BACKUP.setAlphabet(alpha);
    	
    	Set<Character> alphaAll=SCOCR.SET_COMMON_PUCTUATION();
    	alphaAll.addAll(SCOCR.SET_ALPHANUMERIC());
    	OCR_ALL.setAlphabet(alphaAll);
    }
    
	private Bitmap bitmap; // original bitmap
    private Bitmap thin; // thinned bitmap
    
    
    private List<Shape> polygons;
    private List<Line2D> lines;
    private List<Line2D> linesJoined;
    private List<Tuple<Line2D,Integer>> linesOrder;    
    private Map<Shape,List<Tuple<Character,Number>>> ocrAttmept = new HashMap<>();
    private Map<Shape,String> bestGuessOCR = new HashMap<>();
    private Map<Shape,ShapeInfo> shapeTypes = new HashMap<>();
    
    
    private ConnectionTable ctab;
    private ConnectionTable ctabRaw;
    
    
    private final double MAX_REPS = 10;
    private final double INITIAL_MAX_BOND_LENGTH=Double.MAX_VALUE;
    private final double MIN_BOND_TO_AVG_BOND_RATIO = 1/2.7;
    private final double MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL = 1.3;
    private final double MIN_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL = 0.5;
    

    private final double MAX_TOLERANCE_FOR_DASH_BONDS = 3.5;
    private final double MAX_TOLERANCE_FOR_SINGLE_BONDS = 0.4;
    
	private final double OCRcutoffCosine=0.65;
	private final double OCRcutoffCosineRescue=0.50;
	
	private final double MAX_BOND_TO_AVG_BOND_RATIO_TO_KEEP = 1.5;
	private final double MIN_DISTANCE_BEFORE_MERGING_NODES = 4.0;
	private final double maxRatioForIntersection = 1.2;
	private final double maxPerLineDistanceRatioForIntersection = 2;
	private final double minPerLineDistanceRatioForIntersection = 0.7;
	private final double OCR_TO_BOND_MAX_DISTANCE=2.0;
	private final double maxCandidateRatioForIntersection = 1.9;        
	private final double MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS_THIN = 1;
	private final double MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS_FULL = 0.5;
	private final double MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS = 6;
	private final double MIN_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS = 2;
	
	private final double MAX_ANGLE_FOR_JOINING_SEGMENTS=25 * Math.PI/180.0;
	private final double MIN_SIZE_FOR_ANGLE_COMPARE_JOINING_SEGMENTS=8.0;
	
	private final double MAX_DISTANCE_TO_MERGE_PARALLEL_LINES=2;
	private final double MAX_POINT_DISTANCE_TO_BE_PART_OF_MULTI_NODE= 1;
	
	private final double MAX_BOND_RATIO_FOR_OCR_CHAR_SPACING=0.3;
	private final double MAX_THETA_FOR_OCR_SEPERATION=45 * Math.PI/180.0;
	private final double MAX_BOND_RATIO_FOR_MERGING_TO_OCR=0.5;
	private final double MIN_LONGEST_WIDTH_RATIO_FOR_OCR_TO_AVERAGE=0.7;
	
	//For finding high order bonds
	private final double MIN_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS=.5;
	private final double MIN_BIGGER_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS=.25;
	private final double MAX_ANGLE_FOR_PARALLEL=10.0 * Math.PI/180.0;
	
	private final double MIN_ST_DEV_FOR_KEEPING_DASHED_LINES=0.07;
	
	
	
	//Parallel lines
	
	private final double MAX_DELTA_LENGTH_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC=7.0;
	private final double MAX_ANGLE_FOR_PARALLEL_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC=15.0 * Math.PI/180.0;
	
    
	
	
	
	
	public StructureImageExtractor(){
    	
    }
    
    private static Tuple<Character,Number> adjustConfidence(Tuple<Character,Number> tup){
    	String ch=tup.k()+"";
    	double invScore=1-tup.v().doubleValue();
    	if(ch.equals("K") || ch.equals("k") || ch.equals("f")){
    		invScore=invScore*3.5; // penalize "K"
    	}
    	if(ch.equals("R")||
    	   ch.equalsIgnoreCase("A")||
    	   ch.equalsIgnoreCase("Z")||
    	   ch.equalsIgnoreCase("-")||
    	   ch.equalsIgnoreCase("m")||
    	   ch.equals("n")){
    		invScore=invScore*3; // penalize
    	}
    	if(ch.equalsIgnoreCase("X")){
    		invScore=invScore*1.5; // penalize
    	}else if(ch.equalsIgnoreCase("N")){
    		invScore=invScore*(0.9); // promote
    	}
    	
    	return Tuple.of(tup.k(),Math.max(0,1-invScore));
    	
    }
    
    private static boolean OCRIsLikely(Tuple<Character,Number> tup){
    	String t=tup.k()+"";
     	if(!"I".equalsIgnoreCase(t) && 
     	   !"L".equalsIgnoreCase(t) &&
     	   !"1".equalsIgnoreCase(t) &&
     	   !"-".equalsIgnoreCase(t) &&
     	   !"/".equalsIgnoreCase(t) &&
     	   !"K".equalsIgnoreCase(t) &&
     	   !"Y".equalsIgnoreCase(t) &&
     	   !"W".equalsIgnoreCase(t) &&
     	   !"f".equals(t) &&
     	   !"\\".equalsIgnoreCase(t)){
     		return true;
     	}
     	return false;
    }
    
    public static enum ShapeType{
    	TEXT,
    	TEXT_NUMERIC,
    	TEXT_CHEMICAL,
    	LINE,
    	LINES,
    	SIMILAR_ANGLE_LINES,
    	NOISE
    }
    
    
    public static class ShapeInfo{
    	private Shape s;
    	private ShapeType t;
    	private double densityOriginal;
    	private double densityThin;
    	private double area;
    	private int lineCount=0;
    	private double lineLengthAverage=0;
    	
    	public ShapeInfo(Shape s, ShapeType t, double densityO, double densityT, int lcount, double lineLengthAverage){
    		this.s=s;
    		this.t=t;
    		this.densityOriginal=densityO;
    		this.densityThin=densityT;
    		this.lineCount=lcount;
    		this.area=GeomUtil.area(s);
    		this.lineLengthAverage=lineLengthAverage;
    		
    	}
    	
    	
    	public double avgPixelsPerLineUnit(){
    		Rectangle2D rect=s.getBounds2D();
    		double bboxarea= rect.getWidth()*rect.getHeight();
    		double pixelsOn=bboxarea*this.densityOriginal;
    		
    		double totalLineLength=(this.lineCount*this.lineLengthAverage);
    		return pixelsOn/totalLineLength;
    	}
    	
    	
    	public String toString(){
    		return "ShapeInfo:" + s + ", type=" + t + ", densityRaw=" + this.densityOriginal + ", densityThin=" + this.densityThin + ", lines=" + lineCount + ", area=" + this.area + ", lineDensity=" + this.lineCount/area + ", thickness=" + this.avgPixelsPerLineUnit();
    	}
    	
    }
    
    private static ShapeInfo computeShapeType(Shape s, List<Line2D> lines, double cosCutoff, Bitmap bitmap, Bitmap thin){
    	 List<Line2D> containedLines=lines.stream()
				  .filter(l->s.contains(l.getP1()))
				  .collect(Collectors.toList());
    	 Bitmap bmcrop=bitmap.crop(s);
    	 Bitmap thcrop=thin.crop(s);
    	 Supplier<ShapeType> stypeGetter=()->{
	    	 if(s.getBounds2D().getWidth()<=0 || s.getBounds2D().getHeight()<=0)return ShapeType.NOISE;
	    	 
//	    	 
//	    	 {
//	          	List<Tuple<Character,Number>> potential = OCR_DEFAULT.getNBestMatches(4,
//	          			bmcrop.createRaster(),
//	          			thcrop.createRaster()
//	                      )
//	          			.stream()
//	          			.map(Tuple::of)
//	          			.map(t->adjustConfidence(t))
//	          			.collect(Collectors.toList());
//	          	
//	              if(potential.stream().filter(e->e.v().doubleValue()>cosCutoff).findAny().isPresent()){
//	            	  	Tuple<Character,Number> tchar=potential.get(0);
//	             	 	if(OCRIsLikely(tchar)){
//	             	 			if(Character.isDigit(tchar.k())){
//	             	 				return ShapeType.TEXT_NUMERIC;
//	             	 			}
//		                 		return ShapeType.TEXT_CHEMICAL;
//		                }
//		                return ShapeType.TEXT;
//	              }
//	    	 }
//	    	 {         
//	              List<Tuple<Character,Number>> potential = OCR_ALL.getNBestMatches(4,
//	                      bitmap.crop(s).createRaster(),
//	                      thin.crop(s).createRaster()
//	                      )
//	          			.stream()
//	          			.map(Tuple::of)
//	          			.collect(Collectors.toList());
//	          	
//	              if(potential.stream().filter(e->e.v().doubleValue()>cosCutoff).findAny().isPresent()){
//	            	  	Tuple<Character,Number> tchar=potential.get(0);
//	             	 	if(OCRIsLikely(tchar)){
//	             	 			if(Character.isDigit(tchar.k())){
//	             	 				return ShapeType.TEXT_NUMERIC;
//	             	 			}
//		                 		return ShapeType.TEXT;
//		                }
//		                return ShapeType.TEXT;
//	              }
//	        
//	    	 }
	    	
	
	    	 if(containedLines.size()==0)return ShapeType.NOISE;
	    	 if(containedLines.size()==1)return ShapeType.LINE;
	    	 
	    	 
	    	 List<List<Line2D>> parLines=GeomUtil.groupMultipleBonds(containedLines, 3*Math.PI/180, Double.POSITIVE_INFINITY, 0, 0);
	    	 
	    	 int ngroups = parLines.size();
	    	 
	    	 Tuple<Double,Double> expected=REFACTOR_ME.get().get(containedLines.size());
	    	 if(expected!=null){
		    	 if(ngroups<expected.k()-3*expected.v()){
		    		 //3 sigma better than expected
		    		 return ShapeType.SIMILAR_ANGLE_LINES;
		    	 }
	    	 }
	            
	             
	    	return ShapeType.LINES;
    	};
    	
    	
    	double dense1=Optional.ofNullable(bmcrop).map(c->c.fractionPixelsOn()).orElse(0.0);
    	double dense2=Optional.ofNullable(thcrop).map(c->c.fractionPixelsOn()).orElse(0.0);
    	
    	double avg=containedLines.stream().mapToDouble(l->GeomUtil.length(l)).average().orElse(0);
   	 
    	return new ShapeInfo(s,stypeGetter.get(),dense1,dense2,containedLines.size(),avg);
    	
    }
    
    private static CachedSupplier<Map<Integer,Tuple<Double,Double>>> REFACTOR_ME=CachedSupplier.of(()->{
    	 String raw="1	1.0	0.0\n" +
    	    		"2	1.973	0.1620833119108796\n" + 
    	    		"3	2.897	0.31998593719099716\n" + 
    	    		"4	3.803	0.41975111673466614\n" + 
    	    		"5	4.677	0.544675132533149\n" + 
    	    		"6	5.523	0.636766048089879\n" + 
    	    		"7	6.318	0.7105462687256918\n" + 
    	    		"8	7.089	0.8396898236849115\n" + 
    	    		"9	7.855	0.9380698268252714\n" + 
    	    		"10	8.609	1.02670297554843\n" + 
    	    		"11	9.298	1.1194623709620626\n" + 
    	    		"12	9.98	1.145251064177632\n" + 
    	    		"13	10.656	1.2343678544096903\n" + 
    	    		"14	11.249	1.2275988758548042\n" + 
    	    		"15	11.87	1.4153091535067586\n" + 
    	    		"16	12.439	1.5589352135351824\n" + 
    	    		"17	13.037	1.553586495821839\n" + 
    	    		"18	13.65	1.586663165262238\n" + 
    	    		"19	14.106	1.668761217190771\n" + 
    	    		"20	14.573	1.7060688731701248\n" + 
    	    		"21	14.972	1.752488516367517\n" + 
    	    		"22	15.501	1.8176905677259902\n" + 
    	    		"23	15.915	1.7691170113929704\n" + 
    	    		"24	16.386	1.8566108908438494\n" + 
    	    		"25	16.824	1.9470552123655849\n" + 
    	    		"26	17.043	1.9796845708344548\n" + 
    	    		"27	17.412	2.0199643561211715\n" + 
    	    		"28	17.854	2.0944412142621722\n" + 
    	    		"29	18.267	2.0370839452511587\n" + 
    	    		"30	18.381	2.187198893562266\n" + 
    	    		"31	18.615	2.073348740564426\n" + 
    	    		"32	18.857	2.200579696352762\n" + 
    	    		"33	19.292	2.152379148756084\n" + 
    	    		"34	19.537	2.1933150708459688\n" + 
    	    		"35	19.86	2.1859551687992087\n" + 
    	    		"36	19.918	2.218394915248418\n" + 
    	    		"37	20.137	2.3000502168430907\n" + 
    	    		"38	20.355	2.2730981061098046\n" + 
    	    		"39	20.61	2.2781352023091217\n" + 
    	    		"40	20.897	2.3324645763655427\n" + 
    	    		"41	20.86	2.310930548502063\n" + 
    	    		"42	21.18	2.315080992103732\n" + 
    	    		"43	21.192	2.2367690984989896\n" + 
    	    		"44	21.368	2.4491990527517387\n" + 
    	    		"45	21.472	2.4280065897768788\n" + 
    	    		"46	21.51	2.337498663101203\n" + 
    	    		"47	21.704	2.383775157182396\n" + 
    	    		"48	21.815	2.3273966142451776\n" + 
    	    		"49	21.843	2.3787288622287286\n" + 
    	    		"50	21.94	2.3665164271561525\n" + 
    	    		"51	21.981	2.374160693803157\n" + 
    	    		"52	22.117	2.3705085952174865\n" + 
    	    		"53	22.038	2.395110853384454\n" + 
    	    		"54	22.171	2.3349002120005125\n" + 
    	    		"55	22.234	2.339496527033094\n" + 
    	    		"56	22.144	2.4692638579139436\n" + 
    	    		"57	22.307	2.402654989797762\n" + 
    	    		"58	22.282	2.4381296109928217\n" + 
    	    		"59	22.303	2.4444203811946785\n" + 
    	    		"60	22.256	2.3950916475158013\n" + 
    	    		"61	22.31	2.4285592436669217\n" + 
    	    		"62	22.242	2.4862493841125355\n" + 
    	    		"63	22.267	2.433045622260299\n" + 
    	    		"64	22.354	2.4985363715583677\n" + 
    	    		"65	22.301	2.2992170406466794\n" + 
    	    		"66	22.122	2.4110404393124636\n" + 
    	    		"67	22.068	2.4292747889030526\n" + 
    	    		"68	22.055	2.40540537124203\n" + 
    	    		"69	21.945	2.416604022176571\n" + 
    	    		"70	21.958	2.3942088463624294\n" + 
    	    		"71	21.925	2.374736827524257\n" + 
    	    		"72	21.869	2.4420972544106445\n" + 
    	    		"73	21.753	2.3987477983314514\n" + 
    	    		"74	21.725	2.4230920329199104\n" + 
    	    		"75	21.735	2.3817588039094257\n" + 
    	    		"76	21.443	2.427910830323043\n" + 
    	    		"77	21.493	2.3600743632352055\n" + 
    	    		"78	21.345	2.365581323903283\n" + 
    	    		"79	21.306	2.2935483426341774\n" + 
    	    		"80	21.056	2.364077832898063\n" + 
    	    		"81	21.08	2.4231384607570523\n" + 
    	    		"82	21.084	2.3162348758275746\n" + 
    	    		"83	20.935	2.3771358816862085\n" + 
    	    		"84	20.639	2.255810054060419\n" + 
    	    		"85	20.801	2.320215291734811\n" + 
    	    		"86	20.786	2.24860045361552\n" + 
    	    		"87	20.427	2.372060496699029\n" + 
    	    		"88	20.46	2.4783865719455487\n" + 
    	    		"89	20.293	2.282794559306647\n" + 
    	    		"90	20.188	2.393878860761352\n" + 
    	    		"91	20.05	2.3806511714234744\n" + 
    	    		"92	19.916	2.32571365391356\n" + 
    	    		"93	19.784	2.4243234107684666\n" + 
    	    		"94	19.655	2.246769903661688\n" + 
    	    		"95	19.613	2.3730214916852344\n" + 
    	    		"96	19.312	2.2668603838789703\n" + 
    	    		"97	19.249	2.27134299479406\n" + 
    	    		"98	19.188	2.338515768601971\n" + 
    	    		"99	19.109	2.350982560547805";
    	 return Arrays.stream(raw.split("\n"))
		    	       .map(s->s.split("\t"))
		    	       .map(s->Tuple.of(Integer.parseInt(s[0]),Tuple.of(Double.parseDouble(s[1]),Double.parseDouble(s[2]))))
		    	       .collect(Tuple.toMap());
    });
    		
    
    
    private void processOCR(SCOCR socr, List<Shape> polygons,Bitmap bitmap, Bitmap thin, BiConsumer<Shape,List<Tuple<Character,Number>>> onFind){
    	 /*
	     * Looks at each polygon, and gets the likely OCR chars.
	     */   
	    for (Shape s : polygons) {
	         if(s.getBounds2D().getWidth()>0 && s.getBounds2D().getHeight()>0){
	         	List<Tuple<Character,Number>> potential = socr.getNBestMatches(4,
	                     bitmap.crop(s).createRaster(),
	                     thin.crop(s).createRaster()
	                     )
	         			.stream()
	         			.map(Tuple::of)
	         			.map(t->adjustConfidence(t))
	         			.map(t->t.withVComparator())
	         			.sorted(Comparator.reverseOrder())
	         			.collect(Collectors.toList());
	         	 onFind.accept(s, potential);
	         }
	    }
    }
    
    
     
    public StructureImageExtractor load(File file) throws IOException{
    	SCOCR[] socr=new SCOCR[]{OCR_DEFAULT};
    	
    	double[] maxBondLength=new double[]{INITIAL_MAX_BOND_LENGTH};    
        
        
        bitmap = Bitmap.read(file).clean();

        polygons = bitmap.connectedComponents(Bitmap.Bbox.Polygon);

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

        List<Shape> likelyOCR=new ArrayList<Shape>();
        List<Shape> likelyOCRAll=new ArrayList<Shape>();
        /*
         * Looks at each polygon, and gets the likely OCR chars.
         */   
        
    	processOCR(socr[0],polygons,bitmap,thin,(s,potential)->{
    		 ocrAttmept.put(s, potential);
             if(potential.stream().filter(e->e.v().doubleValue()>OCRcutoffCosine).findAny().isPresent()){
            	 	if(OCRIsLikely(potential.get(0))){
                 		likelyOCR.add(s);
                 	}
                 	likelyOCRAll.add(s);
             }
    	});
        
        if(likelyOCR.isEmpty()){
        	socr[0]=OCR_BACKUP;
        	likelyOCRAll.clear();
        	processOCR(socr[0],polygons,bitmap,thin,(s,potential)->{
       		 ocrAttmept.put(s, potential);
                if(potential.stream().filter(e->e.v().doubleValue()>OCRcutoffCosine).findAny().isPresent()){
               	 	if(OCRIsLikely(potential.get(0))){
                    		likelyOCR.add(s);
                    	}
                    	likelyOCRAll.add(s);
                }
        	});
        }
        
        double averageLargestOCR=likelyOCR.stream()
							              .map(s->GeomUtil.getPairOfFarthestPoints(s))
							              .mapToDouble(p->p[0].distance(p[1]))
							              .average()
							              .orElse(0);
        //System.out.println("avg ocr:" + averageLargestOCR);
        
                		 
        likelyOCRAll.retainAll(likelyOCRAll.stream()
                    .map(s->Tuple.of(s,GeomUtil.getPairOfFarthestPoints(s)))
                    .filter(t->t.v()[0].distance(t.v()[1]) > averageLargestOCR*MIN_LONGEST_WIDTH_RATIO_FOR_OCR_TO_AVERAGE)
                    .map(t->t.k())
                    .collect(Collectors.toList()));
        
        lines= GeomUtil.asLines(thin.segments());
        
        
        Predicate<Line2D> isInOCRShape = (l)->{
        	 if(likelyOCR.isEmpty())return false;
        	 Tuple<Shape,Double> shape1=GeomUtil.findClosestShapeTo(likelyOCR, l.getP1());
	       	  if(shape1.v()>OCR_TO_BOND_MAX_DISTANCE){
	       		  return false;
	       	  }
	       	  Tuple<Shape,Double> shape2=GeomUtil.findClosestShapeTo(likelyOCR, l.getP2());
	       	  if(shape2.v()>OCR_TO_BOND_MAX_DISTANCE){
	       		  return false;
	       	  }
	       	  if(shape1.k()==shape2.k()){
	       		  return true;
	       	  }
	       	  return false;
        };
        
        Predicate<Line2D> tryToMerge = isInOCRShape.negate().and((l)->{
        	return true;
        	
        	//return LineUtil.length(l)<largestBond;
        });
        
        
        
        List<Line2D> smallLines=lines.stream()
        						     .filter(tryToMerge)
        						     .collect(Collectors.toList());
        
        List<Line2D> bigLines=lines.stream()
        		                   .filter(tryToMerge.negate())
        		                   .collect(Collectors.toList());
        
        smallLines=smallLines.stream()
			                  .flatMap(l->{
			                	  double len=GeomUtil.length(l);
			                	  if(len<MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS && len>MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS){
			                		  return GeomUtil.splitLineIn2(l).stream();  
			                	  }
			                	  return Stream.of(l);
			                  })
			                  .collect(Collectors.toList());
        
       // smallLines= thin.combineLines(smallLines, MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS, MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS_THIN, MAX_POINT_DISTANCE_TO_BE_PART_OF_MULTI_NODE,MAX_ANGLE_FOR_JOINING_SEGMENTS,MIN_SIZE_FOR_ANGLE_COMPARE_JOINING_SEGMENTS);
        
        smallLines= bitmap.combineLines(smallLines, MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS, MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS_FULL, MAX_POINT_DISTANCE_TO_BE_PART_OF_MULTI_NODE,MAX_ANGLE_FOR_JOINING_SEGMENTS,MIN_SIZE_FOR_ANGLE_COMPARE_JOINING_SEGMENTS);
        smallLines= thin.combineLines(smallLines, MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS, MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS_THIN, MAX_POINT_DISTANCE_TO_BE_PART_OF_MULTI_NODE,MAX_ANGLE_FOR_JOINING_SEGMENTS,MIN_SIZE_FOR_ANGLE_COMPARE_JOINING_SEGMENTS);
        
        smallLines=smallLines.stream()
                .filter(l->GeomUtil.length(l)>MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS)
                .collect(Collectors.toList());
        
        linesJoined=Stream.concat(bigLines.stream(),
        		            smallLines.stream())
        		    .collect(Collectors.toList());
        shapeTypes=new HashMap<>();
        for(Shape s: polygons){
        	
        	ShapeInfo st=computeShapeType(s, linesJoined, OCRcutoffCosine, bitmap, thin);
        	shapeTypes.put(s, st);
        	
        	if(st.t.equals(ShapeType.NOISE)){
        		likelyOCR.remove(s);
        		likelyOCRAll.remove(s);
        	}
        	//System.out.println(st);
        }
//        
//        GeomUtil.groupThings(shapeTypes.values().stream().collect(Collectors.toList()), (tm)->{
//        	ShapeInfo st1=tm.k();
//        	ShapeInfo st2=tm.v();
//        	
//        	double my1=st1.s.getBounds2D().getMinY();
//        	double my2=st2.s.getBounds2D().getMinY();
//        	
//        	
//        	
//        	
//        	return true;
//        });
        
        
        
        
        
        
       // if(true)return this;
        
        
        double largestBond=linesJoined.stream()
		           .mapToDouble(l->GeomUtil.length(l))
		           .max()
		           .orElse(0);
       

        
        List<Line2D> preprocess= GeomUtil.reduceMultiBonds(linesJoined, MAX_ANGLE_FOR_PARALLEL, MAX_DISTANCE_TO_MERGE_PARALLEL_LINES, MIN_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS,0,MAX_DELTA_LENGTH_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC)
        		                         .stream()
        		                         .map(t->t.k())
        		                         .collect(Collectors.toList());
        
        preprocess=GeomUtil.stitchEdgesInMultiBonds(preprocess, MAX_ANGLE_FOR_PARALLEL_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC, largestBond/3, MIN_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS, 0,MAX_DELTA_LENGTH_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC);
        
        linesOrder=GeomUtil.reduceMultiBonds(preprocess, MAX_ANGLE_FOR_PARALLEL, largestBond/3, MIN_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS, MIN_BIGGER_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS,MAX_DELTA_LENGTH_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC);
        
              
        
        
        
     
        int reps=0;
        boolean tooLongBond=true;
        
        while(tooLongBond){
        	List<Tuple<Line2D,Integer>> linesOrderRestricted =linesOrder.stream()
        	          .filter(t->{
        	        	  Line2D l=t.k();
        	        	  return isInOCRShape.negate().test(l);
        	          })
        	          .collect(Collectors.toList());
        	
	        ctab = GeomUtil.getConnectionTable(linesOrderRestricted, likelyOCR, 
	        		maxRatioForIntersection, 
	        		maxCandidateRatioForIntersection,
	        		maxPerLineDistanceRatioForIntersection,
	        		minPerLineDistanceRatioForIntersection,
	        		l-> (GeomUtil.length(l) < maxBondLength[0]))
	        		       .mergeNodesCloserThan(MIN_DISTANCE_BEFORE_MERGING_NODES);
	        
	        
//	        ctab.getTolerancesForAllEdges(bitmap)
//	    	    .forEach(t->{
//	    	    	//System.out.println("Score for edge:" + t.v());
//	    	    	if(t.v()>MAX_TOLERANCE_FOR_DASH_BONDS && t.k().getOrder()==1){	    	    		
//	    	    		ctab.removeEdge(t.k());
//	    	    	}
//	    	    });
	        
	        ctabRaw=ctab.cloneTab();
	        
	        for(Shape s: likelyOCR){
	        	ctab.mergeAllNodesInsideCenter(s, OCR_TO_BOND_MAX_DISTANCE);
	        }
	        //double bl1=ctab.getAverageBondLength();
	        ctab.mergeNodesCloserThan(ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO);
	        ctab.mergeAllNodesOnParLines();
	        ctab.cleanMeaninglessEdges();
	        ctab.cleanDuplicateEdges((e1,e2)->{
	        	if(e1.getOrder()>e2.getOrder()){
	        		return e1;
	        	}
	        	return e2;
	        });
	        
	        
	        ctab.createNodesOnIntersectingLines();
	        ctab.mergeNodesCloserThan(ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO);
	        ctab.cleanMeaninglessEdges();
	        ctab.cleanDuplicateEdges((e1,e2)->{
	        	if(e1.getOrder()>e2.getOrder()){
	        		return e1;
	        	}
	        	return e2;
	        });
	        ctab.mergeNodesExtendingTo(likelyOCR,0.5);
	        
	        ctab.removeOrphanNodes();
	        ctab.mergeNodesCloserThan(ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO);
	        
	        
	        ctab.makeMissingNodesForShapes(likelyOCR,MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL,MIN_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL);
	        
	        Set<Node> toRemove = new HashSet<Node>();
	        
	        ctab.makeMissingBondsToNeighbors(bitmap,MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL,MAX_TOLERANCE_FOR_DASH_BONDS,likelyOCR,OCR_TO_BOND_MAX_DISTANCE, (t)->{
//	        	System.out.println("Tol found:" + t.k());
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
		        		double ddelta=Math.abs(sumd-e.getBondLength());
		        		if(ddelta<MAX_DELTA_LENGTH_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC){
		        			toRemove.add(cn);
		        		}
		        		alreadyExists=true;
	        		}
	        	}
	        	
	        	
	        	//System.out.println("Tol found for add:" + t.k());
	        	if(t.k()>MAX_TOLERANCE_FOR_SINGLE_BONDS){
	        		t.v().setDashed(true);
	        	}
	        });
	        toRemove.forEach(n->ctab.removeNodeAndEdges(n));
	        ctab.removeOrphanNodes();
	        
	     	        
	        double avgBondLength=ctab.getAverageBondLength();
	        maxBondLength[0]=avgBondLength*MAX_BOND_TO_AVG_BOND_RATIO_TO_KEEP;
	        
	        
	        
	        //System.out.println("Average bond length:" + avgBondLength);
	        
	        tooLongBond = ctab.getEdges()
	        		          .stream()
	        	//	          .peek(e->System.out.println(e.getBondDistance()))
	        		          .filter(e->e.getBondLength()>maxBondLength[0])
	        		          .findAny()
	        		          .isPresent();
	        if(tooLongBond){
	        	//System.out.println("No good, try again");
	        	reps++;
	        }
	        if(reps>MAX_REPS)break;
        }
        
        ctab.getEdges()
        .forEach(e->{
        	double wl=bitmap.getWedgeLikeScore(e.getLine());
        	//System.out.println("WL:" + wl);
        	if(wl>.8){
        		e.setWedge(true);            		
        	}
        	else if(wl<-.8){
        		e.setWedge(true);
        		e.switchNodes();
        	}            	
        });
        
        double shortestRealBondRatio = .3;
        ctab.fixBondOrders(likelyOCR,shortestRealBondRatio, e->{
        	e.setOrder(1);
        });
        
        
        double avgTol=ctab.getTolerancesForAllEdges(bitmap,likelyOCR)
				          .stream()
				          .mapToDouble(t->t.v())
				          .average()
				          .orElse(0);

        
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
        double SEED_BOND_RATIO_FOR_OCR_WIDTH=0.20;
        double PROBLEM_BOND_LENGTH_RATIO=0.8;
        double PROBLEM_TOLERANCE_RATIO=1.0;
        
        double SEED_BOND_RATIO_FOR_OCR_WIDTH_FOR_CENTROID=0.5;
        
        List<Node> unmatchedNodes=ctab.getNodesNotInShapes(likelyOCR, OCR_TO_BOND_MAX_DISTANCE + avgbond*SEED_BOND_RATIO_FOR_OCR_WIDTH);
        
        
        
        
        
        //List<Point2D> vertices=lines.stream().flatMap(l->Stream.of(l.getP1(),l.getP2())).collect(Collectors.toList());
        List<Point2D> verticesJ=linesJoined.stream().flatMap(l->Stream.of(l.getP1(),l.getP2())).collect(Collectors.toList());
        
        
        List<Shape> toAddAllOCR=new ArrayList<Shape>();
        
        
        unmatchedNodes.forEach(n->{
        	Tuple<Edge,Double> worstEdge=ctab.getWorstToleranceForNode(n,bitmap,likelyOCR);
        	
//        	System.out.println("Avg tol:" + avgTol);
//        	System.out.println("Worst tol here:" + worstEdge.v());
        	
        	
        	//The worst bond isn't so bad, maybe it's not worth looking
        	if(worstEdge.v()<avgTol*PROBLEM_TOLERANCE_RATIO){
        		
        		//double check to see if the shortest bond is too short
        		//if it isn't, then skip this one, it's probably fine.
        		boolean tooSmallBond=n.getEdges()
        				              .stream()
						              .mapToDouble(e->e.getBondLength())
						              .filter(d->d<avgbond*PROBLEM_BOND_LENGTH_RATIO)
						              .findAny()
						              .isPresent();
        		if(tooSmallBond){
        			//continue
        		}else{
        			return;
        		}
        	}
        	
        	
        	
        	Point2D cpt=n.getPoint();
        	
        	
        	int numEdges=n.getEdges().size();
        	Shape[] area=new Shape[]{null};
        	Shape nshape=null;
        	double radius=Math.max(avgbond*SEED_BOND_RATIO_FOR_OCR_WIDTH_FOR_CENTROID,averageLargestOCR/2);
        	
        	boolean keep=true;
        	for(int i=0;i<3;i++){
        		keep=true;
        		area[0]=GeomUtil.makeShapeAround(cpt,radius);
        		
        		//
            	List<Point2D> insideVertices=verticesJ.stream()
    							        	        .filter(v->area[0].contains(v))
    							        	        .collect(Collectors.toList());
	        	if(insideVertices.size()<=numEdges+1){
	        		keep=false;
	        	}
	        	List<Point2D> insideVertices2=verticesJ.stream()
	        	        .filter(v->area[0].contains(v))
	        	        .collect(Collectors.toList());
        		nshape = GeomUtil.convexHull(insideVertices2.toArray(new Point2D[0]));
        		
        		Point2D[] far=GeomUtil.getPairOfFarthestPoints(nshape);
        		double r=0;
        		if(far!=null){
        			r=far[0].distance(far[1]);	
        		}
        		if(r < averageLargestOCR*MIN_LONGEST_WIDTH_RATIO_FOR_OCR_TO_AVERAGE){
        			keep=false;
        			//break;
                }
            	radius=Math.max(averageLargestOCR/2,r/2);
            	//cpt=GeomUtil.findCenterOfVertices(Arrays.asList(GeomUtil.vertices(nshape)));
            	cpt=GeomUtil.findCenterOfShape(nshape);
        	}
        	
        	if(keep){
        		
        		Bitmap nmap=bitmap.crop(nshape);
                Bitmap nthinmap=thin.crop(nshape);
                if(nmap!=null && nthinmap!=null){
                	//System.out.println("And it looks promising");
                	nshape=GeomUtil.growShape(nshape, 2);
                	nmap=bitmap.crop(nshape);
                    nthinmap=thin.crop(nshape);
                    
                    List<Shape> slist=nmap.connectedComponents(Bitmap.Bbox.Polygon);
                    
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
                   
                    if(nmap!=null && nthinmap!=null){
	                    
	                	List<Tuple<Character,Number>> potential = socr[0].getNBestMatches(4,
		                		nmap.createRaster(),
		                		nthinmap.createRaster()
		                        )
		            			.stream()
		            			.map(Tuple::of)
		            			.map(t->adjustConfidence(t))
		            			.collect(Collectors.toList());
	                	ocrAttmept.put(nshape, potential);
		                polygons.add(nshape);
		                if(potential.stream().findFirst().filter(e->e.v().doubleValue()>OCRcutoffCosineRescue).isPresent()){
		                	toAddAllOCR.add(nshape);
//		                	System.out.println("Found another pot:" +potential.get(0).k());
//		                	System.out.println("Found another pot score:" +potential.get(0).v());
//		                	
//			                if(OCRIsLikely(potential.get(0))){
//			                	likelyOCR.add(nshape);
//			                }
//			                likelyOCRAll.add(nshape);
		                }
                    }
                }
        	}
        	
        	//}
        	
        });
       
        
        
        GeomUtil.mergeOverlappingShapes(toAddAllOCR)
                .forEach(nshape->{
                	boolean sigOverlap = 
                	likelyOCR.stream()
								                    .map(s->Tuple.of(s,GeomUtil.getIntersectionShape(nshape, s)))
								                    .filter(os->os.v().isPresent())
//								                    .peek(t->{
//								                    	System.out.println("Overlap:" + ocrAttmept.get(t.k()).stream().findFirst().get().k());
//								                    })
								                    .map(Tuple.vmap(os->os.get()))
								                    .map(Tuple.vmap(s->GeomUtil.area(s)))
								                    .map(Tuple.kmap(s->GeomUtil.area(s)))
								                    .mapToDouble(t->t.v()/t.k())
//								                    .peek(area->System.out.println(area))
								                    .filter(areaFraction->areaFraction>0.5)
								                    .findAny()
								                    .isPresent();
                	
                	if(sigOverlap){
                		
                		return;
                	}
                	//if(ctab.getNodesInsideShape(nshape, 0).isEmpty())return;
                	Bitmap nmap=bitmap.crop(nshape);
                	Bitmap nthinmap=thin.crop(nshape);
                	if(nmap!=null && nthinmap!=null){
	                    
	                	List<Tuple<Character,Number>> potential = socr[0].getNBestMatches(4,
		                		nmap.createRaster(),
		                		nthinmap.createRaster()
		                        )
		            			.stream()
		            			.map(Tuple::of)
		            			.map(t->adjustConfidence(t))
		            			.collect(Collectors.toList());
	                	ocrAttmept.put(nshape, potential);
		                polygons.add(nshape);
		                if(potential.stream().findFirst().filter(e->e.v().doubleValue()>OCRcutoffCosineRescue).isPresent()){
		                	//toAddAllOCR.add(nshape);
//		                	
			                if(OCRIsLikely(potential.get(0))){
			                	likelyOCR.add(nshape);
			                }
			                likelyOCRAll.add(nshape);
		                }
                    }
                });
        
        
        ctab.mergeNodesExtendingTo(likelyOCR,0.5);
        
        double cosThetaOCRShape =Math.cos(MAX_THETA_FOR_OCR_SEPERATION);
        
        List<List<Shape>> ocrGroupList=GeomUtil.groupShapesIfClosestPointsMatchCriteria(likelyOCRAll, t->{
        	Point2D[] pts=t.v();
        	Shape[] shapes =t.k();
        	Line2D l2 = new Line2D.Double(pts[0],pts[1]);
        	double dist=GeomUtil.length(l2);
        	double cutoff=ctab.getAverageBondLength()*MAX_BOND_RATIO_FOR_OCR_CHAR_SPACING;
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
	        	
	        	for(Shape s: sorted){
	        		String v=(ocrAttmept.get(s).get(0).k() + "");
	        		if(v.equals("-")){
	        			if(making!=null){
	        				bestGuessOCR.put(making, soFar);
	        			}
	        			soFar="";
	        			making=null;
	        			continue;
	        		}
	        		if(making==null){
	        			making=s;
	        		}else{
	        			making=GeomUtil.add(making,s);
	        		}
	        		soFar+=v;
	        	}
	        	
	        	if(making!=null){
    				bestGuessOCR.put(making, soFar);
    			}
	        	
	        	
	        });
        

        
        List<Shape> ocrMeaningful=bestGuessOCR.keySet()
				   .stream()
//				   .peek(t->System.out.println(bestGuessOCR.get(t)))
				   .filter(s->BranchNode.interpretOCRStringAsAtom(bestGuessOCR.get(s))!=null)
				   .collect(Collectors.toList());
        
        for(Shape s: bestGuessOCR.keySet()){
        	String sym=bestGuessOCR.get(s);
        	BranchNode actual=BranchNode.interpretOCRStringAsAtom(sym);
        	if(actual!=null && actual.isRealNode()){
        		Point2D center = GeomUtil.findCenterOfShape(s);
        		ctab.mergeAllNodesInside(s, MAX_BOND_RATIO_FOR_MERGING_TO_OCR*ctab.getAverageBondLength(),(n)->{
        			if(sym.equals("H")){
        				if(GeomUtil.findClosestShapeTo(ocrMeaningful, n.getPoint()).k() !=s){
        					return false;
        				}
        			}
        			
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
        	}
        }
        ctab.cleanMeaninglessEdges();
        ctab.cleanDuplicateEdges((e1,e2)->{
        	if(e1.getOrder()>e2.getOrder()){
        		return e1;
        	}else if(e1.getOrder()>e2.getOrder()){
        		return e2;
        	}else{
        		if(e1.getDashed())return e2;
        		if(e2.getDashed())return e1;
        		if(e1.getWedge())return e1;
        		if(e2.getWedge())return e2;
        	}
        	return e1;
        });
        
        
        
        
        ctab.makeMissingBondsToNeighbors(bitmap,MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL,MAX_TOLERANCE_FOR_SINGLE_BONDS,likelyOCR,OCR_TO_BOND_MAX_DISTANCE, (t)->{
//        	System.out.println("Tol found for add:" + t.k());
        	if(t.k()>MAX_TOLERANCE_FOR_SINGLE_BONDS){
        		t.v().setDashed(true);
        	}
        });
        
        List<Line2D> lj =linesJoined.stream()
                   .flatMap(l->GeomUtil.getLinesNotInside(l, likelyOCR).stream())
                   .collect(Collectors.toList());
                   
        
        GeomUtil.groupShapesIfClosestPointsMatchCriteria(likelyOCR, (t)->{
        	
        	Point2D[] pts=t.v();
        	if(pts[0].distance(pts[1])<ctab.getAverageBondLength()){
        		//It should also have at least 1 line segment between the two
        		Shape cshape = GeomUtil.add(t.k()[0], t.k()[1]);
        		boolean containsLine=lj.stream()
        							   .filter(l1->cshape.contains(l1.getP1()))
        							   .findAny()
        							   .isPresent();
        		if(containsLine){
	        		return true;
        		}else{
        			return false;
        		}
        	}
        	return false;
        })
        .stream()
        .filter(ls->ls.size()>=2)
        .forEach(ls->{
        	List<Node> nodes = ls.stream()
        	  .map(s->ctab.getNodesInsideShape(s, 0))
        	  .flatMap(nds->nds.stream())
        	  .distinct()
        	  .collect(Collectors.toList());
        	//System.out.println("nodes:" + nodes.size());
        	if(nodes.size()==2){
        		Node n1=nodes.get(0);
        		Node n2=nodes.get(1);
        		boolean already=ctab.getEdgeBetweenNodes(n1, n2).isPresent();
        		if(!already){
        			ctab.addEdge(nodes.get(0).getIndex(), nodes.get(1).getIndex(), 1);
        		}
        	}
        	          
        });
        
        //remove triangles that are obviously wrong
        {
        	double avgL = ctab.getAverageBondLength();
        	Set<Edge> skip= new HashSet<Edge>();
        	ctab.getEdges()
        	    .stream()
        	    .filter(e->e.getBondLength()>avgL)
        	    .map(l->Tuple.of(l,GeomUtil.length(l.getLine())).withVComparator())
        	    .sorted(Comparator.reverseOrder())
        	    .map(t->t.k())
        	    .filter(t->!skip.contains(t))
        	    .forEach(e->{
        	    	Node n1= e.getRealNode1();
        	    	Node n2= e.getRealNode2();
        	    	List<KEqualityTuple<Node,Edge>> neigh1=n1.getNeighborNodes();
        	    	List<KEqualityTuple<Node,Edge>> neigh2=n2.getNeighborNodes();
        	    	List<KEqualityTuple<Node,Edge>> things=neigh1.stream()
        	    												 .filter(ne->neigh2.contains(ne))
        	    												 .collect(Collectors.toList());
        	    	
        	    	if(things.size()>0){
        	    		//System.out.println("Triangle found");
        	    		Point2D p1=n1.getPoint();
	        	    	Point2D p2=n2.getPoint();
	        	    	Point2D p3=things.get(0).k().getPoint();
	        	    	
	        	    	double tarea=Math.abs(GeomUtil.areaTriangle(p1,p2,p3));
	        	    	
	        	    	double expected = Math.sqrt(3)/4*Math.pow(e.getBondLength(),2);
	        	    	if(tarea<expected*0.5){
	        	    		//System.out.println("It's a bad one");
	        	    		ctab.removeEdge(e);
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
        
        
        //final cleanup
        {
        	ctab.getDashLikeScoreForAllEdges(bitmap)
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
        
        double fbondlength=ctab.getAverageBondLength();
        
        for(Shape s: bestGuessOCR.keySet()){
        	String sym=bestGuessOCR.get(s);
        	BranchNode actual=BranchNode.interpretOCRStringAsAtom(sym);
        	if(actual!=null && actual.isRealNode()){
        		//System.out.println(actual.toString());
        		List<Node> nlist=ctab.setNodeToSymbol(s, actual.getSymbol());
        		
        		if(nlist.size()==1){
        			Node pnode=nlist.get(0);
        			Point2D ppoint=pnode.getPoint();
        			if(actual.hasChildren()){
        				actual.generateCoordinates();
        				AffineTransform at = new AffineTransform();
        				at.translate(ppoint.getX(), ppoint.getY());
        				at.scale(fbondlength, fbondlength);
        				if(pnode.getEdges().size()>0){
        					Edge edge1= pnode.getEdges().get(0);
        					Point2D otherPoint = edge1.getPoint2();
        					if(!edge1.getRealNode1().equals(pnode)){
        						otherPoint = edge1.getPoint1();
        					}
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
	        					         .setSymbol(curN.getSymbol());
        					 ctab.addEdge(mpnode.getIndex(), n.getIndex(), curN.getOrderToParent());
        					 parentNodes.put(curN, n);
        				});
        				
        				
        			}
        		}
        	}
        	
        }
        
        
        //ctab=ctabRaw;
        

        return this;
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
		return ocrAttmept;
	}

	public ConnectionTable getCtab() {
		return ctab;
	}



	public ConnectionTable getCtabRaw() {
		return this.ctabRaw;
	}

	public Map<Shape,ShapeInfo> getShapeTypes() {
		return shapeTypes;
	}

	
	
	
}
