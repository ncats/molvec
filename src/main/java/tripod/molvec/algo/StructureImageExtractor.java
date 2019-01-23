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
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import gov.nih.ncats.chemkit.api.Chemical;
import tripod.molvec.Bitmap;
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
    private Map<Shape,List<Tuple<Character,Number>>> ocrAttmept = new HashMap<>();
    private Map<Shape,String> bestGuessOCR = new HashMap<>();
    private Map<Shape,ShapeInfo> shapeTypes = new HashMap<>();
    
    
    private ConnectionTable ctab;
    private List<ConnectionTable> ctabRaw = new ArrayList<ConnectionTable>();
    
    
    private final double MAX_REPS = 10;
    private final double INITIAL_MAX_BOND_LENGTH=Double.MAX_VALUE;
    private final double MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE = 1/3.0;
    private final double MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE_INITIAL = 1/2.7;
    
    private final double MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL = 1.3;
    private final double MIN_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL = 0.5;
    

    private final double MAX_TOLERANCE_FOR_DASH_BONDS = 3.0;
    private final double MAX_TOLERANCE_FOR_SINGLE_BONDS = 0.4;
    
	private final double OCRcutoffCosine=0.65;
	private final double OCRcutoffCosineRescue=0.50;
	
	private final double WEDGE_LIKE_PEARSON_SCORE_CUTOFF=.80;
	
	private final double MAX_BOND_TO_AVG_BOND_RATIO_TO_KEEP = 1.8;
	private final double MAX_DISTANCE_BEFORE_MERGING_NODES = 4.0;
	private final double maxRatioForIntersection = 1.2;
	private final double maxPerLineDistanceRatioForIntersection = 2;
	private final double minPerLineDistanceRatioForIntersection = 0.7;
	private final double OCR_TO_BOND_MAX_DISTANCE=3.0;
	private final double maxCandidateRatioForIntersection = 1.8;        
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
	private final double MAX_BOND_RATIO_FOR_MERGING_TO_OCR=0.3165;
	
	
	private final double MIN_LONGEST_WIDTH_RATIO_FOR_OCR_TO_AVERAGE=0.5;
	private final double MIN_AREA_RATIO_FOR_OCR_TO_AVERAGE=0.6;
	
	
	
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
    	   ch.equalsIgnoreCase("W")||
    	   
    	   ch.equals("n")){
    		invScore=invScore*3; // penalize
    	}
    	if(ch.equalsIgnoreCase("X") || ch.equalsIgnoreCase("+") || ch.equals("h")||ch.equalsIgnoreCase("D")){
    		invScore=invScore*1.5; // penalize
    	}else if(ch.equalsIgnoreCase("N") 
    		  || ch.equalsIgnoreCase("C") 
    		  || ch.equalsIgnoreCase("O")){
    		invScore=invScore*(0.85); // promote
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
     	   !"t".equalsIgnoreCase(t) &&
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
    	
    	public double getLineDensity(){
    		return this.lineCount/area;
    	}
    	
    	public String toString(){
    		return "ShapeInfo:" + s + ", type=" + t + ", densityRaw=" + this.densityOriginal + ", densityThin=" + this.densityThin + ", lines=" + lineCount + ", area=" + this.area + ", lineDensity=" + this.lineCount/area + ", thickness=" + this.avgPixelsPerLineUnit();
    	}
    	
    }
    
    private static ShapeInfo computeShapeType(Shape s, List<Line2D> lines, double cosCutoff, Bitmap bitmap, Bitmap thin){
    	 List<Line2D> containedLines=lines.stream()
				  .filter(l->s.contains(l.getP1())||s.contains(l.getP2()))
				  .collect(Collectors.toList());
    	 Bitmap bmcrop=bitmap.crop(s);
    	 Bitmap thcrop=thin.crop(s);
    	 Supplier<ShapeType> stypeGetter=()->{
	    	 if(s.getBounds2D().getWidth()<=0 || s.getBounds2D().getHeight()<=0)return ShapeType.NOISE;
	    	 	
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
    	
	    polygons.stream()
	    		.parallel()
	    	    .flatMap(s->{
	    	    	List<Tuple<Shape,List<Tuple<Character,Number>>>> got= new ArrayList<>();
	    	    	 Rectangle2D bounds2d = s.getBounds2D();
					if(bounds2d.getWidth()>0 && bounds2d.getHeight()>0){
	    	        	 processOCRShape(socr,s,bitmap,thin,(sf,lf)->{
	    	        		 got.add(Tuple.of(sf,lf));
	    	        	 });
	    	         }
	    	    	 return got.stream();
	    	     })
//	    	    .collect(Collectors.toList())
	    	    .forEach(t->{
	    	    	onFind.accept(t.k(), t.v());
	    	    });
    }
    
    private void processOCRShape(SCOCR socr, Shape s, Bitmap bitmap, Bitmap thin,BiConsumer<Shape,List<Tuple<Character,Number>>> onFind){
    	if(s.getBounds2D().getWidth()>0 && s.getBounds2D().getHeight()>0){
    		 double areareal=GeomUtil.area(s);
    		 
    		 if(areareal<=5)return;
    		 
	       	 LinkedHashSet<String> chars = new LinkedHashSet<String>();
	       	 
	       	 List<Tuple<Character,Number>> potential = socr.getNBestMatches(4,
	                    bitmap.crop(s),
	                    thin.crop(s))
	        			.stream()
	        			.map(Tuple::of)
	        			.map(t->adjustConfidence(t))
	        			.map(t->t.withVComparator())
	        			.sorted(Comparator.reverseOrder())
	        			.peek(t->{
	        				chars.add(t.k().toString());
	        			})
	        			.collect(Collectors.toList());
	        	
	       	 if(chars.contains("N") || chars.contains("S")|| chars.contains("s")){
	       		
	       		double areabox=GeomUtil.area(s.getBounds2D());
	       		 
		       		//this usually means it's not a real "N" or S
		       		
		       		if(areareal/areabox <0.5){
		       			if(chars.contains("\\") ||chars.contains("X")||chars.contains("K")||chars.contains("k")||chars.contains("-")){
				       		potential = potential.stream()
				       				             .filter(t->!t.k().toString().equals("N"))
				       				             .filter(t->!t.k().toString().equalsIgnoreCase("S"))
				       				             .collect(Collectors.toList());
		       			}
		       		}
		       		Rectangle2D rbox = s.getBounds2D();
		       		//probably not an N or S
	       		 	if(rbox.getWidth()>rbox.getHeight()*1.3){
		       		 	potential = potential.stream()
	  				             .filter(t->!t.k().toString().equals("N"))
	  				             .filter(t->!t.k().toString().equalsIgnoreCase("S"))
	  				             .collect(Collectors.toList());
	       		 	}
	       	 }
	       	 
	       	 if(chars.iterator().next().equals("L")){
	       		 	Rectangle2D rbox = s.getBounds2D();
	       		 	if(rbox.getWidth()>rbox.getHeight()*0.6){
	       		 		//Too wide for an L, but since it's the highest confidence, it's
	       		 		//probably just part of a bond system. 
	       		 		potential = potential.stream()
      				             .map(Tuple.vmap(n->(Number)0.0))
      				             .collect(Collectors.toList());
	       		 	}
		       	 }
	       	if(chars.contains("K") && chars.contains("X")){
	       		double areabox=GeomUtil.area(s.getBounds2D());
	       		if(areareal/areabox <0.5){
	       			
			       		potential = potential.stream()
			       							 .map(Tuple.vmap(n->(Number)0.0))
			       				             .collect(Collectors.toList());
	       			
	       		}
	       	 }
	       	
	       	if(chars.contains("S") && chars.contains("s")&& chars.contains("8")){
	       		//It's probably an S in this case, slightly adjust numbers for those
	       		potential = potential.stream()
				             .map(t->{
				            	 if(t.k().toString().equalsIgnoreCase("S")){
				            		 return Tuple.of(t.k(),(Number)Math.max(0, (1-(1-t.v().doubleValue())*0.8)));
				            	 }
				            	 return t;
				             })
				             .map(t->t.withVComparator())
				             .sorted(Comparator.reverseOrder())
				             .collect(Collectors.toList());
	       	 }
	       	
	       	if(chars.contains("D") && (chars.contains("U") || chars.contains("u"))){
	       		String best=chars.iterator().next();
	       		if(best.equals("D") || best.equalsIgnoreCase("U")){
	       		//It's probably an O, just got flagged wrong
	       		potential = potential.stream()
				             .map(Tuple.kmap(c->'O'))
				             .collect(Collectors.toList());
	       		}
	       	 }
        	
        	 onFind.accept(s, potential);
        }
    }
    
    
     
    public StructureImageExtractor load(File file) throws IOException{
    	ctabRaw.clear();
    	SCOCR[] socr=new SCOCR[]{OCR_DEFAULT.orElse(OCR_BACKUP, OCRcutoffCosine)};
    	
    	double[] maxBondLength=new double[]{INITIAL_MAX_BOND_LENGTH};    
        
        
        bitmap = Bitmap.read(file).clean();

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

        List<Shape> likelyOCR=new ArrayList<Shape>();
        List<Shape> likelyOCRAll=new ArrayList<Shape>();
        /*
         * Looks at each polygon, and gets the likely OCR chars.
         */   
        
        List<Shape> initialDebug = new ArrayList<>();
        
    	processOCR(socr[0],polygons,bitmap,thin,(s,potential)->{
    		 ocrAttmept.put(s, potential);
    		 initialDebug.add(s);
//    		 System.out.println("initial found:" +potential.get(0).k());
             if(potential.stream().filter(e->e.v().doubleValue()>OCRcutoffCosine).findAny().isPresent()){
            	 	if(OCRIsLikely(potential.get(0))){
                 		likelyOCR.add(s);
                 		
//                 		System.out.println("initial found:" +potential.get(0).k());
                 	}
                 	likelyOCRAll.add(s);
                 	
             }
    	});
        
//        if(likelyOCR.isEmpty()){
//        	socr[0]=OCR_BACKUP;
//        	likelyOCRAll.clear();
//        	processOCR(socr[0],polygons,bitmap,thin,(s,potential)->{
//       		 ocrAttmept.put(s, potential);
//                if(potential.stream().filter(e->e.v().doubleValue()>OCRcutoffCosine).findAny().isPresent()){
//               	 	if(OCRIsLikely(potential.get(0))){
//                    		likelyOCR.add(s);
//                    	}
//                    	likelyOCRAll.add(s);
//                }
//        	});
//        }
        
        double averageLargestOCR=likelyOCR.stream()
							              .map(s->GeomUtil.getPairOfFarthestPoints(s))
							              .mapToDouble(p->p[0].distance(p[1]))
							              .average()
							              .orElse(0);
        double averageAreaOCR=likelyOCR.stream()
	              .mapToDouble(s->GeomUtil.area(s))
	              .average()
	              .orElse(0);
        
        double averageWidthOCR=likelyOCR.stream()
	              .mapToDouble(s->s.getBounds2D().getWidth())
	              .average()
	              .orElse(0);
        //System.out.println("avg ocr:" + averageLargestOCR);
        
                		 
        likelyOCRAll.retainAll(likelyOCRAll.stream()
        			.filter(Objects::nonNull)
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
        
        OptionalDouble avgDistOCRToLine = Optional.of(0)
        		.filter(d->lDistOCRToLine.length>0)
        		.map(d->lDistOCRToLine[lDistOCRToLine.length/2])
        		.map(d->OptionalDouble.of(d))
        		.orElse(OptionalDouble.empty());
        		                    
        		                                  
        
        linesJoined=Stream.concat(bigLines.stream(),
        		            smallLines.stream())
        		    .collect(Collectors.toList());
        shapeTypes=new HashMap<>();
        
        
        
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
        
        List<Line2D> preprocess= GeomUtil.reduceMultiBonds(linesJoined, MAX_ANGLE_FOR_PARALLEL, MAX_DISTANCE_TO_MERGE_PARALLEL_LINES, MIN_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS,0,MAX_DELTA_LENGTH_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC)
        		                         .stream()
        		                         .map(t->t.k())
        		                         .collect(Collectors.toList());
        
       
        linesOrder=GeomUtil.reduceMultiBonds(preprocess, MAX_ANGLE_FOR_PARALLEL, largestBond/3, MIN_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS, MIN_BIGGER_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS,MAX_DELTA_LENGTH_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC);
        
              
        List<Shape> rescueOCRCandidates = new ArrayList<>();
        

        List<Shape> connectedComponents = polygons.stream()
        		                                  .map(s->GeomUtil.growShape(s, 2))
        		                                  .collect(Collectors.toList());
     
        int reps=0;
        boolean tooLongBond=true;
        ctabRaw.clear();
        while(tooLongBond){
        	rescueOCRCandidates.clear();
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
//	        		System.out.println("Wait just a gosh darn minute");
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
	        		double area=GeomUtil.area(candidate);
//	        		System.out.println("Area is:" + area);
	        		if(GeomUtil.area(candidate)>0.5*averageAreaOCR){
	        			candidate=GeomUtil.growShape(candidate,4);
	        			rescueOCRCandidates.add(candidate);
	        			//polygons.add(candidate);
//	        			System.out.println("Candidate");
	        			return GeomUtil.findCenterOfVertices(missingPoints);
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
//	        			System.out.println("Candidate");
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
	        
	        ctab.createNodesOnIntersectingLines(2, elist->{
	        	return true;
	        });
	        
//	        System.out.println("First intersection split:" + ctabRaw.size());
	        ctabRaw.add(ctab.cloneTab());
	        
	        
	        ctab.mergeNodesCloserThan(ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE);
	        ctab.standardCleanEdges();
	        
	        ctabRaw.add(ctab.cloneTab());
	        
	        double maxRatio=0.5;
	        double maxTotalRatio=1.4;
	        
	        if(avgDistOCRToLine.isPresent()){
	        	double nmaxRatio=(avgDistOCRToLine.getAsDouble())/ctab.getAverageBondLength();
	        	
	        	if(nmaxRatio>maxRatio){
		        	//TODO: sometimes there are only pairs of bonds, in such a case the heuristics don't
		        	//work well
	        		maxRatio=nmaxRatio;
		        	
		        	double maxlen=ctab.getEdges()
		        	    .stream()
		        	    .mapToDouble(e->e.getBondLength())
		        	    .max()
		        	    .orElse(1);
		        	
		        	maxlen=Math.max(maxlen, averageWidthOCR);
		        	
		        	maxTotalRatio = Math.max(maxTotalRatio, maxlen/ctab.getAverageBondLength());
	        	}
	        }
	        
	        ctab.mergeNodesExtendingTo(likelyOCR,maxRatio,maxTotalRatio);
//	        System.out.println("Extended the first time:" + ctabRaw.size());
	        ctabRaw.add(ctab.cloneTab());
	        
	        ctab.removeOrphanNodes();
	        ctab.mergeNodesCloserThan(ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO_FOR_MERGE);
	        
	        ctab.mergeNodesExtendingTo(likelyOCR,maxRatio,maxTotalRatio);
	        
	        //basic stuff
	        ctabRaw.add(ctab.cloneTab());
	        for(Shape s: likelyOCR){
	        	ctab.mergeAllNodesInsideCenter(s, OCR_TO_BOND_MAX_DISTANCE);
	        }
	        
//	        System.out.println("Adding back missing OCR:" + ctabRaw.size());
	        ctabRaw.add(ctab.cloneTab());
	        
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
		        		List<Edge> edges=cn.getEdges();
		        		if(edges.size()==2){
			        		if(ddelta<MAX_DELTA_LENGTH_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC){
			        			toRemove.add(cn);
			        			
			        			double o2=edges.stream().map(et->Tuple.of(et,et.getBondLength()))
			        					            .mapToDouble(e1->(e1.k().getOrder() * e1.v()))
			        					            .sum();
			        			int o=(int)Math.round(((o2/sumd)+0.2));
			        			t.v().setOrder(o);
			        		}
			        		if(!edges.stream().anyMatch(e2->e2.getDashed())){
			        			alreadyExists=true;	
			        		}
			        		
			        		
		        		}
		        		
	        		}
	        	}
	        	
	        	
	        	//System.out.println("Tol found for add:" + t.k());
	        	if(t.k()>MAX_TOLERANCE_FOR_SINGLE_BONDS){
	        		if(!alreadyExists){
	        			t.v().setDashed(true);
	        		}
	        	}
	        });
	        
	        
	        
//	        System.out.println("Made new bonds:" + ctabRaw.size());
	      //fuzzy adding missing stuff
	        ctabRaw.add(ctab.cloneTab());
	        toRemove.forEach(n->ctab.removeNodeAndEdges(n));
	        //ctab.removeOrphanNodes();
	        
	        
	        
	        
	        //fuzzy adding missing stuff
	        ctabRaw.add(ctab.cloneTab());
	     	        
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
        

        AtomicBoolean anyOtherIntersections = new AtomicBoolean(false);

        ctab.createNodesOnIntersectingLines(3, elist->{
        	
        	long longEnoughBonds=elist.stream()
        		//	.peek(e->System.out.println("BL:" + e.getBondLength()))
        	     .filter(e->e.getBondLength()>MAX_BOND_TO_AVG_BOND_RATIO_FOR_INTERSECTION*ctab.getAverageBondLength())
        	     .count();
        	if(longEnoughBonds<3)return false;
        	anyOtherIntersections.set(true);
        	return true;
        });
        
        if(anyOtherIntersections.get()){
//        	System.out.println("Second intersection split:" + ctabRaw.size());
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
        	                   //.peek(sc->System.out.println("Found candidate from earlier"))
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
	        	
	        	
	        	
        		nshape = GeomUtil.convexHull2(insideVertices2.toArray(new Point2D[0]));
        		
        		Point2D[] far=GeomUtil.getPairOfFarthestPoints(nshape);
        		
        		double arean = GeomUtil.area(nshape);
        		
        		double r=0;
        		if(far!=null){
        			r=far[0].distance(far[1]);	
        		}
        		if(r < averageLargestOCR*MIN_LONGEST_WIDTH_RATIO_FOR_OCR_TO_AVERAGE){
        			keep=false;
                }
        		if(arean < averageAreaOCR*MIN_AREA_RATIO_FOR_OCR_TO_AVERAGE){
        			keep=false;
                }
        		//polygons.add(nshape);
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
                   
                    if(nmap!=null && nthinmap!=null){
                    	processOCRShape(socr[0],nshape,bitmap,thin,(s,potential)->{
                    		
                    		if(potential.get(0).v().doubleValue()>OCRcutoffCosineRescue){
                    			
                    			String st=potential.get(0).k().toString();
                    			if(BranchNode.interpretOCRStringAsAtom(st)!=null){
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
                	likelyOCR.stream()
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
                	ocrAttmept.put(nshape, matches);
                	//System.out.println("rescue found:" +potential.get(0).k());
                	polygons.add(nshape);
					if (matches.get(0).v().doubleValue() > OCRcutoffCosineRescue) {
						if (OCRIsLikely(matches.get(0))) {
							likelyOCR.add(nshape);
						}
						likelyOCRAll.add(nshape);
					}
                });

        
        
        ctab.mergeNodesExtendingTo(likelyOCR,0.5,1.4);
        
        ctabRaw.add(ctab.cloneTab());
        
        double cosThetaOCRShape =Math.cos(MAX_THETA_FOR_OCR_SEPERATION);
        
        
        
        List<List<Shape>> ocrGroupList=GeomUtil.groupShapesIfClosestPointsMatchCriteria(likelyOCRAll, t->{
        	Point2D[] pts=t.v();
        	Shape[] shapes =t.k();
        	
        	List<Tuple<Character, Number>> attempt0 = ocrAttmept.get(shapes[0]);
        	List<Tuple<Character, Number>> attempt1 = ocrAttmept.get(shapes[1]);
			String v1= (attempt0 ==null || attempt0.isEmpty())? "" : attempt0.get(0).k().toString();
        	String v2= (attempt1 ==null || attempt1.isEmpty())? "" : attempt1.get(0).k().toString();
        	if(v1.equals("\\") || v1.equals("/") || 
        	   v2.equals("\\") || v2.equals("/")){
        		return false;
        	}
        	
        	Line2D l2 = new Line2D.Double(pts[0],pts[1]);
        	double dist=GeomUtil.length(l2);
        	double cutoff=Math.max(ctab.getAverageBondLength()*MAX_BOND_RATIO_FOR_OCR_CHAR_SPACING,averageWidthOCR);
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
//        	System.out.println("Angle was wrong:" + cosTheta + " vs " + cosThetaOCRShape);
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
	        		List<Tuple<Character, Number>> list = ocrAttmept.get(s);
					String v= (list ==null || list.isEmpty())? "":list.get(0).k().toString();
					
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

        ctab.standardCleanEdges();

        
        List<Shape> ocrMeaningful=bestGuessOCR.keySet()
				   .stream()
//				   .peek(t->System.out.println(bestGuessOCR.get(t)))
				   .filter(s->BranchNode.interpretOCRStringAsAtom(bestGuessOCR.get(s))!=null)
				   .collect(Collectors.toList());
        
//        System.out.println("Cleaned edges:" + ctabRaw.size());
        ctabRaw.add(ctab.cloneTab());
        //ctab.removeOrphanNodes();
        
        List<Node> alreadyFixedNodes = new ArrayList<Node>();
        
        
        bestGuessOCR.entrySet()
        		 .stream()
        		 .map(Tuple::of)
        		 .map(Tuple.vmap(s->Tuple.of(s,(s.equals("H"))?1:0).withVComparator()))
        		 .map(t->t.withVComparator())
        		 .sorted()
        		 .map(Tuple.vmap(t->t.k()))
        		 .forEach(shapeString->{
        			 Shape s= shapeString.k();
        			 String sym = shapeString.v();
        			 BranchNode actual = BranchNode.interpretOCRStringAsAtom(sym);
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
             		}
             		Point2D center = centert;
             		    
             		//This is likely the source of lots of problems
             		ctab.mergeAllNodesInside(s, MAX_BOND_RATIO_FOR_MERGING_TO_OCR*ctab.getAverageBondLength(),(n)->{
             			if(sym.equals("H")){
             				if(GeomUtil.findClosestShapeTo(ocrMeaningful, n.getPoint()).k() !=s){
             					return false;
             				}
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
        

//        System.out.println("Merged into atoms:" + ctabRaw.size());
        ctabRaw.add(ctab.cloneTab());

        List<Shape> growLikelyOCR=likelyOCR.stream().map(s->GeomUtil.growShape(s, 2)).collect(Collectors.toList());

        
        
        List<Line2D> lj =Stream.concat(removedTinyLines.stream(), linesJoined.stream())
                   .flatMap(l->GeomUtil.getLinesNotInside(l, growLikelyOCR).stream())
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
//        		System.out.println("Edge already exists? " +  already + " has bond:" + lst.v());
        		
        		
        		if(!already && haspossibleLine){
        			//sometimes adds wrong bonds
        			edgesToMake.add(Tuple.of(lst.v(),Tuple.of(nodes.get(0),nodes.get(1))));
        		}else if(!already && !haspossibleLine){
        			//do nothing
        		}else if(already && !haspossibleLine){
        			ctab.removeEdge(alreadyEdge);
        		}else{
        			//System.out.println("keep, due to line:" + possibleLine.getP1() + "," + possibleLine.getP2());
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
        
//        System.out.println("Made/removed bonds to OCR:" + ctabRaw.size());
        
        ctabRaw.add(ctab.cloneTab());
        
        ctab.makeMissingBondsToNeighbors(bitmap,MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL,MAX_TOLERANCE_FOR_SINGLE_BONDS,likelyOCR,OCR_TO_BOND_MAX_DISTANCE, (t)->{
        	
        	if(t.k()>MAX_TOLERANCE_FOR_SINGLE_BONDS){
        		t.v().setDashed(true);
        	}
        });
        
//        System.out.println("Made new bonds again:" + ctabRaw.size());
        ctabRaw.add(ctab.cloneTab());
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
        	    	List<KEqualityTuple<Node,Edge>> things1=neigh1.stream()
        	    												 .filter(ne->neigh2.contains(ne))
        	    												 .collect(Collectors.toList());
        	    	List<KEqualityTuple<Node,Edge>> things2=neigh2.stream()
																  .filter(ne->neigh1.contains(ne))
																  .collect(Collectors.toList());
        	    	List<KEqualityTuple<Node,Edge>> things = Stream.concat(things1.stream(), things2.stream())
        	    			                                       .collect(Collectors.toList());
        	    	
        	    	if(things.size()>0){
        	    		//System.out.println("Triangle found");
        	    		Point2D p1=n1.getPoint();
	        	    	Point2D p2=n2.getPoint();
	        	    	Point2D p3=things.get(0).k().getPoint();
	        	    	
	        	    	Edge oedge1=things.get(0).v();
	        	    	Edge oedge2=things.get(1).v();
	        	    	
	        	    	double tarea=Math.abs(GeomUtil.areaTriangle(p1,p2,p3));
	        	    	
	        	    	double expected = Math.sqrt(3)/4*Math.pow(e.getBondLength(),2);
	        	    	if(tarea<expected*0.5){
	        	    		//System.out.println("It's a bad one");
	        	    		if(!e.getDashed() && (oedge1.getDashed() && oedge2.getDashed()) ||
	        	    				(oedge1.getDashed() || oedge2.getDashed() && e.getOrder()>1)
	        	    				){
	        	    			ctab.removeEdge(oedge1);
	        	    			ctab.removeEdge(oedge2);
	        	    		}else{
	        	    			ctab.removeEdge(e);	
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
//        System.out.println("Removed Triangles:" + ctabRaw.size());
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
        
        List<Shape> appliedOCR = new ArrayList<Shape>();
        
        for(Shape s: bestGuessOCR.keySet()){
        	String sym=bestGuessOCR.get(s);
//        	System.out.println(sym);
        	BranchNode actual=BranchNode.interpretOCRStringAsAtom(sym);
        	if(actual!=null && actual.isRealNode()){
        		appliedOCR.add(s);
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
//        System.out.println("Added Branch nodes:" + ctabRaw.size());
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
            .filter(t->t.v().get(0).v().getBondLength()<ctab.getAverageBondLength())
            .filter(t->t.v().get(1).v().getBondLength()<ctab.getAverageBondLength())
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
//	                	System.out.println("Adding edge:" + n2.getIndex() + "->" + n3.getIndex());
//	            		System.out.println("Then remove:" + n1.getIndex());
	            		ctab.addEdge(n2.getIndex(), n3.getIndex(), 1);
	            		toRemove.add(n1);
                	}
            	}
            });        
        toRemove.forEach(n->ctab.removeNodeAndEdges(n));
        ctab.standardCleanEdges();
        }while(!toRemove.isEmpty());

        
        
        //Cleanup "duplicate" lines that are probably problems. 
        //criteria is:
        //1. 
        //2. 
        
        ctab.getNodes()
            .stream()
            .filter(n->n.getEdgeCount()>=2)
            .map(n->{
            	return Tuple.of(n,GeomUtil.eachCombination(n.getEdges())
            	        .map(t->{
            	        	if(t.k().getBondLength()>t.v().getBondLength()){
            	        		t=t.swap();
            	        	}
            	        	return t;
            	        })
            	        .filter(t->t.v().getBondLength()>=ctab.getAverageBondLength())
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
	     	        		if(sd1<t.v().getBondLength()){
	     	        		//remove long bond
	     	        		//change point
	     	        		//add edge to long bond other node
	     	        		
//	     	        		System.out.println("Gonna remove edge:" + t.v().toString());
	     	        		tnode.setPoint(ppnt);
	     	        		ctab.addEdge(tnode.getIndex(), t.v().getOtherNode(n).getIndex(), t.v().getOrder());
	     	        		ctab.removeEdge(t.v());
	     	        		}
	     	        	}
            	    });
     	        
            });;
        
        ctab.standardCleanEdges();
        
        
//        System.out.println("Finished:" + ctabRaw.size());
        ctabRaw.add(ctab.cloneTab());
        rescueOCRShapes=rescueOCRCandidates;
        
        

        

        
        ctab.getEdges()
        .forEach(e->{
        	Line2D useLine=GeomUtil.getLinesNotInside(e.getLine(), growLikelyOCR)
        	        .stream()
        	        .map(l->Tuple.of(l, GeomUtil.length(l)).withVComparator())
        	        .max(Comparator.naturalOrder())
        	        .map(t->t.k())
        	        .orElse(null);
        	if(useLine!=null){
        		int mult=1;
        		if(e.getPoint1().distance(useLine.getP1())< e.getPoint2().distance(useLine.getP1())){
        			mult=-1;
        		}
	        	double wl=mult*bitmap.getWedgeLikeScore(useLine);
	        	if(wl>WEDGE_LIKE_PEARSON_SCORE_CUTOFF){
	        		e.setWedge(true);	
	        		e.switchNodes();
	        	}else if(wl<-WEDGE_LIKE_PEARSON_SCORE_CUTOFF){
	        		e.setWedge(true);
	        	}           
        	}
        });
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



	public List<ConnectionTable> getCtabRaw() {
		return this.ctabRaw;
	}

	public Map<Shape,ShapeInfo> getShapeTypes() {
		return shapeTypes;
	}

	
	
	
}
