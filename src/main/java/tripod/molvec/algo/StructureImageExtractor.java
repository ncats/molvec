package tripod.molvec.algo;

import java.awt.Shape;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import gov.nih.ncats.chemkit.api.Chemical;
import tripod.molvec.Bitmap;
import tripod.molvec.CachedSupplier;
import tripod.molvec.ui.RasterCosineSCOCR;
import tripod.molvec.ui.SCOCR;
import tripod.molvec.util.ConnectionTable;
import tripod.molvec.util.GeomUtil;

public class StructureImageExtractor {
	static final SCOCR OCR=new RasterCosineSCOCR();
   
	static{
		Set<Character> alpha=SCOCR.SET_COMMON_CHEM_ALL();
		alpha.add(Character.valueOf('/'));
		alpha.add(Character.valueOf('\\'));
    	OCR.setAlphabet(alpha);
    }
    
	private Bitmap bitmap; // original bitmap
    private Bitmap thin; // thinned bitmap
    
    
    private List<Shape> polygons;
    private List<Line2D> lines;
    private List<Line2D> linesJoined;
    private List<Tuple<Line2D,Integer>> linesOrder;    
    private Map<Shape,List<Entry<Character,Number>>> ocrAttmept = new HashMap<Shape,List<Entry<Character,Number>>>();
    private Map<Shape,String> bestGuessOCR = new HashMap<>();
    
    private ConnectionTable ctab;
    
    
    private final double MAX_REPS = 10;
    private final double INITIAL_MAX_BOND_LENGTH=Double.MAX_VALUE;
    private final double MIN_BOND_TO_AVG_BOND_RATIO = 1/2.5;
    private final double MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL = 1.3;
    
    private final double MAX_TOLERANCE_FOR_DASH_BONDS = 2.0;
    private final double MAX_TOLERANCE_FOR_SINGLE_BONDS = 0.5;
    
	private final double OCRcutoffCosine=0.65;
	private final double MAX_BOND_TO_AVG_BOND_RATIO_TO_KEEP = 1.5;
	private final double MIN_DISTANCE_BEFORE_MERGING_NODES = 4.0;
	private final double maxRatioForIntersection = 1.2;
	private final double maxPerLineDistanceRatioForIntersection = 2;
	private final double OCR_TO_BOND_MAX_DISTANCE=2.0;
	private final double maxCandidateRatioForIntersection = 1.7;        
	private final double MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS_THIN = 1;
	private final double MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS_FULL = 0.5;
	private final double MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS = 6;
	
	private final double MAX_ANGLE_FOR_JOINING_SEGMENTS=25 * Math.PI/180.0;
	private final double MIN_SIZE_FOR_ANGLE_COMPARE_JOINING_SEGMENTS=8.0;
	
	private final double MAX_DISTANCE_TO_MERGE_PARALLEL_LINES=2;
	private final double MAX_POINT_DISTANCE_TO_BE_PART_OF_MULTI_NODE= 1;
	
	private final double MAX_BOND_RATIO_FOR_OCR_CHAR_SPACING=0.3;
	private final double MAX_THETA_FOR_OCR_SEPERATION=45 * Math.PI/180.0;
	private final double MAX_BOND_RATIO_FOR_MERGING_TO_OCR=0.5;
	
	//For finding high order bonds
	private final double MIN_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS=.5;
	private final double MIN_BIGGER_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS=.25;
	private final double MAX_ANGLE_FOR_PARALLEL=10.0 * Math.PI/180.0;
	
	//Parallel lines
	
	
	
    
	private final HashSet<String> accept = CachedSupplier.of(()->{
		HashSet<String> keepers=new HashSet<String>();
		keepers.add("C");
		keepers.add("N");
		keepers.add("O");
		keepers.add("H");
		keepers.add("S");
		keepers.add("P");
		keepers.add("B");
		keepers.add("Br");
		keepers.add("Cl");
		keepers.add("F");
		
	    return keepers;
	}).get();
	
	
	private final HashSet<String> reject = CachedSupplier.of(()->{
		HashSet<String> keepers=new HashSet<String>();
		
		keepers.add("I");
		keepers.add("i");
		keepers.add("l");
		keepers.add("-");
		keepers.add("/");
		keepers.add("\\");
		
	    return keepers;
	}).get();
	
	
	private String interpretOCRStringAsAtom(String s){
		if(accept.contains(s)){
			return s;
		}else if(accept.contains(s.toUpperCase())){
			return s.toUpperCase();
		}
		if(s.contains("H")){
			return interpretOCRStringAsAtom(s.replace("H", ""));
		}
		if(s.contains("I")){
			return interpretOCRStringAsAtom(s.replace("I", "l"));
		}
		if(s.contains("c")){
			return interpretOCRStringAsAtom(s.replace("c", "C"));
		}
		
		return null;
			
	}
	
    
	
    public StructureImageExtractor(){
    	
    }
     
    public StructureImageExtractor load(File file) throws IOException{
    	
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
        for (Shape s : polygons) {
             if(s.getBounds2D().getWidth()>0 && s.getBounds2D().getHeight()>0){
             	List<Entry<Character,Number>> potential = OCR.getNBestMatches(4,
                         bitmap.crop(s).createRaster(),
                         thin.crop(s).createRaster()
                         );
                 ocrAttmept.put(s, potential);
                 if(potential.stream().filter(e->e.getValue().doubleValue()>OCRcutoffCosine).findAny().isPresent()){
                	 String t=potential.get(0).getKey()+"";
                 	if(!"I".equalsIgnoreCase(t) && 
                 	   !"L".equalsIgnoreCase(t) &&
                 	   !"1".equalsIgnoreCase(t) &&
                 	   !"-".equalsIgnoreCase(t) &&
                 	   !"/".equalsIgnoreCase(t) &&
                 	   !"\\".equalsIgnoreCase(t)){
                 		likelyOCR.add(s);
                 	}
                 	likelyOCRAll.add(s);
                 	
                 }
             }
        }
        
        
        
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
        	if(true)return true;
        	return false;
        	//return LineUtil.length(l)<largestBond;
        });
        
        
        
        List<Line2D> smallLines=lines.stream()
        						     .filter(tryToMerge)
        						     .collect(Collectors.toList());
        
        List<Line2D> bigLines=lines.stream()
        		                   .filter(tryToMerge.negate())
        		                   .collect(Collectors.toList());
        
        smallLines= thin.combineLines(smallLines, MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS, MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS_THIN, MAX_POINT_DISTANCE_TO_BE_PART_OF_MULTI_NODE,MAX_ANGLE_FOR_JOINING_SEGMENTS,MIN_SIZE_FOR_ANGLE_COMPARE_JOINING_SEGMENTS);
        
        smallLines= bitmap.combineLines(smallLines, MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS, MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS_FULL, MAX_POINT_DISTANCE_TO_BE_PART_OF_MULTI_NODE,MAX_ANGLE_FOR_JOINING_SEGMENTS,MIN_SIZE_FOR_ANGLE_COMPARE_JOINING_SEGMENTS);
        
        linesJoined=Stream.concat(bigLines.stream(),
        		            smallLines.stream())
        		    .collect(Collectors.toList());
        
        double largestBond=linesJoined.stream()
		           .mapToDouble(l->GeomUtil.length(l))
		           .max()
		           .orElse(0);
       
        System.out.println("Angle:");
        System.out.println(MAX_ANGLE_FOR_PARALLEL);
        System.out.println(Math.cos(MAX_ANGLE_FOR_PARALLEL));
        
        List<Line2D> preprocess= GeomUtil.reduceMultiBonds(linesJoined, MAX_ANGLE_FOR_PARALLEL, MAX_DISTANCE_TO_MERGE_PARALLEL_LINES, MIN_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS,0)
        		                         .stream()
        		                         .map(t->t.k())
        		                         .collect(Collectors.toList());
        
        System.out.println("Detecting multi bonds");
        linesOrder=GeomUtil.reduceMultiBonds(preprocess, MAX_ANGLE_FOR_PARALLEL, largestBond/3, MIN_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS, MIN_BIGGER_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS);
        
              
        
        
        
     
        int reps=0;
        boolean tooLongBond=true;
        
        while(tooLongBond){
        	List<Tuple<Line2D,Integer>> linesOrderRestricted =linesOrder.stream()
        	          .filter(t->{
        	        	  Line2D l=t.k();
        	        	  return isInOCRShape.negate().test(l);
        	          })
        	          .collect(Collectors.toList());
        	
	        ctab = GeomUtil.getConnectionTable(linesOrderRestricted, likelyOCR, maxRatioForIntersection, maxCandidateRatioForIntersection,maxPerLineDistanceRatioForIntersection,l-> (GeomUtil.length(l) < maxBondLength[0]))
	        		       .mergeNodesCloserThan(MIN_DISTANCE_BEFORE_MERGING_NODES);
	        
	        for(Shape s: likelyOCR){
	        	ctab.mergeAllNodesInsideCenter(s, OCR_TO_BOND_MAX_DISTANCE);
	        }
	        double bl1=ctab.getAverageBondLength();
	        ctab.mergeNodesCloserThan(bl1*MIN_BOND_TO_AVG_BOND_RATIO);
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
	        ctab.mergeNodesExtendingTo(likelyOCR);
	        
	        ctab.removeOrphanNodes();
	        ctab.mergeNodesCloserThan(ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO);
	        
	        ctab.makeMissingBondsToNeighbors(bitmap,MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL,MAX_TOLERANCE_FOR_DASH_BONDS,likelyOCR,OCR_TO_BOND_MAX_DISTANCE, (t)->{
//	        	System.out.println("Tol found:" + t.k());
	        	if(t.k()>MAX_TOLERANCE_FOR_SINGLE_BONDS){
	        		t.v().setDashed(true);
	        	}
	        });
	        
	     	        
	        double avgBondLength=ctab.getAverageBondLength();
	        maxBondLength[0]=avgBondLength*MAX_BOND_TO_AVG_BOND_RATIO_TO_KEEP;
	        
	        
	        
	        //System.out.println("Average bond length:" + avgBondLength);
	        
	        tooLongBond = ctab.getEdges()
	        		          .stream()
	        	//	          .peek(e->System.out.println(e.getBondDistance()))
	        		          .filter(e->e.getBondDistance()>maxBondLength[0])
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
        	System.out.println("Fixing bond order");
        	e.setOrder(1);
        });
        
        
        
        
        
        double cosThetaOCRShape =Math.cos(MAX_THETA_FOR_OCR_SEPERATION);
        
        List<List<Shape>> ocrGroupList=GeomUtil.groupShapesIfClosestPointsMatchCriteria(likelyOCRAll, pts->{
        	Line2D l2 = new Line2D.Double(pts[0],pts[1]);
        	double dist=GeomUtil.length(l2);
        	if(dist>ctab.getAverageBondLength()*MAX_BOND_RATIO_FOR_OCR_CHAR_SPACING){
        		return false;
        	}
        	double[] vec=GeomUtil.asVector(l2);
        	double cosTheta=Math.abs(vec[0]/dist);
        	
        	if(cosTheta>cosThetaOCRShape){
        		return true;
        	}
        	return false;
        });
        

        bestGuessOCR= ocrGroupList.stream()
        	//.filter(l->l.size()>1)
	        .map(g->{
	        	String st=g.stream()
	        			.map(s->Tuple.of(s,s))
	        			.map(Tuple.vmap(s->s.getBounds2D().getMinX()))
	        			.map(t->t.withVComparator())
	        			.sorted()
	        			.map(t->t.k())
	        			.map(s->(ocrAttmept.get(s).get(0).getKey() + ""))
	        			.collect(Collectors.joining());
	        	Shape s = g.get(0);
	        	
	        	Shape bigShape=g.stream()
	        					.reduce((s1,s2)->GeomUtil.add(s1, s2))
	        					.orElse(s);
	        	System.out.println(st);
	        	bestGuessOCR.put(bigShape, st);
	        	return Tuple.of(bigShape,st);
	        })
	        .collect(Tuple.toMap());
        
        List<Shape> ocrMeaningful=bestGuessOCR.keySet()
				   .stream()
				   .filter(s->interpretOCRStringAsAtom(bestGuessOCR.get(s))!=null)
				   .collect(Collectors.toList());
        
        for(Shape s: bestGuessOCR.keySet()){
        	String sym=bestGuessOCR.get(s);
        	String actual=this.interpretOCRStringAsAtom(sym);
        	if(actual!=null){
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
        
        
        for(Shape s: bestGuessOCR.keySet()){
        	String sym=bestGuessOCR.get(s);
        	System.out.println("Trying to add:" + sym);
        	String actual=this.interpretOCRStringAsAtom(sym);
        	if(actual!=null){
        		ctab.setNodeToSymbol(s, actual);
        	}
        	
        }
        

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

	public Map<Shape, List<Entry<Character, Number>>> getOcrAttmept() {
		return ocrAttmept;
	}

	public ConnectionTable getCtab() {
		return ctab;
	}

	
	
}
