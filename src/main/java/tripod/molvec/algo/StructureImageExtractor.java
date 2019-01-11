package tripod.molvec.algo;

import java.awt.Shape;
import java.awt.geom.Line2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import gov.nih.ncats.chemkit.api.Chemical;
import tripod.molvec.Bitmap;
import tripod.molvec.CachedSupplier;
import tripod.molvec.algo.LineUtil.ConnectionTable;
import tripod.molvec.ui.RasterCosineSCOCR;
import tripod.molvec.ui.SCOCR;

public class StructureImageExtractor {
	static final SCOCR OCR=new RasterCosineSCOCR();
   
	static{
    	OCR.setAlphabet(SCOCR.SET_COMMON_CHEM_ALL());
    }
    
	private Bitmap bitmap; // original bitmap
    private Bitmap thin; // thinned bitmap
    
    
    private List<Shape> polygons;
    private List<Line2D> lines;
    private List<Tuple<Line2D,Integer>> linesOrder;    
    private Map<Shape,List<Entry<Character,Number>>> ocrAttmept = new HashMap<Shape,List<Entry<Character,Number>>>();
    
    private ConnectionTable ctab;
    
    
    private final double MAX_REPS = 10;
    private final double INITIAL_MAX_BOND_LENGTH=Double.MAX_VALUE;
    private final double MIN_BOND_TO_AVG_BOND_RATIO = 1/2.5;
    private final double MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL = 1.3;
    private final double MAX_TOLERANCE_FOR_DASH_BONDS = 1.0;
	private final double cutoffCosine=0.6;
	private final double MAX_BOND_TO_AVG_BOND_RATIO_TO_KEEP = 1.5;
	private final double MIN_DISTANCE_BEFORE_MERGING_NODES = 4.0;
	private final double maxRatioForIntersection = 1.2;
	private final double OCR_TO_BOND_MAX_DISTANCE=2.0;
	private final double maxCandidateRatioForIntersection = 1.7;        
	private final double MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS = 2;
	private final double MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS = 3;
    
	private final HashSet<String> accept = CachedSupplier.of(()->{
		HashSet<String> keepers=new HashSet<String>();
		keepers.add("C");
		keepers.add("N");
		keepers.add("O");
		keepers.add("H");
		keepers.add("S");
		keepers.add("P");
		keepers.add("B");
		keepers.add("F");
		
	    return keepers;
	}).get();
    
	
    public StructureImageExtractor(){
    	
    }
     
    public StructureImageExtractor load(File file) throws IOException{
    	
    	double[] maxBondLength=new double[]{INITIAL_MAX_BOND_LENGTH};    
        
        
        bitmap = Bitmap.read(file);

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
        
        lines= LineUtil.asLines(thin.segments());
        
        List<Line2D> smallLines=lines.stream()
        						     .filter(l->LineUtil.length(l)<4)
        						     .collect(Collectors.toList());
        
        List<Line2D> bigLines=lines.stream()
        		                   .filter(l->LineUtil.length(l)>=4)
        		                   .collect(Collectors.toList());
        
        smallLines= thin.combineLines(smallLines, MAX_DISTANCE_FOR_STITCHING_SMALL_SEGMENTS, MAX_TOLERANCE_FOR_STITCHING_SMALL_SEGMENTS);
        lines=Stream.concat(bigLines.stream(),
        		            smallLines.stream())
        		    .collect(Collectors.toList());
        
        
        double largestBond=lines.stream()
		           .mapToDouble(l->l.getP1().distance(l.getP2()))
		           .max()
		           .getAsDouble();
        
        List<Line2D> preprocess= LineUtil.reduceMultiBonds(lines, 5 * Math.PI/180.0, 2, .5)
        		                         .stream()
        		                         .map(t->t.k())
        		                         .collect(Collectors.toList());
        
        
        linesOrder=LineUtil.reduceMultiBonds(preprocess, 5 * Math.PI/180.0, largestBond/3, .5);
        
              
        
        
        List<Shape> likelyOCR=new ArrayList<Shape>();
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
                 if(potential.stream().filter(e->e.getValue().doubleValue()>cutoffCosine).findAny().isPresent()){
                	 String t=potential.get(0).getKey()+"";
                 	if(!"I".equalsIgnoreCase(t) && 
                 	   !"L".equalsIgnoreCase(t) &&
                 	   !"1".equalsIgnoreCase(t)){
                 		likelyOCR.add(s);
                 	}
                 }
             }
         }
        
     
        int reps=0;
        boolean tooLongBond=true;
        while(tooLongBond){
	        ctab = LineUtil.getConnectionTable(linesOrder, likelyOCR, maxRatioForIntersection, maxCandidateRatioForIntersection,l-> (LineUtil.length(l) < maxBondLength[0]))
	        		       .mergeNodesCloserThan(MIN_DISTANCE_BEFORE_MERGING_NODES);
	        
	        for(Shape s: likelyOCR){
	        	ctab.mergeAllNodesInside(s, OCR_TO_BOND_MAX_DISTANCE);
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
	        ctab.makeDashBondsToNeighbors(bitmap,MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL,MAX_TOLERANCE_FOR_DASH_BONDS);
	        
	        
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
        	e.order=1;
        });
        
        
        
        
        
        for(Shape s: likelyOCR){
        	String sym=(ocrAttmept.get(s).get(0).getKey() + "").toUpperCase();
        	if(accept.contains(sym)){
        		ctab.setNodeToSymbol(s, sym);
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

	public List<Line2D> getLines() {
		return lines;
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
