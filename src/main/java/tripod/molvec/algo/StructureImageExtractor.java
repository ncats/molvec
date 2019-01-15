package tripod.molvec.algo;

import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import gov.nih.ncats.chemkit.api.Chemical;
import tripod.molvec.Bitmap;
import tripod.molvec.CachedSupplier;
import tripod.molvec.ui.RasterCosineSCOCR;
import tripod.molvec.ui.SCOCR;
import tripod.molvec.util.ConnectionTable;
import tripod.molvec.util.ConnectionTable.Edge;
import tripod.molvec.util.ConnectionTable.Node;
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
    private Map<Shape,List<Tuple<Character,Number>>> ocrAttmept = new HashMap<>();
    private Map<Shape,String> bestGuessOCR = new HashMap<>();
    
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
	private final double maxCandidateRatioForIntersection = 1.7;        
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
	private final double MIN_LONGEST_WIDTH_RATIO_FOR_OCR_TO_AVERAGE=0.5;
	
	//For finding high order bonds
	private final double MIN_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS=.5;
	private final double MIN_BIGGER_PROJECTION_RATIO_FOR_HIGH_ORDER_BONDS=.25;
	private final double MAX_ANGLE_FOR_PARALLEL=10.0 * Math.PI/180.0;
	
	private final double MIN_ST_DEV_FOR_KEEPING_DASHED_LINES=0.07;
	
	
	
	//Parallel lines
	
	private final double MAX_DELTA_LENGTH_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC=5.0;
	private final double MAX_ANGLE_FOR_PARALLEL_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC=15.0 * Math.PI/180.0;
	
    
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
	
	private BranchNode interpretOCRStringAsAtom(String s, boolean tokenOnly){
		if(accept.contains(s)){
			return new BranchNode(s);
		}else if(accept.contains(s.toUpperCase())){
			return new BranchNode(s.toUpperCase());
		}
		if(s.contains("H")){
			return interpretOCRStringAsAtom(s.replace("H", ""),tokenOnly);
		}
		if(s.contains("I")){
			return interpretOCRStringAsAtom(s.replace("I", "l"),tokenOnly);
		}
		if(s.contains("c")){
			return interpretOCRStringAsAtom(s.replace("c", "C"),tokenOnly);
		}
		if(s.equals("0")){
			return interpretOCRStringAsAtom(s.replace("0", "O"),tokenOnly);
		}
		
		try{
			int r=Integer.parseInt(s);
			if(r>0){
				BranchNode repNode=new BranchNode("?");
				repNode.setRepeat(r);
				return repNode;
			}
		}catch(Exception e){
			
		}
		if(tokenOnly){
			return null;
		}
		
		//start breaking up
		if(s.length()>1){
			BranchNode parent=null;
			String rem=null;
			for(int i=s.length()-1;i>=1;i--){
				String substr=s.substring(0, i);
				BranchNode bn1=interpretOCRStringAsAtom(substr,true);
				if(bn1!=null){
					parent=bn1;
					rem=s.substring(i);
					break;
				}
			}
			if(parent!=null){
				BranchNode child=interpretOCRStringAsAtom(rem);
				if(child!=null){
					BranchNode bnFinal= parent.addChild(child);
					return bnFinal;
				}
			}
		}
		
		return null;		
	}
	private BranchNode interpretOCRStringAsAtom(String s){
		return interpretOCRStringAsAtom(s,false);
	}
	
	private static class BranchNode{
		String symbol = "C";
		boolean pseudo=false;
		boolean isRepeat=false;
		int rep=0;
		
		List<BranchNode> children = new ArrayList<BranchNode>();
		
		Point2D suggestedPoint = new Point2D.Double(0, 0);
		
		
		public BranchNode generateCoordinates(){
			suggestedPoint = new Point2D.Double(0, 0);
			int totalChildren = children.size();
			
			double dtheta = (2*Math.PI)/(totalChildren+1.0);
			
			for(int i=0;i<totalChildren;i++){
				BranchNode child = children.get(i);
				child.generateCoordinates();
				AffineTransform at = new AffineTransform();
				
				double ntheta = dtheta*(i+1);
				at.rotate(-ntheta);
				at.translate(-1, 0);
				child.applyTransform(at);
			}
			return this;
			
		}
		
		public BranchNode applyTransform(AffineTransform at){
			
			suggestedPoint=at.transform(suggestedPoint, null);
			for(BranchNode c:this.getChildren()){
				c.applyTransform(at);
			}
			return this;
		}
		
		public String getSymbol(){
			return this.symbol;
		}
		
		public boolean hasChildren(){
			return !children.isEmpty();
		}
		public List<BranchNode> getChildren(){
			return this.children;
		}
		public boolean isRealNode(){
			return !pseudo && !isRepeat;
		}
		public boolean isPseudoNode(){
			return pseudo;
		}
		public boolean isRepeatNode(){
			return this.isRepeat;
		}
		public int getRepeat(){
			return this.rep;
		}
		public BranchNode setRepeat(int r){
			this.isRepeat=true;
			this.rep=r;
			return this;
		}
		
		public BranchNode addChild(BranchNode bn){
			if(bn.isRepeatNode()){
				int rep=bn.getRepeat();
				this.setPseudoNode(true);
				for(int i=0;i<rep;i++){
					this.addChild(new BranchNode(this.symbol));
				}
				if(bn.hasChildren()){
					if(bn.children.size()>1){
						throw new IllegalStateException("Repeat branch nodes should not have more than one child");
					}
					BranchNode realNode=bn.children.get(0);
					return realNode.addChild(this);					
				}
			}else if(bn.isPseudoNode()){
				this.children.addAll(bn.children);
			}else{
				this.children.add(bn);
			}
			return this;
		}
		
		public BranchNode setPseudoNode(boolean p){
			this.pseudo=p;
			return this;
		}
		
		
		public BranchNode(String s){
			this.symbol=s;
		}
		
		public String toString(){
			return this.symbol + "[" + this.children.stream().map(b->b.toString()).collect(Collectors.joining(",")) + "]";
		}
		
		private void forEachBranchNode(BiConsumer<BranchNode,BranchNode> cons, BranchNode parent){
			cons.accept(parent, this);
			for(BranchNode child:this.getChildren()){
				child.forEachBranchNode(cons,this);
			}
		}
		
		public void forEachBranchNode(BiConsumer<BranchNode,BranchNode> parentAndChild){
			forEachBranchNode(parentAndChild,null);
		}
		
		
	}
	
    
	
    public StructureImageExtractor(){
    	
    }
    
    private static Tuple<Character,Number> adjustConfidence(Tuple<Character,Number> tup){
    	String ch=tup.k()+"";
    	double invScore=1-tup.v().doubleValue();
    	if(ch.equals("K") || ch.equals("k") || ch.equals("f")){
    		invScore=invScore*3.5; // penalize "K"
    	}
    	if(ch.equalsIgnoreCase("R")||
    	   ch.equalsIgnoreCase("A")||
    	   ch.equalsIgnoreCase("Z")||
    	   ch.equalsIgnoreCase("-")||
    	   ch.equals("n")){
    		invScore=invScore*3; // penalize
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
             	List<Tuple<Character,Number>> potential = OCR.getNBestMatches(4,
                         bitmap.crop(s).createRaster(),
                         thin.crop(s).createRaster()
                         )
             			.stream()
             			.map(Tuple::of)
             			.map(t->adjustConfidence(t))
             			.collect(Collectors.toList());
                 ocrAttmept.put(s, potential);
                 if(potential.stream().filter(e->e.v().doubleValue()>OCRcutoffCosine).findAny().isPresent()){
                	 
                	 	if(OCRIsLikely(potential.get(0))){
	                 		likelyOCR.add(s);
	                 	}
	                 	likelyOCRAll.add(s);
                 	
                 }
             }
        }
        
        double averageLargestOCR=likelyOCR.stream()
							              .map(s->GeomUtil.getPairOfFarthestPoints(s))
							              .mapToDouble(p->p[0].distance(p[1]))
							              .average()
							              .orElse(0);
        System.out.println("avg ocr:" + averageLargestOCR);
        
                		 
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
	        
	        
	        ctab.getTolerancesForAllEdges(bitmap)
	    	    .forEach(t->{
	    	    	//System.out.println("Score for edge:" + t.v());
	    	    	if(t.v()>MAX_TOLERANCE_FOR_DASH_BONDS && t.k().getOrder()==1){	    	    		
	    	    		ctab.removeEdge(t.k());
	    	    	}
	    	    });
	        
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
	        ctab.mergeNodesExtendingTo(likelyOCR);
	        
	        ctab.removeOrphanNodes();
	        ctab.mergeNodesCloserThan(ctab.getAverageBondLength()*MIN_BOND_TO_AVG_BOND_RATIO);
	        
	        
	        ctab.makeMissingNodesForShapes(likelyOCR,MAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL,MIN_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL);
	        
	        List<Node> toRemove = new ArrayList<Node>();
	        
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
		        		double ddelta=Math.abs(sumd-e.getBondDistance());
		        		if(ddelta<MAX_DELTA_LENGTH_FOR_STITCHING_LINES_ON_BOND_ORDER_CALC){
		        			System.out.println("Should remove this one");
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
        
        
        double avgTol=ctab.getTolerancesForAllEdges(bitmap)
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
        
        
        
        
        System.out.println("Nodes not in shapes:"+ unmatchedNodes.size());
        
        List<Point2D> vertices=lines.stream().flatMap(l->Stream.of(l.getP1(),l.getP2())).collect(Collectors.toList());
        List<Point2D> verticesJ=linesJoined.stream().flatMap(l->Stream.of(l.getP1(),l.getP2())).collect(Collectors.toList());
        
        
        List<Shape> toAddAllOCR=new ArrayList<Shape>();
        
        
        unmatchedNodes.forEach(n->{
        	Tuple<Edge,Double> worstEdge=ctab.getWorstToleranceForNode(n,bitmap);
        	
//        	System.out.println("Avg tol:" + avgTol);
//        	System.out.println("Worst tol here:" + worstEdge.v());
        	
        	
        	//The worst bond isn't so bad, maybe it's not worth looking
        	if(worstEdge.v()<avgTol*PROBLEM_TOLERANCE_RATIO){
        		
        		//double check to see if the shortest bond is too short
        		//if it isn't, then skip this one, it's probably fine.
        		boolean tooSmallBond=n.getEdges()
        				              .stream()
						              .mapToDouble(e->e.getBondDistance())
						              .filter(d->d<avgbond*PROBLEM_BOND_LENGTH_RATIO)
						              .findAny()
						              .isPresent();
        		if(tooSmallBond){
        			//continue
        		}else{
        			System.out.println("Skip");
        			return;
        		}
        	}
        	System.out.println("Keep");
        	
        	
        	
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
        	polygons.add(nshape);
        	if(keep){
        		Bitmap nmap=bitmap.crop(nshape);
                Bitmap nthinmap=thin.crop(nshape);
                if(nmap!=null && nthinmap!=null){
                	System.out.println("And it looks promising");
                	nshape=GeomUtil.growShape(nshape, 2);
                	nmap=bitmap.crop(nshape);
                    nthinmap=thin.crop(nshape);
                    
                    List<Shape> slist=nmap.connectedComponents(Bitmap.Bbox.Polygon);
                    
                    System.out.println("Shapes are:" + slist);
                    Shape bshape=slist.stream()
			                            .map(s->Tuple.of(s,-s.getBounds2D().getWidth()*s.getBounds2D().getHeight()).withVComparator())
			                            .sorted()
			                            .map(t->t.k())
			                            .findFirst()
			                            .orElse(nshape);
                    Rectangle2D rect1 = nshape.getBounds2D();
                    AffineTransform at = new AffineTransform();
                    at.translate(rect1.getX(),rect1.getY());
                    nshape=at.createTransformedShape(bshape).getBounds2D();
//                    
                    
                    nmap=bitmap.crop(nshape);
                    nthinmap=thin.crop(nshape);
                   
                    if(nmap!=null && nthinmap!=null){
	                    
	                	List<Tuple<Character,Number>> potential = OCR.getNBestMatches(4,
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
		                	System.out.println("Found another pot:" +potential.get(0).k());
		                	System.out.println("Found another pot score:" +potential.get(0).v());
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
                	//if(ctab.getNodesInsideShape(nshape, 0).isEmpty())return;
                	Bitmap nmap=bitmap.crop(nshape);
                	Bitmap nthinmap=thin.crop(nshape);
                	if(nmap!=null && nthinmap!=null){
	                    
	                	List<Tuple<Character,Number>> potential = OCR.getNBestMatches(4,
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
		                	System.out.println("Found another:" +potential.get(0).k());
		                	System.out.println("Found another:" +potential.get(0).v());
//		                	
			                if(OCRIsLikely(potential.get(0))){
			                	likelyOCR.add(nshape);
			                }
			                likelyOCRAll.add(nshape);
		                }
                    }
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
	        	
//	        	String st=sorted.stream()
//	        			.map(s->(ocrAttmept.get(s).get(0).k() + ""))
//	        			.collect(Collectors.joining());
	        	
	        	String soFar="";
	        	Shape making=null;
	        	
	        	for(Shape s: sorted){
	        		String v=(ocrAttmept.get(s).get(0).k() + "");
	        		if(v.equals("-")){
	        			if(making!=null){
	        				System.out.println("Found:" + soFar);
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
	        		System.out.println("Found:" + soFar);
    				bestGuessOCR.put(making, soFar);
    			}
	        	
	        	
	        });
        
        List<Shape> ocrMeaningful=bestGuessOCR.keySet()
				   .stream()
				   .filter(s->interpretOCRStringAsAtom(bestGuessOCR.get(s))!=null)
				   .collect(Collectors.toList());
        
        for(Shape s: bestGuessOCR.keySet()){
        	String sym=bestGuessOCR.get(s);
        	BranchNode actual=this.interpretOCRStringAsAtom(sym);
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
        
        
        //final cleanup
        {
        	ctab.getDashLikeScoreForAllEdges(bitmap)
        	    .forEach(t->{
        	    	if(t.v()<MIN_ST_DEV_FOR_KEEPING_DASHED_LINES && t.k().getDashed()){
        	    		System.out.println("Dash score:" + t.v());
        	    		t.k().setDashed(false);
        	    	}
        	    });
        }
        
        double fbondlength=ctab.getAverageBondLength();
        
        for(Shape s: bestGuessOCR.keySet()){
        	String sym=bestGuessOCR.get(s);
        	BranchNode actual=this.interpretOCRStringAsAtom(sym);
        	if(actual!=null && actual.isRealNode()){
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
        				
        				
        				actual.forEachBranchNode((parN,curN)->{
        					 Node mpnode=pnode;
        					 if(parN!=null){
        						 mpnode=parentNodes.get(parN);
        					 }
        					
        					 Node n= ctab.addNode(curN.suggestedPoint)
	        					     .setSymbol(curN.getSymbol());
        					 ctab.addEdge(mpnode.getIndex(), n.getIndex(), 1);
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

	
	
	
}
