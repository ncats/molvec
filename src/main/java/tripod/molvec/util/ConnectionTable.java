package tripod.molvec.util;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.BinaryOperator;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import gov.nih.ncats.chemkit.api.Atom;
import gov.nih.ncats.chemkit.api.Bond;
import gov.nih.ncats.chemkit.api.Bond.BondType;
import gov.nih.ncats.chemkit.api.Bond.Stereo;
import gov.nih.ncats.chemkit.api.Chemical;
import gov.nih.ncats.chemkit.api.ChemicalBuilder;
import tripod.molvec.Bitmap;
import tripod.molvec.CachedSupplier;
import tripod.molvec.algo.Tuple;
import tripod.molvec.algo.Tuple.KEqualityTuple;
import tripod.molvec.util.ConnectionTable.Node;

public class ConnectionTable{
	private List<Node> nodes = new ArrayList<Node>();
	private List<Edge> edges = new ArrayList<Edge>();

	
	private static final int AROMATIC_ORDER=0xDE10CA1;
	
	private CachedSupplier<Map<Integer,List<Edge>>> _bondMap = CachedSupplier.of(()->_getEdgeMap());
	private CachedSupplier<Map<Node,Integer>> _nodeMap = CachedSupplier.of(()->_getNodeMap());
	
	
	public Node addNode(Point2D p){
		Node n=new Node(p,"C");
		nodes.add(n);
		resetCaches();
		return n;
	}
	
	public List<Node> setNodeToSymbol(Shape s, String sym){
		return nodes.stream()
		.filter(n->s.contains(n.point.getX(),n.point.getY()))
		.peek(n->{
			n.symbol=sym;
		})
		.collect(Collectors.toList());
		
	}
	public Shape getConvexHull(){
		return this.nodes.stream().map(n->n.getPoint()).collect(GeomUtil.convexHull());
	}
	
	public Shape getAreaAround(double rad){
		return this.edges.stream()
		          .map(e->e.getLine())
		          .map(l->GeomUtil.growShapeNPoly(l, rad, 10))
		          .collect(GeomUtil.union())
		          .orElse(null);
	}
	
	public List<ConnectionTable> getDisconnectedComponents(){
		
		//note: very inefficient
		return GeomUtil.groupThings(this.nodes, t->t.k().connectsTo(t.v()))
		.stream()
		.map(ln->{
			int[] map = new int[this.getNodes().size()];
			ConnectionTable ct  = new ConnectionTable();
			for(int i=0;i<ln.size();i++){
				Node oldNode=ln.get(i);
				
				Node newNode=ct.addNode(oldNode.getPoint()).setSymbol(oldNode.getSymbol()).setInvented(oldNode.isInvented());
				map[oldNode.getIndex()] = newNode.getIndex();
			}
			
			ln.stream()
			  .flatMap(n->n.getEdges().stream())
			  .distinct()
			  .forEach(e->{
				  int of=map[e.getRealNode1().getIndex()];
				  int os=map[e.getRealNode2().getIndex()];
				  Edge nedge=ct.addEdge(of,os,e.getOrder());
				  nedge.setDashed(e.getDashed());
				  nedge.setWedge(e.getWedge());
			  });
			return ct;
		})
		.collect(Collectors.toList());
		
	}
	
	
	public Chemical toChemical(){
		//return toChemical(this.getAverageBondLength(), false);
		return toChemical(1,true);
	}
	
	public Chemical toChemical(double averageBondLength, boolean center){
		
		AffineTransform at = new AffineTransform();
		double blcur = Math.max(this.getAverageBondLength(),1);
		
		double scale = averageBondLength/blcur;
		
		

		at.scale(scale, scale);
		if(center){
			if(this.getNodes().size()>0){
				Rectangle2D rect=this.getNodes().stream().map(n->n.getPoint()).collect(GeomUtil.convexHull()).getBounds2D();
				at.translate(-rect.getCenterX(), -rect.getCenterY());
			}
		}
		
		ChemicalBuilder cb = new ChemicalBuilder();
		Atom[] atoms = new Atom[nodes.size()];
		
		for(int i=0;i<nodes.size();i++){
			Node n = nodes.get(i);
			Point2D np = at.transform(n.getPoint(), null);
			
			atoms[i]=cb.addAtom(n.symbol,np.getX(),-np.getY());
			if(n.getCharge()!=0){
				atoms[i].setCharge(n.getCharge());
			}
		}
		cb.aromatize(true);
		
		//List<Tuple<String,String>> changeBonds = new ArrayList<>();
		
		List<Tuple<Bond,Tuple<String,String>>> addedAromaticBonds = new ArrayList<>();
		
		for(Edge e : edges){
			if(e.getOrder()==1){
				Bond b=cb.addBond(atoms[e.n1],atoms[e.n2],BondType.SINGLE);
				if(e.getDashed()){
					b.setStereo(Stereo.DOWN);
				}
				if(e.getWedge()){
					b.setStereo(Stereo.UP);
				}
				
			}else if(e.getOrder()==2){
				cb.addBond(atoms[e.n1],atoms[e.n2],BondType.DOUBLE);
			}else if(e.getOrder()==3){
				cb.addBond(atoms[e.n1],atoms[e.n2],BondType.TRIPLE);				
			}else if(e.isAromatic()){
				//doesn't work for some reason, so fix in the molfile and reload
				int n1=e.n1+1;
				int n2=e.n2+1;
				String s1=("   " + n1);
				String s2=("   " + n2);
				s1=s1.substring(s1.length()-3);
				s2=s2.substring(s2.length()-3);
				
				
				Tuple<String,String> trans=Tuple.of(s1 + s2 + "  1", s1 + s2 + "  4");
				
				
				
				Bond b=cb.addBond(atoms[e.n1],atoms[e.n2],BondType.AROMATIC);
				
				addedAromaticBonds.add(Tuple.of(b,trans));
			}else{
				//fall back to single
				cb.addBond(atoms[e.n1],atoms[e.n2],BondType.SINGLE);
				
			}
		}
		
		
		Chemical tc=cb.build();
		if(!addedAromaticBonds.isEmpty()){
			try {
				List<Bond> keepBond= GeomUtil.groupThings(addedAromaticBonds, t->{
					Bond b1=t.k().k();
					Bond b2=t.v().k();
					Atom a1i=b1.getAtom1();
					Atom a2i=b1.getAtom2();
					Atom b1i=b2.getAtom1();
					Atom b2i=b2.getAtom2();
					if(a1i==b1i||a1i==b2i||a2i==b1i||a2i==b2i){
						return true;
					}
					return false;
				})
				.stream()
				.filter(bl->bl.size()>4)
				.map(bl->{
					Map<Atom,AtomicInteger> acounts=new HashMap<>();
					
					bl.stream().forEach(bt->{
						acounts.computeIfAbsent(bt.k().getAtom1(), (k)->new AtomicInteger(0)).incrementAndGet();
						acounts.computeIfAbsent(bt.k().getAtom2(), (k)->new AtomicInteger(0)).incrementAndGet();
					});
					
					for(int i=bl.size()-1;i>=0;i--){
						Bond b=bl.get(i).k();
						Atom ai1=b.getAtom1();
						Atom ai2=b.getAtom2();
						if(acounts.getOrDefault(ai1, new AtomicInteger(0)).get()<2 || acounts.getOrDefault(ai2, new AtomicInteger(0)).get()<2){
							//not really a ring
							bl.remove(i);
							acounts.getOrDefault(ai1, new AtomicInteger(0)).decrementAndGet();
							acounts.getOrDefault(ai2, new AtomicInteger(0)).decrementAndGet();
							i=bl.size();
						}
					}
					
					return bl;
				})
				.peek(bl->System.out.println("Found bonds:" + bl.size()))
				.flatMap(bl->bl.stream())
				.map(t->t.k())
				.collect(Collectors.toList());
				
				List<Tuple<String,String>> realTransForm = addedAromaticBonds.stream()
																			 .filter(t->keepBond.contains(t.k()))
																			 .map(t->t.v())
																			 .collect(Collectors.toList());
				
				
				
				String nmol=Arrays.stream(tc.toMol().split("\n"))
				      .map(l->{
				    	  return realTransForm.stream()
				    	  			 .filter(lc->l.startsWith(lc.k()))
				    	  			 .findFirst()
				    	  			 .map(t->l.replace(t.k(), t.v()))
				    	  			 .orElse(l);
				      })
				      .collect(Collectors.joining("\n"));
				//System.out.println(nmol);
				return ChemicalBuilder.createFromMol(nmol, Charset.defaultCharset()).build();
			} catch (Exception e1) {

				e1.printStackTrace();
			}
		}
		return tc;
	}
	
	public ConnectionTable mergeNodesAverage(int n1, int n2){
		return mergeNodes(n1,n2,(node1,node2)->new Point2D.Double((node1.getX()+node2.getX())/2,(node1.getY()+node2.getY())/2));
	}
	
	public ConnectionTable mergeNodes(int n1, int n2, BinaryOperator<Point2D> op){
		mergeNodesGetTransform(n1,n2,op);
		return this;
	}
	
	public Map<Integer,Integer> mergeNodesGetTransform(int n1, int n2, BinaryOperator<Point2D> op){
		Point2D node1=nodes.get(n1).point;
		Point2D node2=nodes.get(n2).point;
		Point2D np = op.apply(node1, node2);
		int oldMax=nodes.size();
		
		int remNode=Math.max(n1, n2);
		int keepNode=Math.min(n1, n2);
		nodes.set(keepNode,new Node(np,nodes.get(n1).symbol));
		
		nodes.remove(remNode);
		for(Edge e: edges){
			if(e.n1==remNode){
				e.n1=keepNode;
			}
			if(e.n2==remNode){
				e.n2=keepNode;
			}
			if(e.n1>remNode){
				e.n1=e.n1-1;
			}
			if(e.n2>remNode){
				e.n2=e.n2-1;
			}
		}
		Map<Integer,Integer> oldToNew =new HashMap<>();
		for(int i=0;i<oldMax;i++){
			if(i<remNode){
				oldToNew.put(i, i);
			}
			if(i==remNode){
				oldToNew.put(i, keepNode);
			}
			if(i>remNode){
				oldToNew.put(i, i-1);
			}
		}
		resetCaches();
		return oldToNew;
	}
	
	public ConnectionTable mergeNodes(List<Integer> nlist, Function<List<Point2D>, Point2D> op){
		Point2D p = op.apply(nlist.stream().map(i->nodes.get(i).point).collect(Collectors.toList()));
		
		nlist=nlist.stream().sorted().collect(Collectors.toList());
		
		if(nlist.size()==1){
			int ni=nlist.get(0);
			nodes.get(ni).point=p;
		}else{
			for(int i=nlist.size()-1;i>=1;i--){
				int ni1=nlist.get(i);
				int ni2=nlist.get(i-1);
				mergeNodes(ni1,ni2,(p1,p2)->p);
			}
		}
		
		resetCaches();
		return this;
	}
	
	public ConnectionTable removeNode(int remNode){
		
		for(Edge e: edges){
			if(e.n1==remNode || e.n2==remNode){
				throw new IllegalStateException("Can't remove a node that is used in an edge");
			}
			if(e.n1>remNode){
				e.n1=e.n1-1;
			}
			if(e.n2>remNode){
				e.n2=e.n2-1;
			}
		}
		this.nodes.remove(remNode);
		resetCaches();
		return this;
	}
	
	public ConnectionTable removeEdge(Edge e) {
		this.edges.remove(e);
		resetCaches();
		return this;
	}
	
	public Map<Integer,Integer> mergeNodesGetTransform(List<Integer> nlist, Function<List<Point2D>, Point2D> op){
		Point2D p = op.apply(nlist.stream().map(i->nodes.get(i).point).collect(Collectors.toList()));
		
		nlist=nlist.stream().sorted().collect(Collectors.toList());
		
		Map<Integer,Integer> oldToNewMap=new HashMap<>();
		for(int i=0;i<this.nodes.size();i++){
			oldToNewMap.put(i, i);
		}
		
		
		if(nlist.size()==1){
			nodes.get(nlist.get(0)).point=p;
		}else{
			for(int i=nlist.size()-1;i>=1;i--){
				
				int ni1=nlist.get(i);
				int ni2=nlist.get(i-1);
				Map<Integer,Integer> oldNewNew=mergeNodesGetTransform(ni1,ni2,(p1,p2)->p);
				oldToNewMap=oldToNewMap.entrySet().stream()
				   .map(Tuple::of)
		//		   .peek(t->System.out.println("PFrom:"+t.k() + " to "+t.v()))
				   .map(Tuple.vmap(oi->oldNewNew.get(oi)))
		//		   .peek(t->System.out.println("From:"+t.k() + " to "+t.v()))
				   .collect(Tuple.toMap());
				
			}
		}
		return oldToNewMap;
	}
	
	public ConnectionTable cleanMeaninglessEdges(){
		this.edges=edges.stream()
				 		.filter(e->e.n1!=e.n2)
						.collect(Collectors.toList());

		resetCaches();
		
		return this;
	}
	
	public ConnectionTable standardCleanEdges(){
		this.cleanMeaninglessEdges();
        this.cleanDuplicateEdges((e1,e2)->{
        	if(e1.getOrder()>e2.getOrder()){
        		return e1;
        	}else if(e1.getOrder()<e2.getOrder()){
        		return e2;
        	}else{
        		if(e1.getDashed())return e2;
        		if(e2.getDashed())return e1;
        		if(e1.getWedge())return e1;
        		if(e2.getWedge())return e2;
        	}
        	return e1;
        });
        return this;
	}
	
	public ConnectionTable cleanDuplicateEdges(BinaryOperator<Edge> combiner){
		edges=edges.stream().map(e->e.standardize())
		              .map(e->Tuple.of(e,e.n1+"_" + e.n2))
		              .collect(Collectors.toMap(t->t.v(),t-> t.k(), combiner))
		              .values()
		              .stream()
		              .collect(Collectors.toList());
		resetCaches();
		return this;
	}
	public ConnectionTable mergeNodesCloserThan(double maxDistance){
			return mergeFilteredNodesCloserThan(maxDistance,(n->true));
	}
	public ConnectionTable mergeFilteredNodesCloserThan(double maxDistance, Predicate<Node> includeNode){
		
//		return mergeNodesCloserThan(maxDistance, n->true,nl->nl.stream().map(n->n.getPoint()).collect(GeomUtil.averagePoint()));
		
		boolean mergedOne = true;
		
		while(mergedOne){
			mergedOne=false;
			for(int i=0;i<nodes.size();i++){
				
				Node node1=nodes.get(i);
				if(!includeNode.test(node1))continue;
				Point2D pnti=node1.point;
				for(int j=i+1;j<nodes.size();j++){
					Point2D pntj = nodes.get(j).point;
					if(!includeNode.test(nodes.get(j)))continue;
					if(pnti.distance(pntj)<maxDistance){
						mergeNodesAverage(i,j);
						mergedOne=true;
						break;
					}
				}
				if(mergedOne)break;
			}
		}
		resetCaches();

		return this;
	}
	
	
	
	public ConnectionTable mergeNodesCloserThan(double maxDistance, Function<List<Node>, Point2D> combiner){
		return mergeNodesCloserThan(maxDistance,n->true,combiner);
	}
	
	public ConnectionTable mergeNodesCloserThan(double maxDistance, Predicate<Node> incudeNode,Function<List<Node>, Point2D> combiner){
		List<List<Node>> mergeList=GeomUtil.groupThings(nodes.stream().filter(incudeNode).collect(Collectors.toList()), (t)->{
								Node n1=t.k();
								Node n2=t.v();
								return n1.point.distance(n2.point)<maxDistance;
							})
							.stream()
							.filter(nl->nl.size()>=2)
							.collect(Collectors.toList());
		
		mergeList.forEach(ml->{
			List<Integer> nindexs = ml.stream().map(n->n.getIndex()).collect(Collectors.toList());
			//System.out.println("NIA:" + nindexs);
			Point2D centerpt = combiner.apply(ml);
			if(centerpt!=null){
				this.mergeNodes(nindexs, (pl)->{
					return centerpt;
				});			
			}
		});
		return this;
	}
	
	public Edge addEdge(int n1, int n2, int o){
		this.edges.add(new Edge(n1,n2,o));
		resetCaches();
		return this.edges.get(this.edges.size()-1);
	}
	
	public List<Node> getNodesInsideShape(Shape s, double tol){
		List<Node> mnodes= new ArrayList<>();
		
		for(int i=nodes.size()-1;i>=0;i--){
			Point2D pn = nodes.get(i).point;
			if(GeomUtil.distanceTo(s,pn)<tol || s.contains(pn)){
				mnodes.add(nodes.get(i));
			}
		}
		return mnodes;
	}
	
	public Tuple<Node,Double> getClosestNodeToShape(Shape s){
		return nodes.stream()
		     .map(n->Tuple.of(n,GeomUtil.distanceTo(s, n.point)).withVComparator())
		     .min(CompareUtil.naturalOrder())
		     .orElse(null);
		
	}
	
	public List<Tuple<Line2D,Integer>> asBondOrderLines(){
		return this.edges.stream().map(e->Tuple.of(e.getLine(),e.getOrder())).collect(Collectors.toList());
	}
	
	
	public static ConnectionTable fromLinesAndOrders(List<Tuple<Line2D,Integer>> lines){
		ConnectionTable ct=new ConnectionTable();
		
		//Collections.shuffle(lines);
		
		for(int i=0;i<lines.size();i++){
			Tuple<Line2D, Integer> lineOrder1=lines.get(i);
			ct.addNode(lineOrder1.k().getP1());
			ct.addNode(lineOrder1.k().getP2());
			ct.addEdge(i*2, i*2+1, lineOrder1.v());
		}
		return ct;
	}
	

	public ConnectionTable getNewConnectionTable(List<Shape> likelyNodes,
			double maxDistanceRatioNonLikely, 
			double maxDistanceRatioLikely, 
			double maxDistanceRatioPerLine,
			double minPerLineDistanceRatioForIntersection,
			double maxCandidateRatioForIntersectionWithNeighbor,
			Predicate<Line2D> acceptNewLine) {
		
			ConnectionTable copy= this.cloneTab();
			
			//might make things ugly, consider changing
			copy.mergeNodesCloserThan(3.0);
		
			
			
			List<Tuple<Line2D,Integer>> lines = copy.asBondOrderLines();
			List<Tuple<Integer,Integer>> ncount = copy.edges.stream().map(e->Tuple.of(e.getRealNode1().getEdgeCount(),e.getRealNode2().getEdgeCount())).collect(Collectors.toList());
				
				for(int i=0;i<lines.size();i++){
					Tuple<Line2D, Integer> lineOrder1=lines.get(i);
					double distance1 = GeomUtil.length(lineOrder1.k());
					for(int j=i+1;j<lines.size();j++){
						
						
						Tuple<Line2D, Integer> lineOrder2=lines.get(j);
						double distance2 = GeomUtil.length(lineOrder2.k());
						double totalDistance = distance1+distance2;
						
						Point2D intersect = GeomUtil.intersection(lineOrder1.k(),lineOrder2.k());
						if(intersect==null)continue;
						double ndistance1 = Math.max(intersect.distance(lineOrder1.k().getP1()), intersect.distance(lineOrder1.k().getP2()));
						double ndistance2 = Math.max(intersect.distance(lineOrder2.k().getP1()), intersect.distance(lineOrder2.k().getP2()));
						double totalDistanceAfter = ndistance1+ndistance2;
						
						double ratioTotal = Math.max(totalDistanceAfter,totalDistance)/Math.min(totalDistanceAfter, totalDistance);
						
						double ratioLine1 = Math.max(ndistance1,distance1)/Math.min(ndistance1, distance1);
						double ratioLine2 = Math.max(ndistance2,distance2)/Math.min(ndistance2, distance2);
						
						double ratioOldToNew1 = ndistance1/distance1;
						double ratioOldToNew2 = ndistance2/distance2;
						
						boolean merge = false;
						
					
								
						if(ratioTotal<maxDistanceRatioLikely){
							if(ratioTotal<maxDistanceRatioNonLikely){
								if(ratioLine1<maxDistanceRatioPerLine && ratioLine2<maxDistanceRatioPerLine &&
								   ratioOldToNew1>minPerLineDistanceRatioForIntersection && ratioOldToNew2>minPerLineDistanceRatioForIntersection){
									merge=true;
								}
							}else{
								boolean inLikelyNode=likelyNodes.stream().filter(s->s.contains(intersect)).findAny().isPresent();
								if(inLikelyNode){
									
//									System.out.println("Something");
//									System.out.println("Something2");
//									
									if(ratioTotal>maxCandidateRatioForIntersectionWithNeighbor){
										
										Tuple<Integer,Integer> n1count=ncount.get(i);
										Tuple<Integer,Integer> n2count=ncount.get(j);
										int line1MergeSide=0;
										int line2MergeSide=0;
										if(intersect.distance(lineOrder1.k().getP1())<intersect.distance(lineOrder1.k().getP2())){
											line1MergeSide=n1count.k();
										}else{
											line1MergeSide=n1count.v();
										}
										if(intersect.distance(lineOrder2.k().getP1())<intersect.distance(lineOrder2.k().getP2())){
											line2MergeSide=n2count.k();
										}else{
											line2MergeSide=n2count.v();
										}
										if(line1MergeSide>1 || line2MergeSide>1){
											merge=false;
										}else{
											merge=true;
										}
									}else{
										merge=true;	
									}
									
									
								}
							}
						}
						
						if(merge){
							Line2D newLine1 = GeomUtil.longestLineFromOneVertexToPoint(lineOrder1.k(),intersect);
							Line2D newLine2 = GeomUtil.longestLineFromOneVertexToPoint(lineOrder2.k(),intersect);
							if(acceptNewLine.test(newLine1) && acceptNewLine.test(newLine2)){
								Tuple<Line2D,Integer> norder1 = Tuple.of(newLine1,lineOrder1.v());
								Tuple<Line2D,Integer> norder2 = Tuple.of(newLine2,lineOrder2.v());
								lines.set(i, norder1);
								lines.set(j, norder2);
								lineOrder1 = norder1;
							}
						}				
					}
				}
					
				return fromLinesAndOrders(lines);
	}
	
	
	public ConnectionTable mergeAllNodesInside(Shape s, double tol,Predicate<Node> allow, Function<List<Point2D>, Point2D> merger){
		
		List<Integer> toMerge = getAllNodeIndexesInsideShape(s,tol);
		toMerge=toMerge.stream()
		       .filter(i->allow.test(nodes.get(i)))
		        .collect(Collectors.toList());
			
		return this.mergeNodes(toMerge, merger);
	}
	
	private List<Integer> getAllNodeIndexesInsideShape(Shape s,double tol){
		List<Integer> toMerge = new ArrayList<Integer>();
		for(int i=nodes.size()-1;i>=0;i--){
			Point2D pn = nodes.get(i).point;
			if(GeomUtil.distanceTo(s,pn)<=tol){
				toMerge.add(i);
			}
		}
		return toMerge;
	}
	public List<Node> getAllNodesInsideShape(Shape s,double tol){
		return getAllNodeIndexesInsideShape(s,tol).stream().map(i->this.nodes.get(i)).collect(Collectors.toList());
	}
	
	
	public List<Tuple<Edge,Tuple<Node,Node>>> getAllEdgesEntering(Shape s, double tol){
		
		List<Node> toMerge = new ArrayList<Node>();
		for(int i=nodes.size()-1;i>=0;i--){
			Point2D pn = nodes.get(i).point;
			if(GeomUtil.distanceTo(s,pn)<tol){
				toMerge.add(nodes.get(i));
			}
		}
		return toMerge.stream()
				.flatMap(n->n.getEdges().stream())
				.distinct()
				.map(e->{
					Node n1=e.getRealNode1();
					Node n2=e.getRealNode2();
					boolean containsN1=toMerge.contains(n1);
					boolean containsN2=toMerge.contains(n2);
					if(containsN1 && containsN2){
						return Tuple.of(e,Tuple.of(n1,(Node)null));
					}else if(containsN1){
						return Tuple.of(e,Tuple.of(n1,n2));
					}else if(containsN2){
						return Tuple.of(e,Tuple.of(n2,n1));
					}
					return Tuple.of(e,Tuple.of((Node)null,(Node)null));
				})
				.filter(t->t.v().k()!=null)
				.filter(t->t.v().v()!=null)
		        .collect(Collectors.toList());
			
	}
	
	public ConnectionTable mergeAllNodesInsideCenter(Shape s, double tol){
		Rectangle2D r=s.getBounds2D();
		Point2D p = new Point2D.Double(r.getCenterX(),r.getCenterY());
		return mergeAllNodesInside(s,tol,n->true,(l)->p);
	}
	
	public ConnectionTable mergeNodesExtendingTo(Collection<Shape> shapes,double maxAvgBondRatio, double maxTotalAvgBondRatio){
		double avg = this.getAverageBondLength();
		edges.stream()
		     //.filter(e->e.getBondLength()<avg)
		     .forEach(e->{
		    	 boolean term = (e.getRealNode1().getEdgeCount()==1) ||
		    			        (e.getRealNode2().getEdgeCount()==1);
		    	 
		    	 double maxR=maxAvgBondRatio;
		    	 double maxTR=maxTotalAvgBondRatio;
		    	 
		    	 if(!term)maxR = 1;
		    	 if(!term)maxTR = 1;
		    	 
		    	 Point2D p1=e.getPoint1();
		    	 Point2D p2=e.getPoint2();
		    	 double minDistp1=Double.MAX_VALUE;
		    	 double minDistp2=Double.MAX_VALUE;
		    	 Shape closest1=null;
		    	 Shape closest2=null;
		    	 
		    	 for(Shape s: shapes){
		    		 double d1=GeomUtil.distanceTo(s, p1);
		    		 double d2=GeomUtil.distanceTo(s, p2);
		    		 if(d1<minDistp1){
		    			 minDistp1=d1;
		    			 closest1=s;
		    		 }
		    		 if(d2<minDistp2){
		    			 minDistp2=d2;
		    			 closest2=s;
		    		 }
		    	 }
		    	 if(minDistp1>avg*maxR && minDistp2>avg*maxR){
		    		 return;
		    	 }
		    	 Point2D newPoint1=e.getPoint1();
		    	 Point2D newPoint2=e.getPoint2();
		    	 Line2D newLine = null;
		    	 boolean onlyOne=true;
		    	 if(minDistp1<avg*maxR && minDistp2<avg*maxR && 
		    		closest1!=closest2){
		    		 onlyOne=false;
		    		newPoint1=GeomUtil.findCenterOfShape(closest1);
		    		newPoint2=GeomUtil.findCenterOfShape(closest2);
		    		newLine=new Line2D.Double(newPoint1,newPoint2);
		    		double nl=GeomUtil.length(newLine);
		    		if(nl> maxTR*avg ){
		    			
			    		 onlyOne=true;
			    	 }
			    	 double cosTheta = Math.abs(GeomUtil.cosTheta(newLine,e.getLine()));
			    	 if(cosTheta<Math.cos(12.0*Math.PI/180.0)){
			    		 onlyOne=true;
			    	 }
			    	 
			    	 Point2D np=GeomUtil.findCenterOfShape(newLine);
			    	 Point2D op=GeomUtil.findCenterOfShape(e.getLine());
			    	 if(np.distance(op)>nl*0.2){
			    		 onlyOne=true;
			    	 }
		    				    		 
		    	 }
		    	 
		    	 if(onlyOne){
			    	 //Node closestNode = this.nodes.get(e.n1);
			    	 //Point2D newPoint = null;
			    	// Line2D newLine = null;
			    	 if(minDistp1>minDistp2){
			    	//	 closestNode = this.nodes.get(e.n2);
			    		 newPoint2=GeomUtil.findCenterOfShape(closest2);
			    		 newPoint1=e.getPoint1();
			    		 //newLine=new Line2D.Double(p1,newPoint);
			    	 }else{
			    		 newPoint1=GeomUtil.findCenterOfShape(closest1);
			    		 newPoint2=e.getPoint2();
			    		 //newLine=new Line2D.Double(p2,newPoint);
			    	 }
		    	 }
		    	 newLine=new Line2D.Double(newPoint1,newPoint2);
		    	 if(GeomUtil.length(newLine)> maxTR*avg){
		    		 return;
		    	 }
		    	 double cosTheta = Math.abs(GeomUtil.cosTheta(newLine,e.getLine()));
		    	 if(cosTheta<Math.cos(12.0*Math.PI/180.0)){
		    		 return;
		    	 }
		    	 //closestNode.point=newPoint;
		    	 
		    	 this.nodes.get(e.n1).setPoint(newPoint1);
		    	 this.nodes.get(e.n2).setPoint(newPoint2);
		    	 
		     });
		this.mergeNodesCloserThan(avg/40);
		
		return this;
	}
	
	public ConnectionTable mergeAllNodesOnParLines(){
		Map<Line2D,Edge> edgeMap = this.edges.stream().collect(Collectors.toMap(e->e.getLine(),e->e ));
		
		Map<Integer,Integer> oldToNewMap = IntStream.range(0,this.nodes.size())
		         .mapToObj(i->i)
		         .collect(Collectors.toMap(i->i, i->i));
		
		List<LinkedHashSet<Integer>> mergeNodes=GeomUtil.groupMultipleBonds(edgeMap.keySet().stream().collect(Collectors.toList()),5*Math.PI/180, 2, .8, 0)
		.stream()
		.filter(l->l.size()>1)
		.map(l->{
			Line2D keep=l.stream()
			 .map(l1->Tuple.of(GeomUtil.length(l1),l1).withKComparator())
			 .max(CompareUtil.naturalOrder())
			 .get()
			 .v();
			Edge keepEdge=edgeMap.get(keep);
			LinkedHashSet<Integer> mergeNode1 = new LinkedHashSet<Integer>();
			LinkedHashSet<Integer> mergeNode2 = new LinkedHashSet<Integer>();
			mergeNode1.add(keepEdge.n1);
			mergeNode2.add(keepEdge.n2);
			
			Point2D node1Point = keepEdge.getPoint1();
			Point2D node2Point = keepEdge.getPoint2();
			
			l.stream()
			 .map(l1->edgeMap.get(l1))
			 .filter(e->e!=keepEdge)
			 .forEach(me->{
				 double distance1to1=me.getPoint1().distance(node1Point);
				 double distance1to2=me.getPoint1().distance(node2Point);
				 double distance2to1=me.getPoint2().distance(node1Point);
				 double distance2to2=me.getPoint2().distance(node2Point);
				 if(distance1to1<distance1to2){
					 mergeNode1.add(me.n1);
				 }else{
					 mergeNode2.add(me.n1);
				 }
				 if(distance2to1<distance2to2){
					 mergeNode1.add(me.n2);
				 }else{
					 mergeNode2.add(me.n2);
				 }					 
			 });				
			List<LinkedHashSet<Integer>> toMerge=new ArrayList<LinkedHashSet<Integer>>();
			toMerge.add(mergeNode1);
			toMerge.add(mergeNode2);
			return toMerge;
		})
		.flatMap(l->l.stream())
		.filter(l->l.size()>1)
		.collect(Collectors.toList());
		
		mergeNodes.forEach(ls->{
			List<Integer> toMerge=ls.stream().map(i->oldToNewMap.get(i)).collect(Collectors.toList());
			
			Point2D keeper=this.nodes.get(toMerge.get(0)).point;
			toMerge=toMerge.stream().distinct().collect(Collectors.toList());
			if(toMerge.size()<=1)return;

			Map<Integer,Integer> newTrans=mergeNodesGetTransform(toMerge,(pts)->keeper);
			for(int i : oldToNewMap.keySet()){
				int oldMap=oldToNewMap.get(i);

				//System.out.println(newTrans.toString());
				int newMap=newTrans.get(oldMap);
				oldToNewMap.put(i, newMap);
			}
		});
		return this;
	}
	
	public ConnectionTable createNodesOnIntersectingLines(double tol, Predicate<List<Edge>> shouldSplitEdges, Consumer<Node> newNodeConsumer){
		
		boolean splitOne =true;
		while(splitOne){
			splitOne=false;
		
			for(int i=0;i<edges.size();i++){
				Edge e1=edges.get(i);
				Set<Integer> i1set=e1.getNodeSet();
				for(int j=0;j<edges.size();j++){
					Edge e2=edges.get(j);
					//can't share a neighbor
					if(i1set.contains(e2.n1) || i1set.contains(e2.n2)){
						continue;
					}
					Point2D np=GeomUtil.intersection(e1.getLine(), e2.getLine());
					if(np==null)continue;
					
					
					double dist1=e1.getLine().ptSegDist(np);
					double dist2=e2.getLine().ptSegDist(np);
					if(dist1<tol && dist2<tol){
						if(e1.getPoint1().distance(np) <tol ||e1.getPoint2().distance(np) <tol){
							if(e2.getPoint1().distance(np) <tol ||e2.getPoint2().distance(np) <tol){
								continue;
							}
						}
						
							//Point2D np=GeomUtil.intersection(e1.getLine(),e2.getLine());
							int nodeNew = this.nodes.size();
							Node rnode=this.addNode(np);
							Edge nedge1 = new Edge(e1.n1, nodeNew, e1.getOrder());
							Edge nedge2 = new Edge(e1.n2, nodeNew, e1.getOrder());
							Edge nedge3 = new Edge(e2.n1, nodeNew, e2.getOrder());
							Edge nedge4 = new Edge(e2.n2, nodeNew, e2.getOrder());
							edges.set(i, nedge1);
							edges.set(j, nedge2);
							edges.add(nedge3);
							edges.add(nedge4);
							if(!shouldSplitEdges.test(Stream.of(nedge1,nedge2,nedge3,nedge4).collect(Collectors.toList()))){
								edges.set(i, e1);
								edges.set(j, e2);
								edges.remove(edges.size()-1);
								edges.remove(edges.size()-1);
								this.removeNode(nodeNew);
								//this.nodes.remove(nodeNew);
							}else{
								splitOne=true;
								newNodeConsumer.accept(rnode);
							}
							break;
						
					}
				}
				if(splitOne)break;
			}
		}
		resetCaches();
		return this;
	}
	
	private void resetCaches(){
		_bondMap.resetCache();
		_nodeMap.resetCache();
		_averageBondLength.resetCache();
		
	}
	
	public double getAverageBondLength(){
		return _averageBondLength.get();
	}
	
	private CachedSupplier<Double> _averageBondLength=CachedSupplier.of(()->{
			return getMeanBondLength();
	});
	
	public double getMeanBondLength(){
		return edges.stream().mapToDouble(e->e.getEdgeLength()).average().orElse(0);
	}
	
	public double getMedianBondLength(){
		
		double[] lens= edges.stream()
				    	.mapToDouble(e->e.getEdgeLength())
				    	.sorted()
				    	.toArray();
		
		if(lens.length==0)return 0;
		if(lens.length %2==1){
			//1 -> 0
			//3 -> 1
			//5 -> 2
			return lens[lens.length/2];
		}else{
			//2 -> 1,0
			//4 -> 2,
			
			return (lens[lens.length/2] + lens[lens.length/2 - 1] / 2);
		}
	}
	
	public class Node{
		private Point2D point;
		private String symbol="C";
		private int charge=0;
		private boolean invented=false;
		
		
		public List<KEqualityTuple<Node,Edge>> getNeighborNodes(){
			return this.getEdges()
			    .stream()
				.map(ne->Tuple.of(ne.getOtherNode(this),ne).withKEquality())
				.collect(Collectors.toList());
		}
		
		public boolean connectsTo(Node v) {
			return this.getNeighborNodes().stream().filter(t->t.k()==v).findAny().isPresent();
		}

		public int getCharge(){
			return this.charge;
		}
		
		public boolean isInRing(int maxRing){
			Set<Node> neighborsStart =getNeighborNodes().stream().map(t->t.k()).collect(Collectors.toSet());
			Set<Node> neighborsNext =new HashSet<>();
			List<Node> dontInclude =new ArrayList<>();
			dontInclude.add(this);
			
			for(int i=0;i<maxRing;i++){
				for(Node nn: neighborsStart){
					nn.getNeighborNodes().stream()
										 .map(t->t.k())
					                     .filter(n2->!dontInclude.contains(n2))
					                     .forEach(n2->{
					                    	 neighborsNext.add(n2);
					                     });
				}
				
				if(neighborsNext.contains(this)){
					return true;
				}
				dontInclude.clear();
				dontInclude.addAll(neighborsStart);
				neighborsStart= new HashSet<>(neighborsNext);
				neighborsNext.clear();
			}
			return false;
			
			
			
			
		}
		
		public Node setCharge(int c){
			this.charge=c;
			return this;
		}
		
		
		public Node(Point2D p, String s){
			this.point=p;
			this.symbol=s;				
		}
		public double distanceTo(Node n2){
			return this.point.distance(n2.point);
		}
		
		public List<Edge> getEdges(){
			return getEdgeMap().getOrDefault(getIndex(), new ArrayList<>());
		}
		
		public int getIndex(){
			Integer ind=getNodeMap().get(this);
			if(ind==null){
				throw new IllegalStateException("Can't find node index for node");
				//return -1;
			}
			return ind;
		}
		
		
		public Point2D getPoint() {
			return this.point;
		}
		public String getSymbol() {
			return this.symbol;
		}
		public Node setSymbol(String symbol2) {
			this.symbol=symbol2;
			return this;
			
		}

		public int getEdgeCount() {
			return this.getEdges().size();
		}

		public Node setPoint(Point2D ppnt) {
			this.point=ppnt;
			return this;
			
		}

		public Optional<Edge> getBondTo(Node v) {
			return this.getEdges()
			    .stream()
			    .filter(e->e.getOtherNode(this) == v)
			    .findFirst();
		}

		public Node setInvented(boolean b) {
			this.invented=true;
			return this;
		}
		
		public boolean isInvented() {
			return this.invented;
		}

		public int getValanceTotal() {
			return getEdges().stream().mapToInt(e->e.getOrder()).sum();
		}
		
	}
	public class Edge{
		int n1;
		int n2;
		private int order;
		boolean isWedge=false;
		boolean isDash=false;
		public Edge(int n1, int n2, int o){
			this.n1=n1;
			this.n2=n2;
			this.setOrder(o);
		}
		
		public boolean isAromatic() {
			return this.order==AROMATIC_ORDER;
		}

		public Edge setDashed(boolean d){
			this.isDash=d;
			return this;
		}
		public Edge setWedge(boolean d){
			this.isWedge=d;
			return this;
		}
		
		public boolean getDashed(){
			return this.isDash;
		}
		public boolean getWedge(){
			return this.isWedge;
		}
		
		
		
		public double getEdgeLength(){
			return ConnectionTable.this.nodes.get(n1).point.distance(ConnectionTable.this.nodes.get(n2).point);
		}
		public Edge standardize(){
			if(n2<n1){
				int t=this.n1;
				this.n1=this.n2;
				this.n2=t;
			}
			return this;
		}
		public Set<Integer> getNodeSet(){
			Set<Integer> iset = new HashSet<Integer>();
			iset.add(n1);
			iset.add(n2);
			return iset;
			
		}
		public Line2D getLine(){
			Point2D p1 = ConnectionTable.this.nodes.get(n1).point;
			Point2D p2 = ConnectionTable.this.nodes.get(n2).point;
			return new Line2D.Double(p1, p2);
		}
		public Point2D getPoint1(){
			return ConnectionTable.this.nodes.get(n1).point;
		}
		public Point2D getPoint2(){
			return ConnectionTable.this.nodes.get(n2).point;
		}
		
		public Node getOtherNode(Node n){
			if(n.getIndex() == this.n1)return getRealNode2();
			if(n.getIndex() == this.n2)return getRealNode1();
			return null;
		}
		
		public Node getRealNode1(){
			return ConnectionTable.this.nodes.get(n1);
		}
		public Node getRealNode2(){
			return ConnectionTable.this.nodes.get(n2);
		}
		
		
		public int getOrder() {
			
			return this.order;
		}

		public Edge switchNodes() {
			int t=this.n1;
			this.n1=this.n2;
			this.n2=t;
			return this;
			
		}

		public Edge setOrder(int order) {
			this.order = order;
			return this;
		}
		
		public Edge setToAromatic(){
			return setOrder(AROMATIC_ORDER);
		}
		
		public String toString(){
			return "Edge: " + this.n1 + " to " + this.n2 + ", order=" + this.order + ", dash=" + this.isDash + ", wedge=" + this.isWedge;
		}

		public boolean isInventedBond() {
			return this.getRealNode1().isInvented() || this.getRealNode2().isInvented();
		}

		public List<Edge> getNeighborEdges() {
			
			return Stream.concat(this.getRealNode1().getEdges().stream(),
					      this.getRealNode2().getEdges().stream())
			      .filter(e->e!=this)
			      .collect(Collectors.toList());
			      
			
		}

		public boolean hasNode(Node n) {
			if(this.getRealNode1()==n||this.getRealNode2()==n)return true;
			return false;
		}
		
	}
	public Optional<Edge> getEdgeBetweenNodes(Node n1, Node n2){
		List<Edge> edges1=n1.getEdges();
		List<Edge> edges2=n2.getEdges();
		
		return edges1.stream()
		      .filter(e->edges2.contains(e))
		      .findFirst();
		      
		
	}
	public Optional<Edge> getEdgeBetweenNodes(int n1, int n2){
		return getEdgeBetweenNodes(this.nodes.get(n1),this.nodes.get(n2));
	}
	public void draw(Graphics2D g2) {
		int sx=1;
		Stroke old = g2.getStroke();
		
		Stroke dashed = new BasicStroke(3, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0, new float[]{3}, 0);
		Stroke wedge = new BasicStroke(5f);
        
		
		nodes.stream().map(n->n.point)
		.forEach(p->{
			g2.draw(new Ellipse2D.Double((p.getX()-2f/sx), (p.getY()-2f/sx), 4f/sx, 4f/sx));
		});
		
		edges.forEach(l->{
			if(l.getOrder()==1){
				g2.setPaint(Color.BLACK);
				if(l.isDash){
					g2.setStroke(dashed);						
				}
				if(l.isWedge){
					g2.setStroke(wedge);
				}
				
			}else if(l.getOrder()==2){
				g2.setPaint(Color.RED);
			}else if(l.getOrder()==3){
				g2.setPaint(Color.GREEN);
			}
			
			g2.draw(l.getLine());
			g2.setStroke(old);
		});
		
	}

	public List<Edge> getEdges() {
		return this.edges;
	}
	
	public List<Tuple<Edge,Double>> getTolerancesForAllEdges(Bitmap bm,List<Shape> shapes){
		return this.edges.stream()
		          .map(e->Tuple.of(e,getToleranceForEdge(e,bm,shapes)))
		          .collect(Collectors.toList());
	}
	
	
	public double getToleranceForEdge(Edge e, Bitmap bm, Collection<Shape> shapes){
		return GeomUtil.getLongestLineNotInside(e.getLine(), shapes)
		        	   .map(ll->bm.getLineLikeScore(e.getLine()))
		               .orElse(0.0);
	}
	
	public double getAverageToleranceForNode(Node n, Bitmap bm,List<Shape> shapes){
		return n.getEdges().stream()
		          .map(e->Tuple.of(e,getToleranceForEdge(e,bm,shapes)))
		          .mapToDouble(t->t.v())
		          .average()
		          .orElse(0);
	}
	
	public Tuple<Edge,Double> getWorstToleranceForNode(Node n, Bitmap bm,Collection<Shape> shapes){
		return n.getEdges().stream()
		          .map(e->Tuple.of(e,getToleranceForEdge(e,bm,shapes)))
		          .map(t->t.withVComparator())
		          .max(CompareUtil.naturalOrder())
		          .orElse(null);
	}
	
	public List<Tuple<Edge,Double>> getDashLikeScoreForAllEdges(Bitmap bm,Collection<Shape> shapes){
		return this.edges.stream()
		          .map(e->Tuple.of(e,getToleranceForEdge(e,bm,shapes)))
		          .collect(Collectors.toList());
	}

	public ConnectionTable makeMissingBondsToNeighbors(Bitmap bm, double d, double tol, Collection<Shape> OCRSet, double ocrTol,Consumer<Tuple<Double,Edge>> econs) {
		double avg=this.getAverageBondLength();
		Map<Integer,Set<Integer>> nmap = new HashMap<>();
		
		Set<Integer> nullSet = new HashSet<Integer>();
		
		
		edges.forEach(e->{
			nmap.computeIfAbsent(e.n1,k->new HashSet<>()).add(e.n2);
			nmap.computeIfAbsent(e.n2,k->new HashSet<>()).add(e.n1);
		});
		for(int i =0 ;i<this.nodes.size();i++){
			Node n1=nodes.get(i);
			for(int j =i+1 ;j<this.nodes.size();j++){
				if(!nmap.getOrDefault(i,nullSet).contains(j)){
					Node n2=nodes.get(j);
					if(n1.distanceTo(n2)<avg*d){
						Tuple<Shape,Double> t1=GeomUtil.findClosestShapeTo(OCRSet, n1.point);
						Tuple<Shape,Double> t2=GeomUtil.findClosestShapeTo(OCRSet, n2.point);
						
						Point2D pt1=n1.point;
						Point2D pt2=n2.point;
						
						List<Line2D> lines = new ArrayList<Line2D>();
						lines.add(new Line2D.Double(pt1,pt2));
						
						
						if(t1!=null){
							lines=lines.stream()
									   .map(l->GeomUtil.getLinesNotInside(l, t1.k()))
									   .flatMap(ll->ll.stream())
									   .collect(Collectors.toList());
						}
						if(t2!=null){
							lines=lines.stream()
									   .map(l->GeomUtil.getLinesNotInside(l, t2.k()))
									   .flatMap(ll->ll.stream())
									   .collect(Collectors.toList());
						}
						
						
						Line2D tline=lines.stream().map(l->Tuple.of(l,GeomUtil.length(l)).withVComparator())
						                .max(Comparator.naturalOrder())
						                .map(t->t.k())
						                .orElse(null);
						
						
//						if(t1!=null && t1.v()<ocrTol){
//							pt1=GeomUtil.closestPointOnShape(t1.k(), pt2);
//						}
//						if(t2!=null && t2.v()<ocrTol){
//							pt2=GeomUtil.closestPointOnShape(t2.k(), pt1);
//						}
//						
//						tline=new Line2D.Double(pt1,pt2);
						if(tline!=null){
							Edge e= new Edge(i,j,1);
							double score=bm.getLineLikeScore(tline);
							//System.out.println("Score:" + score);
							if(score<tol){
								this.edges.add(e);
								econs.accept(Tuple.of(score,e));
							}
						}
					}
				}
			}
		}
		resetCaches();
		return this;
	}

	public ConnectionTable removeOrphanNodes() {
		Map<Integer,List<Edge>> nmap = getEdgeMap();
		
		for(int i =nodes.size()-1;i>=0;i--){
			if(nmap.get(i).isEmpty()){
				removeNode(i);
			}
		}
		return this;
	}
	
	public ConnectionTable removeNodeAndEdges(Node n){
		List<Edge> nedges=n.getEdges();
		this.edges.removeAll(nedges);
		resetCaches();
		return this.removeNode(n.getIndex());
	}
	
	private Map<Integer,List<Edge>> _getEdgeMap(){
		Map<Integer,List<Edge>> nmap = new HashMap<>();
		IntStream.range(0, nodes.size())
		         .forEach(i->{
		        	 nmap.put(i, new ArrayList<>());
		         });
		edges.forEach(e->{
			nmap.computeIfAbsent(e.n1,k->new ArrayList<>()).add(e);
			nmap.computeIfAbsent(e.n2,k->new ArrayList<>()).add(e);
		});
		return nmap;
	}
	public Map<Integer,List<Edge>> getEdgeMap(){
		return _bondMap.get();
	}
	public Map<Node,Integer> getNodeMap(){
		return _nodeMap.get();
	}
	public Map<Node,Integer> _getNodeMap(){
		return IntStream.range(0, nodes.size())
				 .mapToObj(i->i)
		         .collect(Collectors.toMap(i->nodes.get(i), i->i));
	}
	
	

	public ConnectionTable fixBondOrders(Collection<Shape> likelyOCR, double shortestRealBondRatio, Consumer<Edge> edgeCons) {
		// TODO Auto-generated method stub
		this.edges
		    .stream()
		    .forEach(e->{
		    	Shape s1=GeomUtil.getClosestShapeTo(likelyOCR,e.getPoint1());
		    	Shape s2=GeomUtil.getClosestShapeTo(likelyOCR,e.getPoint2());
		    	if(s1!=s2){
		    		if(s1.contains(e.getPoint1())){
		    			if(s2.contains(e.getPoint2())){
		    				Line2D line = e.getLine();
		    				
		    				Point2D pn1=GeomUtil.getIntersection(s1,line).orElse(null);
		    				Point2D pn2=GeomUtil.getIntersection(s2,line).orElse(null);
		    				if(pn1!=null && pn2!=null){
		    					double realDistance=pn1.distance(pn2);
		    					if(realDistance/e.getEdgeLength()<shortestRealBondRatio){
		    						edgeCons.accept(e);
		    					}
		    				}
		    			}
		    		}
		    	}
		    });
		return this;
	}

	public List<Node> getNodes() {
		return this.nodes;
	}
	
	public List<Node> getNodesNotInShapes(Collection<Shape> shapes, double tol){
		if(shapes.isEmpty())return nodes;
		return nodes.stream()
		     .map(n->Tuple.of(n,GeomUtil.getClosestShapeTo(shapes, n.point)))
		     .filter(t->GeomUtil.distanceTo(t.v(),t.k().point)>tol)
		     .map(t->t.k())
		     .collect(Collectors.toList());
	}

	public ConnectionTable makeMissingNodesForShapes(Collection<Shape> likelyOCR, double mAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL,
			double mIN_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL) {
		double avg=this.getAverageBondLength();
		List<Shape> addShapes=likelyOCR.stream() 
		         .map(oc->Tuple.of(getClosestNodeToShape(oc),oc))
		         .filter(t->t.k().v()>avg*mIN_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL)
				 .filter(t->t.k().v()<avg*mAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL)
				 .map(t->t.v())
				 .collect(Collectors.toList());
		
		for(Shape s:addShapes){
			this.addNode(GeomUtil.findCenterOfShape(s));
		}        
		
		return this;
	}

	public ConnectionTable cloneTab() {
		ConnectionTable ctab2 = new ConnectionTable();
		this.nodes.forEach(n->{
			ctab2.addNode(n.point);
			Node nnode=ctab2.nodes.get(ctab2.nodes.size()-1);
			nnode.symbol=n.symbol;
		});
		this.edges.forEach(e->{
			ctab2.addEdge(e.n1, e.n2,e.order);
			Edge nedge=ctab2.edges.get(ctab2.edges.size()-1);
			nedge.setDashed(e.getDashed());
			nedge.setWedge(e.getWedge());
		});
		
		
		return ctab2;
	}

	public List<Edge> getBondsThatCross(Node n1, Node n2) {
		
		Line2D tline = new Line2D.Double(n1.getPoint(),n2.getPoint());
		
		return this.edges.stream()
		          .filter(e->e.getRealNode1()!=n1 && e.getRealNode1()!=n2)
		          .filter(e->e.getRealNode2()!=n1 && e.getRealNode2()!=n2)
		          .filter(e->GeomUtil.segmentIntersection(e.getLine(), tline).isPresent())
		          .collect(Collectors.toList());
		
		
	}
	
	
	public ConnectionTable simpleClean(){
		this.nodes.stream()
		          .filter(n->n.getEdgeCount()==3)
		          .filter(n->n.getEdges().stream().filter(e->e.getOrder()>1).map(e->e.getOtherNode(n).getSymbol()).anyMatch(s->!"C".equals(s)))
		          .forEach(n->{
		        	  Point2D cpoint=n.getPoint();
		        	  Point2D[] pts=n.getNeighborNodes().stream()
		        	                      .map(nn->nn.k().getPoint())
		        	                      .toArray(i->new Point2D[i]);
		        	  Shape tri=GeomUtil.convexHull2(pts);
		        	  
		        	  if(tri.contains(cpoint)){
		        		 Point2D cpt=GeomUtil.findCenterOfVertices(Arrays.asList(pts));
		        		 Point2D pp1=GeomUtil.projectPointOntoLine(new Line2D.Double(pts[0],cpt), cpoint);
		        		 Point2D pp2=GeomUtil.projectPointOntoLine(new Line2D.Double(pts[1],cpt), cpoint);
		        		 Point2D pp3=GeomUtil.projectPointOntoLine(new Line2D.Double(pts[2],cpt), cpoint);
		        		 Point2D np=Stream.of(pp1,pp2,pp3)
		        		       .map(p->Tuple.of(p,cpt.distance(p)).withVComparator())
		        		       .min(Comparator.naturalOrder())
		        		       .map(t->t.k())
		        		       .orElse(cpt);
		        		 n.setPoint(np);
		        		 
		        	  }
		          });
		return this;
		
		
	}

	public List<Edge> getEdgesWithCenterWithin(Point2D c, double d) {
		return this.edges.stream()
		          .filter(e->GeomUtil.findCenterOfShape(e.getLine()).distance(c)<d)
		          .collect(Collectors.toList());
		
	}

	public int getSumCharge() {
		return this.nodes.stream().mapToInt(n->n.getCharge()).sum();
	}

	
	
	
}