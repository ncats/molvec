package gov.nih.ncats.molvec.internal.util;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.util.*;
import java.util.function.BinaryOperator;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import gov.nih.ncats.molvec.MolvecOptions;
import gov.nih.ncats.molvec.internal.image.Bitmap;
import gov.nih.ncats.molvec.internal.algo.Tuple;
import gov.nih.ncats.molvec.internal.algo.Tuple.KEqualityTuple;
import gov.nih.ncats.molvec.internal.util.GeomUtil.ShapeWrapper;

public class ConnectionTable{
	private List<Node> nodes = new ArrayList<Node>();
	private List<Edge> edges = new ArrayList<Edge>();

	
	private static final int AROMATIC_ORDER=0xDE10CA1;
	
	private CachedSupplier<Map<Integer,List<Edge>>> _bondMap = CachedSupplier.of(()->_getEdgeMap());
	private CachedSupplier<Map<Node,Integer>> _nodeMap = CachedSupplier.of(()->_getNodeMap());
	private CachedSupplier<List<Ring>> _ring = CachedSupplier.of(()->_getRingMap());
	
	
	public static class Ring{
		private List<Node> nodes;
		private List<Edge> edges;
		
		private ConnectionTable parent;
		
		public Ring(List<Node> nodes, List<Edge> edges, ConnectionTable par){
			this.parent=par;
			this.setNodes(nodes);
			this.setEdges(edges);
		}
		
		public static Ring of(List<Node> ring, ConnectionTable ct){
			Set<Node> contain = ring.stream().collect(Collectors.toSet());
			
			List<Edge> edges=ring.stream()
			    .flatMap(n->n.getEdges().stream())
			    .distinct()
			    .filter(e->contain.contains(e.getRealNode1()) && contain.contains(e.getRealNode2()))
			    .collect(Collectors.toList());
			return new Ring(ring,edges,ct);
			
		}
		
		public Shape getConvexHull(){
			return this.getNodes().stream().map(n->n.getPoint()).collect(GeomUtil.convexHull());
		}

		public int size() {
			return getNodes().size();
		}
		
		public boolean isConjugated(){
			long c=this.nodes.stream()
					  .filter(n->n.getEdges().stream().filter(e->e.getOrder()==2).findAny().isPresent())
					  .count();
			return c==this.size();
		}

		public List<Edge> getEdges() {
			return edges;
		}

		public void setEdges(List<Edge> edges) {
			this.edges = edges;
		}

		public List<Node> getNodes() {
			return nodes;
		}

		public void setNodes(List<Node> nodes) {
			this.nodes = nodes;
		}
		
		
	}
	
	public List<Ring> getSmallestSetOfSmallestRings(){
		return null;
	}
	
	
	private static int LARGEST_RING=8;
	
	public List<Ring> getRings(){
		return this._ring.get();
	}
	
	
	private void consumePathsUntilRing(Stack<Node> soFar, Set<Edge> used,Set<Node> ignoreNodes, Consumer<Stack<Node>> found, int MAX_DEPTH) throws InterruptedException{
		Node p=soFar.peek();
		if(soFar.size()>MAX_DEPTH)return;
		
		if(Thread.interrupted()){
			throw new InterruptedException();
			
		}
		
		for(Tuple<Node,Edge> en : p.getNeighborNodes()){
			if(used.contains(en.v()) || ignoreNodes.contains(en.k())){
				continue;
			}
			used.add(en.v());
			
			if(soFar.get(0).equals(en.k())){
				soFar.push(en.k());
				found.accept(soFar);
				soFar.pop();
			}
			soFar.push(en.k());		
			
			consumePathsUntilRing(soFar,used,ignoreNodes,found, MAX_DEPTH);
			
			used.remove(en.v());
			soFar.pop();
		}
	}
	
	
	private List<Ring> _getRingMap(){
		int MAX_RING_TOTAL = 10;
		int MAX_RING_AFTER_INITIAL = 6;
		
		int maxNumberRings = 100;
		
		List<Ring> rings= getDisconnectedNodeSets().stream()
				.flatMap(nl1->{
					try{
						
						long minSSSR= 
								     nl1.stream()
								        .flatMap(n->n.getEdges().stream())
								        .distinct()
								        .count() - nl1.size() +1;

						if(minSSSR>maxNumberRings)return Stream.empty();
												
						Map<Node,List<List<Node>>> nrings = new HashMap<>();
						
						Set<Node> terms = new HashSet<Node>();
						List<Node> check = new ArrayList<Node>(nl1);
						int tbefore=0;
						
						
						
						while(true){
							for(Node nn: check){
								long keepCount = nn.getNeighborNodes().stream().filter(t->!terms.contains(t.k())).count();
								if(keepCount==1){
									terms.add(nn);
								}else{
									if(nn.getEdgeCount()>7){
										terms.add(nn);
									}
								}
								
							}
							check.removeAll(terms);
							
							if(tbefore==terms.size()){
								break;
							}
							tbefore=terms.size();
						}
						
						
						
						for(Node nn: check){
							
							
							Stack<Node> st=new Stack<Node>();
							Set<Edge> nadda=new HashSet<>();
							
							boolean[] foundRing = new boolean[]{false};
							
							st.push(nn);
							
							consumePathsUntilRing(st,nadda,terms,(nst)->{
								Node term = nst.peek();
								List<Node> mlist=new ArrayList<Node>();
								boolean started = false;
								for(Node n1:nst){
									if(started){
										mlist.add(n1);
									}else{
										if(n1==term){
											started=true;
										}
									}
								}
								for(Node n:mlist){
									nrings.computeIfAbsent(n, k->{
										return new ArrayList<List<Node>>();
									}).add(mlist);
								}
								foundRing[0]=true;
							},MAX_RING_AFTER_INITIAL);
							
							if(!foundRing[0]){
								st.clear();
								nadda.clear();
								st.push(nn);
								consumePathsUntilRing(st,nadda,terms,(nst)->{
									Node term = nst.peek();
									List<Node> mlist=new ArrayList<Node>();
									boolean started = false;
									for(Node n1:nst){
										if(started){
											mlist.add(n1);
										}else{
											if(n1==term){
												started=true;
											}
										}
									}
									for(Node n:mlist){
										nrings.computeIfAbsent(n, k->{
											return new ArrayList<List<Node>>();
										}).add(mlist);
									}
									foundRing[0]=true;
								},MAX_RING_TOTAL);
							}			
						
					}
					
					return nrings.entrySet()
					      .stream()
					      .map(Tuple::of)
					      .map(Tuple.vmap(nl->{
					    	  List<Node> best = nl.stream()
					    	    .map(nn->Tuple.of(nn,nn.size()).withVComparator())
					    	    .min(Comparator.naturalOrder())
					    	    .map(t->t.k())
					    	    .orElse(null);
					    	  return nl.stream().filter(nn->nn.size()==best.size()).collect(Collectors.toList());
					      }))
					      .flatMap(t->t.v().stream())
					      .map(t->Tuple.of(t,t.stream().mapToInt(nn->nn.getIndex()).sorted().mapToObj(i->i+"").collect(Collectors.joining())))
					      .map(t->t.swap())
					      .map(t->t.withKEquality())
					      .distinct()
					      .map(t->t.swap())
					      .map(t->t.withVComparator())
					      .sorted(Comparator.reverseOrder())
					      .map(t->t.swap())
					      .map(t->t.v())
					      .map(t->Tuple.of(t,t.size()).withVComparator())
					      .sorted()
					      .map(t->t.k())
					      .map(n->Ring.of(n, this));
					}catch(InterruptedException te){
						te.printStackTrace();
						throw new RuntimeException(te);
					}
				})
				.collect(Collectors.toList());
		

		return rings;
		
		
		
//		
//		Set<Node> notRings = new HashSet<Node>();
//		List<Node> possible=this.getNodes();
//		
//		boolean[] found=new boolean[]{true};
//		
//		while(found[0]){
//			found[0]=false;
//			possible=possible
//			    .stream()
//			    .filter(n->{
//			    	long c=n.getNeighborNodes()
//			    	 .stream()
//			    	 .map(n1->n1.k())
//			    	 .filter(n1->!notRings.contains(n1))
//			    	 .count();
//			    	
//			    	if(c<2){
//			    		notRings.add(n);
//			    		found[0]=true;
//			    		return false;
//			    	}
//			    	return true;
//			    })
//			    .collect(Collectors.toList());
//		}
//		
//		
//		for(Node n:possible){
//		
//			Set<Node> neighborsStart =n.getNeighborNodes().stream().map(t->t.k()).collect(Collectors.toSet());
//			Set<Node> neighborsNext =new HashSet<>();
//			List<Node> dontInclude =new ArrayList<>();
//			dontInclude.addAll(notRings);
//			dontInclude.add(n);
//			
//			for(int i=0;i<LARGEST_RING;i++){
//				for(Node nn: neighborsStart){
//					nn.getNeighborNodes().stream()
//										 .map(t->t.k())
//					                     .filter(n2->!dontInclude.contains(n2))
//					                     .forEach(n2->{
//					                    	 neighborsNext.add(n2);
//					                     });
//				}
//				
//				if(neighborsNext.contains(n)){
//					//ring
//				}
//				dontInclude.clear();
//				dontInclude.addAll(neighborsStart);
//				neighborsStart= new HashSet<>(neighborsNext);
//				neighborsNext.clear();
//			}
//		}
	}
	
	
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
	public List<Tuple<GeomUtil.LineWrapper,List<Edge>>> getEdgesWhichMightBeWiggleLines(){
		double maxRatio = 0.679;

		List<List<Edge>> elist= this.getEdges()
		    .stream()
		    .map(e->Tuple.of(e,e.getLine()))
		    .map(Tuple.vmap(l-> GeomUtil.LineWrapper.of(l)))
		    .collect(GeomUtil.groupThings(t->{
		    	GeomUtil.LineWrapper lw1=t.k().v();
		    	GeomUtil.LineWrapper lw2=t.v().v();
		    	
		    	double rat = lw1.length()*lw2.recipLength();
		    	if(rat>maxRatio && rat<=(1/maxRatio)){
		    		double lim = (lw1.length()+lw2.length())/2;
		    		
		    		if(lw1.centerPoint().distance(lw2.centerPoint()) <=lim){
		    			return true;
		    		}
		    	}
		    	
		    	
		    	return false;
		    }))
		    .stream()
		    .map(l->l.stream().map(t->t.k()).collect(Collectors.toList()))
		    .collect(Collectors.toList());
		
		double maxLent=elist.stream()
		     .filter(ll->ll.size()>1)
		     .mapToDouble(ll->ll.stream().mapToDouble(l->l.getEdgeLength()).average().orElse(0))
		     .max()
		     .orElse(1);
		
		if(maxLent>this.getAverageBondLength()*1.8){
			maxLent=this.getAverageBondLength()*1.8;
		}
		double maxLen = maxLent;
		
		
		List<List<Edge>> elist2= this.getEdges()
			    .stream()
			    .filter(e->e.getEdgeLength()<maxLen*maxRatio)
			    .map(e->Tuple.of(e,e.getLine()))
			    .map(Tuple.vmap(l-> GeomUtil.LineWrapper.of(l)))
			    .collect(GeomUtil.groupThings(t->{
			    	GeomUtil.LineWrapper lw1=t.k().v();
			    	GeomUtil.LineWrapper lw2=t.v().v();
			    	double lim = (lw1.length()+lw2.length())/2;
		    		
		    		if(lw1.centerPoint().distance(lw2.centerPoint()) <=lim){
		    				return true;
		    		}
			    	
			    	
			    	return false;
			    }))
			    .stream()
			    .map(l->l.stream().map(t->t.k()).collect(Collectors.toList()))
			    .collect(Collectors.toList());
		
		//there are some issues here
		
		return elist2.stream()
		     .filter(ll->ll.size()>=5)
		     .filter(ll->ll.stream().allMatch(l1->l1.getOrder()==1))
		     .map(ll->Tuple.of(ll,ll.stream()
		    		                .flatMap(e->e.streamNodes())
		    		                .map(n->n.getPoint())
		    		                .collect(GeomUtil.convexHull())
		    		                ))
		     .map(Tuple.vmap(s->GeomUtil.findLongestSplittingLine(s)))
		     .filter(t->t.v().length()<=maxLen*(1/maxRatio))
		     .filter(t->t.v().length()>=maxLen*(maxRatio))
		     .filter(t->{
		    	 int[] cc=new int[]{0,0};
		    	 t.k().stream()
		    	  .flatMap(e->e.streamNodes()).distinct()
		    	  .map(n->n.getPoint())
		    	  .map(p->GeomUtil.projectPointOntoLineWithRejection(t.v().getLine(), p))
		    	  .map(t1->t1.v())
		    	  .forEach(d->{
		    		  if(d>0){
		    			  cc[0]++;
		    		  }else{
		    			  cc[1]++;
		    		  }
		    	  });
		    	 
		    	 if(cc[0]>0 && cc[1]>0){
		    		 double ratAbove=((double)cc[0])/(double)((cc[0]+cc[1]));
		    		 if(ratAbove>=0.3 && ratAbove<=0.7){
		    			 return true;
		    		 }
		    	 }
		    	 
		    	 return false;
		     })
		     .map(t->t.swap())
		     .collect(Collectors.toList())
		     ; 
		
		
		    
	}
	
	public List<List<Edge>> getEdgesWhichMightBeDottedLines(){
//		double longestEdge1=this.getEdges()
//		    .stream()
//		    .map(e->e.getEdgeLength())
//		    .max(Comparator.naturalOrder())
//		    .orElse(1.0);
		double longestEdge =this.getEdges()
								.stream()
								.filter(e->e.getNeighborEdges().size()>0)
								.mapToDouble(e->e.getEdgeLength())
								.average()
								.orElse(this.getAverageBondLength());
		
		
		
		return this.getEdges()
			    .stream()
			    .filter(e->e.getEdgeLength()<longestEdge*0.3)
			    .map(e->Tuple.of(e, GeomUtil.LineWrapper.of(e.getLine())))
			    .collect(GeomUtil.groupThings(t->{
			    	GeomUtil.LineWrapper l1=t.k().v();
			    	GeomUtil.LineWrapper l2=t.v().v();
			    	
			    	if(l1.absCosTheta(l2)<0.8){
			    		return false;
			    	}
			    	if(l1.centerPoint().distance(l2.centerPoint())> longestEdge*0.3){
			    		return false;
			    	}
			    	return true;
			    }))
			    .stream()
			    .map(l->l.stream().map(t->t.k()).collect(Collectors.toList()))
			    .filter(nl->nl.size()>1)
			    .collect(Collectors.toList());
		    
		    
	}
	
	
	public List<List<Node>> getDisconnectedNodeSets(){
		//note: very inefficient
		return GeomUtil.groupThings(this.nodes, t->t.k().connectsTo(t.v()))
		.stream()
		.collect(Collectors.toList());
	}
	public String toMol(){
		return DEFAULT_OPTIONS.computeResult(this).getMolfile().get();
	}

	private static final MolvecOptions DEFAULT_OPTIONS = new MolvecOptions();



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
		double maxDSq=maxDistance*maxDistance;
		GeomUtil.groupThings(nodes.stream().filter(incudeNode).collect(Collectors.toList()), (t)->{
								Node n1=t.k();
								Node n2=t.v();
								return n1.point.distanceSq(n2.point)<maxDSq;
							})
							.stream()
							.filter(nl->nl.size()>=2)
//							.collect(Collectors.toList())
							.forEach(ml->{
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
	
	public List<Node> getNodesInsideShape(ShapeWrapper s, double tol){
		List<Node> mnodes= new ArrayList<>();
		
		for(int i=nodes.size()-1;i>=0;i--){
			Point2D pn = nodes.get(i).point;
			if(s.distanceTo(pn)<tol || s.contains(pn)){
				mnodes.add(nodes.get(i));
			}
		}
		return mnodes;
	}
	
	public Tuple<Node,Double> getClosestNodeToShape(ShapeWrapper s){
		return nodes.stream()
		     .map(n->Tuple.of(n,s.distanceTo(n.point)).withVComparator())
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
	

	public ConnectionTable getNewConnectionTable(List<ShapeWrapper> likelyNodes,
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
	
	
	public ConnectionTable mergeAllNodesInside(ShapeWrapper s, double tol,Predicate<Node> allow, Function<List<Point2D>, Point2D> merger){
		
		List<Integer> toMerge = getAllNodeIndexesInsideShape(s,tol);
		toMerge=toMerge.stream()
		       .filter(i->allow.test(nodes.get(i)))
		        .collect(Collectors.toList());
			
		return this.mergeNodes(toMerge, merger);
	}
	
	private List<Integer> getAllNodeIndexesInsideShape(ShapeWrapper s,double tol){
		List<Integer> toMerge = new ArrayList<Integer>();
		for(int i=nodes.size()-1;i>=0;i--){
			Point2D pn = nodes.get(i).point;
			if(s.distanceTo(pn)<=tol){
				toMerge.add(i);
			}
		}
		return toMerge;
	}
	public List<Node> getAllNodesInsideShape(ShapeWrapper s,double tol){
		return getAllNodeIndexesInsideShape(s,tol).stream().map(i->this.nodes.get(i)).collect(Collectors.toList());
	}
	
	
	public List<Tuple<Edge,Tuple<Node,Node>>> getAllEdgesEntering(ShapeWrapper s, double tol){
		
		List<Node> toMerge = new ArrayList<Node>();
		for(int i=nodes.size()-1;i>=0;i--){
			Point2D pn = nodes.get(i).point;
			if(s.distanceTo(pn)<tol){
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
	
	public ConnectionTable mergeAllNodesInsideCenter(ShapeWrapper s, double tol){
		Point2D p = s.centerOfBounds();
		return mergeAllNodesInside(s,tol,n->true,(l)->p);
	}
	
	public ConnectionTable mergeNodesExtendingTo(Collection<ShapeWrapper> shapes,double maxAvgBondRatio, double maxTotalAvgBondRatio){
		double avg = this.getAverageBondLength();
		edges
//                .stream()
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
		    	 ShapeWrapper closest1=null;
		    	 ShapeWrapper closest2=null;
		    	 
		    	 for(ShapeWrapper s: shapes){
		    		 double d1=s.distanceTo(p1);
		    		 double d2=s.distanceTo(p2);
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
		    		newPoint1=closest1.centerOfBounds();
		    		newPoint2=closest2.centerOfBounds();
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
			    		 newPoint2=closest2.centerOfBounds();
			    		 newPoint1=e.getPoint1();
			    		 //newLine=new Line2D.Double(p1,newPoint);
			    	 }else{
			    		 newPoint1=closest1.centerOfBounds();
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
		
//		List<LinkedHashSet<Integer>> mergeNodes=
                GeomUtil.groupMultipleBonds(edgeMap.keySet().stream().map(l-> GeomUtil.LineWrapper.of(l)).collect(Collectors.toList()),5*Math.PI/180, 2, .8, 0)
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
//		.collect(Collectors.toList());
//
//		mergeNodes
                .forEach(ls->{
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
		_ring.resetCache();
		_averageBondLength.resetCache();
		
	}
	
	public double getAverageBondLength(){
		return _averageBondLength.get();
	}
	
	public double getAverageBondLengthSquared(){
		double bl = getAverageBondLength();
		return bl*bl;
	}
	public double getLargestBondLength(){
		return this.getEdges().stream().mapToDouble(e->e.getEdgeLength()).max().orElse(1);
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
		private int group=0;
		private String alias = null;
		private boolean tooClose = false;
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
			List<Ring> rings = ConnectionTable.this.getRings();
			
			return rings.stream()
					.filter(r->r.size()<=maxRing)
			     .filter(r->r.getNodes().contains(this))
			     .findAny()
			     .isPresent();
			
			
			
		}
		
		public List<Ring> getAllRings(){
			List<Ring> rings = ConnectionTable.this.getRings();
			
			return rings.stream()
					.filter(r->r.getNodes().contains(this))
					.collect(Collectors.toList());
		}

		public Optional<Ring> getSmallestRing(){
			List<Ring> rings = ConnectionTable.this.getRings();
			
			return rings.stream()
					.filter(r->r.getNodes().contains(this))
				    .map(r->Tuple.of(r,r.size()).withVComparator())
				    .min(Comparator.naturalOrder())
				    .map(t->t.k());
		}
		
		public OptionalInt getSmallestRingSize(){
			List<Ring> rings = ConnectionTable.this.getRings();
			
			return rings.stream()
					.filter(r->r.getNodes().contains(this))
					.mapToInt(r->r.size())
					.min();
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

		public boolean hasInventedBond() {
			return this.getEdges().stream().anyMatch(b->b.isInventedBond());
		}

		public Node markGroup(int gnum) {
			this.group=gnum;
			return this;
		}
		
		public int getGroup() {
			return this.group;
		}
		
		public String getAlias(){
			return this.alias;
		}
		public Node setAlias(String ali){
			this.alias=ali;
			return this;
		}

		public Node markTooClose(boolean b) {
			this.tooClose=true;
			return this;
		}
		public boolean isTooClose(){
			return this.tooClose;
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
		public int getNode1Offset(){
			return n1;
		}
		public int getNode2Offset(){
			return n2;
		}
		public Stream<Node> streamNodes(){
			return Stream.of(this.getRealNode1(),this.getRealNode2());
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

		public boolean isRingEdge() {
			return !getAllRings().isEmpty();

		}

		public List<Ring> getAllRings() {
		
			List<Ring> rings = ConnectionTable.this.getRings();
			
			return rings.stream()
					.filter(r->r.getEdges().contains(this))
					.collect(Collectors.toList());
			
		}

		public Point2D getCenterPoint() {
			return GeomUtil.findCenterOfShape(this.getLine());
		}

		public double getEdgeLengthSquared() {
			return ConnectionTable.this.nodes.get(n1).point.distanceSq(ConnectionTable.this.nodes.get(n2).point);
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
		Stroke old = g2.getStroke();
		
		Stroke dashed = new BasicStroke(3, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0, new float[]{3}, 0);
		Stroke wedge = new BasicStroke(5f);
        
		
		nodes.stream().map(n->n.point)
		.forEach(p->{
			g2.draw(new Ellipse2D.Double((p.getX()-2D), (p.getY()-2D), 4D, 4D));
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

	public ConnectionTable makeMissingBondsToNeighbors(Bitmap bm, double d, double tol, Collection<ShapeWrapper> OCRSet, double ocrTol,Consumer<Tuple<Double,Edge>> econs) {
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
						Tuple<ShapeWrapper,Double> t1=GeomUtil.findClosestShapeWTo(OCRSet, n1.point).orElse(null);
						Tuple<ShapeWrapper,Double> t2=GeomUtil.findClosestShapeWTo(OCRSet, n2.point).orElse(null);
						
						Point2D pt1=n1.point;
						Point2D pt2=n2.point;
						
						List<Line2D> lines = new ArrayList<Line2D>();
						lines.add(new Line2D.Double(pt1,pt2));
						
						
						if(t1!=null){
							lines=lines.stream()
									   .map(l->t1.k().getLinesNotInside(l))
									   .flatMap(ll->ll.stream())
									   .collect(Collectors.toList());
						}
						if(t2!=null){
							lines=lines.stream()
									   .map(l->t2.k().getLinesNotInside(l))
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
		this.edges
//		    .stream()
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
	
	public List<Node> getNodesNotInShapes(Collection<ShapeWrapper> shapes, double tol){
		if(shapes.isEmpty())return nodes;
		return nodes.stream()
		     .map(n->Tuple.of(n,GeomUtil.getClosestShapeToSW(shapes, n.point)))
		     .filter(t->t.v().distanceTo(t.k().point)>tol)
		     .map(t->t.k())
		     .collect(Collectors.toList());
	}

	public ConnectionTable makeMissingNodesForShapes(Collection<ShapeWrapper> likelyOCR, double mAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL,
			double mIN_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL) {
		double avg=this.getAverageBondLength();
		List<ShapeWrapper> addShapes=likelyOCR.stream() 
		         .map(oc->Tuple.of(getClosestNodeToShape(oc),oc))
		        .filter(t->t.k()!=null && t.k().v() !=null)
				.filter(t->t.k().v()>avg*mIN_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL)
				 .filter(t->t.k().v()<avg*mAX_BOND_TO_AVG_BOND_RATIO_FOR_NOVEL)
				 .map(t->t.v())
				 .collect(Collectors.toList());
		
		for(ShapeWrapper s:addShapes){
			this.addNode(s.centerOfBounds());
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
		          .filter(n->!n.isInvented())
		          .filter(n->!n.hasInventedBond())
		          .filter(n->n.getEdges().stream().filter(e->e.getOrder()>1).map(e->e.getOtherNode(n).getSymbol()).anyMatch(s->!"C".equals(s)))
		          .filter(n->n.getSmallestRingSize().orElse(999)>4)
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


	public Node getClosestNodeToPoint(Point2D point) {
		return this.getNodes().stream()
				.map(n->Tuple.of(n,n.getPoint().distanceSq(point)).withVComparator())
				.min(Comparator.naturalOrder())
				.map(t->t.k())
				.orElse(null);
				
	}

	
	
	
}