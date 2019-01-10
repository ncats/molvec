package tripod.molvec.algo;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BinaryOperator;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import gov.nih.ncats.chemkit.api.Atom;
import gov.nih.ncats.chemkit.api.Bond.BondType;
import gov.nih.ncats.chemkit.api.Chemical;
import gov.nih.ncats.chemkit.api.ChemicalBuilder;
import tripod.molvec.util.GeomUtil;

public class LineUtil {

	public static List<Line2D> asLines(Collection<Path2D> segments){
		List<Line2D> lines= new ArrayList<Line2D>();
        for(Path2D p2:segments){
            PathIterator pi=p2.getPathIterator(null);
            double[] prevPt=null;
            while(!pi.isDone()){
        		
                double[] coord= new double[2];
                pi.currentSegment(coord);
                if(prevPt!=null){
                    Line2D line = new Line2D.Double
                        (coord[0], coord[1], prevPt[0], prevPt[1]);
//                    double lineLength=Math.sqrt
//                        ((coord[0]-prevPt[0])*(coord[0]-prevPt[0])
//                         +(coord[1]-prevPt[1])*(coord[1]-prevPt[1]));
                    lines.add(line);
                    //System.out.println(lineLength);
                }
                prevPt=coord;
                pi.next();
            }
        }
        return lines;
	}
	
	public static Path2D fromLine(Line2D line){
		GeneralPath gp = new GeneralPath();
		gp.moveTo(line.getX1(),line.getY1());
		gp.lineTo(line.getX2(),line.getY2());
		return gp;
	}
	
	public static List<Path2D> fromLines(List<Line2D> lines){
		return lines.stream().map(l->fromLine(l)).collect(Collectors.toList());
	}
	
	
	
	public static double length(Line2D l){
		return l.getP1().distance(l.getP2());
	}
	
	/**
	 * <p>Currently, this method takes in a set of lines and attempts to group the lines based on being</p>
	 * <ol>
	 * <li>Sufficiently parallel</li>
	 * <li>Sufficiently close together</li>
	 * <li>Sufficiently "overlapping"</li> 
	 * </ol>
	 * 
	 * 
	 * @param lines
	 * @param maxDeltaTheta
	 * @param maxDeltaOffset
	 * @param minIntersectionRatio
	 * @return
	 */
	public static List<List<Line2D>> groupMultipleBonds(List<Line2D> lines, double maxDeltaTheta, double maxDeltaOffset, double minIntersectionRatio){
		double[] recipLengths=new double[lines.size()];
		Map<Integer,Integer> groups = new HashMap<>();
		for(int i=0;i<lines.size();i++){
			Line2D line=lines.get(i);
			recipLengths[i]=1.0/line.getP1().distance(line.getP2());
			groups.put(i, i);
		}
		double cosThetaCutoff=Math.cos(maxDeltaTheta);
		
		
		//So, this will be an n^2 operation
		for(int i=0;i<lines.size();i++){
			Line2D line1=lines.get(i);
			double[] vec1=asVector(line1);
			for(int j=i+1;j<lines.size();j++){
				Line2D line2=lines.get(j);
				double[] vec2=asVector(line2);
				double dot=vec1[0]*vec2[0]+vec1[1]*vec2[1];
				double cosTheta=Math.abs(dot*recipLengths[i]*recipLengths[j]);
				//Sufficiently parallel
				if(cosTheta>cosThetaCutoff){
					double[] vec2Off = new double[2];
					vec2Off[0] = (line2.getX1()+line2.getX2())/2 -line1.getX1();
					vec2Off[1] = (line2.getY1()+line2.getY2())/2 -line1.getY1();
					double antiDot=(vec2Off[0]*vec1[1]-vec2Off[1]*vec1[0]);
					double rej=Math.abs(antiDot*recipLengths[i]);
					//Sufficiently close to line
					if(rej<maxDeltaOffset){
						double len1=1/recipLengths[i];
						double proj1=0;
						double proj2=0;
						vec2Off[0] = line2.getX1() -line1.getX1();
						vec2Off[1] = line2.getY1() -line1.getY1();
						proj1=(vec2Off[0]*vec1[0]+vec2Off[1]*vec1[1])*recipLengths[i];
						vec2Off[0] = line2.getX2() -line1.getX1();
						vec2Off[1] = line2.getY2() -line1.getY1();
						proj2=(vec2Off[0]*vec1[0]+vec2Off[1]*vec1[1])*recipLengths[i];
						
						if(proj1<0)proj1=0;
						if(proj1>len1){
							proj1=len1;
						}
						if(proj2<0)proj2=0;
						if(proj2>len1){
							proj2=len1;
						}
						double intersectLength=Math.abs(proj1-proj2);
						double cutoffRat=Math.max(recipLengths[i], recipLengths[j]);
						if(intersectLength*cutoffRat>minIntersectionRatio){
							int g1=groups.get(i);
							int g2=groups.get(j);
							groups.put(i, Math.min(g1, g2));
							groups.put(j, Math.min(g1, g2));
						}
					}
				}
			}
		}
		
		return groups.entrySet().stream()
		      .map(Tuple::of)
		      .map(Tuple.kmap(i->lines.get(i)))
		      .map(t->t.swap())
		      .collect(Tuple.toGroupedMap())
		      .entrySet()
		      .stream()
		      .map(Tuple::of)
		      .map(t->t.v())
		      .collect(Collectors.toList());
		
		
	}
	
	/**
	 * This method, like {@link #groupMultipleBonds(List, double, double, double)}, attempts to group bonds, but also 
	 * goes one step further to select the longest line from each group while also returning the count of lines
	 * found in the group. This is an estimate of the bond order.
	 * 
	 * @param lines
	 * @param maxDeltaTheta
	 * @param maxDeltaOffset
	 * @param minIntersectionRatio
	 * @return
	 */
	public static List<Tuple<Line2D, Integer>> reduceMultiBonds(List<Line2D> lines, double maxDeltaTheta, double maxDeltaOffset, double minIntersectionRatio){
		
		return groupMultipleBonds(lines,maxDeltaTheta,maxDeltaOffset,minIntersectionRatio)
				.stream()
				.map(l->Tuple.of(l,l.size()))
				.map(Tuple.kmap(l->l.stream()
								   .map(s->Tuple.of(s,-length(s)).withVComparator())
								   .sorted()
								   .findFirst()
						           .get()
						           .k()
						           ))
				.collect(Collectors.toList());
		
		
		
	}
	
	public static double[] asVector(Line2D line){
		double dx=line.getX2()-line.getX1();
		double dy=line.getY2()-line.getY1();
		return new double[]{dx,dy};
		
	}

	
	
	public static class ConnectionTable{
		List<Node> nodes = new ArrayList<Node>();
		List<Edge> edges = new ArrayList<Edge>();
		
		
		public ConnectionTable addNode(Point2D p){
			nodes.add(new Node(p,"C"));
			return this;
		}
		
		public ConnectionTable setNodeToSymbol(Shape s, String sym){
			nodes.stream()
			.filter(n->s.contains(n.point.getX(),n.point.getY()))
			.forEach(n->{
				n.symbol=sym;
			});
			return this;
			
		}
		
		
		public Chemical toChemical(){
			ChemicalBuilder cb = new ChemicalBuilder();
			Atom[] atoms = new Atom[nodes.size()];
			
			for(int i=0;i<nodes.size();i++){
				Node n = nodes.get(i);
				atoms[i]=cb.addAtom(n.symbol,n.point.getX(),-n.point.getY());
			}
			for(Edge e : edges){
				if(e.order==1){
					cb.addBond(atoms[e.n1],atoms[e.n2],BondType.SINGLE);
				}else if(e.order==2){
					cb.addBond(atoms[e.n1],atoms[e.n2],BondType.DOUBLE);
				}else if(e.order==3){
					cb.addBond(atoms[e.n1],atoms[e.n2],BondType.TRIPLE);
				}else{
					cb.addBond(atoms[e.n1],atoms[e.n2],BondType.SINGLE);
				}
			}
			return cb.build();
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
			
			
			return this;
		}
		
		public Map<Integer,Integer> mergeNodesGetTransform(List<Integer> nlist, Function<List<Point2D>, Point2D> op){
			Point2D p = op.apply(nlist.stream().map(i->nodes.get(i).point).collect(Collectors.toList()));
			
			nlist=nlist.stream().sorted().collect(Collectors.toList());
			
			Map<Integer,Integer> oldToNewMap=new HashMap<>();
			for(int i=0;i<this.nodes.size();i++){
				oldToNewMap.put(i, i);
			}
			
			System.out.println("Start merge");
			if(nlist.size()==1){
				nodes.get(nlist.get(0)).point=p;
			}else{
				for(int i=nlist.size()-1;i>=1;i--){
					
					int ni1=nlist.get(i);
					int ni2=nlist.get(i-1);
					System.out.println("Merging:" + ni1 + " and " + ni2);
					Map<Integer,Integer> oldNewNew=mergeNodesGetTransform(ni1,ni2,(p1,p2)->p);
					oldToNewMap=oldToNewMap.entrySet().stream()
					   .map(Tuple::of)
					   .peek(t->System.out.println("PFrom:"+t.k() + " to "+t.v()))
					   .map(Tuple.vmap(oi->oldNewNew.get(oi)))
					   .peek(t->System.out.println("From:"+t.k() + " to "+t.v()))
					   .collect(Tuple.toMap());
					
				}
			}
			System.out.println("End merge");
			return oldToNewMap;
		}
		
		public ConnectionTable cleanMeaninglessEdges(){
			this.edges=edges.stream().filter(e->e.n1!=e.n2).collect(Collectors.toList());
			return this;
		}
		
		public ConnectionTable cleanDuplicateEdges(BinaryOperator<Edge> combiner){
			edges=edges.stream().map(e->e.standardize())
			              .map(e->Tuple.of(e,e.n1+"_" + e.n2))
			              .collect(Collectors.toMap(t->t.v(),t-> t.k(), combiner))
			              .values()
			              .stream()
			              .collect(Collectors.toList());
			return this;
		}
		
		public ConnectionTable mergeNodesCloserThan(double maxDistance){
			boolean mergedOne = true;
			
			while(mergedOne){
				mergedOne=false;
				for(int i=0;i<nodes.size();i++){
					Point2D pnti=nodes.get(i).point;
					for(int j=i+1;j<nodes.size();j++){
						Point2D pntj = nodes.get(j).point;
						if(pnti.distance(pntj)<maxDistance){
							mergeNodesAverage(i,j);
							System.out.println("Merged");
							mergedOne=true;
							break;
						}
					}
					if(mergedOne)break;
				}
			}
			return this;
		}
		
		public ConnectionTable addEdge(int n1, int n2, int o){
			this.edges.add(new Edge(n1,n2,o));
			return this;
		}
		
		public ConnectionTable mergeAllNodesInside(Shape s, double tol){
			Rectangle2D r=s.getBounds2D();
			Point2D p = new Point2D.Double(r.getCenterX(),r.getCenterY());
			
			List<Integer> toMerge = new ArrayList<Integer>();
			for(int i=nodes.size()-1;i>=0;i--){
				Point2D pn = nodes.get(i).point;
				if(GeomUtil.distanceTo(s,pn)<tol){
					toMerge.add(i);
				}
			}
			return this.mergeNodes(toMerge, (l)->p);
		}
		
		public ConnectionTable mergeAllNodesOnParLines(){
			Map<Line2D,Edge> edgeMap = this.edges.stream().collect(Collectors.toMap(e->e.getLine(),e->e ));
			
			Map<Integer,Integer> oldToNewMap = IntStream.range(0,this.nodes.size())
			         .mapToObj(i->i)
			         .collect(Collectors.toMap(i->i, i->i));
			
			List<LinkedHashSet<Integer>> mergeNodes=groupMultipleBonds(edgeMap.keySet().stream().collect(Collectors.toList()),5*Math.PI/180, 2, .8)
			.stream()
			.filter(l->l.size()>1)
			.map(l->{
				Line2D keep=l.stream()
				 .map(l1->Tuple.of(-LineUtil.length(l1),l1).withKComparator())
				 .sorted()
				 .findFirst()
				 .get()
				 .v();
				Edge keepEdge=edgeMap.get(keep);
				LinkedHashSet<Integer> mergeNode1 = new LinkedHashSet<Integer>();
				LinkedHashSet<Integer> mergeNode2 = new LinkedHashSet<Integer>();
				mergeNode1.add(keepEdge.n1);
				mergeNode2.add(keepEdge.n2);
				
				Point2D node1Point = keepEdge.getNode1();
				Point2D node2Point = keepEdge.getNode2();
				
				l.stream()
				 .map(l1->edgeMap.get(l1))
				 .filter(e->e!=keepEdge)
				 .forEach(me->{
					 double distance1to1=me.getNode1().distance(node1Point);
					 double distance1to2=me.getNode1().distance(node2Point);
					 double distance2to1=me.getNode2().distance(node1Point);
					 double distance2to2=me.getNode2().distance(node2Point);
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
				System.out.println(toMerge.toString());
				Map<Integer,Integer> newTrans=mergeNodesGetTransform(toMerge,(pts)->keeper);
				for(int i : oldToNewMap.keySet()){
					int oldMap=oldToNewMap.get(i);
					System.out.println(i+"=>" + oldMap);
					//System.out.println(newTrans.toString());
					int newMap=newTrans.get(oldMap);
					oldToNewMap.put(i, newMap);
				}
			});
			return this;
		}
		
		public ConnectionTable createNodesOnIntersectingLines(){
			
			boolean splitOne =true;
			while(splitOne){
				splitOne=false;
			
				for(int i=0;i<edges.size();i++){
					Edge e1=edges.get(i);
					Set<Integer> i1set=e1.getNodeSet();
					for(int j=0;j<edges.size();j++){
						Edge e2=edges.get(j);
						if(i1set.contains(e2.n1) || i1set.contains(e2.n2)){
							continue;
						}
						if(e1.getLine().intersectsLine(e2.getLine())){
							Point2D np=GeomUtil.intersection(e1.getLine(),e2.getLine());
							int nodeNew = this.nodes.size();
							this.addNode(np);
							Edge nedge1 = new Edge(e1.n1, nodeNew, 1);
							Edge nedge2 = new Edge(e1.n2, nodeNew, 1);
							Edge nedge3 = new Edge(e2.n1, nodeNew, 1);
							Edge nedge4 = new Edge(e2.n2, nodeNew, 1);
							edges.set(i, nedge1);
							edges.set(j, nedge2);
							edges.add(nedge3);
							edges.add(nedge4);
							splitOne=true;
							break;
						}
					}
					if(splitOne)break;
				}
			}
			
			return this;
		}
		
		public double getAverageBondLength(){
			return edges.stream().mapToDouble(e->e.getBondDistance()).average().getAsDouble();
		}
		public class Node{
			Point2D point;
			String symbol="C";
			public Node(Point2D p, String s){
				this.point=p;
				this.symbol=s;				
			}
			
			
		}
		public class Edge{
			int n1;
			int n2;
			int order;
			public Edge(int n1, int n2, int o){
				this.n1=n1;
				this.n2=n2;
				this.order=o;
			}
			public double getBondDistance(){
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
			public Point2D getNode1(){
				return ConnectionTable.this.nodes.get(n1).point;
			}
			public Point2D getNode2(){
				return ConnectionTable.this.nodes.get(n2).point;
			}
			public int getOrder() {
				
				return this.order;
			}
			
		}
		public void draw(Graphics2D g2) {
			int sx=1;
			
			nodes.stream().map(n->n.point)
			.forEach(p->{
				g2.draw(new Ellipse2D.Double((p.getX()-2f/sx), (p.getY()-2f/sx), 4f/sx, 4f/sx));
			});
			
			edges.forEach(l->{
				if(l.order==1){
					g2.setPaint(Color.BLACK);
				}else if(l.order==2){
					g2.setPaint(Color.RED);
				}else if(l.order==3){
					g2.setPaint(Color.GREEN);
				}
				g2.draw(l.getLine());
			});
			
		}
		
		
	}
	
	
	public static Line2D longestLineFromOneVertexToPoint(Line2D line, Point2D pnt){
		double d1= pnt.distance(line.getP1());
		double d2= pnt.distance(line.getP2());
		if(d1>d2){
			return new Line2D.Double(line.getP1(), pnt);
		}else{
			return new Line2D.Double(line.getP2(), pnt);
		}
	}
	public static ConnectionTable getConnectionTable(List<Tuple<Line2D,Integer>> linest, List<Shape> likelyNodes, double maxDistanceRatioNonLikely,double maxDistanceRatioLikely, Consumer<List<Tuple<Line2D,Integer>>> whenDone){
		//This will find the likely intersection points between lines.
				//If 2 lines would intersect 
				List<Tuple<Line2D,Integer>> lines = new ArrayList<>(linest);
				
				for(int i=0;i<lines.size();i++){
					Tuple<Line2D, Integer> lineOrder1=lines.get(i);
					double distance1 = LineUtil.length(lineOrder1.k());
					for(int j=i+1;j<lines.size();j++){
						
						Tuple<Line2D, Integer> lineOrder2=lines.get(j);
						double distance2 = LineUtil.length(lineOrder2.k());
						double totalDistance = distance1+distance2;
						
						Point2D intersect = GeomUtil.intersection(lineOrder1.k(),lineOrder2.k());
						if(intersect==null)continue;
						double ndistance1 = Math.max(intersect.distance(lineOrder1.k().getP1()), intersect.distance(lineOrder1.k().getP2()));
						double ndistance2 = Math.max(intersect.distance(lineOrder2.k().getP1()), intersect.distance(lineOrder2.k().getP2()));
						double totalDistanceAfter = ndistance1+ndistance2;
						
						double ratio = Math.max(totalDistanceAfter,totalDistance)/Math.min(totalDistanceAfter, totalDistance);
						
						boolean merge = false;
						
						if(ratio<maxDistanceRatioLikely){
							if(ratio<maxDistanceRatioNonLikely){
								System.out.println(totalDistanceAfter +  "vs" + totalDistance);
								merge=true;
							}else{
								boolean inLikelyNode=likelyNodes.stream().filter(s->s.contains(intersect)).findAny().isPresent();
								if(inLikelyNode){
									merge=true;
								}
							}
						}
						
						if(merge){
							Line2D newLine1 = longestLineFromOneVertexToPoint(lineOrder1.k(),intersect);
							Line2D newLine2 = longestLineFromOneVertexToPoint(lineOrder2.k(),intersect);
							Tuple<Line2D,Integer> norder1 = Tuple.of(newLine1,lineOrder1.v());
							Tuple<Line2D,Integer> norder2 = Tuple.of(newLine2,lineOrder2.v());
							lines.set(i, norder1);
							lines.set(j, norder2);
							lineOrder1 = norder1;
						}				
					}
				}
					
				ConnectionTable ct=new ConnectionTable();
				
				for(int i=0;i<lines.size();i++){
					Tuple<Line2D, Integer> lineOrder1=lines.get(i);
					ct.addNode(lineOrder1.k().getP1());
					ct.addNode(lineOrder1.k().getP2());
					ct.addEdge(i*2, i*2+1, lineOrder1.v());
				}
				whenDone.accept(lines);
				return ct;
				
				
	}
	public static ConnectionTable getConnectionTable(List<Tuple<Line2D,Integer>> linest, List<Shape> likelyNodes, double maxDistanceRatioNonLikely,double maxDistanceRatioLikely){
		return getConnectionTable(linest,likelyNodes,maxDistanceRatioNonLikely, maxDistanceRatioLikely);
		
	}
	
}
