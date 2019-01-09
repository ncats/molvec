package tripod.molvec.algo;

import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

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
	
//	public static ConnectionTable getConnectionTable(List<Tuple<Line2D,Integer>> lines, List<Point2D> likelyNodes, double){
//		//This will find the likely intersection points between lines.
//		//If 2 lines would intersect 
//	}
	
	public static class ConnectionTable{
		Point2D[] nodes;
		List<Edge> edges;
		
		
	}
	public static class Edge{
		int n1;
		int n2;
		int order;
	}
}
