package tripod.molvec.util;

import static org.junit.Assert.assertEquals;

import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import static tripod.molvec.util.GeomUtil.*;

import org.junit.Test;

import tripod.molvec.algo.Tuple;


public class GeomUtilTest {

    @Test
    public void projectingPointOntoHorizontalLineShouldReturnSameXCoordinate() throws Exception {
    	
    	Line2D line1 = new Line2D.Double(0,0,100,0);
    	
    	for(int i=-100;i<100;i++){
    		Point2D p = new Point2D.Double(i, i);
    		Point2D proj=GeomUtil.projectPointOntoLine(line1,p);
    		assertEquals(0,proj.getY(), 0.0001);
    		assertEquals(p.getX(),proj.getX(), 0.0001);
    	}
    	
    	for(int i=-100;i<100;i++){
    		Point2D p = new Point2D.Double(i, i*4);
    		Point2D proj=GeomUtil.projectPointOntoLine(line1,p);
    		assertEquals(0,proj.getY(), 0.0001);
    		assertEquals(p.getX(),proj.getX(), 0.0001);
    	}
    }
    @Test
    public void cosThetaOfOrthoLinesShouldBeZero() throws Exception {
    	
    	Line2D line1 = new Line2D.Double(0,0,100,0);
    	Line2D line2 = new Line2D.Double(0,0,0,100);
    	double cos=GeomUtil.cosTheta(line1, line2);
    	assertEquals(0,cos, 0.0001);
    }
    
    @Test
    public void testLineGroupingsBasedOnRejectionWorks(){
    	
    	List<Line2D> lines = new ArrayList<>();
    	
    	lines.add(new Line2D.Double(0,0,1000,0));
    	lines.add(new Line2D.Double(0,1,10,1));
    	lines.add(new Line2D.Double(100,2,0,2));
    	lines.add(new Line2D.Double(0,3,1,3));
    	lines.add(new Line2D.Double(0,4,0.5,4));
    	lines.add(new Line2D.Double(10,5,500,5));
    	lines.add(new Line2D.Double(10,6,500,6));
    	
    	int[] rejections=GeomUtil.getLineOffsetsToLongestLine(lines)
    	.stream()
    	.map(t->t.v())
    	.sorted()
    	.mapToInt(d->(int)Math.round(d))
    	.toArray();
    	
    	int[] expected=new int[]{0,1,2,3,4,5,6};
    	for(int i=0;i<expected.length;i++){
    		assertEquals(expected[i],rejections[i]);	
    	}
    	
    	
    }
    
    @Test
    public void testLineGroupingsBasedOnRejectionRotatedWorks(){
    	
    	int nslices=360;
    	for(int n=0;n<nslices;n++){
    		double theta= n*(2*Math.PI)/((double)nslices);
	    	List<Line2D> lines = new ArrayList<>();
	    	
	    	lines.add(new Line2D.Double(0,0,1000,0));
	    	lines.add(new Line2D.Double(0,1,10,1));
	    	lines.add(new Line2D.Double(100,2,0,2));
	    	lines.add(new Line2D.Double(0,3,1,3));
	    	lines.add(new Line2D.Double(0,4,0.5,4));
	    	lines.add(new Line2D.Double(10,5,500,5));
	    	lines.add(new Line2D.Double(10,6,500,6));
	    	
	    	AffineTransform at = new AffineTransform();
	    	at.rotate(theta);
	    	lines=lines.stream()
	    	     .map(l->{
	    	    	Point2D[] pts = new Point2D[]{l.getP1(),l.getP2()};
	    	    	pts[0]=at.transform(pts[0], null);
	    	    	pts[1]=at.transform(pts[1], null);
	    	    	return new Line2D.Double(pts[0],pts[1]);
	    	     })
	    	     .collect(Collectors.toList());
	    	
	    	
	    	int[] rejections=GeomUtil.getLineOffsetsToLongestLine(lines)
	    	.stream()
	    	.map(t->t.v())
	    	.sorted()
	    	.mapToInt(d->(int)Math.round(d))
	    	.toArray();
	    	
	    	int[] expected=new int[]{0,1,2,3,4,5,6};
	    	assertEquals(Arrays.toString(expected), Arrays.toString(rejections));
	    	
    	}
    }
    @Test
    public void testLineGroupingsBasedOnRejectionTranslateWorks(){
    	
    	int nxslices=10;
    	int nyslices=10;
    	for(int x=0;x<nxslices;x++){
    		for(int y=0;y<nyslices;y++){
	    		List<Line2D> lines = new ArrayList<>();
		    	
		    	lines.add(new Line2D.Double(0,0,1000,0));
		    	lines.add(new Line2D.Double(0,1,10,1));
		    	lines.add(new Line2D.Double(100,2,0,2));
		    	lines.add(new Line2D.Double(0,3,1,3));
		    	lines.add(new Line2D.Double(0,4,0.5,4));
		    	lines.add(new Line2D.Double(10,5,500,5));
		    	lines.add(new Line2D.Double(10,6,500,6));
		    	
		    	AffineTransform at = new AffineTransform();
		    	at.translate(x, y);
		    	lines=lines.stream()
		    	     .map(l->{
		    	    	Point2D[] pts = new Point2D[]{l.getP1(),l.getP2()};
		    	    	pts[0]=at.transform(pts[0], null);
		    	    	pts[1]=at.transform(pts[1], null);
		    	    	return new Line2D.Double(pts[0],pts[1]);
		    	     })
		    	     .collect(Collectors.toList());
		    	
		    	
		    	int[] rejections=GeomUtil.getLineOffsetsToLongestLine(lines)
		    	.stream()
		    	.map(t->t.v())
		    	.sorted()
		    	.mapToInt(d->(int)Math.round(d))
		    	.toArray();
		    	
		    	int[] expected=new int[]{0,1,2,3,4,5,6};
		    	assertEquals(Arrays.toString(expected), Arrays.toString(rejections));
    		}
    	}
    }
    
    @Test
    public void testLineGroupingsBasedOnScaleWorks(){
    	
    	int nxslices=100;
    	for(int x=1;x<nxslices;x++){
	    		List<Line2D> lines = new ArrayList<>();
		    	
		    	lines.add(new Line2D.Double(0,0,1000,0));
		    	lines.add(new Line2D.Double(0,1,10,1));
		    	lines.add(new Line2D.Double(100,2,0,2));
		    	lines.add(new Line2D.Double(0,3,1,3));
		    	lines.add(new Line2D.Double(0,4,0.5,4));
		    	lines.add(new Line2D.Double(10,5,500,5));
		    	lines.add(new Line2D.Double(10,6,500,6));
		    	
		    	AffineTransform at = new AffineTransform();
		    	at.scale(nxslices, nxslices);
		    	lines=lines.stream()
		    	     .map(l->{
		    	    	Point2D[] pts = new Point2D[]{l.getP1(),l.getP2()};
		    	    	pts[0]=at.transform(pts[0], null);
		    	    	pts[1]=at.transform(pts[1], null);
		    	    	return new Line2D.Double(pts[0],pts[1]);
		    	     })
		    	     .collect(Collectors.toList());
		    	
		    	
		    	int[] rejections=GeomUtil.getLineOffsetsToLongestLine(lines)
		    	.stream()
		    	.map(t->t.v())
		    	.sorted()
		    	.mapToInt(d->(int)Math.round(d))
		    	.toArray();
		    	
		    	int[] expected=new int[]{0,1*nxslices,2*nxslices,3*nxslices,4*nxslices,5*nxslices,6*nxslices};
		    	assertEquals(Arrays.toString(expected), Arrays.toString(rejections));
    	}
    }
    @Test
    public void testLineGroupingsBasedOnReflectionWorks(){
    	List<Line2D> lines = new ArrayList<>();
    	
    	lines.add(new Line2D.Double(0,0,1000,0));
    	lines.add(new Line2D.Double(0,1,10,1));
    	lines.add(new Line2D.Double(100,2,0,2));
    	lines.add(new Line2D.Double(0,3,1,3));
    	lines.add(new Line2D.Double(0,4,0.5,4));
    	lines.add(new Line2D.Double(10,5,500,5));
    	lines.add(new Line2D.Double(10,6,500,6));
    	
    	AffineTransform at = new AffineTransform();
    	at.scale(1, -1); // reflect
    	
    	lines=lines.stream()
    	     .map(l->{
    	    	Point2D[] pts = new Point2D[]{l.getP1(),l.getP2()};
    	    	pts[0]=at.transform(pts[0], null);
    	    	pts[1]=at.transform(pts[1], null);
    	    	return new Line2D.Double(pts[0],pts[1]);
    	     })
    	     .collect(Collectors.toList());
    	
    	
    	int[] rejections=GeomUtil.getLineOffsetsToLongestLine(lines)
    	.stream()
    	.map(t->t.v())
    	.sorted()
    	.mapToInt(d->(int)Math.round(d))
    	.toArray();
    	
    	int[] expected=new int[]{-6,-5,-4,-3,-2,-1,0};
    	assertEquals(Arrays.toString(expected), Arrays.toString(rejections));
    }
    
    @Test
    public void testGroupPairsWorks(){
    	
    	List<Integer> ilist1 = IntStream.range(0,100).mapToObj(i->i).collect(Collectors.toList());
    	List<Integer> ilist2 = IntStream.range(102,150).mapToObj(i->i).collect(Collectors.toList());
    	List<Integer> ilist3 = IntStream.range(155,157).mapToObj(i->i).collect(Collectors.toList());
    	
    	List<Integer> ilist = Stream.concat(ilist1.stream(),
    								        Stream.concat(ilist2.stream(),
    								        			  ilist3.stream()))
							    	.collect(Collectors.toList());
    	
    	
    	Collections.shuffle(ilist, new Random(1234l));

    	List<List<Integer>> ilistg = GeomUtil.groupThings(ilist, t->(Math.abs(t.k()-t.v())<=1))
    	        .stream()
    	        .map(l->l.stream().sorted().collect(Collectors.toList()))
    	        .map(l->Tuple.of(l,l.get(0)).withVComparator())
    	        .sorted()
    	        .map(t->t.k())
    	        .collect(Collectors.toList());
    	
    	
    	assertEquals(ilist1,ilistg.get(0));
    	assertEquals(ilist2,ilistg.get(1));
    	assertEquals(ilist3,ilistg.get(2));
    	
    	
    }
    
    @Test
    public void testIntersectingShape(){
    	Shape s = GeomUtil.convexHull(new Point2D.Double(0, 0),new Point2D.Double(1, 0),new Point2D.Double(0, 1), new Point2D.Double(1,1));
    	
    	AffineTransform at = new AffineTransform();
    	at.scale(10, 10);
    	
    	Shape rs=at.createTransformedShape(s);
    	at = new AffineTransform();
    	at.translate(5, 0);
    	Shape ns = at.createTransformedShape(rs);
    	
    	
    	Shape is=GeomUtil.getIntersectionShape(rs, ns).get();
    	assertEquals(50,0, area(is));
    	
    	
    	
    }
    @Test
    public void traingleAreaTest(){
    	double expected= 0.5;
    	Shape si= convexHull(new Point2D[]{new Point2D.Double(2,1),new Point2D.Double(1,0),new Point2D.Double(2,0)});
    	int nrot=100;
    	for(int k=0;k<2;k++){
    		AffineTransform atflip =new AffineTransform();
    		int inv=((k%2)==0)?1:-1;
    		atflip.scale(inv, 1);
    		
    		Shape s=atflip.createTransformedShape(si);
	    	for(int j=0;j<nrot;j++){
	    		AffineTransform atrotate =new AffineTransform();
	    		atrotate.rotate(j*Math.PI/nrot);
	    		Shape rs=atrotate.createTransformedShape(s);
		    	for(int i=1;i<100;i++){
			    	AffineTransform at =new AffineTransform();
			    	at.scale(i, i);
			    	Shape sn=at.createTransformedShape(rs);
			    	double area=areaTriangle(sn);
			    	assertEquals(inv*expected*i*i,area,0.00001);
		    	}
	    	}
    	}
    }
    
    @Test
    public void traingleAreaTest2(){
    	double expected= 0.5;
    	Shape si= convexHull(new Point2D[]{new Point2D.Double(2,1),new Point2D.Double(1,0),new Point2D.Double(2,0)});
    	int nrot=100;
    	for(int k=0;k<2;k++){
    		AffineTransform atflip =new AffineTransform();
    		int inv=((k%2)==0)?1:-1;
    		atflip.scale(inv, 1);
    		
    		Shape s=atflip.createTransformedShape(si);
	    	for(int j=0;j<nrot;j++){
	    		AffineTransform atrotate =new AffineTransform();
	    		atrotate.rotate(j*Math.PI/nrot);
	    		Shape rs=atrotate.createTransformedShape(s);
		    	for(int i=1;i<100;i++){
			    	AffineTransform at =new AffineTransform();
			    	at.scale(i, i);
			    	Shape sn=at.createTransformedShape(rs);
			    	double area=GeomUtil.areaVerticesCW(vertices(sn));
			    	assertEquals(inv*expected*i*i,area,0.00001);
		    	}
	    	}
    	}
    }
    
    @Test
    public void rectangleAreaTestWithShear(){
    	double expected= 1;
    	Random rshear = new Random(1234l);
    	Shape si= convexHull(new Point2D[]{new Point2D.Double(2,1),new Point2D.Double(1,1),new Point2D.Double(1,0),new Point2D.Double(2,0)});
    	int nrot=100;
    	for(int k=0;k<2;k++){
    		AffineTransform atflip =new AffineTransform();
    		int inv=((k%2)==0)?1:-1;
    		atflip.shear(rshear.nextDouble()*1000, 0);
    		atflip.scale(inv, 1);
    		
    		Shape s=atflip.createTransformedShape(si);
	    	for(int j=0;j<nrot;j++){
	    		AffineTransform atrotate =new AffineTransform();
	    		atrotate.rotate(j*Math.PI/nrot);
	    		Shape rs=atrotate.createTransformedShape(s);
		    	for(int i=1;i<100;i++){
			    	AffineTransform at =new AffineTransform();
			    	at.scale(i, i);
			    	Shape sn=at.createTransformedShape(rs);
			    	double area=GeomUtil.areaVerticesCW(vertices(sn));
			    	assertEquals(inv*expected*i*i,area,0.00001);
		    	}
	    	}
    	}
    }
    
    
    
    
    
//    
//    @Test
//    public void f(){
//    	for(int runs=1;runs<100;runs++){
//    		int reps=1000;
//    		List<Double> counts = new ArrayList<Double>();
//    		for(int j=0;j<reps;j++){
//		    	List<Line2D> lines=new ArrayList<Line2D>();
//		    	for(int i=0;i<runs;i++){
//		    		double t=Math.PI*2*Math.random();
//		    		lines.add(new Line2D.Double(0,0,Math.cos(t),
//		    				                    Math.sin(t)));
//		    		
//		    	}
//		    	List<List<Line2D>> parLines=GeomUtil.groupMultipleBonds(lines, 3*Math.PI/180, Double.POSITIVE_INFINITY, 0, 0);
//		    	//System.out.println(runs + "\t" + parLines.size());
//		    	counts.add(parLines.size()+0.0);
//    		}
//    		double avg=counts.stream().mapToDouble(d->d).average().getAsDouble();
//    		double sumsq=counts.stream().mapToDouble(d->d*d).average().getAsDouble();
//    		double stDev=Math.sqrt(sumsq-avg*avg);
//    		System.out.println(runs + "\t"+avg + "\t" + stDev);
//    		//sumsq-
//    	}
//    	
//    }
    
//    @Test
//    public void testAreaForShapeWorks(){
//    	Shape s = GeomUtil.convexHull(new Point2D.Double(0, 0),new Point2D.Double(1, 0),new Point2D.Double(0, 1));
//    	System.out.println("Area:");
//    	System.out.println(GeomUtil.poorMansArea(s));
//    }
//    
//
//    @Test
//    public void testSplitInHalf(){
//    	int segs=700;
//    	double size=1000;
//    	Point2D[] circle = IntStream.range(0,segs)
//    			                    .mapToDouble(i->2*Math.PI*(((double)i)/((double)segs)))
//    			                    .mapToObj(t->new Point2D.Double(size*Math.cos(t),size*Math.sin(t)))
//    			                    .toArray(i->new Point2D[i]);
//    	
//    	Collections.shuffle(Arrays.asList(circle));
//    	System.out.println(Arrays.toString(circle));
//    	
//    	Shape s = GeomUtil.convexHull(circle);
//    	Point2D[] pts= GeomUtil.vertices(s);
//    	
//    	Arrays.stream(pts)
//    	.map(p->Math.atan2(p.getY(), p.getX()))
//    	.forEach(d->{
//    		System.out.println(d);
//    	});
    	
    	
    	//System.out.println(Arrays.toString(GeomUtil.vertices(s)));
//    	Shape[] nshapes = GeomUtil.splitInHalf(s);
//    	
//    	System.out.println(Arrays.toString(GeomUtil.vertices(nshapes[0])));
//    	System.out.println(Arrays.toString(GeomUtil.vertices(nshapes[1])));
//   }
}
