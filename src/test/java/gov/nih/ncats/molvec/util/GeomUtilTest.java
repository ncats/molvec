package gov.nih.ncats.molvec.util;

import static org.junit.Assert.*;
import static gov.nih.ncats.molvec.util.GeomUtil.area;
import static gov.nih.ncats.molvec.util.GeomUtil.areaTriangle;
import static gov.nih.ncats.molvec.util.GeomUtil.shapeFromVertices;
import static gov.nih.ncats.molvec.util.GeomUtil.vertices;

import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import gov.nih.ncats.molvec.util.GeomUtil;
import org.junit.Assert;
import org.junit.Test;

import gov.nih.ncats.molvec.algo.Tuple;


public class GeomUtilTest {

	
	@Test
	public void testMaxThetaForLineWorks(){
		double theta1 = 3*Math.PI/9;
		AffineTransform at= new AffineTransform();
		at.rotate(theta1);
		
		List<Point2D> pts=IntStream.range(0, 100)
		         .mapToObj(i->new Point2D.Double(i,0))
		         .map(p->at.transform(p, null))
		         .collect(Collectors.toList());
		
		Line2D max= GeomUtil.findMaxSeparationLine(pts);
		
		double theta=GeomUtil.angle(max.getP1(),max.getP2());
		
		assertEquals(theta1, theta,0.01);
	}
	
	
	@Test
	public void testPCAMaxThetaForLineWorks(){
		double noise = 0.0;
		int totSlice = 100;
		int npoints = 100;
		for(int i=0;i<totSlice;i++){
			double theta1 = 2*Math.PI*i/((double)totSlice);
			AffineTransform at= new AffineTransform();
			at.rotate(theta1);
			
			
			
			
			List<Point2D> pts=IntStream.range(0, npoints)
			         .mapToObj(j->new Point2D.Double(j-npoints/2,noise*(Math.random()-0.5)))
			         .map(p->at.transform(p, null))
			         .collect(Collectors.toList());
			
			double[] xs=pts.stream()
			.mapToDouble(p->p.getX())
			.toArray();
			double[] ys=pts.stream()
					.mapToDouble(p->p.getY())
					.toArray();
			
			double[] ln=GeomUtil.getPCALikeUnitVector(xs, ys);
			
			double tTheta = theta1;
			if(theta1>Math.PI){
				tTheta = theta1-2*Math.PI; // restrict between -PI and PI
			}
			if(tTheta<0){
				tTheta +=Math.PI;
			}
			
			
			double thetap1=GeomUtil.angle(new Point2D.Double(0,0),new Point2D.Double(ln[0],ln[1]));
			
			
			//can be off by PI in either direction
			try{
				assertEquals(tTheta, thetap1,0.01);
			}catch(Error e2){
				if(tTheta>0.01){
					assertEquals(tTheta-Math.PI, thetap1,0.01);
				}else{
					assertEquals(tTheta+Math.PI, thetap1,0.01);
				}
			}			
		}
	}
	
		
	
	
	@Test 
	public void testSplitShapeIntoNSides(){
		Rectangle2D rect = new Rectangle2D.Double(0,0,10,10);
		
		List<int[]> plist=GeomUtil.getNEquallySpacedPointsAroundShape(rect, 8)
		.stream()
		.map(p->new int[]{(int)p.getX(),(int)p.getY()})
		.collect(Collectors.toList());
		
		
		assertEquals(0,plist.get(0)[0]);
		assertEquals(0,plist.get(0)[1]);
		
		assertEquals(5,plist.get(1)[0]);
		assertEquals(0,plist.get(1)[1]);
		
		assertEquals(10,plist.get(2)[0]);
		assertEquals(0,plist.get(2)[1]);
		
		assertEquals(10,plist.get(3)[0]);
		assertEquals(5,plist.get(3)[1]);
		
		assertEquals(10,plist.get(4)[0]);
		assertEquals(10,plist.get(4)[1]);
		
		assertEquals(5,plist.get(5)[0]);
		assertEquals(10,plist.get(5)[1]);
		
		assertEquals(0,plist.get(6)[0]);
		assertEquals(10,plist.get(6)[1]);
		
		assertEquals(0,plist.get(7)[0]);
		assertEquals(5,plist.get(7)[1]);
		
		
	}
	
	
	
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
    public void transformLineShouldTurnSourceIntoDestination(){
    	Line2D o = new Line2D.Double(-2, -10, 500, 5);
    	Line2D d = new Line2D.Double(100, 1, 0, 0);
    	
    	AffineTransform at=GeomUtil.getTransformFromLineToLine(o, d,true);
    	
    	Line2D nline=GeomUtil.lines(at.createTransformedShape(o))[0];
    	
    	assertEquals(d.getX1(),nline.getX1(),0.001);
    	assertEquals(d.getY1(),nline.getY1(),0.001);
    	assertEquals(d.getX2(),nline.getX2(),0.001);
    	assertEquals(d.getY2(),nline.getY2(),0.001);
    	
    	
    }
    
//    @Test
//    public void pointJustOUtsideOfRectangleShouldHaveCorrectDistance(){
//    	Rectangle2D rect= new Rectangle2D.Double(0, 0, 10, 10).getBounds2D();
//    	Shape s=GeomUtil.convexHull2(vertices(rect));
//    	Arrays.stream(GeomUtil.lines(rect)).map(l->l.getP1() + "," + l.getP2()).forEach(System.out::println);
//    	
//    	//System.out.println(Arrays.toString());//
//    }
    
    @Test
    public void centerOfMassOfSquareShouldBeCenter(){
    	Rectangle2D rect = new Rectangle2D.Double(0,0,10,10);
    	
    	Point2D originalCenter = new Point2D.Double(5.0, 5.0);
    	
    	for(int j=1;j<10;j++){
	    	for(int i=1;i<10;i++){
		    	for(int x=0;x<10;x++){
		    		for(int y=0;y<10;y++){
		    			AffineTransform at = new AffineTransform();
		    			at.translate(x, y);
		    			at.scale(i, j);
		    			Shape trans=at.createTransformedShape(rect);
		    			Point2D transCenter=at.transform(originalCenter, null);
		    	
		
		    	    	Point2D center=GeomUtil.centerOfMass(trans);
		    	    	assertEquals(transCenter.getX(),center.getX(),0.0001);
		    	    	assertEquals(transCenter.getY(),center.getY(),0.0001);
		    			
		    		}
		    	}
	    	}
    	}
    }
    
    @Test
    public void veryRegularPolygonShouldHaveStrongCircleLikeScore(){
    	int totSections = 100;
    	Point2D[] pts = IntStream.range(0, totSections)
    			                 .mapToDouble(i->2*i*Math.PI/totSections)
    			                 .mapToObj(t->new Point2D.Double(100*Math.cos(t),100*Math.sin(t)))
    			                 .toArray(i->new Point2D[i]);
    	
    	Shape regPoly=GeomUtil.convexHull2(pts);
    	
    	double cscore=GeomUtil.getCircleLikeScore(regPoly);
    	
    	assertTrue(cscore>0.99);
    	
    	
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
    	Shape s = GeomUtil.convexHull2(new Point2D.Double(0, 0),new Point2D.Double(1, 0),new Point2D.Double(0, 1), new Point2D.Double(1,1));
    	
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
    public void convexHullShouldEliminateInnerPointsTest(){
		double[] coords = new double[] { 	2,  2,
											2,  2,
											10, 0,
											0, 0, 
			 								0, 0,
			 								0, 10,
			 								0, 10.002,
			 								0, 10};
		
		int[] icoords = Arrays.stream(coords).mapToInt(d->(int)(d+0.5)).toArray(); 
		
		Point2D[] pts = new Point2D[coords.length/2];
		for(int i=0;i<icoords.length;i+=2){
			pts[i/2]=new Point2D.Double(icoords[i], icoords[i+1]);
		}
		Shape s= GeomUtil.convexHull2(pts);
		
		for(Point2D p : pts){
			assertTrue(GeomUtil.distanceTo(s,p)<0.001);
		}
		assertTrue(s.contains(3, 3));
    }
    
    @Test
    public void convexHullDoublePrecissionShouldEliminateInnerPointsTest(){
		double[] coords = new double[] { 	.2,  .2,
											.2,  .2,
											 1, .0,
											.0, .0, 
			 								.0, .0,
			 								.0,  1,
			 								.0,  1.0002,
			 								.0,  1.0};
		
		double[] icoords = Arrays.stream(coords).toArray(); 
		
		Point2D[] pts = new Point2D[coords.length/2];
		for(int i=0;i<icoords.length;i+=2){
			pts[i/2]=new Point2D.Double(icoords[i], icoords[i+1]);
		}
		Shape s= GeomUtil.convexHull2(pts);
		System.out.println(Arrays.toString(vertices(s)));
		
		for(Point2D p : pts){
			assertTrue(GeomUtil.distanceTo(s,p)<0.001);
		}
		assertTrue(s.contains(.3, .3));
    }
    
    
    @Test
    public void traingleAreaTest(){
    	double expected= 0.5;
    	Shape si= GeomUtil.convexHull2(new Point2D[]{new Point2D.Double(2,1),new Point2D.Double(1,0),new Point2D.Double(2,0)});
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
    	Shape si= GeomUtil.convexHull2(new Point2D[]{new Point2D.Double(2,1),new Point2D.Double(1,0),new Point2D.Double(2,0)});
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
    public void convexHullShapeAreaTestWithShear(){
    	double expected= 5;
    	Random rshear = new Random(1234l);
    	
    	Shape si= shapeFromVertices(new Point2D[]{new Point2D.Double(4,2),new Point2D.Double(3,3),new Point2D.Double(2,2),new Point2D.Double(2,0),new Point2D.Double(4,0)});
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
    
    @Test
    public void nonConvexHullShapeAreaTestWithShear(){
    	double expected= 3;
    	Random rshear = new Random(1234l);
    	
    	Shape si= shapeFromVertices(new Point2D[]{new Point2D.Double(4,2),new Point2D.Double(3,1),new Point2D.Double(2,2),new Point2D.Double(2,0),new Point2D.Double(3.0,0),new Point2D.Double(4,0)});
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
