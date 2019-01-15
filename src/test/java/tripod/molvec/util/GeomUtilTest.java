package tripod.molvec.util;

import static org.junit.Assert.assertEquals;

import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.Test;


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
}
