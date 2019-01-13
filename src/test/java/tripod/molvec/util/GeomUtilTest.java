package tripod.molvec.util;

import static org.junit.Assert.*;

import java.awt.geom.Line2D;
import java.awt.geom.Point2D;

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
}
