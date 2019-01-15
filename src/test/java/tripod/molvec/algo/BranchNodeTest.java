package tripod.molvec.algo;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;
public class BranchNodeTest {
	
	@Test
	public void commonAtomsShouldGetCorrectReading(){
		testReadingSingle("C","C");
		testReadingSingle("O","O");
		testReadingSingle("H","H");
		testReadingSingle("N","N");
		testReadingSingle("Br","Br");
		testReadingSingle("F","F");
		testReadingSingle("S","S");
		testReadingSingle("P","P");
		testReadingSingle("Cl","Cl");
	}
	
	@Test
	public void carboxylicAcidShouldHave3NodesWithCorrectBonds(){
		assertEquals("-C(=O,-O)",BranchNode.interpretOCRStringAsAtom("CO2H").toString());
	}
	@Test
	public void terminagedEsterShouldHave3NodesWithCorrectBonds(){
		assertEquals("-C(=O,-O)",BranchNode.interpretOCRStringAsAtom("CO2").toString());
	}
	@Test
	public void methylEsterShouldHave4NodesWithCorrectBonds(){
		String s=BranchNode.interpretOCRStringAsAtom("CO2C").toString();
		
		assertEquals("-C(=O,-O(-C))",s);
	}
	
	@Test
	public void methyoxyEsterShouldHave5NodesWithCorrectBonds(){
		String s=BranchNode.interpretOCRStringAsAtom("CO2CO").toString();
		assertEquals("-C(=O,-O(-C(-O)))",s);
	}
	
	@Test
	public void floridatedMethaneShouldHave4NodesWithCorrectBonds(){
		String s=BranchNode.interpretOCRStringAsAtom("CF3").toString();
		assertEquals("-C(-F,-F,-F)",s);
	}
	
	
	@Test
	public void chloridatedMethaneShouldHave4NodesWithCorrectBonds(){
		String s=BranchNode.interpretOCRStringAsAtom("CCl3").toString();
		assertEquals("-C(-Cl,-Cl,-Cl)",s);
	}
	
	@Test
	public void floridatedInverseMethaneShouldHave4NodesWithCorrectBonds(){
		String s=BranchNode.interpretOCRStringAsAtom("F3C").toString();
		assertEquals("-C(-F,-F,-F)",s);
	}
	
	@Test
	public void floridatedmethylEsterShouldHave4NodesWithCorrectBonds(){
		String s=BranchNode.interpretOCRStringAsAtom("CO2CF3").toString();
		assertEquals("-C(=O,-O(-C(-F,-F,-F)))",s);
	}
	
	@Test
	public void shouldLookNice(){
		BranchNode s=BranchNode.interpretOCRStringAsAtom("CCF3");
		s.generateCoordinates();
		List<String> coordStrings = new ArrayList<>();
		
		s.forEachBranchNode((a,b)->{
			Point2D p=b.suggestedPoint;
			coordStrings.add(b.getSymbol() + " :"+ p.getX() + "," + p.getY());
		});
		assertEquals("C :0.0,0.0",coordStrings.get(0));
		assertEquals("C :0.5000000000000001,-0.8660254037844386",coordStrings.get(1));
		assertEquals("F :3.3306690738754696E-16,-1.7320508075688774",coordStrings.get(2));
		assertEquals("F :1.5,-0.8660254037844386",coordStrings.get(3));
		assertEquals("F :1.0000000000000002,-1.7320508075688772",coordStrings.get(4));
	}
	
	public void testReadingSingle(String input, String expected){
		BranchNode bn= BranchNode.interpretOCRStringAsAtom(input);
		assertTrue("Single input '" + input + "' should be a real node",bn.isRealNode());
		assertEquals(expected,bn.getSymbol());
	}

}
