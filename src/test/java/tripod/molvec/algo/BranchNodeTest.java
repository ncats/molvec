package tripod.molvec.algo;
import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
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
		assertEquals("-C(-O,=O)",BranchNode.interpretOCRStringAsAtom("CO2H").toString());
	}
	@Test
	public void terminagedEsterShouldHave3NodesWithCorrectBonds(){
		assertEquals("-C(-O,=O)",BranchNode.interpretOCRStringAsAtom("CO2").toString());
	}
	@Test
	public void methylEsterShouldHave4NodesWithCorrectBonds(){
		String s=BranchNode.interpretOCRStringAsAtom("CO2C").toString();
		
		assertEquals("-C(-O(-C),=O)",s);
	}
	
	@Test
	public void methyoxyEsterShouldHave5NodesWithCorrectBonds(){
		String s=BranchNode.interpretOCRStringAsAtom("CO2CO").toString();
		assertEquals("-C(-O(-C(-O)),=O)",s);
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
		assertEquals("-C(-O(-C(-F,-F,-F)),=O)",s);
	}
	
	public void testReadingSingle(String input, String expected){
		BranchNode bn= BranchNode.interpretOCRStringAsAtom(input);
		assertTrue("Single input '" + input + "' should be a real node",bn.isRealNode());
		assertEquals(expected,bn.getSymbol());
	}

}
