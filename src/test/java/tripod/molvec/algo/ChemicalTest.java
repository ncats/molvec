package tripod.molvec.algo;

import java.io.IOException;

import org.junit.Test;

import gov.nih.ncats.chemkit.api.Chemical;
import gov.nih.ncats.chemkit.api.ChemicalBuilder;

public class ChemicalTest {

	@Test
	public void testReadAromaticSmiles() throws Exception{
		String s="CC(C)N1CCN(CC1)C(=O)[C@H](Cc2ccc3nc(=O)oc3c2)NC(=O)N4CC[C@@]5(CC4)NC(=O)Nc6ccccc56";
		ChemicalBuilder.createFromSmiles(s).build();
	}
}
