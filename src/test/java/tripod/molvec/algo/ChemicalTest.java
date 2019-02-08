package tripod.molvec.algo;

import java.io.IOException;

import org.junit.Test;

import gov.nih.ncats.chemkit.api.Chemical;
import gov.nih.ncats.chemkit.api.ChemicalBuilder;
import static org.junit.Assert.*;
public class ChemicalTest {

	@Test
	public void testReadAromaticSmiles() throws Exception{
		String s="CC(C)N1CCN(CC1)C(=O)[C@H](Cc2ccc3nc(=O)oc3c2)NC(=O)N4CC[C@@]5(CC4)NC(=O)Nc6ccccc56";
		Chemical c= Chemical.createFromSmiles(s);
		assertTrue(c.getAtomCount() >0);
	}
	
	@Test
	public void testChemicalDoesntAddStereo() throws IOException{
		String mol= "\n" + 
				"  CDK     02081909203D\n" + 
				"\n" + 
				" 20 23  0  0  0  0  0  0  0  0999 V2000\n" + 
				"   -3.4411   -0.8871    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   -1.6363    0.1230    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"    1.7654   -0.9142    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"    0.8607    0.6105    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"    1.7899    0.0910    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"    3.5083    0.1096    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"    0.9975   -1.5245    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   -2.4625   -1.5103    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"    2.5687    0.5905    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"    2.5225   -1.4327    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"    0.0000    0.1068    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"    4.2834   -0.4470    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"    3.4887   -0.9815    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   -1.6511   -0.8723    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   -0.9187   -1.5279    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   -0.8303    0.5858    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   -0.8540    1.5279    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   -3.4224    0.0551    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   -2.5454    0.6312    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   -4.2834   -0.3703    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  1 20  1  0  0  0  0 \n" + 
				"  1 18  1  0  0  0  0 \n" + 
				"  2 19  1  0  0  0  0 \n" + 
				"  2 16  1  0  0  0  0 \n" + 
				"  2 14  1  0  0  0  0 \n" + 
				"  3  5  1  0  0  0  0 \n" + 
				"  4  5  1  0  0  0  0 \n" + 
				"  1  8  1  0  0  0  0 \n" + 
				"  3  7  1  0  0  0  0 \n" + 
				" 12 13  1  0  0  0  0 \n" + 
				" 14 15  1  0  0  0  0 \n" + 
				" 16 17  2  0  0  0  0 \n" + 
				"  3 10  1  0  0  0  0 \n" + 
				" 18 19  1  0  0  0  0 \n" + 
				" 11 16  1  0  0  0  0 \n" + 
				"  5  9  1  0  0  0  0 \n" + 
				" 18 20  1  0  0  0  0 \n" + 
				"  6  9  1  0  0  0  0 \n" + 
				"  8 14  1  0  0  0  0 \n" + 
				"  6 13  1  0  0  0  0 \n" + 
				"  6 12  1  0  0  0  0 \n" + 
				"  4 11  1  0  0  0  0 \n" + 
				" 10 13  1  0  0  0  0 \n" + 
				"M  END";
		//Chemical c1= ChemicalBuilder.createFromMol(mol).build();
		Chemical c= Chemical.parseMol(mol);
		
		assertTrue(!c.toSmiles().contains("@"));
	}
}
