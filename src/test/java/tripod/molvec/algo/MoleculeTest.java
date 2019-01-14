package tripod.molvec.algo;

import static org.junit.Assert.assertEquals;

import java.io.File;

import org.junit.Test;

import gov.nih.ncats.chemkit.api.Chemical;
import gov.nih.ncats.chemkit.api.ChemicalBuilder;

public class MoleculeTest {

	
	private File getFile(String fname){
		ClassLoader classLoader = getClass().getClassLoader();
		return new File(classLoader.getResource(fname).getFile());
		
	}
	
	@Test
	public void fluoxetineWikiTest() throws Exception {
		File f=getFile("moleculeTest/fluoxetine.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals("C17H18F3NO",form);
	}
	
	@Test
	public void gleevecWikiTest() throws Exception {
		File f=getFile("moleculeTest/gleevec.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals("C29H31N7O",form);
	}
	
	@Test
	public void tylenolWikiTest() throws Exception {
		File f=getFile("moleculeTest/tylenol.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals("C8H9NO2",form);
	}
	
	@Test
	public void paxilWikiTest() throws Exception {
		File f=getFile("moleculeTest/paxil.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals("C19H20FNO3",form);
	}
	
	@Test
	public void complexStructure1Test() throws Exception {
		File f=getFile("moleculeTest/complex.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("[H][C@@]12CN(C[C@]1([H])[C@H]2NCc3ccc4cc(F)ccc4n3)c5ncc(cn5)C(=O)NO").build();
		
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
	}
	
	
	
}