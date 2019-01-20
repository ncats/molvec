package tripod.molvec.algo;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.nio.charset.Charset;

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
		System.out.println("Complex1");
		File f=getFile("moleculeTest/complex.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("[H][C@@]12CN(C[C@]1([H])[C@H]2NCc3ccc4cc(F)ccc4n3)c5ncc(cn5)C(=O)NO").build();
		
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
	}
	
	@Test
	public void fuzzyStructure1Test() throws Exception {
		File f=getFile("moleculeTest/fuzzy1.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("ClC(=O)c1ccc(Oc2ccc(cc2)C(Cl)=O)cc1").build();
		
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
	}
	
	@Test
	public void fuzzyStructure2Test() throws Exception {
		File f=getFile("moleculeTest/fuzzy2.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("CC(=C)C(=O)OCCOC(=O)c1ccc(C(=O)OCCCOC(=O)C=C)c(c1)C(=O)OCC(O)COc2ccc(Cc3ccc(OCC4CO4)cc3)cc2").build();
		
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
	}
	
	@Test
	public void zerosForOxygensAndSmallInnerBondTest() throws Exception {
		File f=getFile("moleculeTest/withZerosAsOxygens.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("CC1C2C=CC1C(C2C(=O)OCC(C)=C)C(=O)OCC(C)=C").build();
		
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
	}
	
	@Test
	public void subscriptImplicitAtomsCl3Test() throws Exception {
		File f=getFile("moleculeTest/withSubscriptForCl.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("ClC(Cl)(Cl)c1nc(nc(n1)C(Cl)(Cl)Cl)-c2ccc3OCOc3c2").build();
		
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
	}
	
	@Test
	public void subscriptImplicitAtomsF3Test() throws Exception {
		File f=getFile("moleculeTest/withSubscriptForF.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("FC(F)(F)C1(N=N1)c2ccc(CN3C(=O)C=CC3=O)cc2").build();
		
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
	}
	
	@Test
	public void subscriptImplicitAtomsH2Test() throws Exception {
		File f=getFile("moleculeTest/withSubscriptForH.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("COc1cccc(C(O)c2cc(F)ccc2-N)c1C").build();
		
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
	}
	
	@Test
	public void moleculeWithCloseNitrogensInRingTest() throws Exception {
		File f=getFile("moleculeTest/moleculeWithCloseNitrogensInRing.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("FC(F)(F)CNc1nc(Nc2ccc(cc2)N3CCOCC3)nc4ccsc14").build();
		
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
	}
	@Test
	public void moleculeWith2CarboxyShortHandsTest() throws Exception {
		File f=getFile("moleculeTest/carboxylicShorthandNotation.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("OC(=O)Cc1ccc(OCc2ccccc2C(O)=O)cc1").build();
		
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
	}
	
	@Test
	public void nCarbonChainTest() throws Exception {
		File f=getFile("moleculeTest/carbonChainShorthand.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("CCCc1ccc(CCC)c2cc3c(-c4ccccc4)c5cc6c(CCC)ccc(CCC)c6cc5c(-c7ccccc7)c3cc12").build();
		
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
		
	}
	
	@Test
	public void serifFontTest() throws Exception {
		File f=getFile("moleculeTest/serifFont.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
				"  CDK     01171907453D\n" + 
				"\n" + 
				" 36 40  0  0  0  0  0  0  0  0999 V2000\n" + 
				"   61.6266 -114.0482    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   41.2500  -71.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   39.0000 -156.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   41.9100 -237.9100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   62.8606 -195.7677    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  155.0000 -272.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  115.5000 -249.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  239.7369 -269.2631    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  197.3154 -248.7265    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  271.9032  -71.9032    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  251.1223 -114.4840    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  197.3238  -59.3843    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  157.0000  -36.5000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  203.5000 -106.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  264.0000  -10.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  243.0968  -15.2258    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   49.0000 -299.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   69.0000 -294.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  107.5000 -202.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  301.0000 -262.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  294.0000 -241.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   19.0000 -241.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  294.0000  -69.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   19.0000  -69.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  243.0000 -292.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  272.5204 -237.0864    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"    7.0000  -49.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   70.0000  -16.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   73.2500  -39.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"   74.0610 -269.9695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  239.7773  -39.6384    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  116.4441  -59.7473    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  118.5000 -116.5000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  250.7303 -195.5577    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  192.0000 -193.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  274.0000 -154.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
				"  2 29  1  0  0  0  0 \n" + 
				" 26 34  1  0  0  0  0 \n" + 
				" 34 36  2  0  0  0  0 \n" + 
				" 32 33  1  0  0  0  0 \n" + 
				" 34 35  1  0  0  0  0 \n" + 
				"  7 19  2  0  0  0  0 \n" + 
				" 11 36  1  0  0  0  0 \n" + 
				"  4 30  2  0  0  0  0 \n" + 
				"  5 19  1  0  0  0  0 \n" + 
				" 15 16  1  0  0  0  0 \n" + 
				" 11 14  2  0  0  0  0 \n" + 
				" 18 30  1  0  0  0  0 \n" + 
				" 17 18  1  0  0  0  0 \n" + 
				"  8 25  1  0  0  0  0 \n" + 
				"  8 26  2  0  0  0  0 \n" + 
				"  9 35  1  0  0  0  0 \n" + 
				"  4 22  1  0  0  0  0 \n" + 
				" 16 31  1  0  0  0  0 \n" + 
				" 12 31  1  0  0  0  0 \n" + 
				"  2 24  1  0  0  0  0 \n" + 
				"  1 33  1  0  0  0  0 \n" + 
				" 28 29  1  0  0  0  0 \n" + 
				" 24 27  1  0  0  0  0 \n" + 
				" 10 23  1  0  0  0  0 \n" + 
				"  1  2  2  0  0  0  0 \n" + 
				"  1  3  1  0  0  0  0 \n" + 
				"  3  5  2  0  0  0  0 \n" + 
				"  4  5  1  0  0  0  0 \n" + 
				"  7 30  1  0  0  0  0 \n" + 
				" 12 13  2  0  0  0  0 \n" + 
				" 12 14  1  0  0  0  0 \n" + 
				" 29 32  2  0  0  0  0 \n" + 
				"  6  7  1  0  0  0  0 \n" + 
				"  6  9  2  0  0  0  0 \n" + 
				"  8  9  1  0  0  0  0 \n" + 
				" 21 26  1  0  0  0  0 \n" + 
				" 13 32  1  0  0  0  0 \n" + 
				" 10 31  2  0  0  0  0 \n" + 
				" 20 21  1  0  0  0  0 \n" + 
				" 10 11  1  0  0  0  0 \n" + 
				"M  END", Charset.defaultCharset()).build();
		//This smiles, which should be the same, doesn't seem to give the same formula. Don't know why
		//CCc1c(-C)c2cc3nc(nc4nc(cc5nc(cc1n2)c(-C)c5CC)c(-C)c4CC)c(-C)c3CC
		
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
		//
	}
	//
	@Test
	public void serifFont2Test() throws Exception {
		File f=getFile("moleculeTest/serifFont2.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("C=CC(C(c1ccc(OCCN2CCOCC2)cc1)c3ccc(OCCN4CCOCC4)cc3)c5ccccc5").build();
		
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
		//
	}
	
	//c1c(nn(c1-c2ccccc2)-c3cccc(c3)-n4c5ccc(cc5c6cc(ccc46)-c7ccccc7)-c8ccccc8)-c9ccccc9
	@Test
	public void structureWithWeirdAngleBondToOtherOCRAtomTest() throws Exception {
		File f=getFile("moleculeTest/weirdAngleBetweenNitrogens.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("c1c(nn(c1-c2ccccc2)-c3cccc(c3)-n4c5ccc(cc5c6cc(ccc46)-c7ccccc7)-c8ccccc8)-c9ccccc9").build();
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
		//
	}
	
	@Test
	public void structureWhichIsAlmostARingAndHas2NonBondedOxygensVeryCloseTogether() throws Exception {
		File f=getFile("moleculeTest/almostRing.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("CC(=C)C(=O)OCCOCC(COC(=O)CCC(=O)OCCCCOC(=O)C=C)OC(=O)c1ccccc1C(=O)OCC(O)COc2ccc(cc2)C(C)(C)c3ccc(OCC4CO4)cc3").build();
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
		//
	}
	@Test
	public void structureWhichWithOxygenOffCenterInRing() throws Exception {
		File f=getFile("moleculeTest/connectedOxygen.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("COc1cc2CN(C)c3c(ccc4cc5OCOc5cc34)-c2cc1OC").build();
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
		//
	}
	//
	//CC(=O)CC(C)=O
	
	@Test
	public void structureWithExplicitCarbons() throws Exception {
		File f=getFile("moleculeTest/explicitCarbons.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("CC(=O)CC(C)=O").build();
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
		//
	}
	//ClC(Cl)(Cl)c1nc(nc(n1)C(Cl)(Cl)Cl)-c2ccc3OCOc3c2
	@Test
	public void structureWithCarbonsThatAreSometimesMistakenForOxygens() throws Exception {
		File f=getFile("moleculeTest/carbonVsOxygen.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("ClC(Cl)(Cl)c1nc(nc(n1)C(Cl)(Cl)Cl)-c2ccc3OCOc3c2").build();
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
		//
	}
	
	@Test
	public void structureWithOxygensThatAreSometimesMistakenForCarbons() throws Exception {
		File f=getFile("moleculeTest/carbonVsOxygen2.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("CCOC(=O)CC1CCN(CC1)C(=O)CC2OC(c3cc(Cl)ccc3-n4cccc24)c5cccc6ccccc56").build();
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
		//
	}
	//
	@Test
	public void structureWithVeryCloseExplicitLinearAtoms() throws Exception {
		File f=getFile("moleculeTest/closeOCRShapesVeryExplicit.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("CCC(C)(C)COC(=O)c1ccc(cc1)C(=O)OC").build();
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
		//
	}
	//OC(c1cccc2OCCOc12)c3cc(Cl)ccc3-n4cccc4\C=C\C(=O)OCc5ccccc5
	@Test
	public void structureWithDoubelBondAngleSimilarToN() throws Exception {
		File f=getFile("moleculeTest/doubleBondSometimesSeenAsNAtom.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("OC(c1cccc2OCCOc12)c3cc(Cl)ccc3-n4cccc4\\C=C\\C(=O)OCc5ccccc5").build();
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
		//
	}
	//OCCCCCNc1ncccc1C(=O)Nc2ccc(cc2)-n3nc(cc3C(F)(F)F)-c4cccnc4
	@Test
	public void structureWithCloseNitrogens() throws Exception {
		File f=getFile("moleculeTest/nitrogensSometimesBondedByMistake.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("OCCCCCNc1ncccc1C(=O)Nc2ccc(cc2)-n3nc(cc3C(F)(F)F)-c4cccnc4").build();
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
		//
	}
	@Test
	public void structureWithVeryShortSingleBondBetweenCarbons() throws Exception {
		File f=getFile("moleculeTest/verySmallSingleBondBetweenExplicitCarbons.png");
		
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("C(C=Cc1ccc(cc1)N(c2ccccc2)c3ccc(cc3)-c4ccc(cc4)N(c5ccccc5)c6ccc(C=CC=Cc7ccccc7)cc6)=Cc8ccccc8").build();
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
		//
	}
	//CC(=C)C(=O)OCCOC(=O)c1ccc(C(=O)OCCOC(=O)C=C)c(c1)C(=O)OCC(O)COc2ccc(cc2)C(C)(C)c3ccc(OCC4CO4)cc3
	@Test
	public void structureWithGapInSmallRing() throws Exception {
		File f=getFile("moleculeTest/moleculeWithGapInSmallRing.png");
		
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		
		Chemical cReal=ChemicalBuilder.createFromSmiles("CC(=C)C(=O)OCCOC(=O)c1ccc(C(=O)OCCOC(=O)C=C)c(c1)C(=O)OCC(O)COc2ccc(cc2)C(C)(C)c3ccc(OCC4CO4)cc3").build();
		
		Chemical c=sie.getChemical();
		String form=c.getFormula();
		assertEquals(cReal.getFormula(),form);
		//
	}
}
