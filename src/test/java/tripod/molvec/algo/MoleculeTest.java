package tripod.molvec.algo;

import static org.junit.Assert.assertEquals;

import java.io.*;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.function.Consumer;

import gov.nih.ncats.chemkit.api.util.stream.ThrowingStream;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import gov.nih.ncats.chemkit.api.Chemical;
import gov.nih.ncats.chemkit.api.ChemicalBuilder;
import gov.nih.ncats.chemkit.api.inchi.Inchi;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

@RunWith(Parameterized.class)
public class MoleculeTest {


	private File getFile(String fname){
		ClassLoader classLoader = getClass().getClassLoader();
		return new File(classLoader.getResource(fname).getFile());

	}

	public static class TestSpec{
		private String filePath;
		private ThrowingStream.ThrowingConsumer<Chemical, Exception> assertionConsumer;

		public TestSpec(String filePath, ThrowingStream.ThrowingConsumer<Chemical, Exception> assertionConsumer) {
			this.filePath = filePath;
			this.assertionConsumer = assertionConsumer;
		}
	}





	private TestSpec spec;

	public MoleculeTest(String ignored, TestSpec spec){
		this.spec = spec;
	}


	@Test
	public void testAsFile() throws Exception {
		File f=getFile(spec.filePath);

		StructureImageExtractor sie = new StructureImageExtractor(f);
		spec.assertionConsumer.accept(sie.getChemical());
	}
	@Test
	public void testAsByteArray() throws Exception {
		File f=getFile(spec.filePath);


		ByteArrayOutputStream out = new ByteArrayOutputStream();
		try(InputStream in = new BufferedInputStream(new FileInputStream(f))){
			byte[] buf = new byte[1024];

			int bytesRead=0;
			while( (bytesRead = in.read(buf)) > 0){
				out.write(buf, 0, bytesRead);
			}
		}

		StructureImageExtractor sie = new StructureImageExtractor(out.toByteArray());
		spec.assertionConsumer.accept(sie.getChemical());
	}

	@Parameterized.Parameters(name = "{0}")
	public static List<Object[]> data(){
		List<Object[]> list = new ArrayList<>();

		list.add(new Object[]{"fluoxetineWikiTest", new TestSpec("moleculeTest/fluoxetine.png", c-> assertEquals("C17H18F3NO",c.getFormula()))});
		list.add(new Object[]{"reallyHardPeptideTest", new TestSpec("moleculeTest/peptideFragment.png", c-> {
			String key=Inchi.asStdInchi(c).getKey();
			//String form=c.getFormula();
			assertEquals("PLIFXMBNBJXWIM-MUGJNUQGSA-N",key);
		})});

		list.add(new Object[]{"lipitorWikiTest", new TestSpec("moleculeTest/lipitor.png", c-> {
			String key=Inchi.asStdInchi(c).getKey();
			//String form=c.getFormula();
			assertEquals("XUKUURHRXDUEBC-KAYWLYCHSA-N",key);
		})});
//ringSystemWithHInCenter.png
		list.add(new Object[]{"gleevecWikiTest", new TestSpec("moleculeTest/gleevec.png", c-> assertEquals("C29H31N7O",c.getFormula()))});
		list.add(new Object[]{"tylenolWikiTest", new TestSpec("moleculeTest/tylenol.png", c-> assertEquals("C8H9NO2",c.getFormula()))});

		list.add(new Object[]{"paxilWikiTest", new TestSpec("moleculeTest/paxil.png", c-> assertEquals("C19H20FNO3",c.getFormula()))});

		//phenylShorhand.png
		list.add(new Object[]{"phenylShorthand", new TestSpec("moleculeTest/phenylShorhand.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("CC(C)(C)NS(=O)(=O)C2(Cc1ccccc1)CC2").build();
			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		//cagedStructure.png
		list.add(new Object[]{"cagedStructure", new TestSpec("moleculeTest/cagedStructure.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("CC(=C)C(=O)OC1C2CC3C1OC(=O)C3C2").build();
			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		//NConnectedToDoubleBond.png
		list.add(new Object[]{"NConnectedToDoubleBond", new TestSpec("moleculeTest/NConnectedToDoubleBond.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"  MJ150420                      \n" + 
					"\n" + 
					" 31 34  0  0  0  0  0  0  0  0999 V2000\n" + 
					"    1.4593    2.1617    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.2559    2.4019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3316    1.3347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.0248    0.9364    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.7844    0.9242    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.4776    1.3377    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.0431   -2.3106    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3195   -1.9093    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.6080   -2.7241    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.1398   -2.0684    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3377   -1.0945    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.5229   -1.0945    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.2649    1.1557    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.7302    1.7391    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.1525   -1.0945    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.0710    1.3289    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.6079   -2.3660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.1520   -3.3565    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.0309    2.5478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.7302    2.1404    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.6451    3.1011    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.7559    2.5598    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.0790    2.1404    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3377    2.1404    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.7661    3.3565    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3425   -0.2756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.6451    0.3770    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.6263   -3.1133    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.6501    0.9334    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.0248    0.1337    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.6263    0.1337    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					" 26 31  2  0  0  0  0\n" + 
					" 11 12  2  0  0  0  0\n" + 
					" 13 14  2  0  0  0  0\n" + 
					" 11 15  2  0  0  0  0\n" + 
					"  4 30  2  0  0  0  0\n" + 
					" 19 20  1  0  0  0  0\n" + 
					" 16 29  1  0  0  0  0\n" + 
					" 17 28  1  0  0  0  0\n" + 
					" 19 24  1  0  0  0  0\n" + 
					"  5 16  2  0  0  0  0\n" + 
					"  9 18  1  0  0  0  0\n" + 
					" 22 23  1  0  0  0  0\n" + 
					" 26 30  1  0  0  0  0\n" + 
					" 22 25  2  0  0  0  0\n" + 
					"  2 21  1  0  0  0  0\n" + 
					"  1  2  2  0  0  0  0\n" + 
					"  1  6  1  0  0  0  0\n" + 
					"  3  4  1  0  0  0  0\n" + 
					"  2 14  1  0  0  0  0\n" + 
					" 13 27  1  0  0  0  0\n" + 
					" 11 26  1  0  0  0  0\n" + 
					"  5  6  1  0  0  0  0\n" + 
					"  3 29  2  0  0  0  0\n" + 
					" 29 31  1  0  0  0  0\n" + 
					"  7  8  1  0  0  0  0\n" + 
					" 10 17  2  0  0  0  0\n" + 
					"  9 10  1  0  0  0  0\n" + 
					"  6 13  1  0  0  0  0\n" + 
					" 18 28  2  0  0  0  0\n" + 
					"  3 24  1  0  0  0  0\n" + 
					"  8 17  1  0  0  0  0\n" + 
					" 16 23  1  0  0  0  0\n" + 
					"  1 22  1  0  0  0  0\n" + 
					"  8 11  1  0  0  0  0\n" + 
					"M  END",Charset.defaultCharset()).build();
			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		//withAtomNumbers.png
		list.add(new Object[]{"withAtomNumbers", new TestSpec("moleculeTest/withAtomNumbers.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("[#6]Oc1ccc(cc1)C2(Oc3cc(cc(C(=O)NC\\C=C\\[C@H]4O[C@H]([C@@H]5OC([#6])([#6])O[C@H]45)n6cnc7c(-[#7])ncnc67)c3O2)C(F)(F)F)c8ccc(O[#6])cc8").build();
			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		
		list.add(new Object[]{"colinearTripleBond", new TestSpec("moleculeTest/colinearTripleBond.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"  CDK     01281921493D\n" + 
					"\n" + 
					" 36 39  0  0  0  0  0  0  0  0999 V2000\n" + 
					"  241.0000 -311.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  671.0000   -3.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  616.0000  -42.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  779.5000   -9.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  953.0000 -107.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  422.0000 -291.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  362.7500 -323.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  536.5000 -163.0000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  535.5000  -97.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  535.5000 -230.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  469.0000 -163.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  604.0000 -163.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  205.0155 -257.8492    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  241.0000 -204.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  724.0000  -42.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  420.0000 -204.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  636.0000 -219.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  304.0000 -290.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  703.0000 -105.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  636.2500 -105.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  534.5000 -424.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  592.0000 -391.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  304.7500 -223.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  478.0000 -390.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   36.0000 -314.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  101.0000 -314.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  362.0000 -389.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  478.0000 -324.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  420.0000 -423.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  135.2500 -255.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  101.0000 -199.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   35.5000 -199.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.0000 -255.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  354.0000 -182.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  837.3333  -42.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  895.1667  -74.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1 18  1  0  0  0  0 \n" + 
					" 32 33  2  0  0  0  0 \n" + 
					"  1 13  1  0  0  0  0 \n" + 
					" 15 19  1  0  0  0  0 \n" + 
					" 13 14  2  0  0  0  0 \n" + 
					"  3 20  1  0  0  0  0 \n" + 
					" 11 16  1  0  0  0  0 \n" + 
					" 19 20  1  0  0  0  0 \n" + 
					" 30 31  2  0  0  0  0 \n" + 
					"  7 18  1  0  0  0  0 \n" + 
					"  6 28  2  0  0  0  0 \n" + 
					" 16 34  1  0  0  0  0 \n" + 
					" 26 30  1  0  0  0  0 \n" + 
					" 24 28  1  0  0  0  0 \n" + 
					" 24 29  2  0  0  0  0 \n" + 
					"  2  3  1  0  0  0  0 \n" + 
					" 31 32  1  0  0  0  0 \n" + 
					" 25 33  1  0  0  0  0 \n" + 
					"  2 15  1  0  0  0  0 \n" + 
					" 23 34  1  0  0  0  0 \n" + 
					"  6  7  1  0  0  0  0 \n" + 
					" 12 20  1  0  0  0  0 \n" + 
					" 12 17  1  0  0  0  0 \n" + 
					" 13 30  1  0  0  0  0 \n" + 
					"  8  9  2  0  0  0  0 \n" + 
					"  8 10  2  0  0  0  0 \n" + 
					"  8 12  1  0  0  0  0 \n" + 
					"  7 27  2  0  0  0  0 \n" + 
					" 18 23  2  0  0  0  0 \n" + 
					" 21 24  1  0  0  0  0 \n" + 
					" 14 23  1  0  0  0  0 \n" + 
					"  4 15  1  0  0  0  0 \n" + 
					" 21 22  1  0  0  0  0 \n" + 
					" 25 26  2  0  0  0  0 \n" + 
					" 27 29  1  0  0  0  0 \n" + 
					"  8 11  1  0  0  0  0 \n" + 
					"  4 35  1  0  0  0  0 \n" + 
					" 35 36  3  0  0  0  0 \n" + 
					" 36  5  1  0  0  0  0 \n" + 
					"M  END",Charset.defaultCharset()).build();
			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		//ringSometimesExtendsToOCR.png
		list.add(new Object[]{"ringSometimesExtendsToOCR", new TestSpec("moleculeTest/ringSometimesExtendsToOCR.png", c->{
			
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"  CDK     01271921203D\n" + 
					"\n" + 
					" 42 47  0  0  0  0  0  0  0  0999 V2000\n" + 
					"  117.0000 -256.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  527.0000 -502.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  469.0000 -602.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  355.0000 -135.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  296.0000 -166.7169    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   80.3273 -201.4956    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  115.5000 -148.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  240.0000 -135.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  180.0000 -166.9106    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  468.3914 -268.2704    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  411.0000 -235.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  527.8240 -567.8992    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  409.5016 -568.7134    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  533.0000 -118.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  580.0000 -164.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  180.2809 -234.1626    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  726.7500  -66.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  662.0000  -83.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  459.3488 -120.2028    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  420.5000  -88.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  411.0000 -368.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  469.0000 -335.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  615.0000  -36.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  527.0000 -368.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  644.0000 -148.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  550.2000  -52.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  583.5000 -602.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  525.5000 -235.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  411.0000 -168.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    9.5000 -202.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  296.0000 -501.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  353.0000 -468.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  411.7532 -500.1417    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  237.0000 -267.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  296.0000 -234.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  353.0000 -601.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  296.0000 -567.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  772.0000 -111.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  411.0000 -434.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  527.0000 -434.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  744.0000   -2.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  468.3854 -467.6397    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1 16  1  0  0  0  0 \n" + 
					" 32 33  1  0  0  0  0 \n" + 
					" 34 35  2  0  0  0  0 \n" + 
					"  6 30  2  0  0  0  0 \n" + 
					" 36 37  2  0  0  0  0 \n" + 
					" 12 27  2  0  0  0  0 \n" + 
					"  4 29  1  0  0  0  0 \n" + 
					" 13 36  1  0  0  0  0 \n" + 
					" 14 26  1  0  0  0  0 \n" + 
					" 17 18  1  0  0  0  0 \n" + 
					" 19 20  2  0  0  0  0 \n" + 
					" 10 28  2  0  0  0  0 \n" + 
					" 21 39  1  0  0  0  0 \n" + 
					"  9 16  2  0  0  0  0 \n" + 
					" 19 29  1  6  0  0  0 \n" + 
					"  3 12  1  0  0  0  0 \n" + 
					" 15 25  1  0  0  0  0 \n" + 
					"  5 35  1  0  0  0  0 \n" + 
					"  3 13  1  0  0  0  0 \n" + 
					" 16 34  1  0  0  0  0 \n" + 
					" 22 24  1  0  0  0  0 \n" + 
					" 17 41  1  0  0  0  0 \n" + 
					" 40 42  1  0  0  0  0 \n" + 
					"  2 42  1  1  0  0  0 \n" + 
					" 10 22  1  0  0  0  0 \n" + 
					" 31 37  1  0  0  0  0 \n" + 
					" 31 32  2  0  0  0  0 \n" + 
					"  1  6  1  0  0  0  0 \n" + 
					"  4  5  1  0  0  0  0 \n" + 
					" 14 15  1  0  0  0  0 \n" + 
					"  5  8  2  0  0  0  0 \n" + 
					"  6  7  1  0  0  0  0 \n" + 
					" 11 29  1  0  0  0  0 \n" + 
					" 14 19  1  0  0  0  0 \n" + 
					"  7  9  1  0  0  0  0 \n" + 
					"  8  9  1  0  0  0  0 \n" + 
					" 17 38  1  0  0  0  0 \n" + 
					" 18 23  1  0  0  0  0 \n" + 
					" 13 33  2  0  0  0  0 \n" + 
					" 18 25  1  0  0  0  0 \n" + 
					" 21 22  1  0  0  0  0 \n" + 
					"  2 12  1  0  0  0  0 \n" + 
					" 23 26  1  0  0  0  0 \n" + 
					" 33 42  1  0  0  0  0 \n" + 
					" 24 40  1  0  0  0  0 \n" + 
					" 10 11  1  0  0  0  0 \n" + 
					" 39 42  1  0  0  0  0 \n" + 
					"M  END",Charset.defaultCharset()).build();
			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		//symbolsNMeVeryCloseTogether.png
		list.add(new Object[]{"NMeVeryCloseTogether", new TestSpec("moleculeTest/symbolsNMeVeryCloseTogether.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"  MJ150420                      \n" + 
					"\n" + 
					" 27 29  0  0  0  0  0  0  0  0999 V2000\n" + 
					"   -3.2678   -2.9322    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.0658    0.7611    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.6792   -0.8219    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.8048    3.3439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.5191    2.9497    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.7082   -1.6403    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.3546   -2.1461    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.0519   -2.9029    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.9566    0.4875    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.2521    0.8484    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.9758    2.0962    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.2505    1.6847    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.4045    0.4367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.4232   -0.3992    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.5231   -3.6559    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.0195   -2.1211    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.6760    3.3412    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.4045    2.9322    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.5222    2.1087    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.9680   -0.3858    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.6760    0.8457    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.1282    1.6844    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.9696    2.9322    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.4045    2.1087    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.9758   -1.2352    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.2459    3.3439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.7054   -3.5358    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1 27  1  0  0  0  0\n" + 
					"  9 10  1  1  0  0  0\n" + 
					" 11 12  2  0  0  0  0\n" + 
					" 13 14  1  0  0  0  0\n" + 
					" 17 18  2  0  0  0  0\n" + 
					"  2 22  1  0  0  0  0\n" + 
					"  9 20  1  0  0  0  0\n" + 
					"  3 14  1  0  0  0  0\n" + 
					" 13 21  1  0  0  0  0\n" + 
					" 17 23  1  0  0  0  0\n" + 
					" 22 24  1  0  0  0  0\n" + 
					"  1 16  1  0  0  0  0\n" + 
					"  9 21  1  0  0  0  0\n" + 
					" 10 12  1  0  0  0  0\n" + 
					"  1  8  1  0  0  0  0\n" + 
					"  4  5  1  0  0  0  0\n" + 
					"  3  6  1  0  0  0  0\n" + 
					"  8 15  2  0  0  0  0\n" + 
					" 12 19  1  0  0  0  0\n" + 
					"  6  7  1  0  0  0  0\n" + 
					"  7  8  1  0  0  0  0\n" + 
					"  3 20  1  0  0  0  0\n" + 
					"  6 16  2  0  0  0  0\n" + 
					"  3 25  1  1  0  0  0\n" + 
					"  5 26  1  0  0  0  0\n" + 
					" 13  2  1  6  0  0  0\n" + 
					" 18 24  1  0  0  0  0\n" + 
					" 23 26  1  0  0  0  0\n" + 
					"  5 19  2  0  0  0  0\n" + 
					"M  END\n" + 
					"",Charset.defaultCharset()).build();
			
			Chemical c1=c.connectedComponentsAsStream()
					.map(ct->Tuple.of(ct,ct.getAtomCount()).withVComparator())
					.max(Comparator.naturalOrder())
					.map(t->t.k())
					.orElse(c);
			String form=c1.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		//PMBNShorthand.png
		list.add(new Object[]{"PMBNShorthand", new TestSpec("moleculeTest/PMBNShorthand.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("[H][C@]1(CSC(=O)N1Cc2ccc(OC)cc2)[C@]3(C[C@@H](C[C@@H](CC)O3)OC(=O)c4ccccc4)OC").build();
			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		//
		list.add(new Object[]{"dashedMethylWithNoLabel", new TestSpec("moleculeTest/dashedMethylNoLabel.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("COc1ccc(cc1N2CCCC2)C(=O)N[C@@H](C)C(=O)N(Cc3cccs3)Cc4ccccc4").build();
			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		//terminalGroupCloseToOtherNode.png
		list.add(new Object[]{"terminalGroupCoseToOtherNode", new TestSpec("moleculeTest/terminalGroupCloseToOtherNode.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("COc1ccc(Sc2ccccc2C(C)C)cc1").build();
			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		
		list.add(new Object[]{"ringSystemWithHInMiddle", new TestSpec("moleculeTest/ringSystemWithHInCenter.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("[H][C@@]12CCCC(NC(=O)c3nccc(OC)c3O)[C@]1([H])CC=C(CC\\C=C(/C)C)C2").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		
		list.add(new Object[]{"complexStructure1Test", new TestSpec("moleculeTest/complex.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("[H][C@@]12CN(C[C@]1([H])[C@H]2NCc3ccc4cc(F)ccc4n3)c5ncc(cn5)C(=O)NO").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"fuzzyStructure1Test", new TestSpec("moleculeTest/fuzzy1.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("ClC(=O)c1ccc(Oc2ccc(cc2)C(Cl)=O)cc1").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"fuzzyStructure2Test", new TestSpec("moleculeTest/fuzzy2.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("CC(=C)C(=O)OCCOC(=O)c1ccc(C(=O)OCCCOC(=O)C=C)c(c1)C(=O)OCC(O)COc2ccc(Cc3ccc(OCC4CO4)cc3)cc2").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		
		list.add(new Object[]{"slimmerSulferTest", new TestSpec("moleculeTest/slimmerSulfur.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"  CDK     01251918013D\n" + 
					"\n" + 
					" 42 45  0  0  0  0  0  0  0  0999 V2000\n" + 
					"   98.4258 -543.2349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  156.7500 -574.5789    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  271.5000 -576.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  213.0000 -543.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  486.0000 -212.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  543.2196 -178.1228    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  715.5000  -12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  715.9093  -77.2111    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  388.4502 -442.1299    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  717.0000 -212.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  774.0000 -179.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  445.0000 -676.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  387.5000 -643.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  609.5000 -379.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  542.5000 -379.0000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  475.5000 -379.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  877.0000 -146.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  839.6644 -193.7010    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  504.0000 -443.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  659.0000 -112.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  659.0000 -178.7108    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   42.0000 -576.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  837.5939  -92.4048    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  446.0000 -474.7289    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   99.0000 -676.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  157.0000 -641.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  485.0000  -12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  486.0000  -77.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  870.0000  -34.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  542.5000 -112.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  542.2500 -311.1518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  388.7500 -574.8578    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  331.2292 -542.1304    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  773.2500 -111.1538    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  968.0000 -332.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  934.0000 -273.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  486.0000 -277.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  331.0000 -476.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  601.0000 -211.2965    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  601.0000 -277.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  446.0000 -541.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  872.0000 -257.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  2 26  1  0  0  0  0 \n" + 
					" 32 33  1  0  0  0  0 \n" + 
					" 15 19  1  0  0  0  0 \n" + 
					"  6 30  1  0  0  0  0 \n" + 
					" 11 34  1  0  0  0  0 \n" + 
					" 15 16  2  0  0  0  0 \n" + 
					" 17 18  2  0  0  0  0 \n" + 
					" 11 18  1  0  0  0  0 \n" + 
					"  9 38  1  0  0  0  0 \n" + 
					" 21 39  1  0  0  0  0 \n" + 
					" 36 42  1  0  0  0  0 \n" + 
					" 19 24  1  0  0  0  0 \n" + 
					"  3 33  1  0  0  0  0 \n" + 
					" 17 23  1  0  0  0  0 \n" + 
					"  5 37  2  0  0  0  0 \n" + 
					" 28 30  1  0  0  0  0 \n" + 
					" 32 41  2  0  0  0  0 \n" + 
					" 10 21  2  0  0  0  0 \n" + 
					" 31 37  1  0  0  0  0 \n" + 
					"  1  2  1  0  0  0  0 \n" + 
					" 33 38  2  0  0  0  0 \n" + 
					" 24 41  1  0  0  0  0 \n" + 
					"  2  4  1  0  0  0  0 \n" + 
					"  3  4  1  0  0  0  0 \n" + 
					" 23 34  2  0  0  0  0 \n" + 
					" 35 36  1  0  0  0  0 \n" + 
					" 39 40  1  0  0  0  0 \n" + 
					"  5  6  1  0  0  0  0 \n" + 
					" 12 13  1  0  0  0  0 \n" + 
					"  8 20  1  0  0  0  0 \n" + 
					" 14 15  2  0  0  0  0 \n" + 
					" 31 40  2  0  0  0  0 \n" + 
					"  7  8  2  0  0  0  0 \n" + 
					" 15 31  1  0  0  0  0 \n" + 
					"  9 24  2  0  0  0  0 \n" + 
					"  8 34  1  0  0  0  0 \n" + 
					" 13 32  1  0  0  0  0 \n" + 
					"  6 39  2  0  0  0  0 \n" + 
					" 20 21  1  0  0  0  0 \n" + 
					"  1 22  1  0  0  0  0 \n" + 
					" 25 26  1  0  0  0  0 \n" + 
					" 27 28  1  0  0  0  0 \n" + 
					" 18 42  1  0  0  0  0 \n" + 
					" 23 29  1  0  0  0  0 \n" + 
					" 10 11  1  0  0  0  0 \n" + 
					"M  END", Charset.defaultCharset()).build();
			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
//
		list.add(new Object[]{"nHConnectedTogetherTest", new TestSpec("moleculeTest/NHConnectedTogether.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("[H]C(O)(C(O)=O)C([H])(O)C(O)=O.Brc1c(ccc2nccnc12)\\N=C3\\NCCN3").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		list.add(new Object[]{"iodineContainingTest", new TestSpec("moleculeTest/iodineContaining.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("Ic1ccccc1C2CCCCC2").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		//closeRings.png
		list.add(new Object[]{"closeRingsNotJoinedTest", new TestSpec("moleculeTest/closeRings.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("O=C(NC(Cc1ccc(OCc2ccccc2)nc1)C(=O)N3CCC(CC3)N4CCCCC4)N5CCC(CC5)N6Cc7ccccc7NC6=O").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		
		list.add(new Object[]{"circleAromaticTest", new TestSpec("moleculeTest/circleAromatic.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("CC1CCCc2c(O)ccc(O)c12").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		
		list.add(new Object[]{"problematicIntersectionTest", new TestSpec("moleculeTest/problematicIntersection.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"  CDK     01261912223D\n" + 
					"\n" + 
					" 37 40  0  0  0  0  0  0  0  0999 V2000\n" + 
					"   10.0000 -174.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  330.0000 -227.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  511.0000 -207.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  451.2500 -238.5805    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  624.5000  -13.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  625.5000  -80.0000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  624.5000 -146.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  559.0000  -79.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  330.0000 -120.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  294.4027 -174.0280    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  509.0000 -120.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  393.2947 -206.1649    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  681.0000 -307.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  623.5000 -340.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  393.2500 -139.6275    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  567.0000 -306.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  692.8546  -77.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  190.0000 -230.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  125.0000 -230.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  451.0000 -305.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  725.5818  -21.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  792.0000  -21.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  726.0000 -136.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  792.0000 -136.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  826.0746  -78.1250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  567.0000 -240.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  509.0000 -339.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  190.0000 -115.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  224.9107 -172.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  124.6000 -115.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   91.4322 -171.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  443.0000  -98.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					" 1090.0000 -194.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					" 1026.0000 -194.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  891.0000  -78.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  992.0000 -136.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  926.0000 -136.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					" 34 36  1  0  0  0  0 \n" + 
					" 36 37  1  0  0  0  0 \n" + 
					" 13 14  1  0  0  0  0 \n" + 
					" 30 31  2  0  0  0  0 \n" + 
					"  9 15  1  0  0  0  0 \n" + 
					" 17 21  1  0  0  0  0 \n" + 
					" 17 23  2  0  0  0  0 \n" + 
					" 10 29  1  0  0  0  0 \n" + 
					" 24 25  2  0  0  0  0 \n" + 
					" 28 30  1  0  0  0  0 \n" + 
					" 22 25  1  0  0  0  0 \n" + 
					"  1 31  1  0  0  0  0 \n" + 
					" 28 29  2  0  0  0  0 \n" + 
					" 33 34  1  0  0  0  0 \n" + 
					" 35 37  1  0  0  0  0 \n" + 
					"  3  4  1  0  0  0  0 \n" + 
					" 25 35  1  0  0  0  0 \n" + 
					" 12 15  2  0  0  0  0 \n" + 
					"  3 26  2  0  0  0  0 \n" + 
					" 14 16  1  0  0  0  0 \n" + 
					"  5  6  2  0  0  0  0 \n" + 
					"  2 10  1  0  0  0  0 \n" + 
					" 18 19  2  0  0  0  0 \n" + 
					"  6  7  2  0  0  0  0 \n" + 
					"  4 20  2  0  0  0  0 \n" + 
					"  6  8  1  0  0  0  0 \n" + 
					" 15 32  1  0  0  0  0 \n" + 
					" 16 26  1  0  0  0  0 \n" + 
					"  9 10  2  0  0  0  0 \n" + 
					" 18 29  1  0  0  0  0 \n" + 
					" 16 27  2  0  0  0  0 \n" + 
					"  4 12  1  0  0  0  0 \n" + 
					"  6 17  1  0  0  0  0 \n" + 
					" 20 27  1  0  0  0  0 \n" + 
					" 11 32  1  0  0  0  0 \n" + 
					" 21 22  2  0  0  0  0 \n" + 
					" 19 31  1  0  0  0  0 \n" + 
					" 23 24  1  0  0  0  0 \n" + 
					"  2 12  1  0  0  0  0 \n" + 
					"  8 11  1  0  0  0  0 \n" + 
					"M  END", Charset.defaultCharset()).build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		list.add(new Object[]{"zerosForOxygensAndSmallInnerBondTest", new TestSpec("moleculeTest/withZerosAsOxygens.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("CC1C2C=CC1C(C2C(=O)OCC(C)=C)C(=O)OCC(C)=C").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"subscriptImplicitAtomsCl3Test", new TestSpec("moleculeTest/withSubscriptForCl.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("ClC(Cl)(Cl)c1nc(nc(n1)C(Cl)(Cl)Cl)-c2ccc3OCOc3c2").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"subscriptImplicitAtomsH2Test", new TestSpec("moleculeTest/withSubscriptForH.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("COc1cccc(C(O)c2cc(F)ccc2-N)c1C").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"moleculeWithCloseNitrogensInRingTest", new TestSpec("moleculeTest/moleculeWithCloseNitrogensInRing.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("FC(F)(F)CNc1nc(Nc2ccc(cc2)N3CCOCC3)nc4ccsc14").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"moleculeWith2CarboxyShortHandsTest", new TestSpec("moleculeTest/carboxylicShorthandNotation.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("OC(=O)Cc1ccc(OCc2ccccc2C(O)=O)cc1").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"nCarbonChainTest", new TestSpec("moleculeTest/carbonChainShorthand.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("CCCc1ccc(CCC)c2cc3c(-c4ccccc4)c5cc6c(CCC)ccc(CCC)c6cc5c(-c7ccccc7)c3cc12").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		//CCOC(=O)C[C@H](NS(=O)(=O)c1ccccc1)c2ccc3N(CC)C(C)Cc3c2
		
		list.add(new Object[]{"nCarbonChainTest", new TestSpec("moleculeTest/NHOnTopOfEachOther.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("CCOC(=O)C[C@H](NS(=O)(=O)c1ccccc1)c2ccc3N(CC)C(C)Cc3c2").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"serifFontTest", new TestSpec("moleculeTest/serifFont.png", c->{

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


			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"serifFont2Test", new TestSpec("moleculeTest/serifFont2.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("C=CC(C(c1ccc(OCCN2CCOCC2)cc1)c3ccc(OCCN4CCOCC4)cc3)c5ccccc5").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		list.add(new Object[]{"structureWithWeirdAngleBondToOtherOCRAtomTest", new TestSpec("moleculeTest/weirdAngleBetweenNitrogens.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("c1c(nn(c1-c2ccccc2)-c3cccc(c3)-n4c5ccc(cc5c6cc(ccc46)-c7ccccc7)-c8ccccc8)-c9ccccc9").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWhichIsAlmostARingAndHas2NonBondedOxygensVeryCloseTogether", new TestSpec("moleculeTest/almostRing.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("CC(=C)C(=O)OCCOCC(COC(=O)CCC(=O)OCCCCOC(=O)C=C)OC(=O)c1ccccc1C(=O)OCC(O)COc2ccc(cc2)C(C)(C)c3ccc(OCC4CO4)cc3").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithTightBondToNitrogens", new TestSpec("moleculeTest/tightBondsToNitrogens.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("C-c1cc(CN2C(=O)C3=C(CCCC3)C2=O)c(O)c(c1)-n4nc5ccccc5n4").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithLongLookingBondInAromaticRing", new TestSpec("moleculeTest/longLookingBondInAromaticRing.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("C(\\C=C\\c1ccc(cc1)N(c2ccccc2)c3ccc(\\C=N\\N(c4ccccc4)c5ccccc5)cc3)=C/c6ccccc6").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithAromaticSystemWhichSometimesChoosesWrongDoubleBond", new TestSpec("moleculeTest/aromaticSystemSometimesWrongDoubleBondChosen.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1").build();

			Chemical c1=c.connectedComponentsAsStream()
					.map(ct->Tuple.of(ct,ct.getAtomCount()).withVComparator())
					.max(Comparator.naturalOrder())
					.map(t->t.k())
					.orElse(c);

			String form=c1.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithSmallLineForCl", new TestSpec("moleculeTest/smallLineCl.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("[H]n1c(Cl)nc2n(CCC3CC3)c(=O)n([H])c(=O)c12").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithOxygenOffCenterInRing", new TestSpec("moleculeTest/connectedOxygen.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("COc1cc2CN(C)c3c(ccc4cc5OCOc5cc34)-c2cc1OC").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithOxygenConnectedToBonds", new TestSpec("moleculeTest/connectedOxygen2.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("Oc1ccc-2c(Cc3c-2c4Cc5cc(O)ccc5-c4c6Cc7cc(O)ccc7-c36)c1").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithOxygenConnectedToBonds2", new TestSpec("moleculeTest/connectedOxygen3.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("O=C(Oc1ccc(OC(=O)C2CCC3C(C2)C(=O)OC3=O)cc1)C4CCC5C(C4)C(=O)OC5=O").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithExplicitCarbons", new TestSpec("moleculeTest/explicitCarbons.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("CC(=O)CC(C)=O").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithCarbonsThatAreSometimesMistakenForOxygens", new TestSpec("moleculeTest/carbonVsOxygen.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("ClC(Cl)(Cl)c1nc(nc(n1)C(Cl)(Cl)Cl)-c2ccc3OCOc3c2").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithBridgeHeadInsideRing", new TestSpec("moleculeTest/bridgeHeadMolecule.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("CC1C2C3OC3C1C(C2C(=O)OCC=C)C(=O)OCC(C)=C").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithOxygensThatAreSometimesMistakenForCarbons", new TestSpec("moleculeTest/carbonVsOxygen2.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("CCOC(=O)CC1CCN(CC1)C(=O)CC2OC(c3cc(Cl)ccc3-n4cccc24)c5cccc6ccccc56").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithVeryCloseExplicitLinearAtoms", new TestSpec("moleculeTest/closeOCRShapesVeryExplicit.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("CCC(C)(C)COC(=O)c1ccc(cc1)C(=O)OC").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithDoubelBondAngleSimilarToN", new TestSpec("moleculeTest/doubleBondSometimesSeenAsNAtom.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("OC(c1cccc2OCCOc12)c3cc(Cl)ccc3-n4cccc4\\C=C\\C(=O)OCc5ccccc5").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithCloseNitrogens", new TestSpec("moleculeTest/nitrogensSometimesBondedByMistake.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("OCCCCCNc1ncccc1C(=O)Nc2ccc(cc2)-n3nc(cc3C(F)(F)F)-c4cccnc4").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithVeryLongErroneousLineAndNoisyLine", new TestSpec("moleculeTest/longBadLineNoisyLine.png", c1->{

			Chemical cReal=ChemicalBuilder.createFromSmiles("O=C(Oc1ccc(OC(=O)c2ccc3C(=O)OC(=O)c3c2)cc1)c4ccc5C(=O)OC(=O)c5c4").build();

			Chemical c=c1.connectedComponentsAsStream()
					.map(ct->Tuple.of(ct,ct.getAtomCount()).withVComparator())
					.max(Comparator.naturalOrder())
					.map(t->t.k())
					.orElse(c1);

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithVeryCloseAromaticRings", new TestSpec("moleculeTest/veryCloseAromaticRings.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("COc1ccc(cc1)C2(c3ccccc3-c4ccccc24)c5ccc(OC)cc5").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithVeryShortSingleBondBetweenCarbons", new TestSpec("moleculeTest/verySmallSingleBondBetweenExplicitCarbons.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("C(C=Cc1ccc(cc1)N(c2ccccc2)c3ccc(cc3)-c4ccc(cc4)N(c5ccccc5)c6ccc(C=CC=Cc7ccccc7)cc6)=Cc8ccccc8").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"structureWithGapInSmallRing", new TestSpec("moleculeTest/moleculeWithGapInSmallRing.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("CC(=C)C(=O)OCCOC(=O)c1ccc(C(=O)OCCOC(=O)C=C)c(c1)C(=O)OCC(O)COc2ccc(cc2)C(C)(C)c3ccc(OCC4CO4)cc3").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"aromaticRingSystemSometimesDoubleCounted", new TestSpec("moleculeTest/ringSystemProblem.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("c1ccc(cc1)-c2c3c4ccc5c6cccc7cccc(c8ccc(c3c(-c9ccccc9)c%10ccccc2%10)c4c58)c67").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		list.add(new Object[]{"alphaChannel", new TestSpec("moleculeTest/alphaChannel.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("CCCc1ccc(CCC)c2cc3c(-c4ccccc4)c5cc6c(CCC)ccc(CCC)c6cc5c(-c7ccccc7)c3cc12").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		//This one needs work, it's an outlier
		/*
		list.add(new Object[]{"subscriptImplicitAtomsF3Test", new TestSpec("moleculeTest/withSubscriptForF.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("FC(F)(F)C1(N=N1)c2ccc(CN3C(=O)C=CC3=O)cc2").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		*/
		return list;
	}


}
