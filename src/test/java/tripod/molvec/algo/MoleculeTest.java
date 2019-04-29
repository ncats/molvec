package tripod.molvec.algo;

import static org.junit.Assert.assertEquals;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.*;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import gov.nih.ncats.chemkit.api.io.ChemFormat;
import gov.nih.ncats.chemkit.api.util.stream.ThrowingStream;
import gov.nih.ncats.chemkit.renderer.ChemicalRenderer;
import gov.nih.ncats.chemkit.renderer.RendererOptions;
import gov.nih.ncats.chemkit.renderer.RendererOptions.DrawOptions;
import org.junit.*;
import tripod.molvec.Molvec;

import org.junit.rules.TemporaryFolder;

import gov.nih.ncats.chemkit.api.Chemical;
import gov.nih.ncats.chemkit.api.ChemicalBuilder;
import gov.nih.ncats.chemkit.api.inchi.Inchi;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import javax.imageio.ImageIO;

@RunWith(Parameterized.class)
public class MoleculeTest {

//	static File writeToFolder = new File("testResults");

	static File writeToFolder = null;
	@BeforeClass
	public static void emptyFolder(){
		if(writeToFolder !=null){
			File[] fs= writeToFolder.listFiles();
			if(fs !=null){
				for(File f : fs){
					f.delete();
				}
			}
		}
	}
	@Rule
	public TemporaryFolder tmpDir = new TemporaryFolder();
	
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

	private static final 	ChemFormat.SmilesFormatWriterSpecification smilesSpec = new ChemFormat.SmilesFormatWriterSpecification()
			.setCanonization(ChemFormat.SmilesFormatWriterSpecification.CanonicalizationEncoding.CANONICAL)
			.setKekulization(ChemFormat.KekulizationEncoding.KEKULE);

	ChemicalRenderer renderer = new ChemicalRenderer(RendererOptions.createINNLike()
			.setDrawOption(DrawOptions.DRAW_TERMINAL_CARBON, true)
			.setDrawOption(DrawOptions.DRAW_CARBON, true)
			.setDrawOption(DrawOptions.DRAW_STEREO_LABELS, false)
			.setDrawOption(RendererOptions.DrawOptions.DRAW_GREYSCALE, true)
					).setBackgroundColor(Color.white).setShadowVisible(false);


	private TestSpec spec;

	public MoleculeTest(String ignored, TestSpec spec){
		this.spec = spec;
	}


	@Test
	public void testAsFile() throws Exception {
		File f=getFile(spec.filePath);

		StructureImageExtractor sie = new StructureImageExtractor(f);

		//Chemical c =Chemical.parseMol(sie.getCtab().toMol());
		 Chemical c =Chemical.parseMol(sie.getCtab().toMol()).toBuilder().aromatize(false).build();
		spec.assertionConsumer.accept(c);
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
		String s = new StructureImageExtractor(out.toByteArray()).getCtab().toMol();

		StructureImageExtractor sie = new StructureImageExtractor(out.toByteArray());

		Chemical c = Chemical.parseMol(sie.getCtab().toMol()).toBuilder().aromatize(false).build();

		if(writeToFolder !=null){
			writeToFolder.mkdirs();
			File molvec = new File(writeToFolder, f.getName()+".molvec.png");
			BufferedImage img = renderer.createImage(c, 1000, 1000, false);

			ImageIO.write(img, "png", molvec);
			Files.copy(f.toPath(), new File(writeToFolder, f.getName()+".expected.png").toPath());
		}
		spec.assertionConsumer.accept(c);
	}

	@Test
	@Ignore
	public void rendererRoundTrip() throws Exception {

		File f = getFile(spec.filePath);
		StructureImageExtractor sie = new StructureImageExtractor(f);
		Chemical expected = Chemical.parseMol(sie.getCtab().toMol());

		expected.kekulize();

		BufferedImage img = renderer.createImage(expected, 1000, 1000, false);

		File newFile = tmpDir.newFile("molvec.png");

		ImageIO.write(img, "png", newFile);

		if(writeToFolder !=null){
			Files.copy(newFile.toPath(), new File(writeToFolder, f.getName()+".roundTrip.png").toPath());
		}

		spec.assertionConsumer.accept(Chemical.parseMol(Molvec.ocr(newFile)));
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

		//spiroCarbon.png
				list.add(new Object[]{"lotsOfRings", new TestSpec("moleculeTest/lotsOfRings.png", c->{
					Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
							"  Molvec0104291915382D\n" + 
							"\n" + 
							" 92119  0  0  0  0  0  0  0  0999 V2000\n" + 
							"  -10.0092    2.7588    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -9.1656    3.2224    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -8.2232    2.7512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    0.4180    2.7588    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.2920    3.2528    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    3.8760    1.7328    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    4.8488    1.2768    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"  -10.0092   -0.2508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -9.1352    0.2280    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -8.2232   -0.2584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -6.4904   -0.2584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -7.3720    0.2508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -4.7576   -0.2584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -5.6392    0.2584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -3.0096   -0.2584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -3.9064    0.2584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    0.4104   -0.2432    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.2920    0.2508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    2.1432   -0.2432    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    3.0552    0.2432    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -4.7272   -3.2528    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -5.6468   -2.7740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"  -10.0168   -3.2528    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -9.0744   -2.7816    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -8.1928   -3.2528    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -6.5512   -3.2528    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    6.5816    1.2768    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    5.6696    1.7328    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -5.6696    1.2616    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -4.7728    1.7480    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -3.0856   -3.2528    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -2.1736   -2.7664    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    8.3448    1.2768    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    7.4024    1.7328    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   10.0776    1.2768    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    9.1352    1.7480    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -3.9368    1.2616    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -3.0248    1.7480    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -2.2040    1.2616    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -1.2920    1.7480    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.4408    1.2455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -5.6696   -1.7480    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -4.7728   -1.2616    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -3.9368   -1.7480    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -3.0248   -1.2616    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   10.8376    1.7328    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -6.4980    1.7404    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -6.4980   -1.2692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -2.2027   -1.7919    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    3.0780    1.2388    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    0.4408    1.7708    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -1.2920   -1.2464    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"  -12.6008    0.2736    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"  -12.6008    1.2160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"  -10.8680    0.2508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"  -10.8452    1.2388    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -7.3948    1.2388    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.3148    1.2388    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"  -10.8680   -2.7360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"  -10.8680   -1.7708    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -7.3720   -2.7588    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -7.3492   -1.7708    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -7.3720    3.2528    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -5.6392    3.2528    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -6.4752    2.7588    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -3.9064    3.2528    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -4.7424    2.7588    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -2.1736    3.2528    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -3.0020    2.7588    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.4408    3.2528    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -1.2692    2.7588    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"  -11.7192    2.7360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"  -10.8832    3.2376    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -9.1276    1.2388    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -9.9636    1.7404    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -8.2308    1.7404    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"  -11.6964   -0.2508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -2.1736    0.2508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.4408    0.2508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -1.2692   -0.2508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -8.2536   -1.2692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -9.0820   -1.7708    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -3.9064   -2.7512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    2.1736    1.7556    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    2.1736    2.7360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    9.1352    2.7360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"  -11.7192   -1.2464    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -9.9864   -1.2692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   12.6008    1.7328    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   11.7496    1.2464    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    8.2840    3.2528    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"  -11.7268    1.7328    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							" 48 62  1  0\n" + 
							" 68 71  1  0\n" + 
							" 15 16  1  0\n" + 
							" 50 84  1  0\n" + 
							" 25 61  1  0\n" + 
							" 36 86  1  0\n" + 
							" 56 92  1  0\n" + 
							" 39 78  1  0\n" + 
							"  3 76  1  0\n" + 
							" 68 69  1  0\n" + 
							" 43 44  1  0\n" + 
							"  8 88  1  0\n" + 
							" 40 41  1  0\n" + 
							" 79 80  1  0\n" + 
							" 37 38  1  0\n" + 
							"  4  5  1  0\n" + 
							" 82 88  1  0\n" + 
							" 18 19  1  0\n" + 
							" 63 65  1  0\n" + 
							"  8  9  1  0\n" + 
							" 20 50  1  0\n" + 
							" 21 22  1  0\n" + 
							" 57 76  1  0\n" + 
							" 29 30  1  0\n" + 
							" 74 76  1  0\n" + 
							" 74 75  1  0\n" + 
							"  4 51  1  0\n" + 
							" 55 56  1  0\n" + 
							" 14 29  1  0\n" + 
							" 21 83  1  0\n" + 
							" 17 18  1  0\n" + 
							" 44 83  1  0\n" + 
							" 66 69  1  0\n" + 
							" 66 67  1  0\n" + 
							" 45 49  1  0\n" + 
							"  4 70  1  0\n" + 
							"  9 74  1  0\n" + 
							" 54 92  1  0\n" + 
							"  1 75  1  0\n" + 
							"  1 73  1  0\n" + 
							" 30 67  1  0\n" + 
							" 31 32  1  0\n" + 
							" 52 80  1  0\n" + 
							" 39 40  1  0\n" + 
							" 11 48  1  0\n" + 
							" 47 65  1  0\n" + 
							" 24 82  1  0\n" + 
							" 41 79  1  0\n" + 
							" 61 62  1  0\n" + 
							" 42 48  1  0\n" + 
							" 42 43  1  0\n" + 
							" 23 24  1  0\n" + 
							" 60 88  1  0\n" + 
							" 55 77  1  0\n" + 
							" 72 73  1  0\n" + 
							" 60 87  1  0\n" + 
							" 77 87  1  0\n" + 
							" 78 80  1  0\n" + 
							" 30 37  1  0\n" + 
							" 72 92  1  0\n" + 
							" 53 54  1  0\n" + 
							" 47 57  1  0\n" + 
							" 46 90  1  0\n" + 
							" 11 12  1  0\n" + 
							" 11 14  1  0\n" + 
							" 81 82  1  0\n" + 
							" 19 20  1  0\n" + 
							" 89 90  1  0\n" + 
							" 64 65  1  0\n" + 
							" 16 37  1  0\n" + 
							" 64 67  1  0\n" + 
							" 38 69  1  0\n" + 
							" 13 43  1  0\n" + 
							" 18 58  1  0\n" + 
							" 22 26  1  0\n" + 
							" 33 36  1  0\n" + 
							"  2  3  1  0\n" + 
							" 22 42  1  0\n" + 
							" 10 81  1  0\n" + 
							" 33 34  1  0\n" + 
							"  6  7  1  0\n" + 
							" 44 45  1  0\n" + 
							" 86 91  1  0\n" + 
							" 15 78  1  0\n" + 
							" 40 71  1  0\n" + 
							"  8 55  1  0\n" + 
							" 58 84  1  0\n" + 
							" 49 52  1  0\n" + 
							" 29 47  1  0\n" + 
							" 70 71  1  0\n" + 
							" 53 77  1  0\n" + 
							" 51 58  1  0\n" + 
							" 59 60  1  0\n" + 
							" 13 16  1  0\n" + 
							" 38 39  1  0\n" + 
							" 41 51  1  0\n" + 
							" 13 14  1  0\n" + 
							"  6 50  1  0\n" + 
							" 15 45  1  0\n" + 
							" 56 75  1  0\n" + 
							" 24 25  1  0\n" + 
							" 23 59  1  0\n" + 
							" 27 34  1  0\n" + 
							"  1  2  1  0\n" + 
							" 32 49  1  0\n" + 
							" 35 36  1  0\n" + 
							" 12 57  1  0\n" + 
							" 62 81  1  0\n" + 
							" 26 61  1  0\n" + 
							" 84 85  1  0\n" + 
							"  3 63  1  0\n" + 
							"  9 10  1  0\n" + 
							" 17 79  1  0\n" + 
							"  7 28  1  0\n" + 
							" 35 46  1  0\n" + 
							" 31 83  1  0\n" + 
							" 27 28  1  0\n" + 
							" 10 12  1  0\n" + 
							"  5 85  1  0\n" + 
							"M  END", Charset.defaultCharset()).build();

					String keyReal=Inchi.asStdInchi(cReal).getKey();
					String keyGot=Inchi.asStdInchi(c).getKey();
					assertEquals(keyReal,keyGot);
				} )});
		
		
		
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
		
		//closeFs.png
				list.add(new Object[]{"closeFs", new TestSpec("moleculeTest/closeFs.png", c->{
					Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
							"  CDK     02051920443D\n" + 
							"\n" + 
							" 40 44  0  0  0  0  0  0  0  0999 V2000\n" + 
							"   -2.5943   -2.7489    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -2.0941   -1.8917    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    0.0270    1.7396    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.8724    1.2412    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    2.5902   -2.7574    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    3.4671   -3.2408    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.7433   -0.2578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    0.8739    0.2368    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.7495   -1.2401    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.7282   -3.2769    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    0.8739   -2.7611    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -3.4339    2.7223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -4.3350    3.2348    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -3.4581   -3.2558    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    0.8560   -1.7446    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.0030   -3.2483    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -2.5887   -1.7568    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    0.8851    1.2524    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.0030   -1.2621    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -3.4656    1.7358    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.7208   -4.2601    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -2.5887    3.2460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -2.5887    4.2391    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.0067    2.7514    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    2.6240   -1.7456    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -1.7193   -4.2301    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -1.7081   -3.2445    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.7433    2.7551    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    0.8739    3.2498    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.7433    1.7658    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.8499    4.2541    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -1.7268    4.7413    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.8499    3.2460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -1.7193    2.7551    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    4.3365   -1.7418    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    3.4671   -1.2471    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    4.3365   -2.7311    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.8499   -2.7461    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    0.0045   -4.2301    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.8499   -4.7398    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"  1 17  1  0  0  0  0 \n" + 
							" 24 33  1  0  0  0  0 \n" + 
							"  1 14  1  0  0  0  0 \n" + 
							" 15 19  2  0  0  0  0 \n" + 
							" 11 15  1  0  0  0  0 \n" + 
							"  3 18  1  0  0  0  0 \n" + 
							" 18 30  2  0  0  0  0 \n" + 
							" 11 16  1  0  0  0  0 \n" + 
							"  9 15  1  0  0  0  0 \n" + 
							" 16 38  2  0  0  0  0 \n" + 
							" 22 23  2  0  0  0  0 \n" + 
							" 28 30  1  0  0  0  0 \n" + 
							" 26 27  2  0  0  0  0 \n" + 
							" 28 29  2  0  0  0  0 \n" + 
							" 24 29  1  0  0  0  0 \n" + 
							" 10 21  2  0  0  0  0 \n" + 
							"  1  2  1  0  0  0  0 \n" + 
							" 31 33  1  0  0  0  0 \n" + 
							"  1 27  1  0  0  0  0 \n" + 
							" 27 38  1  0  0  0  0 \n" + 
							" 23 32  1  0  0  0  0 \n" + 
							" 31 32  2  0  0  0  0 \n" + 
							" 25 36  1  0  0  0  0 \n" + 
							" 33 34  2  0  0  0  0 \n" + 
							" 35 37  1  0  0  0  0 \n" + 
							"  3  4  1  0  0  0  0 \n" + 
							" 35 36  2  0  0  0  0 \n" + 
							" 39 40  2  0  0  0  0 \n" + 
							" 12 13  1  0  0  0  0 \n" + 
							"  5  6  1  0  0  0  0 \n" + 
							" 16 39  1  0  0  0  0 \n" + 
							" 12 20  2  0  0  0  0 \n" + 
							"  5 10  1  0  0  0  0 \n" + 
							"  7  8  1  0  0  0  0 \n" + 
							"  7  9  2  0  0  0  0 \n" + 
							"  9 25  1  0  0  0  0 \n" + 
							"  8 18  1  0  0  0  0 \n" + 
							" 12 22  1  0  0  0  0 \n" + 
							"  6 37  2  0  0  0  0 \n" + 
							"  5 25  2  0  0  0  0 \n" + 
							"  3 24  2  0  0  0  0 \n" + 
							" 26 40  1  0  0  0  0 \n" + 
							" 22 34  1  0  0  0  0 \n" + 
							" 10 11  1  0  0  0  0 \n" + 
							"M  END",Charset.defaultCharset()).build();
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
			Chemical cReal= Chemical.parseMol(

					"\n" +
							"  CDK     01311923073D\n" +
							"\n" +
							" 55 62  0  0  0  0  0  0  0  0999 V2000\n" +
							"   -4.3270    0.4544    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -5.2877    0.2122    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    3.9929   -0.7819    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    4.5609   -0.0134    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    6.3485   -1.8010    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    7.2006   -1.2831    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    7.2006   -0.3308    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    6.3513    0.1738    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    4.0514   -4.0397    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    3.2912   -4.4908    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    4.5400   -2.5737    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    4.0848   -3.0289    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    2.5144   -3.0206    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    2.0424   -2.5444    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    2.4977   -4.0397    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    4.5609   -1.6172    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    5.4965   -1.3040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    1.1862   -0.1136    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    0.3508    0.3709    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    6.3485    1.1561    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -2.3724   -4.4406    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -1.6957   -3.7891    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    0.3398    1.3527    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    1.1778    1.8411    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -1.3532    3.3405    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -1.3699    4.3137    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    5.5090   -0.3016    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -5.5508   -0.7485    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -6.5240   -0.9991    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -1.3783    0.3709    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -1.3407    1.3691    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -2.2053    1.8536    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -3.0824    1.3566    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -7.2006   -0.2974    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -1.9422   -2.8076    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -0.5221    1.8536    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -0.5179    2.8435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -3.6588   -0.2180    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -2.4434   -0.9322    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -3.1409   -1.6005    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -1.2530   -2.1184    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -1.4869   -1.1828    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -2.2053    2.8435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -3.9261   -1.1828    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -4.8616   -1.4167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -2.8736   -2.5695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    2.0382   -1.5838    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    1.1862   -1.0826    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    2.8234   -5.3094    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    3.7590   -5.3094    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    3.2829   -2.3022    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -2.7029    0.0347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -2.2444    4.7987    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -0.5126    4.8285    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -1.3871    5.3135    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"  2 28  2  0  0  0  0 \n" +
							" 32 33  1  0  0  0  0 \n" +
							" 36 37  2  0  0  0  0 \n" +
							" 11 12  1  1  0  0  0 \n" +
							"  5 17  1  0  0  0  0 \n" +
							" 11 16  1  0  0  0  0 \n" +
							" 30 31  1  0  0  0  0 \n" +
							" 40 46  1  0  0  0  0 \n" +
							" 14 51  1  0  0  0  0 \n" +
							" 19 23  1  0  0  0  0 \n" +
							" 47 48  2  0  0  0  0 \n" +
							" 23 36  1  0  0  0  0 \n" +
							" 32 43  1  0  0  0  0 \n" +
							" 33 52  1  0  0  0  0 \n" +
							"  9 12  1  0  0  0  0 \n" +
							" 29 34  1  0  0  0  0 \n" +
							"  8 20  1  0  0  0  0 \n" +
							" 18 19  1  0  0  0  0 \n" +
							"  6  7  1  0  0  0  0 \n" +
							" 44 45  2  0  0  0  0 \n" +
							" 21 22  1  0  0  0  0 \n" +
							" 25 26  1  0  0  0  0 \n" +
							" 22 35  1  0  0  0  0 \n" +
							" 38 52  1  0  0  0  0 \n" +
							" 26 54  1  0  0  0  0 \n" +
							" 26 55  1  0  0  0  0 \n" +
							" 25 43  2  0  0  0  0 \n" +
							" 26 53  1  0  0  0  0 \n" +
							"  1 38  2  0  0  0  0 \n" +
							"  4 27  1  0  0  0  0 \n" +
							"  3 16  1  0  0  0  0 \n" +
							" 14 13  1  1  0  0  0 \n" +
							" 13 15  1  0  0  0  0 \n" +
							" 30 52  1  0  0  0  0 \n" +
							" 10 50  1  0  0  0  0 \n" +
							" 10 49  1  0  0  0  0 \n" +
							" 17 27  2  0  0  0  0 \n" +
							" 38 44  1  0  0  0  0 \n" +
							" 41 42  2  0  0  0  0 \n" +
							"  8 27  1  0  0  0  0 \n" +
							" 25 37  1  0  0  0  0 \n" +
							" 28 29  1  0  0  0  0 \n" +
							" 39 52  1  0  0  0  0 \n" +
							" 31 36  1  0  0  0  0 \n" +
							" 28 45  1  0  0  0  0 \n" +
							"  1  2  1  0  0  0  0 \n" +
							" 31 32  2  0  0  0  0 \n" +
							"  3  4  2  0  0  0  0 \n" +
							" 39 40  2  0  0  0  0 \n" +
							"  5  6  2  0  0  0  0 \n" +
							" 12 13  1  1  0  0  0 \n" +
							" 16 17  1  0  0  0  0 \n" +
							"  7  8  2  0  0  0  0 \n" +
							" 18 48  1  0  0  0  0 \n" +
							"  9 10  1  0  0  0  0 \n" +
							" 10 15  1  0  0  0  0 \n" +
							" 11 51  1  0  0  0  0 \n" +
							" 14 47  1  0  0  0  0 \n" +
							" 23 24  2  0  0  0  0 \n" +
							" 35 46  2  0  0  0  0 \n" +
							" 35 41  1  0  0  0  0 \n" +
							" 39 42  1  0  0  0  0 \n" +
							"M  END"
			);
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

			cReal.kekulize();
			c.kekulize();
			String form=c.getFormula();

			c.bonds().filter(b-> b.isInRing()).forEach(b-> System.out.println(b.getBondType()));

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

		//spiroCarbon.png
		list.add(new Object[]{"spiroCarbon", new TestSpec("moleculeTest/spiroCarbon.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" +
					"  CDK     02021921413D\n" +
					"\n" +
					" 53 64  0  0  0  0  0  0  0  0999 V2000\n" +
					"    0.7795    0.5709    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.7796    0.3937    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    3.9528    0.2520    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    3.4566    1.1325    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -3.9686   -0.2835    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -3.4606   -1.1391    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    4.9568    0.2598    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -4.9725   -0.2795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    2.4449    1.1378    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    2.4567   -1.1024    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.7953   -0.3937    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -1.1811   -2.2992    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -2.1260   -2.1103    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    4.9450    2.0158    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    3.9607    2.0079    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -5.4686    0.6102    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -6.4686    0.5866    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -0.4961   -1.5591    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -0.8071   -0.5984    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -2.1654    2.0394    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -1.1811    2.2520    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    0.4803    1.5433    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.1811   -2.2835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    2.1418   -2.0945    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -4.9607    1.4646    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -5.4647    2.3307    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -6.9607    1.4488    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -6.4647    2.3229    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    2.1260    2.0945    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.1575    2.2756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    5.4568    1.1575    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    6.4568    1.1496    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    6.9607    0.2677    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -6.4568   -1.1339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -6.9607   -0.2835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -5.4686   -1.1457    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    5.4647   -2.3307    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    4.9607   -1.4803    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    6.9607   -1.4646    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    6.4568   -2.3307    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    6.4529   -0.5866    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    5.4528   -0.6102    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -3.9686   -2.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -4.9607   -2.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -1.8111    0.3780    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -2.4725    1.0866    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -2.4607   -1.1339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -1.8032   -0.4173    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -0.8189    0.5630    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -0.5118    1.5355    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    0.8032   -0.5866    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    0.5000   -1.5709    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -0.0095   -0.0106    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					" 32 33  2  0  0  0  0 \n" +
					" 36 44  1  0  0  0  0 \n" +
					" 43 44  2  0  0  0  0 \n" +
					" 47 48  2  0  0  0  0 \n" +
					" 18 52  1  0  0  0  0 \n" +
					" 26 28  1  0  0  0  0 \n" +
					"  1 53  1  0  0  0  0 \n" +
					" 10 24  1  0  0  0  0 \n" +
					" 22 30  1  0  0  0  0 \n" +
					"  7 31  2  0  0  0  0 \n" +
					" 23 52  1  0  0  0  0 \n" +
					" 37 38  2  0  0  0  0 \n" +
					" 13 47  1  0  0  0  0 \n" +
					" 37 40  1  0  0  0  0 \n" +
					"  2  9  2  0  0  0  0 \n" +
					" 14 15  2  0  0  0  0 \n" +
					" 18 19  2  0  0  0  0 \n" +
					"  4  9  1  0  0  0  0 \n" +
					" 20 46  1  0  0  0  0 \n" +
					"  8 16  1  0  0  0  0 \n" +
					"  9 29  1  0  0  0  0 \n" +
					"  4 15  1  0  0  0  0 \n" +
					"  1 22  2  0  0  0  0 \n" +
					" 25 26  2  0  0  0  0 \n" +
					" 33 41  1  0  0  0  0 \n" +
					" 29 30  2  0  0  0  0 \n" +
					"  7 42  1  0  0  0  0 \n" +
					" 51 52  2  0  0  0  0 \n" +
					" 51 53  1  0  0  0  0 \n" +
					" 34 36  1  0  0  0  0 \n" +
					" 34 35  2  0  0  0  0 \n" +
					" 49 53  1  0  0  0  0 \n" +
					"  6 43  1  0  0  0  0 \n" +
					" 14 31  1  0  0  0  0 \n" +
					" 38 42  1  0  0  0  0 \n" +
					" 17 27  1  0  0  0  0 \n" +
					" 41 42  2  0  0  0  0 \n" +
					"  6 47  1  0  0  0  0 \n" +
					" 45 46  2  0  0  0  0 \n" +
					" 19 48  1  0  0  0  0 \n" +
					" 45 49  1  0  0  0  0 \n" +
					" 49 50  2  0  0  0  0 \n" +
					" 22 50  1  0  0  0  0 \n" +
					"  1  2  1  0  0  0  0 \n" +
					" 31 32  1  0  0  0  0 \n" +
					"  3  4  2  0  0  0  0 \n" +
					" 39 40  2  0  0  0  0 \n" +
					"  3  7  1  0  0  0  0 \n" +
					"  5  6  2  0  0  0  0 \n" +
					" 12 13  2  0  0  0  0 \n" +
					"  5  8  1  0  0  0  0 \n" +
					" 16 17  2  0  0  0  0 \n" +
					" 12 18  1  0  0  0  0 \n" +
					"  8 36  2  0  0  0  0 \n" +
					" 21 50  1  0  0  0  0 \n" +
					" 16 25  1  0  0  0  0 \n" +
					" 11 51  1  0  0  0  0 \n" +
					" 17 35  1  0  0  0  0 \n" +
					" 20 21  2  0  0  0  0 \n" +
					" 23 24  2  0  0  0  0 \n" +
					" 19 53  1  0  0  0  0 \n" +
					" 27 28  2  0  0  0  0 \n" +
					" 10 11  2  0  0  0  0 \n" +
					" 39 41  1  0  0  0  0 \n" +
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		//thinThickBondAndLongDashedLines.png
		list.add(new Object[]{"thinThickBondAndLongDashedLines", new TestSpec("moleculeTest/thinThickBondAndLongDashedLines.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"  CDK     02271911163D\n" + 
					"\n" + 
					" 34 38  0  0  0  0  0  0  0  0999 V2000\n" + 
					"    5.7962    1.9430    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    6.6422    2.4897    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    6.6645    0.4643    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.7823    0.9604    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -6.4986    0.7246    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.9750    1.5084    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.1180   -0.2346    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.9976    0.2378    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3488   -0.2346    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.2549   -0.6836    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.2355    0.5399    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.3978    0.9871    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.5551   -0.7108    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.9750   -0.0717    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.5100    0.4822    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.2286   -0.4789    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.0058    1.1989    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -7.5388    0.7265    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.9528    2.4470    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.0940    1.9694    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.0791    0.9851    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.9224    0.4811    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    7.5360    0.9795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    7.5244    1.9936    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.2610   -1.8635    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.2320    1.6713    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.2883   -0.1938    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.5571    0.6939    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.3158    0.2401    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.3035    1.1663    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.2737   -2.4825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.1859   -2.4825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.2576   -0.4300    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.3284   -1.1822    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1 19  1  0  0  0  0 \n" + 
					"  5 18  1  0  0  0  0 \n" + 
					" 11 12  1  0  0  0  0 \n" + 
					" 11 16  2  0  0  0  0 \n" + 
					" 19 20  2  0  0  0  0 \n" + 
					" 17 26  1  0  0  0  0 \n" + 
					"  9 13  1  0  0  0  0 \n" + 
					"  4 22  1  0  0  0  0 \n" + 
					" 10 25  1  0  0  0  0 \n" + 
					"  5 14  1  0  0  0  0 \n" + 
					" 11 21  1  0  0  0  0 \n" + 
					"  2 24  1  0  0  0  0 \n" + 
					" 26 30  1  0  0  0  0 \n" + 
					" 27 33  1  0  0  0  0 \n" + 
					" 25 32  1  0  0  0  0 \n" + 
					"  1  2  2  0  0  0  0 \n" + 
					"  1  4  1  0  0  0  0 \n" + 
					" 31 32  1  0  0  0  0 \n" + 
					"  3  4  2  0  0  0  0 \n" + 
					" 15 12  1  1  0  0  0 \n" + 
					" 27 13  1  1  0  0  0 \n" + 
					"  5  6  2  0  0  0  0 \n" + 
					" 15 28  1  0  0  0  0 \n" + 
					"  7  8  1  0  0  0  0 \n" + 
					"  8 14  1  0  0  0  0 \n" + 
					" 15 33  1  0  0  0  0 \n" + 
					"  9 10  1  0  0  0  0 \n" + 
					"  3 23  1  0  0  0  0 \n" + 
					"  7 29  1  0  0  0  0 \n" + 
					"  6 17  1  0  0  0  0 \n" + 
					" 21 22  2  0  0  0  0 \n" + 
					"  8 17  2  0  0  0  0 \n" + 
					" 23 24  2  0  0  0  0 \n" + 
					" 20 21  1  0  0  0  0 \n" + 
					" 27 28  1  0  0  0  0 \n" + 
					" 29 30  1  0  0  0  0 \n" + 
					" 29 10  1  6  0  0  0 \n" + 
					" 27 34  1  6  0  0  0 \n" + 
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		//moreDottedLinesThanNormalLines.png
		list.add(new Object[]{"moreDottedLinesThanNormalLines", new TestSpec("moleculeTest/moreDottedLinesThanNormalLines.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"\n" + 
					"\n" + 
					" 31 33  0  0  0  0  0  0  0  0999 V2000\n" + 
					"    0.1697    2.2164    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.6172   -0.2997    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.8846   -0.7227    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.0615   -0.3042    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.7899   -1.5601    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.0643   -1.9786    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.6150   -1.9729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.1266    0.5372    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.5480    1.2586    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.5637    0.9593    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.1623    0.5594    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.3373   -1.5586    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.3529   -0.7241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.3811    1.2712    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.7899    0.5302    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.8876    0.9439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.7899    1.9815    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.3769    2.7043    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.8819    1.8001    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.6131    2.2076    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.2893    2.2076    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.1590   -0.3042    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.5316   -0.1842    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.5637    1.8149    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.1973   -2.7043    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.2921    0.5330    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.6103    0.5330    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.0317   -2.7043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.1210   -0.9181    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.8875   -1.5517    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.3769   -0.1842    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1 19  1  0  0  0  0\n" + 
					"  2 27  1  6  0  0  0\n" + 
					" 17 18  1  0  0  0  0\n" + 
					" 11 16  1  0  0  0  0\n" + 
					" 19 20  1  1  0  0  0\n" + 
					"  7 12  1  0  0  0  0\n" + 
					"  8 26  1  6  0  0  0\n" + 
					" 10 26  1  6  0  0  0\n" + 
					"  8 23  1  0  0  0  0\n" + 
					"  9 14  1  0  0  0  0\n" + 
					" 11 22  1  1  0  0  0\n" + 
					" 10 24  1  0  0  0  0\n" + 
					" 23 31  1  0  0  0  0\n" + 
					"  2  3  1  0  0  0  0\n" + 
					"  1 24  1  0  0  0  0\n" + 
					" 14 17  1  0  0  0  0\n" + 
					"  7 30  1  0  0  0  0\n" + 
					" 12 13  1  0  0  0  0\n" + 
					"  5  6  1  0  0  0  0\n" + 
					" 14 15  2  0  0  0  0\n" + 
					"  3 30  1  0  0  0  0\n" + 
					" 16 19  1  0  0  0  0\n" + 
					"  8  9  1  0  0  0  0\n" + 
					"  7 25  1  1  0  0  0\n" + 
					" 15 31  1  0  0  0  0\n" + 
					" 16 27  1  6  0  0  0\n" + 
					" 12  6  1  1  0  0  0\n" + 
					" 13  4  1  6  0  0  0\n" + 
					" 24 21  1  1  0  0  0\n" + 
					"  7 28  1  6  0  0  0\n" + 
					"  2 13  1  0  0  0  0\n" + 
					" 23 29  1  6  0  0  0\n" + 
					" 10 11  1  0  0  0  0\n" + 
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		//HCloseToDoubleBond.png
		list.add(new Object[]{"HCloseToDoubleBond", new TestSpec("moleculeTest/HCloseToDoubleBond.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"\n" + 
					"\n" + 
					" 33 36  0  0  0  0  0  0  0  0999 V2000\n" + 
					"   -3.4685    1.2569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.3349    0.7583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.7301    1.3019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.8648    0.7594    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.0007    1.2566    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.8678    0.7732    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.4701   -0.7431    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.6027   -0.2429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.7340   -1.7507    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.8667   -2.2511    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.8859   -0.2583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.7470   -0.7349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.7367    1.3101    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.6032    0.7568    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.5947    0.7646    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.4653    1.3113    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.0193   -0.7332    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.7390   -0.7419    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.8473   -0.2592    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.8633    1.7520    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.8680   -1.2503    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.8667   -1.2503    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.3360    0.7512    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.7387   -1.7507    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.6013   -2.2511    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.0040    2.2524    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.3354   -0.2431    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.6013   -0.2445    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.4733    2.2524    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.4653   -0.7499    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.7387    2.2524    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.4653    2.2524    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.7307    2.2524    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  2 27  2  0  0  0  0 \n" + 
					"  1 14  2  0  0  0  0 \n" + 
					"  3 15  1  0  0  0  0 \n" + 
					" 11 12  1  0  0  0  0 \n" + 
					" 12 24  1  6  0  0  0 \n" + 
					" 13 14  1  0  0  0  0 \n" + 
					" 15 16  1  0  0  0  0 \n" + 
					" 11 17  1  0  0  0  0 \n" + 
					" 17 19  1  0  0  0  0 \n" + 
					"  3 33  2  0  0  0  0 \n" + 
					" 16 32  2  0  0  0  0 \n" + 
					" 11 21  1  1  0  0  0 \n" + 
					" 18  9  1  1  0  0  0 \n" + 
					" 19 22  1  1  0  0  0 \n" + 
					" 28 30  1  0  0  0  0 \n" + 
					"  1 29  1  0  0  0  0 \n" + 
					"  1  2  1  0  0  0  0 \n" + 
					"  3  4  1  0  0  0  0 \n" + 
					"  4  5  1  0  0  0  0 \n" + 
					" 15 28  2  0  0  0  0 \n" + 
					"  5  6  2  0  0  0  0 \n" + 
					"  4 19  1  0  0  0  0 \n" + 
					"  4 20  1  1  0  0  0 \n" + 
					" 18 19  1  0  0  0  0 \n" + 
					"  7  8  2  0  0  0  0 \n" + 
					"  6 11  1  0  0  0  0 \n" + 
					"  8 14  1  0  0  0  0 \n" + 
					" 13 31  2  0  0  0  0 \n" + 
					"  6 13  1  0  0  0  0 \n" + 
					"  9 10  1  0  0  0  0 \n" + 
					"  8 12  1  0  0  0  0 \n" + 
					" 18 28  1  0  0  0  0 \n" + 
					"  9 25  1  0  0  0  0 \n" + 
					"  7 27  1  0  0  0  0 \n" + 
					" 16 23  1  0  0  0  0 \n" + 
					"  5 26  1  0  0  0  0 \n" + 
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		list.add(new Object[]{"multiIodine", new TestSpec("moleculeTest/multiIodine.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"  CDK     02141919243D\n" + 
					"\n" + 
					" 21 21  0  0  0  0  0  0  0  0999 V2000\n" + 
					"   -0.4405    0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3195    0.1981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4462   -1.3216    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4318   -0.8187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.9975   -0.8513    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.9177   -1.3256    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3210   -0.8137    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.1744   -1.3256    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.3121   -1.3256    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.1758   -0.7622    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4945    1.7292    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3061    2.1945    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.0420   -1.3088    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.9170   -0.8227    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4342    0.1948    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.3054    0.6831    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4305    2.1945    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.1811    0.6936    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4378   -2.1938    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.1737    0.1830    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.0494    0.1830    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1  2  1  0  0  0  0 \n" + 
					"  2 18  1  0  0  0  0 \n" + 
					"  1 15  2  0  0  0  0 \n" + 
					"  3  4  2  0  0  0  0 \n" + 
					"  2  7  2  0  0  0  0 \n" + 
					" 11 12  2  0  0  0  0 \n" + 
					"  3  7  1  0  0  0  0 \n" + 
					"  5  6  1  0  0  0  0 \n" + 
					" 13 14  1  0  0  0  0 \n" + 
					" 15 16  1  0  0  0  0 \n" + 
					"  5  8  1  0  0  0  0 \n" + 
					" 11 17  1  0  0  0  0 \n" + 
					"  4  9  1  0  0  0  0 \n" + 
					"  3 19  1  0  0  0  0 \n" + 
					"  7  8  1  0  0  0  0 \n" + 
					"  5 21  2  0  0  0  0 \n" + 
					"  9 10  1  0  0  0  0 \n" + 
					" 10 20  2  0  0  0  0 \n" + 
					"  4 15  1  0  0  0  0 \n" + 
					"  1 11  1  0  0  0  0 \n" + 
					" 10 13  1  0  0  0  0 \n" + 
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});

		list.add(new Object[]{"noisySulfur", new TestSpec("moleculeTest/noisySulfur.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"\n" + 
					"\n" + 
					" 27 30  0  0  0  0  0  0  0  0999 V2000\n" + 
					"    4.2123   -0.2320    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.4460    0.1822    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.8508   -0.9695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.2721   -0.2200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.6901    0.4838    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.9912   -0.6948    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.4609    1.0670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.1530    1.4823    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.2700   -0.2412    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.9330   -1.0486    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.7296   -0.8806    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.3652   -0.2181    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.0104    0.5035    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.0401   -0.9695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.9044    1.0670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.2070   -1.0570    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.7807   -0.1885    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.9877    0.1673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.5304   -0.6716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.5205    0.1948    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.2286   -1.0961    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.2720    0.6221    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.9739    0.1970    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.9044    0.1970    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.8705    0.5134    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.9330   -1.4823    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.4616   -1.4645    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1 16  1  0  0  0  0\n" + 
					"  5 20  1  0  0  0  0\n" + 
					" 17 18  1  0  0  0  0\n" + 
					" 19 20  1  0  0  0  0\n" + 
					" 11 19  1  0  0  0  0\n" + 
					"  6 21  1  0  0  0  0\n" + 
					"  6 23  2  0  0  0  0\n" + 
					" 15 24  2  0  0  0  0\n" + 
					"  9 18  1  0  0  0  0\n" + 
					"  3 14  1  0  0  0  0\n" + 
					" 22 23  1  0  0  0  0\n" + 
					" 19 21  2  0  0  0  0\n" + 
					"  1  2  2  0  0  0  0\n" + 
					"  2 17  1  0  0  0  0\n" + 
					"  1 24  1  0  0  0  0\n" + 
					"  2  7  1  0  0  0  0\n" + 
					"  4  5  2  0  0  0  0\n" + 
					" 12 13  1  0  0  0  0\n" + 
					" 13 25  1  0  0  0  0\n" + 
					" 12 14  1  0  0  0  0\n" + 
					"  3  9  1  0  0  0  0\n" + 
					"  7  8  2  0  0  0  0\n" + 
					"  6 10  1  0  0  0  0\n" + 
					"  8 15  1  0  0  0  0\n" + 
					" 16 26  2  0  0  0  0\n" + 
					"  9 25  1  0  0  0  0\n" + 
					"  4 11  1  0  0  0  0\n" + 
					" 16 27  1  0  0  0  0\n" + 
					"  4 12  1  0  0  0  0\n" + 
					" 20 22  2  0  0  0  0\n" + 
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		//cagedStructure4.png
		list.add(new Object[]{"cagedStructure4", new TestSpec("moleculeTest/cagedStructure4.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"          \n" + 
					"\n" + 
					" 45 51  0  0  0  0  0  0  0  0999 V2000\n" + 
					"   -2.4884   -1.5794    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.7193    0.5561    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.5386    0.0037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.0623   -1.0158    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.8164   -1.4941    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.2427    0.5709    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.2788    0.2595    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.3899    0.5621    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.7786    2.3060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.5736    2.8867    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.1533    1.2086    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    6.1425    1.0603    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    6.6170    1.9353    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.6987   -1.0121    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.6772    1.2882    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.9057    1.9353    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4627    0.8972    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    6.0832    2.7657    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.2818   -0.9639    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.2754   -0.9565    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.1172    0.9713    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.5675   -0.4078    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -6.6185    1.6239    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.1162    2.2096    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.9434   -0.0556    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.9481   -0.0556    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.8275   -2.4840    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.0821    2.7138    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.1785    1.0158    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.6261    1.8315    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.9002   -1.4978    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.8927   -2.4840    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.0437   -2.9771    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.0060    2.5507    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.7232   -1.0158    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.9872    3.7742    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.0326   -3.9373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.4864    0.3230    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -6.2403   -0.3040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.6400    1.8129    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.8147    2.4766    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.6007    3.0401    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.0159    3.9373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.3061   -0.6154    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3845    1.5609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1 19  1  0  0  0  0 \n" + 
					" 30 40  1  0  0  0  0 \n" + 
					" 32 33  2  0  0  0  0 \n" + 
					"  1 14  1  0  0  0  0 \n" + 
					"  2 26  1  0  0  0  0 \n" + 
					"  3 38  1  0  0  0  0 \n" + 
					" 15 16  1  0  0  0  0 \n" + 
					" 19 20  1  0  0  0  0 \n" + 
					" 12 29  1  0  0  0  0 \n" + 
					" 36 43  1  0  0  0  0 \n" + 
					"  5 14  1  0  0  0  0 \n" + 
					"  6 25  1  0  0  0  0 \n" + 
					" 19 26  2  0  0  0  0 \n" + 
					" 11 21  1  0  0  0  0 \n" + 
					" 15 23  2  0  0  0  0 \n" + 
					"  2  3  2  0  0  0  0 \n" + 
					" 33 37  1  0  0  0  0 \n" + 
					"  4  5  1  0  0  0  0 \n" + 
					" 11 24  1  0  0  0  0 \n" + 
					"  6  7  1  0  0  0  0 \n" + 
					" 18 28  2  0  0  0  0 \n" + 
					" 10 36  1  0  0  0  0 \n" + 
					" 16 42  1  0  0  0  0 \n" + 
					" 14 25  2  0  0  0  0 \n" + 
					"  5 27  2  0  0  0  0 \n" + 
					"  8 17  1  0  0  0  0 \n" + 
					" 25 26  1  0  0  0  0 \n" + 
					" 29 30  2  0  0  0  0 \n" + 
					"  8 11  1  0  0  0  0 \n" + 
					" 34 40  1  0  0  0  0 \n" + 
					" 24 34  1  0  0  0  0 \n" + 
					" 38 39  1  0  0  0  0 \n" + 
					" 15 38  1  0  0  0  0 \n" + 
					"  3 20  1  0  0  0  0 \n" + 
					" 13 18  1  0  0  0  0 \n" + 
					" 21 40  1  0  0  0  0 \n" + 
					" 38 44  1  0  0  0  0 \n" + 
					" 41 42  1  0  0  0  0 \n" + 
					"  7 17  1  0  0  0  0 \n" + 
					" 28 30  1  0  0  0  0 \n" + 
					"  8 22  2  0  0  0  0 \n" + 
					" 27 33  1  0  0  0  0 \n" + 
					" 31 32  1  0  0  0  0 \n" + 
					" 31 35  1  0  0  0  0 \n" + 
					" 12 13  2  0  0  0  0 \n" + 
					"  4 31  2  0  0  0  0 \n" + 
					" 10 16  1  0  0  0  0 \n" + 
					"  9 10  1  0  0  0  0 \n" + 
					" 42 43  1  0  0  0  0 \n" + 
					"  9 41  1  0  0  0  0 \n" + 
					"  6 45  1  6  0  0  0 \n" + 
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		
		//lineDash.png
		list.add(new Object[]{"lineDash", new TestSpec("moleculeTest/lineDash.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"\n" + 
					"\n" + 
					" 29 31  0  0  0  0  0  0  0  0999 V2000\n" + 
					"    1.2725    3.4561    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.1686    2.9818    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.1686    0.4819    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3026   -0.0376    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4066    0.4518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4795   -2.5091    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.2499   -3.1474    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3177   -2.0481    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.8814   -4.0763    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4367    3.4637    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4367    2.9818    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3102    2.9667    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4292    4.9696    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4735    4.4858    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3177    4.9696    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.2901    1.4738    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4141    1.9577    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4528    0.0378    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.2876    0.4518    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.0527   -4.0962    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.3765   -3.1474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3177   -1.0542    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.4532   -4.9696    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4368   -1.5433    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4367   -1.0240    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3102    1.9878    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.1686    1.4608    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.1686    1.9878    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4141   -2.0481    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  2 28  1  0  0  0  0 \n" + 
					" 12 26  1  0  0  0  0 \n" + 
					"  5 18  1  0  0  0  0 \n" + 
					" 13 14  1  0  0  0  0 \n" + 
					" 16 28  1  0  0  0  0 \n" + 
					"  6 21  1  0  0  0  0 \n" + 
					"  6 24  1  0  0  0  0 \n" + 
					"  4 22  1  0  0  0  0 \n" + 
					"  9 20  1  0  0  0  0 \n" + 
					" 22 24  1  0  0  0  0 \n" + 
					"  1 11  1  0  0  0  0 \n" + 
					" 24 25  1  0  0  0  0 \n" + 
					" 26 27  1  0  0  0  0 \n" + 
					" 24 29  1  1  0  0  0 \n" + 
					"  1  2  1  0  0  0  0 \n" + 
					"  4  3  1  6  0  0  0 \n" + 
					"  4  5  1  0  0  0  0 \n" + 
					"  3 27  1  0  0  0  0 \n" + 
					" 14 15  2  0  0  0  0 \n" + 
					" 18 19  1  1  0  0  0 \n" + 
					" 16 17  2  0  0  0  0 \n" + 
					"  6  7  1  0  0  0  0 \n" + 
					"  6  8  1  1  0  0  0 \n" + 
					" 16 19  1  0  0  0  0 \n" + 
					"  7  9  1  0  0  0  0 \n" + 
					" 10 14  1  0  0  0  0 \n" + 
					" 18 25  1  0  0  0  0 \n" + 
					" 20 21  1  0  0  0  0 \n" + 
					"  9 23  2  0  0  0  0 \n" + 
					" 10 12  1  0  0  0  0 \n" + 
					" 10 11  1  0  0  0  0 \n" + 
					"M  END", Charset.defaultCharset()).build();

			System.out.println(cReal.toMol());
			System.out.println(c.toMol());
			
			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		//OConnectedToDash.png
		list.add(new Object[]{"OConnectedToDash", new TestSpec("moleculeTest/OConnectedToDash.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"\n" + 
					"\n" + 
					" 47 48  0  0  0  0  0  0  0  0999 V2000\n" + 
					"    2.3473    3.3156    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.3515    2.9020    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3604    3.2977    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4623    1.9122    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.4175    1.4332    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.5333   -1.0054    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.4153   -1.4553    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.2974   -2.3815    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4506   -3.2870    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.3223   -0.9340    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.2172   -0.9341    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.2172    2.8674    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.4150    0.4803    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4906   -0.0214    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.2857   -1.4297    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.3633   -0.9626    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.4015   -0.0025    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4103   -2.3799    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.3147   -2.8633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4906   -2.8635    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.3918   -3.8147    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.0888    1.8467    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.3241    0.8342    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.3704    1.9109    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.5173   -2.4671    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.3267    0.4544    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.3704   -0.0214    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.4415    1.4642    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.2380   -0.0449    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4530    1.4580    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.3444   -0.9697    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.2439    2.3744    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.1573   -0.3494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.2385   -2.8592    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.2838    0.8803    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4947    2.8191    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.3704    2.8450    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.2332   -4.2425    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.3373    1.9038    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4578   -1.4403    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.2029   -3.7861    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4506    0.4706    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3391    4.2425    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.3197    0.2414    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.2510    2.2745    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.4667   -1.7826    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.2385    3.6150    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1 37  1  0  0  0  0 \n" + 
					" 13 17  1  0  0  0  0 \n" + 
					"  5 39  1  0  0  0  0 \n" + 
					" 13 14  2  0  0  0  0 \n" + 
					" 11 15  1  0  0  0  0 \n" + 
					"  4 30  1  0  0  0  0 \n" + 
					" 15 16  2  0  0  0  0 \n" + 
					" 12 47  1  0  0  0  0 \n" + 
					" 21 38  1  0  0  0  0 \n" + 
					" 18 40  2  0  0  0  0 \n" + 
					"  4 45  1  6  0  0  0 \n" + 
					" 19  9  1  1  0  0  0 \n" + 
					"  5 13  1  0  0  0  0 \n" + 
					" 12 32  1  1  0  0  0 \n" + 
					"  3 36  1  0  0  0  0 \n" + 
					"  1 12  1  0  0  0  0 \n" + 
					" 26 29  1  0  0  0  0 \n" + 
					" 24 28  1  0  0  0  0 \n" + 
					" 26 27  2  0  0  0  0 \n" + 
					" 26 28  1  0  0  0  0 \n" + 
					" 19 21  1  6  0  0  0 \n" + 
					" 34 41  1  0  0  0  0 \n" + 
					" 24 30  1  0  0  0  0 \n" + 
					" 22 28  1  0  0  0  0 \n" + 
					"  7 31  1  0  0  0  0 \n" + 
					"  2  3  1  0  0  0  0 \n" + 
					" 11 46  1  0  0  0  0 \n" + 
					"  4  5  1  0  0  0  0 \n" + 
					" 18 19  1  0  0  0  0 \n" + 
					"  6  7  2  0  0  0  0 \n" + 
					" 18 20  1  0  0  0  0 \n" + 
					" 30 42  2  0  0  0  0 \n" + 
					" 11 29  1  0  0  0  0 \n" + 
					"  8 15  1  0  0  0  0 \n" + 
					"  8 34  1  0  0  0  0 \n" + 
					"  7 25  1  0  0  0  0 \n" + 
					"  3 43  1  0  0  0  0 \n" + 
					" 12 22  1  0  0  0  0 \n" + 
					" 20 25  1  0  0  0  0 \n" + 
					"  8 19  1  0  0  0  0 \n" + 
					"  4 36  1  1  0  0  0 \n" + 
					" 17 35  1  1  0  0  0 \n" + 
					" 24 37  1  0  0  0  0 \n" + 
					" 33 44  1  0  0  0  0 \n" + 
					" 17 31  1  0  0  0  0 \n" + 
					" 28 23  1  1  0  0  0 \n" + 
					" 11 10  1  1  0  0  0 \n" + 
					" 17 33  1  6  0  0  0 \n" + 
					"M  END", Charset.defaultCharset()).build();

			System.out.println(cReal.toMol());
			System.out.println(c.toMol());
			
			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		//smallRingWithDash.png
		list.add(new Object[]{"smallRingWithDash", new TestSpec("moleculeTest/smallRingWithDash.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"\n" + 
					"\n" + 
					" 53 59  0  0  0  0  0  0  0  0999 V2000\n" + 
					"    1.2977    4.1055    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.8143    3.2452    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.1410    3.1563    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.6728    2.7592    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.5481    2.2468    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.8169    5.0015    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.7670    4.9182    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.9741    3.2221    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.3876    1.4273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.6970    2.3940    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.5287    3.0651    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.6947    2.3738    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.5802    3.8762    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.5615    2.9081    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.7986    2.0337    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.6764    2.7494    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.1722    3.6058    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.8125    3.2452    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.8299   -5.0030    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.3483   -4.1316    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.6613    1.7428    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.5718    2.8846    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.5484    3.1400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.6791   -2.4038    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.3239   -2.3963    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.5736    3.8762    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.2064   -0.1653    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.1989    0.8451    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.9952    1.4273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.1824   -1.5099    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.7984    2.2386    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.3275   -0.6498    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.5336    4.9729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.6773    4.4922    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.7106    4.1166    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.8224   -3.2602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.6839   -0.6460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.6848   -0.6535    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.8233   -3.2602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.1777   -1.5287    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.6923   -2.3888    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.6839   -2.3888    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.1797   -1.5174    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.3257   -2.3738    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.8224   -1.5249    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.6791   -0.6460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.8308   -1.5249    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.5243    3.2377    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.6845    1.4123    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.9567    1.4123    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.9567    0.4357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.6845    0.4357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.8131   -0.0751    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					" 30 40  1  0  0  0  0 \n" + 
					" 11 12  1  0  0  0  0 \n" + 
					" 36 39  2  0  0  0  0 \n" + 
					" 11 14  1  0  0  0  0 \n" + 
					" 19 20  1  0  0  0  0 \n" + 
					" 12 29  1  0  0  0  0 \n" + 
					"  7 34  1  0  0  0  0 \n" + 
					" 20 36  1  0  0  0  0 \n" + 
					" 36 44  1  0  0  0  0 \n" + 
					"  3 35  2  0  0  0  0 \n" + 
					"  4 48  1  0  0  0  0 \n" + 
					" 22 23  1  0  0  0  0 \n" + 
					" 40 41  2  0  0  0  0 \n" + 
					" 52 53  1  0  0  0  0 \n" + 
					" 30 46  2  0  0  0  0 \n" + 
					"  2  3  1  0  0  0  0 \n" + 
					" 33 34  1  0  0  0  0 \n" + 
					"  4 16  1  0  0  0  0 \n" + 
					"  4  5  1  0  0  0  0 \n" + 
					" 37 38  2  0  0  0  0 \n" + 
					" 11 26  2  0  0  0  0 \n" + 
					" 14 15  1  6  0  0  0 \n" + 
					"  6  7  1  0  0  0  0 \n" + 
					" 37 43  1  0  0  0  0 \n" + 
					"  8 18  1  0  0  0  0 \n" + 
					" 44 45  2  0  0  0  0 \n" + 
					"  9 28  1  0  0  0  0 \n" + 
					" 14 22  1  0  0  0  0 \n" + 
					" 49 52  1  0  0  0  0 \n" + 
					" 51 53  1  0  0  0  0 \n" + 
					" 25 47  1  0  0  0  0 \n" + 
					" 38 40  1  0  0  0  0 \n" + 
					" 14 13  1  1  0  0  0 \n" + 
					" 41 42  1  0  0  0  0 \n" + 
					" 45 47  1  0  0  0  0 \n" + 
					" 24 25  2  0  0  0  0 \n" + 
					" 25 39  1  0  0  0  0 \n" + 
					" 28 29  1  0  0  0  0 \n" + 
					" 24 30  1  0  0  0  0 \n" + 
					"  1  2  1  0  0  0  0 \n" + 
					" 27 32  1  0  0  0  0 \n" + 
					" 32 46  1  0  0  0  0 \n" + 
					" 32 47  2  0  0  0  0 \n" + 
					"  1  6  1  1  0  0  0 \n" + 
					"  1  8  1  0  0  0  0 \n" + 
					" 15 50  1  0  0  0  0 \n" + 
					" 15 49  1  0  0  0  0 \n" + 
					" 10  3  1  6  0  0  0 \n" + 
					" 16 17  2  0  0  0  0 \n" + 
					" 50 51  1  0  0  0  0 \n" + 
					" 16 18  1  0  0  0  0 \n" + 
					"  9 10  1  0  0  0  0 \n" + 
					" 42 43  2  0  0  0  0 \n" + 
					" 16 21  2  0  0  0  0 \n" + 
					"  5 48  1  0  0  0  0 \n" + 
					" 28 27  1  1  0  0  0 \n" + 
					" 10 12  1  0  0  0  0 \n" + 
					"  8 31  2  0  0  0  0 \n" + 
					"  1  7  1  6  0  0  0 \n" + 
					"M  END", Charset.defaultCharset()).build();

			System.out.println(cReal.toMol());
			System.out.println(c.toMol());
			
			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		//CO2AsEster
		list.add(new Object[]{"CO2AsEster", new TestSpec("moleculeTest/CO2AsEster.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"\n" + 
					"\n" + 
					" 12 11  0  0  0  0  0  0  0  0999 V2000\n" + 
					"   -2.9387   -0.8013    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4275   -0.7382    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4275    0.2203    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.8907   -0.7361    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3029   -0.7287    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.0650   -0.8013    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.3257   -0.8013    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.3545   -0.7506    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.8902    0.2928    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.4432    0.8067    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.1605   -0.4172    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.4952   -0.8013    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  8  9  1  0  0  0  0 \n" + 
					"  5 12  1  0  0  0  0 \n" + 
					"  2  3  1  0  0  0  0 \n" + 
					"  1  4  1  0  0  0  0 \n" + 
					"  2  5  1  0  0  0  0 \n" + 
					"  2  6  1  0  0  0  0 \n" + 
					"  1 12  1  0  0  0  0 \n" + 
					"  1 11  2  0  0  0  0 \n" + 
					"  6  7  1  0  0  0  0 \n" + 
					"  3 10  1  0  0  0  0 \n" + 
					"  7  8  1  0  0  0  0 \n" + 
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		//parentheticalEthyl.png
		list.add(new Object[]{"parentheticalEthyl", new TestSpec("moleculeTest/parentheticalEthyl.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"  CDK     02091912193D\n" + 
					"\n" + 
					" 55 60  0  0  0  0  0  0  0  0999 V2000\n" + 
					"    9.1612    2.8167    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.0790    2.8046    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.8515    2.8147    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.1094   -4.0882    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.0946    2.8167    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.8598    2.8425    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -8.1417    2.8425    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -9.1139    2.8167    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.8923    2.8425    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    8.1742    2.8425    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.0946    1.8334    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.3347    3.6668    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.8448    4.5170    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.8662    2.8425    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.1057   -3.0865    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.3672    3.6668    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.8477    4.5170    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3662    1.9739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.3273    1.9739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.3987    1.9739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.3598    1.9739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.9522    1.3381    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.9522    0.3622    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.9522   -1.6190    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.9522   -2.5949    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -6.1678    2.8425    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -6.6594    3.6742    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    6.6994    3.6668    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    6.2004    2.8425    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    7.6826    3.6742    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.7629    1.3381    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.1020   -0.1331    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.7629    0.3622    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.7629   -1.6190    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.0946   -1.1237    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.7629   -2.5949    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.8337    2.8425    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -7.6353    1.9739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -6.6594    1.9739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    6.7068    1.9739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    7.6826    1.9739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3662    3.6742    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -7.6427    3.6668    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.3987    3.6742    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.8633    2.8046    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    9.7017    3.6987    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    9.2081    4.6078    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    9.6548    1.9075    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   10.6889    1.8804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.7845   -4.6088    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.0672    2.8147    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -9.6071    1.9073    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -9.0662    1.0255    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -9.6547    3.6985    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  -10.6889    3.6709    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					" 27 43  1  0  0  0  0 \n" + 
					" 32 35  1  0  0  0  0 \n" + 
					" 32 33  2  0  0  0  0 \n" + 
					"  2 26  1  0  0  0  0 \n" + 
					"  4 50  1  0  0  0  0 \n" + 
					"  5 11  1  0  0  0  0 \n" + 
					"  7 38  1  0  0  0  0 \n" + 
					" 15 25  1  0  0  0  0 \n" + 
					"  9 20  1  0  0  0  0 \n" + 
					" 11 22  1  0  0  0  0 \n" + 
					" 10 30  2  0  0  0  0 \n" + 
					" 22 23  2  0  0  0  0 \n" + 
					" 26 27  2  0  0  0  0 \n" + 
					" 40 41  2  0  0  0  0 \n" + 
					" 52 53  1  0  0  0  0 \n" + 
					" 29 40  1  0  0  0  0 \n" + 
					" 12 37  2  0  0  0  0 \n" + 
					" 14 16  2  0  0  0  0 \n" + 
					" 18 19  2  0  0  0  0 \n" + 
					" 14 21  1  0  0  0  0 \n" + 
					" 11 31  2  0  0  0  0 \n" + 
					" 37 45  1  0  0  0  0 \n" + 
					" 15 36  2  0  0  0  0 \n" + 
					" 12 42  1  0  0  0  0 \n" + 
					" 19 37  1  0  0  0  0 \n" + 
					" 16 44  1  0  0  0  0 \n" + 
					"  4 15  1  0  0  0  0 \n" + 
					" 48 49  1  0  0  0  0 \n" + 
					" 26 39  1  0  0  0  0 \n" + 
					"  8 54  1  0  0  0  0 \n" + 
					"  7 43  2  0  0  0  0 \n" + 
					"  8 52  1  0  0  0  0 \n" + 
					" 34 36  1  0  0  0  0 \n" + 
					" 24 35  1  0  0  0  0 \n" + 
					" 34 35  2  0  0  0  0 \n" + 
					" 38 39  2  0  0  0  0 \n" + 
					"  3 51  1  0  0  0  0 \n" + 
					"  3 14  1  0  0  0  0 \n" + 
					" 24 25  2  0  0  0  0 \n" + 
					" 10 41  1  0  0  0  0 \n" + 
					"  2 45  1  0  0  0  0 \n" + 
					" 28 30  1  0  0  0  0 \n" + 
					" 28 29  2  0  0  0  0 \n" + 
					"  6 42  2  0  0  0  0 \n" + 
					" 54 55  1  0  0  0  0 \n" + 
					" 31 33  1  0  0  0  0 \n" + 
					" 23 32  1  0  0  0  0 \n" + 
					"  1 48  1  0  0  0  0 \n" + 
					"  1 46  1  0  0  0  0 \n" + 
					"  6 18  1  0  0  0  0 \n" + 
					" 12 13  1  0  0  0  0 \n" + 
					"  1 10  1  0  0  0  0 \n" + 
					"  5  6  1  0  0  0  0 \n" + 
					" 16 17  1  0  0  0  0 \n" + 
					"  5  9  1  0  0  0  0 \n" + 
					"  7  8  1  0  0  0  0 \n" + 
					" 29 51  1  0  0  0  0 \n" + 
					" 46 47  1  0  0  0  0 \n" + 
					" 20 21  2  0  0  0  0 \n" + 
					"  9 44  2  0  0  0  0 \n" + 
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		//NNring.png
		list.add(new Object[]{"NNring", new TestSpec("moleculeTest/NNring.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"\n" + 
					"\n" + 
					" 28 31  0  0  0  0  0  0  0  0999 V2000\n" + 
					"    0.8397    0.2564    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.6669    0.0289    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4962    0.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.2284   -0.2771    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.6524    0.8851    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.0733    1.4642    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.0825    1.5884    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.0247    0.4384    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.4300    1.1333    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.7838   -0.7859    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.6184    0.0372    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.1964   -1.3401    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4756    0.8686    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.8675    1.7372    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.0247    2.5562    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.1146    1.3981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.0567    0.3143    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.2377    0.5625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.4300    0.6535    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3981   -1.0961    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.2222    0.5873    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.8407    0.2730    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.2532   -1.3318    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.8323   -0.7859    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.8810   -1.7455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.1033   -1.7455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4425   -1.0961    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.2771   -0.2771    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  5 17  1  0  0  0  0\n" + 
					"  9 14  1  0  0  0  0\n" + 
					" 17 18  1  0  0  0  0\n" + 
					"  4 20  2  0  0  0  0\n" + 
					"  4 22  1  0  0  0  0\n" + 
					" 11 21  1  0  0  0  0\n" + 
					" 11 22  2  0  0  0  0\n" + 
					"  2 24  1  0  0  0  0\n" + 
					"  3 13  2  0  0  0  0\n" + 
					" 26 27  1  0  0  0  0\n" + 
					" 20 25  1  0  0  0  0\n" + 
					"  8 21  1  0  0  0  0\n" + 
					"  1  2  2  0  0  0  0\n" + 
					" 10 11  1  0  0  0  0\n" + 
					"  3  4  1  0  0  0  0\n" + 
					" 10 12  2  0  0  0  0\n" + 
					"  3 28  1  0  0  0  0\n" + 
					" 14 15  1  0  0  0  0\n" + 
					"  2 18  1  0  0  0  0\n" + 
					"  1 28  1  0  0  0  0\n" + 
					" 14 16  1  0  0  0  0\n" + 
					"  5  6  2  0  0  0  0\n" + 
					"  5  7  2  0  0  0  0\n" + 
					" 16 21  2  0  0  0  0\n" + 
					" 12 20  1  0  0  0  0\n" + 
					" 23 24  2  0  0  0  0\n" + 
					" 23 27  1  0  0  0  0\n" + 
					" 25 26  2  0  0  0  0\n" + 
					" 27 28  2  0  0  0  0\n" + 
					"  8  9  2  0  0  0  0\n" + 
					"  5 19  1  0  0  0  0\n" + 
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			Chemical c1=c.connectedComponentsAsStream()
						 .map(t->Tuple.of(t,t.getAtomCount()).withVComparator())
						 .max(Comparator.naturalOrder())
						 .map(t->t.k())
						 .orElse(c);
			String keyGot=Inchi.asStdInchi(c1).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		
		//methoxyTerm.png
		list.add(new Object[]{"methoxyTerm", new TestSpec("moleculeTest/methoxyTerm.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"\n" + 
					"\n" + 
					" 23 24  0  0  0  0  0  0  0  0999 V2000\n" + 
					"    3.9696    1.2913    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.9759    1.2782    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.9217   -0.4385    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4674    0.3725    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.5106    0.3853    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.4978    1.3936    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.5106   -0.6143    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.2751   -1.2081    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.5024    0.3811    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.0001   -0.5203    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.0025   -0.5075    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.4889    0.4152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.4879    0.4067    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.4549    0.4002    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.0139   -0.4904    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.9799   -0.4520    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.9802   -0.4605    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.5052   -1.3577    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.4877   -1.3662    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.5176    0.3640    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.5052    0.3683    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.6908   -0.9493    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.4889   -1.3927    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					" 11 20  2  0  0  0  0 \n" + 
					"  1  2  2  0  0  0  0 \n" + 
					"  1 14  1  0  0  0  0 \n" + 
					"  3  4  2  0  0  0  0 \n" + 
					"  3 15  1  0  0  0  0 \n" + 
					" 15 18  1  0  0  0  0 \n" + 
					"  5  6  2  0  0  0  0 \n" + 
					" 13 14  1  0  0  0  0 \n" + 
					" 14 16  2  0  0  0  0 \n" + 
					"  3  8  1  0  0  0  0 \n" + 
					"  5  7  2  0  0  0  0 \n" + 
					" 16 17  1  0  0  0  0 \n" + 
					" 18 19  2  0  0  0  0 \n" + 
					"  4  9  1  0  0  0  0 \n" + 
					"  5  9  1  0  0  0  0 \n" + 
					" 12 17  2  0  0  0  0 \n" + 
					" 11 19  1  0  0  0  0 \n" + 
					" 15 21  2  0  0  0  0 \n" + 
					"  5 12  1  0  0  0  0 \n" + 
					" 20 21  1  0  0  0  0 \n" + 
					"  2 12  1  0  0  0  0 \n" + 
					"  8 22  3  0  0  0  0 \n" + 
					" 10 23  1  0  0  0  0 \n" + 
					" 10 11  1  0  0  0  0 \n" + 
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		
		
		//NHCloseToChain.png
		list.add(new Object[]{"NHCloseToChain", new TestSpec("moleculeTest/NHCloseToChain.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"\n" + 
					"\n" + 
					" 25 25  0  0  0  0  0  0  0  0999 V2000\n" + 
					"   -0.6629    0.5343    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.6579    1.4952    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.3524   -1.4952    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.4593   -1.0121    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.7159    0.5061    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.5781    0.0211    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.0175   -1.4722    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.8801   -0.9968    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.1817    0.0153    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.5096    0.4907    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.5128    0.0444    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.0060    0.5176    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.8840    0.0230    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.3951    0.5176    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.2616    0.0153    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.1626   -0.9930    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.5517   -0.9930    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.4181   -1.4952    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.7351   -1.4722    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.6905   -0.1131    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.9895   -1.0658    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.6617   -1.4952    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.4687    0.4447    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.2623   -0.1150    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.9556   -1.0658    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1  2  2  0  0  0  0 \n" + 
					"  3  4  2  0  0  0  0 \n" + 
					"  1  9  1  0  0  0  0 \n" + 
					" 12 13  1  0  0  0  0 \n" + 
					"  5  6  1  0  0  0  0 \n" + 
					" 14 15  1  0  0  0  0 \n" + 
					" 17 18  1  0  0  0  0 \n" + 
					" 17 19  1  0  0  0  0 \n" + 
					"  7  8  1  0  0  0  0 \n" + 
					"  4 11  1  0  0  0  0 \n" + 
					"  4 22  1  0  0  0  0 \n" + 
					"  7 16  1  0  0  0  0 \n" + 
					" 21 25  1  0  0  0  0 \n" + 
					"  8 19  1  0  0  0  0 \n" + 
					"  5 13  1  0  0  0  0 \n" + 
					"  6 14  1  0  0  0  0 \n" + 
					" 10 20  1  0  0  0  0 \n" + 
					" 16 22  1  0  0  0  0 \n" + 
					" 20 21  2  0  0  0  0 \n" + 
					" 23 24  1  0  0  0  0 \n" + 
					"  1 11  1  0  0  0  0 \n" + 
					" 24 25  2  0  0  0  0 \n" + 
					" 20 23  1  0  0  0  0 \n" + 
					"  9 12  1  0  0  0  0 \n" + 
					" 10 11  2  0  0  0  0 \n" + 
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		
		//colinearAromaticBond.png
		list.add(new Object[]{"colinearAromaticBond", new TestSpec("moleculeTest/colinearAromaticBond.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"\n" + 
					"\n" + 
					" 39 41  0  0  0  0  0  0  0  0999 V2000\n" + 
					"   -1.2417    1.5482    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.6595    2.3517    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.3186    0.7781    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.7592    1.5720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.9530    0.0476    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.4016   -0.7542    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.5993   -0.7860    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.0519   -0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.2576   -1.5164    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.6746   -2.3183    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.9030   -2.3500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.6883   -2.4294    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.7669    2.6139    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.8636    3.2472    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.7419   -2.0752    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.8319   -1.4450    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.9347    2.3421    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.7528    2.3342    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.7370   -2.3500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    6.6420   -2.4136    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    6.6659    2.2786    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.7041    2.2627    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.0845    2.2706    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4637    2.2786    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.0607   -2.4215    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4954   -2.4215    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.7528    1.4291    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.7290   -3.2472    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.5914    0.7979    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.7552    0.0357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -6.6643    0.0159    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.0360    1.5958    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.5676    2.3659    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.0519   -1.5323    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.4016    0.8019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.7022   -0.7701    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.3066   -0.7542    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.5914   -2.3262    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.6863    0.7860    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					" 30 37  2  0  0  0  0 \n" + 
					"  1 39  1  0  0  0  0 \n" + 
					" 34 38  1  0  0  0  0 \n" + 
					" 32 33  2  0  0  0  0 \n" + 
					" 15 19  1  0  0  0  0 \n" + 
					" 11 12  1  0  0  0  0 \n" + 
					" 13 17  1  0  0  0  0 \n" + 
					" 13 14  2  0  0  0  0 \n" + 
					" 11 15  1  0  0  0  0 \n" + 
					" 15 16  2  0  0  0  0 \n" + 
					" 12 25  1  0  0  0  0 \n" + 
					" 19 20  2  0  0  0  0 \n" + 
					" 30 31  1  0  0  0  0 \n" + 
					" 13 18  1  0  0  0  0 \n" + 
					"  7 34  2  0  0  0  0 \n" + 
					" 19 28  1  0  0  0  0 \n" + 
					"  7 36  1  0  0  0  0 \n" + 
					" 10 26  1  0  0  0  0 \n" + 
					"  9 36  2  0  0  0  0 \n" + 
					"  8 29  1  0  0  0  0 \n" + 
					"  5 35  1  0  0  0  0 \n" + 
					" 17 22  1  0  0  0  0 \n" + 
					"  3 35  2  0  0  0  0 \n" + 
					"  2 24  1  0  0  0  0 \n" + 
					" 22 23  1  0  0  0  0 \n" + 
					"  1  2  2  0  0  0  0 \n" + 
					" 29 39  2  0  0  0  0 \n" + 
					"  3  4  1  0  0  0  0 \n" + 
					"  5  6  2  0  0  0  0 \n" + 
					"  3 30  1  0  0  0  0 \n" + 
					" 29 32  1  0  0  0  0 \n" + 
					"  5  8  1  0  0  0  0 \n" + 
					"  7  8  1  0  0  0  0 \n" + 
					" 10 38  2  0  0  0  0 \n" + 
					" 18 27  1  0  0  0  0 \n" + 
					"  9 10  1  0  0  0  0 \n" + 
					"  6 37  1  0  0  0  0 \n" + 
					" 23 24  1  0  0  0  0 \n" + 
					" 25 26  1  0  0  0  0 \n" + 
					"  2 33  1  0  0  0  0 \n" + 
					" 18 21  2  0  0  0  0 \n" + 
					"M  END", Charset.defaultCharset()).build();

			System.out.println(cReal.toMol());
			System.out.println(c.toMol());
			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		
		//terminalGroupWithPhAtEnd.png
		list.add(new Object[]{"terminalGroupWithPhAtEnd", new TestSpec("moleculeTest/terminalGroupWithPhAtEnd.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"  CDK     02091912223D\n" + 
					"\n" + 
					" 31 33  0  0  0  0  0  0  0  0999 V2000\n" + 
					"    4.7276   -1.3075    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.2408   -1.3005    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.7719   -0.4348    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.7731   -2.1802    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.2690   -1.3040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.7662   -2.1652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.7458   -1.2812    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.7768   -0.4311    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.2841    0.4228    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.7774   -0.4085    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.3003    0.4228    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.2362   -0.4348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.7328    0.4254    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.7794   -2.1802    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.2525   -0.4198    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.2305   -0.4198    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.2777   -1.3075    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.2525   -2.1652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.7409   -1.2850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.2305   -2.1652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.7478   -1.2925    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.2338    1.2894    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.7328    0.4216    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.2507   -0.4552    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.7741    0.4239    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.2971    1.2762    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.8205    2.1553    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.8209    2.1821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.2978    1.3298    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.7744    0.4507    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.2996    0.3850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					" 15 19  1  0  0  0  0 \n" + 
					"  7 20  2  0  0  0  0 \n" + 
					" 15 16  2  0  0  0  0 \n" + 
					"  7 16  1  0  0  0  0 \n" + 
					" 13 23  1  0  0  0  0 \n" + 
					" 12  3  1  6  0  0  0 \n" + 
					" 13 22  2  0  0  0  0 \n" + 
					" 24 25  1  0  0  0  0 \n" + 
					" 26 27  1  0  0  0  0 \n" + 
					" 19 21  1  0  0  0  0 \n" + 
					" 28 29  1  0  0  0  0 \n" + 
					"  2  3  1  0  0  0  0 \n" + 
					"  2  5  1  0  0  0  0 \n" + 
					"  2 14  2  0  0  0  0 \n" + 
					"  1  7  1  0  0  0  0 \n" + 
					"  1 24  1  0  0  0  0 \n" + 
					"  4  5  2  0  0  0  0 \n" + 
					"  4  6  1  0  0  0  0 \n" + 
					" 12 13  1  0  0  0  0 \n" + 
					"  5  8  1  0  0  0  0 \n" + 
					" 18 19  2  0  0  0  0 \n" + 
					" 18 20  1  0  0  0  0 \n" + 
					" 10 17  1  0  0  0  0 \n" + 
					"  8  9  1  0  0  0  0 \n" + 
					"  8 10  2  0  0  0  0 \n" + 
					" 11 31  1  0  0  0  0 \n" + 
					"  6 17  2  0  0  0  0 \n" + 
					" 12 21  1  0  0  0  0 \n" + 
					" 25 26  2  0  0  0  0 \n" + 
					" 27 28  2  0  0  0  0 \n" + 
					" 29 30  2  0  0  0  0 \n" + 
					" 25 30  1  0  0  0  0 \n" + 
					" 10 11  1  0  0  0  0 \n" + 
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		//smallImage.png
		list.add(new Object[]{"smallImage", new TestSpec("moleculeTest/smallImage.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"  MJ150420                      \n" + 
					"\n" + 
					" 57 58  0  0  0  0  0  0  0  0999 V2000\n" + 
					"   -1.5225   -2.0160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.7370   -2.4088    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.3633    0.5367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.1095    0.1178    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -7.8718    0.3272    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -6.6020    0.9425    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    7.3764   -1.0558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.1400   -2.0030    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    7.3436   -1.9279    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.6617   -2.3826    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.7714   -2.0160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.7510    1.3483    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.8658   -1.9226    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.8610   -1.0735    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    6.5941   -2.3303    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.5466   -0.4058    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.6694   -0.7069    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.6580    0.1385    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.1561   -2.3433    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.5015   -1.9925    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.9102   -1.9715    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.6739   -2.4094    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -7.3351   -0.3011    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -8.0027   -0.6938    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.9416    0.5629    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -6.6150    1.7673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -6.5719   -2.4313    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -7.3351   -2.0160    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.4209   -1.9375    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    8.0865   -3.1811    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    8.1192   -2.3237    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -8.8033   -0.3493    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -8.7882    0.4425    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.8426   -1.1520    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -6.5627   -1.5578    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    8.0930   -0.6153    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    8.0865    0.2225    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.6694   -3.2728    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -9.5265   -0.7462    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.1116   -2.3695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.1357    1.8328    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.8499    2.2928    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    9.5265   -2.3303    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    8.8066   -1.9113    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.0703    0.9556    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    8.8066   -1.0342    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.8361   -2.0226    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    6.6334   -0.6283    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -7.3220   -1.1389    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -6.5758   -3.2728    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.9917   -1.3222    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.7902    3.2728    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -7.4922    1.0997    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.7423   -2.3433    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.1297   -2.3570    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.4419   -2.0030    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.0703   -2.4219    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					" 32 39  1  0  0  0  0\n" + 
					" 32 33  2  0  0  0  0\n" + 
					" 27 47  1  0  0  0  0\n" + 
					" 36 37  1  0  0  0  0\n" + 
					" 47 57  1  0  0  0  0\n" + 
					" 19 20  1  0  0  0  0\n" + 
					" 11 19  1  0  0  0  0\n" + 
					" 30 31  1  0  0  0  0\n" + 
					"  9 15  1  0  0  0  0\n" + 
					"  5 53  1  0  0  0  0\n" + 
					"  7 36  2  0  0  0  0\n" + 
					" 43 44  1  0  0  0  0\n" + 
					"  6 26  1  0  0  0  0\n" + 
					"  4 45  1  0  0  0  0\n" + 
					" 36 46  1  0  0  0  0\n" + 
					" 10 29  1  0  0  0  0\n" + 
					" 26 42  1  0  0  0  0\n" + 
					" 56 57  1  0  0  0  0\n" + 
					" 34 47  1  0  0  0  0\n" + 
					"  4 16  1  0  0  0  0\n" + 
					"  2  8  1  0  0  0  0\n" + 
					"  5 23  1  0  0  0  0\n" + 
					" 44 46  2  0  0  0  0\n" + 
					"  7 48  1  0  0  0  0\n" + 
					" 12 45  1  0  0  0  0\n" + 
					" 18 25  1  0  0  0  0\n" + 
					" 21 22  1  0  0  0  0\n" + 
					" 13 55  1  0  0  0  0\n" + 
					" 27 50  1  0  0  0  0\n" + 
					" 23 49  1  0  0  0  0\n" + 
					"  8 54  1  0  0  0  0\n" + 
					" 22 38  2  0  0  0  0\n" + 
					" 34 51  1  0  0  0  0\n" + 
					"  1 40  1  0  0  0  0\n" + 
					" 31 44  1  0  0  0  0\n" + 
					" 35 49  1  0  0  0  0\n" + 
					" 24 32  1  0  0  0  0\n" + 
					" 34 35  1  0  0  0  0\n" + 
					" 13 14  2  0  0  0  0\n" + 
					"  3 18  1  0  0  0  0\n" + 
					" 13 15  1  0  0  0  0\n" + 
					" 17 18  2  0  0  0  0\n" + 
					" 21 40  1  0  0  0  0\n" + 
					" 41 42  2  0  0  0  0\n" + 
					" 20 54  1  0  0  0  0\n" + 
					" 16 51  1  0  0  0  0\n" + 
					"  9 31  2  0  0  0  0\n" + 
					"  1  2  1  0  0  0  0\n" + 
					" 29 55  1  0  0  0  0\n" + 
					"  3  4  1  0  0  0  0\n" + 
					" 42 52  1  0  0  0  0\n" + 
					"  7  9  1  0  0  0  0\n" + 
					"  6 12  1  0  0  0  0\n" + 
					" 23 24  1  0  0  0  0\n" + 
					" 22 56  1  0  0  0  0\n" + 
					" 27 28  2  0  0  0  0\n" + 
					" 10 11  1  0  0  0  0\n" + 
					"  6 53  1  0  0  0  0\n" + 
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});


		//dashedMethylNoLabelWithAtomNumbers.png
		list.add(new Object[]{"dashedToPhenyl", new TestSpec("moleculeTest/dashedToPhenyl.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" +
					"  CDK     02021921273D\n" +
					"\n" +
					" 12 13  0  0  0  0  0  0  0  0999 V2000\n" +
					"   -0.3758   -0.0902    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -0.8410    1.5370    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -0.8794    2.5706    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -1.6686    0.9922    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -0.0601    0.9922    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -1.3529    0.0451    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    0.1203   -0.8418    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -0.3271   -1.7362    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    0.2237   -2.5708    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.2219   -2.5111    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.6693   -1.6167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.1185   -0.7821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"  8  9  1  0  0  0  0 \n" +
					"  7 12  1  0  0  0  0 \n" +
					"  9 10  2  0  0  0  0 \n" +
					"  2  3  2  0  0  0  0 \n" +
					"  1  5  1  0  0  0  0 \n" +
					"  2  4  1  0  0  0  0 \n" +
					"  1  6  1  0  0  0  0 \n" +
					"  2  5  1  0  0  0  0 \n" +
					"  1  7  1  6  0  0  0 \n" +
					"  4  6  1  0  0  0  0 \n" +
					" 11 12  2  0  0  0  0 \n" +
					"  7  8  2  0  0  0  0 \n" +
					" 10 11  1  0  0  0  0 \n" +
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});

		//dashedMethylNoLabelWithAtomNumbers.png
		list.add(new Object[]{"dashedMethylNoLabelWithAtomNumbers", new TestSpec("moleculeTest/dashedMethylNoLabelWithAtomNumbers.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" +
					"  MJ150420                      \n" +
					"\n" +
					" 38 38  0  0  0  0  0  0  0  0999 V2000\n" +
					"   -3.8527    0.8463    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -1.7491    0.5050    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    2.5128    0.5050    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.0654   -0.3818    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -3.1503   -2.8454    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -0.3695    2.0817    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    0.3633    1.6813    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -2.4513   -1.6136    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -3.1903   -2.0078    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -1.0347   -1.6136    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -1.7522   -2.0108    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    0.3941   -1.6136    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -0.3448   -2.0078    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    2.5050    2.0932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -0.2709    0.4926    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -1.7614    2.8454    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -1.7799   -2.8454    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.0654    2.0817    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.7815   -0.7840    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    2.4943   -0.3818    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -3.9048   -1.5890    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -3.9048   -0.7883    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.7861    1.6752    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.7861   -1.5890    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.0839   -2.0078    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    3.9048    2.0817    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    3.2149    1.6875    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -2.5375   -0.7390    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -3.1903   -0.3757    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -1.7614    2.0878    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    0.3818    0.8561    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -1.0593    1.6875    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.0839    0.4557    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -2.4759    0.8745    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -2.4759    1.6752    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.7861    0.8745    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -3.1903    0.4434    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -1.0732   -0.8393    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					" 22 29  1  0  0  0  0\n" +
					" 30 35  1  0  0  0  0\n" +
					" 37  1  1  6  0  0  0\n" +
					" 30 32  1  0  0  0  0\n" +
					" 34 35  1  0  0  0  0\n" +
					" 12 25  1  0  0  0  0\n" +
					" 34 37  1  0  0  0  0\n" +
					" 14 27  1  0  0  0  0\n" +
					" 29 28  1  1  0  0  0\n" +
					" 11 17  1  1  0  0  0\n" +
					"  9  5  1  6  0  0  0\n" +
					"  4 19  1  0  0  0  0\n" +
					"  7 18  1  1  0  0  0\n" +
					" 23 36  1  0  0  0  0\n" +
					" 24 25  2  0  0  0  0\n" +
					" 30 16  1  1  0  0  0\n" +
					" 26 27  2  0  0  0  0\n" +
					"  7 31  1  0  0  0  0\n" +
					" 29 37  1  0  0  0  0\n" +
					" 33 36  1  0  0  0  0\n" +
					" 31 33  1  0  0  0  0\n" +
					" 10 13  1  0  0  0  0\n" +
					"  9 21  1  0  0  0  0\n" +
					" 10 11  1  0  0  0  0\n" +
					" 12 13  2  0  0  0  0\n" +
					"  6  7  1  0  0  0  0\n" +
					" 14 23  2  0  0  0  0\n" +
					"  6 32  1  0  0  0  0\n" +
					" 33  4  1  6  0  0  0\n" +
					" 19 24  1  0  0  0  0\n" +
					" 31 15  1  1  0  0  0\n" +
					" 19 20  2  0  0  0  0\n" +
					" 21 22  2  0  0  0  0\n" +
					"  8 11  1  0  0  0  0\n" +
					"  8  9  1  0  0  0  0\n" +
					" 34  2  1  6  0  0  0\n" +
					" 36  3  1  6  0  0  0\n" +
					" 10 38  1  6  0  0  0\n" +
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			
			Chemical c1=c.connectedComponentsAsStream()
					.map(ct->Tuple.of(ct,ct.getAtomCount()).withVComparator())
					.max(Comparator.naturalOrder())
					.map(t->t.k())
					.orElse(c);
			String keyGot=Inchi.asStdInchi(c1).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		
		//dashBondVerticalLinesToH.png
		list.add(new Object[]{"dashBondVerticalLinesToH", new TestSpec("moleculeTest/dashBondVerticalLinesToH.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"  CDK     02061900183D\n" + 
					"\n" + 
					" 50 53  0  0  0  0  0  0  0  0999 V2000\n" + 
					"   -1.8041   -2.4756    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.2442   -1.7496    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.7271   -1.9649    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.3175   -1.5923    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.8020   -2.4756    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.8085   -2.4606    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.8202    2.6769    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.8116    2.7445    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.3224    1.8657    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.3138    1.8732    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.8383    1.0211    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.3288    0.1232    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.8246    2.7445    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.8160   -0.7105    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.3268    0.1232    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.8267   -0.7331    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.6199   -3.4834    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.7271   -2.9713    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.8095   -0.7256    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.3138   -1.5818    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.7301   -0.7932    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.8321   -0.6013    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.3258    0.2584    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.2118   -3.2567    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.3258   -1.4691    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.5833   -4.4735    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.8310    0.9944    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.6058   -1.4579    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.3042   -1.6043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.2582   -3.4070    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.3248   -3.4670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.4771   -2.9638    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.7977   -2.4756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.3313   -1.4391    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.8431   -0.5904    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.3313    0.2584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.4771   -1.9649    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.5307   -0.4477    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.3203    1.8657    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -6.3264    0.0971    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.7753   -4.3598    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.7737   -4.4161    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    6.3216   -3.5795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.8711   -2.6868    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.8727   -2.6305    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.3303    3.6138    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.8360    4.4765    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.8360    4.4699    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.3303    3.6006    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.8245    2.7379    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					" 32 37  1  0  0  0  0 \n" + 
					"  6 29  1  0  0  0  0 \n" + 
					" 11 12  1  0  0  0  0 \n" + 
					"  5 20  1  0  0  0  0 \n" + 
					" 32 30  1  1  0  0  0 \n" + 
					" 19 20  2  0  0  0  0 \n" + 
					" 30 31  1  0  0  0  0 \n" + 
					" 43 44  1  0  0  0  0 \n" + 
					" 13 46  2  0  0  0  0 \n" + 
					" 47 48  2  0  0  0  0 \n" + 
					" 22 23  2  0  0  0  0 \n" + 
					" 22 25  1  0  0  0  0 \n" + 
					"  2 21  2  0  0  0  0 \n" + 
					"  2  3  1  0  0  0  0 \n" + 
					"  4 16  2  0  0  0  0 \n" + 
					" 13 50  1  0  0  0  0 \n" + 
					" 14 15  1  0  0  0  0 \n" + 
					"  4 20  1  0  0  0  0 \n" + 
					" 11 27  2  0  0  0  0 \n" + 
					" 14 19  1  0  0  0  0 \n" + 
					" 10 39  1  6  0  0  0 \n" + 
					"  8  9  1  0  0  0  0 \n" + 
					" 44 45  2  0  0  0  0 \n" + 
					" 18 24  1  0  0  0  0 \n" + 
					" 48 49  1  0  0  0  0 \n" + 
					" 22 35  1  0  0  0  0 \n" + 
					"  2 33  1  0  0  0  0 \n" + 
					" 24 33  1  0  0  0  0 \n" + 
					" 28 38  1  1  0  0  0 \n" + 
					" 28 37  1  0  0  0  0 \n" + 
					" 34 35  1  0  0  0  0 \n" + 
					" 31 45  1  0  0  0  0 \n" + 
					" 14 29  2  0  0  0  0 \n" + 
					" 15 40  1  0  0  0  0 \n" + 
					" 17 18  1  0  0  0  0 \n" + 
					" 17 26  2  0  0  0  0 \n" + 
					" 41 42  1  0  0  0  0 \n" + 
					" 49 50  2  0  0  0  0 \n" + 
					" 33  1  1  6  0  0  0 \n" + 
					" 37 25  1  6  0  0  0 \n" + 
					"  1  4  1  0  0  0  0 \n" + 
					" 35 36  1  0  0  0  0 \n" + 
					"  3 28  1  0  0  0  0 \n" + 
					"  5  6  2  0  0  0  0 \n" + 
					" 12 19  1  0  0  0  0 \n" + 
					"  7  8  2  0  0  0  0 \n" + 
					" 10  9  1  1  0  0  0 \n" + 
					" 17 32  1  0  0  0  0 \n" + 
					" 42 43  2  0  0  0  0 \n" + 
					" 46 47  1  0  0  0  0 \n" + 
					" 31 41  2  0  0  0  0 \n" + 
					" 10 13  1  0  0  0  0 \n" + 
					" 10 11  1  0  0  0  0 \n" + 
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			
			Chemical c1=c.connectedComponentsAsStream()
					.map(ct->Tuple.of(ct,ct.getAtomCount()).withVComparator())
					.max(Comparator.naturalOrder())
					.map(t->t.k())
					.orElse(c);
			String keyGot=Inchi.asStdInchi(c1).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		//wedgeBondWeirdAngle.png
		list.add(new Object[]{"wedgeBondWeirdAngle", new TestSpec("moleculeTest/wedgeBondWeirdAngle.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" +
					"  CDK     02021919453D\n" +
					"\n" +
					" 29 32  0  0  0  0  0  0  0  0999 V2000\n" +
					"   -1.6302    2.1505    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -2.4763    2.7022    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.8919   -0.6634    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -0.6388   -1.4498    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -1.4551   -1.8859    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    0.0748    1.1338    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -0.7767    0.6336    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    2.4146    0.1253    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    0.0880   -0.8572    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -0.7794   -0.3587    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.8695    0.9467    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -1.6414   -0.8572    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    0.9377   -0.3541    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -3.2739    2.1245    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.8943   -1.8629    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.0645   -1.3045    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    0.2780   -1.8710    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -2.8995    1.2179    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    0.0656    2.1245    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -0.7767    2.6313    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -1.9620    1.1852    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -3.4528    0.3205    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    0.9377    0.6448    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    3.4498    0.1416    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    1.5788   -2.8252    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    0.5948   -2.8252    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"   -2.4465    3.7048    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    0.0656   -3.7048    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"    2.0633   -3.7048    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
					"  1 20  1  0  0  0  0 \n" +
					"  2 27  1  0  0  0  0 \n" +
					" 15  3  1  1  0  0  0 \n" +
					" 15 16  1  0  0  0  0 \n" +
					" 19 20  1  0  0  0  0 \n" +
					" 17 26  1  0  0  0  0 \n" +
					"  9 13  1  0  0  0  0 \n" +
					"  8 24  2  0  0  0  0 \n" +
					"  6 23  1  0  0  0  0 \n" +
					" 13 23  2  0  0  0  0 \n" +
					" 15 25  1  0  0  0  0 \n" +
					"  3 13  1  0  0  0  0 \n" +
					" 26 28  1  6  0  0  0 \n" +
					"  1  2  2  0  0  0  0 \n" +
					"  2 14  1  0  0  0  0 \n" +
					"  4  5  1  0  0  0  0 \n" +
					" 14 18  1  0  0  0  0 \n" +
					"  6 19  1  0  0  0  0 \n" +
					" 11 23  1  0  0  0  0 \n" +
					" 17  4  1  1  0  0  0 \n" +
					"  3  8  1  0  0  0  0 \n" +
					"  6  7  2  0  0  0  0 \n" +
					" 16 17  1  0  0  0  0 \n" +
					"  7 10  1  0  0  0  0 \n" +
					"  9 10  2  0  0  0  0 \n" +
					" 18 22  2  0  0  0  0 \n" +
					" 25 26  1  0  0  0  0 \n" +
					" 25 29  1  6  0  0  0 \n" +
					"  1 21  1  0  0  0  0 \n" +
					" 18 21  1  0  0  0  0 \n" +
					"  8 11  1  0  0  0  0 \n" +
					" 10 12  1  0  0  0  0 \n" +
					"M  END", Charset.defaultCharset()).build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		
		

		//closeNonBondedAtoms.png
		list.add(new Object[]{"closeNonBondedAtoms", new TestSpec("moleculeTest/closeNonBondedAtoms.png", c->{
			Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
					"  CDK     02071919263D\n" + 
					"\n" + 
					" 68 71  0  0  0  0  0  0  0  0999 V2000\n" + 
					"    4.6945   -0.8455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.0594   -3.8374    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.6286    5.2578    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -6.5345    4.7747    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4341    2.2457    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4469    1.7550    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.1500    2.2382    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.0438    1.7429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.0680    2.2382    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.1380    1.7197    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3075    0.7374    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.1742    0.2452    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.2225   -0.7675    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.0735   -1.2571    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.6377   -2.2943    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    6.2930   -3.0523    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.8073   -2.6295    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.0589    3.7510    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.1621    3.2649    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.9135    3.2649    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4359    0.2760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4348   -0.7634    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.8073    4.7747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.9135    5.2760    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.9256    0.7002    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.7952    0.2459    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.1621   -4.2843    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.0801   -3.8012    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.9135   -0.7634    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.6219   -1.2774    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    6.5345   -0.7634    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.0438   -2.2642    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.1742   -2.7805    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.9074   -2.7805    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.7228   -2.3517    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.0317   -4.2843    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.9135   -3.7891    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.2803    2.2382    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.3287   -3.8012    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4228   -4.2903    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4590   -3.8012    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.2164    5.2518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3045    4.7626    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.8098   -3.9461    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.9376   -4.2843    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4831    3.3011    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.7590   -3.8012    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4348    0.7465    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.6649    0.7465    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.4348   -2.7805    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4107   -2.2913    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.9135    2.2744    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3045   -1.2646    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.9195   -1.2646    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.0680   -0.7815    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.1802    0.7465    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.0438    0.7404    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.9256    1.7309    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.0680    0.2452    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3075    3.7571    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.0438    4.7626    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.2924   -4.2843    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.4348   -5.2747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.3045    1.7429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.9135   -1.7659    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    6.2359   -4.8503    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.0332    3.2371    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.2953    2.7530    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					" 32 33  1  0  0  0  0 \n" + 
					" 32 34  1  0  0  0  0 \n" + 
					" 36 37  1  0  0  0  0 \n" + 
					" 11 12  1  0  0  0  0 \n" + 
					" 15 17  1  1  0  0  0 \n" + 
					" 15 16  1  0  0  0  0 \n" + 
					" 30 31  2  0  0  0  0 \n" + 
					"  9 58  1  0  0  0  0 \n" + 
					" 11 64  1  0  0  0  0 \n" + 
					" 11 21  1  1  0  0  0 \n" + 
					" 25 59  1  0  0  0  0 \n" + 
					"  1 54  1  6  0  0  0 \n" + 
					" 27 39  1  0  0  0  0 \n" + 
					"  2 62  1  6  0  0  0 \n" + 
					" 40 41  1  0  0  0  0 \n" + 
					" 61 24  1  6  0  0  0 \n" + 
					"  1 30  1  0  0  0  0 \n" + 
					"  1 26  1  0  0  0  0 \n" + 
					"  2 36  1  0  0  0  0 \n" + 
					" 15 30  1  0  0  0  0 \n" + 
					" 41 62  1  0  0  0  0 \n" + 
					" 18 19  1  0  0  0  0 \n" + 
					" 18 20  1  1  0  0  0 \n" + 
					" 13 53  2  0  0  0  0 \n" + 
					" 10 38  1  0  0  0  0 \n" + 
					"  8 57  1  0  0  0  0 \n" + 
					" 16 44  1  0  0  0  0 \n" + 
					"  6 38  1  0  0  0  0 \n" + 
					" 21 22  2  0  0  0  0 \n" + 
					" 44 47  1  0  0  0  0 \n" + 
					" 25 26  1  0  0  0  0 \n" + 
					"  2 33  1  0  0  0  0 \n" + 
					" 18 61  1  0  0  0  0 \n" + 
					" 34 37  1  0  0  0  0 \n" + 
					" 34 35  1  1  0  0  0 \n" + 
					" 14 29  2  0  0  0  0 \n" + 
					" 42 61  1  0  0  0  0 \n" + 
					" 13 14  1  0  0  0  0 \n" + 
					" 14 32  1  0  0  0  0 \n" + 
					" 60 46  1  1  0  0  0 \n" + 
					" 45 47  2  0  0  0  0 \n" + 
					" 32 65  1  1  0  0  0 \n" + 
					" 26 49  1  1  0  0  0 \n" + 
					" 20 52  1  0  0  0  0 \n" + 
					" 40 63  1  0  0  0  0 \n" + 
					" 28 45  1  0  0  0  0 \n" + 
					" 54 55  1  0  0  0  0 \n" + 
					"  3  4  1  0  0  0  0 \n" + 
					" 25 58  2  0  0  0  0 \n" + 
					" 39 40  2  0  0  0  0 \n" + 
					"  5  6  1  0  0  0  0 \n" + 
					" 12 13  1  0  0  0  0 \n" + 
					" 12 57  1  0  0  0  0 \n" + 
					" 50 51  1  0  0  0  0 \n" + 
					"  7  8  1  0  0  0  0 \n" + 
					" 43 60  1  0  0  0  0 \n" + 
					" 19 60  1  0  0  0  0 \n" + 
					"  9 10  1  0  0  0  0 \n" + 
					" 21 48  1  0  0  0  0 \n" + 
					" 42 43  1  0  0  0  0 \n" + 
					"  5 46  1  0  0  0  0 \n" + 
					"  3 23  1  0  0  0  0 \n" + 
					" 23 24  1  0  0  0  0 \n" + 
					" 27 28  2  0  0  0  0 \n" + 
					" 10 56  2  0  0  0  0 \n" + 
					"  7 64  1  0  0  0  0 \n" + 
					"  5 68  1  6  0  0  0 \n" + 
					"  6 48  1  6  0  0  0 \n" + 
					"  9 67  1  6  0  0  0 \n" + 
					" 41 50  1  6  0  0  0 \n" + 
					" 44 66  1  6  0  0  0 \n" + 
					"M  END").build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		//explicitCarbonStructureWithTallerFont.png
		list.add(new Object[]{"explicitCarbonStructureWithTallerFont", new TestSpec("moleculeTest/explicitCarbonStructureWithTallerFont.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("CC1CC2OC2CC1COC(=O)C3CC4OC4CC3C").build();

			
			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		//chainOnEdge.png
		list.add(new Object[]{"chainOnEdge", new TestSpec("moleculeTest/chainOnEdge.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("COc1ccc(CCN2CCCc3cc(O)c(OC)cc23)cc1O").build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		
		/*
		//SandOCloseTogether.png
				list.add(new Object[]{"SandOCloseTogether", new TestSpec("moleculeTest/SandOCloseTogether.png", c->{
					Chemical cReal=ChemicalBuilder.createFromMol("\n" + 
							"  CDK     02281917583D\n" + 
							"\n" + 
							" 82 91  0  0  0  0  0  0  0  0999 V2000\n" + 
							"    5.9870    5.8640    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    5.1213    5.3267    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.8274    5.8043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -1.6934    5.3028    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    4.2597    5.8042    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    3.4337    5.2794    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -1.6921    4.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -2.6250    3.8217    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.8044    2.7957    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -1.6838    2.3028    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -5.9698    0.9818    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -6.8576    0.5204    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.8182    3.8017    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    0.0420    4.3031    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.2395   -2.2215    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    0.0807   -1.3067    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    3.3419   -2.9532    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    4.2088   -2.4508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -4.5010    2.0469    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -3.7610    2.7242    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.8731   -4.7895    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.1379   -5.4397    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -1.8237   -5.0766    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -5.3379   -0.5014    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -6.3179   -0.5990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.6278   -2.5047    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    2.6056   -2.2881    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -2.7632    2.8192    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    0.0120    0.3671    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.4911   -0.5191    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -4.7621   -1.3786    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    3.3435   -3.9572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    4.2128   -4.4560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -3.7693   -3.0110    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -4.6754   -3.1188    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -4.1091   -4.5735    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -4.9853   -4.0838    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -4.3157    1.0268    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.0015    0.3638    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.5045   -0.5191    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    2.5017   -0.5180    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    3.0104   -1.3819    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -5.1319    0.4538    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -3.3762   -3.9172    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.7519   -3.7651    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    0.1219   -3.1643    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    4.0065   -1.4852    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    2.2075    6.7604    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    2.5206    5.8042    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    0.0417    5.3042    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    0.9984    5.7849    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -4.1558    3.6420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -2.3118   -4.1408    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.5298    1.2475    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.0780    2.1262    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.5427   -5.8072    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    2.4773   -5.4631    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    0.8273   -5.1466    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.3313   -6.8090    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    6.8583    5.3111    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -1.2273   -2.2714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.7211    5.2112    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -5.2618   -2.2315    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    2.4807   -4.4603    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.0381   -4.1838    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -1.4939   -3.2509    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -2.0736   -3.3575    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -6.2613   -2.2281    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -6.7677   -1.3786    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -5.9215    1.9729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -0.8276    6.7970    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    5.0798   -3.9536    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    5.0791   -2.9513    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -3.4262    4.2984    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    3.0104    0.3471    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    3.3902    4.3050    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -1.4872    1.3200    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.2113    6.7437    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    1.1181   -1.1587    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"   -5.3984    2.3794    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    5.1259    4.3050    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							"    5.9821    6.8103    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
							" 27 42  2  0  0  0  0 \n" + 
							"  1 60  1  0  0  0  0 \n" + 
							" 32 33  2  0  0  0  0 \n" + 
							" 51 78  1  0  0  0  0 \n" + 
							" 36 37  1  0  0  0  0 \n" + 
							" 14 50  1  0  0  0  0 \n" + 
							" 11 12  1  0  0  0  0 \n" + 
							" 15 16  1  0  0  0  0 \n" + 
							" 19 20  1  0  0  0  0 \n" + 
							" 16 30  1  0  0  0  0 \n" + 
							" 16 79  2  0  0  0  0 \n" + 
							"  9 13  1  0  0  0  0 \n" + 
							" 68 69  1  0  0  0  0 \n" + 
							" 36 44  2  0  0  0  0 \n" + 
							" 26 27  1  0  0  0  0 \n" + 
							"  9 55  1  0  0  0  0 \n" + 
							" 15 61  2  0  0  0  0 \n" + 
							"  3 71  1  0  0  0  0 \n" + 
							" 40 41  1  0  0  0  0 \n" + 
							" 57 64  1  0  0  0  0 \n" + 
							" 49 62  1  0  0  0  0 \n" + 
							" 53 67  2  0  0  0  0 \n" + 
							" 54 77  1  0  0  0  0 \n" + 
							" 56 57  1  0  0  0  0 \n" + 
							" 23 53  1  0  0  0  0 \n" + 
							" 29 39  1  0  0  0  0 \n" + 
							" 56 58  1  0  0  0  0 \n" + 
							"  2  5  1  0  0  0  0 \n" + 
							" 56 59  2  0  0  0  0 \n" + 
							" 44 53  1  0  0  0  0 \n" + 
							"  3 50  2  0  0  0  0 \n" + 
							" 11 70  1  0  0  0  0 \n" + 
							"  4  7  2  0  0  0  0 \n" + 
							" 45 66  1  0  0  0  0 \n" + 
							" 19 38  1  0  0  0  0 \n" + 
							"  6 76  2  0  0  0  0 \n" + 
							" 63 68  2  0  0  0  0 \n" + 
							" 21 23  1  0  0  0  0 \n" + 
							" 19 80  2  0  0  0  0 \n" + 
							" 48 49  2  0  0  0  0 \n" + 
							" 21 22  1  0  0  0  0 \n" + 
							" 29 30  1  0  0  0  0 \n" + 
							" 31 63  1  0  0  0  0 \n" + 
							" 35 63  1  0  0  0  0 \n" + 
							" 10 77  1  0  0  0  0 \n" + 
							" 24 31  1  0  0  0  0 \n" + 
							" 34 35  2  0  0  0  0 \n" + 
							" 13 14  2  0  0  0  0 \n" + 
							" 17 18  2  0  0  0  0 \n" + 
							" 38 43  1  0  0  0  0 \n" + 
							"  2 81  2  0  0  0  0 \n" + 
							" 17 27  1  0  0  0  0 \n" + 
							"  7 13  1  0  0  0  0 \n" + 
							" 41 42  1  0  0  0  0 \n" + 
							" 45 46  2  0  0  0  0 \n" + 
							" 15 46  1  0  0  0  0 \n" + 
							" 11 43  1  0  0  0  0 \n" + 
							"  8 28  2  0  0  0  0 \n" + 
							"  6 49  1  0  0  0  0 \n" + 
							" 18 73  1  0  0  0  0 \n" + 
							" 32 64  1  0  0  0  0 \n" + 
							" 34 44  1  0  0  0  0 \n" + 
							" 24 25  2  0  0  0  0 \n" + 
							" 33 72  1  0  0  0  0 \n" + 
							" 52 74  1  0  0  0  0 \n" + 
							" 51 62  2  0  0  0  0 \n" + 
							" 20 52  2  0  0  0  0 \n" + 
							"  1  2  1  0  0  0  0 \n" + 
							" 54 55  2  0  0  0  0 \n" + 
							" 24 43  1  0  0  0  0 \n" + 
							"  3  4  1  0  0  0  0 \n" + 
							" 35 37  1  0  0  0  0 \n" + 
							" 39 40  1  0  0  0  0 \n" + 
							"  5  6  1  0  0  0  0 \n" + 
							" 48 78  1  0  0  0  0 \n" + 
							" 29 54  1  0  0  0  0 \n" + 
							" 50 51  1  0  0  0  0 \n" + 
							"  7  8  1  0  0  0  0 \n" + 
							" 42 47  1  0  0  0  0 \n" + 
							" 20 28  1  0  0  0  0 \n" + 
							"  9 10  2  0  0  0  0 \n" + 
							" 17 32  1  0  0  0  0 \n" + 
							" 41 75  2  0  0  0  0 \n" + 
							" 61 66  1  0  0  0  0 \n" + 
							" 21 45  1  0  0  0  0 \n" + 
							" 18 47  1  0  0  0  0 \n" + 
							" 22 58  1  0  0  0  0 \n" + 
							" 58 65  1  0  0  0  0 \n" + 
							" 72 73  2  0  0  0  0 \n" + 
							"  8 74  1  0  0  0  0 \n" + 
							"  1 82  2  0  0  0  0 \n" + 
							"M  END").build();

					String keyReal=Inchi.asStdInchi(cReal).getKey();
					String keyGot=Inchi.asStdInchi(c).getKey();
					assertEquals(keyReal,keyGot);
				} )});
				*/
		//cagedStructure2.png
		list.add(new Object[]{"cagedStructure2", new TestSpec("moleculeTest/cagedStructure2.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("COc1ccnc(C(=O)NC2CC3CCC2C3)c1O").build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});

		//cagedStructure3.png
		list.add(new Object[]{"cagedStructure3", new TestSpec("moleculeTest/cagedStructure3.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("COc1ccnc(C(=O)N[C@H]2CC3CC2C[C@@H]3c4ccccc4)c1O").build();

			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		
		//cagedStructure5.png
				list.add(new Object[]{"cagedStructure5", new TestSpec("moleculeTest/cagedStructure5.png", c->{
					Chemical cReal=ChemicalBuilder.createFromSmiles("COC1C2CC3=CC=C(O)C=C3C1(C)CCN2CC4CC4").build();

					String keyReal=Inchi.asStdInchi(cReal).getKey();
					String keyGot=Inchi.asStdInchi(c).getKey();
					assertEquals(keyReal,keyGot);
				} )});
		

		
		list.add(new Object[]{"nhOnTopOfEachOther", new TestSpec("moleculeTest/NHOnTopOfEachOther.png", c->{
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
//
		list.add(new Object[]{"crossColinearBonds", new TestSpec("moleculeTest/crossBonds.png", c->{
			//C1(CO1)COc2ccc(C(C)(C)c3ccc(OCC4CO4)cc3)cc2
			Chemical cReal=ChemicalBuilder.createFromSmiles("C1(CO1)COc2ccc(C(C)(C)c3ccc(OCC4CO4)cc3)cc2").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});
		//NHConnectedTogetherInChain.png
		list.add(new Object[]{"NHConnectedTogetherInChain", new TestSpec("moleculeTest/NHConnectedTogetherInChain.png", c->{

//			System.out.println("HERE!!!!!\n"+c.toMol());
			String mol =
					"\n" +
							"  CDK     01311923023D\n" +
							"\n" +
							" 50 53  0  0  0  0  0  0  0  0999 V2000\n" +
							"    2.3159    3.0415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    1.3853    2.5641    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    0.5280    3.0519    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -1.1629    3.0415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -2.0625    2.5594    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -3.0117    1.2126    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -3.6081    2.0277    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    4.0242    3.0598    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -6.0297   -0.5234    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -5.0394   -0.4970    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    5.7408   -3.2865    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    6.7291   -3.2602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    3.1409    0.5169    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    3.1509    1.5307    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    2.2662    1.0337    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    3.9957    0.0199    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    4.8704    2.5445    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -3.0316    2.8427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    5.0316   -2.5865    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -3.5286    0.3777    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -6.5154    0.3429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -2.0774    1.5456    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -6.2719    1.3120    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -6.9677    2.0277    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    2.0873   -3.3596    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    2.8129   -4.0355    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    5.2680   -1.6102    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    0.5169    4.0355    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    5.8544    2.5644    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    2.8825   -0.4373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    3.9957    4.0355    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    2.3696   -2.3795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    3.1359    2.5594    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    8.2201   -2.1470    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    7.2360   -2.3954    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    3.5783   -1.1232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -0.3181    2.5545    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -0.3280    1.5506    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -1.2027    1.0536    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    3.3496   -2.0873    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -7.5243    0.1193    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    1.6599   -1.6699    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    1.8885   -0.6660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    1.2226   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    4.0454   -2.8228    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    3.7671   -3.7969    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -8.2201    0.8150    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -7.9517    1.7792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -4.5225    0.3777    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    4.5424   -0.8946    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							" 30 36  2  0  0  0  0 \n" +
							" 32 40  2  0  0  0  0 \n" +
							"  3 37  1  0  0  0  0 \n" +
							"  5 18  1  0  0  0  0 \n" +
							" 11 12  1  0  0  0  0 \n" +
							" 36 40  1  0  0  0  0 \n" +
							" 11 19  1  0  0  0  0 \n" +
							"  6 22  1  0  0  0  0 \n" +
							" 19 27  2  0  0  0  0 \n" +
							" 43 44  1  0  0  0  0 \n" +
							" 47 48  1  0  0  0  0 \n" +
							" 24 48  1  0  0  0  0 \n" +
							" 32 42  1  0  0  0  0 \n" +
							" 25 32  1  0  0  0  0 \n" +
							"  2  3  1  0  0  0  0 \n" +
							" 26 46  1  0  0  0  0 \n" +
							" 37 38  1  0  0  0  0 \n" +
							"  4  5  1  0  0  0  0 \n" +
							" 12 35  1  0  0  0  0 \n" +
							"  6  7  1  0  0  0  0 \n" +
							" 30 43  1  0  0  0  0 \n" +
							"  5 22  2  0  0  0  0 \n" +
							" 20 49  1  0  0  0  0 \n" +
							" 21 23  1  0  0  0  0 \n" +
							"  8 17  1  0  0  0  0 \n" +
							" 25 26  1  0  0  0  0 \n" +
							" 22 39  1  0  0  0  0 \n" +
							"  9 21  1  0  0  0  0 \n" +
							" 34 35  1  0  0  0  0 \n" +
							" 13 16  2  0  0  0  0 \n" +
							" 38 39  2  0  0  0  0 \n" +
							" 13 14  1  0  0  0  0 \n" +
							" 21 41  2  0  0  0  0 \n" +
							" 13 15  2  0  0  0  0 \n" +
							" 41 47  1  0  0  0  0 \n" +
							" 10 49  1  0  0  0  0 \n" +
							" 45 46  1  0  0  0  0 \n" +
							" 19 45  1  0  0  0  0 \n" +
							"  7 18  2  0  0  0  0 \n" +
							" 14 33  1  0  0  0  0 \n" +
							"  1 33  1  0  0  0  0 \n" +
							"  1  2  1  0  0  0  0 \n" +
							"  3 28  2  0  0  0  0 \n" +
							" 36 50  1  0  0  0  0 \n" +
							"  6 20  1  0  0  0  0 \n" +
							" 17 29  1  0  0  0  0 \n" +
							" 13 30  1  0  0  0  0 \n" +
							"  9 10  1  0  0  0  0 \n" +
							" 42 43  2  0  0  0  0 \n" +
							"  4 37  2  0  0  0  0 \n" +
							" 23 24  1  0  0  0  0 \n" +
							" 33  8  1  1  0  0  0 \n" +
							"  8 31  2  0  0  0  0 \n" +
							"M  END";
			Chemical cReal=Chemical.parseMol(mol);
			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});
		//wiggleBond.png
		list.add(new Object[]{"wiggleBond", new TestSpec("moleculeTest/wiggleBond.png", c->{
			String mol =
					"\n" + 
					"\n" + 
					"\n" + 
					" 19 20  0  0  0  0  0  0  0  0999 V2000\n" + 
					"   -0.2234    0.4263    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.2859    0.4702    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -4.2501    0.4263    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.2510    0.4429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.7563   -0.4216    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.2376   -1.3028    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.2526   -1.3028    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.7538   -0.4183    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.2359    1.2998    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.2411    1.3040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.7392    0.4429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.2339   -0.4204    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.7364    1.2908    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.7306    0.4346    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.2287   -0.4245    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.7669    0.4702    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.2488    1.2908    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -2.7538    1.2975    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    4.2488   -0.4382    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1  2  1  0  0  0  0 \n" + 
					"  2  5  1  0  0  0  0 \n" + 
					"  4  3  1  6  0  0  0 \n" + 
					" 12 15  1  0  0  0  0 \n" + 
					" 11 12  2  0  0  0  0 \n" + 
					"  5  6  1  6  0  0  0 \n" + 
					" 14 16  1  0  0  0  0 \n" + 
					"  4 18  1  0  0  0  0 \n" + 
					" 14 15  2  0  0  0  0 \n" + 
					"  4  8  1  0  0  0  0 \n" + 
					"  5  8  1  0  0  0  0 \n" + 
					" 16 17  2  0  0  0  0 \n" + 
					"  8  7  1  1  0  0  0 \n" + 
					" 16 19  1  0  0  0  0 \n" + 
					" 13 18  1  0  0  0  0 \n" + 
					"  9 10  2  0  0  0  0 \n" + 
					"  9 14  1  0  0  0  0 \n" + 
					"  2 13  1  0  0  0  0 \n" + 
					"  1 11  1  0  0  0  0 \n" + 
					" 10 11  1  0  0  0  0 \n" + 
					"M  END";
			Chemical cReal=Chemical.parseMol(mol);
			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
		} )});

		//nitrogenAttachedToBondAndPlus.png
		list.add(new Object[]{"nitrogenAttachedToBondAndPlus", new TestSpec("moleculeTest/nitrogenAttachedToBondAndPlus.png", c->{

//			System.out.println("HERE!!!!!\n"+c.toMol());
			String mol =
					"\n" +
							"  CDK     02041901303D\n" +
							"\n" +
							" 24 24  0  0  0  0  0  0  0  0999 V2000\n" +
							"    3.2688    0.0989    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    3.2853    1.0966    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -4.2046   -1.0215    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0\n" +
							"    0.7916   -0.0976    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0\n" +
							"    1.4526    0.5708    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -5.0368   -1.5172    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n" +
							"    1.7080   -0.3906    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -3.3694   -1.6073    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -2.5357   -1.0215    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -3.7973   -0.0753    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    2.3840   -0.4056    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n" +
							"   -2.8661   -0.1650    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    0.0481   -0.8262    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -0.4702    1.0891    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -1.1687    0.3755    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -2.1517    0.6027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    0.5137    0.8262    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -0.8983   -0.5858    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -4.3533    0.7962    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"   -2.4230    1.6073    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    4.1346   -0.3928    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    4.1640   -1.4120    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    3.2928   -0.9163    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"    5.0353    0.0901    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n" +
							"  1  2  2  0  0  0  0 \n" +
							"  3  6  1  0  0  0  0 \n" +
							"  4  5  1  0  0  0  0 \n" +
							" 12 16  1  0  0  0  0 \n" +
							" 14 17  1  0  0  0  0 \n" +
							" 15 18  1  0  0  0  0 \n" +
							"  4 17  1  0  0  0  0 \n" +
							"  3  8  1  0  0  0  0 \n" +
							"  4  7  1  0  0  0  0 \n" +
							" 15 16  1  0  0  0  0 \n" +
							" 14 15  1  0  0  0  0 \n" +
							"  3 10  2  0  0  0  0 \n" +
							" 13 18  1  0  0  0  0 \n" +
							" 16 20  2  0  0  0  0 \n" +
							"  8  9  1  0  0  0  0 \n" +
							"  4 13  1  0  0  0  0 \n" +
							" 21 24  1  0  0  0  0 \n" +
							" 10 19  1  0  0  0  0 \n" +
							" 21 23  1  0  0  0  0 \n" +
							" 21 22  1  0  0  0  0 \n" +
							"  1 11  1  0  0  0  0 \n" +
							"  1 21  1  0  0  0  0 \n" +
							"  9 12  2  0  0  0  0 \n" +
							" 10 12  1  0  0  0  0 \n" +
							"M  CHG  1   3   1\n" +
							"M  CHG  1   4   1\n" +
							"M  CHG  1   6  -1\n" +
							"M  CHG  1  11  -1\n" +
							"M  END";
			Chemical cReal=Chemical.parseMol(mol);
			String keyReal=Inchi.asStdInchi(cReal).getKey();
			String keyGot=Inchi.asStdInchi(c).getKey();
			assertEquals(keyReal,keyGot);
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

		list.add(new Object[]{"subscriptImplicitAtomsF3Test", new TestSpec("moleculeTest/withSubscriptForF.png", c->{
			Chemical cReal=ChemicalBuilder.createFromSmiles("FC(F)(F)C1(N=N1)c2ccc(CN3C(=O)C=CC3=O)cc2").build();

			String form=c.getFormula();
			assertEquals(cReal.getFormula(),form);
		} )});

		return list;
	}


}
