package tripod.molvec.algo;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.Test;

public class SegmentationTest {

	
	private File getFile(String fname){
		ClassLoader classLoader = getClass().getClassLoader();
		return new File(classLoader.getResource(fname).getFile());
		
	}
	
	@Test
	public void parallelLinesAtProblemAngleShouldBeJoined() throws Exception {
		File f=getFile("segmentTest/badparlines.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		assertEquals(6,sie.getLinesJoined().size());	
	}
	
	@Test
	public void thickWigglyLineShouldBeJoined() throws Exception {
		File f=getFile("segmentTest/angled_wiggle_line.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		assertEquals(1,sie.getLinesJoined().size());	
	}
	
	@Test
	public void thickWigglyLinesShouldBeJoinedAndStitched() throws Exception {
		File f=getFile("segmentTest/angled_wiggle_line_branch.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		int num=sie.getLinesJoined().size();
		

		assertTrue("Expected number of segments was 3 or 4, found:" + num,num>=3 && num<=4);
		if(num!=3){
			int finalEdges=sie.getCtab().getEdges().size();
			assertTrue("Expected 3 final edges, found:" + finalEdges,finalEdges==3);
		}
	}
	
	@Test
	public void thickWigglyLinesWithDoubleBondsShouldBeJoinedAndStitched() throws Exception {
		File f=getFile("segmentTest/angled_wiggle_line_double.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		int num=sie.getLinesJoined().size();
		

		assertTrue("Expected number of segments was between 6 and 8, found:" + num,num>=6 && num<=8);
		if(num!=6){
			int finalEdges=sie.getCtab().getEdges().size();
			assertTrue("Expected 4 final edges, found:" + finalEdges,finalEdges==4);
		}
	}

	@Test
	public void thickWigglyLinesWithSmallAngleSectionShouldBeStitched() throws Exception {
		File f=getFile("segmentTest/angled_wiggle_line_with_small_part.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		int num=sie.getLinesJoined().size();
		assertTrue("Expected number of segments was between 2 and 3, found:" + num,num>=2 && num<=3);
	}
	
	@Test
	public void thinLinesWithDoubleBondsShouldBeStitchedToCorrectNumber() throws Exception {
		File f=getFile("segmentTest/thin_double.png");
		StructureImageExtractor sie = new StructureImageExtractor();
		sie.load(f);
		int num=sie.getLinesJoined().size();
		assertTrue("Expected number of segments was between 9 and 11, found:" + num,num>=9 && num<=11);
		if(num!=9){
			int finalEdges=sie.getCtab().getEdges().size();
			assertTrue("Expected 7 final edges, found:" + finalEdges,finalEdges==7);
		}
	}
	
}