package gov.nih.ncats.molvec;

import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.*;
public class Java11RegressionTest {

    @Test
    public void canReadFileCorrectly() throws IOException {
        File f = new File(getClass().getResource("/moleculeTest/circleAromatic.png").getFile());
        String result = Molvec.ocr(f);
//        System.out.println(result);

        assertMolFileNotBlank(result);
    }

    private void assertMolFileNotBlank(String result) {
        int numLines = result.split("\n").length;
        assertTrue("blank mol file? only " + numLines + " lines long", numLines > 5);
    }

    @Test
    public void grayscaleImageWithoutAllBands() throws IOException{
        File f = new File(getClass().getResource("/convertedGrayscale.png").getFile());
        String result = Molvec.ocr(f);
        System.out.println(result);
        assertMolFileNotBlank(result);
    }
}
