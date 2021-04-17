package gov.nih.ncats.molvec;

import org.junit.Test;

import java.awt.geom.Area;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import static org.junit.Assert.*;
public class MolvecTest {
    private static final String lineSep = System.lineSeparator();

    @Test
    public void getBoundingBox() throws Exception{
        File f = new File(MolvecTest.class.getResource("/moleculeTest/tylenol.png").getFile());
        MolvecResult result = Molvec.ocr(f, new MolvecOptions());
        Rectangle2D boundingBox = result.getOriginalBoundingBox().get();

        String mol=result.getMolfile().get();
        assertFalse(result.hasError());
        assertTrue(mol, mol.contains("11 11  0  0  0  0  0  0  0  0999 V2000"));
        assertFalse(mol, mol.contains("$$$$"));
        assertTrue(new Area(boundingBox).contains(100,70,300, 100));

    }

    @Test
    public void asSdNoProperties() throws Exception{
        File f = new File(MolvecTest.class.getResource("/moleculeTest/tylenol.png").getFile());
        MolvecResult result = Molvec.ocr(f, new MolvecOptions());

        String mol=result.getSDfile().get();
        assertFalse(result.hasError());
        assertTrue(mol, mol.contains("11 11  0  0  0  0  0  0  0  0999 V2000"));
        assertTrue(mol, mol.contains("$$$$"));

    }

    @Test
    public void asSdWithProperties() throws Exception{
        File f = new File(MolvecTest.class.getResource("/moleculeTest/tylenol.png").getFile());
        MolvecResult result = Molvec.ocr(f, new MolvecOptions());
        Map<String,String> properties = new LinkedHashMap<>();
        properties.put("foo", "my foo value");
        properties.put("bar", "my bar value");

        String mol=result.getSDfile(properties).get();
        assertFalse(result.hasError());
        assertTrue(mol, mol.contains("11 11  0  0  0  0  0  0  0  0999 V2000"));
        assertTrue(mol, mol.contains(">  <foo>"+lineSep + "my foo value" +lineSep + lineSep +
                ">  <bar>"+lineSep + "my bar value" +lineSep + lineSep));
        assertTrue(mol, mol.contains("$$$$"));

    }
}
