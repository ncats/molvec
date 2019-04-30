package gov.nih.ncats.molvec.ui;

import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.font.GlyphVector;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;

public class FontBasedRasterCosineSCOCR extends RasterBasedCosineSCOCR {

	private List<Font> fontList = new ArrayList<Font>();
	
	public static List<Font> DEFAULT_FONTS(){
		return SANS_SERIF_FONTS();
	}
	
	public static List<Font> SERIF_FONTS(){
		return Arrays.asList(new Font[] {
				new Font(Font.SERIF, Font.PLAIN, 8),
				new Font(Font.SERIF, Font.BOLD, 8) 
				});
	}
	public static List<Font> SANS_SERIF_FONTS(){
		return Arrays.asList(new Font[] {
				new Font(Font.SANS_SERIF, Font.BOLD, 8),
				new Font(Font.SANS_SERIF, Font.PLAIN, 8),
				new Font("DejaVu Sans", Font.BOLD, 8),
				new Font("Khmer OS", Font.PLAIN, 8),
				new Font("DejaVu Sans", Font.PLAIN, 8),
				new Font("Khmer OS", Font.BOLD, 8) 
				});
	}
	
	public FontBasedRasterCosineSCOCR(List<Font> flist) {
		fontList.addAll(flist);
	}
	
	public FontBasedRasterCosineSCOCR() {
		this(Arrays.asList(new Font[] {
				new Font(Font.SANS_SERIF, Font.BOLD, 8),
				new Font(Font.SANS_SERIF, Font.PLAIN, 8),
				new Font("DejaVu Sans", Font.BOLD, 8),
				new Font("Khmer OS", Font.PLAIN, 8),
				new Font("DejaVu Sans", Font.PLAIN, 8),
				new Font("Khmer OS", Font.BOLD, 8), 
				new Font(Font.SERIF, Font.PLAIN, 8),
				new Font(Font.SERIF, Font.BOLD, 8) 
				}));
	}
	
	
	@Override
	public void getBitmapsForChar(Character c, Consumer<RasterChar> rconsumer) {
		for (Font f : fontList) {
			int[][] charMat = new int[DEF_WIDTH][DEF_HEIGHT];
			Rectangle2D whrat = getBitmap(c,f,charMat);
			rconsumer.accept(new RasterChar(charMat,whrat));
		}
	}
	
	
	/**
	 * 
	 * @param c
	 *            Character to generate bitmap for
	 * @param f
	 *            Font to use in generation
	 * @param toRet
	 *            Integer array to put result into
	 * @return Returns the Rectangle2D bounding box of the font glyph
	 */
	private static Rectangle2D getBitmap(Character c, Font f, int[][] toRet) {
		if (toRet == null) {
			toRet = new int[DEF_WIDTH][DEF_HEIGHT];
		}

		Graphics2D g = (Graphics2D) bi.getGraphics();
		GlyphVector gv = f.createGlyphVector(g.getFontRenderContext(),
				new char[] { c });
		Shape outline=gv.getOutline();
		Rectangle2D r = gv.getVisualBounds();

		double w = r.getWidth();
		double h = r.getHeight();

		Point2D testPoint = new Point2D.Double();
		
		for (int j = 0; j < DEF_HEIGHT; j++) {
			for (int i = 0; i < DEF_WIDTH; i++) {
				testPoint.setLocation(r.getMinX() + (w * i) / DEF_WIDTH, 
									  r.getMinY() + (h * j) / DEF_HEIGHT);
				toRet[i][j] = outline.contains(testPoint) ? 1
						: 0;
			}
		}
		
		
		        
		blurbmap(toRet);
		return r;
	}

}
