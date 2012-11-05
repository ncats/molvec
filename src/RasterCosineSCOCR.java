import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.font.GlyphVector;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.Raster;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class RasterCosineSCOCR extends SCOCR {
	Set<Character> _alphabet;
	Map<Character, Collection<int[][]>> charVal = new HashMap<Character, Collection<int[][]>>();
	Map<Character, Double> WH_RATIO = new HashMap<Character, Double>();

	static int DEF_WIDTH = 20;
	static int DEF_HEIGHT = 20;
	private static BufferedImage bi = new BufferedImage(DEF_WIDTH, DEF_HEIGHT,
			BufferedImage.TYPE_BYTE_GRAY);
	static List<Font> FONT_LIST = new ArrayList<Font>();
	static {
		FONT_LIST.addAll(Arrays.asList(new Font[] {
				new Font(Font.SANS_SERIF, Font.BOLD, 8),
				new Font(Font.SANS_SERIF, Font.PLAIN, 8),
				new Font("DejaVu Sans", Font.BOLD, 8),
				new Font("Khmer OS", Font.PLAIN, 8),
				new Font("DejaVu Sans", Font.PLAIN, 8),
				new Font("Khmer OS", Font.BOLD, 8) }));

	}

	public RasterCosineSCOCR() {

	}

	/**
	 * Vertical "motion" Blur: This is a one-dimensional blurring of an image,
	 * without renormalizing.
	 * 
	 * @param bmap
	 *            bitmap as 2d int array
	 * @param rad
	 *            radius for vertical blur
	 */
	private static void vblur(int[][] bmap, int rad) {
		int sofar = 0;
		List<Integer> added = new ArrayList<Integer>();
		for (int i = 0; i < bmap.length; i++) {
			for (int j = 0; j < rad * 2 + bmap[0].length; j++) {
				// System.out.println(sofar + ":" + added.size());
				if (j < bmap.length) {
					added.add(bmap[i][j]);
					sofar += bmap[i][j];
				} else {
					if (j < bmap.length + rad) {
						added.add(bmap[i][bmap[0].length - 1]);
						sofar += bmap[i][bmap[0].length - 1];
					}
				}
				if (j >= rad) {
					if (j - rad < bmap[0].length) {
						bmap[i][j - rad] = sofar;
					}
					sofar -= added.get(0);
					added.remove(0);
				}
			}
		}
	}

	/**
	 * Horizontal "motion" Blur: This is a one-dimensional blurring of an image,
	 * without renormalizing.
	 * 
	 * @param bmap
	 *            bitmap as 2d int array
	 * @param rad
	 *            radius for horizontal blur
	 */
	private static void hblur(int[][] bmap, int rad) {
		int sofar = 0;
		List<Integer> added = new ArrayList<Integer>();
		for (int j = 0; j < bmap[0].length; j++) {
			for (int i = 0; i < rad * 2 + bmap.length; i++) {
				// System.out.println(sofar + ":" + added.size());
				if (i < bmap.length) {
					added.add(bmap[i][j]);
					sofar += bmap[i][j];
				} else {
					if (i < bmap.length + rad) {
						added.add(bmap[bmap.length - 1][j]);
						sofar += bmap[bmap.length - 1][j];
					}
				}
				if (i >= rad) {
					if (i - rad < bmap.length) {
						bmap[i - rad][j] = sofar;
					}
					sofar -= added.get(0);
					added.remove(0);
				}
			}
		}
	}

	private static void blurbmap(int[][] bmap) {
		vblur(bmap, 2);
		hblur(bmap, 2);

	}

	/**
	 * 
	 * @param c
	 *            Character to generate bitmap for
	 * @param f
	 *            Font to use in generation
	 * @param toRet
	 *            Integer array to put result into
	 * @return Returns the width/height of the character
	 */
	private static double getBitmap(Character c, Font f, int[][] toRet) {
		if (toRet == null) {
			toRet = new int[DEF_WIDTH][DEF_HEIGHT];
		}

		Graphics2D g = (Graphics2D) bi.getGraphics();
		GlyphVector gv = f.createGlyphVector(g.getFontRenderContext(),
				new char[] { c });
		Rectangle2D r = gv.getVisualBounds();

		double w = r.getWidth();
		double h = r.getHeight();

		for (int j = 0; j < DEF_HEIGHT; j++) {
			for (int i = 0; i < DEF_WIDTH; i++) {
				toRet[i][j] = (gv.getOutline().contains(r.getMinX() + (w * i)
						/ DEF_WIDTH, r.getMinY() + (h * j) / DEF_HEIGHT)) ? 1
						: 0;
			}
		}
		blurbmap(toRet);
		return w / h;
	}

	private void makeAlphabetMaps() {
		for (Font f : FONT_LIST) {
			for (char c : _alphabet) {
				int[][] charMat = new int[DEF_WIDTH][DEF_HEIGHT];
				double wh = getBitmap(c, f, charMat);
				WH_RATIO.put(c, wh);
				Collection<int[][]> bmap = charVal.get(c);
				if (bmap == null) {
					bmap = new ArrayList<int[][]>();
					charVal.put(c, bmap);
				}
				bmap.add(charMat);
			}
		}
	}

	@Override
	public void setAlphabet(Set<Character> charSet) {
		_alphabet = charSet;
		makeAlphabetMaps();
	}

	@Override
	public Set<Character> getAlphabet() {
		return _alphabet;
	}

	@Override
	public Map<Character, Number> getRanking(Raster r) {
		Map<Character, Number> myMap = new HashMap<Character, Number>();
		int[][] bmap = new int[r.getHeight()][r.getWidth()];
		for (int i = 0; i < r.getHeight(); i++) {
			r.getPixels(0, i, r.getWidth(), 1, bmap[i]);
		}
		int[][] bmap2 = new int[r.getWidth()][r.getHeight()];
		int twidth = bmap2.length;
		int theight = bmap2[0].length;
		for (int j = 0; j < theight; j++) {
			for (int i = 0; i < twidth; i++) {
				bmap2[i][j] = (bmap[j][i] > 0) ? 0 : 1;
			}
		}
		// debugPrintBmap(bmap2);
		for (char c2 : _alphabet) {
			myMap.put(c2, correlation(bmap2, c2));
		}
		return myMap;
	}

	public static void debugPrintBmap(int[][] test) {
		int twidth = test.length;
		int theight = test[0].length;

		for (int j = 0; j < theight; j++) {
			for (int i = 0; i < twidth; i++) {
				System.out.print(test[i][j]);
			}
			System.out.println();
		}

	}

	public double correlation(int[][] test, Character c) {
		double maxCor = Double.MIN_VALUE;
		for (int[][] cM : charVal.get(c)) {
			int twidth = test.length;
			int theight = test[0].length;
			int total = 0;
			int totalC = 0;
			double cor = 0;
			for (int i = 0; i < twidth; i++) {
				for (int j = 0; j < theight; j++) {
					int cx = Math.round((i * DEF_WIDTH) / twidth);
					int cy = Math.round((j * DEF_HEIGHT) / theight);
					int val = cM[cx][cy];
					cor += val * test[i][j];
					total += test[i][j] * test[i][j];
					totalC += val * val;
				}
			}
			double whrat = (double) twidth / (double) theight;
			whrat = WH_RATIO.get(c) / whrat;
			if (whrat > 1)
				whrat = 1 / whrat;
			whrat = 1 - Math.pow(1 - whrat, 2);

			double tcor = whrat * cor / (Math.sqrt(total) * Math.sqrt(totalC));
			maxCor = Math.max(tcor, maxCor);
		}
		return maxCor;
	}

}
