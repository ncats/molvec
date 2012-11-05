import java.awt.image.Raster;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 * @author peryeata
 * 
 *         An SCOCR (Single Character Optical Character Recognition) object is
 *         one that can be passed a raster image, and attempt to determine the
 *         most likely characters that the bitmap represents, along with some
 *         numerical confidence level for each guess.
 * 
 *         Extending classes should allow setting the "alphabet", and implement
 *         returning a Map from each character in the alphabet to a confidence
 *         number.
 * 
 * 
 */
public abstract class SCOCR {

	public static Set<Character> SET_ALPHA_UPPER() {
		Set<Character> toRet = new LinkedHashSet<Character>();
		for (char c = 'A'; c <= 'Z'; c++) {
			toRet.add(c);
		}
		return toRet;
	}

	public static Set<Character> SET_ALPHA_LOWER() {
		Set<Character> toRet = new LinkedHashSet<Character>();
		for (char c = 'a'; c <= 'z'; c++) {
			toRet.add(c);
		}
		return toRet;
	}

	public static Set<Character> SET_NUMERIC() {
		Set<Character> toRet = new LinkedHashSet<Character>();
		for (char c = '0'; c <= '9'; c++) {
			toRet.add(c);
		}
		return toRet;
	}

	public static Set<Character> SET_ALPHANUMERIC() {
		Set<Character> toRet = new LinkedHashSet<Character>();
		toRet.addAll(SET_NUMERIC());
		toRet.addAll(SET_ALPHA_UPPER());
		toRet.addAll(SET_ALPHA_LOWER());
		return toRet;
	}

	public static Set<Character> SET_PERIODIC_TABLE_CHARS() {
		String ptable = "AaBbCcDdEeFfGgHhIiKkLlMmNnOoPpRrSsTtUuVvWXYyZ";
		Set<Character> toRet = new LinkedHashSet<Character>();

		for (char c : ptable.toCharArray()) {
			toRet.add(c);
		}
		return toRet;
	}

	public static Set<Character> SET_COMMON_CHEM_EXTRAS() {
		String ptable = "+-";
		Set<Character> toRet = new LinkedHashSet<Character>();
		for (char c : ptable.toCharArray()) {
			toRet.add(c);
		}
		return toRet;
	}

	public static Set<Character> SET_COMMON_CHEM_ALL() {
		Set<Character> toRet = new LinkedHashSet<Character>();
		toRet.addAll(SET_NUMERIC());
		toRet.addAll(SET_COMMON_CHEM_EXTRAS());
		toRet.addAll(SET_PERIODIC_TABLE_CHARS());
		return toRet;
	}

	private static Map<Character, Number> bestOf(Map<Character, Number>... map) {
		for (char c : map[0].keySet()) {
			double largest = map[0].get(c).doubleValue();
			for (int i = 1; i < map.length; i++) {
				largest = Math.max(largest, map[i].get(c).doubleValue());
			}
			map[0].put(c, largest);
		}
		return map[0];
	}
	private static List<Entry<Character, Number>> sortMap(
			Map<Character, Number> rmap) {
		List<Entry<Character, Number>> ranks = new ArrayList<Entry<Character, Number>>(
				rmap.entrySet());
		Collections.sort(ranks, new Comparator<Entry<Character, Number>>() {
			@Override
			public int compare(Entry<Character, Number> arg0,
					Entry<Character, Number> arg1) {
				return -Double.compare(arg0.getValue().doubleValue(), arg1
						.getValue().doubleValue());
			}
		});
		return ranks;
	}

	//Setters:
	public abstract void setAlphabet(Set<Character> charSet);
	public void setAlphabet(String alphaString) {
		Set<Character> toRet = new LinkedHashSet<Character>();
		for (char c : alphaString.toCharArray()) {
			toRet.add(c);
		}
		setAlphabet(toRet);
	}
	
	//Getters:
	public abstract Set<Character> getAlphabet();
	public abstract Map<Character, Number> getRanking(Raster r);
	public Entry<Character, Number> getBestMatch(Raster... r) {
		return getNBestMatches(1, r).get(0);
	}
	public List<Entry<Character, Number>> getNBestMatches(int n, Raster... r) {
		Map[] g = new Map[r.length];
		for (int i = 0; i < r.length; i++) {
			g[i] = getRanking(r[i]);
		}
		List<Entry<Character, Number>> retRanks = new ArrayList<Entry<Character, Number>>();
		List<Entry<Character, Number>> ranks = sortMap(bestOf(g));
		for (int i = 0; i < n; i++) {
			retRanks.add(ranks.get(i));
		}
		return retRanks;
	}


}
