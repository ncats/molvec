package gov.nih.ncats.molvec.ui;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.function.Function;

import gov.nih.ncats.molvec.image.Bitmap;
import gov.nih.ncats.molvec.algo.Tuple;

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
public interface SCOCR {

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
	
	public static Set<Character> SET_COMMON_PUCTUATION() {
		Set<Character> toRet = new LinkedHashSet<Character>();
		toRet.add('.');
		toRet.add('!');
		toRet.add('?');
		toRet.add(',');
		toRet.add('-');
		toRet.add('\'');
		toRet.add('$');
		toRet.add(')');
		toRet.add('(');
		toRet.add('"');
		toRet.add('/');
		toRet.add('\\');
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

	public static Map<Character, Number> bestOf(Map<Character, Number>... map) {
		for (char c : map[0].keySet()) {
			double largest = map[0].get(c).doubleValue();
			for (int i = 1; i < map.length; i++) {
				largest = Math.max(largest, map[i].get(c).doubleValue());
			}
			map[0].put(c, largest);
		}
		return map[0];
	}
	
	public static List<Entry<Character, Number>> sortMap(
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
	public void setAlphabet(Set<Character> charSet);
	public default void setAlphabet(String alphaString) {
		Set<Character> toRet = new LinkedHashSet<Character>();
		for (char c : alphaString.toCharArray()) {
			toRet.add(c);
		}
		setAlphabet(toRet);
	}
	
	//Getters:
	public Set<Character> getAlphabet();
	public Map<Character, Number> getRanking(Bitmap r);
	public default Entry<Character, Number> getBestMatch(Bitmap... r) {
		return getNBestMatches(1, r).get(0);
	}
	public default List<Entry<Character, Number>> getNBestMatches(int n, Bitmap... r) {
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
	
	public static class OrElseSCOCR implements SCOCR{
		private double keepCutoff =0.6;
		
		List<SCOCR> scocrList = new ArrayList<>();
		
		
		public OrElseSCOCR(List<SCOCR> scocrlist, double c){
			this.scocrList=scocrlist;
			keepCutoff=c;
		}
		
		@Override
		public void setAlphabet(Set<Character> charSet) {
			for(SCOCR s:scocrList){
				s.setAlphabet(charSet);
			}
		}

		@Override
		public Set<Character> getAlphabet() {
			for(SCOCR s:scocrList){
				return s.getAlphabet();
			}
			return null;
		}

		@Override
		public Map<Character, Number> getRanking(Bitmap r) {
			Map<Character, Number> res=new HashMap<>();
			for(SCOCR s:scocrList){
				res=s.getRanking(r);
				boolean pres=res.values()
				   .stream()
				   .filter(n->n.doubleValue()>this.keepCutoff)
				   .findAny()
				   .isPresent();
				if(pres)return res;
			}
			return res;
		}
	}
	
	public default OrElseSCOCR orElse(SCOCR backup, double cut){
		List<SCOCR> nlist = new ArrayList<>();
		nlist.add(this);
		nlist.add(backup);
		return new OrElseSCOCR(nlist,cut);
	}
	
	public default SCOCR adjustWeights(Function<Tuple<Character,Number>,Tuple<Character,Number>> transform){
		SCOCR _this=this;
		return new SCOCR(){

			@Override
			public void setAlphabet(Set<Character> charSet) {
				_this.setAlphabet(charSet);
				
			}

			@Override
			public Set<Character> getAlphabet() {
				return _this.getAlphabet();
			}

			@Override
			public Map<Character, Number> getRanking(Bitmap r) {
				Map<Character, Number> map = _this.getRanking(r);
				
				return map.entrySet()
				   .stream()
				   .map(Tuple::of)
				   .map(transform)
				   .collect(Tuple.toMap());
			}
			
		};
	}
	
	
	

}
