package gov.nih.ncats.molvec.util;

import java.util.Comparator;

public class CompareUtil {
	public static <T> Comparator<T> naturalOrder(){
		return (a,b)->{
			if(a instanceof Comparable){
				Comparable c=(Comparable)a;
				return c.compareTo(b);
			}
			if(a==null && b==null)return 0;
			if(a==null)return -1;
			if(b==null)return 1;
			throw new IllegalStateException("Don't know how to compare these things");
		};
	}
}
