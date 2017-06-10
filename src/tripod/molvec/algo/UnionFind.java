// $Id: UnionFind.java 3342 2009-10-12 19:50:39Z nguyenda $
package tripod.molvec.algo;

import java.util.Map;
import java.util.List;
import java.util.ArrayList;
import java.util.TreeMap;

public class UnionFind {
    int[] nodes;
    int[] sizes;

    public UnionFind (int N) {
	nodes = new int[N];
	sizes = new int[N];
	for (int i = 0; i < N; ++i) {
	    nodes[i] = i;
	    sizes[i] = 1;
	}
    }

    protected int getRoot (int n) {
	while (n != nodes[n]) {
	    nodes[n] = nodes[nodes[n]]; // path compression
	    n = nodes[n];
	}
	return n;
    }

    public boolean find (int p, int q) {
	return getRoot (p) == getRoot (q);
    }

    public void union (int p, int q) {
	int i = getRoot (p);
	int j = getRoot (q);
	if (sizes[i] < sizes[j]) {
	    nodes[i] = j;
	    sizes[j] += sizes[i];
	}
	else {
	    nodes[j] = i;
	    sizes[i] += sizes[j];
	}
    }

    public int getComponent (int p) { return getRoot (p); }
    public int[][] getComponents () {
	Map<Integer, List<Integer>> eqmap = 
	    new TreeMap<Integer, List<Integer>>();
	for (int i = 0; i < nodes.length; ++i) {
	    List<Integer> v = eqmap.get(nodes[i]);
	    if (v == null) {
		eqmap.put(nodes[i], v = new ArrayList<Integer>());
	    }
	    v.add(i);
	}
	    
	int[][] eqv = new int[eqmap.size()][];

	int eq = 0;
	for (Map.Entry<Integer, List<Integer>> e : eqmap.entrySet()) {
	    List<Integer> v = e.getValue();
	    eqv[eq] = new int[v.size()];
	    for (int i = 0; i < v.size(); ++i) {
		eqv[eq][i] = v.get(i);
	    }
	    ++eq;
	}

	return eqv;
    }

    public String toString () {
	String s = "";
	for (int i = 0; i < nodes.length; ++i) {
	    s += String.valueOf(i+1) + ":" + nodes[i]+" ";
	}
	return s;
    }
}
