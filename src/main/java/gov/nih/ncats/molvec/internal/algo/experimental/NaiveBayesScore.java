package gov.nih.ncats.molvec.internal.algo.experimental;

import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.util.stream.Collectors;

import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.inchi.Inchi;

public class NaiveBayesScore {
    static final Logger logger = Logger.getLogger(NaiveBayesScore.class.getName());
    static final int MINCOUNT = 20;
    
    /**
     * load in features file of the form:
99994 61985
C(-C,-S,-S)	19	17
C(-C,-Cl,=O)	29	33
C(-N,-P,=N)	1	2
C(-O,-S,=S)	6	3
...
    */
    final int total; // number of samples spanned by the features
    final int matchedCount;
    final double matchedWt;
    final double notmatchedWt;
    final double prior;
    
    final List<String> features = new ArrayList<>(); // selected list of features to use
    final ConcurrentHashMap<String, Integer> matched = new ConcurrentHashMap<>();
    final ConcurrentHashMap<String, Integer> notmatched = new ConcurrentHashMap<>();
    
    public NaiveBayesScore (File file) throws Exception {
	logger.info("Loading features '"+file+"'...");
	try (BufferedReader br = new BufferedReader (new FileReader (file))) {
	    String[] toks = br.readLine().split("\t");
	    total = Integer.parseInt(toks[0]);
	    matchedCount = Integer.parseInt(toks[1]);
	    int lines = 1;
	    for (String line; (line = br.readLine()) != null; ++lines) {
		toks = line.split("\t");
		if (toks.length == 3) {
		    String f = toks[0];
		    int mc = 0;
		    if (!"".equals(toks[1]))
			mc = Integer.parseInt(toks[1]);
		    matched.put(f, mc+1);

		    int nc = 0;
		    if (!"".equals(toks[2])) 
			nc = Integer.parseInt(toks[2]);
		    notmatched.put(f, nc+1);
		    features.add(f);
		}
	    }
	    Collections.sort(features, (a, b) -> (matched.get(b)+notmatched.get(b)) -
			     (matched.get(a)+notmatched.get(a)));
	    /*
	    for (String f : features) {
		System.out.println(f+"\t"+matched.get(f)+"\t"+notmatched.get(f));
	    }
	    */
	    // laplace smoothing
	    matchedWt = matchedCount + 2.0;
	    notmatchedWt = (total - matchedCount) + 2.0;
	    prior = (double)matchedCount / total; // empirical prior
	    logger.info(lines+" features loaded!");
	}
    }

    public double score (Chemical c) {
	Set<String> features = c.atoms()
	    .map(a -> ChemFixer.atFeat(a, 1))
	    .collect(Collectors.toSet());
	return score (features);
    }

    public double score (Set<String> fset) {
	double pm = 0., pn = 0.;
	for (String f : features) {
	    if (fset.contains(f)) {
		Integer mc = matched.get(f);
		Integer nc = notmatched.get(f);
		pm += Math.log(mc / matchedWt);
		pn += Math.log(nc / notmatchedWt);
	    }
	}
	pm = Math.exp(pm);
	pn = Math.exp(pn);

	double posterior = pm * prior / (pm*prior + pn*(1 - prior));

	return posterior;
    }

    public static void main (String[] argv) throws Exception {
	if (argv.length < 2) {
	    logger.info("usage: NaiveBayesScore FILE INCHI...");
	    System.exit(1);
	}
	NaiveBayesScore nbs = new NaiveBayesScore (new File (argv[0]));
	for (int i = 1; i < argv.length; ++i) {
	    logger.info(argv[i]+"\n===> posterior = "+ nbs.score(Inchi.toChemical(argv[i])));
	}
    }
}
