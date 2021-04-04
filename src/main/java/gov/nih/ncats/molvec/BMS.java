package gov.nih.ncats.molvec;

import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.io.*;
import java.nio.file.Files;
import java.util.concurrent.atomic.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

import gov.nih.ncats.molvec.internal.algo.experimental.ModifiedMolvecPipeline;
import gov.nih.ncats.molwitch.*;
import gov.nih.ncats.molwitch.inchi.Inchi;
import gov.nih.ncats.molvec.internal.algo.experimental.ChemFixer;

public class BMS {
    static final Logger logger = Logger.getLogger(BMS.class.getName());
    static final String DEFAULT_INCHI = "InChI=1S/C12H24N2O/c1-13-12-4-2-11(3-5-12)10-14-6-8-15-9-7-14/h11-13H,2-10H2,1H3";

    public static int calculateLD (String x, String y) {
        int[][] dp = new int[x.length() + 1][y.length() + 1];
        
        for (int i = 0; i <= x.length(); i++) {
            for (int j = 0; j <= y.length(); j++) {
                if (i == 0) {
                    dp[i][j] = j;
                }
                else if (j == 0) {
                    dp[i][j] = i;
                }
                else {
                    dp[i][j] = min(dp[i - 1][j - 1]
                                   + costOfSubstitution(x.charAt(i - 1), y.charAt(j - 1)),
                                   dp[i - 1][j] + 1, dp[i][j - 1] + 1);
                }
            }
        }

        return dp[x.length()][y.length()];
    }

    public static int costOfSubstitution (char a, char b) {
        return a == b ? 0 : 1;
    }
    
    public static int min (int... numbers) {
        return Arrays.stream(numbers).min().orElse(Integer.MAX_VALUE);
    }
    
    static class Result {
        String id;
        String inchi;
        String type;
	double score;
        int distance;
        String truth;

        Result (String id, String truth) {
            this.id = id;
            this.truth = "InChI=1S/H2O/h1H2".equals(truth) ? DEFAULT_INCHI : truth;
        }

        String getImage () {
            return id.charAt(0)+File.separator+id.charAt(1)+File.separator+id.charAt(2)
                +File.separator+id+".png";
        }

        int calcLD () {
            return distance = calculateLD (inchi, truth);
        }

        public String toString () {
            return id+"\t"+inchi+"\t"+type+"\t"
		+String.format("%1$.3f", score)+"\t"+distance+"\t"+truth;
        }
    }

    static class Trainer {
	ConcurrentHashMap<String, Integer> matched = new ConcurrentHashMap<>();
	ConcurrentHashMap<String, Integer> notmatched = new ConcurrentHashMap<>();
	AtomicInteger total = new AtomicInteger ();
	AtomicInteger matchedCount = new AtomicInteger ();

	public void update (Chemical c, Result r) {
	    if (r.distance > 0) {
		update (c, notmatched);
	    }
	    else {
		update (c, matched);
		matchedCount.incrementAndGet();
	    }
	    total.incrementAndGet();
	}

	public boolean isEmpty () { return total.get() == 0; }
	public void write (String file) throws Exception {
	    write (new File (file));
	}
	public void write (File file) throws Exception {
	    try (PrintStream ps = new PrintStream (new FileOutputStream (file))) {
		ps.println(total.get()+"\t"+matchedCount.get());
		Set<String> features = new HashSet<>();
		features.addAll(matched.keySet());
		features.addAll(notmatched.keySet());
		for (String f : features) {
		    ps.print(f+"\t");
		    Integer c = matched.get(f);
		    if (c == null) c = 0;
		    ps.print(c+"\t");
		    c = notmatched.get(f);
		    if (c == null) c = 0;
		    ps.println(c);
		}
	    }
	}

	void update (Chemical c, ConcurrentHashMap<String, Integer> accum) {
	    c.atoms()
		.map(a -> ChemFixer.atFeat(a, 1))
		.collect(Collectors.toSet())
		.stream().forEach(f -> accum.merge(f, 1, (x,y)-> x+y));
	}
    }
    
    public static void main (String[] argv) throws Exception {
        if (argv.length < 2) {
            logger.info("Usage: "+BMS.class.getName()+" PATH FILE...");
            System.exit(1);
        }

        File inDir = new File (argv[0]);
        if (!inDir.isDirectory()) {
            logger.log(Level.SEVERE, argv[0]+": Not a directory!");
            System.exit(1);
        }

	Trainer trainer = new Trainer ();

	long now = System.currentTimeMillis();
	ModifiedMolvecPipeline.setInChIDefaultScorer(); //loads the ikeys.txt file from resources
	logger.info("InChI keys loaded in "+String.format("%1$.3fs!",
							  1e-3*(System.currentTimeMillis()-now)));
	//this is optional, but if you want to make it a smaller footprint you can play with this number
	ModifiedMolvecPipeline.setTryLimit(23);
	//ModifiedMolvecPipeline.setTrySkipAndLimit(20, 15);
	ModifiedMolvecPipeline.DO_MULTI_TRIES=true;

        AtomicInteger count = new AtomicInteger ();
	AtomicLong ts = new AtomicLong (System.currentTimeMillis());
        for (int i = 1; i < argv.length; ++i) {
            File f = new File (argv[i]);
            if (!f.exists()) {
                logger.warning(argv[i]+": File not found!");
                continue;
            }
            Files.lines(f.toPath())
                .map(l -> {
                        int pos = l.indexOf(','), offset = 0;
                        String id = l.substring(0, pos);
                        if (l.charAt(++pos) == '"')
                            offset = 1;
                        return new Result (id, l.substring(pos+offset, l.length()-offset));
                    })
		//.limit(10000)                   // TP ADDED
		.collect(Collectors.toList())   // TP ADDED
                .parallelStream()               // TP ADDED
		//.parallel()
                .forEach(r -> {
			int cnt = count.incrementAndGet();
			if (cnt % 100 == 0) {
			    long current = System.currentTimeMillis();
			    logger.info(String.format("%1$ 6d...%2$.3fs ",
						      cnt, 1e-3*(current-ts.get()))+r.id);
			    ts.set(current);
			}
                        try {
                            ModifiedMolvecPipeline.ModifiedMolvecResult mvr =
				ModifiedMolvecPipeline.process
                                (new File (inDir, r.getImage()), new MolvecOptions());
                            Chemical mc = Chemical.parse(mvr.getMolfile().get());
                            r.inchi = mc.toInchi().getInchi();
                            r.type = mvr.getType();
			    r.score = mvr.getScore();
			    r.calcLD();
			    if (argv[0].startsWith("train")) {
				trainer.update(Inchi.toChemical(r.truth), r);
			    }
                        }
                        catch (Exception e1) {
                            r.inchi = DEFAULT_INCHI;
			    r.calcLD();
                        }
                        System.out.println(r);
                    });
        }

	if (!trainer.isEmpty()) {
	    trainer.write("trainer.txt");
	}
    }
}
