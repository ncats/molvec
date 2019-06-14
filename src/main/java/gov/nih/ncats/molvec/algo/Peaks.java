package gov.nih.ncats.molvec.algo;

import java.util.*;
import java.util.logging.Logger;

/**
 * A simple peak detection class
 */
public class Peaks {
    private static final Logger logger = 
        Logger.getLogger(Peaks.class.getName());

    static final int DEFAULT_WINDOW = 5;
    static final double DEFAULT_STD_FACTOR = 1.;

    static class Peak implements Comparable<Peak> {
        int location;
        double value;

        Peak (int location, double value) {
            this.location = location;
            this.value = value;
        }

        public int compareTo (Peak p) {
            if (p.value > value) return 1;
            if (p.value < value) return -1;
            return location - p.location;
        }
    }
    
    private int window; // averaging window
    private double psr; // peak to signal ratio

    public Peaks () {
        this (DEFAULT_WINDOW, DEFAULT_STD_FACTOR);
    }

    public Peaks (int window) {
        this.window = window;
    }

    public Peaks (double psr) {
        this.psr = psr;
    }

    public Peaks (int window, double psr) {
        this.window = window;
        this.psr = psr;
    }

    public int getWindow () { return window; }
    public void setWindow (int window) { this.window = window; }

    public double getPeakSignalRatio () { return psr; }
    public void setPeakSignalRatio (double psr) { this.psr = psr; }

    public int[] detect (int[] signal) {
        if (signal.length <= window || signal.length == 0) {
            return new int[0];
        }
        
        double[] smoothed = new double[signal.length];
        double avg = 0, left, right;

        // first, smooth out the signal
        for (int i = 0, j, k; i < signal.length; ++i) {
            // left
            k = Math.max(i - window, 0);
            left = 0.;
            for (j = k; j < i; ++j) 
                left += signal[j];
            if (j > k)
                left /= (j - k);

            // right
            k = Math.min(i+window+1, signal.length);
            right = 0.;
            for (j = i +1; j < k; ++j)
                right += signal[j];
            if (j > i+1)
                right /= (j - i -1);

            smoothed[i] = ((signal[i] - left) + (signal[i] - right))/2D;
            avg += smoothed[i];
        }

        avg /= signal.length;
        
        double std = 0.;
        for (int i = 0; i < smoothed.length; ++i) {
            double x = smoothed[i] - avg;
            std += x*x;
        }
        std = Math.sqrt(std/smoothed.length);

        double thres = psr * std;

        // now locate peak candidates
        List<Peak> peaks = new ArrayList<>(smoothed.length);
        for (int i = 0; i < smoothed.length; ++i) {
            if ((smoothed[i] - avg) >= thres) {
                peaks.add(new Peak (i, smoothed[i]));
            }
        }

//        try {
//            PrintStream ps = new PrintStream
//                (new FileOutputStream ("smoothed.txt"));
//            for (int i = 0; i < smoothed.length; ++i) {
//                if (smoothed[i] > 0) {
//                    ps.println(i+" "+smoothed[i]);
//                }
//            }
//            ps.close();
//        }
//        catch (Exception ex) {
//            ex.printStackTrace();
//        }

        // sort them in descending order
        Collections.sort(peaks);

        // extract peaks
        BitSet mask = new BitSet (signal.length);

        Iterator<Peak> iter= peaks.iterator();
        while(iter.hasNext()){
            Peak p = iter.next();
            if (mask.get(p.location)) {
               iter.remove();
            }
            else {
                // mask the left
                for (int i = Math.max(p.location-window, 0); 
                     i < p.location; ++i) {
                    mask.set(i);
                }
                // mask the right
                int end = Math.min(p.location+window, signal.length);
                for (int i = p.location; i < end; ++i) {
                    mask.set(i);
                }
            }
        }

        int[] locs = new int[peaks.size()];
        int i=0;
        for(Peak p : peaks){
            locs[i++] = p.location;
        }

        return locs;
    }
}