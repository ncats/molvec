package tripod.molvec.algo;

import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;

public class NearestNeighbors<T> {
    private static final Logger logger = 
        Logger.getLogger(NearestNeighbors.class.getName());

    private static final int PAD = 10;
    private static final int DEFAULT_MAX_NEIGHBORS = 10;

    public class Entry implements Comparable<Entry> {
        double value;
        T neighbor;

        Entry (double value, T neighbor) {
            this.value = value;
            this.neighbor = neighbor;
        }

        public int compareTo (Entry nb) {
            if (value < nb.value) return -1;
            if (value > nb.value) return 1;
            return 0;
        }

        public double getValue () { return value; }
        public T getNeighbor () { return neighbor; }
    }

    protected Map<T, Queue<Entry>> neighbors = new HashMap<T, Queue<Entry>>();
    protected Metric<T> metric;
    protected double threshold = Double.MAX_VALUE;
    protected int maxNb = DEFAULT_MAX_NEIGHBORS;

    public NearestNeighbors (Metric<T> metric) {
        if (metric == null) {
            throw new IllegalArgumentException ("Metric can't be null");
        }

        this.metric = metric;
    }

    public NearestNeighbors (int maxNb, Metric<T> metric) {
        if (maxNb <= 0) {
            throw new IllegalArgumentException ("Max neighbors must > 0");
        }

        if (metric == null) {
            throw new IllegalArgumentException ("Metric can't be null");
        }

        this.maxNb = maxNb;
        this.metric = metric;
    }

    public NearestNeighbors (double threshold, Metric<T> metric) {
        if (threshold <= 0.) {
            throw new IllegalArgumentException ("Threshold must be > 0.");
        }

        if (metric == null) {
            throw new IllegalArgumentException ("Metric can't be null");
        }

        this.metric = metric;
        this.threshold = threshold;
    }

    public void setMaxNeighbors (int maxNb) {
        this.maxNb = maxNb;
    }
    public int getMaxNeighbors () { return maxNb; }

    public void setThreshold (double threshold) {
        if (threshold <= 0.) {
            throw new IllegalArgumentException ("Threshold must be > 0.");
        }
        this.threshold = threshold;
    }
    public double getThreshold () { return threshold; }

    public void addAll (Collection<T> entries) {
        for (T e : entries) {
            add (e);
        }
    }

    public void add (T entry) {
        Queue<Entry> nq = new PriorityQueue<Entry>();
        for (Map.Entry<T, Queue<Entry>> me : neighbors.entrySet()) {
            double xv = metric.evaluate(me.getKey(), entry);
            if (xv < threshold) {
                Queue<Entry> cq = me.getValue();
                cq.add(new Entry (xv, entry));
                nq.add(new Entry (xv, me.getKey()));

                adjustNeighborQueue (cq);
                adjustNeighborQueue (nq);
            }
        }
        neighbors.put(entry, nq);
    }

    protected void adjustNeighborQueue (Queue<Entry> q) {
        // truncate the queue to fit
        if (q.size() > (maxNb+PAD)) {
            List<Entry> entries = new ArrayList<Entry>();
            for (int i = 0; i < maxNb && !q.isEmpty(); ++i) {
                entries.add(q.poll());
            }
            q.clear();
            q.addAll(entries);
        }
    }

    public T neighbor (T entry) {
        Queue<Entry> q = neighbors.get(entry);
        if (q != null) {
            Entry nb = q.peek();
            return nb.getNeighbor();
        }
        return null;
    }

    public List<T> neighbors (T entry) {
        return neighbors (entry, maxNb);
    }

    public List<T> neighbors (T entry, int K) {
        if (K < 1) {
            throw new IllegalArgumentException ("K must > 0");
        }

        List<T> nbs = new ArrayList<T>(K);
        Queue<Entry> q = neighbors.get(entry);
        if (q != null) {
            int size = Math.min(K, q.size());
            Iterator<Entry> iter = q.iterator();
            for (int k = 0; k < size; ++k) {
                Entry nb = iter.next();
                nbs.add(nb.getNeighbor());
            }
        }
        return nbs;
    }

    public List<T> neighbors (T entry, double cutoff) {
        List<T> nbs = new ArrayList<T>();
        Queue<Entry> q = neighbors.get(entry);
        if (q != null) {
            Iterator<Entry> iter = q.iterator();
            while (iter.hasNext()) {
                Entry nb = iter.next();
                if (nb.getValue() <= cutoff) {
                    nbs.add(nb.getNeighbor());
                }
                else {
                    break;
                }
            }
        }
        return nbs;
    }

    public void clear () { neighbors.clear(); }
    public int size () { return neighbors.size(); }
}
