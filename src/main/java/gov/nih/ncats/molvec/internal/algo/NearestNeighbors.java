package gov.nih.ncats.molvec.internal.algo;

import java.util.*;
import java.util.logging.Logger;

/**
 * K-Nearest neighbor class
 */
public class NearestNeighbors<T> {
    private static final Logger logger = 
        Logger.getLogger(NearestNeighbors.class.getName());

    private static final int PAD = 10;
    private static final int DEFAULT_MAX_NEIGHBORS = 10;

    public static class Neighbor<E> implements Comparable<Neighbor<E>> {
        double value;
        E neighbor;

        Neighbor (double value, E neighbor) {
            this.value = value;
            this.neighbor = neighbor;
        }

        public int compareTo (Neighbor<E> nb) {
            if (value < nb.value) return -1;
            if (value > nb.value) return 1;
            return 0;
        }

        public double getValue () { return value; }
        public E getNeighbor () { return neighbor; }
    }

    protected Map<T, Queue<Neighbor<T>>> neighbors = 
        new HashMap<T, Queue<Neighbor<T>>>();
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
        Queue<Neighbor<T>> nq = new PriorityQueue<Neighbor<T>>();
        for (Map.Entry<T, Queue<Neighbor<T>>> me : neighbors.entrySet()) {
            double xv = metric.evaluate(me.getKey(), entry);
            if (xv < threshold) {
                Queue<Neighbor<T>> cq = me.getValue();
                cq.add(new Neighbor<T> (xv, entry));
                nq.add(new Neighbor<T> (xv, me.getKey()));

                adjustNeighborQueue (cq);
                adjustNeighborQueue (nq);
            }
        }
        neighbors.put(entry, nq);
    }

    protected void adjustNeighborQueue (Queue<Neighbor<T>> q) {
        // truncate the queue to fit
        if (q.size() > (maxNb+PAD)) {
            List<Neighbor<T>> entries = new ArrayList<Neighbor<T>>();
            for (int i = 0; i < maxNb && !q.isEmpty(); ++i) {
                entries.add(q.poll());
            }
            q.clear();
            q.addAll(entries);
        }
    }

    public T neighbor (T entry) {
        Queue<Neighbor<T>> q = neighbors.get(entry);
        if (q != null) {
            Neighbor<T> nb = q.peek();
            return nb.getNeighbor();
        }
        return null;
    }

    public T nearest (T entry) {
        Queue<Neighbor<T>> q = neighbors.get(entry);
        if (q != null) {
            Neighbor<T> nb = q.peek();
            return nb != null ? nb.getNeighbor() : null;
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
        Queue<Neighbor<T>> q = neighbors.get(entry);
        if (q != null) {
            int size = Math.min(K, q.size());
            Iterator<Neighbor<T>> iter = q.iterator();
            for (int k = 0; k < size; ++k) {
                Neighbor<T> nb = iter.next();
                nbs.add(nb.getNeighbor());
            }
        }
        return nbs;
    }

    public List<Neighbor<T>> neighborList (T entry) {
        return neighborList (entry, maxNb);
    }

    public List<Neighbor<T>> neighborList (T entry, int K) {
        if (K < 1) {
            throw new IllegalArgumentException ("K must > 0");
        }

        List<Neighbor<T>> nbs = new ArrayList<Neighbor<T>>(K);
        Queue<Neighbor<T>> q = neighbors.get(entry);
        if (q != null) {
            int size = Math.min(K, q.size());
            Iterator<Neighbor<T>> iter = q.iterator();
            for (int k = 0; k < size; ++k) {
                nbs.add(iter.next());
            }
        }
        return nbs;
    }

    public List<T> neighbors (T entry, double cutoff) {
        List<T> nbs = new ArrayList<T>();
        Queue<Neighbor<T>> q = neighbors.get(entry);
        if (q != null) {
            Iterator<Neighbor<T>> iter = q.iterator();
            while (iter.hasNext()) {
                Neighbor<T> nb = iter.next();
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

    public List<Neighbor<T>> neighborList (T entry, double cutoff) {
        List<Neighbor<T>> nbs = new ArrayList<Neighbor<T>>();
        Queue<Neighbor<T>> q = neighbors.get(entry);
        if (q != null) {
            Iterator<Neighbor<T>> iter = q.iterator();
            for (int k = 0; iter.hasNext() && k < maxNb; ++k) {
                Neighbor<T> nb = iter.next();
                if (nb.getValue() <= cutoff) {
                    nbs.add(nb);
                }
                else {
                    break;
                }
            }
        }
        return nbs;
    }

    public Set<T> entries () { return neighbors.keySet(); }

    public void clear () { neighbors.clear(); }
    public int size () { return neighbors.size(); }
}
