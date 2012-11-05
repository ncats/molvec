package tripod.molvec;

import java.io.Serializable;
import java.util.*;
import java.awt.geom.*;

/**
 * A zone encapsulates a collection of connected components
 * as a coherent unit such as word, line, paragraph, graphic.
 */
public class Zone implements Serializable {
    private static final long serialVersionUID = 0xa874cc00c159c356l;

    protected Area bounds = new Area ();
    protected List<Zone> children = new ArrayList<Zone>();

    public Zone () {
    }
    
    public Rectangle2D getBounds () { return bounds.getBounds2D(); }

    public void add (Zone zone) {
        bounds.add(zone.bounds);
        children.add(zone);
    }

    public boolean remove (Zone zone) {
        boolean ok = children.remove(zone);
        if (ok) {
            bounds.subtract(zone.bounds);
        }
        return ok;
    }

    public boolean contains (Zone zone) {
        return children.contains(zone);
    }

    public boolean intersects (Zone zone) {
        return bounds.intersects(zone.getBounds());
    }

    public Iterator<Zone> children () { return children.iterator(); }
    public Collection<Zone> getChildren () { 
        return Collections.unmodifiableCollection(children);
    }

    /**
     * return all children that intersects the give bounds
     */
    public Collection<Zone> getChildren (Rectangle2D bounds) {
        List<Zone> zones = new ArrayList<Zone>();
        for (Zone z : children) {
            if (z.bounds.intersects(bounds)) {
                zones.add(z);
            }
        }
        return zones;
    }

    public int getChildCount () { return children.size(); }

}
