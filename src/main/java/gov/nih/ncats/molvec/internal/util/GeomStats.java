package gov.nih.ncats.molvec.internal.util;

import java.awt.Shape;
import java.awt.Rectangle;

import java.util.*;
import java.util.logging.Logger;


public class GeomStats {
    private static final Logger logger = 
        Logger.getLogger (GeomStats.class.getName ());

    public double avgW, avgH, medW, medH, medA;
    public double stdW, stdH, varW, varH;
    public double avgA, stdA, varA; // area

    public int minW, maxW, minH, maxH, minA, maxA;
    public int[] histW, histH, histA;

    public GeomStats (Collection<Shape> objs) {
        if (objs == null || objs.isEmpty()) {
            throw new IllegalArgumentException ("Objects is null or empty");
        }

        maxW = maxH = maxA = 0;
        minW = minH = minA = Integer.MAX_VALUE;
        avgW = avgH = 0.;
        int size = objs.size();

        int[] ws = new int[size];
        int[] hs = new int[size];
        int[] as = new int[size];

        {  int i = 0;
            for (Shape s : objs) {
                Rectangle r = s.getBounds();
                if (r.width > maxW) maxW = r.width;
                if (r.width < minW) minW = r.width;
                if (r.height > maxH) maxH = r.height;
                if (r.height < minH) minH = r.height;

                avgW += r.width;
                avgH += r.height;
                ws[i] = r.width;
                hs[i] = r.height;
                int area = r.width*r.height;
                if (area > maxA) maxA = area;
                if (area < minA) minA = area;
                as[i] = area;
                avgA += area;
                ++i;
            }
        }
        avgW /= size;
        avgH /= size;
        avgA /= size;

        Arrays.sort(ws);
        Arrays.sort(hs);
        Arrays.sort(as);
        if (size % 2 == 0) {
            medW = (ws[size/2] + ws[size/2 -1])/2.;
            medH = (hs[size/2] + hs[size/2 -1])/2.;
            medA = (as[size/2] + as[size/2 -1])/2.;
        }
        else {
            medW = ws[size/2];
            medH = hs[size/2];
            medA = as[size/2];
        }

        histW = new int[maxW+1];
        histH = new int[maxH+1];
        histA = new int[maxA+1];
        // do another pass to calc std and hist
        varW = varH = varA = 0.;
        for (int i = 0; i < ws.length; ++i) {
            ++histW[ws[i]];
            ++histH[hs[i]];
            ++histA[as[i]];
            
            double w = ws[i] - avgW;
            double h = hs[i] - avgH;
            double a = as[i] - avgA;
            varW += w*w;
            varH += h*h;
            varA += a*a;
        }
        varW /= size;
        varH /= size;
        varA /= size;
        stdW = Math.sqrt(varW);
        stdH = Math.sqrt(varH);
        stdA = Math.sqrt(varA);

        hs = null;
        ws = null;
        as = null;
    }

    public String toString () {
        return "{avgW="+String.format("%1$.1f",avgW)+",avgH="
            +String.format("%1$.1f",avgH)+",medW="
            +String.format("%1$.1f",medW)+",medH="
            +String.format("%1$.1f",medH)+",stdW="
            +String.format("%1$.2f",stdW)+",stdH="
            +String.format("%1$.2f",stdH)+",avgA="
            +String.format("%1$.1f",avgA)+",medA="
            +String.format("%1$.1f",medA)+",stdA="
            +String.format("%1$.1f",stdA)
            +",minW="+minW+",maxW="+maxW
            +",minH="+minH+",maxH="+maxH
            +",minA="+minA+",maxA="+maxA
            +"}";
    }
}
