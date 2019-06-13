package gov.nih.ncats.molvec.ui;

import gov.nih.ncats.molvec.CachedSupplier;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.BiConsumer;
import java.util.function.Consumer;

/**
 * Created by katzelda on 6/13/19.
 */
public abstract class AbstractStupidestPossibleSCOCR extends RasterBasedCosineSCOCR {
    private final Map<String, List<RasterChar>> readAhead = new HashMap<>();

    AbstractStupidestPossibleSCOCR(){
        loadRasters( (k, v)-> readAhead.computeIfAbsent(k, newKey -> new ArrayList<>()).add(v));

    }
    @Override
    public void getBitmapsForChar(Character c, Consumer<RasterChar> rconsumer) {
//        if(readAhead.isEmpty()){
//            loadRasters( (k, v)-> readAhead.computeIfAbsent(k, newKey -> new ArrayList<>()).add(v));
//        }
        List<RasterChar> list = readAhead.get(c.toString());
        if(list !=null){
            list.forEach(rconsumer);
        }
    }

    protected abstract void loadRasters(BiConsumer<String, RasterChar> consumer);
}
