package tripod.molvec.algo;

import java.awt.Shape;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import tripod.molvec.Bitmap;
import tripod.molvec.CachedSupplier;
import tripod.molvec.util.GeomUtil;

public class ImageDecomposition {
	private Bitmap bitmap;
	private Bitmap thin;
	
	private CachedSupplier<Collection<Shape>> polygons = CachedSupplier.of(()->bitmap.connectedComponents(Bitmap.Bbox.Polygon));
	private CachedSupplier<Collection<Path2D>> segments = CachedSupplier.of(()->thin.segments());
	private CachedSupplier<Collection<Line2D>> lines = CachedSupplier.of(()->GeomUtil.asLines(segments.get()));
	
	
	
	
	private Map<Shape,List<Entry<Character,Number>>> ocrAttmept = new HashMap<Shape,List<Entry<Character,Number>>>();
	
	
	public ImageDecomposition(Bitmap bm){
		this.bitmap=bm;
		this.thin=bitmap.thin();
	}
	
	public ImageDecomposition(File f) throws IOException{
		this(Bitmap.read(f));
	}
	
	
	
	
	
	

}
