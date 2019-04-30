package gov.nih.ncats.molvec.ui;

import java.io.IOException;
import java.io.InputStream;
import java.util.function.Consumer;
import java.util.stream.Stream;

public abstract class PrecomputedRasterCosineSCOCR extends RasterBasedCosineSCOCR{

	
	@Override
	public void getBitmapsForChar(Character c, Consumer<RasterChar> rconsumer) {
		try(InputStream ins=getInputStream(c)){
			extractRasters(ins)
			         .forEach(rconsumer);
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	
	public abstract InputStream getInputStream(Character c);
	public abstract Stream<RasterChar> extractRasters(InputStream is);
	
	
	
	
	
	

}
