package gov.nih.ncats.molvec.internal.algo.experimental;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import gov.nih.ncats.molwitch.Chemical;

public class InChIKeySetScorer implements ResultScorer{
	private Set<Long> ikeys;
	
	public InChIKeySetScorer(File iKeysFile){
		try(Stream<String> sf= Files.lines((iKeysFile).toPath())){
			ikeys = sf.map(l->encodeKey(l))
    			.collect(Collectors.toCollection(()->new HashSet<>(100000000)));
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	private boolean hasIKey(String ikey){
		return this.ikeys.contains(encodeKey(ikey));
	}
	
	private static long encodeKey(String s){
		//stereo-insensitive
		s=s.split("-")[0];
		
		int x=s.hashCode();
		int y=(s.substring(2)+"!").hashCode()^x;
		long l = (((long)x) << 32) | (y & 0xffffffffL);
		return l;
	}

	@Override
	public double score(Chemical c) {
		try {
			String ikey=c.toInchi().getKey();
			if(hasIKey(ikey)){
				return 1;
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
		return 0;
	}

}
