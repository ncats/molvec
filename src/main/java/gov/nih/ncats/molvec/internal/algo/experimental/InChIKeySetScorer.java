package gov.nih.ncats.molvec.internal.algo.experimental;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.BinaryOperator;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import gov.nih.ncats.common.Tuple;
import gov.nih.ncats.molwitch.Chemical;

public class InChIKeySetScorer implements ResultScorer{
	private ConcurrentHashMap<Long,Long> ikeys;
	
	public InChIKeySetScorer(File iKeysFile){
		try(Stream<String> sf= Files.lines((iKeysFile).toPath())){
			BinaryOperator<Long> bin = (a,b)->a|b;
			
			
		 
		 
			ikeys = sf.map(l->Tuple.of(encodeKey(l), encodeStereoKey(l)))
					
    			.collect(Collectors.toConcurrentMap(t->t.k(),t->t.v(), bin, ()-> new ConcurrentHashMap<Long,Long>(100000000)));
    			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	private boolean hasIKey(String ikey){
		return this.ikeys.containsKey(encodeKey(ikey));
	}
	
	private boolean hasIKeyStereo(String ikey){
		return (this.ikeys.getOrDefault(encodeKey(ikey),0l) & encodeStereoKey(ikey))!=0;
	}
	
	private static long encodeKey(String s){
		//stereo-insensitive
		s=s.split("-")[0];
		
		int x=s.hashCode();
		int y=(s.substring(2)+"!").hashCode()^x;
		long l = (((long)x) << 32) | (y & 0xffffffffL);
		return l;
	}
	
	private static long encodeStereoKey(String s){
		//stereo-sensitive
		s=s.split("-")[1];
		
		int x=s.hashCode();
		int y=(s.substring(2)+"!").hashCode()^x;
		long l = (((long)x) << 32) | (y & 0xffffffffL);
		l=Math.abs(l);
		
		return (1<<((int) (l/64)));
	}

	@Override
	public double score(Chemical c) {
		try {
			String ikey=c.toInchi().getKey();
			if(hasIKey(ikey)){
				if(hasIKeyStereo(ikey)){
					return 1;
				}
				return 0.95;
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
//			e.printStackTrace();
		}		
		return 0;
	}

}
