package gov.nih.ncats.molvec.internal.algo.experimental;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.BinaryOperator;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;

import gov.nih.ncats.common.Tuple;
import gov.nih.ncats.molwitch.Chemical;

public class InChIKeySetScorer implements ResultScorer{
	private ConcurrentHashMap<Long,Long> ikeys;

	private static BinaryOperator<Long> bin = (a,b)->a|b;

	public InChIKeySetScorer(File iKeysFile, boolean compressed){
		ikeys=new ConcurrentHashMap<Long,Long>(119803351);


		if(compressed){
			try(InputStream fis= new FileInputStream(iKeysFile)){
				try(InputStream in = new GZIPInputStream(fis)){
					new BufferedReader(new InputStreamReader(in,StandardCharsets.UTF_8)).lines().parallel()
					.map(l->Tuple.of(encodeKey(l), encodeStereoKey(l)))

					.forEach(t->{
						ikeys.merge(t.k(), t.v(), bin);
					});	

					;

				}

			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}else{
			try(Stream<String> sf= Files.lines((iKeysFile).toPath())){
				sf.parallel()
				.map(l->Tuple.of(encodeKey(l), encodeStereoKey(l)))

				.forEach(t->{
					ikeys.merge(t.k(), t.v(), bin);
				});	
			}catch(Exception e){
				e.printStackTrace();
			}
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
