package tripod.molvec.algo;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.junit.Ignore;
import org.junit.Test;

import gov.nih.ncats.chemkit.api.Chemical;
import gov.nih.ncats.chemkit.api.ChemicalBuilder;
import gov.nih.ncats.chemkit.api.inchi.Inchi;
import tripod.molvec.CachedSupplier;

public class RegressionTest {

	public static enum Result{
		CORRECT_FULL_INCHI,
		CORRECT_STEREO_INSENSITIVE_INCHI,
		LARGEST_FRAGMENT_CORRECT_FULL_INCHI,
		LARGEST_FRAGMENT_CORRECT_STEREO_INSENSITIVE_INCHI,
		CORRECT,
		LARGEST_FRAGMENT_CORRECT,
		WEIRD_SOURCE,
		ATOMS_RIGHT_WRONG_BONDS,
		LARGEST_FRAGMENT_ATOM_RIGHT_WRONG_BONDS,
		ATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY,
		LARGEST_FRAGMENT_ATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY,
		RIGHT_HEVAY_ATOMS,
		RIGHT_BONDS,
		LARGEST_FRAGMENT_RIGHT_BONDS,
		INCORRECT,
		FOUND_NOTHING,
		ERROR,
		TIMEOUT
	}
	

	private File getFile(String fname){
		ClassLoader classLoader = getClass().getClassLoader();
		return new File(classLoader.getResource(fname).getFile());
		
	}
	
	
	
	public static Result testMolecule(File image, File sdf){
		
		try{
			AtomicBoolean lastWasAlias = new AtomicBoolean(false);
			String rawMol=null;
			try(Stream<String> stream = Files.lines(sdf.toPath())){
				rawMol = stream
			
					.filter(l->!l.contains("SUP"))
					.filter(l->!l.contains("SGROUP"))
					.filter(l->!l.contains("M  S"))
					.filter(l->{
						if(lastWasAlias.get()){
							lastWasAlias.set(false);
							return false;
						}
						return true;
					})
					.filter(l->{
						if(l.startsWith("A")){
							lastWasAlias.set(true);
							return false;
						}
						lastWasAlias.set(false);
						return true;
					})
					.collect(Collectors.joining("\n"));
			}
			
			Chemical c1=ChemicalBuilder.createFromMol(rawMol,Charset.defaultCharset()).build();
			System.out.println("--------------------------------");
			System.out.println(sdf.getAbsolutePath());
			System.out.println("--------------------------------");
			
			CachedSupplier<StructureImageExtractor> cget = CachedSupplier.of(()->{
				try {
					return new StructureImageExtractor(image);
				} catch (IOException e1) {
					throw new RuntimeException(e1);
				}
			});
			
			
			Thread th = new Thread(()->{
				cget.get();
			});
			
			try{
				long start = System.currentTimeMillis();
				th.start();
				while(true){
					if(cget.hasRun())break;
					Thread.sleep(1);
					if(System.currentTimeMillis()-start > 60000){
						System.out.println("OH NO TIMEOUT!\t" + image.getAbsolutePath());
						return Result.TIMEOUT;
					}
				}
			
			
			}catch(Exception e){
				throw e;
			}
			
			StructureImageExtractor sie = cget.get();
			
			Chemical c=sie.getChemical();
			
			//c1.makeHydrogensImplicit();
			//c.makeHydrogensImplicit();
			
			String iinchi=Inchi.asStdInchi(c).getKey();
			String rinchi=Inchi.asStdInchi(c1).getKey();
			
			
					
			
			int ratomCount=c1.getAtomCount();
			int iatomCount=c.getAtomCount();
			
			int rbondCount=c1.getBondCount();
			int ibondCount=c.getBondCount();
			
			
			String smilesReal=c1.toSmiles();
			String smilesFound=c.toSmiles();
			
			String formReal=c1.getFormula();
			String formFound=c.getFormula();
			
			System.out.println("Real:" + c1.toSmiles());
			System.out.println("Image:" + c.toSmiles());
			
			if(smilesFound.equals("")){
				return Result.FOUND_NOTHING;
			}
			
			if(rinchi.equals(iinchi)){
				return Result.CORRECT_FULL_INCHI;
			}
			if(rinchi.split("-")[0].equals(iinchi.split("-")[0])){
				return Result.CORRECT_STEREO_INSENSITIVE_INCHI;
			}
			
			if(formReal.equals(formFound)){
				//System.out.println("Matched!");
				return Result.CORRECT;
			}else{
				String withoutHydrogensReal =formReal.replaceAll("H[0-9]*", "");
				String withoutHydrogensFound=formFound.replaceAll("H[0-9]*", "");
				
				if(withoutHydrogensReal.equals(withoutHydrogensFound)){
					return Result.ATOMS_RIGHT_WRONG_BONDS;
				}
				
				if(smilesReal.contains(".") || smilesFound.contains(".") ){
					String largestR=Arrays.stream(smilesReal.split("[.]"))
					      .map(t->Tuple.of(t,t.length()).withVComparator())
					      .max(Comparator.naturalOrder())
					      .map(t->t.k())
					      .orElse(null);
					
					String largestF=Arrays.stream(smilesFound.split("[.]"))
					      .map(t->Tuple.of(t,t.length()).withVComparator())
					      .max(Comparator.naturalOrder())
					      .map(t->t.k())
					      .orElse(null);
					Chemical clargestR=ChemicalBuilder.createFromSmiles(largestR).build();
					Chemical clargestF=ChemicalBuilder.createFromSmiles(largestF).build();

					
					iinchi=Inchi.asStdInchi(clargestF).getKey();
					rinchi=Inchi.asStdInchi(clargestR).getKey();
					
					if(rinchi.equals(iinchi)){
						return Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI;
					}
					if(rinchi.split("-")[0].equals(iinchi.split("-")[0])){
						return Result.LARGEST_FRAGMENT_CORRECT_STEREO_INSENSITIVE_INCHI;
					}
					
					if(clargestR.getFormula().equals(clargestF.getFormula())){
						return Result.LARGEST_FRAGMENT_CORRECT;
					}
					String fragwHReal =clargestR.getFormula().replaceAll("H[0-9]*", "");
					String fragwHFound=clargestF.getFormula().replaceAll("H[0-9]*", "");
					
					if(fragwHReal.equals(fragwHFound)){
						return Result.LARGEST_FRAGMENT_ATOM_RIGHT_WRONG_BONDS;
					}
					clargestR.makeHydrogensImplicit();
					clargestF.makeHydrogensImplicit();
					
					int fratomCount=clargestR.getAtomCount();
					int fiatomCount=clargestF.getAtomCount();
					
					int frbondCount=clargestR.getBondCount();
					int fibondCount=clargestF.getBondCount();
					if(fratomCount==fiatomCount && frbondCount == fibondCount){
						return Result.LARGEST_FRAGMENT_ATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY;
					}
				}
				if(smilesReal.contains("|")){
					return Result.WEIRD_SOURCE;
				}
				if(ratomCount==iatomCount && rbondCount == ibondCount){
					return Result.ATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY;
				}
				
				//System.out.println("NO MATCH!");
				return Result.INCORRECT;
			}
		}catch(Exception e){
			e.printStackTrace();
		}
		return Result.ERROR;
	}
	
	
	//@Ignore
	@Test
	public void test1(){
		File dir1 = getFile("regressionTest/uspto");
		
		try {
			ChemicalBuilder cb = ChemicalBuilder.createFromSmiles("CCCC");
			String ii = Inchi.asStdInchi(cb.build()).getKey();
			System.out.println(ii);
			Thread.sleep(1000);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		Arrays.stream(dir1.listFiles())
		      .filter(f->f.getName().contains("."))
		      //.filter(f->f.getName().contains("2008058707_41_chem"))
		      .map(f->Tuple.of(f.getName().split("[.]")[0],f))
		      .collect(Tuple.toGroupedMap())
		      .values()
		      .stream()
		      .filter(l->l.size()==2)
		      
		      .map(l->{
		    	  	if(!l.get(1).getName().toLowerCase().endsWith("tif") && !l.get(1).getName().toLowerCase().endsWith("png")){
		    	  		List<File> flist = new ArrayList<File>();
		    	  		flist.add(l.get(1));
		    	  		flist.add(l.get(0));
		    	  		return flist;
					}
		    	  	return l;
		      })
//		      .collect(shuffler(new Random(11111124l)))		      
//		      .limit(100)

//NOTE, I THINK THIS TECHNICALLY WORKS, BUT SINCE THERE IS PARALLEL THINGS GOING ON IN EACH, IT SOMETIMES WILL STARVE A CASE FOR A LONG TIME
//		      .parallel()
		      
		      
		      .map(fl->Tuple.of(fl,testMolecule(fl.get(1),fl.get(0))))
		      .map(t->t.swap())
		      .peek(t->System.out.println(t.k()))
		      
		      .collect(Tuple.toGroupedMap())
		      .entrySet()
		      .stream()
		      .map(Tuple::of)
		      .forEach(t->{
		    	  Result r=t.k();
		    	  List<List<File>> fl = t.v();
		    	  System.out.println("======================================");
		    	  System.out.println(r.toString() + "\t" + fl.size());
		    	  System.out.println("--------------------------------------");
		    	  System.out.println(fl.stream().map(f->f.get(1).getAbsolutePath()).collect(Collectors.joining("\n")));
		      });
	}
	
	public static <T> Collector<T,List<T>,Stream<T>> shuffler(Random r){
		
		return new Collector<T,List<T>,Stream<T>>(){

			@Override
			public BiConsumer<List<T>, T> accumulator() {
				return (l,t)->{
					l.add(t);
				};
			}

			@Override
			public Set<java.util.stream.Collector.Characteristics> characteristics() {
				//java.util.stream.Collector.Characteristics.
				return new HashSet<java.util.stream.Collector.Characteristics>();
			}

			@Override
			public BinaryOperator<List<T>> combiner() {
				return (l1,l2)->{
					l1.addAll(l2);
					return l1;
				};
			}

			@Override
			public Function<List<T>, Stream<T>> finisher() {
				return (u)->{
					Collections.shuffle(u,r);
					return u.stream();
				};
			}

			@Override
			public Supplier<List<T>> supplier() {
				return ()-> new ArrayList<T>();
			}

		};
	}
	
	
}
