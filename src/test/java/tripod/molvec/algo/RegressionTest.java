package tripod.molvec.algo;

import java.io.File;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.Ignore;
import org.junit.Test;

import gov.nih.ncats.chemkit.api.Chemical;
import gov.nih.ncats.chemkit.api.ChemicalBuilder;

public class RegressionTest {

	public static enum Result{
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
		ERROR
	}
	

	private File getFile(String fname){
		ClassLoader classLoader = getClass().getClassLoader();
		return new File(classLoader.getResource(fname).getFile());
		
	}
	
	
	
	public static Result testMolecule(File image, File sdf){
		
		try{
			String rawMol = Files.lines(sdf.toPath())
					.filter(l->!l.contains("SUP"))
					.filter(l->!l.contains("SGROUP"))
					.collect(Collectors.joining("\n"));
			
			Chemical c1=ChemicalBuilder.createFromMol(rawMol,Charset.defaultCharset()).build();
			System.out.println("--------------------------------");
			System.out.println(sdf.getAbsolutePath());
			System.out.println("--------------------------------");
			
			StructureImageExtractor sie = new StructureImageExtractor();
			sie.load(image);
			Chemical c=sie.getChemical();
			
			c1.makeHydrogensImplicit();
			c.makeHydrogensImplicit();
			
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
					if(clargestR.getFormula().equals(clargestF.getFormula())){
						System.out.println("LARGEST MATCHED!");
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
	
	
	@Ignore
	@Test
	public void test1(){
		File dir1 = getFile("regressionTest/testSet1");
		
		Arrays.stream(dir1.listFiles())
		      .filter(f->f.getName().contains("."))
		      //.filter(f->f.getName().contains("2008058707_41_chem"))
		      .map(f->Tuple.of(f.getName().split("[.]")[0],f))
		      .collect(Tuple.toGroupedMap())
		      .values()
		      .stream()
		      .filter(l->l.size()==2)
		      
		      .map(l->{
		    	  	if(!l.get(0).getName().endsWith("sdf")){
		    	  		List<File> flist = new ArrayList<File>();
		    	  		flist.add(l.get(1));
		    	  		flist.add(l.get(0));
		    	  		return flist;
					}
		    	  	return l;
		      })
		      .map(fl->Tuple.of(fl,testMolecule(fl.get(1),fl.get(0))))
		      .peek(t->{
		    	  System.out.println(t.v());
		      })
		      .map(t->t.swap())
		      .collect(Tuple.toGroupedMap())
		      .forEach((r,fl)->{
		    	  System.out.println("======================================");
		    	  System.out.println(r.toString() + "\t" + fl.size());
		    	  System.out.println("--------------------------------------");
		    	  System.out.println(fl.stream().map(f->f.get(1).getAbsolutePath()).collect(Collectors.joining("\n")));
		      });
	}
	
	
}
