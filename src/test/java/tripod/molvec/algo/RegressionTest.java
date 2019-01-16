package tripod.molvec.algo;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import org.junit.Ignore;
import org.junit.Test;

import gov.nih.ncats.chemkit.api.Chemical;
import gov.nih.ncats.chemkit.api.ChemicalBuilder;

public class RegressionTest {

	public static enum Result{
		CORRECT,
		LARGEST_FRAGMENT_CORRECT,
		WEIRD_SOURCE,
		RIGHT_HEVAY_ATOMS,
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
			Chemical c1=ChemicalBuilder.createFromMol(sdf).build();
			System.out.println("--------------------------------");
			System.out.println(sdf.getAbsolutePath());
			System.out.println("--------------------------------");
			
			StructureImageExtractor sie = new StructureImageExtractor();
			sie.load(image);
			Chemical c=sie.getChemical();
			
			String smilesReal=c1.toSmiles();
			String smilesFound=c.toSmiles();
			
			System.out.println("Real:" + c1.toSmiles());
			System.out.println("Image:" + c.toSmiles());
			
			
			if(c1.getFormula().equals(c.getFormula())){
				System.out.println("Matched!");
				return Result.CORRECT;
			}else{
				if(smilesFound.equals("")){
					return Result.FOUND_NOTHING;
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
				}
				if(smilesReal.contains("|")){
					return Result.WEIRD_SOURCE;
				}
				
				System.out.println("NO MATCH!");
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
		    	  System.out.println("--------------------------------------");
		    	  System.out.println("RESULT:" + r.toString());
		    	  System.out.println("COUNT:" + fl.size());
		      });
	}
	
	
}
