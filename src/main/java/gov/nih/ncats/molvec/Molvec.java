package gov.nih.ncats.molvec;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.PrintStream;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;

import gov.nih.ncats.molvec.algo.StructureImageExtractor;

public class Molvec {

	public static String ocr(File image) throws Exception{
		StructureImageExtractor sie = new StructureImageExtractor(image);
		String mol = sie.getCtab().toMol();
		return mol;
		
	}

	public static String ocr(byte[] image) throws Exception{
		StructureImageExtractor sie = new StructureImageExtractor(image);
		return sie.getCtab().toMol();

	}
	public static String ocr(BufferedImage image) throws Exception{
		StructureImageExtractor sie = StructureImageExtractor.createFromImage(image);
		return sie.getCtab().toMol();
		
	}
	public static CompletableFuture<String> ocrAsync(byte[] image){
		return CompletableFuture.supplyAsync(() -> {
			try{
				return ocr(image);
			}catch(Exception e){
				e.printStackTrace();
				return null;
			}
		});
	}
	
	public static CompletableFuture<String> ocrAsync(File image){
		return CompletableFuture.supplyAsync(() -> {
			try{
				return ocr(image);
			}catch(Exception e){
				e.printStackTrace();
				return null;
			}
		});
	}
	public static CompletableFuture<String> ocrAsync(File image, Executor executor){
		return CompletableFuture.supplyAsync(() -> {
			try{
				return ocr(image);
			}catch(Exception e){
				e.printStackTrace();
				return null;
			}
		},executor);
	}
	

	public static CompletableFuture<String> ocrAsync(BufferedImage image){
		return CompletableFuture.supplyAsync(() -> {
			try{
				return ocr(image);
			}catch(Exception e){
				e.printStackTrace();
				return null;
			}
		});
	}
	public static CompletableFuture<String> ocrAsync(BufferedImage image, Executor executor){
		return CompletableFuture.supplyAsync(() -> {
			try{
				return ocr(image);
			}catch(Exception e){
				e.printStackTrace();
				return null;
			}
		},executor);
	}

	public static String getVersion(){
		return "0.9.2";
	}
	
	private static void printUsage(){
		System.err.println("MolVec -- Java library for converting structure images into molfiles\n" + 
				"\n" + 
				"USAGE:\n" + 
				"\n" + 
				"	java -jar molvec.jar <image_filename> [-o <output_filename>]\n" + 
				"\n" + 
				"Note that the output of this command will write to standard out if no -o argument is specified.\n" + 
				"\n" + 
				"");
	}
	
	public static void main(String[] args) throws Exception{
		//TODO: CLI stuff
		if(args.length==0){
			System.err.println("ERROR -- No file specified\n\n");
			printUsage();
			return;
		}
		
		PrintStream os = System.out;
		
		String outf = null;
		String inf = null;
	
		for(int i=0;i<args.length;i++){
			if(args[i].equals("-o")){
				outf=args[i+1];
				i++;
			}else{
				inf = args[i];
			}
		}
		if(inf == null){
			System.err.println("ERROR -- No file specified\n\n");
			printUsage();
			return;
		}
		
		File f = new File(inf);
		if(!f.exists()){
			System.err.println("ERROR -- file \"" + f.getAbsolutePath() + "\" does not exist!\n\n");
			printUsage();
			return;
		}
	
		if(outf!=null){
			os = new PrintStream(new File(outf));
		}
		
		try{
			os.println(ocr(f));	
		}finally{
			os.close();
		}
			
			
		
	}
}
