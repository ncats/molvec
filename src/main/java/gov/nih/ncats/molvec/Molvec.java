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


	

}
