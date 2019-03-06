package tripod.molvec;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;
import java.util.concurrent.FutureTask;

import gov.nih.ncats.chemkit.api.Chemical;
import tripod.molvec.algo.StructureImageExtractor;

public class Molvec {

	public static Chemical ocr(File image) throws Exception{
		StructureImageExtractor sie = new StructureImageExtractor(image);
		return sie.getChemical();
		
	}
	
	public static Chemical ocr(BufferedImage image) throws Exception{
		StructureImageExtractor sie = StructureImageExtractor.createFromImage(image);
		return sie.getChemical();
		
	}
	
	public static CompletableFuture<Chemical> ocrAsync(File image){
		return CompletableFuture.supplyAsync(() -> {
			try{
				return ocr(image);
			}catch(Exception e){
				e.printStackTrace();
				return null;
			}
		});
	}
	public static CompletableFuture<Chemical> ocrAsync(File image, Executor executor){
		return CompletableFuture.supplyAsync(() -> {
			try{
				return ocr(image);
			}catch(Exception e){
				e.printStackTrace();
				return null;
			}
		},executor);
	}
	

	public static CompletableFuture<Chemical> ocrAsync(BufferedImage image){
		return CompletableFuture.supplyAsync(() -> {
			try{
				return ocr(image);
			}catch(Exception e){
				e.printStackTrace();
				return null;
			}
		});
	}
	public static CompletableFuture<Chemical> ocrAsync(BufferedImage image, Executor executor){
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
