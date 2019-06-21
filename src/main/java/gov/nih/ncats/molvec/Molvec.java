package gov.nih.ncats.molvec;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Objects;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;

import gov.nih.ncats.molvec.algo.StructureImageExtractor;

/**
 *
 */
public final class Molvec {

	/**
	 * Analyze the given image and try to recognize a molecular structure.
	 * @param image the image to analyze, can not be null.
	 * @return a String encoded in mol format of the recognized molecular structure.
	 * @throws IOException if there are any problems parsing the images.
	 * @throws NullPointerException if image is null.
	 */
	public static String ocr(File image) throws IOException{
		checkNotNull(image);
		StructureImageExtractor sie = new StructureImageExtractor(image);
		String mol = sie.getCtab().toMol();
		return mol;
		
	}

	private static void checkNotNull(Object obj){
		Objects.requireNonNull(obj, "image can not be null");
	}
	/**
	 * Analyze the given image encoded data as a bytre array and try to recognize a molecular structure.
	 * @param image the image to analyze, can not be null.
	 * @return a String encoded in mol format of the recognized molecular structure.
	 * @throws IOException if there are any problems parsing the images.
	 * @throws NullPointerException if image is null.
	 */
	public static String ocr(byte[] image) throws IOException{
		checkNotNull(image);
		StructureImageExtractor sie = new StructureImageExtractor(image);
			return sie.getCtab().toMol();

	}
	/**
	 * Analyze the given image and try to recognize a molecular structure.
	 * @param image the image to analyze, can not be null.
	 * @return a String encoded in mol format of the recognized molecular structure.
	 * @throws IOException if there are any problems parsing the images.
	 * @throws NullPointerException if image is null.
	 */
	public static String ocr(BufferedImage image) throws IOException{
		checkNotNull(image);
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
