package gov.nih.ncats.molvec;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Objects;
import java.util.Optional;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;

import gov.nih.ncats.molvec.internal.algo.StructureImageExtractor;

/**
 *
 */
public final class Molvec {

	private static MolvecOptions DEFAULT_OPTIONS = new MolvecOptions();
	/**
	 * Analyze the given image and try to recognize a molecular structure.
	 * @param image the image to analyze, can not be null.
	 * @return a String encoded in mol format of the recognized molecular structure.
	 * @throws IOException if there are any problems parsing the images.
	 * @throws NullPointerException if image is null.
	 */
	public static String ocr(File image) throws IOException{
		return ocr(image, DEFAULT_OPTIONS).getMolfile().get();
		
	}

	/**
	 * Analyze the given image and try to recognize a molecular structure
	 * and compute a {@link MolvecResult} using the given {@link MolvecOptions}.
	 * @param image the image to analyze, can not be null.
	 *
	 * @param options the {@link MolvecOptions} to use; if options is null, then the default options are used.
	 *
	 * @return a {@link MolvecResult} which includes a String encoded in mol format of the recognized molecular structure.
	 * @throws IOException if there are any problems parsing the images.
	 * @throws NullPointerException if image is null.
	 *
	 * @since 0.9.8
	 */
	public static MolvecResult ocr(File image, MolvecOptions options) throws IOException{
		checkNotNull(image);
		options = Optional.ofNullable(options).orElse(DEFAULT_OPTIONS);
		StructureImageExtractor sie = new StructureImageExtractor(image);
		return options.computeResult(sie.getCtab());

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
		return ocr(image, DEFAULT_OPTIONS).getMolfile().get();

	}

	/**
	 * Analyze the given image encoded data as a bytre array and try to recognize a molecular structure
	 * and compute a {@link MolvecResult} using the given {@link MolvecOptions}.
	 *
	 * @param image the image to analyze, can not be null.
	 *
	 * @param options the {@link MolvecOptions} to use; if options is null, then the default options are used.
	 *
	 * @return a {@link MolvecResult} which includes a String encoded in mol format of the recognized molecular structure.
	 *
	 * @throws IOException if there are any problems parsing the images.
	 * @throws NullPointerException if image is null.
	 *
	 * @since 0.9.8
	 */
	public static MolvecResult ocr(byte[] image, MolvecOptions options) throws IOException{
		checkNotNull(image);
		options = Optional.ofNullable(options).orElse(DEFAULT_OPTIONS);
		StructureImageExtractor sie = new StructureImageExtractor(image,options.getValues());
		
		return options.computeResult(sie.getCtab());

	}
	/**
	 * Analyze the given image and try to recognize a molecular structure.
	 * @param image the image to analyze, can not be null.
	 * @return a String encoded in mol format of the recognized molecular structure.
	 * @throws IOException if there are any problems parsing the images.
	 * @throws NullPointerException if image is null.
	 */
	public static String ocr(BufferedImage image) throws IOException{
		return ocr(image, DEFAULT_OPTIONS).getMolfile().get();
		
	}
	/**
	 * Analyze the given image and try to recognize a molecular structure
	 * and compute a {@link MolvecResult} using the given {@link MolvecOptions}.
	 *
	 * @param image the image to analyze, can not be null.
	 * @param options the {@link MolvecOptions} to use; if options is null, then the default options are used.
	 *
	 * @return a {@link MolvecResult} which includes a String encoded in mol format of the recognized molecular structure.
	 *
	 * @throws IOException if there are any problems parsing the images.
	 * @throws NullPointerException if image is null.
	 *
	 * @since 0.9.8
	 */
	public static MolvecResult ocr(BufferedImage image, MolvecOptions options) throws IOException{
		checkNotNull(image);
		options = Optional.ofNullable(options).orElse(DEFAULT_OPTIONS);
		StructureImageExtractor sie = StructureImageExtractor.createFromImage(image,options.getValues());
		
		return options.computeResult(sie.getCtab());

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

	public static CompletableFuture<MolvecResult> ocrAsync(byte[] image, MolvecOptions options){
		return CompletableFuture.supplyAsync(() -> {
			try{
				return ocr(image, options);
			}catch(Exception e){
				return MolvecResult.createFromError(e);
			}
		});
	}
	
	public static CompletableFuture<String> ocrAsync(File image){
		return CompletableFuture.supplyAsync(() -> {
			try{
				return ocr(image);
			}catch(Exception e){
//				e.printStackTrace();
				return null;
			}
		});
	}
	public static CompletableFuture<MolvecResult> ocrAsync(File image, MolvecOptions options){
		return CompletableFuture.supplyAsync(() -> {
			try{
				return ocr(image, options);
			}catch(Exception e){
				return MolvecResult.createFromError(e);
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

	public static CompletableFuture<MolvecResult> ocrAsync(File image, MolvecOptions options, Executor executor){
		return CompletableFuture.supplyAsync(() -> {
			try{
				return ocr(image, options);
			}catch(Exception e){
				return MolvecResult.createFromError(e);
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
	public static CompletableFuture<MolvecResult> ocrAsync(BufferedImage image, MolvecOptions options, Executor executor){
		return CompletableFuture.supplyAsync(() -> {
			try{
				return ocr(image, options);
			}catch(Exception e){
				return MolvecResult.createFromError(e);
			}
		},executor);
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
