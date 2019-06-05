package gov.nih.ncats.molvec.image.binarization;

public class ImageStats {
	public double min;
	public double max;
	public double stdev;
	public double mean;
	public double threshold;
	public int[] histogram;
	public int[] histogramRaw;
	public double count;

	public double getPercentageThreshold(){
		return 100*(threshold-min)/(max-min);
	}
}