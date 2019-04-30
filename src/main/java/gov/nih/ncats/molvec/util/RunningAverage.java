package gov.nih.ncats.molvec.util;

/**
 * Created by katzelda on 3/11/19.
 */
public final class RunningAverage {

    private final double defaultValue;
    private double current;
    private int num;

    public RunningAverage(){
        this(0);
    }

    public RunningAverage(double defaultValue){
        this.defaultValue = defaultValue;
    }
    public void add(double value){
        current+=value;
        num++;
    }

    public double computeAvg(){
        if(num==0){
            return defaultValue;
        }
        return current/num;
    }
}
