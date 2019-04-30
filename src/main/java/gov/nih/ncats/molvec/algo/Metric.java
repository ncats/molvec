package gov.nih.ncats.molvec.algo;

public interface Metric<T> {
    double evaluate (T arg0, T arg1);
}
