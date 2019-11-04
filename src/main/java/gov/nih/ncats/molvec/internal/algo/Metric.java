package gov.nih.ncats.molvec.internal.algo;

public interface Metric<T> {
    double evaluate (T arg0, T arg1);
}
