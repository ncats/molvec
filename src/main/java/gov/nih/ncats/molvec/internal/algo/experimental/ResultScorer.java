package gov.nih.ncats.molvec.internal.algo.experimental;

import gov.nih.ncats.molwitch.Chemical;

public interface ResultScorer {
	public double score(Chemical c);
	
	public default ResultScorer ifBelow(double s, ResultScorer alt){
		return (c)->{
			double d=this.score(c);
			if(d<s){
				return alt.score(c);
			}
			return d;
		};
	}
	public default ResultScorer scale(double scale){
		return (c)->{
			double m=this.score(c);
			return m*scale;
		};
	}
}
