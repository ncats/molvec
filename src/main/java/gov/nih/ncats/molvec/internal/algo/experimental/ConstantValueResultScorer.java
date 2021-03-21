package gov.nih.ncats.molvec.internal.algo.experimental;

import gov.nih.ncats.molwitch.Chemical;

public class ConstantValueResultScorer implements ResultScorer{
	double con=0;
	public ConstantValueResultScorer(double v){
		this.con=v;
	}
	@Override
	public double score(Chemical c) {
		return con;
	}

}
