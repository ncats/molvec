package gov.nih.ncats.molvec.internal.algo.experimental;

import gov.nih.ncats.molwitch.Chemical;

public class InChILengthScorer implements ResultScorer{
	
    int minLength = 0;

	public InChILengthScorer(int l){
	    this.minLength=l;
	}

	@Override
	public double score(Chemical c) {
	    try {
    	    String inch=c.toInchi().getInchi();
    	    if(inch.length()<minLength)return 0;
    		return 1;
	    }catch(Exception e) {
	        return 0;
	    }

	}

}
