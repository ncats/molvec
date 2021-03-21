package gov.nih.ncats.molvec.internal.algo.experimental;

import java.util.DoubleSummaryStatistics;

import gov.nih.ncats.molwitch.Chemical;

public class FeatureRarityResultScorer implements ResultScorer{
	

	public FeatureRarityResultScorer(){}

	@Override
	public double score(Chemical c) {
		return Math.random();

	}

}
