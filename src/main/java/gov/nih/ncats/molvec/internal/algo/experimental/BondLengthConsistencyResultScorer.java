package gov.nih.ncats.molvec.internal.algo.experimental;

import java.util.DoubleSummaryStatistics;

import gov.nih.ncats.molwitch.Chemical;

public class BondLengthConsistencyResultScorer implements ResultScorer{
	private static double RATIO_RANGE=0.10;

	public BondLengthConsistencyResultScorer(){}

	@Override
	public double score(Chemical c) {

		DoubleSummaryStatistics dss=c.bonds()
				.mapToDouble(bb->bb.getBondLength())
				.summaryStatistics();

		double avg=dss.getAverage();

		long count=c.bonds()
				.mapToDouble(bb->bb.getBondLength())
				.filter(d->d<avg*(1+RATIO_RANGE))
				.filter(d->d>avg*(1-RATIO_RANGE))
				.count();
		;

		return count/(0.0+c.getBondCount());

	}

}
