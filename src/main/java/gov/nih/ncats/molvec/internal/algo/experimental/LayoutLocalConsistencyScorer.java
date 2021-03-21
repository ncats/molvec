package gov.nih.ncats.molvec.internal.algo.experimental;

import java.util.Optional;
import java.util.stream.IntStream;

import gov.nih.ncats.common.Tuple;
import gov.nih.ncats.common.stream.StreamUtil;
import gov.nih.ncats.molvec.internal.util.GeomUtil;
import gov.nih.ncats.molvec.internal.util.GeomUtil.ShapeWrapper;
import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.MolwitchException;

public class LayoutLocalConsistencyScorer implements ResultScorer{

	public LayoutLocalConsistencyScorer(){}

	public Optional<ShapeWrapper> getNeighborShape(Atom at){
		ShapeWrapper sw= null;
		try{
			sw= ShapeWrapper.of(StreamUtil.with(at.getNeighbors().stream())
					.and(at)
					.stream()
					.map(aa->ChemFixer.getPoint(aa))
					.collect(GeomUtil.convexHull())
					)
					;
			
		}catch(Exception ee){
			ee.printStackTrace();
			return Optional.empty();
		}

		return Optional.ofNullable(sw);

	}

	@Override
	public double score(Chemical c) {
		try{
			Chemical c2=c.copy();
			try {
				c2.generateCoordinates();
			} catch (MolwitchException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			return IntStream.range(0,c.getAtomCount()-1)
					.mapToObj(ii->Tuple.of(c.getAtom(ii), c2.getAtom(ii)))
					.map(Tuple.vmap(a->getNeighborShape(a)))
					.map(Tuple.kmap(a->getNeighborShape(a)))
					.filter(t->t.v().isPresent())
					.filter(t->t.k().isPresent())
					.map(t->Tuple.of(t.k().get(),t.v().get()))
					.mapToDouble(t->t.v().normalize().similarity(t.k().normalize()))
					.filter(d->!Double.isNaN(d))
					.average()
					.orElse(0);
		}catch(Exception e){
			e.printStackTrace();
		}
		return 0;


	}
	
//	public static void main(String[] as) throws Exception{
//		Chemical c= Chemical.parse("\n" + 
//				"  Molvec0106031917102D\n" + 
//				"\n" + 
//				" 43 46  0  0  0  0  0  0  0  0999 V2000\n" + 
//				"   -1.2796    0.2441    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -1.2796    1.2545    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    1.3428    1.7583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -2.9905    3.2683    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -2.1439   -0.2761    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    3.5875    4.0171    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    3.2870    3.0680    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    4.5759    4.2266    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    2.1982    0.2441    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    1.3428    0.7274    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -1.2796   -1.7484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -2.1439   -1.2444    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -2.1439   -3.2477    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -1.2796   -2.7574    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -3.8725    1.7583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -3.8725    2.7523    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -2.1439   -4.2546    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -2.1439    2.7523    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -2.9905    1.2545    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -2.1439    1.7583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -0.3958    1.7583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    0.4723    1.2545    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    3.0717    0.7274    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -2.9905   -1.7484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -4.7264    1.2545    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    2.9197    4.7609    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -4.7264    3.2683    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -2.9171   -2.7890    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    4.7260    2.7523    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    3.8817    2.2568    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    2.3327    2.7523    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    3.2926    1.4416    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    4.7260    1.7583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -2.9905    0.2441    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -1.2796   -4.7607    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -3.8725   -1.2444    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -0.3958    2.7523    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    1.3428    2.7523    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    2.3365    1.7583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    2.1982   -0.7751    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"   -3.8725   -3.2477    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    0.6224    3.4520    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				"    2.3420    0.8942    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
//				" 32 39  1  0\n" + 
//				" 15 19  1  0\n" + 
//				"  3 38  1  0\n" + 
//				" 13 17  1  0\n" + 
//				" 11 12  1  0\n" + 
//				"  3 39  1  0\n" + 
//				" 12 24  1  0\n" + 
//				" 13 14  1  0\n" + 
//				" 11 14  1  0\n" + 
//				" 15 16  2  0\n" + 
//				" 30 33  1  0\n" + 
//				" 30 32  1  0\n" + 
//				" 19 20  2  0\n" + 
//				" 38 42  2  0\n" + 
//				"  5 12  1  0\n" + 
//				"  5 34  2  0\n" + 
//				" 21 37  2  0\n" + 
//				"  6 26  2  0\n" + 
//				" 15 25  1  0\n" + 
//				"  9 40  2  0\n" + 
//				" 24 28  1  0\n" + 
//				"  2 21  1  0\n" + 
//				"  7 31  1  0\n" + 
//				"  2 20  1  0\n" + 
//				" 31 39  1  0\n" + 
//				"  2  1  1  6\n" + 
//				" 31 38  1  0\n" + 
//				"  1  5  1  0\n" + 
//				"  4 16  1  0\n" + 
//				" 13 28  1  0\n" + 
//				"  7 30  1  0\n" + 
//				"  4 18  2  0\n" + 
//				"  7  6  1  1\n" + 
//				"  3 10  1  1\n" + 
//				" 28 41  2  0\n" + 
//				" 18 20  1  0\n" + 
//				"  6  8  1  0\n" + 
//				" 39 43  1  1\n" + 
//				"  9 10  1  0\n" + 
//				" 16 27  1  0\n" + 
//				" 17 35  1  0\n" + 
//				" 21 22  1  0\n" + 
//				" 24 36  2  0\n" + 
//				" 29 30  1  0\n" + 
//				"  9 23  1  0\n" + 
//				"  3 22  1  6\n" + 
//				"M  END");
////		c.generateCoordinates();
//		ResultScorer rs=new LayoutLocalConsistencyScorer();
//		System.out.println(rs.score(c));
//	}

}
