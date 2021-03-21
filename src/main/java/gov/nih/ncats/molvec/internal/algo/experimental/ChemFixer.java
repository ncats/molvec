package gov.nih.ncats.molvec.internal.algo.experimental;

import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import gov.nih.ncats.common.stream.StreamUtil;
import gov.nih.ncats.molvec.internal.algo.Tuple;
import gov.nih.ncats.molvec.internal.util.GeomUtil;
import gov.nih.ncats.molvec.internal.util.GeomUtil.LineWrapper;
import gov.nih.ncats.molvec.internal.util.GeomUtil.ShapeWrapper;
import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.Bond;
import gov.nih.ncats.molwitch.Bond.BondType;
import gov.nih.ncats.molwitch.Bond.Stereo;
import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.MolwitchException;
import gov.nih.ncats.molwitch.isotopes.Elements;

public class ChemFixer {

	public static enum KnownMissingBond{
		LEFT,
		TOP,
		RIGHT,
		BOTTOM
	}
	
	public static enum FixType{
		NORMAL,
		MERGED,
		FAKE,
		NULL
	}
	public static class ChemFixResult{
		public FixType type=FixType.NORMAL;
		public Chemical c;
	}

	private static final int CLEAN_ADD_MISSING_BONDS = 35;
	private static final int CLEAN_FORCE_CARBON_ON_RARE_ATOM_RING = 36;
	private static final int CLEAN_FORCE_CARBON_ON_TRIVALENT_OXYGEN = 37;
	private static final int CLEAN_FORCE_CARBON_ON_BIVALENT_F = 38;
	private static final int CLEAN_FORCE_NITROGEN_ON_BIVALENT_CL = 39;
	private static final int CLEAN_REMOVE_O_CL_BOND = 40;
	private static final int CLEAN_FORCE_F_ON_TERM = 41;
	private static final int CLEAN_REMOVE_DASH_BOND_ON_PENT_C = 42;
	private static final int CLEAN_FORCE_VERTICAL_I_TO_C = 43;
	private static final int CLEAN_FORCE_IODINES_TO_BE_F = 44;
	private static final int CLEAN_REMOVE_H_WITH_DOUBLE_BONDS = 45;
	private static final int CLEAN_FORCE_PEPTIDE_BOND = 46;
	private static final int CLEAN_FORCE_DASH_N_TO_DOUBLE = 47;
	private static final int CLEAN_FORCE_BORON_TO_C = 48;
	private static final int CLEAN_FORCE_B_TO_BR = 49;
	private static final int CLEAN_FORCE_DASH_RING_AS_DOUBLE = 50;
	private static final int CLEAN_REMOVE_EXTRA_RING_LINK = 51;
	private static final int CLEAN_REMOVE_EXTRA_RING_LINK_MAKE_DOUBLE = 52;
	private static final int CLEAN_FORCE_REMOVE_SMALL_BOND = 53;
	private static final int CLEAN_FORCE_5_MEMBER_RING_DOUBLE = 54;
	private static final int CLEAN_FORCE_STEREO_ON_H = 55;
	private static final int CLEAN_FORCE_DOUBLE_BOND_FROM_RING_TO_N = 56;
	private static final int CLEAN_REMOVE_TINY_BOND_OFF_RING = 57;
	private static final int CLEAN_FORCE_DOUBLE_BOND_ON_CARDINAL_6RING = 58;  //heavy lifter
	private static final int CLEAN_FORCE_DOUBLE_BOND_ON_CARDINAL_5RING = 59; //???
	private static final int CLEAN_REMOVE_DASHED_BONDS_ON_TRIANGLES = 60;
	private static final int CLEAN_FORCE_N_ON_VERY_CLOSE_H_BOND = 61;
	private static final int CLEAN_REMOVE_DISTORTED_TRANGLES = 62;
	private static final int CLEAN_ADD_ATOM_ON_5_MEMBER_RING_WITH_LONG_BOND = 63;
	private static final int CLEAN_DEMOTE_DOUBLE_BOND_ON_PENT_C = 64;
	private static final int CLEAN_FORCE_QUAT_N_TO_C = 65;
	private static final int CLEAN_FORCE_N_ON_SHORT_N_DOUBLE_C = 66;
	private static final int CLEAN_FORCE_N_ON_SHORT_RING_C_SINGLE_C = 67;
	private static final int CLEAN_EXTEND_CLOSE_BONDS = 68;
	private static final int CLEAN_STITCH_CLOSE_COMPONENTS = 69; //heavy lifter
	private static final int CLEAN_SECOND_CLEAN = 70;
	private static final int CLEAN_S_TO_SI_WHEN_ORG = 71;
	private static final int CLEAN_S_TO_S2O_WHEN_ORG = 72;
	private static final BitSet DO_ALL = new BitSet();
	private static final int CLEAN_INFER_MISSING_FROM_MARGIN = 73;
	static{
		DO_ALL.set(0, 100);
	}
	
	
	
	public static Chemical makeMissingBonds(Chemical c){
		double avg= c.bonds().mapToDouble(b->b.getBondLength())
				.average().getAsDouble();
		GeomUtil.eachCombination(c.atoms()
				.filter(at->at.getSymbol().equals("N") || at.getSymbol().equals("C") || at.getSymbol().equals("O"))
				.filter(at->{
					if(at.getSymbol().equals("N")){
						return at.getBonds().stream().mapToInt(b->b.getBondType().getOrder()).sum() < 3;
					}else if(at.getSymbol().equals("C")){
						return at.getBonds().stream().mapToInt(b->b.getBondType().getOrder()).sum() < 4;
					}else if(at.getSymbol().equals("O")){

						if(at.getBonds().stream().mapToInt(b->b.getBondType().getOrder()).sum() < 2){
							long cc=at.getNeighbors().stream()
									.flatMap(aa->aa.getBonds().stream())
									.filter(bb->bb.getBondType().equals(BondType.DOUBLE))
									.filter(bb->bb.getAtom1().getSymbol().equals("O") || bb.getAtom2().getSymbol().equals("O") || 
											bb.getAtom1().getSymbol().equals("N") || bb.getAtom2().getSymbol().equals("N"))
									.count();
							if(cc==0){
								return true;
							}else{
								return false;
							}
						}
						return false;
					}else{
						return true;
					}
				})
				.collect(Collectors.toList()))
		.filter(t->{

			return !t.k().bondTo(t.v()).isPresent();
		})
		.filter(t->{
			double ds= t.k().getAtomCoordinates().distanceSquaredTo(t.v().getAtomCoordinates());

			if(ds<avg*avg*(1.3*1.3) && ds> avg*avg*0.75*0.75){
				return true;
			}
			return false;
		})
		//		.collect(Collectors.toList())
		//		.stream()
		.filter(t->{
			double dx= t.k().getAtomCoordinates().getX() - t.v().getAtomCoordinates().getX();
			double dy= t.k().getAtomCoordinates().getY() - t.v().getAtomCoordinates().getY();
			if(Math.abs(dx)< avg*.10 || Math.abs(dy)< avg*.10){
				return true;
			}else if(Math.abs(dx)< avg*.30 || Math.abs(dy)< avg*.30){
				//				System.out.println("Almost");

				Chemical cop= c.copy();

				Atom ca1=cop.getAtom(t.k().getAtomIndexInParent());
				Atom ca2=cop.getAtom(t.v().getAtomIndexInParent());

				Bond nb=cop.addBond(ca1, ca2, BondType.SINGLE);
								
				if(nb.isInRing() && ca1.getSmallestRingSize()>=5 || ca2.getSmallestRingSize()>=5){
					//					System.out.println("Is in ring");
					return true;
				}

			}
			return false;
		})
		//		.filter(t->false)
		.filter(t->{
			double sd = t.k().getAtomCoordinates().distanceSquaredTo(t.v().getAtomCoordinates());
			for(Atom na:t.k().getNeighbors()){
				if(na.getAtomCoordinates().distanceSquaredTo(t.v().getAtomCoordinates())<sd){
					return false;
				}
			}
			for(Atom na:t.v().getNeighbors()){
				if(na.getAtomCoordinates().distanceSquaredTo(t.k().getAtomCoordinates())<sd){
					return false;
				}
			}
			boolean shared=t.k().getNeighbors().stream().filter(nn->t.v().getNeighbors().contains(nn)).findAny().isPresent();

			return !shared;
		})
		.filter(t->{
			double dx= t.k().getAtomCoordinates().getX()-t.v().getAtomCoordinates().getX();
			double dy= t.k().getAtomCoordinates().getY()-t.v().getAtomCoordinates().getY();

			if(Math.abs(dx)<Math.abs(dy)){
				//vertical

				if(t.k().getNeighbors()
						.stream()
						.filter(nn->Math.abs(t.k().getAtomCoordinates().getX()-nn.getAtomCoordinates().getX())<0.3*avg)
						.count()>0){
					return false;
				}
				if(t.v().getNeighbors()
						.stream()
						.filter(nn->Math.abs(t.v().getAtomCoordinates().getX()-nn.getAtomCoordinates().getX())<0.3*avg)
						.count()>0){
					return false;
				}
			}else{
				//horizontal
				if(t.k().getNeighbors()
						.stream()
						.filter(nn->Math.abs(t.k().getAtomCoordinates().getY()-nn.getAtomCoordinates().getY())<0.3*avg)
						.count()>0){
					return false;
				}
				if(t.v().getNeighbors()
						.stream()
						.filter(nn->Math.abs(t.v().getAtomCoordinates().getY()-nn.getAtomCoordinates().getY())<0.3*avg)
						.count()>0){
					return false;
				}
			}
			return true;
		})
		.filter(t->{
			return !c.atoms()
					.filter(aa->!aa.bondTo(t.k()).isPresent())
					.filter(aa->!aa.equals(t.k()))
					.filter(aa->aa.getAtomCoordinates().distanceSquaredTo(t.k().getAtomCoordinates())<avg*avg*0.5*0.5)
					.findAny()
					.isPresent() &&
					!c.atoms()
					.filter(aa->!aa.bondTo(t.v()).isPresent())
					.filter(aa->!aa.equals(t.v()))
					.filter(aa->aa.getAtomCoordinates().distanceSquaredTo(t.v().getAtomCoordinates())<avg*avg*0.5*0.5)
					.findAny()
					.isPresent();
		})
		.collect(Collectors.toList())
		.forEach(t->{
			//			
			Chemical cop= c.copy();

			Atom ca1=cop.getAtom(t.k().getAtomIndexInParent());
			Atom ca2=cop.getAtom(t.v().getAtomIndexInParent());

			Bond nb=cop.addBond(ca1, ca2, BondType.QUADRUPLE);

			if(nb.isInRing() && ( (ca1.getSmallestRingSize()==3)||
					(ca1.getSmallestRingSize()==5 && ca2.getSmallestRingSize()==5))){
				if((ca1.getSmallestRingSize()==5 && ca2.getSmallestRingSize()==5)){
					cop.setAtomMapToPosition();
					cop.bonds().filter(bbb->!bbb.isInRing())
					.forEach(bbb->cop.removeBond(bbb));
					Chemical pent= StreamUtil.forIterable(cop.getConnectedComponents())
							.filter(ccc->ccc.atoms().filter(aaa->aaa.getAtomToAtomMap().getAsInt() == ca1.getAtomToAtomMap().getAsInt()).findAny().isPresent())
							.findFirst()
							.get();
					if(pent.getAtomCount()==5||pent.getAtomCount()==6||pent.getAtomCount()==9){
						c.addBond(t.k(),t.v(),BondType.SINGLE);
						t.k().setImplicitHCount(null);
						t.v().setImplicitHCount(null);
					}else{

					}


				}

			}else{
				c.addBond(t.k(),t.v(),BondType.SINGLE);
				t.k().setImplicitHCount(null);
				t.v().setImplicitHCount(null);
			}


		});
		return c;
	}
	
	public static double correlationToClean(Chemical c) throws MolwitchException{
		
		Chemical cop = c.copy();
		cop.setAtomMapToPosition();
		Chemical cop2= cop.copy();
		cop2.generateCoordinates();
		
		double avg1 = cop.bonds().mapToDouble(b->b.getBondLength()).average().getAsDouble();
		double avg2 = cop2.bonds().mapToDouble(b->b.getBondLength()).average().getAsDouble();
		
		double cx1= cop.atoms().mapToDouble(at->at.getAtomCoordinates().getX()).average().getAsDouble();
		double cy1= cop.atoms().mapToDouble(at->at.getAtomCoordinates().getY()).average().getAsDouble();
		double cx2= cop2.atoms().mapToDouble(at->at.getAtomCoordinates().getX()).average().getAsDouble();
		double cy2= cop2.atoms().mapToDouble(at->at.getAtomCoordinates().getY()).average().getAsDouble();
		
		
		
		List<Point2D> pts1 = cop.atoms().map(at->getPoint(at)).map(p->new Point2D.Double((p.getX()-cx1)/avg1, (p.getY()-cy1)/avg1))
				.collect(Collectors.toList());
		
		List<Point2D> pts2 = cop.atoms().map(at->getPoint(at)).map(p->new Point2D.Double((p.getX()-cx2)/avg1, (p.getY()-cy2)/avg2))
				.collect(Collectors.toList());
		
		Point2D cent = new Point2D.Double(0,0);
		double smallestTot = Double.POSITIVE_INFINITY; 
		for(int i=0;i<pts1.size();i++){
			Point2D anch1 = pts1.get(i);
			Point2D anch2 = pts2.get(i);
			AffineTransform afft1=GeomUtil.getTransformFromLineToLine(new Line2D.Double(cent, anch1), new Line2D.Double(cent, anch2), false);
			
			double totSq1=0;
			
			for(int j=0;j< pts1.size();j++){
				Point2D oa1=pts1.get(j);
				Point2D oa2=pts2.get(j);
				
				Point2D npt= afft1.transform(oa1, null);
				totSq1+=npt.distanceSq(oa2);
			}
			AffineTransform afft2=GeomUtil.getTransformFromLineToLine(new Line2D.Double(cent, anch1), new Line2D.Double(cent, anch2), true);
			
			double totSq2=0;
			
			for(int j=0;j< pts1.size();j++){
				Point2D oa1=pts1.get(j);
				Point2D oa2=pts2.get(j);
				
				Point2D npt= afft2.transform(oa1, null);
				totSq2+=npt.distanceSq(oa2);
			}
			smallestTot = Math.min(smallestTot,Math.min(totSq1, totSq2));
		}
		
		
		
		
		return Math.sqrt(smallestTot/pts1.size());
	}

	public static Chemical combineCloseBonds(Chemical c){
		//		if(true)return c;
		double avg= c.bonds().mapToDouble(b->b.getBondLength())
				.average().getAsDouble();
		GeomUtil.eachCombination(c.atoms()
				.filter(at->at.getSymbol().equals("N") || at.getSymbol().equals("C") || at.getSymbol().equals("O"))
				.filter(at->{
					if(at.getSymbol().equals("N")){
						return at.getBonds().stream().mapToInt(b->b.getBondType().getOrder()).sum() < 3;
					}else if(at.getSymbol().equals("C")){
						return at.getBonds().stream().mapToInt(b->b.getBondType().getOrder()).sum() < 4;
					}else if(at.getSymbol().equals("O")){
						return at.getBonds().stream().mapToInt(b->b.getBondType().getOrder()).sum() < 2;
					}else{
						return true;
					}
				})
				.collect(Collectors.toList()))
		.filter(t->{

			return !t.k().bondTo(t.v()).isPresent();
		})
		.filter(t->{
			double ds= t.k().getAtomCoordinates().distanceSquaredTo(t.v().getAtomCoordinates());

			if(ds<avg*avg*0.5*0.5){
				return true;
			}
			return false;
		})
		.collect(Collectors.toList())
		.forEach(t->{
			boolean doit=false;
			int b1=t.k().getBondCount();
			int b2=t.v().getBondCount();
			
			if(t.k().getSymbol().equals("C") && t.k().getBondCount()==1 && t.v().getBondCount()>=1){
				doit=true;
			}else{
				t=t.swap();
				if(t.k().getSymbol().equals("C") && t.k().getBondCount()==1 && t.v().getBondCount()>=1){
				doit=true;
				}
			}
			if(b1==b2 && t.k().getBonds().get(0).getBondLength()>t.v().getBonds().get(0).getBondLength()){
				t=t.swap();
			}
			
			
			if(doit){
				
				
				Atom oa=t.k().getNeighbors().get(0);
				LineWrapper lw=LineWrapper.of(new Line2D.Double(getPoint(oa),getPoint(t.v())));
				LineWrapper lw2=LineWrapper.of(new Line2D.Double(getPoint(oa),getPoint(t.k())));
				if(lw.absCosTheta(lw2)>Math.cos(35*Math.PI/180.0)){

					c.addBond(oa,t.v(),BondType.SINGLE);
					c.removeAtom(t.k());	
				}
				
				
				
			}
		});
		return c;
	}


	public static Tuple<int[],Double> closestAtoms(Chemical c1, Chemical c2){
		int mini=-1;
		int minj=-1;
		double minsq=99999;
		for(int i=0;i<c1.getAtomCount();i++){
			Atom a1=c1.getAtom(i);
			if(a1.getSymbol().equals("H"))continue;
			if(a1.getSymbol().equals("O") && sumOrder(a1)>=2)continue;
			for(int j=0;j<c2.getAtomCount();j++){
				Atom a2=c2.getAtom(j);
				if(a2.getSymbol().equals("H"))continue;
				if(a2.getSymbol().equals("O") && sumOrder(a2)>=2)continue;
				double[] xy1=a1.getAtomCoordinates().xy();
				double[] xy2=a2.getAtomCoordinates().xy();

				double l=0.8;
				double sq = Math.pow(Math.pow(Math.abs(xy1[0]-xy2[0]),l) + Math.pow(Math.abs(xy1[1]-xy2[1]),l),1/l) ;

				if(sq<minsq){
					minsq=sq;
					mini=i;
					minj=j;
				}
			}
		}
		return Tuple.of(new int[]{mini,minj},Math.sqrt(minsq));

	}
	public static Chemical cleanDupeBonds(Chemical c){
		Chemical c2=c.copy();
		Set<Bond> sb=new HashSet<Bond>();
		for(Bond b: c2.getBonds()){
			if(sb.contains(b)){
				c2.removeBond(b);
			}
			sb.add(b);
		}
		return c2;
	}
	
	
	
	
	public static Chemical dumbClean(Chemical c, BitSet ops){
		double avg= c.bonds().mapToDouble(b->b.getBondLength()).average().getAsDouble();
		
		c=cleanDupeBonds(c);
		
		if(ops.get(CLEAN_ADD_MISSING_BONDS)){
			makeMissingBonds(c);
		}

		if(ops.get(CLEAN_FORCE_CARBON_ON_RARE_ATOM_RING)){
			c.atoms()
			.filter(at->at.getSymbol().equals("H")  || at.getSymbol().equals("Cl") 
					||(at.getSymbol().equals("P") && sumOrder(at)<=4)
					||(at.getSymbol().equals("S") && sumOrder(at)>2)
					
					)
			.filter(at->at.getSmallestRingSize()==6)
			.forEach(at->{
				at.setAtomicNumber(6);
				at.setImplicitHCount(null);
			});
		}
		if(ops.get(CLEAN_FORCE_CARBON_ON_TRIVALENT_OXYGEN)){
			
			c.atoms()
			.filter(at->at.getSymbol().equals("O"))
			.filter(at->{
				return at.getBonds().stream().mapToInt(b->b.getBondType().getOrder()).sum()>=3;
			})
			.forEach(at->{
				at.setAtomicNumber(6);;
				at.setImplicitHCount(null);
			});
		}

		if(ops.get(CLEAN_FORCE_CARBON_ON_BIVALENT_F)){
			
			c.atoms()
			.filter(at->at.getSymbol().equals("F"))
			.filter(at->{
				return at.getBonds().stream().mapToInt(b->b.getBondType().getOrder()).sum()>=2;
			})
			.forEach(at->{
				at.setAtomicNumber(6);;
				at.setImplicitHCount(null);
			});
		}


		if(ops.get(CLEAN_FORCE_NITROGEN_ON_BIVALENT_CL)){
			c.atoms()
			.filter(at->at.getSymbol().equals("Cl"))
			.filter(at->{
				return at.getBonds().stream().mapToInt(b->b.getBondType().getOrder()).sum()>=2;
			})
			.forEach(at->{
				at.setAtomicNumber(7);;
				at.setImplicitHCount(null);
			});
		}

		if(ops.get(CLEAN_REMOVE_O_CL_BOND)){
			Chemical tc=c;
			c.atoms()
			.filter(at->at.getSymbol().equals("Cl"))
			.filter(at->at.getBondCount()==1)
			.filter(at->at.getNeighbors().get(0).getSymbol().equals("O"))
			.collect(Collectors.toList())
			.forEach(at->{
				tc.removeAtom(at);
			});
		}

		if(ops.get(CLEAN_FORCE_F_ON_TERM)){
			Chemical tc=c;
			c.atoms()
			.filter(at->at.getSymbol().equals("F"))
			.filter(at->at.getBondCount()>0)
			.forEach(at->{
				List<Atom> fs=at.getNeighbors().get(0).getNeighbors().stream().filter(nn->nn.getBondCount()==1)
						.filter(a2->a2.getSymbol().equals("C") || a2.getSymbol().equals("F") || a2.getSymbol().equals("I")|| a2.getSymbol().equals("B")|| 
								a2.getSymbol().equals("N") ||
								a2.getSymbol().equals("O") //riskier
								
								
								)
						.collect(Collectors.toList())
						;
				if(fs.size()==3){
					fs.forEach(att->{
						att.setAtomicNumber(9);	
						att.setImplicitHCount(null);
					});

				}
				if(fs.size()==2){
					if((at.getNeighbors().get(0).getImplicitHCount()==1 && at.getNeighbors().get(0).getNeighbors().stream().filter(nn->nn.getSymbol().equals("O")).count()==0)){
						Atom ca=at.getNeighbors().get(0);
						fs.forEach(att->{
							att.setAtomicNumber(9);	
							att.setImplicitHCount(null);
						});
						ShapeWrapper s=ShapeWrapper.of(ca.getNeighbors().stream().map(aa->getPoint(aa))
						.collect(GeomUtil.convexHull()));
						ShapeWrapper s2=s.normalize();
						ShapeWrapper s3=GeomUtil.makeNPolygonOriginCenter(3, 1);
						double sim=s2.similarity(s3);
						
						if(sim<0.5){
						//
							Atom an=tc.addAtom("F", ca.getAtomCoordinates().getX()+1,ca.getAtomCoordinates().getY());
							tc.addBond(ca,an, BondType.SINGLE);
						}
					}else{
						fs.forEach(att->{
							att.setAtomicNumber(9);	
							att.setImplicitHCount(null);
						});	
					}

				}


			});
		}

		if(!ops.get(CLEAN_REMOVE_DASH_BOND_ON_PENT_C)){
			Chemical tc=c;
			c.atoms()
			.filter(at->at.getSymbol().equals("C"))
			.filter(at->{
				return at.getBonds().stream().mapToInt(b->b.getBondType().getOrder()).sum()>=5;
			})
			.forEach(at->{
				at.getBonds().stream().filter(b->b.getBondType().equals(BondType.SINGLE))
				.filter(b->b.getStereo().equals(Stereo.DOWN)||b.getStereo().equals(Stereo.DOWN_INVERTED))
				.findFirst()
				.ifPresent(b->{
					tc.removeBond(b);
					b.getAtom1().setImplicitHCount(null);
					b.getAtom2().setImplicitHCount(null);
				});;

			});
		}

		if(ops.get(CLEAN_FORCE_VERTICAL_I_TO_C)){

			List<Atom> rem = c.atoms()
					.filter(at->at.getSymbol().equals("I"))
					.filter(at->at.getBondCount()>0)
					.filter(at->{
						Atom n1=at.getNeighbors().get(0);
						if(Math.abs(n1.getAtomCoordinates().getX()-at.getAtomCoordinates().getX())<0.3){
							return true;
						}else if(at.getBondCount()>1){
							return true;
						}
						return false;
					})
					.collect(Collectors.toList());
			rem.forEach(at->{
				at.setAtomicNumber(6);
				at.setImplicitHCount(null);
			});
		}
		
		if(ops.get(CLEAN_FORCE_IODINES_TO_BE_F)){
			long fcount=c.atoms().filter(at->at.getSymbol().equals("F")).count();
			if(fcount>0){
				c.atoms()
				.filter(at->at.getSymbol().equals("I"))
				.filter(at->at.getBondCount()==0 || at.getNeighbors().get(0)
					.getBonds()
					.stream()
					.filter(bb->bb.getBondType().getOrder()==2)
					.count()>0)
				.forEach(at->{
					at.setAtomicNumber(9);
				});
			}
		}

		if(ops.get(CLEAN_REMOVE_H_WITH_DOUBLE_BONDS)){

			Chemical tc=c;
			c.atoms()
			.filter(at->at.getSymbol().equals("H"))
			.filter(at->at.getBondCount()>0)
			.filter(at->{
				if(at.getBonds().get(0).getBondType().equals(BondType.DOUBLE)){
					return true;
				}
				return false;
			})
			.collect(Collectors.toList())
			.forEach(at->{
				tc.removeAtom(at);
			});
		}


		if(ops.get(CLEAN_FORCE_PEPTIDE_BOND)){


			c.atoms()
			.filter(at->at.getSymbol().equals("O"))
			.filter(at->at.getBondCount()==1)
			.filter(at->at.getBonds().get(0).getBondType().equals(BondType.SINGLE))
			//					.filter(at->Stereo.DOWN.equals(at.getBonds().get(0).getStereo()) ||at.getBonds().get(0).getStereo().equals(Stereo.DOWN_INVERTED))
			.filter(at->sumOrder(at.getNeighbors().get(0))
					<4)
			.filter(at->at.getNeighbors().get(0).getNeighbors().stream()
					.filter(nn->nn.getSymbol().equals("N"))
					.filter(nn->nn.getBonds().stream().allMatch(bn->bn.getBondType().equals(BondType.SINGLE))).count()>0)
			.forEach(at->{
				at.getBonds().get(0).setBondType(BondType.DOUBLE);
				at.setImplicitHCount(null);
				at.getNeighbors().get(0).setImplicitHCount(null);
			});
		}

		if(!ops.get(CLEAN_FORCE_DASH_N_TO_DOUBLE)){


			c.atoms()
			.filter(at->at.getSymbol().equals("N"))
			.filter(at->at.getBonds().stream().filter(bb->bb.getBondType().equals(BondType.DOUBLE)).count()==0)
			.filter(at->at.getBondCount()<=2)
			.forEach(na->{
				List<Bond> stereoBond = na.getBonds().stream().filter(bb->bb.getBondType().equals(BondType.SINGLE))
						.filter(bb->bb.getStereo().equals(Stereo.DOWN)||bb.getStereo().equals(Stereo.DOWN_INVERTED))
						.collect(Collectors.toList());
				for(Bond b:stereoBond){
					if(b.getOtherAtom(na).getBonds().stream().filter(bb->bb.getBondType().equals(BondType.DOUBLE)).count()>0){
						continue;
					}

					if(b.getAtom1().equals(na) && b.getStereo().equals(Stereo.DOWN)){
						b.setBondType(BondType.DOUBLE);
						b.getAtom1().setImplicitHCount(null);
						b.getAtom2().setImplicitHCount(null);
					}else if(b.getAtom2().equals(na) && b.getStereo().equals(Stereo.DOWN_INVERTED)){
						b.setBondType(BondType.DOUBLE);
						b.getAtom1().setImplicitHCount(null);
						b.getAtom2().setImplicitHCount(null);
					}
				}
			});
		}



		
		if(!ops.get(CLEAN_FORCE_BORON_TO_C)){
			Set<String> allowb=Arrays.stream("C,O,F,Cl".split(",")).collect(Collectors.toSet());
			c.atoms()
			.filter(at->at.getSymbol().equals("B"))
			.forEach(na->{

				if(na.getNeighbors().stream().anyMatch(nn->!allowb.contains(nn.getSymbol()))){
					na.setAtomicNumber(6);
					na.setImplicitHCount(null);
				}
			});
		}
		
		
		if(ops.get(CLEAN_FORCE_B_TO_BR)){

			//If there's a little MORE it's Br, if there's a little LESS, it's probably F
			//Br=35
			c.atoms()
			.filter(at->at.getSymbol().equals("B"))
			.filter(at->at.getBondCount()==1)
			.filter(at->at.getBonds().get(0).getBondType().getOrder()==1)
			.filter(at->at.getNeighbors().get(0).isInRing())
			.forEach(na->{
				na.setAtomicNumber(35);
			});
		}
		
		//!@@@@@@@@@@@@@@
		if(ops.get(CLEAN_FORCE_DASH_RING_AS_DOUBLE)){


			c.bonds()
			.filter(bb->bb.getBondType().getOrder()==1)
			.filter(bb->bb.getStereo().equals(Stereo.DOWN)||bb.getStereo().equals(Stereo.DOWN_INVERTED))
			.filter(bb->bb.getAtom1().getSymbol().equals("C") && bb.getAtom2().getSymbol().equals("C"))
			.filter(bb->bb.isInRing())
			.filter(bb->bb.getAtom1().getSmallestRingSize()==6 && bb.getAtom2().getSmallestRingSize()==6)
			.filter(bb->bb.getAtom1().getBonds().stream().mapToInt(b->b.getBondType().getOrder()).sum()<=3)
			.filter(bb->bb.getAtom2().getBonds().stream().mapToInt(b->b.getBondType().getOrder()).sum()<=3)
			.filter(bb->bb.getAtom1().getBonds().stream().mapToInt(b->b.getBondType().getOrder()).filter(k->k==2).count()==0)
			.filter(bb->bb.getAtom2().getBonds().stream().mapToInt(b->b.getBondType().getOrder()).filter(k->k==2).count()==0)
			.forEach(bb->{
				bb.setBondType(BondType.DOUBLE);

			});
		}

		c=cleanDupeBonds(c);
		
		if(ops.get(CLEAN_REMOVE_EXTRA_RING_LINK_MAKE_DOUBLE)){

			Chemical tc=c;

			
			c.bonds()
			.filter(bb->bb.getBondLength()<avg*0.4)
			.forEach(bb->{
				Atom tlatom=null;
				Atom tpatom=null;
				if(bb.getAtom1().getBondCount()==2 && bb.getAtom2().getBondCount()==3){
					tlatom=bb.getAtom1();
					tpatom=bb.getAtom2();
				}else if(bb.getAtom1().getBondCount()==3 && bb.getAtom2().getBondCount()==2){
					tlatom=bb.getAtom2();
					tpatom=bb.getAtom1();
				}
				Atom latom=tlatom;
				Atom patom=tpatom;

				if(latom!=null){
					List<Bond> ringBonds=patom.getBonds().stream().filter(b->!bb.equals(b))
							.filter(b->b.isInRing())
							.collect(Collectors.toList());

					if(ringBonds.size()==2){

						ringBonds.stream()
						.filter(rb->rb.getOtherAtom(patom).getBonds().stream()
								.filter(bn->bn.getBondType().equals(BondType.DOUBLE))
								.count()<=0)
						.findFirst()
						.ifPresent(bbb->{
							bbb.setBondType(BondType.DOUBLE);
							bbb.getAtom1().setImplicitHCount(null);
							bbb.getAtom2().setImplicitHCount(null);
						});
						Atom nlatom=latom.getBonds().stream().filter(lb->!lb.equals(bb))
								.findFirst().get().getOtherAtom(latom);
						tc.removeAtom(latom);
						tc.addBond(patom, nlatom, BondType.SINGLE);

					}

				}

			});
		}
		
		
		if(ops.get(CLEAN_FORCE_REMOVE_SMALL_BOND)){

			Chemical tc=c;
			c.bonds()
			.filter(bb->bb.getBondLength()<avg*0.3)
			.filter(bb->bb.getAtom1().getSymbol().equals("C") && bb.getAtom2().getSymbol().equals("C") )
			.forEach(bb->{
				Atom toD=null;
				if(bb.getAtom1().getBondCount()==2 && bb.getAtom2().getBondCount()!=2){
					toD=bb.getAtom1();
				}else if(bb.getAtom2().getBondCount()==2 && bb.getAtom1().getBondCount()!=2){
					toD=bb.getAtom2();
				}
				if(toD!=null){
					List<Atom> nat=toD.getNeighbors();
					tc.removeAtom(toD);
					tc.addBond(nat.get(0),nat.get(1), BondType.SINGLE);
				}

			});
		}

		c=cleanDupeBonds(c);
		
		if(ops.get(CLEAN_FORCE_5_MEMBER_RING_DOUBLE)){
			try{
						//			double avg= c.bonds().mapToDouble(b->b.getBondLength()).average().getAsDouble();
						c.bonds()
						.filter(bb->bb.getBondType().equals(BondType.SINGLE))
						.filter(bb->bb.getStereo().equals(Stereo.DOWN) || bb.getStereo().equals(Stereo.DOWN_INVERTED)  )
						.filter(bb->bb.isInRing())
						.filter(bb->bb.getAtom1().getSmallestRingSize()==5)
						.filter(bb->bb.getAtom2().getSmallestRingSize()==5)
						.filter(bb->!bb.getAtom2().getSymbol().equals("O") && !bb.getAtom2().getSymbol().equals("S"))
						.filter(bb->!bb.getAtom1().getSymbol().equals("O") && !bb.getAtom1().getSymbol().equals("S"))
						.filter(bb->Stream.of(bb.getAtom1(),bb.getAtom2()).flatMap(at->at.getBonds().stream()).filter(b1->b1.getBondType().equals(BondType.DOUBLE)).count()==0)
						.filter(b->Stream.of(b.getAtom1(),b.getAtom2())
								.filter(at->at.getSymbol().equals("N"))
								.filter(at->sumOrder(at)>=3)
								.count() ==0
								)
						.collect(Collectors.toList())
						
						.forEach(bb->{
							
							bb.setBondType(BondType.DOUBLE);
			
						});
			}catch(Exception e){
				throw e;
			}
		}

		c.atoms().forEach(at->{
			if(at.getCharge()!=0){
				at.setCharge(0);
				at.setImplicitHCount(null);
			}

		});

		if(ops.get(CLEAN_FORCE_STEREO_ON_H)){

			Chemical tc=c;
			c.atoms().filter(at->at.getSymbol().equals("H"))
	
			.forEach(at->{
				if(at.getBondCount()>0){
					Bond b = at.getBonds().get(0);
					if(b.getBondType().equals(BondType.DOUBLE)){
						b.setBondType(BondType.SINGLE);
					}
					if(b.getStereo().equals(Stereo.NONE)){
						if(b.getAtom1().equals(at)){
							b.setStereo(Stereo.DOWN_INVERTED);
						}else{
							b.setStereo(Stereo.DOWN);
						}
					}}else{
						tc.removeAtom(at);
					}
			});
		}



		c=cleanDupeBonds(c);
		
		
		
		if(ops.get(CLEAN_FORCE_DOUBLE_BOND_FROM_RING_TO_N)){

			Chemical tc=c;
			//remove dumb double bonds
			c.atoms()
			.filter(at->at.getSmallestRingSize()==6)
			.filter(at->at.getSymbol().equals("C"))

			.filter(b->b.getBonds().stream()
					.filter(bn->bn.getBondType().equals(BondType.DOUBLE))
					.count()==2)
			.flatMap(a->a.getBonds().stream().filter(b->b.getBondType().getOrder()==2).map(b->b.getOtherAtom(a)))
			.filter(a->a.getBondCount()==1)
			.collect(Collectors.toList())
			.forEach(a->{
				a.getNeighbors().get(0).setAtomicNumber(7);
				tc.removeAtom(a);
			});
		}

		if(ops.get(CLEAN_ADD_MISSING_BONDS)){

			makeMissingBonds(c);
		}
		
		
		
		c=cleanDupeBonds(c);
		
		if(ops.get(CLEAN_REMOVE_TINY_BOND_OFF_RING)){
			Chemical tc=c;
			//			
			c.bonds()
			.filter(b->b.getBondLength()<avg*0.35)
			//			.peek(b->System.out.println(b))
			.filter(b->b.getAtom1().getSymbol().equals("C") && b.getAtom2().getSymbol().equals("C") )
			.filter(b->(b.getAtom1().getSmallestRingSize()==6 && b.getAtom2().getBondCount()==1) ||
					(b.getAtom2().getSmallestRingSize()==6 && b.getAtom1().getBondCount()==1)
					)
			.forEach(b->{

				Atom at1=(b.getAtom1().getBondCount()==1)?b.getAtom1():b.getAtom2();
				Atom at2=b.getOtherAtom(at1);
				at2.setAtomCoordinates(at1.getAtomCoordinates());
				tc.removeAtom(at1);
				//				b.setBondType(BondType.DOUBLE);
			});
			;
		}

		if(ops.get(CLEAN_EXTEND_CLOSE_BONDS)){
			combineCloseBonds(c);
		}
		
		if(ops.get(CLEAN_FORCE_DOUBLE_BOND_ON_CARDINAL_6RING)){

			boolean[] changed = new boolean[]{false};
			
			Set<String> createDoublesFor = Stream.of("C","N").collect(Collectors.toSet());
			do{
			changed[0]=false;
			c.atoms()
			.filter(at->at.getSmallestRingSize()==6)
			.flatMap(at->at.getBonds().stream())
			.filter(b->b.getBondType().equals(BondType.SINGLE))
			.filter(b->Math.abs(b.getAtom1().getAtomCoordinates().getX() - b.getAtom2().getAtomCoordinates().getX())<avg*0.2 ||
					Math.abs(b.getAtom1().getAtomCoordinates().getY() - b.getAtom2().getAtomCoordinates().getY())<avg*0.2
					)
			.filter(b->b.getAtom1().getSmallestRingSize()>=5 &&b.getAtom2().getSmallestRingSize()>=5)
			.filter(b->b.getAtom1().getBonds().stream().filter(bb->!bb.getBondType().equals(BondType.SINGLE)).count()==0)
			.filter(b->b.getAtom2().getBonds().stream().filter(bb->!bb.getBondType().equals(BondType.SINGLE)).count()==0)
			.filter(b->createDoublesFor.contains(b.getAtom1().getSymbol()))
			.filter(b->createDoublesFor.contains(b.getAtom2().getSymbol()))
			//			.peek(b->System.out.println(b))
			.filter(b->Stream.of(b.getAtom1(),b.getAtom2())
					.filter(at->at.getSymbol().equals("N"))
					.filter(at->sumOrder(at)>=3)
					.count() ==0
					)
			.filter(b->Stream.of(b.getAtom1(),b.getAtom2())
					.filter(at->at.getSymbol().equals("C"))
					.filter(at->sumOrder(at)>=4)
					.count() ==0
					)
			.filter(b->b.getAtom1().getNeighbors().stream().filter(at->at.bondTo(b.getAtom1()).get().isInRing()).filter(at->at.getSmallestRingSize()>=5).flatMap(aa->aa.getBonds().stream()).filter(bb->bb.getBondType().equals(BondType.DOUBLE)).count()>=1)
			.filter(b->b.getAtom2().getNeighbors().stream().filter(at->at.bondTo(b.getAtom2()).get().isInRing()).filter(at->at.getSmallestRingSize()>=5).flatMap(aa->aa.getBonds().stream()).filter(bb->bb.getBondType().equals(BondType.DOUBLE)).count()>=1)

			.forEach(b->{
				
				changed[0]=true;
				b.setBondType(BondType.DOUBLE);
			});
			}while(changed[0]);
		}
		
		if(!ops.get(CLEAN_FORCE_DOUBLE_BOND_ON_CARDINAL_5RING)){

			Set<String> createDoublesFor = Stream.of("C","N").collect(Collectors.toSet());
			c.atoms()
			.filter(at->at.getSmallestRingSize()==5)
			.flatMap(at->at.getBonds().stream())
			.filter(b->b.getBondType().equals(BondType.SINGLE))
			.filter(b->Math.abs(b.getAtom1().getAtomCoordinates().getX() - b.getAtom2().getAtomCoordinates().getX())<avg*0.2 ||
					Math.abs(b.getAtom1().getAtomCoordinates().getY() - b.getAtom2().getAtomCoordinates().getY())<avg*0.2
					)
			.filter(b->b.getAtom1().getSmallestRingSize()==5 &&b.getAtom2().getSmallestRingSize()==5)
			.filter(b->b.getAtom1().getBonds().stream().filter(bb->!bb.getBondType().equals(BondType.SINGLE)).count()==0)
			.filter(b->b.getAtom2().getBonds().stream().filter(bb->!bb.getBondType().equals(BondType.SINGLE)).count()==0)
			.filter(b->createDoublesFor.contains(b.getAtom1().getSymbol()))
			.filter(b->createDoublesFor.contains(b.getAtom2().getSymbol()))
			.filter(b->Stream.of(b.getAtom1(),b.getAtom2())
					.filter(at->at.getSymbol().equals("N"))
					.filter(at->sumOrder(at)>=3)
					.count() ==0
					)
			.filter(b->Stream.of(b.getAtom1(),b.getAtom2())
					.filter(at->at.getSymbol().equals("C"))
					.filter(at->sumOrder(at)>=4)
					.count() ==0
					)
			//			.peek(b->System.out.println(b))
			.filter(b->{
				boolean at1dub=b.getAtom1().getNeighbors().stream().filter(at->at.bondTo(b.getAtom1()).get().isInRing()).filter(at->at.getSmallestRingSize()>=5).flatMap(aa->aa.getBonds().stream()).filter(b1->b1.isInRing()).filter(bb->bb.getBondType().equals(BondType.DOUBLE)).count()>=1;
				boolean at2dub=b.getAtom2().getNeighbors().stream().filter(at->at.bondTo(b.getAtom2()).get().isInRing()).filter(at->at.getSmallestRingSize()>=5).flatMap(aa->aa.getBonds().stream()).filter(b1->b1.isInRing()).filter(bb->bb.getBondType().equals(BondType.DOUBLE)).count()>=1;
				return at1dub^at2dub;
			})
			.forEach(b->{
				b.setBondType(BondType.DOUBLE);
			});
			;
		}
		
		if(ops.get(CLEAN_REMOVE_DASHED_BONDS_ON_TRIANGLES)){
			Chemical tc=c;
			//remove dumb triangles
			c.atoms()
			.filter(at->at.getSmallestRingSize()==3)
			.flatMap(at->at.getBonds().stream())
			.filter(b->b.getBondType().equals(BondType.SINGLE))
			.filter(b->b.getStereo().equals(Stereo.DOWN)||b.getStereo().equals(Stereo.DOWN_INVERTED))
			.distinct()
			.filter(b->b.getAtom1().getSmallestRingSize()==3 &&b.getAtom2().getSmallestRingSize()==3)
			.collect(Collectors.toList())
			.forEach(b->{
				boolean addback=true;
				if(b.getAtom1().getNeighbors().stream().filter(nn->nn.getSmallestRingSize()==4).count()>0 ||
						b.getAtom2().getNeighbors().stream().filter(nn->nn.getSmallestRingSize()==4).count()>0
						){
					addback=false;
				}
				tc.removeBond(b);
				
				if(b.getAtom1().getSmallestRingSize()==6 ||b.getAtom2().getSmallestRingSize()==6 ){
					addback=false;
				}
				
				if(addback){
					tc.addBond(b);
				}
			});
			c=c.copy();
		}
		
		if(ops.get(CLEAN_FORCE_N_ON_VERY_CLOSE_H_BOND)){
			Chemical tc=c;
			//probably N
			c.atoms()
			.filter(at->at.getSymbol().equals("H"))
			.filter(at->at.getNeighbors().stream().filter(n->n.getSymbol().equals("C")).count()>0)
			.filter(at->at.getBonds().stream().filter(b->b.getBondLength()<avg*0.4).count()>0)
			.collect(Collectors.toList())
			.forEach(b->{
				Atom a1=b.getNeighbors().get(0);
				if(sumOrder(a1)==3){
						a1.setAtomicNumber(7);
						a1.setImplicitHCount(null);
				}else if(sumOrder(a1)==2){
					a1.getNeighbors().stream()
					.filter(nn->!b.equals(nn))
					.findFirst()
					.ifPresent(oo->{
						tc.removeAtom(a1);
						Bond bn=tc.addBond(oo,b, BondType.SINGLE);
						bn.setStereo(Stereo.DOWN);
					});
					
					
//					c.addBond(b)
				}
			});
		}

		if(ops.get(CLEAN_REMOVE_DISTORTED_TRANGLES)){
			Chemical tc=c;
			c.atoms()
			.filter(at->at.getSymbol().equals("C"))
			.filter(at->at.getBondCount()==2)
			.forEach(at->{
				Atom a1=at.getNeighbors().get(0);
				Atom a2=at.getNeighbors().get(1);
				if(a1.bondTo(a2).isPresent()){
					double hyp = a1.bondTo(a2).get().getBondLength();
					//triangle
					if(
							//at.getBonds().get(0).getBondLength()+ at.getBonds().get(1).getBondLength() < avg*1.7 || 
							at.getBonds().get(0).getBondLength()+ at.getBonds().get(1).getBondLength() < hyp*1.7 ||
							at.getBonds().get(0).getBondLength() < hyp*0.6 || at.getBonds().get(1).getBondLength() < hyp*0.6){
						Bond bb=a1.getBonds().stream()
								.filter(b->b.getOtherAtom(a1).equals(a2))
								.findFirst().get();
						if(bb.getBondType().getOrder() == 1 && !bb.getStereo().equals(Stereo.NONE)){
							bb.setBondType(BondType.DOUBLE);
							bb.getAtom1().setImplicitHCount(null);
							bb.getAtom2().setImplicitHCount(null);
						}
						tc.removeAtom(at);
						a1.setImplicitHCount(null);
						a2.setImplicitHCount(null);
					}
				}
			});
		}

		if(ops.get(CLEAN_ADD_ATOM_ON_5_MEMBER_RING_WITH_LONG_BOND)){
			Chemical tc=c;
			Set<Bond> obond = getOverlappingBonds(c).stream().flatMap(t->Stream.of(t.k(),t.v())).collect(Collectors.toSet());
			c.bonds()
			.filter(b->b.getBondLength()>avg*1.3)
			.filter(b->b.isInRing())
			.filter(b->b.getAtom1().getSmallestRingSize()==5 && b.getAtom2().getSmallestRingSize()==5 )
			.filter(b->b.getAtom1().getSymbol().equals("C") && b.getAtom2().getSymbol().equals("C") )
			.filter(b->b.getAtom1().getSymbol().equals("C") && b.getAtom2().getSymbol().equals("C") )
			.filter(b->!obond.contains(b))
			.collect(Collectors.toList())
			.forEach(b->{
				Atom a1=b.getAtom1();
				Atom a2=b.getAtom2();
				double fy=(a1.getAtomCoordinates().getY() +a2.getAtomCoordinates().getY())/2;
				double fx=(a1.getAtomCoordinates().getX()+ a2.getAtomCoordinates().getX())/2;
				Atom mat=tc.addAtom("C", fx,fy);
				
				tc.removeBond(b);
//				Bond b1=tc.addBond(b.getAtom1(), mat, BondType.SINGLE);
//				Bond b2=tc.addBond(b.getAtom2(), mat, BondType.SINGLE);
//				if(b.getAtom1().getBonds().stream().filter(bb->bb.getBondType().getOrder()==2).count()==1){
//					b2.setBondType(BondType.DOUBLE);
//				}else if(b.getAtom2().getBonds().stream().filter(bb->bb.getBondType().getOrder()==2).count()==1){
//					b1.setBondType(BondType.DOUBLE);
//				}
				mat.setImplicitHCount(null);
				b.getAtom1().setImplicitHCount(null);
				b.getAtom2().setImplicitHCount(null);
			});
		}


		if(!ops.get(CLEAN_DEMOTE_DOUBLE_BOND_ON_PENT_C)){
			
			c.atoms()
			.filter(at->at.getSymbol().equals("C"))
			.filter(at->sumOrder(at)>4)
			.filter(at->at.getBondCount(BondType.DOUBLE)>=2)
			.forEach(at->{
				List<Bond> tbds=at.getBonds().stream().filter(bb->bb.getBondType().getOrder()==2)
				.filter(b->b.getOtherAtom(at).getSymbol().equals("O"))
				.collect(Collectors.toList());
				;
				if(tbds.size()==2){
					tbds.get(0).setBondType(BondType.SINGLE);
				}
			});
		}
		
		
		if(ops.get(CLEAN_FORCE_QUAT_N_TO_C)){

			c.atoms()
			.filter(at->at.getSymbol().equals("N"))
			.forEach(na->{
				if(sumOrder(na)>3){
					na.setAtomicNumber(6);
				}
			});
		}
		

		if(ops.get(CLEAN_FORCE_N_ON_SHORT_N_DOUBLE_C)){
			Chemical ct=c;

			c.atoms()
			.filter(at->at.getSymbol().equals("N"))
			.filter(at->at.getBondCount()==1)
			.filter(at->at.getBonds().get(0).getBondType().getOrder()==2)
			.filter(at->at.getBonds().get(0).getOtherAtom(at).getSymbol().equals("C"))
			.map(a->Tuple.of(a,a.getBonds().get(0)))
			.filter(b->b.v().getBondLength()<avg*0.5)
			.collect(Collectors.toList())
			.forEach(nb->{
				Atom aa=nb.v().getOtherAtom(nb.k());
				ct.removeAtom(nb.k());
				aa.setAtomicNumber(7);
			});
		}
		
		if(ops.get(CLEAN_FORCE_N_ON_SHORT_RING_C_SINGLE_C)){
			Chemical ct=c;

			c.atoms()
			.filter(at->at.getSymbol().equals("C"))
			.filter(at->at.getBondCount()==1)
			.filter(at->at.getBonds().get(0).getBondType().getOrder()==1)
			.filter(at->at.getBonds().get(0).getOtherAtom(at).getSymbol().equals("C"))
			.filter(at->at.getBonds().get(0).getOtherAtom(at).isInRing())
			.filter(at->sumOrder(at.getBonds().get(0).getOtherAtom(at))==3)
			.map(a->Tuple.of(a,a.getBonds().get(0)))
			.filter(b->b.v().getBondLength()<avg*0.5)
			//todo:roughly vertical
			
			.collect(Collectors.toList())
			.forEach(nb->{
				Atom aa=nb.v().getOtherAtom(nb.k());
				ct.removeAtom(nb.k());
				aa.setAtomicNumber(7);
			});
		}
		
		
		if(ops.get(CLEAN_S_TO_SI_WHEN_ORG)){
			c.atoms()
			.filter(at->at.getSymbol().equals("S"))
			.filter(at->sumOrder(at)>2)
			.filter(at->at.getBonds().stream().filter(bb->!bb.getBondType().equals(BondType.SINGLE)).count()==0)
			.forEach(nb->{
				nb.setAtomicNumber(14);
			});
		}
		
		
		c.atoms().forEach(at->{
			at.setImplicitHCount(null);		
		});
		
		
		
		

		cleanDupeBonds(c);

		return c;
	}
	
	public static int sumOrder(Atom a){
		return (int) a.getBonds().stream().mapToDouble(b->b.getBondType().getOrder()).sum();
	}

	public static Rectangle2D bounds(Chemical c){
		return c.atoms()
				.map(a->a.getAtomCoordinates())
				.map(a->new Point2D.Double(a.getX(), a.getY()))
				.collect(GeomUtil.convexHull())
				.getBounds2D();
	}

	public static boolean isBound(Bond b1, Bond b2){
		return Stream.of(b1.getAtom1(),b1.getAtom2(), 
				b2.getAtom1(),b2.getAtom2())
				.distinct()
				.count()!=4;

	}
	public static Point2D getPoint(Atom a){
		return new Point2D.Double(a.getAtomCoordinates().getX(),a.getAtomCoordinates().getY());
	}
	public static Tuple<Bond,Line2D> wrapB(Bond b){
		return Tuple.of(b, new Line2D.Double(getPoint(b.getAtom1()),getPoint(b.getAtom2())));
	}

	public static List<Tuple<Bond,Bond>> getOverlappingBonds(Chemical c){
		return GeomUtil.eachCombination(c.bonds().collect(Collectors.toList()))
				.filter(t->!isBound(t.k(),t.v()))
				.map(Tuple.vmap(b->wrapB(b)))
				.map(Tuple.kmap(b->wrapB(b)))
				.filter(t->GeomUtil.segmentIntersection(t.k().v(), t.v().v()).isPresent())
				.map(Tuple.vmap(b->b.k()))
				.map(Tuple.kmap(b->b.k()))
				.collect(Collectors.toList());
	}

	public static void setHStereo(Chemical c){
		c.atoms().filter(at->at.getSymbol().equals("H"))
		.forEach(at->{
					if (at.getBondCount() > 0 && !at.getNeighbors().get(0).getSymbol().equals("H")) {
						Bond b = at.getBonds().get(0);
						b.setBondType(BondType.SINGLE);
						
						if (at.getNeighbors().get(0).getBonds().stream()
								.filter(bb -> bb.getBondType().equals(BondType.DOUBLE)).count() == 0) {
							if (b.getStereo().equals(Stereo.NONE)) {
								if (b.getAtom1().equals(at)) {
									b.setStereo(Stereo.DOWN_INVERTED);
								} else {
									b.setStereo(Stereo.DOWN);
								}
							}
						}else{
							b.setStereo(Stereo.NONE);
						}
					} else {
						c.removeAtom(at);
					}
		});
	}

	public static Tuple<Chemical,Boolean> stitchChemical(Chemical cf){
		
		boolean stitched=false;
		for(int li=0;li<20;li++){
			double avg= cf.bonds().mapToDouble(b->b.getBondLength())
					.average().getAsDouble();
			Chemical c=cf.copy();
			c.clearAtomMaps();
			c.setAtomMapToPosition();

			//			System.out.println("loop");
			List<Chemical> comps = StreamUtil.forIterable(c.getConnectedComponents())
					.collect(Collectors.toList());
			if(comps.size()>1){
				boolean tryit=true;
				if(li==0){
					List<Tuple<Bond,Bond>> intesecting = getOverlappingBonds(c);
					if(intesecting.size()>0 && intesecting.size()<3){

						Chemical cc=c;
						intesecting.forEach(t->{
							Point2D pi=GeomUtil.segmentIntersection(wrapB(t.k()).v(),wrapB(t.v()).v()).get();
							Atom na=cc.addAtom("C", pi.getX(), pi.getY());
							cc.removeBond(t.k());
							cc.removeBond(t.v());
							cc.addBond(t.k().getAtom1(),na, BondType.SINGLE);
							cc.addBond(t.k().getAtom2(),na, BondType.SINGLE);
							cc.addBond(t.v().getAtom1(),na, BondType.SINGLE);
							cc.addBond(t.v().getAtom2(),na, BondType.SINGLE);
						});
						comps = StreamUtil.forIterable(c.getConnectedComponents())
								.collect(Collectors.toList());
						if(comps.size()==1){
							tryit=false;
						}
						cf=c;
						continue;
					}
				}
				if(tryit){
					int mini=-1;
					int minj=-1;
					double mind=9999;
					for(int i=0;i<comps.size();i++){
						Chemical c1 = comps.get(i);

						for(int j=i+1;j<comps.size();j++){
							Chemical c2 = comps.get(j);
							Tuple<int[],Double> pair=closestAtoms(c1,c2);

							if(pair.v()<mind){
								try{
									int mi=c1.getAtom(pair.k()[0]).getAtomToAtomMap().getAsInt()-1;
									int mj=c2.getAtom(pair.k()[1]).getAtomToAtomMap().getAsInt()-1;
									mind=pair.v();
									mini=mi;
									minj=mj;
								}catch(Exception e){
									
								}
							}
						}
					}
					//							System.out.println("mini:" + mini);
					//							System.out.println("minj:" + minj);
					if(mini>=0 && minj>=0){
						Atom a1=c.getAtom(mini);
						Atom a2=c.getAtom(minj);
						boolean addBond=true;
						if(a1.getAtomCoordinates().distanceSquaredTo(a2.getAtomCoordinates()) < avg*avg*0.5*0.5){
							//probably incomplete extension
							if(!a1.getSymbol().equals("C") && a2.getSymbol().equals("C")){
								Atom t=a1;
								a1=a2;
								a2=t;
							}
							if(a1.getSymbol().equals("C") && !a2.getSymbol().equals("C")){
								addBond=false;
							}
						}
						if(a1.getBondCount()==0){
							addBond=true;
						}

						if(addBond){
							c.addBond(a1,a2,BondType.SINGLE);
							a1.setImplicitHCount(null);
							a2.setImplicitHCount(null);
						}else{
							
							Bond ob=a1.getBonds().get(0);
							Atom na =ob.getOtherAtom(a1);
							c.removeAtom(a1);

							c.addBond(a2,na,ob.getBondType());
							na.setImplicitHCount(null);
							a2.setImplicitHCount(null);
						}
						stitched=true;
					}
					//					res.type=FixType.MERGED;
				}
				cf=c;
				
			}else{


				break;
			}
		}
		return Tuple.of(cf,stitched);

	}
	
	public static ChemFixResult fixChemical(Chemical cf){
		return fixChemical(cf, new HashSet<>());
	}
	public static Chemical secondClean(Chemical ct, BitSet ops){
		Chemical c=ct.copy();
		c.atoms().filter(ca->ca.getSymbol().equals("F"))
		.flatMap(ca->ca.getNeighbors().stream())
		.distinct()
		.filter(ca->sumOrder(ca)>4)
		.forEach(ca->{
			ca.getNeighbors().stream().filter(nn->nn.getSymbol().equals("F"))
			.map(tt->Tuple.of(tt,-tt.getAtomCoordinates().distanceSquaredTo(ca.getAtomCoordinates())).withVComparator())
			.sorted()
			.map(tt->tt.k())
			.findFirst().ifPresent(nn->{
				c.removeAtom(nn);
			});
		});
		
		c.bonds().filter(b->b.getBondType().equals(BondType.AROMATIC) || b.getBondType().equals(BondType.SINGLE_OR_DOUBLE) || b.getBondType().equals(BondType.QUADRUPLE))
//		.peek(b->System.out.println("Found"))
//		.peek(b->b.setBondType(BondType.SINGLE))
//		.filter(b->b.getBondType().equals(BondType.AROMATIC) || b.getBondType().equals(BondType.SINGLE_OR_DOUBLE) || b.getBondType().equals(BondType.QUADRUPLE))
//		.peek(b->System.out.println("Found AGAIN"))
		.forEach(b->{
//			c.removeBond(b);
			b.setBondType(BondType.SINGLE);
		});
		;
//		try {
//			System.out.println(c.toMol());
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
		
		if(ops.get(CLEAN_S_TO_S2O_WHEN_ORG)){
			c.atoms()
			.filter(at->at.getSymbol().equals("S"))
			.filter(at->sumOrder(at)>=4)
			.filter(at->sumOrder(at)<6)
			.filter(at->at.getNeighbors().stream().filter(nn->nn.getSymbol().equals("O")).count()>=2)
			.forEach(nb->{
				nb.getNeighbors().stream().filter(nn->nn.getSymbol().equals("O"))
				.filter(b->b.bondTo(nb).get().getBondType().getOrder()==1)
				.limit(2)
				.forEach(at->at.bondTo(nb).get().setBondType(BondType.DOUBLE));;
			});
		}
		return c;
		
	}

	public static ChemFixResult fixChemical(Chemical cf, BitSet bs){
		return fixChemical(cf,new HashSet<>(),bs);
	}
	
	public static ChemFixResult fixChemical(Chemical cf, Set<KnownMissingBond> mset ){
		return fixChemical(cf,mset,DO_ALL);
	}

	public static ChemFixResult fixChemical(Chemical cf, Set<KnownMissingBond> mset , BitSet bs){
		ChemFixResult res = new ChemFixResult();

		res.c=cf;
		try {


			cf=dumbClean(cf, bs);	
			if(bs.get(CLEAN_STITCH_CLOSE_COMPONENTS)){
				Tuple<Chemical, Boolean> tup=stitchChemical(cf);
				cf=tup.k();
				if(tup.v()){
					res.type=FixType.MERGED;
				}
			}
			if(bs.get(CLEAN_SECOND_CLEAN)){
				cf=secondClean(cf,bs);
			}
			setHStereo(cf);

			String[] lines=cf.toMol().split("\n");
			for(int ii=4;ii<4+Integer.parseInt(lines[3].substring(0,3).trim());ii++){
				lines[ii]=lines[ii].substring(0,33) + "  0  0  0  0  0  0  0  0  0  0  0  0";
			}
			String mm = Arrays.stream(lines).collect(Collectors.joining("\n"));
			cf= Chemical.parse(mm);
			
			if(!bs.get(CLEAN_INFER_MISSING_FROM_MARGIN)){
			if(mset.size()>0){
				double avg = cf.bonds().mapToDouble(b->b.getBondLength()).average().getAsDouble();
				double tol = avg*0.2;
				List<Atom> byleft = cf.atoms().map(a->Tuple.of(-a.getAtomCoordinates().getX(),a).withKComparator()).sorted().map(t->t.v()).collect(Collectors.toList());
				List<Atom> bytop = cf.atoms().map(a->Tuple.of(-a.getAtomCoordinates().getY(),a).withKComparator()).sorted().map(t->t.v()).collect(Collectors.toList());
				Atom bottomMost=bytop.get(bytop.size()-1);
				Atom topMost=bytop.get(0);
				Atom rightMost=byleft.get(byleft.size()-1);
				Atom leftMost=byleft.get(0);
				
				List<Atom> leftContend = byleft.stream().filter(ba->ba.getAtomCoordinates().getX()<leftMost.getAtomCoordinates().getX()+tol).collect(Collectors.toList());
				List<Atom> rightContend = byleft.stream().filter(ba->ba.getAtomCoordinates().getX()>rightMost.getAtomCoordinates().getX()-tol).collect(Collectors.toList());
				List<Atom> topContend = bytop.stream().filter(ba->ba.getAtomCoordinates().getY()<topMost.getAtomCoordinates().getY()+tol).collect(Collectors.toList());
				List<Atom> bottomContend = bytop.stream().filter(ba->ba.getAtomCoordinates().getY()>bottomMost.getAtomCoordinates().getY()-tol).collect(Collectors.toList());
				
				Chemical ch=cf;
				mset.forEach(ss->{
					List<Atom> attempt;
					double dx=0;
					double dy=0;
					switch(ss){
						case BOTTOM:
							attempt=bottomContend;
							dy=avg;
							break;
						case LEFT:
							dx=-avg;
							attempt=leftContend;
							break;
						case RIGHT:
							dx=avg;
							attempt=rightContend;
							break;
						case TOP:
							attempt=topContend;
							dy=-avg;
							break;
						default:
							attempt=new ArrayList<>();
							break;
					}
					Atom neigh=null;
					if(attempt.size()==1){
						neigh=attempt.get(0);
					}else if(attempt.size()>1){
						neigh=attempt.get(0);
					}
					
					if(neigh!=null){
						Atom nat=ch.addAtom("C",neigh.getAtomCoordinates().getX()+dx,neigh.getAtomCoordinates().getY()+dy);
						ch.addBond(nat,neigh,BondType.SINGLE);
						nat.setImplicitHCount(null);
						neigh.setImplicitHCount(null);
					}
				});
			}
			}
			res.c=cf;
			cf.clearAtomMaps();


//			System.out.println(correlationToClean(cf));


		} catch (Exception e) {
			res.type=FixType.NULL;
			
			// TODO Auto-generated catch block
//									e.printStackTrace();
			return res;
		}
		return res;
	}

	public static String atFeat(Atom a, int r){
		if(r==0){
			return a.getSymbol();
		}
		return a.getSymbol() + a.getNeighbors().stream().map(n->n.bondTo(a).get().getBondType().getSymbol() + atFeat(n,r-1)).sorted().collect(Collectors.joining(","));
	}

	public static List<String> getAllFeats(Chemical c, int r){
		return c.atoms().map(a->atFeat(a,r)).sorted().collect(Collectors.toList());
	}
	
	public static List<String> getAllBondFeats(Chemical c){
		return c.bonds()
				.map(b->Stream.of(b.getAtom1(),b.getAtom2())
						      .map(at->at.getSymbol())
						      .sorted()
						      .collect(Collectors.joining(b.getBondType().getSymbol()+""
						    		                     )
						    		  )
			        )
				.sorted()
				.collect(Collectors.toList());
	}

	public static List<String> compareFeats(Chemical c1, Chemical c2, int r){
		List<String> f1= getAllFeats(c1,r);
		List<String> f2= getAllFeats(c2,r);


		return getDiffs(f1,f2);

	}
	
	public static List<String> compareBondFeats(Chemical c1, Chemical c2){
		List<String> f1= getAllBondFeats(c1);
		List<String> f2= getAllBondFeats(c2);


		return getDiffs(f1,f2);

	}

	public static List<String> getDiffs(List<String> f1, List<String> f2){
		List<String> diff = new ArrayList<String>();

		int i=0;
		int j=0;
		while(true){
			if(i>=f1.size() && j<f2.size()){
				diff.add("+" + f2.get(j));
				j++;
			}else if(j>=f2.size() && i<f1.size()){
				diff.add("-" + f1.get(i));
				i++;
			}else if(j<f2.size() && i<f1.size()){
				String sf1=f1.get(i);
				String sf2=f2.get(j);
				int kk=sf1.compareTo(sf2);
				if(kk<0){
					diff.add("-" + sf1);
					i++;
				}else if(kk>0){
					diff.add("+" + sf2);
					j++;
				}else{
					j++;
					i++;
				}
			}else{
				break;
			}

		}
		return diff;

	}

	public static void main(String[] h) throws Exception{
		Chemical c= Chemical.parse("CCCCC");
		
		System.out.println(correlationToClean(c));
		System.out.println(c.toMol());
	}

}
