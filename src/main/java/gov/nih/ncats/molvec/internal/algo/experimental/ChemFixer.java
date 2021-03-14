package gov.nih.ncats.molvec.internal.algo.experimental;

import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import gov.nih.ncats.common.stream.StreamUtil;
import gov.nih.ncats.molvec.internal.algo.Tuple;
import gov.nih.ncats.molvec.internal.util.GeomUtil;
import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.Bond;
import gov.nih.ncats.molwitch.Bond.BondType;
import gov.nih.ncats.molwitch.Bond.Stereo;
import gov.nih.ncats.molwitch.Chemical;

public class ChemFixer {

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

			if(ds<avg*avg*(1.3*1.3) && ds> avg*avg*0.68*0.68){
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

				Bond nb=cop.addBond(ca1, ca2, BondType.QUADRUPLE);

				if(nb.isInRing() && ca1.getSmallestRingSize()>=6 && ca2.getSmallestRingSize()>=6){
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
					if(pent.getAtomCount()==5){
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

			if(ds<avg*avg*0.33*0.33){
				return true;
			}
			return false;
		})
		.collect(Collectors.toList())
		.forEach(t->{

			if(t.k().getSymbol().equals("C") && t.k().getBondCount()==1 && t.v().getBondCount()>1){
				Atom oa=t.k().getNeighbors().get(0);
				c.addBond(oa,t.v(),BondType.SINGLE);
				c.removeAtom(t.k());
			}else{
				t=t.swap();
				if(t.k().getSymbol().equals("C") && t.k().getBondCount()==1 && t.v().getBondCount()>1){
					Atom oa=t.k().getNeighbors().get(0);
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
			for(int j=0;j<c2.getAtomCount();j++){
				Atom a2=c2.getAtom(j);
				if(a2.getSymbol().equals("H"))continue;
				double[] xy1=a1.getAtomCoordinates().xy();
				double[] xy2=a2.getAtomCoordinates().xy();

				double sq = Math.abs(xy1[0]-xy2[0]) + Math.abs(xy1[1]-xy2[1]) ;

				if(sq<minsq){
					minsq=sq;
					mini=i;
					minj=j;
				}
			}
		}
		return Tuple.of(new int[]{mini,minj},Math.sqrt(minsq));

	}
	public static Chemical dumbClean(Chemical c){
		makeMissingBonds(c);
		;

		{
			c.atoms()
			.filter(at->at.getSymbol().equals("H") || at.getSymbol().equals("P")  || at.getSymbol().equals("Cl"))
			.filter(at->at.getSmallestRingSize()==6)
			.forEach(at->{
				at.setAtomicNumber(6);
				at.setImplicitHCount(null);
			});
		}
		{
			List<Atom> rem = c.atoms()
					.filter(at->at.getSymbol().equals("I"))
					.filter(a->a.getBondCount()>=1)
					.filter(at->{
						Atom n1=at.getNeighbors().get(0);
						if(!n1.getSymbol().equals("C"))return true;

						return false;
					})
					.collect(Collectors.toList());
			rem.forEach(at->{
				List<Atom> nl=at.getNeighbors();
				c.removeAtom(at);
				nl.forEach(na->{
					na.setImplicitHCount(null);
				});
			});
		}
		{
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

		{
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

		{
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

		{
			c.atoms()
			.filter(at->at.getSymbol().equals("F"))
			.filter(at->at.getBondCount()>0)
			.forEach(at->{
				List<Atom> fs=at.getNeighbors().get(0).getNeighbors().stream().filter(nn->nn.getBondCount()==1)
						.filter(a2->a2.getSymbol().equals("C") || a2.getSymbol().equals("F") || a2.getSymbol().equals("I")|| a2.getSymbol().equals("B"))
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
						Atom an=c.addAtom("F", ca.getAtomCoordinates().getX()+1,ca.getAtomCoordinates().getY());
						c.addBond(ca,an, BondType.SINGLE);
					}else{
						fs.forEach(att->{
							att.setAtomicNumber(9);	
							att.setImplicitHCount(null);
						});	
					}

				}


			});
		}

		{
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
					c.removeBond(b);
					b.getAtom1().setImplicitHCount(null);
					b.getAtom2().setImplicitHCount(null);
				});;

			});
		}

		{
			List<Atom> rem = c.atoms()
					.filter(at->at.getSymbol().equals("I"))
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
		{
			long fcount=c.atoms().filter(at->at.getSymbol().equals("F")).count();
			if(fcount>0){
				c.atoms()
				.filter(at->at.getSymbol().equals("I"))
				.filter(at->at.getBondCount()==0 || at.getNeighbors().get(0).getBonds().stream().filter(bb->bb.getBondType().getOrder()==2).count()>0)
				.forEach(at->{
					at.setAtomicNumber(9);
				});
			}

		}

		{
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
				List<Atom> nl=at.getNeighbors();
				c.removeAtom(at);
				nl.forEach(na->{
					na.setImplicitHCount(null);
				});
			});
		}


		{

			c.atoms()
			.filter(at->at.getSymbol().equals("O"))
			.filter(at->at.getBondCount()==1)
			.filter(at->at.getBonds().get(0).getBondType().equals(BondType.SINGLE))
			//					.filter(at->Stereo.DOWN.equals(at.getBonds().get(0).getStereo()) ||at.getBonds().get(0).getStereo().equals(Stereo.DOWN_INVERTED))
			.filter(at->at.getNeighbors().get(0).getBonds().stream().mapToInt(b->b.getBondType().getOrder()).sum()<4)
			.filter(at->at.getNeighbors().get(0).getNeighbors().stream()
					.filter(nn->nn.getSymbol().equals("N"))
					.filter(nn->nn.getBonds().stream().allMatch(bn->bn.getBondType().equals(BondType.SINGLE))).count()>0)
			.forEach(at->{
				at.getBonds().get(0).setBondType(BondType.DOUBLE);
				at.setImplicitHCount(null);
				at.getNeighbors().get(0).setImplicitHCount(null);
			});
		}

		{

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



		{

			c.atoms()
			.filter(at->at.getSymbol().equals("N"))
			.forEach(na->{
				if(na.getBonds().stream().mapToInt(b->b.getBondType().getOrder()).sum()>3){
					na.setAtomicNumber(6);
					na.setImplicitHCount(null);
				}
			});
		}
		{
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
		{
			c.atoms()
			.filter(at->at.getSymbol().equals("B"))
			.filter(at->at.getBondCount()==1)
			.filter(at->at.getNeighbors().get(0).isInRing() && at.getNeighbors().get(0).getSmallestRingSize()==6)
			.forEach(na->{
				na.setAtomicNumber(9);
			});
		}
		//!@@@@@@@@@@@@@@
		{


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
		{

			double avg= c.bonds().mapToDouble(b->b.getBondLength()).average().getAsDouble();
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
						c.removeAtom(latom);
						c.addBond(patom, nlatom, BondType.SINGLE);

					}

				}

			});
		}
		{

			double avg= c.bonds().mapToDouble(b->b.getBondLength()).average().getAsDouble();
			c.bonds()
			.filter(bb->bb.getBondLength()<avg*0.3)
			//					.filter(bb->bb.isInRing())
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
					c.removeAtom(toD);
					c.addBond(nat.get(0),nat.get(1), BondType.SINGLE);
					nat.get(0).setImplicitHCount(null);
					nat.get(1).setImplicitHCount(null);
				}

			});
		}

		{

			//			double avg= c.bonds().mapToDouble(b->b.getBondLength()).average().getAsDouble();
			c.bonds()
			.filter(bb->bb.getBondType().equals(BondType.SINGLE))
			.filter(bb->bb.getStereo().equals(Stereo.DOWN) || bb.getStereo().equals(Stereo.DOWN_INVERTED)  )
			.filter(bb->bb.isInRing())
			.filter(bb->bb.getAtom1().getSmallestRingSize()==5)
			.filter(bb->Stream.of(bb.getAtom1(),bb.getAtom2()).flatMap(at->at.getBonds().stream()).filter(b1->b1.getBondType().equals(BondType.DOUBLE)).count()==0)
			.forEach(bb->{
				bb.setBondType(BondType.DOUBLE);

			});
		}

		c.atoms().forEach(at->{
			if(at.getCharge()!=0){
				at.setCharge(0);
				at.setImplicitHCount(null);
			}

		});


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
					c.removeAtom(at);
				}
		});



		{
			//remove dumb triangles
			c.atoms()
			.filter(at->at.getSmallestRingSize()==3)
			.filter(at->at.getSymbol().equals("C"))
			.flatMap(at->at.getBonds().stream())
			.filter(b->b.getBondType().equals(BondType.SINGLE))
			.filter(b->b.getStereo().equals(Stereo.DOWN)||b.getStereo().equals(Stereo.DOWN_INVERTED))
			.distinct()
			.filter(b->b.getAtom1().getSmallestRingSize()==3 &&b.getAtom2().getSmallestRingSize()==3)
			.forEach(b->{
				c.removeBond(b);
				if(b.getAtom1().getSmallestRingSize()==6 &&b.getAtom2().getSmallestRingSize()==6 ){

				}else{
					c.addBond(b);
				}
			});
		}
		{
			//remove dumb double bonds
			c.atoms()
			.filter(at->at.getSmallestRingSize()==6)
			.filter(at->at.getSymbol().equals("C"))

			.filter(b->b.getBonds().stream()
					.filter(bn->bn.getBondType().equals(BondType.DOUBLE))
					.count()==2)
			.flatMap(a->a.getBonds().stream().filter(b->b.getBondType().getOrder()==2).map(b->b.getOtherAtom(a)))
			.filter(a->a.getBondCount()==1)
			.forEach(a->{
				a.getNeighbors().get(0).setAtomicNumber(7);
				c.removeAtom(a);
			});
		}


		makeMissingBonds(c);

		{
			//			
			double avg= c.bonds().mapToDouble(b->b.getBondLength()).average().getAsDouble();
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
				c.removeAtom(at1);
				//				b.setBondType(BondType.DOUBLE);
			});
			;
		}
		{
			double avg= c.bonds().mapToDouble(b->b.getBondLength()).average().getAsDouble();
			c.atoms()
			.filter(at->at.getSmallestRingSize()==6)
			.flatMap(at->at.getBonds().stream())
			.filter(b->b.getBondType().equals(BondType.SINGLE))
			.filter(b->Math.abs(b.getAtom1().getAtomCoordinates().getX() - b.getAtom2().getAtomCoordinates().getX())<avg*0.2 ||
					Math.abs(b.getAtom1().getAtomCoordinates().getY() - b.getAtom2().getAtomCoordinates().getY())<avg*0.2
					)
			.filter(b->b.getAtom1().getSmallestRingSize()==6 &&b.getAtom2().getSmallestRingSize()==6)
			.filter(b->b.getAtom1().getBonds().stream().filter(bb->!bb.getBondType().equals(BondType.SINGLE)).count()==0)
			.filter(b->b.getAtom2().getBonds().stream().filter(bb->!bb.getBondType().equals(BondType.SINGLE)).count()==0)
			//			.peek(b->System.out.println(b))
			.filter(b->b.getAtom1().getNeighbors().stream().filter(at->at.getSmallestRingSize()==6).flatMap(aa->aa.getBonds().stream()).filter(bb->bb.getBondType().equals(BondType.DOUBLE)).count()>=1)
			.filter(b->b.getAtom2().getNeighbors().stream().filter(at->at.getSmallestRingSize()==6).flatMap(aa->aa.getBonds().stream()).filter(bb->bb.getBondType().equals(BondType.DOUBLE)).count()>=1)

			.forEach(b->{
				b.setBondType(BondType.DOUBLE);
			});
			;
		}

		{
			double avg= c.bonds().mapToDouble(b->b.getBondLength()).average().getAsDouble();
			c.atoms()
			.filter(at->at.getSymbol().equals("C"))
			.filter(at->at.getBondCount()==2)
			.forEach(at->{
				Atom a1=at.getNeighbors().get(0);
				Atom a2=at.getNeighbors().get(1);
				if(a1.bondTo(a2).isPresent()){
					double hyp = a1.bondTo(a2).get().getBondLength();
					//triangle
					if(at.getBonds().get(0).getBondLength()+ at.getBonds().get(1).getBondLength() < avg*1.7 || at.getBonds().get(0).getBondLength()+ at.getBonds().get(1).getBondLength() < hyp*1.7 ||
							at.getBonds().get(0).getBondLength() < hyp*0.6 || at.getBonds().get(1).getBondLength() < hyp*0.6){
						Bond bb=a1.getBonds().stream()
								.filter(b->b.getOtherAtom(a1).equals(a2))
								.findFirst().get();
						if(bb.getBondType().getOrder() == 1 && !bb.getStereo().equals(Stereo.NONE)){
							bb.setBondType(BondType.DOUBLE);
							bb.getAtom1().setImplicitHCount(null);
							bb.getAtom2().setImplicitHCount(null);
						}
						c.removeAtom(at);
						a1.setImplicitHCount(null);
						a2.setImplicitHCount(null);
					}
				}
			});
		}

		{
			double avg= c.bonds().mapToDouble(b->b.getBondLength()).average().getAsDouble();
			Set<Bond> obond = getOverlappingBonds(c).stream().flatMap(t->Stream.of(t.k(),t.v())).collect(Collectors.toSet());
			c.bonds()
			.filter(b->b.getBondLength()>avg*1.3)
			.filter(b->b.isInRing())
			.filter(b->b.getAtom1().getSmallestRingSize()==5 && b.getAtom2().getSmallestRingSize()==5 )
			.filter(b->b.getAtom1().getSymbol().equals("C") && b.getAtom2().getSymbol().equals("C") )
			.filter(b->b.getAtom1().getSymbol().equals("C") && b.getAtom2().getSymbol().equals("C") )
			.filter(b->!obond.contains(b))
			.forEach(b->{
				double fy=(b.getAtom1().getAtomCoordinates().getY()+b.getAtom2().getAtomCoordinates().getY())/2;
				double fx=(b.getAtom1().getAtomCoordinates().getX()+b.getAtom2().getAtomCoordinates().getX())/2;
				Atom mat=c.addAtom("C", fx,fy);
				c.removeBond(b);
				Bond b1=c.addBond(b.getAtom1(), mat, BondType.SINGLE);
				Bond b2=c.addBond(b.getAtom2(), mat, BondType.SINGLE);
				if(b.getAtom1().getBonds().stream().filter(bb->bb.getBondType().getOrder()==2).count()==1){
					b2.setBondType(BondType.DOUBLE);
				}else if(b.getAtom2().getBonds().stream().filter(bb->bb.getBondType().getOrder()==2).count()==1){
					b1.setBondType(BondType.DOUBLE);
				}
				mat.setImplicitHCount(null);
				b.getAtom1().setImplicitHCount(null);
				b.getAtom2().setImplicitHCount(null);
			});
		}


		combineCloseBonds(c);
		c.atoms().forEach(at->{
			at.setImplicitHCount(null);		
		});

		Set<Bond> sb=new HashSet<Bond>();
		for(Bond b: c.getBonds()){
			if(sb.contains(b)){
				c.removeBond(b);
			}
			sb.add(b);
		}



		//				;
		//		c.bonds().filter(b->!sb.contains(b))
		//		.collect(Collectors.toList()).forEach(b->c.removeBond(b));

		return c;
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
					c.removeAtom(at);
				}
		});
	}

	public static Tuple<Chemical,Boolean> stitchChemical(Chemical cf){
		double avg= cf.bonds().mapToDouble(b->b.getBondLength())
				.average().getAsDouble();
		boolean stitched=false;
		for(int li=0;li<20;li++){
			Chemical c=cf.copy();
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
								mini=c1.getAtom(pair.k()[0]).getAtomToAtomMap().getAsInt()-1;
								minj=c2.getAtom(pair.k()[1]).getAtomToAtomMap().getAsInt()-1;
								mind=pair.v();
							}
						}
					}
					//							System.out.println("mini:" + mini);
					//							System.out.println("minj:" + minj);
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
		ChemFixResult res = new ChemFixResult();

		try {


			dumbClean(cf);	
			Tuple<Chemical, Boolean> tup=stitchChemical(cf);
			cf=tup.k();
			if(tup.v()){
				res.type=FixType.MERGED;
			}

			setHStereo(cf);

			String[] lines=cf.toMol().split("\n");
			for(int ii=4;ii<4+Integer.parseInt(lines[3].substring(0,3).trim());ii++){
				lines[ii]=lines[ii].substring(0,33) + "  0  0  0  0  0  0  0  0  0  0  0  0";
			}
			String mm = Arrays.stream(lines).collect(Collectors.joining("\n"));
			cf= Chemical.parse(mm);
			res.c=cf;
			cf.clearAtomMaps();




		} catch (Exception e) {
			// TODO Auto-generated catch block
			//						e.printStackTrace();
			return null;
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

	public static List<String> compareFeats(Chemical c1, Chemical c2, int r){
		List<String> f1= getAllFeats(c1,r);
		List<String> f2= getAllFeats(c2,r);


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

	public static void main(String[] h){
		List<String> l1=Arrays.stream("A,B,C,D,E,F,F,G".split(",")).sorted().collect(Collectors.toList());
		List<String> l2=Arrays.stream("A,B,C,D,D,E,F,G,G".split(",")).sorted().collect(Collectors.toList());
		getDiffs(l1,l2)
		.forEach(c->System.out.println(c));
	}

}
