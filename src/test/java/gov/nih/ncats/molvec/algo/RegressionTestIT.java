package gov.nih.ncats.molvec.algo;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.RenderedImage;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Random;
import java.util.Set;
import java.util.UUID;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageWriteParam;
import javax.imageio.ImageWriter;
import javax.imageio.plugins.jpeg.JPEGImageWriteParam;
import javax.imageio.stream.FileImageOutputStream;

import org.junit.Test;

import com.mortennobel.imagescaling.ResampleOp;

import gov.nih.ncats.chemkit.api.Atom;
import gov.nih.ncats.chemkit.api.AtomCoordinates;
import gov.nih.ncats.chemkit.api.Bond;
import gov.nih.ncats.chemkit.api.Bond.Stereo;
import gov.nih.ncats.chemkit.api.Chemical;
import gov.nih.ncats.chemkit.api.ChemicalBuilder;
import gov.nih.ncats.chemkit.api.inchi.Inchi;
import gov.nih.ncats.molvec.Molvec;
import gov.nih.ncats.molvec.algo.ShellCommandRunner.Monitor;
import gov.nih.ncats.molvec.image.ImageUtil;
import gov.nih.ncats.molvec.util.GeomUtil;

public class RegressionTestIT {
	
	
	
	private static boolean DO_ALIGN = false;
	private static boolean EXPORT_CORRECT = false;
	
	private static String exportDir = "/home/tyler/workspace/results";
	

	public static enum Result{
		CORRECT_FULL_INCHI,
		CORRECT_STEREO_INSENSITIVE_INCHI,
		LARGEST_FRAGMENT_CORRECT_FULL_INCHI,
		LARGEST_FRAGMENT_CORRECT_STEREO_INSENSITIVE_INCHI,
		FORMULA_CORRECT,
		LARGEST_FRAGMENT_FORMULA_CORRECT,

		ATOMS_RIGHT_WRONG_BONDS,
		LARGEST_FRAGMENT_ATOM_RIGHT_WRONG_BONDS,
		ATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY,
		LARGEST_FRAGMENT_ATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY,
		WEIRD_SOURCE,
		RIGHT_HEVAY_ATOMS,
		RIGHT_BONDS,
		LARGEST_FRAGMENT_RIGHT_BONDS,
		INCORRECT,
		FOUND_NOTHING,
		ERROR,
		TIMEOUT
	}
	public static class TestResult{
		public Result result;
		public long time;
		public double RMSE=Double.POSITIVE_INFINITY;
		
		public static TestResult of(Result r, long ms, double rmse){
			TestResult tr = new TestResult();
			tr.result=r;
			tr.time=ms;
			tr.RMSE=rmse;
			return tr;
		}
		public static TestResult of(Result r, long ms){
			TestResult tr = new TestResult();
			tr.result=r;
			tr.time=ms;
			return tr;
		}
	}
	

	private File getFile(String fname){
		ClassLoader classLoader = getClass().getClassLoader();
		return new File(classLoader.getResource(fname).getFile());
		
	}
	

	public static Chemical combineChemicals(Chemical c1, Chemical c2){
		ChemicalBuilder nc = c1.copy().toBuilder();
		
		Map<Atom,Atom> oldToNew = new HashMap<>();
		
		for(int i=0;i<c1.getAtomCount();i++){
			Atom aa=nc.atomAt(i);
			AtomCoordinates ac=aa.getAtomCoordinates();
			aa.setAtomCoordinates(AtomCoordinates.valueOf(ac.getX(), ac.getY(), ac.getZ().orElse(0)));
		}
		
		c2.atoms()
		  .forEach(a->{
			  AtomCoordinates ac=a.getAtomCoordinates();
			  
			  Atom na=nc.addAtom(a.getSymbol(), ac.getX(), ac.getY(), ac.getZ().orElse(0));
			  oldToNew.put(a, na);
			  na.setCharge(a.getCharge());
			  na.setMassNumber(a.getMassNumber());
			  na.setAtomCoordinates(AtomCoordinates.valueOf(ac.getX(), ac.getY(), ac.getZ().orElse(0)));
		  });
		
		c2.bonds()
		  .forEach(b->{
			  Atom na1=oldToNew.get(b.getAtom1());
			  Atom na2=oldToNew.get(b.getAtom2());
			  Bond nb=nc.addBond(na1,na2, b.getBondType());
			  if(b.getStereo()!=null && b.getStereo()!=Stereo.NONE){
				  nb.setStereo(b.getStereo());  
			  }
			  //
		  });
//		
		
		return nc.build();
	}
	

	

	public static Chemical getOSRAChemical(File f) throws IOException, InterruptedException{
		StringBuilder resp = new StringBuilder();
		AtomicBoolean done=new AtomicBoolean(false);
		
		Monitor m=(new ShellCommandRunner.Builder()).activeDir("./")
	               .command("osra", "-f sdf", f.getAbsolutePath())
	               .build()
	               .run();
		m.onError(l->{
			try{
				System.err.println(l);
			m.kill();
			}catch(Exception e){
				e.printStackTrace();
			}
		});
		m.onInput(l->{resp.append(l + "\n");});
		m.onKilled(k->{done.set(true);});
	
		while(!done.get()){
			Thread.sleep(5);
		}
		
		StringBuilder sbnew = new StringBuilder();
		
		Chemical fc= Arrays.stream(resp.toString().split("\n"))
		      .map(l->{
		    	  sbnew.append(l+"\n");
		    	  if(l.equals("$$$$")){
		    		  String f1=removeSGroupStuff(sbnew.toString().replace(" *   ", " C   "));
		    		  sbnew.setLength(0);
		    		  try {
						return Chemical.parseMol(f1);
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
		    	  }
		    	  return null;
		      })
		      .filter(c->c!=null)
		      .reduce((c1,c2)->combineChemicals(c1,c2))
		      .orElse(new ChemicalBuilder().build());
		
		return fc;
	}

	public static Chemical getImagoChemical(File f) throws IOException, InterruptedException{
		AtomicBoolean done=new AtomicBoolean(false);
		
		
		
		
		String fname = UUID.randomUUID().toString();
		
		File imageFile = File.createTempFile(fname, ".png");
		File molFile = File.createTempFile(fname, ".mol");
		imageFile=stdResize(f, imageFile, 1, Interpolation.NEAREST_NEIGHBOR,1);
				
		String tmpFileNameImage = imageFile.getAbsolutePath();
		String tmpFileNameMol = molFile.getAbsolutePath();
		
		
		//System.out.println(raw1);
		//if(true)return new ChemicalBuilder().build();
		
		Monitor m=(new ShellCommandRunner.Builder())
	               .command("./imago_console", tmpFileNameImage, "-o", tmpFileNameMol)
	               .build()
	               .run();
		m.onError(l->{
			try{
				System.err.println("err:" + l);
			//m.kill();
			}catch(Exception e){
				e.printStackTrace();
			}
		});
		m.onKilled(k->{done.set(true);});
	
		while(!done.get()){
			Thread.sleep(5);
		}

		StringBuilder sbnew = new StringBuilder();
		try(Stream<String> slines=Files.lines(Paths.get(tmpFileNameMol))){
			Chemical fc=Stream.concat(slines, Stream.of("$$$$"))
				 .map(l->{
			    	  sbnew.append(l+"\n");
			    	  if(l.equals("$$$$")){
			    		  String f1=removeSGroupStuff(sbnew.toString().replace(" *   ", " C   "));
			    		  sbnew.setLength(0);
			    		  try {
							return Chemical.parseMol(f1);
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
			    	  }
			    	  return null;
			      })
			      .filter(c->c!=null)
			      .reduce((c1,c2)->combineChemicals(c1,c2))
			      .orElse(new ChemicalBuilder().build());
			
			fc.atoms()
			 .filter(a->a.hasAromaticBond())
			 .filter(a->a.getSymbol().equals("C"))
			 .forEach(aa->aa.setImplicitHCount(Math.max(0,3-aa.getBonds().size())));
			fc.kekulize();
			
			return fc;
		}
		
	}
	
	private static String removeSGroupStuff(String mol){
		AtomicBoolean lastWasAlias = new AtomicBoolean(false);
		Stream<String> stream = Arrays.stream(mol.split("\n"));
		
			return stream
		
				.filter(l->!l.contains("SUP"))
				.filter(l->!l.contains("SGROUP"))
				.filter(l->!l.contains("M  S"))
				.filter(l->{
					if(lastWasAlias.get()){
						lastWasAlias.set(false);
						return false;
					}
					return true;
				})
				.filter(l->{
					if(l.startsWith("A")){
						lastWasAlias.set(true);
						return false;
					}
					lastWasAlias.set(false);
					return true;
				})
				.collect(Collectors.joining("\n"));
		
	}
	private static Chemical getCleanChemical(String mol) throws IOException{
		Chemical nc= Chemical.parseMol(mol);
		
		
		Set<String> metals = Stream.of("Na","K","Li","Mg", "Pt").collect(Collectors.toSet());
		
		Set<String> bondsToRemove = nc.bonds()
		.filter(b->{
			if(metals.contains(b.getAtom1().getSymbol())){
				
				b.getAtom1().setCharge(b.getAtom1().getCharge()+1);
				b.getAtom2().setCharge(b.getAtom2().getCharge()-1);
				return true;
			}else if(metals.contains(b.getAtom2().getSymbol())){
				
				b.getAtom2().setCharge(b.getAtom2().getCharge()+1);
				b.getAtom1().setCharge(b.getAtom1().getCharge()-1);
				return true;
			}
			return false;
		})
		.map(b-> Tuple.of(b.getAtom1(),b.getAtom2()))
		.map(Tuple.vmap(a->a.getAtomIndexInParent()+1))
		.map(Tuple.kmap(a->a.getAtomIndexInParent()+1))
		.map(Tuple.vmap(i->("   " + i)))
		.map(Tuple.kmap(i->("   " + i)))
		.map(Tuple.vmap(i->i.toString().substring(i.toString().length()-3)))
		.map(Tuple.kmap(i->i.toString().substring(i.toString().length()-3)))
		.map(t->t.k()+t.v())
		//.peek(s->System.out.println(s))
		.collect(Collectors.toSet());
		
		if(!bondsToRemove.isEmpty()){
			String[] lines2 = nc.toMol().split("\n");
			//System.out.println("OLD:" + nc.toMol());
			String padBonds = "   " + (nc.getBondCount()-bondsToRemove.size());
			padBonds=padBonds.substring(padBonds.length()-3);
			
			lines2[3]=lines2[3].substring(0, 3) + padBonds + lines2[3].substring(6);
			String mol2=Arrays.stream(lines2)
					          .filter(bl->!bondsToRemove.contains((bl+"      ").substring(0, 6)))
					          .collect(Collectors.joining("\n"));
			
			nc= Chemical.parseMol(mol2);
			//System.out.println("NEW:" + nc.toMol());
		}
		  
		
		return nc;
	}
	
	 /**
     * Simple morgan's algorithm for graph invariants. This requires k*N operations
     * where k is a constant that is large enough to "absorb" the whole graph (13 here).
     * 
     * @param m
     * @return
     */
    public static long[] morgans(Chemical m){
    	int MAX_ROUND = 13;
    	int[] atno = new int[m.getAtomCount()];
        for (int i = 0; i < atno.length; ++i) {
            Atom a = m.getAtom(i);
            atno[i]=a.getAtomicNumber();
        }
        long[] rank;
        {
            
            long[][] hash = new long[MAX_ROUND][atno.length];
            for (int i = 0; i < atno.length; ++i)
                hash[0][i] = atno[i];

            int round = 1;
            for (; round < MAX_ROUND; ++round) {
                int p = round - 1;
                for (int i = 0; i < atno.length; ++i) {
                    Atom a = m.getAtom(i);
                    long ha = hash[p][i];
                    ha+=a.getBonds()
                     .stream()
                     .mapToLong(b->{
                    	 Atom oa=b.getOtherAtom(a);
                    	 int k = oa.getAtomIndexInParent();
                    	 return (1 << oa.getImplicitHCount()) + hash[p][k];
                     })
                     .sum();
                    if (ha < 0) {
                        ha = hash[round-1][i];
                    }
                    hash[round][i] = ha;
                }
            }
            rank = hash[round-1];
        }
        return rank;
    }

	
	//need to align them:
	// 1. Resize to same size (easy enough based on bond length probably)
	// 2. Center?
	// 3. Something else ...
	
	public static double align(Chemical c1, Chemical c2) throws IOException{
		if(!DO_ALIGN)return Double.POSITIVE_INFINITY;
		c1.makeHydrogensImplicit();
		c2.makeHydrogensImplicit();
		
		List<Chemical> c1chems = c1.connectedComponentsAsStream().collect(Collectors.toList());
		
		
		
		if(c1chems.size()>1){
			
			Map<String,Chemical> c2chems = c2.connectedComponentsAsStream()
					.map(cc->{
						String ikey;
						try{
							ikey=Inchi.asStdInchi(cc).getKey().split("-")[0];
						}catch(Exception e){
							e.printStackTrace();
							ikey="ERROR";
						}
						return Tuple.of(ikey,cc);	
					})
					.collect(Tuple.toMap());
			
			return c1chems.stream().map(cc->{
						String ikey;
						try{
							ikey=Inchi.asStdInchi(cc).getKey().split("-")[0];
						}catch(Exception e){
							ikey="ERROR";
						}
						return Tuple.of(ikey,cc);	
					})
					.map(Tuple.kmap(ik->c2chems.get(ik)))
					.map(t->{
						try{
							double rmse= align(t.v(),t.k());
							double unrmse= rmse*rmse*t.v().getAtomCount();
							return Tuple.of(t.v().getAtomCount(),unrmse);
						}catch(Exception e){
							throw new RuntimeException(e);
						}
					})
					.reduce((a,b)->Tuple.of(a.k()+b.k(),a.v()+b.v()))
					.map(t->{
						return Math.sqrt(t.v()/t.k());
					})
					.orElse(0.0);
			
		}
		

		double b1avg = c1.bonds().mapToDouble(b->b.getBondLength()).average().orElse(1);

		AffineTransform at2= new AffineTransform();
		at2.scale(1/b1avg, 1/b1avg);
		c1.atoms()
		  .map(a->Tuple.of(a,a.getAtomCoordinates()))
		  .map(Tuple.vmap(ac->asPoint(ac)))
		  .map(Tuple.vmap(p->at2.transform(p, null)))
		  .map(Tuple.vmap(p->fromPoint(p)))
		  .forEach(t->{
			  t.k().setAtomCoordinates(t.v());
		  });
				
		long[] l1= morgans(c1);
		long[] l2= morgans(c2);
		
		Set<Long> morganIgnore = new HashSet<Long>();
		
		
	
		
		Map<Long,List<Integer>> omap =IntStream.range(0,l1.length)
				 .mapToObj(i->Tuple.of(i, c1.getAtom(i)))
				 .map(Tuple.kmap(i->Tuple.of(i,l1[i])))
				 .peek(t->{
					 if(t.v().getAtomToAtomMap().isPresent()){
						 morganIgnore.add(t.k().v());
					 }
				 })
				 .map(t->t.k())
				 .map(t->t.swap())
				 .collect(Tuple.toGroupedMap());
		
		
		
		
		
		if(!morganIgnore.isEmpty()){
			IntStream.range(0,l1.length)
					.mapToObj(i->Tuple.of(i,c1.getAtom(i)))
					.collect(Collectors.toList())
					.stream()
					.filter(t->morganIgnore.contains(l1[t.k()]))
					.forEach(t->{
						c1.removeAtom(t.v());
					});
			
			IntStream.range(0,l2.length)
						.mapToObj(i->Tuple.of(i,c2.getAtom(i)))
						.collect(Collectors.toList())
						.stream()
						.filter(t->morganIgnore.contains(l2[t.k()]))
						.forEach(t->{
							c2.removeAtom(t.v());
						});
			System.out.println("IGNORE IT");
			return align(c1,c2);
		}
		
		List<Tuple<Integer,Integer>> sameIndex =				
				IntStream.range(0,l2.length)
				 		 .mapToObj(i->Tuple.of(i,l2[i]))
				 		 .filter(t->!morganIgnore.contains(t.k()))
				 		 .map(Tuple.vmap(m->omap.get(m)))
				 		 .filter(t->t.v()!=null)
				 		 .filter(t->t.v().size()==1)
				 		 .map(Tuple.vmap(ol->ol.get(0)))
				 		 .map(t->t.swap())
				 		 .collect(Collectors.toList());
		
		Point2D cpt1=
				IntStream.range(0,c1.getAtomCount())
				.mapToObj(i->Tuple.of(i,c1.getAtom(i)))
				.map(t->t.v())
				  .map(a->a.getAtomCoordinates())
				  .map(ac->asPoint(ac))
				  .collect(GeomUtil.averagePoint());
				
		Point2D cpt2=IntStream.range(0,l2.length)
				.mapToObj(i->Tuple.of(i,c2.getAtom(i)))
				.map(t->t.v())
				.map(a->a.getAtomCoordinates())
				  .map(ac->asPoint(ac))
				  .collect(GeomUtil.averagePoint());
		
		
		double tdx1 = cpt1.getX();
		double tdy1 = cpt1.getY();
		double tdx2 = cpt2.getX();
		double tdy2 = cpt2.getY();
		
		double s=1;
		
		
		if(sameIndex.size()>=2){
			
			System.out.println("Matched:" + sameIndex.size());
			
			AtomCoordinates center2=fromPoint(cpt2);
			double cscale = sameIndex.stream()
					.mapToDouble(t->c2.getAtom(t.v()).getAtomCoordinates().distanceSquaredTo(center2))
					.sum()
					;
			
			
			s = sameIndex.stream()
					            .map(Tuple.kmap(i->c1.getAtom(i).getAtomCoordinates()))
					            .map(Tuple.vmap(i->c2.getAtom(i).getAtomCoordinates()))
					            .mapToDouble(t->{
					            	return ((t.k().getX()-tdx1)*(t.v().getX()-tdx2) +(t.k().getY()-tdy1)*(t.v().getY()-tdy2))/ cscale;
					            })
					            .sum() ;
			
//			System.out.println("Scale:" + s);
//			double b1 = c1.bonds().mapToDouble(b->b.getBondLength()).average().orElse(1);
//			double b2 = c2.bonds().mapToDouble(b->b.getBondLength()).average().orElse(1);
//			System.out.println("Scale2:" + b1/b2);
			
		}else{
			//scale c2 to c1:
			double b1 = c1.bonds().mapToDouble(b->b.getBondLength()).average().orElse(1);
			double b2 = c2.bonds().mapToDouble(b->b.getBondLength()).average().orElse(1);
			s = b1/b2;
		}
		
		
		
		
		
		

		
	
		
		AffineTransform at= new AffineTransform();
		at.translate(cpt1.getX(), cpt1.getY());
		at.scale(s, s);
		at.translate(-cpt2.getX(), -cpt2.getY());
		
		
		c2.atoms()
		  .map(a->Tuple.of(a,a.getAtomCoordinates()))
		  .map(Tuple.vmap(ac->asPoint(ac)))
		  .map(Tuple.vmap(p->at.transform(p, null)))
		  .map(Tuple.vmap(p->fromPoint(p)))
		  .forEach(t->{
			  t.k().setAtomCoordinates(t.v());
		  });
		
		System.out.println("COMBINED");
		
		System.out.println(combineChemicals(c1, c2).toMol());
		
		
		long c1TotCount = c1.atoms().count();
		
		double dd=		Math.sqrt(c1.atoms()
				  .map(a->a.getAtomCoordinates())
				  .mapToDouble(p->c2.atoms().mapToDouble(a->a.getAtomCoordinates().distanceSquaredTo(p))
						  					   .min()
						  					   .orElse(0))
				  .sum()/c1TotCount);
		//center, I guess
		//System.out.println("RMSE=" + dd);
		return dd;
		
		
		
	}
	
	private static Point2D asPoint(AtomCoordinates ac){
		return new Point2D.Double(ac.getX(),ac.getY());
	}
	private static AtomCoordinates fromPoint(Point2D pt){
		return AtomCoordinates.valueOf(pt.getX(),pt.getY());
	}
	
	private static Chemical wiggleNoise(Chemical cc) throws IOException{
		Chemical c1 = Chemical.parseMol(cc.toMol());
		double b1avg = c1.bonds().mapToDouble(b->b.getBondLength()).average().orElse(1);

		AffineTransform at2= new AffineTransform();
		at2.scale(1/b1avg, 1/b1avg);
		c1.atoms()
		  .map(a->Tuple.of(a,a.getAtomCoordinates()))
		  .map(Tuple.vmap(ac->asPoint(ac)))
		  .map(Tuple.vmap(p->at2.transform(p, null)))
		  .map(Tuple.vmap(p->fromPoint(p)))
		  .forEach(t->{
			  t.k().setAtomCoordinates(t.v());
		  });
		
		c1.atoms()
		  .map(a->Tuple.of(a,a.getAtomCoordinates()))
		  .map(Tuple.vmap(ac->asPoint(ac)))
		  .map(Tuple.vmap(p->new Point2D.Double(p.getX() + (Math.random()-0.5)*1/(35.0), p.getY() + (Math.random()-0.5)*1/(35.0))))
		  .map(Tuple.vmap(p->fromPoint(p)))
		  .forEach(t->{
			  t.k().setAtomCoordinates(t.v());
		  });
		return c1;
	}
	private static File stdResize(File f, File imageFile, double scale, Interpolation terp, double quality) throws IOException{
		

		


		RenderedImage ri = ImageUtil.decode(f);
		
		int nwidth=(int) (ri.getWidth() *scale);
		int nheight=(int) (ri.getHeight() *scale);
		BufferedImage outputImage=null;
		
		if(Interpolation.SINC.equals(terp)){
			
			ResampleOp resizeOp = new ResampleOp(nwidth, nheight);
			outputImage = resizeOp.filter(convertRenderedImage(ri), null);
			
		}else{
			
	        // creates output image
	        outputImage = new BufferedImage(nwidth,
	                nheight,BufferedImage.TYPE_3BYTE_BGR);
	 
	        // scales the input image to the output image
	        Graphics2D g2d = outputImage.createGraphics();
	        
	        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
	        g2d.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
	        
	        if(Interpolation.BICUBIC.equals(terp)){
	        	g2d.setRenderingHint(RenderingHints.KEY_INTERPOLATION,
	        	       RenderingHints.VALUE_INTERPOLATION_BICUBIC);
	        }else if(Interpolation.BILINEAR.equals(terp)){
	        	g2d.setRenderingHint(RenderingHints.KEY_INTERPOLATION,
		        	       RenderingHints.VALUE_INTERPOLATION_BILINEAR);
	        }else if(Interpolation.NEAREST_NEIGHBOR.equals(terp)){
	        	g2d.setRenderingHint(RenderingHints.KEY_INTERPOLATION,
		        	       RenderingHints.VALUE_INTERPOLATION_NEAREST_NEIGHBOR);
	        }
	        g2d.scale(scale, scale);
	        g2d.drawImage(convertRenderedImage(ri), 0, 0,null);
	        g2d.dispose();

	        for (int x = 0; x < outputImage.getWidth(); x++) {
	            for (int y = 0; y < outputImage.getHeight(); y++) {
	                int rgba = outputImage.getRGB(x, y);
	                Color col = new Color(rgba, false);
	                col = new Color(col.getRed(),
	                		col.getRed(),
	                		col.getRed());
	                outputImage.setRGB(x, y, col.getRGB());
	                
	            }
	        }
		}
  
		if(quality<1 && quality>0){
			JPEGImageWriteParam jpegParams = new JPEGImageWriteParam(null);
			jpegParams.setCompressionMode(ImageWriteParam.MODE_EXPLICIT);
			jpegParams.setCompressionQuality((float) quality);
			
			final ImageWriter writer = ImageIO.getImageWritersByFormatName("jpg").next();
			// specifies where the jpg image has to be written
			
			File tfile = new File(imageFile.getAbsolutePath() + ".jpg");
			
			
			try(FileImageOutputStream fos = new FileImageOutputStream(
					tfile)){
				writer.setOutput(fos);
				writer.write(null, new IIOImage(outputImage, null, null), jpegParams);
				return tfile;
			}
		}
		
		ImageIO.write(outputImage, "png", imageFile);
		
		return imageFile;
	}
	public static BufferedImage convertRenderedImage(RenderedImage img) {
	    if (img instanceof BufferedImage) {
	        return (BufferedImage)img;  
	    }   
	    ColorModel cm = img.getColorModel();
	    int width = img.getWidth();
	    int height = img.getHeight();
	    WritableRaster raster = cm.createCompatibleWritableRaster(width, height);
	    boolean isAlphaPremultiplied = cm.isAlphaPremultiplied();
	    Hashtable properties = new Hashtable();
	    String[] keys = img.getPropertyNames();
	    if (keys!=null) {
	        for (int i = 0; i < keys.length; i++) {
	            properties.put(keys[i], img.getProperty(keys[i]));
	        }
	    }
	    BufferedImage result = new BufferedImage(cm, raster, isAlphaPremultiplied, properties);
	    img.copyData(raster);
	    return result;
	}
	
	
	public static TestResult testMolecule(File image, File sdf){
		return testMolecule(image, sdf, 60, Method.MOLVEC.adapt());
	}
	public static TestResult testMolecule(File image, File sdf, long timeoutInSeconds, MethodAdapted meth){
		long start= System.currentTimeMillis();
		try{
			AtomicBoolean lastWasAlias = new AtomicBoolean(false);
			String rawMol=null;
			
			boolean assmiles = sdf.getAbsolutePath().endsWith("smi");
			
			boolean[] isSgroupMol = new boolean[]{false};
			
			List<Integer> indexesInSGroup = new ArrayList<Integer>();
			
			try(Stream<String> stream = Files.lines(sdf.toPath())){
				String rawMol2 = stream.collect(Collectors.joining("\n"));
				
				
				
				
				if(rawMol2.contains("SUP") || rawMol2.contains("SGROUP") || rawMol2.contains("M  S")){
					isSgroupMol[0]=true;
				}
				rawMol =Arrays.stream(rawMol2.split("\n"))
					.filter(l->!l.contains("SUP"))
					.filter(l->!l.contains("SGROUP"))
					.filter(l->{
						boolean ignore = l.contains("M  S");
						
						if(ignore){
							if(l.startsWith("M  SAL")){
								String  ll = l.substring("M  SAL   1  2".length());
								Arrays.stream(ll.split(" "))
								      .filter(mm->mm.length()>0)
								      .map(mm->Integer.parseInt(mm))
								      .forEach(ii->{
								    	  indexesInSGroup.add(ii);
								      });
							}
						}
						
						
						return !ignore;
					})
					.filter(l->{
						if(lastWasAlias.get()){
							lastWasAlias.set(false);
							return false;
						}
						return true;
					})
					.filter(l->{
						if(l.startsWith("A")){
							lastWasAlias.set(true);
							return false;
						}
						lastWasAlias.set(false);
						return true;
					})
					.collect(Collectors.joining("\n"));
			}
			Chemical c1=null;
			if(assmiles){
				c1=ChemicalBuilder.createFromSmiles(rawMol).build();
			}else{
				c1=ChemicalBuilder.createFromMol(rawMol,Charset.defaultCharset()).build();
			}
			System.out.println("--------------------------------");
			System.out.println(sdf.getAbsolutePath());
			System.out.println("--------------------------------");
			
			for(int ii : indexesInSGroup){
				c1.getAtom(ii-1).setAtomToAtomMap(ii);
			}
			

			
			Chemical c=null;
			
			if(Method.EXACT.equals(meth.method)){
				c= getCleanChemical(c1.toMol());
			}else{
				c=meth.getChem(image);
			}
			
			c.atoms()
			 .filter(a->a.hasAromaticBond())
			 .filter(a->a.getSymbol().equals("C"))
			 .forEach(aa->aa.setImplicitHCount(Math.max(0,3-aa.getBonds().size())));
			c.atoms()
			 .filter(a->a.hasAromaticBond())
			 .filter(a->a.getSymbol().equals("N"))
			 .forEach(aa->aa.setImplicitHCount(0));
			c.kekulize();
			
//			c =Chemical.parseMol(c.toMol()).toBuilder().aromatize(false).build();
			
			long total = System.currentTimeMillis()-start;
			
			if(c.getAtomCount()==0){
				return TestResult.of(Result.FOUND_NOTHING,total);
			}
			
			String tiiinchi=null;
			try{
				tiiinchi=Inchi.asStdInchi(c).getKey();
			}catch(Exception e){
				System.out.println(c.toMol());
				throw e;
			}
			String rinchi=Inchi.asStdInchi(c1).getKey();
			
			String iinchi= tiiinchi;
			
					
			
			int ratomCount=c1.getAtomCount();
			int iatomCount=c.getAtomCount();
			
			int rbondCount=c1.getBondCount();
			int ibondCount=c.getBondCount();
			
			
			String smilesReal=c1.toSmiles();
			String smilesFound=c.toSmiles();
			
			String formReal=c1.getFormula();
			String formFound=c.getFormula();
			
			System.out.println("Real:" + c1.toSmiles());
			System.out.println("Image:" + c.toSmiles());
			
			if(smilesFound.equals("")){
				return TestResult.of(Result.FOUND_NOTHING,total);
			}
			
			if(rinchi.equals(iinchi)){
				//if(isSgroupMol[0]){
				
				if(EXPORT_CORRECT){
					String exp1=exportDir + "/" + meth + "/correct";
					File dir = new File(exp1);
					dir.mkdirs();
					
					String exFile = exp1 + "/" + image.getName() + ".mol";
					
					String mfile = c.toMol();
					try(PrintWriter pw = new PrintWriter(exFile)){
						pw.print(mfile);
					}
					try(FileOutputStream fos = new FileOutputStream(exp1 + "/" + image.getName())){
						Files.copy(image.toPath(), fos);	
					}
				}
				
				return TestResult.of(Result.CORRECT_FULL_INCHI,total,align(c1,c));
				//}else{
				//	return TestResult.of(Result.CORRECT_FULL_INCHI,total);	
				//}
				
			}
			if(rinchi.split("-")[0].equals(iinchi.split("-")[0])){
				//if(isSgroupMol[0]){
					return TestResult.of(Result.CORRECT_STEREO_INSENSITIVE_INCHI,total,align(c1,c));
				//}else{
				//	return TestResult.of(Result.CORRECT_STEREO_INSENSITIVE_INCHI,total);	
				//}
				
			}
			
			if(formReal.equals(formFound)){
				//System.out.println("Matched!");
				return TestResult.of(Result.FORMULA_CORRECT,total);
			}else{
				String withoutHydrogensReal =formReal.replaceAll("H[0-9]*", "");
				String withoutHydrogensFound=formFound.replaceAll("H[0-9]*", "");
				
				if(withoutHydrogensReal.equals(withoutHydrogensFound)){
					return TestResult.of(Result.ATOMS_RIGHT_WRONG_BONDS,total);
				}
				
				if(smilesReal.contains(".") || smilesFound.contains(".") ){
					String largestR=Arrays.stream(smilesReal.split("[.]"))
					      .map(t->Tuple.of(t,t.length()).withVComparator())
					      .max(Comparator.naturalOrder())
					      .map(t->t.k())
					      .orElse(null);
					
					String largestF=Arrays.stream(smilesFound.split("[.]"))
					      .map(t->Tuple.of(t,t.length()).withVComparator())
					      .max(Comparator.naturalOrder())
					      .map(t->t.k())
					      .orElse(null);
					Chemical clargestR=ChemicalBuilder.createFromSmiles(largestR).build();
					Chemical clargestF=ChemicalBuilder.createFromSmiles(largestF).build();

					
					iinchi=Inchi.asStdInchi(clargestF).getKey();
					rinchi=Inchi.asStdInchi(clargestR).getKey();
					
					if(rinchi.equals(iinchi)){
						return TestResult.of(Result.LARGEST_FRAGMENT_CORRECT_FULL_INCHI,total);
					}
					if(rinchi.split("-")[0].equals(iinchi.split("-")[0])){
						return TestResult.of(Result.LARGEST_FRAGMENT_CORRECT_STEREO_INSENSITIVE_INCHI,total);
					}
					
					if(clargestR.getFormula().equals(clargestF.getFormula())){
						return TestResult.of(Result.LARGEST_FRAGMENT_FORMULA_CORRECT,total);
					}
					String fragwHReal =clargestR.getFormula().replaceAll("H[0-9]*", "");
					String fragwHFound=clargestF.getFormula().replaceAll("H[0-9]*", "");
					
					if(fragwHReal.equals(fragwHFound)){
						return TestResult.of(Result.LARGEST_FRAGMENT_ATOM_RIGHT_WRONG_BONDS,total);
					}
					clargestR.makeHydrogensImplicit();
					clargestF.makeHydrogensImplicit();
					
					int fratomCount=clargestR.getAtomCount();
					int fiatomCount=clargestF.getAtomCount();
					
					int frbondCount=clargestR.getBondCount();
					int fibondCount=clargestF.getBondCount();
					if(fratomCount==fiatomCount && frbondCount == fibondCount){
						return TestResult.of(Result.LARGEST_FRAGMENT_ATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY,total);
					}
				}
				if(smilesReal.contains("|")){
					return TestResult.of(Result.WEIRD_SOURCE,total);
				}
				if(ratomCount==iatomCount && rbondCount == ibondCount){
					return TestResult.of(Result.ATOM_COUNT_BOND_COUNT_RIGHT_WRONG_LABELS_OR_CONNECTIVITY,total);
				}
				
				//System.out.println("NO MATCH!");
				return TestResult.of(Result.INCORRECT,total);
			}
		}catch(Exception e){
			e.printStackTrace();
		}
		return TestResult.of(Result.ERROR,System.currentTimeMillis()-start);
	}
	
	
	public void doAllIdentityDataSetTests(MethodAdapted adapted) throws FileNotFoundException{
		testSet("uspto", adapted);
		testSet("usan", adapted);
		testSet("testSet1", adapted);
		testSet("maybridge", adapted);
		testSet("trec", adapted);		
	}
	
	
	public void doAllCompressionQualityDataSetTests(MethodAdapted adapted) throws FileNotFoundException{
		for(int i=1;i<=10;i++){
			double q=i/10.0;
			adapted= adapted.suffix("_jpg[" + q + "]")
						    .quality(q);
			testSet("uspto", adapted);
			testSet("usan", adapted);
			testSet("testSet1", adapted);
			testSet("maybridge", adapted);
			testSet("trec", adapted);
		}		
	}
	
	public void doAllBicubicScaleQualityDataSetTests(MethodAdapted adapted) throws FileNotFoundException{
		for(int i=3;i<=20;i++){
			adapted= adapted.suffix("bicubic")
							.interpolation(Interpolation.BICUBIC)
							.scale(i/10.0);
			testSet("uspto", adapted);
			testSet("usan", adapted);
			testSet("testSet1", adapted);
			testSet("maybridge", adapted);
			testSet("trec", adapted);
		}
	}
	
	public void doAllBilinearScaleQualityDataSetTests(MethodAdapted adapted) throws FileNotFoundException{
		for(int i=3;i<=20;i++){
			adapted= adapted.suffix("bilinear")
							.interpolation(Interpolation.BILINEAR)
							.scale(i/10.0);
			testSet("uspto", adapted);
			testSet("usan", adapted);
			testSet("testSet1", adapted);
			testSet("maybridge", adapted);
			testSet("trec", adapted);
		}
	}
	
	public void doAllSincScaleQualityDataSetTests(MethodAdapted adapted) throws FileNotFoundException{
		for(int i=3;i<=20;i++){
			adapted= adapted.suffix("sinc")
							.interpolation(Interpolation.SINC)
							.scale(i/10.0);
			testSet("uspto", adapted);
			testSet("usan", adapted);
			testSet("testSet1", adapted);
			testSet("maybridge", adapted);
			testSet("trec", adapted);
		}
	}
	
	
	public void doAllRMSEDataSetTests(MethodAdapted adapted) throws FileNotFoundException{
			adapted= adapted.rmse(true);
			testSet("uspto", adapted);
			testSet("trec", adapted);		
	}
	
	
	
	
	@Test
	public void test1() throws FileNotFoundException{
		
//		
		
		
		
		//*********************		
		//      Full quality analysis molvec
		//*********************		
//		for(int i=1;i<=10;i++){
//			double q=i/10.0;
//			testSet("trec", Method.MOLVEC.adapt().suffix("_jpg[" + q+ "]").quality(q));
//		}		
		//*********************

//		
//		for(int i=1;i<=10;i++){
//			double q=i/10.0;
//			testSet("trec", Method.IMAGO.adapt().suffix("_jpg[" + q+ "]").quality(i/10.0));
//		}
//		for(int i=1;i<=10;i++){
//			double q=i/10.0;
//			testSet("trec", Method.OSRA.adapt().suffix("_jpg[" + q+ "]").quality(i/10.0));
//		}
		
		
		//Full rescale for Molvec
//		for(int i=3;i<=20;i++){
//			testSet("trec", Method.MOLVEC.adapt().suffix("bicubicLim").interpolation(Interpolation.BICUBIC).scale(i/10.0));
//		}
//		for(int i=3;i<=20;i++){
//			testSet("trec", Method.MOLVEC.adapt().suffix("sincLim").interpolation(Interpolation.SINC).scale(i/10.0));
//		}
//		for(int i=3;i<=20;i++){
//			testSet("trec", Method.MOLVEC.adapt().suffix("bilinearLim").interpolation(Interpolation.BILINEAR).scale(i/10.0));
//		}
		
		
		
//		for(int i=3;i<=20;i++){
//			testSet("uspto", Method.MOLVEC.adapt().suffix("bicubicLim").limit(100).interpolation(Interpolation.BICUBIC).scale(i/10.0));
//		}
//		for(int i=3;i<=20;i++){
//			testSet("uspto", Method.MOLVEC.adapt().suffix("sincLim").limit(100).interpolation(Interpolation.SINC).scale(i/10.0));
//		}
//		for(int i=3;i<=20;i++){
//			testSet("uspto", Method.MOLVEC.adapt().suffix("bilinearLim").limit(100).interpolation(Interpolation.BILINEAR).scale(i/10.0));
//		}
		
		
		
		
		testSet("uspto", Method.MOLVEC.adapt());
		testSet("usan", Method.MOLVEC.adapt());
		
//		testSet("maybridge", Method.MOLVEC.adapt());
//		testSet("testSet1", Method.MOLVEC.adapt());
		
//		for(int i=20;i<=20;i++){
//			testSet("trec", Method.MOLVEC.adapt().suffix("bilinearLim2").interpolation(Interpolation.BILINEAR).scale(i/10.0));
//		}
		
		
		
//		
		
		
//		
//		for(int i=17;i<=20;i++){
//			testSet("trec", Method.OSRA.adapt().suffix("_sinc").interpolation(Interpolation.SINC).scale(i/10.0));
//		}
		
//		for(int i=3;i<=20;i++){
//			testSet("trec", Method.IMAGO.adapt().suffix("_sinc").scale(i/10.0));
//		}
//		for(int i=3;i<=20;i++){
//			testSet("trec", Method.OSRA.adapt().suffix("_sinc").scale(i/10.0));
//		}
//		
//		testSet("uspto", Method.MOLVEC.adapt().rmse(true));
		//All data sets stats
		
//		for(int j=0;j<20;j++){
//			double sig = (j-10)/5.0;
//			StructureImageExtractor.THRESH_STDEV=sig;
//			for(int i=10;i<=10;i++){
//				testSet("trec", Method.MOLVEC.adapt().limit(10).suffix("_sig[" + sig + "]").scale(i/10.0));
//			}
//		}
		
//		testSet("usan", Method.OSRA.adapt().scale(0.5));
//		testSet("maybridge", Method.MOLVEC.adapt());
//		testSet("trec", Method.MOLVEC.adapt());
//		testSet("uspto", Method.MOLVEC.adapt());
//		testSet("testSet1", Method.MOLVEC.adapt());
//		testSet("usan", Method.MOLVEC.adapt());
//		testSet("usan", Method.MOLVEC.adapt().scale(0.5));
		
		
//		
////		testSet("uspto", Method.MOLVEC.adapt().rmse(true));
////		
////		
//		for(int i=3;i<=20;i++){
//			testSet("trec", Method.IMAGO.adapt().scale(i/10.0));
//		}
//		for(int i=3;i<=20;i++){
//			testSet("trec", Method.MOLVEC.adapt().scale(i/10.0));
//		}
//		for(int i=4;i<=20;i++){
//			testSet("trec", Method.OSRA.adapt().scale(i/10.0));
//		}
//		
		
//		for(int j=1;j<20;j++){
//			double sig = j*5.0;
//			StructureImageExtractor.DEF_BINARIZATION = new RangeFractionThreshold(sig/100.0);
//			StructureImageExtractor.RESIZE_BINARIZATION = new RangeFractionThreshold(sig/100.0);
//			
//			for(int i=3;i<=10;i++){
//				testSet("trec", Method.MOLVEC.adapt().limit(100).suffix("_sig[" + sig + "]").scale(i/10.0));
//			}
//		}
		
//		testSet("usan", Method.MOLVEC.adapt());
//		testSet("usan", Method.MOLVEC.adapt().scale(0.5));
//		
		
		//testSet("uspto",true);

		//IMAGO_SCALE
		
		//RegressionTestIT.EXPORT_CORRECT=true;

		
//
//		for(int i=15;i<=20;i++){
//			testSet("trec", Method.IMAGO.adapt().scale(i/10.0));
//		}
//		for(int i=14;i<=20;i++){
//			testSet("trec", Method.OSRA.adapt().scale(i/10.0));
//		}
		//
//			
		
	}
	
	private static long DEF_TIMEOUT = 400;
	
	
	
	public static enum Method{
		OSRA(f->()->getOSRAChemical(f).toMol()),
		MOLVEC(f->()->{
			CompletableFuture<String> chemicalCompletableFuture = Molvec.ocrAsync(f);
			//
			try {
				return chemicalCompletableFuture.get(DEF_TIMEOUT, TimeUnit.SECONDS);
			}catch(TimeoutException te) {
				System.out.println("timeout!!");
				chemicalCompletableFuture.cancel(true);
				throw new RuntimeException(te);
			}
		}),
		IMAGO(f->()->getImagoChemical(f).toMol()),
		EXACT(f->null);
		
		Function<File,Callable<String>> fetcher;
		
		
		Method(Function<File,Callable<String>> mfileFetcher ){
			this.fetcher=mfileFetcher;
		}
		
		
		public Chemical getChem(File f) throws Exception{
			return getCleanChemical(fetcher.apply(f).call());
		}
		
		public MethodAdapted adapt(){
			return MethodAdapted.from(this);
		}
		
	}
	
	public static enum Interpolation{
		BICUBIC,
		BILINEAR,
		SINC,
		NEAREST_NEIGHBOR;
	}
	
	
	public static class MethodAdapted{
		public Method method;
		Interpolation terp = Interpolation.BICUBIC;
		double scale = 1;
		double quality = 1;
		long max = Long.MAX_VALUE;
		boolean RMSE = false;
		String suffix =""; 
		
		
		
		public MethodAdapted (Method m){
			this.method=m;
		}
		
		public MethodAdapted quality(double q){
			this.quality=q;
			return this;
		}
		public MethodAdapted scale(double s){
			this.scale=s;
			return this;
		}
		
		public MethodAdapted interpolation(Interpolation terp){
			this.terp=terp;
			return this;
		}
		
		public MethodAdapted suffix(String suf){
			this.suffix=suf;
			return this;
		}
		
		public MethodAdapted limit(long l){
			this.max=l;
			return this;
		}
		public MethodAdapted rmse(boolean rmse){
			this.RMSE=rmse;
			return this;
		}
		
		public boolean rmse(){
			return this.RMSE;
		}
		
		public long getLimit(){
			return this.max;
		}
		
		public Chemical getChem(File f) throws Exception{
			File imageFile = File.createTempFile("tmpStr" + UUID.randomUUID().toString() +"_scale_" +  scale + "x" , ".png");
			imageFile=stdResize(f, imageFile, scale, terp, this.quality);
			return method.getChem(imageFile);
		}
		
		public static MethodAdapted from(Method m){
			return new MethodAdapted(m);
		}
		
		public String toString(){
			String limitAdd = "";
			if(this.max<Long.MAX_VALUE){
				limitAdd = "_limit_"+this.max;
			}
			return this.method.toString() + "_" + scale + "x" + ((RMSE)?"_RMSE":"" + limitAdd + suffix);
		}
		
	}
	
	
	
//	@Ignore
	
	public void testSet(String set, MethodAdapted meth) throws FileNotFoundException{
		
		RegressionTestIT.DO_ALIGN=meth.rmse();
		
		
		
		try(PrintWriter pw1 = new PrintWriter("/home/tyler/workspace/molvec/reports/" + set + "_" +  meth.toString() + "_" +
				LocalDate.now().format( DateTimeFormatter.ofPattern("YYYYMMdd")) + 				
				".txt")){
			
		
		long start = System.currentTimeMillis();
		
		
		File dir1 = getFile("regressionTest/" + set);
		
		List<String> dataMethod = new ArrayList<>();
		dataMethod.add("@Parameterized.Parameters(name=\"{0}\")");
		dataMethod.add("public static List<Object[]> getData(){");
		dataMethod.add("\tFile dir = new File(RegressionTest2.class.getResource(\"/regressionTest/testSet1\").getFile());");

		dataMethod.add("\n\tList<Object[]> list = new ArrayList<>();\n");

		Map<Result, List<Tuple<Tuple<TestResult, List<File>>, List<File>>>> mm=
			Arrays.stream(dir1.listFiles())
		      .filter(f->f.getName().contains("."))
		      //.filter(f->f.getName().contains("cas-382-67-2"))
		      .map(f->Tuple.of(f.getName().split("[.]")[0],f))
		      
		      .collect(Tuple.toGroupedMap())
		      .values()
		      .stream()
		      .filter(l->l.size()==2)
		      
		      .map(l->{
		    	  	if(!l.get(1).getName().toLowerCase().endsWith("tif") && !l.get(1).getName().toLowerCase().endsWith("png")){
		    	  		List<File> flist = new ArrayList<File>();
		    	  		flist.add(l.get(1));
		    	  		flist.add(l.get(0));
		    	  		return flist;
					}
		    	  	return l;
		      })
		      
		      .collect(shuffler(new Random(12440l)))		      
		      .limit(meth.getLimit())

//NOTE, I THINK THIS TECHNICALLY WORKS, BUT SINCE THERE IS PARALLEL THINGS GOING ON IN EACH, IT SOMETIMES WILL STARVE A CASE FOR A LONG TIME
		      .parallel()
		      
		      
		      .map(fl->Tuple.of(fl,testMolecule(fl.get(1),fl.get(0), 400, meth)))
		      .map(t->t.swap())
		      .map(t->Tuple.of(t.k().result,Tuple.of(t,t.v())))
		      .peek(t->System.out.println(t.v().v().get(1).getAbsolutePath() + ":" +t.k()))
		      .collect(Tuple.toGroupedMap());
			
		
		  pw1.println("~~~~~~~~~~~~~~~~~~~~~~~~~~");
		  pw1.println("Time\t" + ((System.currentTimeMillis()-start)/1000.0) + " seconds");
		  pw1.println("Scale\t" + meth.scale);
		  pw1.println("Method\t" + meth.method);
		  pw1.println("!!!!!!!!!!!!!!!!!!!!!!!!!!");
		
		  for(Result r:Result.values()){
			  int c=Optional.ofNullable(mm.get(r))
			          .map(t->t.size())
			          .orElse(0);
			  pw1.println(r.toString() + "\t" + c);
		  }
		
			  mm
			  .entrySet()
		      .stream()
		      .map(Tuple::of)
		      .sorted(Comparator.comparing(t->t.k().ordinal()))
		      .forEach(t->{
		    	  Result r=t.k();
		    	  List<List<File>> fl = t.v().stream().map(t1->t1.v()).collect(Collectors.toList());
		    	  
		    	  pw1.println("======================================");
		    	  pw1.println(r.toString() + "\t" + fl.size());
		    	  pw1.println("--------------------------------------");
		    	  pw1.println(t.v().stream().map(tf->tf.v().get(1).getAbsolutePath() + "\t" + tf.k().k().time + "\t" + tf.k().k().RMSE).collect(Collectors.joining("\n")));

		      });
			}
	}
	
	public static <T> Collector<T,List<T>,Stream<T>> shuffler(Random r){
		
		return new Collector<T,List<T>,Stream<T>>(){

			@Override
			public BiConsumer<List<T>, T> accumulator() {
				return (l,t)->{
					l.add(t);
				};
			}

			@Override
			public Set<java.util.stream.Collector.Characteristics> characteristics() {
				//java.util.stream.Collector.Characteristics.
				return new HashSet<java.util.stream.Collector.Characteristics>();
			}

			@Override
			public BinaryOperator<List<T>> combiner() {
				return (l1,l2)->{
					l1.addAll(l2);
					return l1;
				};
			}

			@Override
			public Function<List<T>, Stream<T>> finisher() {
				return (u)->{
					Collections.shuffle(u,r);
					return u.stream();
				};
			}

			@Override
			public Supplier<List<T>> supplier() {
				return ()-> new ArrayList<T>();
			}

		};
	}
	
	
}
