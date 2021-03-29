package gov.nih.ncats.molvec;

import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.text.NumberFormat;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import gov.nih.ncats.common.stream.StreamUtil;
import gov.nih.ncats.molvec.internal.algo.StructureImageExtractor;
import gov.nih.ncats.molvec.internal.algo.experimental.ConstantValueResultScorer;
import gov.nih.ncats.molvec.internal.algo.experimental.ResultScorer;
import gov.nih.ncats.molvec.internal.util.CachedSupplier;
import gov.nih.ncats.molvec.internal.util.ConnectionTable;
import gov.nih.ncats.molvec.internal.util.GeomUtil;

public class MolvecOptions {
	private static int DEFAULT_LIMIT_TRIES=7;
	private static int DEFAULT_MAX_LIMIT_TRIES=28;
	
	private static int[][] ADDITIONAL_PAIRS = new int[][] {
		new int[]{74,30,75},
		new int[]{74,59},
		new int[]{74,30,1,75},
	};
	
	private static int[] FLAG_TRY_ORDER=new int[]{
			-1,
			74,
			30,
			35,
			1,
			75,
			59,
			5,
			7,
			6,
			76,
			34, //
			13,
			36,
			29, //
			25,
			12, //
			19, //
			43, //
            77,
            84, //Testing
            85, //Testing 2
            86, //Testing 3
//            87,
			58,
			2,  //
			26, //
			28, //
			21,
			4,
			8,
			14,
			20,
			3,
			31,
			44,
			54,
			17,
			23,
			41,
			47,
			50,
			63,
			70,
			71,
			0,
			18,
			46,
			64,
			72,
			9,
			};
	static{
		FLAG_TRY_ORDER=StreamUtil.with(Arrays.stream(FLAG_TRY_ORDER).boxed())
				.and(IntStream.range(0,100).boxed())
				.stream()
				.distinct()
				.mapToInt(i->i)
				.toArray();
		
	}
    private double averageBondLength = 1D;
    private boolean center = true;
    private boolean includeSgroups = true;

    private String name;
    
    private BitSet flags = new BitSet();
    private int[] flagTries = new int[]{-1};
    private ResultScorer scorer = new ConstantValueResultScorer(1.0);
    
    
    private boolean overrideScorer=true;
  
    
    
    public boolean isOverrideScorer() {
		return overrideScorer;
	}



	public MolvecOptions setOverrideScorer(boolean overrideScorer) {
		this.overrideScorer = overrideScorer;
		return this;
	}



	public ResultScorer getScorer() {
		return scorer;
	}



	public MolvecOptions setScorer(ResultScorer scorer) {
		this.scorer = scorer;
		return this;
	}



	public int[] getFlagTries() {
		return flagTries;
	}



	public void setFlagTries(int[] flagTries) {
		this.flagTries = flagTries;
	}

	private StructureImageExtractor.ImageExtractionValues values = StructureImageExtractor.DEFAULT_VALUES;
    
    public StructureImageExtractor.ImageExtractionValues getValues() {
		return values;
	}



	public void setValues(StructureImageExtractor.ImageExtractionValues values) {
		this.values = values;
	}



	public BitSet getFlags() {
		return flags;
	}



	public MolvecOptions setFlags(BitSet bs){
    	flags=bs;
    	return this;
    }
	
	public MolvecOptions setDebug(){
		values= this.values.debug(true);
		return this;
	}
	
	public MolvecOptions modFlags(){
		BitSet bs = new BitSet();
		bs.set(0, 100);
		flagTries = Arrays.copyOf(FLAG_TRY_ORDER, FLAG_TRY_ORDER.length);
		
		return this.setFlags(bs).limitAttempts(DEFAULT_LIMIT_TRIES);
	}
	
	public MolvecOptions limitAttempts(int max){
		this.flagTries=Arrays.stream(Arrays.copyOf(FLAG_TRY_ORDER, FLAG_TRY_ORDER.length))
				.limit(max)
				.toArray();
		return this;
	}
	
	public MolvecOptions clearFlag(int p){
		this.flags.clear(p);
		return this;
	}
    
    
    
    public MolvecOptions setName(String name){
        this.name = name;
        return this;
    }

    public MolvecOptions averageBondLength(double averageBondLength){
        if(averageBondLength <=0){
            throw new IllegalArgumentException("avg bond length must be > 0");
        }
        this.averageBondLength = averageBondLength;
        return this;
    }

    public MolvecOptions center(boolean center){
        this.center = center;
        return this;
    }
    public MolvecOptions includeSgroups(boolean includeSgroups){
        this.includeSgroups = includeSgroups;
        return this;
    }

    public MolvecResult computeResult(ConnectionTable ct){
        String mol= toMol(ct);

        return new Result(mol, name, CachedSupplier.of(()-> ct.getNodes()
                                                    .stream()
                                                    .map(n->n.getPoint())
                                                    .collect(GeomUtil.convexHull())
                                                    .getBounds2D()));
    }

    private String toMol(ConnectionTable ct){
        AffineTransform at = new AffineTransform();
        double blcur = Math.max(ct.getAverageBondLength(),1);

        double scale = averageBondLength/blcur;



        at.scale(scale, scale);
        if(center){
            if(ct.getNodes().size()>0){
                Rectangle2D rect=ct.getNodes().stream().map(n->n.getPoint()).collect(GeomUtil.convexHull()).getBounds2D();
                at.translate(-rect.getCenterX(), -rect.getCenterY());
            }
        }



        String newLine = System.lineSeparator();
        StringBuilder headerBuilder = new StringBuilder(80);
        if(name !=null){
            headerBuilder.append(name);
        }
        String header = headerBuilder
                .append(newLine)
//		IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR
//				(FORTRAN: A2<--A8--><---A10-->A2I2<--F10.5-><---F12.5--><-I6-> )
                .append("  Molvec01")
                //date/time (M/D/Y,H:m)
                .append(MOL_DATETIME_FORMATTER.format(LocalDateTime.now()))
                .append("2D") //always write out 2D coords
                .append(newLine).append(newLine)
                .toString();

        String countsLine = new StringBuilder(80)
                .append(writeMolInt(ct.getNodes().size(), 3))
                .append(writeMolInt(ct.getEdges().size(), 3))
                //TODO for now mark eveything as chiral
                //CDK sets this to 1 only if there's a tetrahedral stereo in molecule
                .append("  0  0  0  0  0  0  0  0999 V2000")
                .append(newLine)
                .toString();

        StringBuilder atomBlockBuilder = makeAtomBlock(ct, at, newLine);

        StringBuilder bondBuilder = makeBondBlock(ct,newLine);

        StringBuilder sgroupBuilder =new StringBuilder(0);
        if(includeSgroups){
            sgroupBuilder = makeSGroupBlock(ct, at, newLine, true);
        }



        String mol= new StringBuilder(header.length() + countsLine.length() + atomBlockBuilder.length() + bondBuilder.length() +sgroupBuilder.length())
                .append(header)
                .append(countsLine)
                .append(atomBlockBuilder)
                .append(bondBuilder)
                .append(sgroupBuilder)
                .append("M  END")
                .toString();
        return mol;
    }

    private StringBuilder makeBondBlock(ConnectionTable ct, String newLine){
        StringBuilder bondBuilder = new StringBuilder(ct.getEdges().size() *6);
        boolean useSingle=true;

        for(ConnectionTable.Edge e : ct.getEdges()){
            int order = e.getOrder();
            if(e.isAromatic()){
                order =4;
//				order =useSingle? 1: 2;
//				useSingle = !useSingle;
            }else if(order <1 || order > 4){
                order =1;
            }
            int bondStereo=0;
            if(e.getOrder() ==1){
                if(e.getWedge()){
                    bondStereo =1; // up
                }else if(e.getDashed()){
                    bondStereo=6; //down
                }
            }

            bondBuilder.append(writeMolInt(e.getNode1Offset()+1, 3))
                    .append(writeMolInt(e.getNode2Offset()+1, 3))
                    .append(writeMolInt(order, 3))
                    .append(writeMolInt(bondStereo, 3))
                    .append(newLine);

        }
        return bondBuilder;
    }
    private StringBuilder makeAtomBlock(ConnectionTable ct, AffineTransform at, String newLine){
        StringBuilder atomBlockBuilder = new StringBuilder(82* ct.getNodes().size());
        for(ConnectionTable.Node n : ct.getNodes()){
            Point2D np = at.transform(n.getPoint(), null);

            String sym = n.getSymbol();
            int massDifference = 0;
            if(n.getSymbol().equals("D")){
                sym="H";
                massDifference=1;
            }


            atomBlockBuilder.append(writeMolDouble(np.getX(), 10))
                    .append(writeMolDouble(-np.getY(), 10))
                    .append(writeMolDouble(Double.NaN, 10))	//only write 2d coords
                    .append(' ')
                    .append(leftPaddWithSpaces(sym, 3))		//TODO should we check to make sure sym is always < 3 chars?
                    .append(writeMolInt(massDifference, 2))
                    .append(writeMolInt(computeMolCharge(n.getCharge()), 3))
                    .append("  0  0  0  0  0  0  0  0  0  0")		//TODO really compute charge
                    .append(newLine);
        }
        return atomBlockBuilder;
    }

    private StringBuilder makeSGroupBlock(ConnectionTable ct, AffineTransform at, String newLine, boolean onlyIfTooClose){
        StringBuilder sgroupBlockBuilder = new StringBuilder();
        Map<Integer, List<ConnectionTable.Node>> groups = new HashMap<>();

        Map<ConnectionTable.Edge,Integer> bindex = new HashMap<ConnectionTable.Edge,Integer>();
        for(int i=0;i<ct.getEdges().size();i++){
            bindex.put(ct.getEdges().get(i), i+1);
        }
        Set<Integer> dontDo = new HashSet<>();

        for(ConnectionTable.Node n : ct.getNodes()){
            if(n.getGroup()!=0){
                groups.computeIfAbsent(n.getGroup(), k->new ArrayList<>()).add(n);
                if(onlyIfTooClose){
                    if(!n.isTooClose())dontDo.add(n.getGroup());
                }
            }
        }


        groups.forEach((g,al)->{
            if(dontDo.contains(g))return;
            ConnectionTable.Node aNode = al.stream()
                    .filter(n->n.getAlias()!=null)
                    .findFirst()
                    .orElse(null);
            if(aNode==null)return;

            List<ConnectionTable.Edge> me = al.stream()
                    .flatMap(n->n.getEdges().stream().filter(e->!e.isInventedBond()))
                    .distinct()
                    .collect(Collectors.toList());
            String l1 = "M  STY  1" +rightPaddWithSpaces(g+"", 4) + " SUP";
            String l2 = "M  SLB  1" +rightPaddWithSpaces(g+"", 4) + rightPaddWithSpaces(g+"", 4);
            String la = al.stream().map(n->n.getIndex()+1)
                    .map(i->rightPaddWithSpaces(i+"", 4))
                    .collect(Collectors.joining());
            String bl = me.stream().map(e->bindex.get(e))
                    .map(i->rightPaddWithSpaces(i+"", 4))
                    .collect(Collectors.joining());

            String l3 = "M  SAL" + rightPaddWithSpaces(g+"", 4) +rightPaddWithSpaces(al.size()+"", 3) + la;
            String l4 = "M  SBL" + rightPaddWithSpaces(g+"", 4) +rightPaddWithSpaces(me.size()+"", 3) + bl;
            String l5 = "M  SMT" + rightPaddWithSpaces(g+"", 4) + " " +aNode.getAlias();
            //String l6 = "M  SBV" + leftPaddWithSpaces(g+"", 4) + " " +aNode.getAlias();
            sgroupBlockBuilder.append(l1).append(newLine);
            sgroupBlockBuilder.append(l2).append(newLine);
            sgroupBlockBuilder.append(l3).append(newLine);
            sgroupBlockBuilder.append(l4).append(newLine);
            sgroupBlockBuilder.append(l5).append(newLine);

        });

        return sgroupBlockBuilder;
    }

    private static int computeMolCharge(int charge){
        switch(charge){
            case -3: return 7;
            case -2: return 6;
            case -1: return 5;
            case 0: return 0;
            case 1: return 3;
            case 2: return 2;
            case 3: return 1;
            default: return 0;
        }
    }
    private static String writeMolInt(int value, int numDigits){
        String s = Integer.toString(value);
        if(s.length()>numDigits){
            s="0";
        }
        return rightPaddWithSpaces(s, numDigits);
    }

    private static String writeMolDouble(double d, int width) {
        String value;
        if (Double.isNaN(d) || Double.isInfinite(d)){
            value = "0.0000";
        }else{
            value = MOL_FLOAT_FORMAT.get().format(d);
        }
        return rightPaddWithSpaces(value, width);
    }

    private static String rightPaddWithSpaces(String value, int numDigits) {
        int padd = numDigits - value.length();
        StringBuilder builder = new StringBuilder(numDigits);
        for(int i=0; i< padd; i++){
            builder.append(' ');
        }
        builder.append(value);
        return builder.toString();
    }
    private static String leftPaddWithSpaces(String value, int numDigits) {
        int padd = numDigits - value.length();
        StringBuilder builder = new StringBuilder(numDigits);
        builder.append(value);
        for(int i=0; i< padd; i++){
            builder.append(' ');
        }

        return builder.toString();
    }

    private static DateTimeFormatter MOL_DATETIME_FORMATTER = DateTimeFormatter.ofPattern("MMddyyHHmm");

    private static StringBuilder EMPTY_STRING_BUILDER = new StringBuilder();
    private static ThreadLocal<NumberFormat> MOL_FLOAT_FORMAT = ThreadLocal.withInitial(()->{
        NumberFormat nf = NumberFormat.getNumberInstance(Locale.ENGLISH);
        nf.setMinimumIntegerDigits(1);
        nf.setMaximumIntegerDigits(4);
        nf.setMinimumFractionDigits(4);
        nf.setMaximumFractionDigits(4);
        nf.setGroupingUsed(false);

        return nf;
    });

    private static class Result implements MolvecResult{
        private final String mol;
        private final CachedSupplier<Rectangle2D> boundsSupplier;
        private String name;

        private static final String lineSep = System.lineSeparator();
        public Result(String mol, String name, CachedSupplier<Rectangle2D> boundsSupplier) {
            this.mol = mol;
            this.name = name;
            this.boundsSupplier = boundsSupplier;
        }

        @Override
        public Optional<String> getMolfile() {
            return Optional.of(mol);
        }
        
        @Override
        public Optional<Map<String, String>> getProperties(){
        	if(name!=null){
        		return Optional.of(Collections.singletonMap("Molecule Name", name));
        	}else{
        		return Optional.empty();
        	}
        }

       

        @Override
        public Optional<String> getSDfile(Map<String, String> properties) {
            StringBuilder propertiesBuilder = formatProperties(properties);
            StringBuilder sdBuilder = new StringBuilder(mol.length()+propertiesBuilder.length() + 5);
            sdBuilder.append(mol).append(lineSep)
                    .append(propertiesBuilder)
                    .append("$$$$");

            return Optional.of(sdBuilder.toString());
        }

        private StringBuilder formatProperties(Map<String, String> properties){
            if(properties==null || properties.isEmpty()){
                return EMPTY_STRING_BUILDER;
            }

            StringBuilder builder = new StringBuilder(2000);
            for(Map.Entry<String, String> entry: properties.entrySet()){
                builder.append(">  <").append(entry.getKey()).append('>').append(lineSep)
                        .append(entry.getValue()).append(lineSep).append(lineSep);
            }
            return builder;
        }

        @Override
        public Optional<Rectangle2D> getOriginalBoundingBox() {
            return Optional.of(boundsSupplier.get());
        }

        @Override
        public boolean hasError() {
            return false;
        }

        @Override
        public Optional<Throwable> getError() {
            return Optional.empty();
        }


    }
}
