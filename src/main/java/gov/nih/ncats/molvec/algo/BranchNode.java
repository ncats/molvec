package gov.nih.ncats.molvec.algo;

import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.Supplier;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import gov.nih.ncats.molvec.util.CachedSupplier;


class BranchNode{
	private final static HashSet<String> accept = CachedSupplier.of(()->{
		HashSet<String> keepers=new HashSet<String>();
		keepers.add("C");
		keepers.add("N");
		keepers.add("O");
		keepers.add("H");
		keepers.add("S");
		keepers.add("P");
		keepers.add("Pt");
		keepers.add("K");
		keepers.add("As");
		keepers.add("Au");
		keepers.add("Al");
		keepers.add("D");
		keepers.add("Hg");
		//keepers.add("I");
		keepers.add("Na");
		keepers.add("B"); //(not sure I want to confirm that yet, it's so rare)
		keepers.add("Br");
		keepers.add("Cl");
		keepers.add("F");
		keepers.add("Si");
		
//		keepers.add("Ca");
		
	    return keepers;
	}).get();
	
	
	
	private static Map<String,Optional<BranchNode>> _cache = new ConcurrentHashMap<>();
	
	//TODO
	//Need clone mechanism, but there is a lot of embedded state / links to other things
	
	
	private String symbol = "C";
	
	
	private boolean pseudo=false;
	private boolean isRepeat=false;
	private int rep=0;
	
	private int charge=0;
	private int orderToParent=1;
	
	
	private int thetaOffset=0;
	private boolean linkable =true;
	private boolean isCombiner=false;
	
	//1 is wedge, -1 is dash
	private int wedgeToParent=0;
	private boolean isTerminal =false;
	private boolean combineLinearly = false;
	
	private Point2D suggestedPoint = new Point2D.Double(0, 0);
	
	
	//These are the difficult parts for cloning
	private List<BranchNode> children = new ArrayList<BranchNode>();
	private Tuple<BranchNode,Integer> ringBond = null; 
	private BranchNode childForCombine=null;
	private BranchNode rightBranch=this;
	private BranchNode leftBranch=this;
	
	private String alias=null;
	
	public String getAlias(){
		return this.alias;
	}
	public BranchNode setAlias(String alias){
		this.alias=alias;
		return this;
	}
	
	public Tuple<BranchNode,Boolean> cloneNode(Map<BranchNode,BranchNode> oldToNew){
		BranchNode bn = new BranchNode(this.symbol);
		bn.pseudo=this.pseudo;
		bn.isRepeat=this.isRepeat;
		bn.rep=this.rep;
		bn.charge=this.charge;
		bn.orderToParent=this.orderToParent;
		bn.thetaOffset=this.thetaOffset;
		bn.linkable=this.linkable;
		bn.isCombiner=this.isCombiner;
		bn.wedgeToParent=this.wedgeToParent;
		bn.alias=this.alias;
		bn.isTerminal=this.isTerminal;
		bn.combineLinearly=this.combineLinearly;
		bn.setSuggestedPoint(this.getSuggestedPoint());
		
		
		
		oldToNew.put(this, bn);
		
		boolean complete=true;
		
		for(BranchNode c: children){
			Tuple<BranchNode,Boolean> nc = c.cloneNode(oldToNew);
			if(!nc.v())complete=false;
			oldToNew.put(c, nc.k());
			bn.children.add(nc.k());
		}
		bn.rightBranch = oldToNew.get(this.rightBranch);
		bn.leftBranch = oldToNew.get(this.leftBranch);
		bn.childForCombine =  oldToNew.get(this.childForCombine);
		
		if(ringBond!=null){
			BranchNode newRing = oldToNew.get(this.ringBond.k());
			if(newRing==null){
				complete=false;
				//This part might not work
				newRing=this.ringBond.k();
			}
			bn.ringBond=Tuple.of(newRing,this.ringBond.v());
		}		
		return Tuple.of(bn,complete);
	}
	
	public BranchNode cloneNode(){
		Map<BranchNode,BranchNode> oldToNewMap = new HashMap<>();
		Tuple<BranchNode,Boolean> clone = cloneNode(oldToNewMap);
		
		if(clone.v()){
			return clone.k();
		}else{
			clone.k().forEachBranchNode((p,c)->{
				if(c.ringBond!=null){
					c.ringBond = Tuple.of(oldToNewMap.get(c.ringBond.k()),c.ringBond.v());
				}
			});
		}
		return clone.k();
	}
	
	
	
	public int getOrderToParent(){
		return this.orderToParent;
	}
	
	
	public BranchNode setLinkable(boolean b){
		this.linkable=b;
		return this;
	}
	public BranchNode addRing(BranchNode to, int order){
		this.ringBond=Tuple.of(to,order);
		return this;
	}
	public Optional<Tuple<BranchNode,Integer>> getRingBond(){
		return Optional.ofNullable(ringBond);
	}
	
	public BranchNode thetaOffset(int o){
		this.thetaOffset=o;
		return this;
	}
	
	public BranchNode setOrderToParent(int o){
		this.orderToParent=o;
		return this;
	}
	
	public BranchNode setCharge(int c){
		this.charge=c;
		return this;
	}
	public int getCharge(){
		return this.charge;
	}
	
	public BranchNode generateCoordinates(){
		setSuggestedPoint(new Point2D.Double(0, 0));
		int totalChildren = children.size();
		
		double[] thetas=new double[]{-Math.PI/3,Math.PI/3,0,-2*Math.PI/3,2*Math.PI/3};
		
		for(int i=0;i<totalChildren;i++){
			BranchNode child = children.get(i);
			child.generateCoordinates();
			AffineTransform at = new AffineTransform();
			
			int ti = i+thetaOffset;
			double ntheta = ti<thetas.length?thetas[ti]:0;
			at.rotate(ntheta);
			at.translate(1, 0);
			
			
			child.applyTransform(at);
			
		}
		return this;
		
	}
	
	public BranchNode applyTransform(AffineTransform at){
		
		setSuggestedPoint(at.transform(getSuggestedPoint(), null));
		for(BranchNode c:this.getChildren()){
			c.applyTransform(at);
		}
		return this;
	}
	
	public String getSymbol(){
		if(this.isPseudoNode()){
			return "?";
		}else if(this.isRepeatNode()){
			return "'" + this.rep + "'";
		}
		return this.symbol;
	}
	
	public boolean shouldComineLinearly(){
		if(combineLinearly){
			return true;
		}
		return false;
	}
	
	public boolean hasChildren(){
		return !children.isEmpty();
	}
	public List<BranchNode> getChildren(){
		return this.children;
	}
	public boolean isRealNode(){
		return !pseudo && !isRepeat;
	}
	public boolean isPseudoNode(){
		return pseudo;
	}
	public boolean isRepeatNode(){
		return this.isRepeat;
	}
	public int getRepeat(){
		return this.rep;
	}
	public BranchNode setRepeat(int r){
		this.isRepeat=true;
		this.rep=r;
		return this;
	}
	
	public BranchNode getNodeForLinearCombine(){
		if(childForCombine!=null){
			return childForCombine;
		}
		return this;
	}
	
	public BranchNode setCombiningNode(BranchNode bn){
		this.childForCombine=bn;
		return this;
	}
	
	public BranchNode flagForCombining(){
		this.isCombiner=true;
		return this;
	}
	
	public BranchNode addChild(BranchNode bn){
		BranchNode ret=_addChild(bn);
		
		return ret;
	}
	
	public static BranchNode of(String s){
		return new BranchNode(s);
	}
	
	
	//Always add to the left
	//
	private BranchNode _addChild(BranchNode bn){
		BranchNode nToAddTo= this.getNodeForLinearCombine();
		
		if(bn.isCombiner){
			this.setCombiningNode(bn);
		}
		
		if(bn.isRepeatNode()){
			int rep=bn.getRepeat();
			if(this.shouldComineLinearly()){
				BranchNode parent=nToAddTo;
				for(int i=1;i<rep;i++){
					BranchNode bnew=new BranchNode(nToAddTo.symbol).thetaOffset(i%2);
					parent.addChild(bnew);
					parent=bnew;
				}
			}else{
				nToAddTo.setPseudoNode(true);
				for(int i=0;i<rep;i++){
					nToAddTo.addChild(new BranchNode(nToAddTo.symbol));
				}
			}
						
			if(bn.hasChildren()){
				//This usually means the text is inverted
				if(bn.getChildren().size()!=1){
					throw new IllegalStateException("Node:[" + bn.toString() + "] needs to be added, but has ambiguous meaning");
				}
				
				BranchNode nchild=bn.getChildren().get(0);
				return nchild.addChild(this);
			}
		}else if(bn.isPseudoNode()){
			nToAddTo.children.addAll(bn.children);
		}else if(!this.isPseudoNode() && nToAddTo.isTerminal() && !bn.isTerminal() && bn.isRealNode()){
			return bn.addChild(nToAddTo);
		}else{
			nToAddTo.children.add(bn);
		}
		return this;
	}
	
	public BranchNode setPseudoNode(boolean p){
		this.pseudo=p;
		return this;
	}
	
	public BranchNode getLeftBranchNode(){
		return rightBranch;
	}
	
	public BranchNode getRightBranchNode(){
		return leftBranch;
	}
	
	public BranchNode setRightBranchNode(BranchNode bn){
		this.rightBranch=bn;
		return this;
	}
	public BranchNode setLeftBranchNode(BranchNode bn){
		this.leftBranch=bn;
		return this;
	}
	
	public boolean canBeChain(){
		return this.leftBranch!=this.rightBranch;
	}
	
	public BranchNode(String s){
		this.symbol=s;
		if(this.symbol.equals("C")){
			combineLinearly=true;
		}
	}
	
	public BranchNode setShouldCombineLinearly(boolean t){
		this.combineLinearly=t;
		return this;
	}
	
	public boolean isLinkable(){
		return linkable;
	}
	
	public String toString(){
		String bondToParent="-";
		if(this.getOrderToParent()==2){
			bondToParent="=";
		}else if(this.getOrderToParent()==3){
			bondToParent="#";
		}
		if(this.hasChildren()){
			return bondToParent + this.getSymbol() + "(" + this.children.stream()
			                               .map(b->b.toString())
					                       .collect(Collectors.joining(",")) + ")";
		}
		return bondToParent+this.getSymbol();
	}
	
	private void forEachBranchNode(BiConsumer<BranchNode,BranchNode> cons, BranchNode parent){
		cons.accept(parent, this);
		this.getChildren()
		.stream()
		.collect(Collectors.toList())
		.forEach(child->{
			child.forEachBranchNode(cons,this);
		});
	}
	
	public void forEachBranchNode(BiConsumer<BranchNode,BranchNode> parentAndChild){
		forEachBranchNode(parentAndChild,null);
	}
	
	public static interface Token{
		public String getTokenName();
		public String getTokenPreferredStyle();
		public List<Token> getToken(String q);
		public String[] getForms();
		
		public static Token groupedToken(String tname,List<Token> tlist){
			Map<String,List<Token>> tokensForForms = new LinkedHashMap<String, List<Token>>();
			String[] forms = tlist.stream()
				 .map(tt->Tuple.of(tt,tt.getForms()))
				 .flatMap(tt->Arrays.stream(tt.v()).map(f->Tuple.of(tt.k(),f)))
				 .peek(tt->tokensForForms.computeIfAbsent(tt.v(), k->new ArrayList<>()).add(tt.k()))
				 .map(tt->tt.v())
			     .distinct()
			     .toArray(i->new String[i]);
			return new Token(){

				@Override
				public String getTokenName() {
					return tname;
				}

				@Override
				public List<Token> getToken(String q) {
					return tlist.stream()
					     .flatMap(tt->tt.getToken(q).stream())
					     .collect(Collectors.toList());
				}
				@Override
				public String[] getForms() {
					return forms;
				}

				@Override
				public String getTokenPreferredStyle() {
					return null;
				}
				@Override
				public List<Tuple<Token, String>> getFirstMatchingTokens(String s){
					
					List<Tuple<Token,String>> found = new ArrayList<>();
					for(String f:forms){
						if(s.startsWith(f)){
							tokensForForms.get(f)
							.forEach(fl->{
								found.add(Tuple.of(fl,f));	
							});
						}
					}
					return found;
				}
				
			};
		}
		public default List<Tuple<Token, String>> getFirstMatchingTokens(String s){
			String[] forms = this.getForms();
			Map<String,List<Token>> tokensForForms = new LinkedHashMap<String, List<Token>>();
			for(String f:forms){
				tokensForForms.computeIfAbsent(f, k->new ArrayList<>())
				.addAll(getToken(f));
			}
			List<Tuple<Token,String>> found = new ArrayList<>();
			for(String f:forms){
				if(s.startsWith(f)){
					tokensForForms.get(f)
					.forEach(fl->{
						found.add(Tuple.of(fl,f));	
					});
				}
			}
			return found;
		}
		
		public default Token compoundToken(Token t2){
			String[] forms1 =this.getForms();
			String[] forms2 =t2.getForms();
			List<String> nformsCombined = new ArrayList<String>();
			for(int i=0;i<forms1.length;i++){
				for(int j=0;j<forms2.length;j++){
					nformsCombined.add(forms1[i] + forms2[j]);
				}
			}
			return SimpleToken.of(this.getTokenName() + t2.getTokenName(), this.getTokenPreferredStyle() + t2.getTokenPreferredStyle(), nformsCombined.toArray(new String[0]));
		}
	}
	
	public static class SimpleToken implements Token{
		String[] forms;
		String name;
		String pref;
		
		
		@Override
		public String getTokenName() {
			return name;
		}
		@Override
		public List<Token> getToken(String q) {
			for(int i=0;i<forms.length;i++){
				if(forms[i].equals(q)){
					return Arrays.asList(this);
				}
			}
			return Arrays.asList();
		}
		
		@Override
		public String[] getForms() {
			return forms;
		}
		@Override
		public String getTokenPreferredStyle() {
			return pref;
		}
		
		public SimpleToken(String name, String pref, String[] forms){
			this.name=name;
			this.pref=pref;
			this.forms=forms;
		}
		
		public static SimpleToken of(String name, String pref, String ...forms){
			return new SimpleToken(name,pref,forms);
		}
		public static SimpleToken of(String name){
			return of(name, name, name);
		}
		
		public boolean equals(Object o){
			if(o instanceof Token){
				if(((Token)(o)).getTokenName().equals(this.getTokenName())){
					return true;
				}
			}
			return false;
			
		}
		public int hashCode(){
			return this.getTokenName().hashCode();
			
		}
	}
	
	
	private static Map<String, Token> tokenList = new LinkedHashMap<String,Token>();
	private static Map<String, Token> masterTokenList = new LinkedHashMap<String,Token>();
	
	private static List<Token> atomicSet = new ArrayList<>();
	private static List<Token> saturatedAtomicSet = new ArrayList<>();
	private static List<Token> numericSet = new ArrayList<>();
	private static void initializeTokenSet(){
		
		atomicSet.add(SimpleToken.of("C", "C", "C","c"));
		atomicSet.add(SimpleToken.of("O", "O", "O", "()", "0", "o"));
		atomicSet.add(SimpleToken.of("N", "N", "N", "1O"));
		atomicSet.add(SimpleToken.of("H", "H", "H", "tt", "1t", "t1", "I1", "t4", "11"));
		atomicSet.add(SimpleToken.of("S", "S", "S", "s"));
		atomicSet.add(SimpleToken.of("P", "P", "P", "p"));
		atomicSet.add(SimpleToken.of("Pt", "Pt", "Pt", "pt"));
		atomicSet.add(SimpleToken.of("K", "K", "K", "k"));
		atomicSet.add(SimpleToken.of("As", "As", "As", "AS", "A8"));
		atomicSet.add(SimpleToken.of("Au", "Au", "Au", "AU"));
		atomicSet.add(SimpleToken.of("Al", "Al", "Al", "A1", "At", "AI"));
		atomicSet.add(SimpleToken.of("F", "F", "F", "f"));
		atomicSet.add(SimpleToken.of("Cl", "Cl", "Cl", "CT", "Ct", "C)", "CI", "C1","cl", "cT", "ct", "c)", "cI", "c1"));
		atomicSet.add(SimpleToken.of("Br", "Br", "Br", "Sr", "sr", "8r", "BT"));
		atomicSet.add(SimpleToken.of("Na", "Na", "Na"));
		atomicSet.add(SimpleToken.of("I", "I", "t","1"));
		atomicSet.add(SimpleToken.of("B", "B", "B"));
		atomicSet.add(SimpleToken.of("Si", "Si", "Si", "SI", "Sl", "S1", "St", "si", "sI", "sl", "s1", "st", "8i", "8I", "8l", "81", "8t"));
		atomicSet.add(SimpleToken.of("Hg", "Hg", "Hg", "ttg", "1tg", "t1g", "I1g", "t4g", "11g"));
		atomicSet.add(SimpleToken.of("D", "D", "D"));
		
		atomicSet.forEach(tt->registerToken(tt,false));
		
		registerToken(Token.groupedToken("Atomic", atomicSet),true);
		
		
		numericSet.add(SimpleToken.of("1", "1", "1","t"));
		numericSet.add(SimpleToken.of("2", "2", "2"));
		numericSet.add(SimpleToken.of("3", "3", "3"));
		numericSet.add(SimpleToken.of("4", "4", "4"));
		numericSet.add(SimpleToken.of("5"));
		numericSet.add(SimpleToken.of("6"));
		numericSet.add(SimpleToken.of("7"));
		numericSet.add(SimpleToken.of("8"));
		numericSet.add(SimpleToken.of("9"));
		
		numericSet.forEach(tt->registerToken(tt,false));
		registerToken(Token.groupedToken("Numeric", numericSet),true);
		
		if(false){
			saturatedAtomicSet.add(getCombinedToken("C", "H", "3"));
			saturatedAtomicSet.add(getCombinedToken("C", "H", "2"));
			saturatedAtomicSet.add(getCombinedToken("C"));
			saturatedAtomicSet.add(getCombinedToken("N", "H", "2"));
			saturatedAtomicSet.add(getCombinedToken("N", "H"));
			saturatedAtomicSet.add(getCombinedToken("O", "H"));
			saturatedAtomicSet.add(getCombinedToken("S", "H"));
			saturatedAtomicSet.add(getCombinedToken("H", "3", "C"));
			saturatedAtomicSet.add(getCombinedToken("H", "2", "C"));
			saturatedAtomicSet.add(getCombinedToken("H", "C"));
			saturatedAtomicSet.add(getCombinedToken("H", "O"));
			saturatedAtomicSet.add(getCombinedToken("H", "2", "N"));
			saturatedAtomicSet.add(getCombinedToken("H", "N"));
			saturatedAtomicSet.add(getCombinedToken("H", "S"));
			saturatedAtomicSet.forEach(tt->registerToken(tt,false));
		}
		
		
		//Special Tokens:
		//Me	Me, Mc, MC
		//Et	Et
		//Pr	Pr
		//Ms	Ms, MS, M8
		//Ac	Ac, AC
		//Bn	Bn, Bt1
		//Cbz	CbZ, Cbz, cbZ, cbz, C6Z, c6Z, C6z, c6z 
		//tBu	t-Bu, tbu, t-Bo, tBo
		//nBu	n-Bu, nbu, n-Bo, nBo
		//Cys	CyS,Cys,cyS,cys,Cy8,cy8
		//Ph	Ph, ph, Pb, pb
		//pTol	p-tol
		//Ts	Ts, TS
		//Tf	Tf, Tr
		//PMBN	pMBN,PMBN
		//PLUS	+
		//OPEN	(
		//CLOSE	)
		List<Token> specialSet = new ArrayList<>();
		specialSet.add(SimpleToken.of("Me","Me","Me","Mc", "MC"));
		specialSet.add(SimpleToken.of("Et"));
		specialSet.add(SimpleToken.of("Pr"));
		specialSet.add(SimpleToken.of("Ms","Ms","Ms","MS","M8"));
		specialSet.add(SimpleToken.of("Ac","Ac","Ac","AC"));
		specialSet.add(SimpleToken.of("Bn","Bn","Bn","Bt1"));
		specialSet.add(SimpleToken.of("Cbz","Cbz", "CbZ", "Cbz", "cbZ", "cbz", "C6Z", "c6Z", "C6z", "c6z"));
		specialSet.add(SimpleToken.of("tBu","tBu", "t-Bu","tbu","t-Bo","tBo"));
		specialSet.add(SimpleToken.of("nBu","nBu", "n-Bu","nbu","n-Bo","nBo"));
		specialSet.add(SimpleToken.of("Cys","Cys", "CyS", "Cys", "cyS", "cys", "Cy8", "cy8"));
		specialSet.add(SimpleToken.of("Ph","Ph", "Ph", "ph", "Pb" , "pb"));
		specialSet.add(SimpleToken.of("pTol","p-tol", "p-tol", "ptol", "pTol"));
		specialSet.add(SimpleToken.of("Ts","Ts", "Ts", "TS"));
		specialSet.add(SimpleToken.of("Tf","Tf", "Tf", "Tr"));
		specialSet.add(SimpleToken.of("PMBN","PMBN", "PMBN", "pMBN"));
		specialSet.add(SimpleToken.of("PLUS","+", "+"));
		specialSet.add(SimpleToken.of("OPEN","(", "("));
		specialSet.add(SimpleToken.of("CLOSE",")", ")"));
		
		specialSet.forEach(tt->registerToken(tt,false));
		registerToken(Token.groupedToken("Special", specialSet),true);
	}
	
	private static void registerToken(Token t, boolean master){
		if(master){
			masterTokenList.put(t.getTokenName(), t);
		}
		tokenList.put(t.getTokenName(), t);
	}
	
	private static Token getCombinedToken(String ... tokNames){
		return Arrays.stream(tokNames)
						.map(s->tokenList.get(s))
						.reduce((t1,t2)->t1.compoundToken(t2))
						.orElse(null);
	}
	
	
	public static class TokenTree{
		private Token root;
		private int invalidChildren = 0;
		private boolean isInvalid=false;
		private List<TokenTree> children= new ArrayList<TokenTree>();
		
		
		public TokenTree(Token t){
			this.root=t;
		}
		
		public TokenTree addChild(TokenTree t){
			if(t.isInvalid()){
				invalidChildren++;
			}else{
				children.add(t);
			}			
			return this;
		}
		
		public List<TokenTree> asEnumeratedList(){
			List<TokenTree> uTrees = new ArrayList<>();
			getAllTrees((tl)->{
				TokenTree rparent=null;
				TokenTree parent = null;
				for(TokenTree c: tl){
					TokenTree clone = new TokenTree(c.root);
					if(parent!=null){
						parent.addChild(clone);
					}else{
						rparent=clone;
					}
					parent=clone;
				}
				uTrees.add(rparent);
			});
			return uTrees;
		}
		
		public void getAllTrees(Consumer<List<TokenTree>> cons){
			if(!children.isEmpty()){
				for(TokenTree tt: children){
					tt.getAllTrees(ll->{
						List<TokenTree> mlist = new ArrayList<>();
						mlist.add(this);
						mlist.addAll(ll);
						cons.accept(mlist);
					});
				}
			}else{
				cons.accept(Arrays.asList(this));
			}
		}
		public void printAlllTrees(){
			this.getAllTrees(ll->{
				String full=ll.stream()
				  .map(tt1->tt1.getDisplay())
				  .collect(Collectors.joining("-"));
				System.out.println(full);
			});
		}
		public Token getCurrentToken(){
			return root;
		}
		
		public String getDisplay(){
			if(this.root==null)return "null";
			return this.root.getTokenPreferredStyle();
		}
		public TokenTree markInvalid(){
			this.isInvalid=true;
			return this;
		}
		public boolean isInvalid(){
			if(isInvalid)return true;
			if(invalidChildren>0 && children.isEmpty())return true;
			return false;
		}

		public boolean hasChildren() {
			return !this.children.isEmpty();
		}

		public int maxLength() {
			int sofar=0;
			if(this.root==null){
				sofar=0;
			}else{
				sofar=1;
			}
			
			return sofar + children.stream()
							.map(tt->tt.maxLength())
							.max(Comparator.naturalOrder())
							.orElse(0);
		}

		public Tuple<TokenTree, TokenTree> splitTreeAtDepth(int i) {
			if(this.root==null)i++;
			TokenTree shallowCopy=depthLimitedCopy(i-1);
			TokenTree deepCopy=new TokenTree(null);
			consumeAllNodesAtDepth((tn)->{
				deepCopy.addChild(tn);
			}, i);
			
			return Tuple.of(shallowCopy,deepCopy);
		}
		public Tuple<TokenTree, TokenTree> splitTreeToRemaingDepth(int i) {
			
			TokenTree deepCopy=new TokenTree(null);
			Set<TokenTree> trimAt = new HashSet<TokenTree>();
			
			consumeAllNodesAtFirstMatch((tn)->{
				for(TokenTree child:tn.children){
					deepCopy.addChild(child);	
				}
				trimAt.add(tn);
			}, (tn)->tn.maxLength()==i+1);
			
			TokenTree shallowCopy = this.copyTruncatedAt(tn->{
				return trimAt.contains(tn);
			});
			
			
			return Tuple.of(shallowCopy,deepCopy);
		}
		public void consumeAllNodesAtDepth(Consumer<TokenTree> cons, int initDepth){
			if(initDepth==0){
				cons.accept(this);
			}else{
				for(TokenTree child:children){
					child.consumeAllNodesAtDepth(cons, initDepth-1);
				}
			}
		}
		
		public void consumeAllNodesAtFirstMatch(Consumer<TokenTree> cons, Predicate<TokenTree> stop){
			if(stop.test(this)){
				cons.accept(this);
			}else{
				for(TokenTree child:children){
					child.consumeAllNodesAtFirstMatch(cons, stop);
				}
			}
		}
		
		public TokenTree depthLimitedCopy(int maxDepth){
			TokenTree tnew = new TokenTree(this.root);
			if(maxDepth>0){
				for(TokenTree child:children){
					tnew.addChild(child.depthLimitedCopy(maxDepth-1));
				}
			}
			return tnew;
		}
		
		public TokenTree copyTruncatedAt(Predicate<TokenTree> truncateAfter){
			TokenTree tnew = new TokenTree(this.root);
			
			if(!truncateAfter.test(this)){
				for(TokenTree child:children){
					tnew.addChild(child.copyTruncatedAt(truncateAfter));
				}
			}
			return tnew;
		}
		
		public List<List<Token>> getAllTokenPathsWhichAllMatch(Predicate<Token> keep){
			List<List<Token>> ttlist = new ArrayList<>();
			
			getAllTrees(lt->{
				List<Token> tlist=lt.stream()
						  .map(tt->tt.getCurrentToken())
						  .filter(tok->tok!=null)
						  .collect(Collectors.toList());
				if(tlist.size()>0){
					if(tlist.stream().allMatch(keep)){
						ttlist.add(tlist);
					}
				}
			});
			return ttlist;
		}
	}
	
	private static TokenTree parseTokenTree(TokenTree starts, String t){
		if(t.length()==0)return starts;
		AtomicBoolean foundany = new AtomicBoolean(false);
		
		masterTokenList.forEach((s,tt)->{
			List<Tuple<Token,String>> parsed = tt.getFirstMatchingTokens(t);
			parsed.stream()
			.forEach(tup->{
				foundany.set(true);
				Token ftok=tup.k();
				
				TokenTree found = new TokenTree(ftok);
				
				String remainder = t.substring(tup.v().length());
//				System.out.println("Found:" + ftok.getTokenName() + " in " + t + " now check:" + remainder);
				parseTokenTree(found,remainder);
				
				starts.addChild(found);				
			});
		});
		if(!foundany.get()){
			starts.addChild(new TokenTree(null).markInvalid());
		}
		return starts;
	}
	
	public static TokenTree parseTokenTree(String t){
		TokenTree start = new TokenTree(null);
		return parseTokenTree(start,t);
	}
	
	static{
		initializeTokenSet();
	}
	
	
	private static List<ParsingRule> parsingRules = new ArrayList<>();
	
	
	
	
	public static interface ParsingRule{
		public String getRuleName();
		public Optional<Tuple<BranchNode, String>> parse(TokenTree tt);
	}
	
	public static class RegexTokenParsingRule implements ParsingRule{
		
		private String name;
		private Pattern p;
		private Function<Matcher,Tuple<BranchNode,String>> nodeMaker;
		private String returnAs=null;
		
		
		
		public RegexTokenParsingRule(String name,String pattern, Function<Matcher, BranchNode> maker){
			this(name,pattern
					 ,(m)->Optional.ofNullable(maker.apply(m)).map(b->Tuple.of(b,(String)null)).orElse(null)
					 ,true);
		}
		
		
		public RegexTokenParsingRule(String name,String pattern, Function<Matcher,Tuple<BranchNode,String>> maker, boolean change){
			p=Pattern.compile(pattern);
			nodeMaker=maker;
			this.name=name;
		}
		
		public RegexTokenParsingRule setReturn(String ret){
			this.returnAs=ret;
			return this;
		}

		@Override
		public String getRuleName() {
			return name;
		}

		@Override
		public Optional<Tuple<BranchNode, String>> parse(TokenTree tt) {
			return tt.getAllTokenPathsWhichAllMatch((t)->true)
			  .stream()
			  .map(tl->Tuple.of(tl,tl.stream().map(t->t.getTokenName()).collect(Collectors.joining())))
			  .map(Tuple.vmap(s->p.matcher(s)))
			  .filter(m->m.v().matches())
			  .findAny()
			  .map(Tuple.vmap(m->nodeMaker.apply(m)))
			  .filter(t->t.v()!=null)
			  .map(b->Tuple.of(b.v(),b.k()))
			  .map(t->{
				  String retcalc = t.k().v();
				  if(returnAs!=null){
					  return Tuple.of(t.k().k(),returnAs);
				  }else if(retcalc!=null){
					  return Tuple.of(t.k().k(),retcalc);
				  }else{
					  return Tuple.of(t.k().k(),t.v().stream().map(tok->tok.getTokenPreferredStyle()).collect(Collectors.joining()));
				  }
			  });
		}
		
	}
	
	public static class TemplateTokenParsingRule implements ParsingRule{
		String name;
		String returnAs;
		List<Token> expectedTokens;
		Supplier<BranchNode> bnMaker;
		@Override
		public String getRuleName() {
			return name;
		}
		
		public TemplateTokenParsingRule (String name, String returnAs, List<Token> expectedTokens, Supplier<BranchNode> bnmaker){
			this.name=name;
			this.returnAs=returnAs;
			this.expectedTokens=expectedTokens;
			this.bnMaker=bnmaker;
		}
		
		public static TemplateTokenParsingRule fromTokenShorthand(String name, String returnAs, String tokenString, Supplier<BranchNode> bnmaker){
			String syntax = tokenString.replace(" ", "").replace("[", "").replace("]", "").replace("><", ",").replace("<","").replace(">", "");
			String[] tokens =  syntax.split(",");
			List<Token> resolved = Arrays.stream(tokens)
					                     .map(t->tokenList.get(t))
					                     .collect(Collectors.toList());
			
			return new TemplateTokenParsingRule(name, returnAs, resolved, bnmaker);
			
		}
		
		
		
		@Override
		public Optional<Tuple<BranchNode, String>> parse(TokenTree tt) {
			
			if(tt.root==null){
				for(TokenTree tchild: tt.children){
					Optional<Tuple<BranchNode,String>> parsed = parse(tchild);
					if(parsed.isPresent())return parsed;
				}
			}else{
				boolean matched = matches(tt,expectedTokens,0);
				if(matched){
					return Optional.of(Tuple.of(bnMaker.get(), returnAs));
				}
			}
			return Optional.empty();
		}
		
		private static boolean matches(TokenTree tt, List<Token> expected, int offset){
			Token expectedToken = expected.get(offset);
//			System.out.println("Token:" + offset + " of " + expected.size());
//			System.out.println("Expecting:" + expectedToken.getTokenName());
			if(tt.getCurrentToken().getTokenName()
					.equals(expectedToken.getTokenName())){
				if(offset+1==expected.size()){
					return !tt.hasChildren();
				}else{
					for(TokenTree tchild: tt.children){
						boolean worked = matches(tchild, expected, offset+1);
						if(worked)return true;
					}
				}
			}
			return false;
		}
		
	}
	
	
/*
Atomic Tokens:
C	C,c
N	N, 1O
O	(),0,O,o
H	H, tt, 1t, t1, I1, t4, 11
S	S, s
P	P, p
Pt	Pt, pt
K	K, k
As	As, AS, A8
Au	Au, AU
Al	Al, A1, At, AI
F	F, f
Cl	Cl, CT, Ct, C), CI, C1
Br	Br, Sr, sr, 8r, BT
Na	Na
I	I*, t*, 1*
B	B
Si	Si, SI, Sl, S1, St, si, sI, sl, s1, st, 8i, 8I, 8l, 81, 8t
Hg	Hg, ttg, 1tg, t1g, I1g, t4g, 11g
D	D


Numeric Tokens:
1	1
2	2
3	3
4	4
5	5
6	6
7	7
8	8
9	9

Special Tokens:
Me	Me, Mc, MC
Et	Et
Pr	Pr
Ms	Ms, MS, M8
Ac	Ac, AC
Bn	Bn, Bt1
Cbz	CbZ, Cbz, cbZ, cbz, C6Z, c6Z, C6z, c6z 
tBu	t-Bu, tbu, t-Bo, tBo
nBu	n-Bu, nbu, n-Bo, nBo
Cys	CyS,Cys,cyS,cys,Cy8,cy8
Ph	Ph, ph, Pb, pb
pTol	p-tol
Ts	Ts, TS
Tf	Tf, Tr
PMBN	pMBN,PMBN
PLUS	+
OPEN	(
CLOSE	)


Parsing Groups:
Alkyl_Chain_Linker	[<C><H><2>]+
Ester_Linker	[<C><O><O><C><H><2>]
Alkyl_Chain_Linker_2	[<C><H><2>]+<C>
N_Ketone_Linker	[<N><H><C><O>]
N_Keto_Methyl	[<C><H><3><C><O><N><H>]
CARBOXYLIC_ACID_F1	[<C><O><2><H>]	CO2H
CARBOXYLIC_ACID_F2	[<C><O><O><H>]	COOH
CARBOXYLIC_ACID_F3	[<H><O><O><C>]	HOOC
CARBOXYLIC_ACID_F4	[<H><O><2><C>]	HO2C
ESTER_LINKER_F1	[<C><O><2>]	CO2
ESTER_LINKER_F2	[<C><O><O>]	COO
ESTER_LINKER_F3	[<O><O><C>] 	OOC
ESTER_LINKER_F4	[<O><O><O>]	COO
CYANIDE_F1	[<C><N>]	CN
CYANIDE_F2	[<N><C>]	NC
SULFATE	[<S><O><3>] || [<O><3><S>]	SO3
NITRATE	[<N><O><2>] || [<O><2><N>]	NO2
METHYL_1	[<C><H><3>]	CH3
METHYL_2	[<Me>]	Me
ETHYL	[<Et>]	Et
PROPYL	[<Pr>]	Pr
N_METHYL_F1	[<N><C><H><3>]	NCH3
N_METHYL_F2	[<N><H><C><H><3>]	NHCH3
SULFONE_F1	[<S><O><2>]	SO2
SULFONE_F2	[<2><O><S>]	O2S
METHYL_SULFONE	<Ms>	Ms
BOC	[<B><O><C>]	Boc
BOC_2N	[<BOC><2><N>]	Boc2N
AZIDE	[<N><3>]	N3
KETO	[<C><O>]	CO
ALCOHOL	[<O><H>]	OH
METHOXY_LINK	[<C><H><2><O>]	CH2O
ACETATE	[<Ac>]	Ac
ETHYL_ESTER	[<Et><O><2><C>]	EtO2C
BENZO	[<Bn>]	Bn
DOUBLE_BOND_NITROGEN	[<B><N>]	HN
BENZYL_CARBAMATE	[<Cbz>]	Cbz
TERT_BUTYL	[<tBu>]	t-Bu
N_BUTYL	[<nBu>]	n-Bu
CYSTEINE	[<Cys>]	Cys
CYSTEINE_ESTER	[<ESTER_LINKER_F1><Cys>]	Cys
PHENYL	[<Ph>]	Ph
PARA_TOLUENE	[<pTol>]	p-tol
PARA_SULFO_TOLUENE	[<Ts>]	Ts
PARA_MBN	[<PMBN>]	PMBN
METHOXY_F1	[<C><H><3><O>]	CH3O
METHOXY_F2	[<H><3><C><O>]	H3CO
N_ETHYL	[<Et><H><N>]	EtHN
ALDEHYDE	[<O><H><C>]	OHC
NOISY_O	[<I>+<O>]	O
NOISY_N	[<I>+<N>]	N
S_PLUS	[<S><PLUS>]	S
N_PLUS	[<N><PLUS>]	N
REVERSE_O_BENZO	[<Bn><O>]	BnO
REVERSE_OME	[<Me><O>]	MeO
REVERSE_N_BENZO	[<Bn><N>]	BnN
REVERSE_N_METHYL	[<Me><N>]	MeN
REVERSE_S_METHYL	[<Me><S>]	MeS
REVERSE_OXY_ETHYL	[<Et><O>]	EtO
REVERSE_OCF3	[<F><3><C><O>]	F3CO
REVERSE_ETHYL_ESTER	[<Et><O><O><C>]	EtOOC
REVERSE_S_METHYL	[<H><3><C><S>]	H3CS
REVERSE_TRIPHEYNYL	[<Ph><3><C>]	Ph3C
REVERSE_TERM_ALCOHOL	[<H><O><H><2><C>]	HOH2C
TERM_ALCOHOL	[<C><H><2><O><H>]	CH2OH
SULFOXY_METHYL	[<S><O><C><H><3>]	SOCH3
N_CHAIN	[<OPEN><!OPEN>+,<CLOSE><Numeric>+]	???
N_CHAIN_N	[<OPEN><!OPEN>+,<CLOSE><Numeric>+<N>]	???
SO2CF3	[<Tf>]	Tf
METHYL_ESTER	[<Me><O><O><C>]	MeOOC
REVERESE_NHCH3	[<H><3><C><H><N>]	H3CHN
REVERESE_OXY_ALKYL	[<C><Numeric>+<O>]	???
N_C_COOH	[<N><H><C><H><2><CARBOXYLIC_ACID_F2>]	NHCH2COOH
REVERSE_SO2CH3	[<H><3><C><O><2><S>]	H3CO2S
REVERSE_ALKYL_OXY_ALKYL	[<H><3><C>(<C><H><2>)*<O>(<C><H><2>)*]	???	
ALKYL_OXY_ALKYL	[<C><H><2>(<C><H><2>)*<O>(<C><H><2>)*<C><H><3>]	???	
N_N_LINKER	[<N><N>]	???
ATOMIC_SYMBOL	[<Atomic>]	???
NUMERIC_SYMBOL	[<Numeric>+]	???

 */
	
	//TODO: implement shorthand
	public static BranchNode fromShortHand(String sh){
		//"-?(-C(-C),-C(-C))"
		//[-=#](<WEDGE>)?<SYMBOL>([:]<LOCUS>)?(!<THETA_OFF>)?([-=#](<LOCUS>)?)?(<CHILD>(,<CHILD>)*)
		//WEDGE	@,@@
		//SYMBOL [A-Z][a-z]
		//BOND_TYPE [-=#]
		//LOCUS [1-9][0-9]+
		//TEHTA_OFF [1-9]
		
		
		return null;
	}
	
	
	
	
	
	static{
		
		
		//Alkyl_Chain_Linker	[<C><H><2>]+
		parsingRules.add(new RegexTokenParsingRule("Alkyl_Chain_Linker","((CH[2]*)+)", (m)->{
			Matcher mm=m;
			
			String cc=mm.group(1);
			
			int c=cc.split("H[2|3]*").length;
			
			BranchNode full = of("C");
			BranchNode parent=full;
			for(int i=1;i<c;i++){
				BranchNode nn=of("C").thetaOffset(i%2);
				parent.addChild(nn);
				parent=nn;
			}
			full.setRightBranchNode(parent);
			return full;
		}));
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("Ester_Linker", "COOCH2", "[<C><O><O><C><H><2>]", ()->{
			BranchNode full = of("C").addChild(of("O").setOrderToParent(2).thetaOffset(1));
			BranchNode OMethyl = of("O").setOrderToParent(1).thetaOffset(1);
			BranchNode meth = of("C").setOrderToParent(1);
			OMethyl.addChild(meth);
			full.addChild(OMethyl);
			full.setRightBranchNode(meth);
			full.setCombiningNode(meth);
			return full;
		}));
		
		
	
		parsingRules.add(new RegexTokenParsingRule("Alkyl_Chain_Linker","((CH[2]*)+)C", (m)->{
			Matcher mm=m;
			String cc=mm.group(1);
			int c=cc.split("H[2|3]*").length;
			BranchNode full = new BranchNode("C");
			BranchNode nfull = new BranchNode("C");
			full.addChild(nfull);
			BranchNode parent=nfull;
			for(int i=1;i<c;i++){
				BranchNode nn=new BranchNode("C").thetaOffset(i%2);
				parent.addChild(nn);
				parent=nn;
			}
			full.setRightBranchNode(parent);
			return full;
		}));
		
		parsingRules.add(new RegexTokenParsingRule("HH_START_FIXER","H(H.*)", (m)->{
			Matcher mm=m;
			String cc=mm.group(1);
			return parseBranchNode(cc).orElse(null);
		},true));
		
	
		
		
		//N_Ketone_Linker
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("N_Ketone_Linker", "NHCO", "[<N><H><C><O>]", ()->{
			BranchNode full = of("N");
			BranchNode carb = of("C").addChild(new BranchNode("O").setOrderToParent(2));
			full.addChild(carb);
			full.setRightBranchNode(carb);
			full.setCombiningNode(carb);
			return full;
		}));
		
		
		//N_Keto_Methyl
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("N_Keto_Methyl", "CH3CONH", "[<C><H><3><C><O><N><H>]", ()->{
			return new BranchNode("N").addChild(parseBranchNode("CO").get().k().addChild(new BranchNode("C")));
		}));
		
		Supplier<BranchNode> carboxylicAcidSupplier = ()->{
			BranchNode bn = new BranchNode("C");
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			BranchNode oalcohol=new BranchNode("O").setOrderToParent(1).flagForCombining();
			bn.addChild(oalcohol);
			bn.setRightBranchNode(oalcohol);
			return bn;
		};
		//CARBOXYLIC_ACID_F1	[<C><O><2><H>]	CO2H
		//CARBOXYLIC_ACID_F2	[<C><O><O><H>]	COOH
		//CARBOXYLIC_ACID_F3	[<H><O><O><C>]	HOOC
		//CARBOXYLIC_ACID_F4	[<H><O><2><C>]	HO2C
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("CARBOXYLIC_ACID_F1", "CO2H", "[<C><O><2><H>]",carboxylicAcidSupplier));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("CARBOXYLIC_ACID_F2", "COOH", "[<C><O><O><H>]",carboxylicAcidSupplier));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("CARBOXYLIC_ACID_F3", "HOOC", "[<H><O><O><C>]",carboxylicAcidSupplier));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("CARBOXYLIC_ACID_F4", "HO2C", "[<H><O><2><C>]",carboxylicAcidSupplier));
		
		//ESTER_LINKER_F1	[<C><O><2>]	CO2
		//ESTER_LINKER_F2	[<C><O><O>]	COO
		//ESTER_LINKER_F3	[<O><O><C>]	OOC
		//ESTER_LINKER_F4	[<O><O><O>]	COO
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("ESTER_LINKER_F1", "CO2", "[<C><O><2>]",carboxylicAcidSupplier));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("ESTER_LINKER_F2", "COO", "[<C><O><O>]",carboxylicAcidSupplier));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("ESTER_LINKER_F3", "OOC", "[<O><O><C>]",carboxylicAcidSupplier));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("ESTER_LINKER_F4", "O2C", "[<O><2><C>]",carboxylicAcidSupplier));
		
		
		Supplier<BranchNode> cyanoSupplier = ()->{
			BranchNode bn = new BranchNode("C").thetaOffset(2);
			bn.addChild(new BranchNode("N").setOrderToParent(3));
			//bn.addChild(new BranchNode("O").setOrderToParent(1).flagForCombining());
			return bn;
		};
		
		//CYANIDE_F1	[<C><N>]	CN
		//CYANIDE_F2	[<N><C>]	NC
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("CYANIDE_F1", "CN", "[<C><N>]",cyanoSupplier));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("CYANIDE_F2", "NC", "[<N><C>]",cyanoSupplier));
		
		Supplier<BranchNode> sulfateSupplier = ()->{
			BranchNode bn = new BranchNode("S");
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			bn.addChild(new BranchNode("O").setOrderToParent(1).flagForCombining());
			return bn;
		};
		
		
		
		//SULFATE	[<S><O><3>] || [<O><3><S>]	SO3
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("SULFATE_F1", "SO3", "[<S><O><3>]",sulfateSupplier));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("SULFATE_F2", "O3S", "[<O><3><S>]",sulfateSupplier));
		
		
		Supplier<BranchNode> nitroSupplier = ()->{
			BranchNode bn = new BranchNode("N").setCharge(1);
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			bn.addChild(new BranchNode("O").setOrderToParent(1).setCharge(-1));
			//bn.addChild(new BranchNode("O").setOrderToParent(1).flagForCombining());
			return bn;
		};
		
		//NITRATE	[<N><O><2>] || [<O><2><N>]	NO2
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("NITRO_F1", "NO2", "[<N><O><2>]",nitroSupplier));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("NITRO_F2", "O2N", "[<O><2><N>]",nitroSupplier));
		
		
		//METHYL_1	[<C><H><3>]	CH3
		//METHYL_2	[<Me>]	Me		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("METHYL_F1", "CH3", "[<C><H><3>]",()->new BranchNode("C").setShouldCombineLinearly(false).setTerminal(true)));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("METHYL_F2", "Me", "[<Me>]",()->new BranchNode("C").setShouldCombineLinearly(false).setTerminal(true)));
		

		//ETHYL	[<Et>]	Et
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("ETHYL", "Et", "[<Et>]",()->new BranchNode("C").addChild(new BranchNode("C")).setTerminal(true)));
		
		//PROPYL	[<Pr>]	Pr
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("PROPYL", "Pr", "[<Pr>]",()->new BranchNode("C").addChild(new BranchNode("C")
				                                                                                                                    .addChild(new BranchNode("C"))
				                                                                                                                    ).setTerminal(true)
				));
		
		
		//N_METHYL_F1	[<N><C><H><3>]	NCH3
		//N_METHYL_F2	[<N><H><C><H><3>]	NHCH3
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("N_METHYL_F1", "NCH3", "[<N><C><H><3>]",()->new BranchNode("N").addChild(new BranchNode("C").setTerminal(true))));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("N_METHYL_F2", "NHCH3", "[<N><H><C><H><3>]",()->new BranchNode("N").addChild(new BranchNode("C").setTerminal(true))));
		
		
		//SULFONE_F1	[<S><O><2>]	SO2
		//SULFONE_F2	[<2><O><S>]	O2S
		
		Supplier<BranchNode> sulfoneMaker = ()->{
			BranchNode bn = new BranchNode("S");
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			return bn;
		};
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("SULFONE_F1", "SO2", "[<S><O><2>]",sulfoneMaker));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("SULFONE_F2", "2OS", "[<2><O><S>]",sulfoneMaker));
		
		
		//METHYL_SULFONE	<Ms>	Ms
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("METHYL_SULFONE", "Ms", "[<Ms>]",()->{
			BranchNode bn = new BranchNode("S");
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			bn.addChild(new BranchNode("C").setOrderToParent(1));
			return bn;
		}));
		
		
		Supplier<BranchNode> bocMaker = ()->{
			BranchNode bn = new BranchNode("C");
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			bn.addChild(new BranchNode("O").setOrderToParent(1)
					.addChild(new BranchNode("C")
							.addChild(new BranchNode("C"))
							.addChild(new BranchNode("C"))
							.addChild(new BranchNode("C"))
							));
			return bn;
		};
		
		//BOC	[<B><O><C>]	Boc
		//BOC_2N	[<BOC><2><N>]	Boc2N
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("BOC", "Boc", "[<B><O><C>]",bocMaker));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("BOC_2N", "Boc2N", "[<B><O><C><2><N>]",()->{
			BranchNode b1=bocMaker.get();
			BranchNode b2=bocMaker.get();
			BranchNode bn = new BranchNode("N");
			bn.addChild(b1);
			bn.addChild(b2);
			return bn;
		}));
		
		//AZIDE	[<N><3>]	N3
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("AZIDE", "N3", "[<N><3>]",()->{
			BranchNode bn = new BranchNode("N");
			bn.addChild(new BranchNode("N").setOrderToParent(2).setCharge(1).thetaOffset(2)
					.addChild(new BranchNode("N").setOrderToParent(2).setCharge(-1)));
			return bn;	
		}));
		
		
		//KETO	[<C><O>]	CO
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("KETO", "CO", "[<C><O>]", ()->{
			return new BranchNode("C").addChild(new BranchNode("O").setOrderToParent(2));
		}));
		
				
		//ALCOHOL	[<O><H>]	OH
//		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("ALCOHOL", "OH", "[<O><H>]", ()->{
//			return new BranchNode("O").setTerminal(true);
//		}));
		
		//METHOXY_LINK	[<C><H><2><O>]	CH2O
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("METHOXY_LINK", "CH2O", "[<C><H><2><O>]", ()->{
			BranchNode bn = new BranchNode("C");
			BranchNode p=new BranchNode("O").setOrderToParent(1);
			bn.addChild(p);
			bn.setRightBranchNode(p);
			return bn;
		}));
		
		//ACETATE	[<Ac>]	Ac
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("ACETATE", "Ac", "[<Ac>]", ()->{
			BranchNode bn = new BranchNode("C").setOrderToParent(1)
					   .addChild(new BranchNode("C").setOrderToParent(1))
					   .addChild(new BranchNode("O").setOrderToParent(2));
			return bn;
		}));
		
		//O_ACETATE	[<Ac><O>]	AcO
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("O_ACETATE", "AcO", "[<Ac><O>]", ()->{
			BranchNode bn = new BranchNode("O");
			BranchNode ac=parseBranchNode("Ac").get().k();
			bn.addChild(ac);
			return bn;
		}));
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("ETHYL_ESTER", "EtO2C", "[<Et><O><2><C>]", ()->{
			BranchNode co2=parseBranchNode("CO2").get().k();
			BranchNode et=parseBranchNode("Et").get().k();
			co2.addChild(et);
			return co2;
		}));
		
		Supplier<BranchNode> phenSupplier = ()->{
			BranchNode bn = new BranchNode("C").thetaOffset(1);
			bn.addChild(new BranchNode("C").setOrderToParent(2)
					    			       .addChild(new BranchNode("C").setOrderToParent(1)
					    			    		   .addChild(new BranchNode("C").setOrderToParent(2)
					    			    				   .addChild(new BranchNode("C").setOrderToParent(1).addChild(new BranchNode("C").setOrderToParent(2).addRing(bn, 1)))
					    			    		   )
					    			    	)
					);
			return bn.setTerminal(true);
		};
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("BENZO", "Bn", "[<Bn>]", ()->{
			BranchNode carb=new BranchNode("C");
			BranchNode ben=phenSupplier.get();			
			carb.addChild(ben);
			return carb;
		}));
		
		//DOUBLE_BOND_NITROGEN	[<B><N>]	HN
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("DOUBLE_BOND_NITROGEN", "HN", "[<B><N>]", ()->{
			return new BranchNode("N");
		}));
		
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("BENZYL_CARBAMATE", "Cbz", "[<Cbz>]", ()->{
			BranchNode carb=new BranchNode("C").addChild(new BranchNode("O").setOrderToParent(2));
			BranchNode ox=new BranchNode("O").setOrderToParent(1);
			carb.addChild(ox);
			BranchNode ben=parseBranchNode("Bn").get().k();
			ox.addChild(ben);
			return carb;
		}));
		
		//TERT_BUTYL
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("TERT_BUTYL", "t-Bu", "[<tBu>]", ()->{
			BranchNode bn = new BranchNode("C");
			bn.addChild(new BranchNode("C"));
			bn.addChild(new BranchNode("C"));
			bn.addChild(new BranchNode("C"));
			return bn;
		}));
		
		//TERT_BUTYL
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("N_BUTYL", "n-Bu", "[<nBu>]", ()->{
			BranchNode bn = new BranchNode("C");
			bn.addChild(new BranchNode("C")
						     .addChild(new BranchNode("C")
						    		 .addChild(new BranchNode("C"))));
			return bn;
		}));
		
		//CYSTEINE	[<Cys>]	Cys
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("CYSTEINE", "Cys", "[<Cys>]", ()->{
			BranchNode bn = new BranchNode("N");
			BranchNode carbonyl=new BranchNode("C").addChild(new BranchNode("O").setOrderToParent(2));
			bn.addChild(new BranchNode("C").setOrderToParent(1)
										   .flagForCombining()
					                       .addChild(new BranchNode("C").thetaOffset(1).setWedgeToParent(1).addChild(new BranchNode("S")))
					                       .addChild(carbonyl)
					);
			bn.setCombiningNode(carbonyl);
			return bn;
		}));
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("CYSTEINE_ESTER", "CO2Cys", "[<C><O><2><Cys>]", ()->{
			BranchNode bn = new BranchNode("N");
			BranchNode carbonyl=new BranchNode("C").addChild(new BranchNode("O").setOrderToParent(2));
			bn.addChild(new BranchNode("C").setOrderToParent(1)
										   .flagForCombining()
					                       .addChild(new BranchNode("C").thetaOffset(1).setWedgeToParent(1).addChild(new BranchNode("S")))
					                       .addChild(carbonyl)
					);
			bn.setCombiningNode(carbonyl);
			return bn;
		}));
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("PHENYL", "Ph", "[<Ph>]", phenSupplier));
		
		//PARA_TOLUENE	[<pTol>]	p-tol
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("PARA_TOLUENE", "p-tol", "[<pTol>]", ()->{
			BranchNode bn = new BranchNode("C").thetaOffset(1);
			bn.addChild(new BranchNode("C").setOrderToParent(2) //ortho
					    			       .addChild(new BranchNode("C").setOrderToParent(1) //meta
					    			    		   .addChild(new BranchNode("C").setOrderToParent(2) //para
					    			    				   .addChild(new BranchNode("C").setOrderToParent(1) //meta
					    			    						   .addChild(new BranchNode("C").setOrderToParent(2).addRing(bn, 1))) //ortho
					    			    				   .addChild(new BranchNode("C"))
					    			    		   )
					    			    	)
					);
			
			
			return bn;
		}));
		
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("PARA_SULFO_TOLUENE", "Ts", "[<Ts>]", ()->{
			BranchNode ben=parseBranchNode("SO2").get().k();
			BranchNode tol=parseBranchNode("p-tol").get().k();
			ben.addChild(tol);
			return ben;
		}));
		

		//PARA_MBN	[<PMBN>]	PMBN
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("PARA_MBN", "PMBN", "[<PMBN>]", ()->{
			BranchNode bno = new BranchNode("N").thetaOffset(1);
			BranchNode ml = new BranchNode("C");
			bno.addChild(ml);
			
			BranchNode bn = new BranchNode("C").thetaOffset(1);
			
			bn.addChild(new BranchNode("C").setOrderToParent(2)
					    			       .addChild(of("C").setOrderToParent(1)
					    			    		   .addChild(of("C").setOrderToParent(2)
					    			    				   .addChild(of("C").setOrderToParent(1).addChild(of("C").setOrderToParent(2).addRing(bn, 1)))
					    			    				   .addChild(of("O").addChild(of("C")))
					    			    		   )
					    			    	)
					);
			ml.addChild(bn);
			return bno;
		}));
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("METHOXY_F1", "CH3O", "[<C><H><3><O>]", ()->{
			return parseBranchNode("OCH3").get().k();
		}));
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("METHOXY_F2", "OCH3", "[<O><C><H><3>]", ()->{
			return new BranchNode("O").addChild(new BranchNode("C"));
		}));
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("N_ETHYL", "EtNH", "[<Et><H><N>]", ()->{
			return parseBranchNode("NHEt").get().k();
		}));

		//ALDEHYDE	[<O><H><C>]	OHC
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("ALDEHYDE_F1", "OHC", "[<O><H><C>]", ()->{
			return parseBranchNode("COH").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("ALDEHYDE_F2", "CHO", "[<C><H><O>]", ()->{
			return parseBranchNode("COH").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("ALDEHYDE_F3", "COH", "[<C><O><H>]", ()->{
			return new BranchNode("C").addChild(new BranchNode("O").setOrderToParent(2));
		}));
		
		
		
		parsingRules.add(new RegexTokenParsingRule("NOISY_O","I+O", (m)->{
			return new BranchNode("O");
		}).setReturn("O"));
		
		parsingRules.add(new RegexTokenParsingRule("NOISY_N","I+N", (m)->{
			return new BranchNode("N");
		}).setReturn("N"));
		
		
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("S_PLUS", "S", "[<S><PLUS>]", ()->{
			return new BranchNode("S").setCharge(1);
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("N_PLUS", "N", "[<N><PLUS>]", ()->{
			return new BranchNode("N").setCharge(1);
		}));
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_O_BENZO", "BnO", "[<Bn><O>]", ()->{
			return parseBranchNode("OBn").get().k();
		}));
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_OME", "MeO", "[<Me><O>]", ()->{
			return parseBranchNode("OMe").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_HO", "HO", "[<H><O>]", ()->{
			return parseBranchNode("OH").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_SO3H", "HO3S", "[<H><O><3><S>]", ()->{
			return parseBranchNode("SO3H").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_HN", "HN", "[<H><N>]", ()->{
			return parseBranchNode("NH").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_HC", "HC", "[<H><C>]", ()->{
			return parseBranchNode("CH").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("NORMAL_NH", "NH", "[<N><H>]", ()->{
			return of("N");
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("NORMAL_NH2", "NH2", "[<N><H><2>]", ()->{
			return of("N");
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_HS", "HS", "[<H><S>]", ()->{
			return parseBranchNode("SH").get().k();
		}));
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_N_BENZO", "BnN", "[<Bn><N>]", ()->{
			return parseBranchNode("NBn").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_N_METHYL", "MeN", "[<Me><N>]", ()->{
			return parseBranchNode("NMe").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_S_METHYL", "MeS", "[<Me><S>]", ()->{
			return parseBranchNode("SMe").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_OXY_ETHYL", "EtO", "[<Et><O>]", ()->{
			return parseBranchNode("OEt").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_OCF3", "F3CO", "[<F><3><C><O>]", ()->{
			return parseBranchNode("OCF3").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_OCF3", "F3CO", "[<F><3><C><O>]", ()->{
			return parseBranchNode("OCF3").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_ETHYL_ESTER", "EtOOC", "[<Et><O><O><C>]", ()->{
			return parseBranchNode("COOEt").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_S_METHYL", "H3CS", "[<H><3><C><S>]", ()->{
			return parseBranchNode("SCH3").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_S_METHYL", "H3CS", "[<H><3><C><S>]", ()->{
			return parseBranchNode("SCH3").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_TRIPHENYL", "Ph3C", "[<P><H><3><C>]", ()->{
			return parseBranchNode("CPh3").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_TERM_ALCOHOL", "HOH2C", "[<H><O><H><2><C>]", ()->{
			return parseBranchNode("C2HOH").get().k();
		}));
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("TERM_ALCOHOL", "CH2OH", "[<C><H><2><O><H>]", ()->{
			return new BranchNode("C").addChild(new BranchNode("O"));
		}));
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("SULFOXY_METHYL", "SOCH3", "[<S><O><C><H><3>]", ()->{
			BranchNode bn = new BranchNode("S");
			bn.setCharge(1);
			bn.addChild(new BranchNode("O").setCharge(-1));
			bn.addChild(new BranchNode("C"));			
			return bn;
		}));
		
		
		
		//N_CHAIN	[<OPEN><!OPEN>+,<CLOSE><Numeric>+]	???
		parsingRules.add(new RegexTokenParsingRule("N_CHAIN","OPEN(.*)CLOSE([0-9]+)", (m)->{
					String s1=m.group(1);
					String s2=m.group(2);
					int n=Integer.parseInt(s2);
					BranchNode ps = new BranchNode("?").setPseudoNode(true);
					if(n<10){
						for(int i=0;i<n;i++){
							BranchNode bnnew=parseBranchNode(s1).get().k();
							if(bnnew!=null && bnnew.isRealNode()){
								ps.addChild(bnnew);
							}else{
								ps=null;
								
							}
						}
					}		
					return ps;
		}));
		
		parsingRules.add(new RegexTokenParsingRule("N_CHAIN_N","OPEN(.*)CLOSE([0-9]+)N", (m)->{
			String s1=m.group(1);
			String s2=m.group(2);
			return parseBranchNode("N(" + s1 + ")" +s2).get().k();
		}));
		
		
		//N_CHAIN_N	[<OPEN><!OPEN>+,<CLOSE><Numeric>+<N>]	???
		//N_CHAIN_N	[<OPEN><!OPEN>+,<CLOSE><Numeric>+<N>]	???
		
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("SO2CF3", "Tf", "[<Tf>]", ()->{
			return parseBranchNode("SO2CF3").get().k();
		}));
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("METHYL_ESTER", "MeOOC", "[<Me><O><O><C>]", ()->{
			return parseBranchNode("COOC").get().k();
		}));

		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERESE_NHCH3", "H3CHN", "[<H><3><C><H><N>]", ()->{
			return parseBranchNode("NHCH3").get().k();
		}));
		
		
		
		//REVERESE_OXY_ALKYL	[(<C><Numeric>)+<O>]	???
		parsingRules.add(new RegexTokenParsingRule("REVERESE_OXY_ALKYL","(C[0-9]+H[0-9]*)O", (m)->{
			String s1=m.group(1);
			return parseBranchNode("O" + s1).get().k();
		}));
		
		
		//N_C_COOH	[<N><H><C><H><2><CARBOXYLIC_ACID_F2>]	NHCH2COOH
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("N_C_COOH", "NHCH2COOH", "[<N><H><C><H><2><C><O><O><H>]", ()->{
			return new BranchNode("N")
					 .addChild(new BranchNode("C")
							     .addChild(parseBranchNode("COOH").get().k())
							 );
		}));
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_SO2CH3", "H3CO2S", "[<H><3><C><O><2><S>]", ()->{
			return parseBranchNode("SO2CH3").get().k();
		}));
		
		//TODO:
		//REVERSE_ALKYL_OXY_ALKYL	[<H><3><C>(<C><H><2>)*<O>(<C><H><2>)*]	???
//		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("REVERSE_ALKYL_OXY_ALKYL", "???", "[<H><3><C>(<C><H><2>)*<O>(<C><H><2>)*]", ()->{
//			return parseBranchNode("SO2CH3").get().k();
//		}));
		//REVERESE_OXY_ALKYL	[(<C><Numeric>)+<O>]	???
		parsingRules.add(new RegexTokenParsingRule("REVERSE_ALKYL_OXY_ALKYL","H3C(CH2)*O(CH2)*", (m)->{
			BranchNode bn1=new BranchNode("O");
			
			BranchNode addTo=bn1;
			String ss1=m.group(1);
			String ss2=m.group(2);
			
			if(ss1!=null && ss1.length()>0){
				BranchNode child=parseBranchNode(ss1).get().k();
				bn1.addChild(child);
				addTo=child;
			}
			BranchNode on=new BranchNode("C");
			addTo.addChild(on);
			addTo=on;
			if(ss2!=null && ss2.length()>0){
				BranchNode parent=parseBranchNode(ss2).get().k();
				parent.addChild(bn1);
				bn1=parent;
			}
			
			return bn1;
		}));		
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("N_N_LINKER", null, "[<N><N>]", ()->{
			BranchNode nn1= new BranchNode("N");
			BranchNode nn2= new BranchNode("N").setOrderToParent(2);
			nn1.addChild(nn2);
			nn1.setRightBranchNode(nn2);
			return nn1;
		}));
		
		
		
		parsingRules.add(TemplateTokenParsingRule.fromTokenShorthand("N_N_LINKER", null, "[<N><N>]", ()->{
			BranchNode nn1= new BranchNode("N");
			BranchNode nn2= new BranchNode("N").setOrderToParent(2);
			nn1.addChild(nn2);
			nn1.setRightBranchNode(nn2);
			return nn1;
		}));
		
		
		
		List<ParsingRule> atomicRules = new ArrayList<ParsingRule>();
		atomicSet.forEach(tok->{
			String SYM = tok.getTokenPreferredStyle();
			boolean termt=false;
			if(SYM.equals("F") ||  SYM.equals("Cl") || SYM.equals("Br")){
				termt=true;
			}
			boolean term=termt;
			atomicRules.add(
			new TemplateTokenParsingRule("ATOMIC_SYMBOL_" + SYM, SYM, Arrays.asList(tok), ()->{
				return new BranchNode(SYM).setTerminal(term);
			})
			);	
		});
		
		parsingRules.add(new ParsingRule(){
			
			@Override
			public String getRuleName() {
				return "ATOMIC_SYMBOL";
			}

			@Override
			public Optional<Tuple<BranchNode, String>> parse(TokenTree tt) {
				for(ParsingRule pr:atomicRules){
					Optional<Tuple<BranchNode, String>> op = pr.parse(tt);
					if(op.isPresent()){
						return op;
					}
				}
				return Optional.empty();
			}
		});
	
		
		
		//
		// getAllTokenPathsWhichAllMatch
		//
		
		parsingRules.add(new ParsingRule(){
			
			@Override
			public String getRuleName() {
				return "NUMERIC";
			}

			@Override
			public Optional<Tuple<BranchNode, String>> parse(TokenTree tt) {
				List<List<Token>> tlist = tt.getAllTokenPathsWhichAllMatch(t->{
					return numericSet.contains(t);
				});
				if(tlist.isEmpty())return Optional.empty();
				if(tlist.size()>0){
					String num=tlist.stream()
					.map(t->t.stream().map(t1->t1.getTokenPreferredStyle()).collect(Collectors.joining()))
					.findFirst()
					.orElse("");
					
					int r= Integer.parseInt(num);
					if(r>0 && r<20){
						BranchNode repNode=new BranchNode("?");
						repNode.setRepeat(r);
						return Optional.of(Tuple.of(repNode, num));
					}
				}
				return Optional.empty();
			}
		});		
		
		parsingRules.add(new RegexTokenParsingRule("METHYL_EXTRACTOR","CH[23]*(.*)", (m)->{
			String ss1=m.group(1);
			BranchNode bn = new BranchNode("C");
			
			
			BranchNode child=parseBranchNode(ss1).map(t->t.k()).orElse(null);
			if(child!=null){
				return bn.addChild(child);
			}
			return null;
		}));
		
		parsingRules.add(new RegexTokenParsingRule("N_CARBON_CHAIN_EXTRACTOR","C([1-9][0-9]*)H[1-9][0-9]*", (m)->{
			Matcher mm=m;
			
			int c=Integer.parseInt(mm.group(1));
			
			BranchNode full = new BranchNode("C");
			BranchNode parent=full;
			for(int i=1;i<c;i++){
				BranchNode nn=new BranchNode("C").thetaOffset(i%2);
				parent.addChild(nn);
				parent=nn;
			}
//			full.setRightBranchNode(parent);
			return full;
		}));
		
		parsingRules.add(new RegexTokenParsingRule("TERM_FH_EXTRACTOR","(.*F[2]*)H+", (m)->{
			Matcher mm=m;
			String cc=mm.group(1);
			return parseBranchNode(cc).map(Tuple.vmap(hp->mm.group(0))).orElse(null);
		},true));
		
//		parsingRules.add(new RegexTokenParsingRule("H_REPLACER",".*H.*", (m)->{
//			String ss1=m.group(0);
//			
//			return parseBranchNode(ss1.replace("H","")).map(t->t.k()).orElse(null);
//		}));

	}
	
	
	public static Optional<Tuple<BranchNode, String>> parseBranchNode(TokenTree tt, boolean breakUp){
//		System.out.println("subparse:");
//		tt.printAlllTrees();
//		System.out.println("-----------------");
		for(ParsingRule pr: parsingRules){
			Optional<Tuple<BranchNode, String>> op = pr.parse(tt);
			if(op.isPresent()){
//				System.out.println("Found rule:" + pr.getRuleName() + " for :");
//				tt.printAlllTrees();
//				System.out.println("#####");
				return op;
			}
		}
		if(!breakUp)return Optional.empty();
		
		List<TokenTree> ttlist = tt.asEnumeratedList()
				.stream()
				.map(tt1->Tuple.of(tt1,tt1.maxLength()).withVComparator())
				.sorted()
				.map(t->t.k())
				.collect(Collectors.toList());
		for(TokenTree mt: ttlist){
//			System.out.println("Gonna try:");
//			mt.printAlllTrees();
			//start breaking up
			int maxLength=mt.maxLength();
			if(maxLength>1){
				Optional<Tuple<BranchNode, String>> parent=null;
				for(int i=maxLength-1;i>=0;i--){
					Tuple<TokenTree,TokenTree> splitTree=mt.splitTreeAtDepth(i);
					TokenTree head=splitTree.k();
					TokenTree tail=splitTree.v();
//					System.out.println("HEAD:");
//					head.printAlllTrees();
//					System.out.println("TAIL:");
//					tail.printAlllTrees();
					Optional<Tuple<BranchNode, String>> bn1=parseBranchNode(head,false);
					if(bn1.isPresent()){
						parent=bn1;
						//System.out.println(bn1.toString());
						Optional<Tuple<BranchNode, String>> child=parseBranchNode(tail,true);
						if(child.isPresent()){
							Optional<Tuple<BranchNode, String>> bnFinal= parent.map(bnn->Tuple.of(bnn.k().addChild(child.get().k()),
									bnn.v() +child.get().v()
									));
							//System.out.println("F:" + bnFinal.toString());
							return bnFinal;
						}
						//break;
					}
				}
				if(parent!=null){
					
				}
			}	
		}
		
		
		
		
		return Optional.empty();
	}
	
	public static Optional<Tuple<BranchNode, String>> parseBranchNode(String p){
		TokenTree tt=parseTokenTree(p);
		//tt.printAlllTrees();
		//System.out.println("-------");
		return parseBranchNode(tt,true);
	}
	
	public static Optional<Tuple<BranchNode, String>> parseBranchNodeInit(String p){
//		System.out.println("Parsing:" + p);
		Optional<Tuple<BranchNode,String>> bns=parseBranchNode(p);
		if(bns.isPresent()){
//			System.out.println("Found:" +bns.get().k().toString() + " as " + bns.get().v());
		}else{
//			System.out.println("Not found:" + p);
		}
		return bns;
	}
	
	
	
	
	

	private static BranchNode interpretOCRStringAsAtom(String s, boolean tokenOnly){
//		System.out.println("wat:" + s);
		if(s.equals("1O")){
			return interpretOCRStringAsAtom("N");
		}
		
		if(s.equals("1")){
			return interpretOCRStringAsAtom("I");
		}
		if(s.contains("()")){
			return interpretOCRStringAsAtom(s.replace("()", "O"));
		}
		if(s.contains("tt")){
			return interpretOCRStringAsAtom(s.replace("tt", "H"));
		}
		if(s.contains("1t")){
			return interpretOCRStringAsAtom(s.replace("1t", "H"));
		}
		if(s.contains("t1")){
			return interpretOCRStringAsAtom(s.replace("t1", "H"));
		}
		if(s.contains("I1")){
			return interpretOCRStringAsAtom(s.replace("I1", "H"));
		}
		if(s.equals("t4")){
			return interpretOCRStringAsAtom("H");
		}
		if(s.equals("11")){
			return interpretOCRStringAsAtom("H");
		}
		
		
		if(s.contains("t1N")){
			return interpretOCRStringAsAtom(s.replace("t1N", "HN"));
		}
		
		if(s.contains("t43")){
			return interpretOCRStringAsAtom(s.replace("t43", "H3"));
		}
		if(s.equals("t42")){
			return interpretOCRStringAsAtom(s.replace("t42", "H2"));
		}
		
		if(s.equals("HCl") || s.equals("HC1") || s.equals("Hcl") || s.equals("Hc1")){
			return new BranchNode("Cl").setLinkable(false);
		}
		
		//first step is to see if this is an aliphatic chain
		if(s.matches("([cC]H[2]*)+")){
			int c=s.split("H[2]*").length;
			
			BranchNode full = new BranchNode("C");
			BranchNode parent=full;
			for(int i=1;i<c;i++){
				BranchNode nn=new BranchNode("C").thetaOffset(i%2);
				parent.addChild(nn);
				parent=nn;
			}
			full.setRightBranchNode(parent);
			return full;
		}else if(s.matches("[Cc][Oo][Oo][Cc][H][2]")){ //Ester linker
			
			BranchNode full = new BranchNode("C").addChild(new BranchNode("O").setOrderToParent(2));
			
			BranchNode OMethyl = new BranchNode("O").setOrderToParent(2);
			BranchNode meth = new BranchNode("C").setOrderToParent(1);
			
			OMethyl.addChild(meth);
			
			full.addChild(OMethyl);
			full.setRightBranchNode(meth);
			return full;
		}else if(s.matches("([cC]H[2]*)+[Cc]")){ //Aliphatic chain ending in Carbon
			int c=s.split("H[2]*").length;
			
			BranchNode full = new BranchNode("C");
			BranchNode parent=full;
			for(int i=1;i<c;i++){
				BranchNode nn=new BranchNode("C").thetaOffset(i%2);
				parent.addChild(nn);
				parent=nn;
			}
			full.setRightBranchNode(parent);
			return full;
		}else if(s.equalsIgnoreCase("NHCO")){ // N Ketone linker
			
			BranchNode full = new BranchNode("N");
			
			BranchNode carb = new BranchNode("C").thetaOffset(1).addChild(new BranchNode("O").setOrderToParent(2));
			full.addChild(carb);
			full.setRightBranchNode(carb);
			
			return full;
		}
		
		if(s.equalsIgnoreCase("CH3CONH")){
			return interpretOCRStringAsAtom("N").addChild(interpretOCRStringAsAtom("CO").addChild(new BranchNode("C")));
		}else if((  s.equalsIgnoreCase("CO2H")
				|| s.equalsIgnoreCase("CO2")
				|| s.equalsIgnoreCase("COOH")
				|| s.equalsIgnoreCase("HOOC")
				|| s.equalsIgnoreCase("HO2C")
				|| s.equalsIgnoreCase("OOC")
				|| s.equalsIgnoreCase("COO")
				
				

				|| s.equalsIgnoreCase("C02H")
				|| s.equalsIgnoreCase("C02")
				|| s.equalsIgnoreCase("C00H")
				|| s.equalsIgnoreCase("H00C")
				|| s.equalsIgnoreCase("H02C")
				|| s.equalsIgnoreCase("O00")
				|| s.equalsIgnoreCase("C00")
				)){
			BranchNode bn = new BranchNode("C");
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			BranchNode oalcohol=new BranchNode("O").setOrderToParent(1).flagForCombining();
			bn.addChild(oalcohol);
			bn.setRightBranchNode(oalcohol);
			return bn;
		}else if((       s.equalsIgnoreCase("CN")
				|| s.equalsIgnoreCase("NC"))){
			BranchNode bn = new BranchNode("C");
			bn.addChild(new BranchNode("N").setOrderToParent(3));
			//bn.addChild(new BranchNode("O").setOrderToParent(1).flagForCombining());
			return bn;
		}else if((s.equalsIgnoreCase("SO3")
				|| s.equalsIgnoreCase("O3S"))){
			BranchNode bn = new BranchNode("S");
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			bn.addChild(new BranchNode("O").setOrderToParent(1).flagForCombining());
			//bn.addChild(new BranchNode("O").setOrderToParent(1).flagForCombining());
			return bn;
		}else if((s.equalsIgnoreCase("NO2")
				|| s.equalsIgnoreCase("O2N"))){
			BranchNode bn = new BranchNode("N").setCharge(1);
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			bn.addChild(new BranchNode("O").setOrderToParent(1).setCharge(-1));
			//bn.addChild(new BranchNode("O").setOrderToParent(1).flagForCombining());
			return bn;
		}else if(s.equalsIgnoreCase("CH3")){
			BranchNode bn = new BranchNode("C").setShouldCombineLinearly(false).setTerminal(true);
			return bn;
		}else if(s.equals("Me")||s.equals("Mc")||s.equals("MC")){
			BranchNode bn = new BranchNode("C").setShouldCombineLinearly(false).setTerminal(true);
			return bn;
		}else if(s.equals("Et")){
			return new BranchNode("C").addChild(new BranchNode("C"));
		}else if(s.equals("Pr")){
			BranchNode bn = new BranchNode("C");
			bn.addChild(new BranchNode("C").addChild(new BranchNode("C")));
			return bn;
		}else if(s.equals("t")){
			BranchNode bn = new BranchNode("I");
			return bn;
		}else if(s.equals("AS") || s.equals("A8")){
			BranchNode bn = new BranchNode("As");
			return bn;
		}else if(s.equals("NCH3") || s.equals("NcH3") || s.equals("NHcH3") || s.equals("NHCH3")){
			BranchNode bn = new BranchNode("N");
			bn.addChild(new BranchNode("C"));
			return bn;
		}else if((s.equalsIgnoreCase("SO2")
				|| s.equalsIgnoreCase("S02")
				|| s.equalsIgnoreCase("O2S")
				|| s.equalsIgnoreCase("02S"))){
			BranchNode bn = new BranchNode("S");
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			//bn.addChild(new BranchNode("O").setOrderToParent(1).flagForCombining());
			return bn;
		}else if(s.equalsIgnoreCase("Ms") || s.equalsIgnoreCase("M8")){
			BranchNode bn = new BranchNode("S");
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			bn.addChild(new BranchNode("C").setOrderToParent(1));
			//bn.addChild(new BranchNode("O").setOrderToParent(1).flagForCombining());
			return bn;
		}else if(s.equalsIgnoreCase("Boc")){
			BranchNode bn = new BranchNode("C");
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			bn.addChild(new BranchNode("O").setOrderToParent(1)
					.addChild(new BranchNode("C")
							.addChild(new BranchNode("C"))
							.addChild(new BranchNode("C"))
							.addChild(new BranchNode("C"))
							));
			//bn.addChild(new BranchNode("O").setOrderToParent(1).flagForCombining());
			return bn;
		}else if(s.equalsIgnoreCase("Boc2N")){
			BranchNode b1=interpretOCRStringAsAtom("Boc");
			BranchNode b2=interpretOCRStringAsAtom("Boc");
			BranchNode bn = new BranchNode("N");
			bn.addChild(b1);
			bn.addChild(b2);
			return bn;
		}else if(s.equalsIgnoreCase("N3")){
			BranchNode bn = new BranchNode("N");
			bn.addChild(new BranchNode("N").setOrderToParent(2).setCharge(1).thetaOffset(2)
					.addChild(new BranchNode("N").setOrderToParent(2).setCharge(-1)));
			//bn.addChild(new BranchNode("O").setOrderToParent(1).flagForCombining());
			return bn;
		}else if(s.equalsIgnoreCase("CO")){
			BranchNode bn = new BranchNode("C");
			bn.addChild(new BranchNode("O").setOrderToParent(2));
			return bn;
		}else if(s.equalsIgnoreCase("OH")){
			return new BranchNode("O").setTerminal(true);
		}else if(s.equalsIgnoreCase("CH2O")){
			BranchNode bn = new BranchNode("C");
			bn.addChild(new BranchNode("O").setOrderToParent(1));
			return bn;
		}else if(s.matches("A[Cc]")){
			BranchNode bn = new BranchNode("C").setOrderToParent(1)
										   .addChild(new BranchNode("C").setOrderToParent(1))
					                       .addChild(new BranchNode("O").setOrderToParent(2));
			return bn;
		}else if(s.matches("A[Cc][Oo]")){
			BranchNode bn = new BranchNode("O");
			BranchNode ac=interpretOCRStringAsAtom("Ac");
			bn.addChild(ac);

			return bn;
		}else if(s.equalsIgnoreCase("EtO2C")){
			
			BranchNode co2=interpretOCRStringAsAtom("CO2");
			BranchNode et=interpretOCRStringAsAtom("Et");
			
			co2.addChild(et);
			return co2;
		}else if(s.equals("Bn") || s.equalsIgnoreCase("Bt1")){
			
			BranchNode carb=interpretOCRStringAsAtom("C");
			BranchNode ben=interpretOCRStringAsAtom("Ph");
			
			
			carb.addChild(ben);
			return carb;
		}else if(s.equals("BN")){
			return interpretOCRStringAsAtom("HN");
		}else if(s.equalsIgnoreCase("CBZ") || s.equalsIgnoreCase("C6Z")){
			BranchNode carb=interpretOCRStringAsAtom("C").addChild(new BranchNode("O").setOrderToParent(2));
			BranchNode ox=new BranchNode("O").setOrderToParent(1);
			carb.addChild(ox);
			BranchNode ben=interpretOCRStringAsAtom("Bn");
			ox.addChild(ben);
			
			return carb;
		}else if(s.equals("t-Bu") || s.equals("tBu") || s.equals("t-Bo") || s.equals("tBo")){
			BranchNode bn = new BranchNode("C");
			bn.addChild(new BranchNode("C"));
			bn.addChild(new BranchNode("C"));
			bn.addChild(new BranchNode("C"));
			return bn;
		}else if(s.equals("n-Bu") || s.equals("nBu") || s.equals("n-Bo") || s.equals("nBo")){
			return interpretOCRStringAsAtom("C4H9");
		}else if(s.equalsIgnoreCase("CO2Cys") || s.equalsIgnoreCase("CO2Cy8")){
			
			BranchNode co2=interpretOCRStringAsAtom("CO2");
			BranchNode cys=interpretOCRStringAsAtom("Cys");
			
			BranchNode combine=cys.getNodeForLinearCombine();
			
			co2.addChild(cys);
			co2.setCombiningNode(combine);
			return co2;
		}else if(s.equalsIgnoreCase("Cys")){
			BranchNode bn = new BranchNode("N");
			BranchNode carbonyl=new BranchNode("C").addChild(new BranchNode("O").setOrderToParent(2));
			bn.addChild(new BranchNode("C").setOrderToParent(1)
										   .flagForCombining()
					                       .addChild(new BranchNode("C").thetaOffset(1).setWedgeToParent(1).addChild(new BranchNode("S")))
					                       .addChild(carbonyl)
					);
			
			bn.setCombiningNode(carbonyl);
			return bn;
		}else if(s.matches("[Pp]h") || s.equals("Pb") || s.equals("pb")){
			BranchNode bn = new BranchNode("C").thetaOffset(1);
			
			bn.addChild(new BranchNode("C").setOrderToParent(2)
					    			       .addChild(new BranchNode("C").setOrderToParent(1)
					    			    		   .addChild(new BranchNode("C").setOrderToParent(2)
					    			    				   .addChild(new BranchNode("C").setOrderToParent(1).addChild(new BranchNode("C").setOrderToParent(2).addRing(bn, 1)))
					    			    		   )
					    			    	)
					);
			
			
			return bn;
		}else if(s.equals("[p-tol]")){
			BranchNode bn = new BranchNode("C").thetaOffset(1);
			
			bn.addChild(new BranchNode("C").setOrderToParent(2) //ortho
					    			       .addChild(new BranchNode("C").setOrderToParent(1) //meta
					    			    		   .addChild(new BranchNode("C").setOrderToParent(2) //para
					    			    				   .addChild(new BranchNode("C").setOrderToParent(1) //meta
					    			    						   .addChild(new BranchNode("C").setOrderToParent(2).addRing(bn, 1))) //ortho
					    			    				   .addChild(new BranchNode("C"))
					    			    		   )
					    			    	)
					);
			
			
			return bn;
		}else if(s.equals("Ts")){
			return interpretOCRStringAsAtom("SO2[p-tol]");
		}else if(s.equalsIgnoreCase("PMBN")){
			BranchNode bno = new BranchNode("N").thetaOffset(1);
			BranchNode ml = new BranchNode("C");
			bno.addChild(ml);
			
			BranchNode bn = new BranchNode("C").thetaOffset(1);
			
			bn.addChild(new BranchNode("C").setOrderToParent(2)
					    			       .addChild(new BranchNode("C").setOrderToParent(1)
					    			    		   .addChild(new BranchNode("C").setOrderToParent(2)
					    			    				   .addChild(new BranchNode("C").setOrderToParent(1).addChild(new BranchNode("C").setOrderToParent(2).addRing(bn, 1)))
					    			    				   .addChild(new BranchNode("O").addChild(new BranchNode("C")))
					    			    		   )
					    			    	)
					);
			ml.addChild(bn);
			
			
			return bno;
		}else if(s.equalsIgnoreCase("CH3O") || s.equalsIgnoreCase("H3CO")){
			return interpretOCRStringAsAtom("OCH3");
		}else if(s.equalsIgnoreCase("OCH3")){
			return interpretOCRStringAsAtom("O").addChild(interpretOCRStringAsAtom("C"));
		}else if(s.equalsIgnoreCase("EtHN")){
			return interpretOCRStringAsAtom("NHEt");
		}else if(s.equalsIgnoreCase("OHC")){
			return interpretOCRStringAsAtom("COH");
		}else if(s.matches("II*O")){ // usually from dash bonds
			return interpretOCRStringAsAtom("O");
		}else if(s.matches("N[lI][lI]*")){ // usually from dash bonds
			return interpretOCRStringAsAtom("N");
		}else if(s.equalsIgnoreCase("S+")){
			BranchNode bn = new BranchNode("S");
			bn.setCharge(1);
			return bn;
		}else if(s.equalsIgnoreCase("N+") || s.equalsIgnoreCase("+N")){
			BranchNode bn = new BranchNode("N");
			bn.setCharge(1);
			return bn;
		}else if(s.equalsIgnoreCase("BnO")){
			return interpretOCRStringAsAtom("OBn");
		}else if(s.equalsIgnoreCase("MeO") || s.equalsIgnoreCase("McO")){
			return interpretOCRStringAsAtom("OMe");
		}else if(s.equals("BnN") || s.equals("Bt1N")){
			return interpretOCRStringAsAtom("NBn");
		}else if(s.equals("MeN") || s.equals("McN")){
			return interpretOCRStringAsAtom("NMe");
		}else if(s.equals("MeS") || s.equals("McS")){
			return interpretOCRStringAsAtom("SMe");
		}else if(s.equalsIgnoreCase("EtO") ){
			return interpretOCRStringAsAtom("OEt");
		}else if(s.equalsIgnoreCase("F3CO") ){
			return interpretOCRStringAsAtom("OCF3");
		}else if(s.equalsIgnoreCase("EtOOC") ){
			return interpretOCRStringAsAtom("COOEt");
		}else if(s.equalsIgnoreCase("H3CS") ){
			return interpretOCRStringAsAtom("SCH3");
		}else if(s.equalsIgnoreCase("PH3C") ){
			return interpretOCRStringAsAtom("CPH3");
		}else if(s.equalsIgnoreCase("HOH2C") ){
			return interpretOCRStringAsAtom("CH2OH");
		}else if(s.equalsIgnoreCase("CH2OH") ){
			return interpretOCRStringAsAtom("C").addChild(interpretOCRStringAsAtom("O"));
		}
//		else if(s.equalsIgnoreCase("CH3CONH") ){
//			return interpretOCRStringAsAtom("NHOCCH3");
//		}
		else if(s.equalsIgnoreCase("SOCH3") ){
			BranchNode bn = new BranchNode("S");
			bn.setCharge(1);
			bn.addChild(new BranchNode("O").setCharge(-1));
			bn.addChild(interpretOCRStringAsAtom("CH3"));			
			return bn;
		}else if(s.matches("[(].*[)][0-9][0-9]*")){
			Pattern p = Pattern.compile("[(](.*)[)]([0-9][0-9]*)");
			Matcher m=p.matcher(s);
			if(m.find()){
				String s1=m.group(1);
				String s2=m.group(2);
				
				int n=Integer.parseInt(s2);
				
				BranchNode ps = new BranchNode("?").setPseudoNode(true);
				
				
				if(n<10){
					for(int i=0;i<n;i++){
						BranchNode bnnew=interpretOCRStringAsAtom(s1);
						if(bnnew!=null && bnnew.isRealNode()){
							ps.addChild(bnnew);
						}else{
							ps=null;
							break;
						}
					}
				}
				if(ps!=null)return ps;
				
				
				//System.out.println("Got:" + s1 + "," + s2 + " from " + s);
			}
			 
		}else if(s.matches("[(].*[)][0-9][0-9]*N")){
			Pattern p = Pattern.compile("[(](.*)[)]([0-9][0-9]*)N");
			Matcher m=p.matcher(s);
			if(m.find()){
				return interpretOCRStringAsAtom("N(" + m.group(1) + ")" + m.group(2));
			}
		}else if(s.equals("Tf") || s.equals("Tr")){
			return interpretOCRStringAsAtom("SO2CF3");
		}else if(s.matches("M[ecC][O0o][O0o][Cc]")){
			return interpretOCRStringAsAtom("COOC");
		}else if(s.equalsIgnoreCase("H3CHN")){
			return interpretOCRStringAsAtom("N").addChild(interpretOCRStringAsAtom("CH3"));
		}else if(s.matches("[Cc][0-9][0-9]*H[0-9][0-9]*O")){
			Pattern p = Pattern.compile("([Cc][0-9][0-9]*H[0-9][0-9]*)O");
			Matcher m=p.matcher(s);
			if(m.find()){
				return interpretOCRStringAsAtom("O" + m.group(1));
			}
		}else if(s.matches("NH[Cc]H2[cC][Oo0][Oo0]H")){
			return interpretOCRStringAsAtom("N")
					 .addChild(interpretOCRStringAsAtom("C")
							     .addChild(interpretOCRStringAsAtom("COOH"))
							 );
		}else if(s.matches("H3[Cc][Oo0]2S")){
			return interpretOCRStringAsAtom("SO2CH3");
		}else if(s.matches("H3[Cc]([Cc]H2)*[Oo0]([Cc]H2)*")){
			Pattern p = Pattern.compile("H3[Cc]([Cc]H2)*[Oo0]([Cc]H2)*");
			Matcher m=p.matcher(s);
			if(m.find()){
				
				BranchNode bn1=new BranchNode("O");
				
				BranchNode addTo=bn1;
				String ss1=m.group(1);
				String ss2=m.group(2);
				
				if(ss1!=null && ss1.length()>0){
					BranchNode child=interpretOCRStringAsAtom(ss1);
					bn1.addChild(child);
					addTo=child;
				}
				BranchNode on=new BranchNode("C");
				addTo.addChild(on);
				addTo=on;
				
				if(ss2!=null && ss2.length()>0){
					BranchNode parent=interpretOCRStringAsAtom(ss2);
					parent.addChild(bn1);
					bn1=parent;
				}
				
				return bn1;
			}
			
		}else if(s.matches("[Cc]H2([Cc]H2)*[Oo0]([Cc]H2)*[Cc]H3")){
			Pattern p = Pattern.compile("[Cc]H2([Cc]H2)*[Oo0]([Cc]H2)*[Cc]H3");
			Matcher m=p.matcher(s);
			if(m.find()){
				
				BranchNode bn1=new BranchNode("C");
				
				BranchNode addTo=bn1;
				String ss1=m.group(1);
				String ss2=m.group(2);
				
				if(ss1!=null && ss1.length()>0){
					BranchNode child=interpretOCRStringAsAtom(ss1);
					bn1.addChild(child);
					addTo=child;
				}
				BranchNode on=new BranchNode("O");
				addTo.addChild(on);
				addTo=on;
				
				if(ss2!=null && ss2.length()>0){
					BranchNode child=interpretOCRStringAsAtom(ss2);
					addTo.addChild(child);
					addTo=child;
				}
				addTo.addChild(new BranchNode("C"));
				
				return bn1;
			}
			
		}else if(s.equals("NN")){
			BranchNode nn1= new BranchNode("N");
			BranchNode nn2= new BranchNode("N").setOrderToParent(2);
			nn1.addChild(nn2);
			nn1.setRightBranchNode(nn2);
			return nn1;
		}else if(s.equals("Sl") || s.equals("SI")){
			return new BranchNode("Si");
		}else if(s.equals("Sr") || s.equals("sr")|| s.equals("8r") || s.equals("BT")){
			return new BranchNode("Br");
		}
		
		
		
		//TODO:
		
		// Et
		// Ph
		
		if(accept.contains(s)){
			BranchNode bn=new BranchNode(s);
			if(s.equals("F") || s.equals("Br") || s.equals("Cl")){
				bn=bn.setTerminal(true);
			}
			return bn;
		}else if(accept.contains(s.toUpperCase()) && !s.equals("b") && !s.equals("n")){
			
			
			return new BranchNode(s.toUpperCase());
		}
		if(s.contains("Ct")){
			return interpretOCRStringAsAtom(s.replaceAll("C[t]", "Cl"),tokenOnly);
		}
		if(s.equals("C)") || s.equals("c)")){
			return interpretOCRStringAsAtom("Cl",tokenOnly);
		}
		if(s.contains("Htt")){
			return interpretOCRStringAsAtom(s.replace("Htt", "H11"),tokenOnly);
		}
		if(s.contains("Ht")){
			return interpretOCRStringAsAtom(s.replace("Ht", "H1"),tokenOnly);
		}
		if(s.contains("H")){
			return interpretOCRStringAsAtom(s.replaceAll("H[1-9][0-9]*", "").replace("H", ""),tokenOnly);
		}
		if(s.contains("I")){
			return interpretOCRStringAsAtom(s.replace("I", "l"),tokenOnly);
		}

		if(s.contains("t1")){
			
			return interpretOCRStringAsAtom(s.replace("t1", "n"),tokenOnly);
		}
		if(s.contains("1")){
			return interpretOCRStringAsAtom(s.replace("1", "l"),tokenOnly);
		}
		if(s.contains("c")){
			return interpretOCRStringAsAtom(s.replace("c", "C"),tokenOnly);
		}
		if(s.equals("0")){
			return interpretOCRStringAsAtom(s.replace("0", "O"),tokenOnly);
		}
		
		
		try{
			int r=Integer.parseInt(s);
			if(r>0 && r<20){
				BranchNode repNode=new BranchNode("?");
				repNode.setRepeat(r);
				return repNode;
			}
		}catch(Exception e){
			
		}
		if(tokenOnly){
			return null;
		}
		
		//start breaking up
		if(s.length()>1){
			BranchNode parent=null;
			String rem=null;
			for(int i=s.length()-1;i>=1;i--){
				String substr=s.substring(0, i);
				BranchNode bn1=interpretOCRStringAsAtom(substr,true);
				if(bn1!=null){
					parent=bn1;
					//System.out.println(bn1.toString());
					rem=s.substring(i);
					BranchNode child=interpretOCRStringAsAtom(rem);
					if(child!=null){
						BranchNode bnFinal= parent.addChild(child);
						//System.out.println("F:" + bnFinal.toString());
						return bnFinal;
					}
					//break;
				}
			}
			if(parent!=null){
				
			}
		}
		
		return null;		
	}
	
	public BranchNode removeHydrogens(){
		BranchNode _this=this;
		this.forEachBranchNode((p,c)->{
			if(p!=null){
				if(c.getSymbol().equals("H")){
					p.children.remove(c);
				}
			}
		});
		
		return this;
	}
	public BranchNode setWedgeToParent(int i) {
		this.wedgeToParent=i;
		return this;
	}
	public int getWedgeType(){
		return this.wedgeToParent;
	}

	private static BranchNode interpretOCRStringAsAtom(String s){
		try{
			Optional<BranchNode> bn= _cache.get(s);
			if(bn==null){
				
				
				//old way
				if(false){
					bn=Optional.ofNullable(interpretOCRStringAsAtom(s,false));
					_cache.put(s, bn);
				}else{
				//new way
					bn=parseBranchNodeInit(s).map(b->b.k()
													.removeHydrogens()
													.setAlias(b.v())
													);
					_cache.put(s, bn);
				}
			}
			return bn.map(b->b.cloneNode()).orElse(null);
		}catch(Exception e){
			return null;
		}
	}
	
	public static BranchNode interpretOCRStringAsAtom2(String s){
		try{
			return Optional.ofNullable(interpretOCRStringAsAtom(s)).map(b->b.removeHydrogens()).orElse(null);
		}catch(Exception e){
			e.printStackTrace();
			return null;
		}
	}
	
	

	public boolean isTerminal() {
		return isTerminal;
	}
	public BranchNode setTerminal(boolean t){
		this.isTerminal=t;
		return this;
	}

	public Point2D getSuggestedPoint() {
		return suggestedPoint;
	}

	public void setSuggestedPoint(Point2D suggestedPoint) {
		this.suggestedPoint = suggestedPoint;
	}
	
}