package gov.nih.ncats.molvec.pub;

import java.io.IOException;
import java.util.Arrays;

import gov.nih.ncats.molvec.image.Bitmap;

/**
 * Aides in creating figures for explaining thinning algorithm
 * @author tyler
 *
 */
public class ThinningAlgorithmFigureAide {

	/*
	 *   
       Step 1:  All contour pixels satisfying the following conditions are
       marked for deletion.
       
       (a) 2 <= N(p1) <= 6
       (b) S(p1) = 1
       (c) p2 * p4 * p6 = 0
       (d) p4 * p6 * p8 = 0
       
       where N(p1) is the number of 8-connected foreground pixels around p1.  
       S(p1) is the number of 0-1 transitions from p2, p3, p4, ..., p8, p9, p2
       according to the following labeling scheme.
       p9  p2  p3
       p8  p1  p4
       p7  p6  p5
       After all contour pixels have been processed, those that were marked for
       deletion are deleted from the image.  Next, step 2 is applied to the 
       image.
       
       Step 2:  This step is almost identical to step 1.  The only difference 
       here is that conditions (c) and (d) are changed to
       
       (c') p2 * p4 * p8 = 0
       (d') p2 * p6 * p8 = 0
       
	 */
	
	public static void main(String[] args) throws IOException{
		
		
		int c1=0;
		
//		Bitmap bm = new Bitmap(3*8, 28*3);
		//bm=bm.invert();
		
		//Both agree on
		// 
		// ##.
		// .X. (all 4 rotations)
		// ...
		//
		// ###
		// .X. (all 4 rotations)
		// ...
		//
		// ##.
		// #X. (all 4 rotations)
		// ...
		//
		// ###
		// #X. (all 4 rotations)
		// ...
		//
		// ###
		// .X# (all 4 rotations)
		// ...
		//
		// .##
		// .X. (all 4 rotations)
		// ...
		//
		// ###
		// .X# (all 4 rotations)
		// ..#
		
		
		
		
		for(int i=0;i<256;i++){
			int c=Integer.bitCount(i);
//			if(c>=2 && c<=6){
				int tc = countTransition(i);
//				if(tc!=1)continue;
				boolean rule1New= rule1New(i);
				boolean rule2New= rule2New(i);
				boolean rule1= rule1(i);
				boolean rule2= rule2(i);
				String criteria = "";
				if(rule1 && rule2)criteria = "both rules";
				if(rule1 && !rule2)criteria = "rule 1 only";
				if(!rule1 && rule2)criteria = "rule 2 only";
				if(rule1New && ! rule2){
					
//					drawShape((byte)i,bm, (c1/7)*3, (c1%7)*3);
					printShape((byte)i,i+":" + tc + ":" + criteria + " " + (c1++));
					System.out.println();
				}
//			}
		}
//		ImageIO.write(bm.createBufferedImage(), "png",
//				new File("tmp.png"));
		
	}
	
//	public static int countTransition(byte b){
//		int tc = 0;
//		int[] on = new int[8];
//		for(int i=0;i<8;i++){
//			on[i] = ((b & (1<<i))!=0)?1:0;
//		}
//		
//		for(int i=0;i<8;i++){
//			if(on[(i+1)%8] >on[i]){
//				tc++;
//			}
//		}
//		return tc;
//		
//	}
//	
	public static boolean rule1(int b){
		int c=Integer.bitCount(b);
		if(c<2 || c>6)return false;
		if(countTransition(b)!=1)return false;
		
		int[] on = new int[8];
		for(int i=0;i<8;i++){
			on[i] = ((b & (1<<i))!=0)?1:0;
		}
		
		return ((on[0] & on[2] & on[4]) == 0) && ((on[6] & on[2] & on[4]) == 0);
		
	}
	

	public static boolean rule2(int b){
		int c=Integer.bitCount(b);
		if(c<2 || c>6)return false;
		if(countTransition(b)!=1)return false;
		
		int[] on = new int[8];
		for(int i=0;i<8;i++){
			on[i] = ((b & (1<<i))!=0)?1:0;
		}
		
		return ((on[0] & on[2] & on[6]) == 0) && ((on[0] & on[4] & on[6]) == 0);
		
	}
//	
	public static void drawShape(byte b, Bitmap bm, int yoff, int xoff){
		int[] on = new int[8];
		
		for(int i=0;i<8;i++){
			on[i] = ((b & (1<<i))!=0)?1:0;
		}
		bm.set(xoff+ 1, yoff+1, true);
		bm.set(xoff+ 0, yoff, on[0]!=0);
		bm.set(xoff+ 1, yoff, on[1]!=0);
		bm.set(xoff+ 2, yoff, on[2]!=0);
		bm.set(xoff+ 2, yoff+1, on[3]!=0);
		bm.set(xoff+ 2, yoff+2, on[4]!=0);
		bm.set(xoff+ 1, yoff+2, on[5]!=0);
		bm.set(xoff+ 0, yoff+2, on[6]!=0);
		bm.set(xoff+ 0, yoff+1, on[7]!=0);
		
	}
	
	public static void printShape(byte b , String and){
		
		
		
		int[] on = new int[8];
		
		for(int i=0;i<8;i++){
			on[i] = ((b & (1<<i))!=0)?1:0;
		}
		
		String[] ob=Arrays.stream(on)
		      .mapToObj(o->((o!=0)?"#":"."))
		      .toArray(j->new String[j]);
		System.out.println(ob[7] + ob[0] + ob[1]);
		System.out.println(ob[6] + "X" + ob[2] + " " + and);
		System.out.println(ob[5] + ob[4] + ob[3]);
		
	}
	

   	
	private static boolean rule1New(int b){
		int tot = Integer.bitCount(b);
		if(tot>=2 && tot<=6){
			boolean[] on = new boolean[8];
			for(int i=0;i<8;i++){
				on[i] = ((b & (1<<i))!=0)?true:false;
			}
			boolean ex1 = !on[0] && !on[1] && !on[2] && !on[5] && on[4] && on[6];
			boolean ex2 = !on[2] && !on[3] && !on[4] && !on[7] && on[6] && on[0];
			
			
			
			if(countTransition(b)==1 || ex1 || ex2){				
				
				return ((!on[2] || !on[0] || !on[6])) && ((!on[4] || !on[0] || !on[6]));
				
			}
		}
		return false;
		
	}
	
	private static boolean rule2New(int b){
		int tot = Integer.bitCount(b);
		if(tot>=2 && tot<=6){
			boolean[] on = new boolean[8];
			for(int i=0;i<8;i++){
				on[i] = ((b & (1<<i))!=0)?true:false;
			}
		     
//          if (nb > 1 && nb < 7 
//              && (ap == 1 || ((1-parity)*cp + parity*dp) == 1)) {
//              int ep = (thin.p2(x, y) + thin.p4(x, y))
//                  *thin.p0(x, y) * thin.p6(x, y);
//              int fp = (thin.p6(x, y) + thin.p0(x, y))
//                  *thin.p4(x, y) * thin.p2(x, y);
			boolean ex1 = !on[1] && !on[4] && !on[5] && !on[6] && on[0] && on[2];
			boolean ex2 = !on[0] && !on[3] && !on[6] && !on[7] && on[2] && on[4];
			
			
			
			if(countTransition(b)==1 || ex1 || ex2){				
				return ((!on[2] || !on[4] || !on[6])) && ((!on[4] || !on[0] || !on[2]));
			}
		}
		return false;
		
	}
	
	public static int countTransition(int b){
		int tc = 0;
		int[] on = new int[8];
		for(int i=0;i<8;i++){
			on[i] = ((b & (1<<i))!=0)?1:0;
		}
		
		for(int i=0;i<8;i++){
			if(on[(i+1)%8] >on[i]){
				tc++;
			}
		}
		return tc;
		
	}
    
}
