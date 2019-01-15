package tripod.molvec.algo;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import org.junit.Test;

import tripod.molvec.util.CompareUtil;
public class UtilTest {

	
	@Test
	public void maxStreamUsingNaturalOrderComparatorShouldGiveMaximumValue(){
		int[] is = new int[]{1,3,6,2,31,602,-10,8};
		
		int mini = -10;
		int maxi = 602;
		
		int maxg=Arrays.stream(is)
					   .mapToObj(i->i)
		               .max(CompareUtil.naturalOrder()).orElse(null);
		
		int ming=Arrays.stream(is)
				   .mapToObj(i->i)
	               .min(CompareUtil.naturalOrder()).orElse(null);
		
		assertEquals(mini,ming);
		assertEquals(maxi,maxg);
		
		
	}
}
