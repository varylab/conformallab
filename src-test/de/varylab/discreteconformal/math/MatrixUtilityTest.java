package de.varylab.discreteconformal.math;

import org.junit.Assert;
import org.junit.Test;

import de.jreality.math.Rn;

public class MatrixUtilityTest {

	@Test
	public void testMakeMappingMatrix() throws Exception {
		double[] s1 = {2,0,2,1};  
		double[] s2 = {1,1,0,0};  
		double[] s3 = {1,0,8,0};  
		double[] s4 = {3,4,0,1};  
		
		double[] t1 = {2,0,1,0};  
		double[] t2 = {0,3,0,2};  
		double[] t3 = {1,0,4,0};  
		double[] t4 = {0,3,0,5};
		
		double[][] S = {s1,s2,s3,s4};
		double[][] T = {t1,t2,t3,t4};
		double[] R = MatrixUtility.makeMappingMatrix(S, T);
		Assert.assertArrayEquals(t1, Rn.matrixTimesVector(null, R, s1), 1E-15);
		Assert.assertArrayEquals(t2, Rn.matrixTimesVector(null, R, s2), 1E-15);
		Assert.assertArrayEquals(t3, Rn.matrixTimesVector(null, R, s3), 1E-15);
		Assert.assertArrayEquals(t4, Rn.matrixTimesVector(null, R, s4), 1E-15);
	}
	
}
