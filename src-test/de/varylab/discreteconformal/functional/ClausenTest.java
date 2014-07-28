package de.varylab.discreteconformal.functional;

import org.junit.Assert;
import org.junit.Test;

import de.jtem.mfc.field.Complex;

public class ClausenTest {

	@Test
	public void testImLi2() throws Exception {
		// compare with results from mathematica
		Complex z = new Complex(0.5, 0.5);
		double r = Clausen.ImLi2(z);
		Assert.assertEquals(0.643767332889268748742017, r, 1E-15);
		
		z = new Complex(1.234, -3.554);
		r = Clausen.ImLi2(z);
		Assert.assertEquals(-2.784709023190200241929031, r, 1E-15);	
		
		z = new Complex(-4.254, 2.965);
		r = Clausen.ImLi2(z);
		Assert.assertEquals(1.104092479314665907914598, r, 1E-15);
		
		z = new Complex(-2.264, -1.05);
		r = Clausen.ImLi2(z);
		Assert.assertEquals(-0.540683148040955780529547, r, 1E-15);
		
		z = new Complex(-4, 2);
		r = Clausen.ImLi2(z);
		Assert.assertEquals(0.7855694340809751306351005, r, 1E-15);
	}
	
}
