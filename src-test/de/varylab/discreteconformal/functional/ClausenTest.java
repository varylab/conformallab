package de.varylab.discreteconformal.functional;

import org.junit.Test;

import de.jtem.mfc.field.Complex;

public class ClausenTest {

	@Test
	public void testImLi2() throws Exception {
		Complex z = new Complex(0.5, 0.5);
		double r = Clausen.ImLi2(z);
		System.out.println(r);
	}
	
}
