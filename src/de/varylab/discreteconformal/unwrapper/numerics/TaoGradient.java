package de.varylab.discreteconformal.unwrapper.numerics;

import de.jtem.halfedgetools.functional.Gradient;
import de.jtem.jpetsc.Vec;

public class TaoGradient extends TaoDomain implements Gradient {

	public TaoGradient(Vec u) {
		super(u);
	}

	@Override
	public void add(double coeff, Gradient g) {
		super.add(coeff, (TaoGradient)g);
	}
	
}