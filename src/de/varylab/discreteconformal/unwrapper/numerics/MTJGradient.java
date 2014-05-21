package de.varylab.discreteconformal.unwrapper.numerics;

import de.jtem.halfedgetools.functional.Gradient;
import no.uib.cipr.matrix.Vector;

public class MTJGradient extends MTJDomain implements Gradient {

	public MTJGradient(Vector u) {
		super(u);
	}

	@Override
	public void add(double coeff, Gradient g) {
		super.add(coeff, (MTJGradient)g);
	}

}