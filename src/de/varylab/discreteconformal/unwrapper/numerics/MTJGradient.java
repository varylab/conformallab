package de.varylab.discreteconformal.unwrapper.numerics;

import no.uib.cipr.matrix.Vector;
import de.jtem.halfedgetools.functional.Gradient;

public class MTJGradient implements Gradient {

	private Vector
		G = null;
	
	public MTJGradient(Vector G) {
		this.G = G;
	}
	
	@Override
	public void add(int i, double value) {
		G.add(i, value);
	}

	@Override
	public void set(int i, double value) {
		G.set(i, value);
	}
	
	@Override
	public void setZero() {
		G.zero();
	}
	
	@Override
	public double get(int i) {
		return G.get(i);
	}
	
}