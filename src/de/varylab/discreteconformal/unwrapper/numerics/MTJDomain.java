package de.varylab.discreteconformal.unwrapper.numerics;

import no.uib.cipr.matrix.Vector;
import de.jtem.halfedgetools.functional.DomainValue;

public class MTJDomain implements DomainValue {

	private Vector
		u = null;
	
	public MTJDomain(Vector u) {
		this.u = u;
	}
	
	@Override
	public void add(int i, double value) {
		u.add(i, value);
	}

	@Override
	public void set(int i, double value) {
		u.set(i, value);
	}

	@Override
	public void setZero() {
		u.zero();
	}
	
	@Override
	public double get(int i) {
		return u.get(i);
	}
	
}