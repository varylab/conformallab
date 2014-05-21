package de.varylab.discreteconformal.unwrapper.numerics;

import de.jtem.halfedgetools.functional.DomainValue;
import de.jtem.jpetsc.InsertMode;
import de.jtem.jpetsc.Vec;

public class TaoDomain implements DomainValue {

	private Vec
		u = null;
	
	public TaoDomain(Vec u) {
		this.u = u;
	}

	@Override
	public void add(int i, double value) {
		u.add(i, value);
	}
	
	@Override
	public void add(double coeff, DomainValue x) {
		u.aXPY(coeff, ((TaoDomain)x).u);
	}

	@Override
	public void set(int i, double value) {
		u.setValue(i, value, InsertMode.INSERT_VALUES);
	}

	@Override
	public void setZero() {
		u.zeroEntries();
	}

	@Override
	public double get(int i) {
		return u.getValue(i);
	}
	
}