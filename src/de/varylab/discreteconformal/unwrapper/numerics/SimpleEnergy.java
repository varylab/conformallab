/**
 * 
 */
package de.varylab.discreteconformal.unwrapper.numerics;

import de.jtem.halfedgetools.functional.Energy;

public class SimpleEnergy implements Energy {

	public double 
		E = 0.0;
	
	@Override
	public double get() {
		return E;
	}
	
	@Override
	public void add(double E) {
		this.E += E;
	}

	@Override
	public void set(double E) {
		this.E = E;
	}

	@Override
	public void setZero() {
		E = 0.0;
	}
	
}