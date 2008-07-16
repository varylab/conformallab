package de.varylab.discreteconformal.math;

import de.varylab.discreteconformal.heds.CHDS;

public class FunctionalVertexImpl extends Functional {

	private CHDS
		S = null;
	
	public FunctionalVertexImpl(CHDS S) {
		this.S = S;
	}
	
	
	@Override
	void evaluate(double[] u, double[] f, double[] grad, double[][] hess) {
		
	}

}
