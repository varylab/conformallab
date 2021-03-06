package de.varylab.discreteconformal.unwrapper.numerics;

import no.uib.cipr.matrix.Matrix;
import de.jtem.halfedgetools.functional.Hessian;

public class MTJHessian implements Hessian {
	
	private Matrix
		H = null;
	
	public MTJHessian(Matrix H) {
		this.H = H;
	}

	@Override
	public void add(int i, int j, double value) {
		H.add(i, j, value);
	}
	
	@Override
	public void add(double coeff, Hessian h) {
		H.add(coeff, ((MTJHessian)h).H);
	}

	@Override
	public void set(int i, int j, double value) {
		H.set(i, j, value);
	}
	
	@Override
	public void setZero() {
		H.zero();
	}
	
	@Override
	public double get(int i, int j) {
		return H.get(i, j);
	}
	
}