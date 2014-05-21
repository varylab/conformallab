package de.varylab.discreteconformal.unwrapper.numerics;

import static de.jtem.jpetsc.MatStructure.SAME_NONZERO_PATTERN;
import de.jtem.halfedgetools.functional.Hessian;
import de.jtem.jpetsc.InsertMode;
import de.jtem.jpetsc.Mat;

public class TaoHessian implements Hessian {
	
	private Mat
		H = null;
	
	public TaoHessian(Mat H) {
		this.H = H;
	}

	@Override
	public void add(int i, int j, double value) {
		H.add(i, j, value);
	}
	@Override
	public void add(double alpha, Hessian h) {
		H.aXPY(alpha, ((TaoHessian)h).H, SAME_NONZERO_PATTERN);
	}
	@Override
	public void setZero() {
		H.zeroEntries();
	}

	@Override
	public void set(int i, int j, double value) {
		H.setValue(i, j, value, InsertMode.INSERT_VALUES);
	}
	
	@Override
	public double get(int i, int j) {
		return H.getValue(i, j);
	}
	
}