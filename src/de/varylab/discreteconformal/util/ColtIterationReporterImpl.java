package de.varylab.discreteconformal.util;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.algo.solver.DoubleIterationReporter;

public class ColtIterationReporterImpl implements DoubleIterationReporter{

	@Override
	public void monitor(double arg0, DoubleMatrix1D arg1, int arg2) {
		monitor(arg0, arg2);
	}

	@Override
	public void monitor(double arg0, int arg1) {
		if (arg1 % 1000 == 0)
			System.err.println("iteration = " + arg1 + ", value = " + arg0);
	}
}
