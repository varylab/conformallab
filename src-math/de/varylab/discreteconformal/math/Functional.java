package de.varylab.discreteconformal.math;

public abstract class Functional {

	/**
	 * Evaluation of the functional
	 * @param x the point to evaluate at
	 * @param f f.length = 1 the value of f(x) is null if only the gradient is needed
	 * @param grad the gradient of f at x if null if only the value is needed
	 * @param hess the hessian or null if not needed
	 */
	abstract void evaluate(double[] x, double[] f, double[] grad, double[][] hess);
	
	
}
