package de.varylab.discreteconformal.math.optimization;


import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;


/**
 * Implementers must provide methods for evaluating the value, gradient
 * and hessian. The mtj vector classes are used.
 * <p>
 * Copyright 2005 <a href="http://www.sechel.de">Stefan Sechelmann</a>
 * <a href="http://www.math.tu-berlin.de/geometrie">TU-Berlin</a> 
 * @author Stefan Sechelmann
 * @see math.optimization.newton.NewtonOptimizer
 * @see mt.Vector
 * @see mt.Matrix
 */
public interface Optimizable {

	/**
	 * Evaluates this Optimizable an the given x and returns the gradient
	 * and hessian
	 * @param x
	 * @param gradient the gradient return
	 * @param hessian the hessian return
	 * @return the value at x
	 */
	public Double evaluate(Vector x, Vector gradient, Matrix hessian);
	
	public Double evaluate(Vector x, Vector gradient);
	
	public Double evaluate(Vector x, Matrix hessian);
	
	public Double evaluate(Vector x);
	
	public Integer getDomainDimension();
	
}
