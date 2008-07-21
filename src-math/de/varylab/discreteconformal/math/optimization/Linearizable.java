package de.varylab.discreteconformal.math.optimization;


import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;


/**
 * Implementers must provide methods for evaluating the value and gradient.
 * <p>
 * Copyright 2005 <a href="http://www.sechel.de">Stefan Sechelmann</a>
 * <a href="http://www.math.tu-berlin.de/geometrie">TU-Berlin</a> 
 * @author Stefan Sechelmann
 * @see math.optimization.newton.NewtonSolver
 */
public interface Linearizable {

	public void evaluate(Vector x, Vector fx, Vector offset) throws FunctionNotDefinedException;
	
	public void evaluate(Vector x, Vector fx, Vector offset, Matrix jacobian) throws FunctionNotDefinedException;
	
	public void evaluate(Vector x, Matrix jacobian) throws FunctionNotDefinedException;
	
	public Integer getDomainDimension();
	
	public Integer getCoDomainDimension();
	
}
