package de.varylab.discreteconformal.math.optimization;


import no.uib.cipr.matrix.Vector;


/**
 * Base interface for an Solver
 * <p>
 * Copyright 2005 <a href="http://www.sechel.de">Stefan Sechelmann</a>
 * <a href="http://www.math.tu-berlin.de/geometrie">TU-Berlin</a> 
 * @author Stefan Sechelmann
 * @see math.optimization.newton.NewtonOptimizer
 */
public interface Solver {

	/**
	 * Solves the problem fun(x)=b and stores the result in guess
	 * @param fun a linearizable function
	 * @param x a guess for x
	 * @param b
	 * @throws NotConvergentException 
	 */
	public void solve(Linearizable fun, Vector x, Vector b) throws FunctionNotDefinedException, NotConvergentException;
	
	public void setIterationMonitor(IterationMonitor monitor);

	public void setMaxIterations(Integer maxIterations);
	
	public Integer getMaxIterations();
	
	public void setError(Double error);
	
	public Double getError();
	
}
