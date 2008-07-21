package de.varylab.discreteconformal.math.optimization;


/**
 * IterationMonitor interface for the NewtonOptimizer and NewtonSolver
 * <p>
 * Copyright 2005 <a href="http://www.sechel.de">Stefan Sechelmann</a>
 * <a href="http://www.math.tu-berlin.de/geometrie">TU-Berlin</a> 
 * @author Stefan Sechelmann
 * @see math.optimization.newton.NewtonOptimizer
 * @see math.optimization.newton.NewtonSolver
 */
public interface IterationMonitor {

	public void setIteration(Integer iterations, Double error);
	
	public void done(Double error);
	
	public void start(Double error);
	
}
