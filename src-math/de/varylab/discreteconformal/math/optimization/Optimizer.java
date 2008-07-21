package de.varylab.discreteconformal.math.optimization;

import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Vector.Norm;
import de.varylab.discreteconformal.math.optimization.stepcontrol.StepController;


/**
 * Base interface for an Optimizer
 * <p>
 * Copyright 2005 <a href="http://www.sechel.de">Stefan Sechelmann</a>
 * <a href="http://www.math.tu-berlin.de/geometrie">TU-Berlin</a> 
 * @author Stefan Sechelmann
 * @see math.optimization.newton.NewtonOptimizer
 */
public interface Optimizer {
	
	/**
	 * Minimizes the given optimizable
	 * @param guess the first guess
	 * @param func the function to optimize
	 */
	public void minimize(Vector guess, Optimizable func) throws NotConvergentException;
	
	public void setIterationMonitor(IterationMonitor monitor);

	public void setMaxIterations(Integer maxIterations);
	
	public Integer getMaxIterations();
	
	public void setError(Double error);
	
	public Double getError();
	
	public Norm getNorm();

	public void setNorm(Norm norm);
	
	public void setStepController(StepController stepController);
	
	public StepController getStepController();
	
}
