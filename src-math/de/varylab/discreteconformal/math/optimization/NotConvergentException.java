package de.varylab.discreteconformal.math.optimization;


/**
 * Thrown if a calculation is not converging 
 * <p>
 * Copyright 2005 <a href="http://www.sechel.de">Stefan Sechelmann</a>
 * <a href="http://www.math.tu-berlin.de/geometrie">TU-Berlin</a> 
 * @author Stefan Sechelmann
 */
@SuppressWarnings("serial")
public class NotConvergentException extends Exception {

	private Double
		error = 0.0;
	
	public NotConvergentException(String message, Double reachedError) {
		super(message);
		error = reachedError;
	}
	
	public String getMessage() {
		return super.getMessage() + "\nReached Error: " + error;
	}
	
}
