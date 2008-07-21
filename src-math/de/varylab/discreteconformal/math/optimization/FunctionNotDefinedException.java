package de.varylab.discreteconformal.math.optimization;


/**
 * Thrown if the domain value of a Linearizable is invalid.
 * <p>
 * Copyright 2005 <a href="http://www.sechel.de">Stefan Sechelmann</a>
 * <a href="http://www.math.tu-berlin.de/geometrie">TU-Berlin</a> 
 * @author Stefan Sechelmann
 * @see math.optimization.Linearizable
 */
@SuppressWarnings("serial")
public class FunctionNotDefinedException extends Exception {

	public FunctionNotDefinedException() {
		super();
	}

	public FunctionNotDefinedException(String message) {
		super(message);
	}

	public FunctionNotDefinedException(String message, Throwable cause) {
		super(message, cause);
	}

	public FunctionNotDefinedException(Throwable cause) {
		super(cause);
	}

}
