package de.varylab.discreteconformal.math.optimization.stepcontrol;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import de.varylab.discreteconformal.math.optimization.Optimizable;


/**
 * A step width controller wich implements the armijo line search
 * <p>
 * Copyright 2005 <a href="http://www.sechel.de">Stefan Sechelmann</a>
 * <a href="http://www.math.tu-berlin.de/geometrie">TU-Berlin</a> 
 * @author Stefan Sechelmann
 * @see <a href="http://www.stanford.edu/~boyd/cvxbook/">Stephen Boyd 
 * and Lieven Vandenberghe - Convex Optimization </a>
 */
public class ArmijoStepController implements StepController {

	private Double
		alpha = 0.2,
		beta = 0.5;
	
	
	public Double step(Vector x, Double value, Vector dx, Optimizable func, Vector grad, Matrix hess) {
		Vector oldX = new DenseVector(x);
		Vector oldGrad = new DenseVector(grad);
		x.add(dx);
		Double result = func.evaluate(x, grad, hess);
//		---
//		if (result < value + alpha * oldGrad.dot(dx)){
//			oldGrad.scale(-1);
//			grad.scale(-1);
//		}
//		----
		while (result > value + alpha * oldGrad.dot(dx)){
			dx.scale(beta);
			x.set(oldX).add(dx);
			result = func.evaluate(x, grad, hess);
		}
		return result;
	}


	public Double getAlpha() {
		return alpha;
	}


	public void setAlpha(Double alpha) {
		this.alpha = alpha;
	}


	public Double getBeta() {
		return beta;
	}


	public void setBeta(Double beta) {
		this.beta = beta;
	}

}
