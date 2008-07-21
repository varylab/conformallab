package de.varylab.discreteconformal.math.optimization.stepcontrol;

import static no.uib.cipr.matrix.Vector.Norm.Two;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import de.varylab.discreteconformal.math.optimization.Optimizable;

public class ShortGradientStepController implements StepController {

	private Double 
		alpha = 0.5;
	
	
	public Double step(Vector x, Double value, Vector dx, Optimizable func, Vector grad, Matrix hess) {
		Vector oldX = new DenseVector(x);
		Vector oldGrad = new DenseVector(grad);
		Double oldGradLength = oldGrad.norm(Two);
		
		x.add(dx);
		Double result = func.evaluate(x, grad, hess);
		Double gradLength = grad.norm(Two);
		int counter = 0; 
		while (oldGradLength <= gradLength || gradLength.equals(Double.NaN)) {
			dx.scale(alpha);
			x.set(oldX).add(dx);
			result = func.evaluate(x, grad, hess);
			gradLength = grad.norm(Two);
			counter++;
			if (counter == 100)
				throw new RuntimeException("No valid step in step controller!");
		}
		return result;
	}


	public Double getAlpha() {
		return alpha;
	}


	public void setAlpha(Double alpha) {
		this.alpha = alpha;
	}

}
