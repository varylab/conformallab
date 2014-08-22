package de.varylab.discreteconformal.functional.node;

import de.jtem.halfedge.Edge;

public class ConformalEdge <
	V extends ConformalVertex<V, E, F>,
	E extends ConformalEdge<V, E, F>,
	F extends ConformalFace<V, E, F>
> extends Edge<V, E, F> {

    protected double
    	lambda = 1.0,
    	alpha = 0.0,
    	theta = 0.0,
    	phi = Math.PI;
	protected Integer
		solverIndex = -1;
	
	public Integer getSolverIndex() {
		return solverIndex;
	}
	public void setSolverIndex(Integer solverIndex) {
		this.solverIndex = solverIndex;
	}

	public double getLambda() {
		return lambda;
	}
	public void setLambda(double lambda) {
		this.lambda = lambda;
	}

	public double getAlpha() {
		return alpha;
	}
	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}
	
	public double getPhi() {
		return phi;
	}
	public void setPhi(double phi) {
		this.phi = phi;
	}
	
	public double getTheta() {
		return theta;
	}
	public void setTheta(double theta) {
		this.theta = theta;
	}
	
	@Override
	public void copyData(E e) {
		lambda = e.lambda;
		alpha = e.alpha;
		theta = e.theta;
		solverIndex = e.solverIndex;
	};
	
}
