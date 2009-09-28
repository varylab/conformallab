package de.jtem.halfedgetools.functional.conformal.node;

import static java.lang.Math.PI;
import de.jtem.halfedge.Vertex;

public class ConformalVertex <
	V extends ConformalVertex<V, E, F>,
	E extends ConformalEdge<V, E, F>,
	F extends ConformalFace<V, E, F>
> extends Vertex<V, E, F> {

	private double 
		theta = 2 * PI;
	private Integer
		solverIndex = -1;

	public double getTheta() {
		return theta;
	}

	public void setTheta(double theta) {
		this.theta = theta;
	}

	public Integer getSolverIndex() {
		return solverIndex;
	}

	public void setSolverIndex(Integer solverIndex) {
		this.solverIndex = solverIndex;
	}
	
	@Override
	public void copyData(V v) {
		setSolverIndex(v.getSolverIndex());
		setTheta(v.getTheta());
	};
	
}
