package de.varylab.discreteconformal.adapter;

import static java.lang.Math.exp;
import no.uib.cipr.matrix.Vector;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.Search.WeightAdapter;

public class EuclideanLengthWeightAdapter implements WeightAdapter<CoEdge> {

	private Vector
		u = null;
	
	public EuclideanLengthWeightAdapter(Vector u) {
		this.u = u;
	}
	
	
	@Override
	public double getWeight(CoEdge e) {
		return getNewLength(e);
	}
	
	/**
	 * Calculate the edge length for the flat metric
	 * @param e
	 * @param u
	 * @return the new edge length
	 */
	public Double getNewLength(CoEdge e) {
		if (u == null) {
			return e.getLength();
		}
		CoVertex v1 = e.getStartVertex();
		CoVertex v2 = e.getTargetVertex();
		Double u1 = v1.getSolverIndex() >= 0 ? u.get(v1.getSolverIndex()) : 0.0; 
		Double u2 = v2.getSolverIndex() >= 0 ? u.get(v2.getSolverIndex()) : 0.0;
		Double lambda = e.getLambda();
		return exp(lambda + u1 + u2);
	}
	
}
