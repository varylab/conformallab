package de.varylab.discreteconformal.plugin.adapter;

import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.sqrt;

import java.util.HashSet;
import java.util.Set;

import no.uib.cipr.matrix.Vector;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.Search.WeightAdapter;

public class HyperbolicLengthWeightAdapter implements WeightAdapter<CoEdge> {

	private Vector
		u = null;
	private Set<CoEdge> 
		infiniteWeightSet = new HashSet<CoEdge>();
	
	public HyperbolicLengthWeightAdapter(Vector u) {
		this.u = u;
	}
	
	
	@Override
	public double getWeight(CoEdge e) {
		if (infiniteWeightSet.contains(e)) {
			return Double.POSITIVE_INFINITY;
		} else {
			return getNewLength(e);
		}
	}

	/**
	 * Calculate the edge length for the flat metric
	 * @param e
	 * @param u
	 * @return the new edge length
	 */
	public Double getNewLength(CoEdge e) {
		CoVertex v1 = e.getStartVertex();
		CoVertex v2 = e.getTargetVertex();
		Double u1 = v1.getSolverIndex() >= 0 ? u.get(v1.getSolverIndex()) : 0.0; 
		Double u2 = v2.getSolverIndex() >= 0 ? u.get(v2.getSolverIndex()) : 0.0;
		Double lambda = e.getLambda();
		Double lambdaNew = lambda + u1 + u2;
		return 2 * arsinh( exp(lambdaNew / 2) );
	}
	
	
	
	private double arsinh(double x) {
		double r = x + sqrt(x*x + 1);
		return log(r);
	}
	
	
	@Override
	public void setInfiniteWeightPaths(Set<CoEdge> paths) {
		infiniteWeightSet = paths;
	}
	
}
