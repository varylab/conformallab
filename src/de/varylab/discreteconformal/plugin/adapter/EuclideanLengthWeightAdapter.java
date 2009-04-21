package de.varylab.discreteconformal.plugin.adapter;

import java.util.HashSet;
import java.util.Set;

import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.util.Search.WeightAdapter;

public class EuclideanLengthWeightAdapter implements WeightAdapter<CoEdge> {

	private Set<CoEdge> 
		infiniteWeightSet = new HashSet<CoEdge>();
	
	@Override
	public double getWeight(CoEdge e) {
		if (infiniteWeightSet.contains(e)) {
			return Double.POSITIVE_INFINITY;
		} else {
			return e.getLength();
		}
	}
	
	@Override
	public void setInfiniteWeightPaths(Set<CoEdge> paths) {
		infiniteWeightSet = paths;
	}
	
}
