package de.varylab.discreteconformal.plugin.adapter;

import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.util.Search.WeightAdapter;

public class EuclideanLengthWeightAdapter implements WeightAdapter<CoEdge> {

	
	public EuclideanLengthWeightAdapter() {
	}
	
	
	@Override
	public double getWeight(CoEdge e) {
		return e.getLength();
	}
	
}
