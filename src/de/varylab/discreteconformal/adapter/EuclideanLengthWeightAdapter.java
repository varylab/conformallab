package de.varylab.discreteconformal.adapter;

import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.util.Search.WeightAdapter;

public class EuclideanLengthWeightAdapter implements WeightAdapter<CoEdge> {

	@Override
	public double getWeight(CoEdge e) {
		return e.getLength();
	}
	
}
