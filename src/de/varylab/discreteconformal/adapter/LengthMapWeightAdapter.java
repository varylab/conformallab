package de.varylab.discreteconformal.adapter;

import java.util.Map;

import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.util.Search.WeightAdapter;

public class LengthMapWeightAdapter implements WeightAdapter<CoEdge> {

	private Map<CoEdge, Double>
		lMap = null;
	
	public LengthMapWeightAdapter(Map<CoEdge, Double> lMap) {
		this.lMap = lMap;
	}
	
	@Override
	public double getWeight(CoEdge e) {
		if (!lMap.containsKey(e)) {
			return Double.POSITIVE_INFINITY;
		} else {
			return lMap.get(e);
		}
	}

}
