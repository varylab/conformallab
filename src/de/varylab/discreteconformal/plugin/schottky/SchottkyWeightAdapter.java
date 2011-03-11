package de.varylab.discreteconformal.plugin.schottky;

import java.util.Set;

import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.util.Search.WeightAdapter;

public class SchottkyWeightAdapter implements WeightAdapter<CoEdge> {

	private Set<Set<CoEdge>>
		cycles = null;
	
	public SchottkyWeightAdapter(Set<Set<CoEdge>> cycles) {
		this.cycles = cycles;
	}
	
	@Override
	public double getWeight(CoEdge e) {
		for (Set<CoEdge> cycle : cycles) {
			if (cycle.contains(e)) {
				return 0;
			}
		}
		return 1;
	}

}
