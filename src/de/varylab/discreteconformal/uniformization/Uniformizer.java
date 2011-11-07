package de.varylab.discreteconformal.uniformization;

import java.util.Map;

import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.EuclideanUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.HyperbolicUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.UnwrapException;

public class Uniformizer {

	public static Map<CoVertex, Double> uniformizationFactors(CoHDS hds, AdapterSet a) throws UnwrapException {
		int g = HalfEdgeUtils.getGenus(hds);
		switch (g) {
		case 0:
			throw new RuntimeException("genus 0 not supported in Uniformizer");
		case 1:
			EuclideanUnwrapperPETSc euclUnwrapper = new EuclideanUnwrapperPETSc();
			return euclUnwrapper.calculateConformalFactors(hds, a);
		default:
			HyperbolicUnwrapperPETSc hypUnwrapper = new HyperbolicUnwrapperPETSc();
			return hypUnwrapper.calculateConformalFactors(hds, a);
		}
	}
	
}
