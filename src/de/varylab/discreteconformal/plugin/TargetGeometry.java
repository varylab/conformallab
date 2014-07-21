package de.varylab.discreteconformal.plugin;

import java.util.List;

import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;

public enum TargetGeometry {
	Automatic,
	Euclidean,
	Spherical,
	Hyperbolic;

	public static TargetGeometry calculateTargetGeometry(CoHDS hds) {
		int g = HalfEdgeUtils.getGenus(hds);
		List<List<CoEdge>> bc = HalfEdgeUtils.boundaryComponents(hds);
		return TargetGeometry.calculateTargetGeometry(g, bc.size());
	}

	public static TargetGeometry calculateTargetGeometry(int g, int numBoundaryComponents) {
		switch (g) {
		case 0:
			if (numBoundaryComponents > 0) {
				return Euclidean;
			} else {
				return Spherical;
			}
		case 1:
			return Euclidean;
		default:
			return Hyperbolic;
		}
	}
}