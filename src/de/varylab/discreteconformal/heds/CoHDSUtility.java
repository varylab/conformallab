package de.varylab.discreteconformal.heds;

import java.util.Map;

import de.jtem.halfedge.util.HalfEdgeUtils;

public final class CoHDSUtility {

	private CoHDSUtility() {
	}

	public static void createSurfaceCopy(CoHDS source, CoHDS target, Map<CoVertex, CoVertex> vertexMap) {
		HalfEdgeUtils.copy(source, target);
		for (int i = 0; i < source.numVertices(); i++) {
			CoVertex v = source.getVertex(i);
			CoVertex vv = target.getVertex(i);
			vv.P = v.P;
			vv.info = v.info;
			vertexMap.put(v, vv);
		}
		for (int i = 0; i < source.numEdges(); i++) {
			CoEdge e = source.getEdge(i);
			CoEdge ee = target.getEdge(i);
			ee.info = e.info; 
		}
	}
	
}
