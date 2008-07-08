package de.varylab.discreteconformal.heds;

import de.jtem.halfedge.HalfEdgeDataStructure;

public class HDS extends HalfEdgeDataStructure<CVertex, CEdge, CFace> {

	public HDS() {
		super(CVertex.class, CEdge.class, CFace.class);
	}

}
