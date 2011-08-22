package de.varylab.discreteconformal.functional.hds;

import de.jtem.halfedge.HalfEdgeDataStructure;

public class ConformalHDS extends HalfEdgeDataStructure<MyConformalVertex, MyConformalEdge, MyConformalFace> {

	public ConformalHDS() {
		super(MyConformalVertex.class, MyConformalEdge.class, MyConformalFace.class);
	}
	
}
