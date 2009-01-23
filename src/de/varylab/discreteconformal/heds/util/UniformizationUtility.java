package de.varylab.discreteconformal.heds.util;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static de.jtem.halfedge.util.HalfEdgeUtils.neighboringVertices;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;

public class UniformizationUtility {

	
	public static interface EdgeFlipAdapter<E extends Edge<?, E, ?>> {
		
		public void flip(E edge);
		
	}
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void reduceToFundamentalPolygon(
		HDS hds, 
		V root,
		EdgeFlipAdapter<E> flip
	) {
		while (hds.numVertices() > 1) {
			boolean removed = false;
			for (V v : hds.getVertices()) {
				if (v == root) {
					continue;
				}
				int starSize = neighboringVertices(v).size();
				if (starSize == 3) {
					removeTrivalentVertex(v);
					removed = true;
				}
			}
			V v = root.getIncomingEdge().getStartVertex();
			if (!removed && hds.numVertices() != 1) {
				while (incomingEdges(v).size() != 3) {
					flip.flip(v.getIncomingEdge());
				}
			}
			removeTrivalentVertex(v);
		}
	}
	
	
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void removeTrivalentVertex(V v) {
		HalfEdgeDataStructure<V, E, F> hds = v.getHalfEdgeDataStructure();
		E in1 = v.getIncomingEdge();
		E out1 = in1.getOppositeEdge();
		E outer1 = in1.getPreviousEdge();
		E in2 = out1.getPreviousEdge();
		E out2 = in2.getOppositeEdge();
		E outer2 = in2.getPreviousEdge();
		E in3 = out2.getPreviousEdge();
		E out3 = in3.getOppositeEdge();
		E outer3 = in3.getPreviousEdge();
		
		F f1 = in1.getLeftFace();
		F f2 = in2.getLeftFace();
		F f3 = in3.getLeftFace();

		hds.removeVertex(v);
		hds.removeEdge(in1);
		hds.removeEdge(in2);
		hds.removeEdge(in3);
		hds.removeEdge(out1);
		hds.removeEdge(out2);
		hds.removeEdge(out3);
		hds.removeFace(f2);
		hds.removeFace(f3);
		
		outer1.setLeftFace(f1);
		outer2.setLeftFace(f1);
		outer3.setLeftFace(f1);
		
		outer1.linkNextEdge(outer2);
		outer2.linkNextEdge(outer3);
		outer3.linkNextEdge(outer1);
	}
	
	
	
}
