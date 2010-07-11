package de.varylab.discreteconformal.util;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;

public class SurgeryUtility {

	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void cutAndGluePaths(List<E> p1, List<E> p2) {
		if (p1.size() != p2.size()) {
			throw new IllegalArgumentException("Paths of different lengths in cutAndGluePaths()");
		}
		V start = p1.get(0).getStartVertex();
		Set<V> allVertices = new HashSet<V>();
		allVertices.addAll(PathUtility.getVerticesOnPath(p1));
		allVertices.addAll(PathUtility.getVerticesOnPath(p2));
		
		// link sheets
		for (int i = 0; i < p1.size(); i++) {
			E e1 = p1.get(i);
			E e2 = p2.get(i);
			E e1Opp = e1.getOppositeEdge();
			E e2Opp = e2.getOppositeEdge();
			V v1 = e1Opp.getTargetVertex();
			V v2 = e2Opp.getTargetVertex();
			e1.linkOppositeEdge(e2Opp);
			e2.linkOppositeEdge(e1Opp);
			e2Opp.setTargetVertex(v1);
			e1Opp.setTargetVertex(v2);
		}
		
		// relocate vertices
		Set<V> isolatedVertices = new HashSet<V>();
		for (V v : allVertices) {
			E e0 = v.getIncomingEdge();
			if (e0 == null) {
				isolatedVertices.add(v);
				continue;
			}
			E e = e0;
			do {
				e = e.getNextEdge();
				e = e.getOppositeEdge();
				e.setTargetVertex(v);
			} while (e != e0);
		}
		
		// remove isolated vertices
		assert isolatedVertices.size() == 2;
		HalfEdgeDataStructure<V, E, F> hds = start.getHalfEdgeDataStructure();
		for (V v : isolatedVertices) {
			hds.removeVertex(v);
		}
	}

}
