package de.varylab.discreteconformal.util;

import static de.jtem.halfedge.util.HalfEdgeUtils.boundaryEdges;
import static de.jtem.halfedge.util.HalfEdgeUtils.boundaryVertices;
import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;

public class SurgeryUtility {

	/**
	 * Glue surface along boundary components starting with
	 * edges e1 and e2.  
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param e1
	 * @param e2
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Set<E> glueAlongBoundaries(E be1, E be2) {
		HalfEdgeDataStructure<V, E, F> hds = be1.getHalfEdgeDataStructure();
		List<E> b1 = boundaryEdges(be1);
		List<E> b2 = boundaryEdges(be2);
		Set<V> vb2 = PathUtility.getVerticesOnPath(b2);
		Set<E> edgeCycle = new HashSet<E>();
		if (b1.size() != b2.size()) {
			throw new RuntimeException("Boundary components have different lengths in glueAlongBoundaries()");
		}
		// collect data to link
		Map<E, E> oppMap = new HashMap<E, E>(); // new opposites
		Map<V, List<E>> targetMap = new HashMap<V, List<E>>(); // new targets
		E e1 = be1;
		E e2 = be2;
		do {
			oppMap.put(e1.getOppositeEdge(), e2.getOppositeEdge());
			List<E> v2Star = incomingEdges(e2.getStartVertex());
			targetMap.put(e1.getTargetVertex(), v2Star);
			e1 = e1.getNextEdge();
			e2 = e2.getPreviousEdge();
		} while (e1 != be1 && e2 != be2);
		assert e1 == be1 && e2 == be2 : "cycles are not completed";
		
		// re-link
		for (E e : oppMap.keySet()) {
			E opp = oppMap.get(e);
			e.linkOppositeEdge(opp);
			edgeCycle.add(e);
			edgeCycle.add(opp);
		}
		// set targets
		for (V v : targetMap.keySet()) {
			for (E e : targetMap.get(v)) {
				e.setTargetVertex(v);
			}
		}
			
		// remove waste
		for (E e : b1) hds.removeEdge(e);
		for (E e : b2) hds.removeEdge(e);
		for (V v : vb2) hds.removeVertex(v);
		return edgeCycle;
	} 

	
	
	/**
	 * Identifies the vertices on the given boundaries. The vertices
	 * of boundary one are retained
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param b1
	 * @param b2
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void glueAlongFaces(F f1, F f2, V v1, V v2) {
		HalfEdgeDataStructure<V, E, F> hds = f1.getHalfEdgeDataStructure();
		List<E> b1 = boundaryEdges(f1);
		List<E> b2 = boundaryEdges(f2);
		List<V> vb1 = boundaryVertices(f1);
		List<V> vb2 = boundaryVertices(f2);
		if (!vb1.contains(v1) || !vb2.contains(v2)) {
			throw new RuntimeException("Vertex not on boundary in glueAlongFaces()");
		}
		if (vb1.size() != vb2.size()) {
			throw new RuntimeException("Face boundaries have different lengths in glueAlongFaces()");
		}
		
		// find first edges to glue
		E be1 = null;
		E be2 = null;
		for (E e : incomingEdges(v1)) {
			if (e.getLeftFace() == f1) {
				be1 = e;
			}
		}
		for (E e : incomingEdges(v2)) {
			if (e.getLeftFace() == f2) {
				be2 = e.getNextEdge();
			}
		}
		assert be1 != null && be2 != null : "bundary corrupt";
		if (be1 == null || be2 == null) {
			throw new RuntimeException("corrupt boundary in glueAlongFaces().");
		}
		assert be1.getTargetVertex() == v1 && be2.getStartVertex() == v2 : "wrong vertices";
		if (be1.getTargetVertex() != v1 || be2.getStartVertex() != v2) {
			throw new RuntimeException("wrong vertices in glueAlongFaces().");
		}
		
		// collect data to link
		Map<E, E> oppMap = new HashMap<E, E>(); // new opposites
		Map<V, List<E>> targetMap = new HashMap<V, List<E>>(); // new targets
		E e1 = be1;
		E e2 = be2;
		do {
			oppMap.put(e1.getOppositeEdge(), e2.getOppositeEdge());
			List<E> v2Star = incomingEdges(e2.getStartVertex());
			targetMap.put(e1.getTargetVertex(), v2Star);
			e1 = e1.getNextEdge();
			e2 = e2.getPreviousEdge();
		} while (e1 != be1 && e2 != be2);
		assert e1 == be1 && e2 == be2 : "cycles are not completed";
		
		// re-link
		for (E e : oppMap.keySet()) {
			E opp = oppMap.get(e);
			e.linkOppositeEdge(opp);
		}
		// set targets
		for (V v : targetMap.keySet()) {
			for (E e : targetMap.get(v)) {
				e.setTargetVertex(v);
			}
		}
			
		// remove waste
		for (E e : b1) hds.removeEdge(e);
		for (E e : b2) hds.removeEdge(e);
		for (V v : vb2) hds.removeVertex(v);
		hds.removeFace(f1);
		hds.removeFace(f2);
	}
	
	
	
	
	
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
