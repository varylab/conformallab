package de.varylab.discreteconformal.heds.util;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static de.jtem.halfedge.util.HalfEdgeUtils.isBoundaryVertex;
import static de.jtem.halfedge.util.HalfEdgeUtils.isManifoldVertex;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;




public class CuttingUtility {

	/**
	 * Cuts an edge and returns the mapping 
	 * oldVertex -> newVertex for split vertices
	 * @param edge
	 * @return
	 * @throws IllegalArgumentException
	 */
	public static 
	<
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Map<V, V> cutAtEdge(E edge) throws IllegalArgumentException{
		Map<V, V> result = new HashMap<V, V>();
		HalfEdgeDataStructure<V, E, F> hds = edge.getHalfEdgeDataStructure();
		V v1 = edge.getStartVertex();
		V v2 = edge.getTargetVertex();
		E opp = edge.getOppositeEdge();
		if (v1 == v2) {
			return cutLoopEdge(edge);
		}
		
		boolean splitV1 = isBoundaryVertex(v1);
		boolean splitV2 = isBoundaryVertex(v2);
		
		List<E> v1Star = incomingEdges(v1);
		List<E> v2Star = incomingEdges(v2);
		
		E new1 = hds.addNewEdge();
		E new2 = hds.addNewEdge();
		
		new1.linkOppositeEdge(opp);
		new2.linkOppositeEdge(edge);
		new1.setTargetVertex(v2);
		new2.setTargetVertex(v1);
		
		new1.linkNextEdge(new2);
		new2.linkNextEdge(new1);
		
		if (splitV1){
			E b = null;
			for (E e : v1Star) {
				if (e.getLeftFace() == null){
					b = e;
					break;
				}
			}
			assert (b != null);
			E b2 = b.getNextEdge();
			List<E> newTargetEdges = new LinkedList<E>();
			E actEdge = opp;
			do {
				newTargetEdges.add(actEdge);
				actEdge = actEdge.getNextEdge().getOppositeEdge();
			} while (actEdge != b);
			newTargetEdges.add(b);
			
			V newV = hds.addNewVertex();
			result.put(v1, newV);
			
			b.linkNextEdge(new1);
			new2.linkNextEdge(b2);
			
			for (E e : newTargetEdges) {
				e.setTargetVertex(newV);
			}
		}
		if (splitV2){
			E b = null;
			for (E e : v2Star) {
				if (e.getLeftFace() == null){
					b = e;
					break;
				}
			}
			assert (b != null);
			E b2 = b.getNextEdge();
			List<E> newTargetEdges = new LinkedList<E>();
			E actEdge = edge;
			do {
				newTargetEdges.add(actEdge);
				actEdge = actEdge.getNextEdge().getOppositeEdge();
			} while (actEdge != b);
			newTargetEdges.add(b);
			
			V newV = hds.addNewVertex();
			result.put(v2, newV);
			
			b.linkNextEdge(new2);
			new1.linkNextEdge(b2);
			
			for (E e : newTargetEdges) {
				e.setTargetVertex(newV);
			}
		}
		return result;
		
	}
	
	
	
	/**
	 * Cuts an edge which is a loop and returns the mapping 
	 * oldVertex -> newVertex for split vertices
	 * @param edge
	 * @return
	 * @throws IllegalArgumentException
	 */
	public static 
	<
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Map<V, V> cutLoopEdge(E edge) throws IllegalArgumentException{
		Map<V, V> result = new HashMap<V, V>();
		HalfEdgeDataStructure<V, E, F> hds = edge.getHalfEdgeDataStructure();
		System.out.println("CutLoop 1: " + HalfEdgeUtils.isValidSurface(hds, true));
		V oldVertex = edge.getStartVertex();
		V v2 = edge.getTargetVertex();
		if (oldVertex != v2) {
			throw new IllegalArgumentException("Start vertex != target vertex in cutLootEdge()");
		}
		
		boolean repairNonManifold = isBoundaryVertex(oldVertex);
		
		E opp = edge.getOppositeEdge();
		V newVertex = hds.addNewVertex();
		result.put(oldVertex, newVertex);
		
		// get edges of the new vertex 
		List<E> newTargetEdges = new LinkedList<E>();
		E actEdge = edge;
		do {
			newTargetEdges.add(actEdge);
			actEdge = actEdge.getNextEdge().getOppositeEdge();
			// if this is no loop we have a dead lock here
		} while (actEdge != edge.getOppositeEdge()); 
		
		// create new edges
		E new1 = hds.addNewEdge();
		E new2 = hds.addNewEdge();
		
		new1.linkOppositeEdge(opp);
		new2.linkOppositeEdge(edge);
		new1.setTargetVertex(oldVertex);
		new2.setTargetVertex(newVertex);
		
		for (E e : newTargetEdges) {
			e.setTargetVertex(newVertex);
		}
		
		// link boundary loops
		new1.linkNextEdge(new1);
		new2.linkNextEdge(new2);

		System.out.println("CutLoop 2: " + HalfEdgeUtils.isValidSurface(hds, true));
		
		// repair the created non manifold vertex
		if (repairNonManifold) {
			System.out.println("repairing non manifold vertex...");
			V nmv = null;
			if (isManifoldVertex(oldVertex)) {
				nmv = newVertex;
			} else {
				nmv = oldVertex;
			}
			assert isManifoldVertex(nmv == oldVertex ? newVertex : oldVertex);
			
			V v = hds.addNewVertex();
			
			E be1 = null;
			List<E> star = incomingEdges(nmv);
			for (E e : star) {
				if (e.getLeftFace() == null) {
					be1 = e;
					break;
				}
			}
			assert be1 != null;
			
			E be2 = be1.getNextEdge();
			actEdge = be1;
			while (actEdge.getLeftFace() != null) {
				actEdge.setTargetVertex(v);
				actEdge = actEdge.getOppositeEdge().getPreviousEdge();
			}
			E be4 = actEdge;
			E be3 = be4.getNextEdge();
			
			be1.linkNextEdge(be3);
			be4.linkNextEdge(be2);
		}
		
		System.out.println("CutLoop 3: " + HalfEdgeUtils.isValidSurface(hds, true));
		
		return result;
	}

	
}
