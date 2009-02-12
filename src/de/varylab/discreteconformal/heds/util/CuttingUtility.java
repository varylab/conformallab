package de.varylab.discreteconformal.heds.util;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static de.jtem.halfedge.util.HalfEdgeUtils.isBoundaryVertex;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;




public class CuttingUtility {

	/**
	 * Cuts an edge and returns the mapping oldVertex <-> newVertex for split vertices
	 * @param <V>
	 * @param <E>
	 * @param <F>
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
		HalfEdgeDataStructure<V, E, F> graph = edge.getHalfEdgeDataStructure();
		V v1 = edge.getStartVertex();
		V v2 = edge.getTargetVertex();
		
		E opp = edge.getOppositeEdge();
	
		boolean splitV1 = isBoundaryVertex(v1);
		boolean splitV2 = isBoundaryVertex(v2);
		
		List<E> v1Star = incomingEdges(v1);
		List<E> v2Star = incomingEdges(v2);
		
		E new1 = graph.addNewEdge();
		E new2 = graph.addNewEdge();
		
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
				System.err.println("do1");
				newTargetEdges.add(actEdge);
				actEdge = actEdge.getNextEdge().getOppositeEdge();
			} while (actEdge != b);
			newTargetEdges.add(b);
			
			V newV = graph.addNewVertex();
			result.put(v1, newV);
			
			b.linkNextEdge(new1);
			new2.linkNextEdge(b2);
			
			for (E e : newTargetEdges)
				e.setTargetVertex(newV);
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
				System.err.println("do2");
				newTargetEdges.add(actEdge);
				actEdge = actEdge.getNextEdge().getOppositeEdge();
			} while (actEdge != b);
			newTargetEdges.add(b);
			
			V newV = graph.addNewVertex();
			result.put(v2, newV);
			
			b.linkNextEdge(new2);
			new1.linkNextEdge(b2);
			
			for (E e : newTargetEdges)
				e.setTargetVertex(newV);
		}
		return result;
		
	}
	
	
}
