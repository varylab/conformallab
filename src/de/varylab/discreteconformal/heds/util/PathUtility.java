package de.varylab.discreteconformal.heds.util;

import java.util.HashSet;
import java.util.Set;
import java.util.Stack;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;

public class PathUtility {

	
	/**
	 * Returns the vertices on a given path
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param path
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Set<V> getVerticesOnPath(Set<E> path) {
		HashSet<V> result = new HashSet<V>();
		for (E e : path) {
			result.add(e.getStartVertex());
			result.add(e.getTargetVertex());
		}
		return result;
	}
	
	
	
	/**
	 * Checks whether a given cycle is a simple path 
	 * i.e. it has no self-intersections
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param cycle a cycle of consistently oriented half-edges
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> boolean isCycleSimple(Set<E> cycle) {
		Set<V> vSet = new HashSet<V>();
		for (E e : cycle) {
			V v = e.getTargetVertex();
			if (vSet.contains(v)) {
				return false;
			}
			vSet.add(v);
		}
		return true;
	}
	
	
	
	/**
	 * Checks whether a cycle is essential
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param cycle must be a cycle of consistently oriented half-edges
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> boolean isCycleEssential(Set<E> cycle) {
		if (cycle.size() == 0) {
			return false;
		}
		if (!isCycleSimple(cycle)) {
			return false;
		}
		E seed = cycle.iterator().next();
		V startLeft = seed.getNextEdge().getTargetVertex();
		V startRight = seed.getOppositeEdge().getNextEdge().getTargetVertex();
		Set<V> leftVisited = new HashSet<V>();
		Set<V> rightVisited = new HashSet<V>();
		Stack<V> leftPath = new Stack<V>();
		Stack<V> rightPath = new Stack<V>();
		leftPath.add(startLeft);
		rightPath.add(startRight);
		leftVisited.add(startLeft);
		rightVisited.add(startRight);
		
		HalfEdgeDataStructure<V, E, F> hds = startLeft.getHalfEdgeDataStructure();
		Set<V> validVertices = new HashSet<V>(hds.getVertices()); 
		validVertices.removeAll(cycle);
		
		
		
		return false;
	}
	


}