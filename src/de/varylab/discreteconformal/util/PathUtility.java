package de.varylab.discreteconformal.util;

import java.util.HashSet;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.varylab.discreteconformal.util.Search.WeightAdapter;

public class PathUtility {

	
	/**
	 * Returns a path which does not include opposite edges of each edge.
	 * All edges of the result will be positive
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
	> Set<E> getHalfPath(Iterable<E> path) {
		Set<E> result = new TreeSet<E>(new NodeIndexComparator<E>());
		Set<E> negSet = new HashSet<E>();
		for (E e : path) {
			if (negSet.contains(e) || result.contains(e)) {
				continue;
			}
			if (e.isPositive()) {
				result.add(e);
				negSet.add(e.getOppositeEdge());
			} else {
				result.add(e.getOppositeEdge());	
				negSet.add(e);
			}
		}
		return result;
	}
	
	
	
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
	> Set<V> getVerticesOnPath(Iterable<E> path) {
		Set<V> result = new TreeSet<V>(new NodeIndexComparator<V>());
		for (E e : path) {
			result.add(e.getStartVertex());
			result.add(e.getTargetVertex());
		}
		return result;
	}
	
	
	/**
	 * Returns a path which includes the opposite edges of each edge
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
	> Set<E> getFullPath(Iterable<E> path) {
		Set<E> result = new TreeSet<E>(new NodeIndexComparator<E>());
		for (E e : path) {
			result.add(e);
			result.add(e.getOppositeEdge());
		}
		return result;
	}
	
	
	
	/**
	 * The path must not contain opposite edges!
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param path
	 * @param wa
	 * @return
	 */
	public static  <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double getTotalPathWeight(Set<E> path, WeightAdapter<E> wa) {
		double result = 0.0;
		for (E e : path) {
			result += wa.getWeight(e);
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
	 * TODO This is not ready yet!
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
		leftPath.push(startLeft);
		rightPath.push(startRight);
		leftVisited.add(startLeft);
		rightVisited.add(startRight);
		
		HalfEdgeDataStructure<V, E, F> hds = startLeft.getHalfEdgeDataStructure();
		Set<V> validVertices = new HashSet<V>(hds.getVertices()); 
		validVertices.removeAll(cycle);
		
		
		
		return false;
	}
	

}
