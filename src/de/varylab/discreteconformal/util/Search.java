package de.varylab.discreteconformal.util;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static de.jtem.halfedge.util.HalfEdgeUtils.isInteriorEdge;
import static java.lang.Double.POSITIVE_INFINITY;
import static java.util.Collections.singleton;
import static java.util.Collections.sort;

import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;


public class Search {

	/**
	 * Breadth-first-search from start to end
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param start
	 * @param end
	 * @param avoidBorderEdges
	 * @return
	 * @throws NoSuchElementException
	 */
	public static 
	<
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<E> bFS(V start, V end, boolean avoidBorderEdges) throws NoSuchElementException{
		Set<E> allEdges = new HashSet<E>(start.getHalfEdgeDataStructure().getEdges());
		return bFS(allEdges, start, singleton(end), avoidBorderEdges, null);
	}

	public static 
	<
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<E> bFS(Collection<E> valid, V start, V end) throws NoSuchElementException{
		return bFS(valid, start, singleton(end), false, null);
	}
	
	public static 
	<
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<E> dualBFS(Set<E> valid, F start, F end) throws NoSuchElementException{
		return dualBFS(valid, start, singleton(end), false, null);
	}
	
	
	/**
	 * Breadth-first-search to the first hit in endPoints
	 * @param start
	 * @param end
	 * @param graph
	 * @return
	 * @throws NoSuchElementException
	 */
	public static 
	<
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<E> bFS(Collection<E> valid, V start, Collection<V> endPoints, boolean avoidBorderEdges, Comparator<E> comp) throws NoSuchElementException{
		if (endPoints.contains(start)) {
			return Collections.emptyList();
		}
		HashMap<V, Stack<E>> pathMap = new HashMap<V, Stack<E>>();
		LinkedList<V> queue = new LinkedList<V>();
		HashSet<V> visited = new HashSet<V>();
		V actVertex = start;
		queue.add(start);
		pathMap.put(start, new Stack<E>());
		while (!queue.isEmpty()){
			actVertex = queue.poll();
			Stack<E> path = pathMap.get(actVertex);
			List<E> star = incomingEdges(actVertex);
			if (comp != null) { 
				sort(star, comp);
			}
			for (E e : star){
				E pathEdge = e.getOppositeEdge();
				V v = pathEdge.getTargetVertex();
				if (visited.contains(v)) {
					continue;
				}
				if (!isInteriorEdge(e) && avoidBorderEdges) {
					continue;
				}
				if (!valid.contains(e)) {
					continue;
				}
				Stack<E> newPath = new Stack<E>();
				newPath.addAll(path);
				newPath.push(pathEdge);
				pathMap.put(v, newPath);
				if (!visited.contains(v)){
					visited.add(v);
					if (endPoints.contains(v))
						return newPath;
					else
						queue.offer(v);
				}
			}
		}
		throw new NoSuchElementException();
	}

	/**
	 * Breadth-first-search to the first hit in endPoints
	 * @param start
	 * @param end
	 * @param graph
	 * @return
	 * @throws NoSuchElementException
	 */
	public static 
	<
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<E> dualBFS(Set<E> valid, F start, Set<F> endPoints, boolean avoidBorder, Comparator<E> comp) throws NoSuchElementException{
		HashMap<F, Stack<E>> pathMap = new HashMap<F, Stack<E>>();
		LinkedList<F> queue = new LinkedList<F>();
		HashSet<F> visited = new HashSet<F>();
		F actFace = start;
		queue.add(start);
		pathMap.put(start, new Stack<E>());
		while (!queue.isEmpty()){
			actFace = queue.poll();
			Stack<E> path = pathMap.get(actFace);
			List<E> star = HalfEdgeUtils.boundaryEdges(actFace);
			if (comp != null) { 
				sort(star, comp);
			}
			for (E e : star){
				E pathEdge = e;//.getOppositeEdge();
				F f = pathEdge.getRightFace();
				if (visited.contains(f)) {
					continue;
				}
				if (!isInteriorEdge(e) && avoidBorder) {
					continue;
				}
				if (!valid.contains(e)) {
					continue;
				}
				Stack<E> newPath = new Stack<E>();
				newPath.addAll(path);
				newPath.push(pathEdge);
				pathMap.put(f, newPath);
				if (!visited.contains(f)){
					visited.add(f);
					if (endPoints.contains(f))
						return newPath;
					else
						queue.offer(f);
				}
			}
		}
		throw new NoSuchElementException();
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<E> bFS(V start, V end, Collection<V> avoid) throws NoSuchElementException {
		HashMap<V, Stack<E>> pathMap = new HashMap<V, Stack<E>>();
		LinkedList<V> queue = new LinkedList<V>();
		HashSet<V> visited = new HashSet<V>(avoid);
		V actVertex = start;
		queue.add(start);
		pathMap.put(start, new Stack<E>());
		while (!queue.isEmpty()){
			actVertex = queue.poll();
			Stack<E> path = pathMap.get(actVertex);
			List<E> star = incomingEdges(actVertex);
			for (E e : star){
				E pathEdge = e.getOppositeEdge();
				V v = pathEdge.getTargetVertex();
				if (visited.contains(v)) {
					continue;
				}
				Stack<E> newPath = new Stack<E>();
				newPath.addAll(path);
				newPath.push(pathEdge);
				pathMap.put(v, newPath);
				visited.add(v);
				if (end == v) {
					return newPath;
				} else {
					queue.offer(v);
				}
			}
		}
		throw new NoSuchElementException();
	}
	
	
	
	
	/**
	 * Depth first search - Untested!!
	 * @param start
	 * @param end
	 * @param graph
	 * @return
	 * @throws NoSuchElementException
	 */
	public static 
	<
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<E> dFS(V start, V end) throws NoSuchElementException{
		Stack<E> path = new Stack<E>();
		HashSet<V> visited = new HashSet<V>();
		visited.add(start);
		if (!dFS_R(start, end, path, visited)) {
			throw new NoSuchElementException();
		} else {
			return path;
		}
	} 


	protected static 
	<
		V extends Vertex<V, E, ?>,
		E extends Edge<V, E, ?>
	> boolean dFS_R(V start, V end, Stack<E> path, HashSet<V> visited) {
		if (start == end)
			return true;
		List<E> star = incomingEdges(start);
		for (E e : star){
			if (!visited.contains(e.getTargetVertex())){
				visited.add(e.getTargetVertex());
				path.push(e);
				if (dFS_R(e.getTargetVertex(), end, path, visited)) {
					return true;
				}
				path.pop();
			}
		}
		return false;
	}

	
	
	/**
	 * The weights for Dijkstra's algorithm 
	 * @param <E>
	 */
	public static interface WeightAdapter<E extends Edge<?, ?, ?>> {
		
		public double getWeight(E e);
		
	}
	
	
	public static class DefaultWeightAdapter <E extends Edge<?, ?, ?>>
		implements WeightAdapter<E> {
		
		@Override
		public double getWeight(E e) {
			return 1;
		};
	}
	
	
	
	/**
	 * Dijkstra's Algorithm
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param start
	 * @param ends
	 * @param w
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	>  TreeSet<E> getAllShortestPathsTree(V start, Set<V> ends, WeightAdapter<E> w, Set<V> avoid) {
		Map<V, List<E>> paths = getAllShortestPaths(start, ends, w, avoid);
		TreeSet<E> tree = new TreeSet<E>(new NodeIndexComparator<E>());
		for (V v : ends) {
			List<E> path = paths.get(v);
			for (E e : path) {
				tree.add(e);
				tree.add(e.getOppositeEdge());
			}
		}
		return tree;
	}
	
	
	
	
	/**
	 * Dijkstra's Algorithm
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param start
	 * @param ends
	 * @param w
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	>  Map<V, List<E>> getAllShortestPaths(V start, Set<V> ends, WeightAdapter<E> w, Set<V> avoid) {
		HashMap<V, List<E>> result = new HashMap<V, List<E>>();
		// init
		HalfEdgeDataStructure<V, E, F> hds = start.getHalfEdgeDataStructure();
		double[] d = new double[hds.numVertices()];
		int[] p = new int[hds.numVertices()];
		
		// run Dijkstra's Algorithm
		initializeSingleSource(hds, start, d, p);
		Set<V> S = new TreeSet<V>(new Comparator<V>() {
			@Override
			public int compare(V o1, V o2) {
				return o1.getIndex() - o2.getIndex();
			}
		});
		PriorityQueue<V> Q = new PriorityQueue<V>(hds.numVertices(), new DComparator<V>(d));
		Q.addAll(hds.getVertices());
		while (!Q.isEmpty()) {
			V u = Q.poll();
			S.add(u);
			for (E e : incomingEdges(u)) {
				if (avoid.contains(e.getStartVertex())) {
					continue;
				}
				relax(e.getOppositeEdge(), d, p, w, Q);
			}
		}
		
		// assemble the vertex path
		for (V v : ends) {
			V end = v;
			Stack<V> vPath = new Stack<V>();
			while (v != start) {
				assert !vPath.contains(v);
				vPath.push(v);
				if (p[v.getIndex()] < 0)
					break;
				v = hds.getVertex(p[v.getIndex()]);
			}
			vPath.push(start);

			Stack<E> ePath = new Stack<E>();
			V lastVertex = vPath.pop();
			V pathVertex = lastVertex;
			while (!vPath.isEmpty()) {
				lastVertex = pathVertex;
				pathVertex = vPath.pop();
				for (E e : HalfEdgeUtils.incomingEdges(lastVertex)) {
					if (e.getStartVertex() == pathVertex) {
						ePath.push(e.getOppositeEdge());
						break;
					}
				}
			}
			result.put(end, ePath);
		}
		
		return result;
	}
	
	
	
	
	/**
	 * Dijkstra's Algorithm
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param start
	 * @param end
	 * @param w
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	>  List<E> getShortestPath(V start, V end, WeightAdapter<E> w) {
		return getShortestPath(start, singleton(end), w);
	}
	
	
	
	/**
	 * Dijkstra's Algorithm
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param start
	 * @param ends
	 * @param w
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	>  List<E> getShortestPath(V start, Set<V> ends, WeightAdapter<E> w) {
		// init
		HalfEdgeDataStructure<V, E, F> hds = start.getHalfEdgeDataStructure();
		double[] d = new double[hds.numVertices()];
		int[] p = new int[hds.numVertices()];
		
		// run Dijkstra's Algorithm
		initializeSingleSource(hds, start, d, p);
		Set<V> S = new TreeSet<V>(new Comparator<V>() {
			@Override
			public int compare(V o1, V o2) {
				return o1.getIndex() - o2.getIndex();
			}
		});
		PriorityQueue<V> Q = new PriorityQueue<V>(hds.numVertices(), new DComparator<V>(d));
		Q.addAll(hds.getVertices());
		while (!Q.isEmpty()) {
			V u = Q.poll();
			S.add(u);
			for (E e : incomingEdges(u)) {
				relax(e.getOppositeEdge(), d, p, w, Q);
			}
		}
		
		// assemble the vertex path
		V v = Collections.min(ends, new DComparator<V>(d));
		Stack<V> vPath = new Stack<V>();
		while (v != start) {
			assert !vPath.contains(v);
			vPath.push(v);
			if (p[v.getIndex()] < 0)
				break;
			v = hds.getVertex(p[v.getIndex()]);
		}
		vPath.push(start);
		
		// get the edge path
		Stack<E> ePath = new Stack<E>();
		V lastVertex = vPath.pop();
		V pathVertex = lastVertex;
		while (!vPath.isEmpty()) {
			lastVertex = pathVertex;
			pathVertex = vPath.pop();
			for (E e : HalfEdgeUtils.incomingEdges(lastVertex)) {
				if (e.getStartVertex() == pathVertex) {
					ePath.push(e.getOppositeEdge());
					break;
				}
			}
		}
		return ePath;
	}
	
	
	protected static class DComparator<T extends Vertex<?,?,?>> implements Comparator<T> {

		private double[]
		    d = null;
		
		public DComparator(double[] d) {
			this.d = d;
		}
		
		@Override
		public int compare(T o1, T o2) {
			double val = d[o1.getIndex()] - d[o2.getIndex()];
			if (val == 0.0)
				return 0;
			else
				return val < 0.0 ? -1 : 1;
		}
		
	}
	
	
	
	/**
	 * Algorithm from the Cormen, Leiserson, Rivest, Stein
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param e
	 * @param d
	 * @param p
	 * @param w
	 * @param Q
	 */
	protected static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void relax(E e, double[] d, int[] p, WeightAdapter<E> w, PriorityQueue<V> Q) {
		V u = e.getStartVertex();
		V v = e.getTargetVertex();
		double weight = 1.0;
		if (w != null) weight = w.getWeight(e);
		double newWeight = d[u.getIndex()] + weight;
		if (d[v.getIndex()] > newWeight) {
			d[v.getIndex()] = newWeight;
			p[v.getIndex()] = u.getIndex();
			if (Q.remove(v)) { // d[v] has changed, we have to update the queue
				Q.add(v);
			}
		}
	}
	
	
	/**
	 * Algorithm from the Cormen, Leiserson, Rivest, Stein
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param G
	 * @param s
	 */
	protected static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void initializeSingleSource(HalfEdgeDataStructure<V, E, F> G, V s, double[] d, int[] p) {
		for (V v : G.getVertices()) {
			d[v.getIndex()] = POSITIVE_INFINITY;
			p[v.getIndex()] = -1;
		}
		d[s.getIndex()] = 0;
	}
	
}
